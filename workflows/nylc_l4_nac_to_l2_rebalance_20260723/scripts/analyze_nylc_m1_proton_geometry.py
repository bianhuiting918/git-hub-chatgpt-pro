#!/usr/bin/env python
import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import calc_angles, calc_bonds, distance_array

FAVORABLE = {"da": 3.5, "ha": 2.5, "ang": 135.0}
STRONG = {"da": 3.0, "ha": 2.2, "ang": 150.0}
NAC = {"dist": 3.5, "amin": 95.0, "amax": 115.0}

def sha256(path, chunk=1024*1024):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

def scalar_bond(a, b, box):
    return float(calc_bonds(np.asarray([a]), np.asarray([b]), box=box)[0])

def scalar_angle(d, h, a, box):
    v = float(calc_angles(np.asarray([d]), np.asarray([h]), np.asarray([a]), box=box)[0])
    return math.degrees(v) if math.isfinite(v) else float("nan")

def passes(da, ha, ang, threshold):
    return bool(da <= threshold["da"] and ha <= threshold["ha"] and ang >= threshold["ang"])

def best_direct(donor_pos, hydrogen_positions, acceptor_positions, box):
    best = None
    records = []
    for hi, hpos in enumerate(hydrogen_positions):
        for ai, apos in enumerate(acceptor_positions):
            da = scalar_bond(donor_pos, apos, box)
            ha = scalar_bond(hpos, apos, box)
            ang = scalar_angle(donor_pos, hpos, apos, box)
            rec = {
                "donor_h_index": hi,
                "acceptor_index": ai,
                "d_a_A": da,
                "h_a_A": ha,
                "dha_deg": ang,
                "favorable": passes(da, ha, ang, FAVORABLE),
                "strong": passes(da, ha, ang, STRONG),
            }
            records.append(rec)
            key = (0 if rec["favorable"] else 1, ha, -ang if math.isfinite(ang) else 999.0)
            if best is None or key < best[0]:
                best = (key, rec)
    out = dict(best[1])
    out["favorable"] = any(x["favorable"] for x in records)
    out["strong"] = any(x["strong"] for x in records)
    return out

def best_water_acceptor(donor_pos, hydrogen_positions, water_o, box):
    n_h = len(hydrogen_positions)
    n_w = len(water_o)
    if n_h == 0 or n_w == 0:
        return None, []
    opos = water_o.positions
    hpos = np.asarray(hydrogen_positions)
    ha = distance_array(hpos, opos, box=box)
    da = distance_array(np.asarray([donor_pos]), opos, box=box)[0]
    drep = np.repeat(np.asarray([donor_pos]), n_h*n_w, axis=0)
    hrep = np.repeat(hpos, n_w, axis=0)
    arep = np.tile(opos, (n_h, 1))
    ang = np.degrees(calc_angles(drep, hrep, arep, box=box)).reshape(n_h, n_w)
    fav = (da[None, :] <= FAVORABLE["da"]) & (ha <= FAVORABLE["ha"]) & (ang >= FAVORABLE["ang"])
    strong = (da[None, :] <= STRONG["da"]) & (ha <= STRONG["ha"]) & (ang >= STRONG["ang"])
    candidates = []
    for hi, wi in np.argwhere(fav):
        candidates.append((int(hi), int(wi)))
    rank = np.lexsort((-np.nan_to_num(ang, nan=-999.0).ravel(), ha.ravel(), (~fav).ravel()))
    k = int(rank[0])
    hi, wi = np.unravel_index(k, ha.shape)
    atom = water_o[wi]
    out = {
        "donor_h_index": int(hi),
        "water_array_index": int(wi),
        "water_atom_index1": int(atom.index+1),
        "water_resid": int(atom.resid),
        "d_a_A": float(da[wi]),
        "h_a_A": float(ha[hi, wi]),
        "dha_deg": float(ang[hi, wi]),
        "favorable": bool(np.any(fav)),
        "strong": bool(np.any(strong)),
        "favorable_water_count": int(np.count_nonzero(np.any(fav, axis=0))),
    }
    return out, candidates

def water_bridge(og, hg, asp_o, water_o, water_h_by_o, box):
    first, candidates = best_water_acceptor(og, [hg], water_o, box)
    if first is None or not candidates:
        return {
            "d_a_A": None, "h_a_A": None, "dha_deg": None,
            "favorable": False, "strong": False,
            "water_atom_index1": None, "water_resid": None,
        }
    records = []
    for _, wi in candidates:
        ow = water_o[wi]
        h_atoms = water_h_by_o.get(int(ow.index), [])
        if not h_atoms:
            continue
        second = best_direct(ow.position, [x.position for x in h_atoms], asp_o.positions, box)
        f_da = scalar_bond(og, ow.position, box)
        f_ha = scalar_bond(hg, ow.position, box)
        f_ang = scalar_angle(og, hg, ow.position, box)
        first_strong = passes(f_da, f_ha, f_ang, STRONG)
        rec = {
            "water_atom_index1": int(ow.index+1),
            "water_resid": int(ow.resid),
            "first_d_a_A": f_da,
            "first_h_a_A": f_ha,
            "first_dha_deg": f_ang,
            "second_d_a_A": second["d_a_A"],
            "second_h_a_A": second["h_a_A"],
            "second_dha_deg": second["dha_deg"],
            "d_a_A": max(f_da, second["d_a_A"]),
            "h_a_A": max(f_ha, second["h_a_A"]),
            "dha_deg": min(f_ang, second["dha_deg"]),
            "favorable": bool(second["favorable"]),
            "strong": bool(first_strong and second["strong"]),
        }
        records.append(rec)
    if not records:
        return {
            "d_a_A": None, "h_a_A": None, "dha_deg": None,
            "favorable": False, "strong": False,
            "water_atom_index1": None, "water_resid": None,
        }
    records.sort(key=lambda x: (0 if x["favorable"] else 1, x["h_a_A"], -x["dha_deg"]))
    out = dict(records[0])
    out["favorable"] = any(x["favorable"] for x in records)
    out["strong"] = any(x["strong"] for x in records)
    return out

def finite_stats(values):
    a = np.asarray([x for x in values if x is not None and math.isfinite(float(x))], dtype=float)
    if len(a) == 0:
        return {"n": 0, "min": None, "p05": None, "median": None, "mean": None, "max": None}
    return {
        "n": int(len(a)),
        "min": float(np.min(a)),
        "p05": float(np.percentile(a, 5)),
        "median": float(np.median(a)),
        "mean": float(np.mean(a)),
        "max": float(np.max(a)),
    }

def longest_run(mask, times):
    mask = np.asarray(mask, dtype=bool)
    times = np.asarray(times, dtype=float)
    best = None
    start = None
    for i, val in enumerate(mask):
        if val and start is None:
            start = i
        if start is not None and ((not val) or i == len(mask)-1):
            end = i if val and i == len(mask)-1 else i-1
            count = end-start+1
            duration = float(times[end]-times[start]) if count > 1 else 0.0
            rec = {"start_frame": int(start), "end_frame": int(end), "frame_count": int(count),
                   "start_ps": float(times[start]), "end_ps": float(times[end]), "duration_ps": duration}
            if best is None or (duration, count) > (best["duration_ps"], best["frame_count"]):
                best = rec
            start = None
    return best or {"start_frame": None, "end_frame": None, "frame_count": 0,
                    "start_ps": None, "end_ps": None, "duration_ps": 0.0}

def route_summary(frames, route_name, mask):
    selected = [f for f, keep in zip(frames, mask) if keep]
    vals = [f["routes"][route_name] for f in selected]
    n = len(vals)
    fav = [bool(v and v.get("favorable")) for v in vals]
    strong = [bool(v and v.get("strong")) for v in vals]
    full_joint_fav = [bool(keep and f["routes"][route_name] and f["routes"][route_name].get("favorable"))
                      for f, keep in zip(frames, mask)]
    all_times = [f["time_ps"] for f in frames]
    water_counts = {}
    for v in vals:
        if v and v.get("favorable") and v.get("water_atom_index1") is not None:
            key = str(v["water_atom_index1"])
            water_counts[key] = water_counts.get(key, 0) + 1
    top_water = None
    if water_counts:
        key, count = sorted(water_counts.items(), key=lambda x: (-x[1], int(x[0])))[0]
        top_water = {
            "atom_index1": int(key),
            "favorable_frame_count": int(count),
            "share_of_favorable_frames": float(count/sum(fav)) if sum(fav) else None,
            "occupancy_over_denominator": float(count/n) if n else None,
            "unique_favorable_water_count": int(len(water_counts)),
        }
    return {
        "denominator_frames": n,
        "favorable_frames": int(sum(fav)),
        "favorable_occupancy": float(sum(fav)/n) if n else None,
        "strong_frames": int(sum(strong)),
        "strong_occupancy": float(sum(strong)/n) if n else None,
        "h_a_A": finite_stats([v.get("h_a_A") if v else None for v in vals]),
        "d_a_A": finite_stats([v.get("d_a_A") if v else None for v in vals]),
        "dha_deg": finite_stats([v.get("dha_deg") if v else None for v in vals]),
        "longest_favorable_episode": longest_run(full_joint_fav, all_times) if frames else None,
        "top_favorable_water": top_water,
    }

def sanitize(x):
    if isinstance(x, dict):
        return {k: sanitize(v) for k, v in x.items()}
    if isinstance(x, list):
        return [sanitize(v) for v in x]
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.floating,)):
        x = float(x)
    if isinstance(x, float) and not math.isfinite(x):
        return None
    return x

def analyze(case, outdir):
    tpr = Path(case["tpr"])
    xtc = Path(case["xtc"])
    if not tpr.is_file() or not xtc.is_file():
        raise FileNotFoundError(f"missing input for {case['id']}: {tpr} {xtc}")
    u = mda.Universe(str(tpr), str(xtc))
    atom = lambda index1: u.atoms[int(index1)-1]
    og = atom(case["thr267_og1_index1"])
    carbon_c = atom(case["carbonyl_c_index1"])
    carbon_o = atom(case["carbonyl_o_index1"])
    if og.name != "OG1" or og.resname != "THR":
        raise ValueError(f"{case['id']}: Thr index maps to {og.resname}:{og.name}")
    thr = og.residue
    nalpha = thr.atoms.select_atoms("name N")
    hg = thr.atoms.select_atoms("name HG1")
    nh = thr.atoms.select_atoms("name H1 H2 H3")
    if len(nalpha) != 1 or len(hg) != 1:
        raise ValueError(f"{case['id']}: incomplete Thr267 donor atoms")
    seg_res = list(og.segment.residues)
    thr_pos = next(i for i,r in enumerate(seg_res) if r.ix == thr.ix)
    asp306 = seg_res[thr_pos+39]
    asp308 = seg_res[thr_pos+41]
    if asp306.resname != "ASH" or asp308.resname != "ASP":
        raise ValueError(f"{case['id']}: expected ASH at Thr+39 and ASP at Thr+41, got {asp306.resname}/{asp308.resname}")
    asp306_o = asp306.atoms.select_atoms("name OD1 OD2")
    asp306_h = asp306.atoms.select_atoms("name HD1 HD2")
    asp308_o = asp308.atoms.select_atoms("name OD1 OD2")
    if len(asp306_o) != 2 or len(asp308_o) != 2:
        raise ValueError(f"{case['id']}: ASP oxygen mapping failed")
    water_o = u.select_atoms("resname SOL and name OW")
    if len(water_o) == 0:
        water_o = u.select_atoms("(resname WAT or resname HOH) and (name O or name OH2)")
    water_h_by_o = {}
    for ow in water_o:
        hs = [a for a in ow.residue.atoms if a.name.upper().startswith("H")]
        water_h_by_o[int(ow.index)] = hs
    frames = []
    route_names = [
        "OgH_to_Asp306", "OgH_to_Asp308", "OgH_to_water",
        "NalphaH_to_Asp306", "NalphaH_to_Asp308", "NalphaH_to_water",
        "OgH_to_Nalpha_preorg", "OgH_via_water_to_Asp306", "OgH_via_water_to_Asp308",
    ]
    for ts in u.trajectory:
        box = ts.dimensions
        ogp = og.position.copy()
        hgp = hg[0].position.copy()
        np_ = nalpha[0].position.copy()
        ccp = carbon_c.position.copy()
        cop = carbon_o.position.copy()
        attack_dist = scalar_bond(ogp, ccp, box)
        attack_ang = scalar_angle(cop, ccp, ogp, box)
        nac = bool(attack_dist <= NAC["dist"] and NAC["amin"] <= attack_ang <= NAC["amax"])
        r = {}
        r["OgH_to_Asp306"] = best_direct(ogp, [hgp], asp306_o.positions, box)
        r["OgH_to_Asp308"] = best_direct(ogp, [hgp], asp308_o.positions, box)
        r["OgH_to_water"], _ = best_water_acceptor(ogp, [hgp], water_o, box)
        if len(nh):
            hpos = [x.position.copy() for x in nh]
            r["NalphaH_to_Asp306"] = best_direct(np_, hpos, asp306_o.positions, box)
            r["NalphaH_to_Asp308"] = best_direct(np_, hpos, asp308_o.positions, box)
            r["NalphaH_to_water"], _ = best_water_acceptor(np_, hpos, water_o, box)
        else:
            r["NalphaH_to_Asp306"] = None
            r["NalphaH_to_Asp308"] = None
            r["NalphaH_to_water"] = None
        r["OgH_to_Nalpha_preorg"] = best_direct(ogp, [hgp], np.asarray([np_]), box)
        r["OgH_via_water_to_Asp306"] = water_bridge(ogp, hgp, asp306_o, water_o, water_h_by_o, box)
        r["OgH_via_water_to_Asp308"] = water_bridge(ogp, hgp, asp308_o, water_o, water_h_by_o, box)
        frames.append({
            "frame": int(ts.frame),
            "time_ps": float(ts.time),
            "attack_distance_A": attack_dist,
            "attack_angle_deg": attack_ang,
            "nac": nac,
            "routes": r,
        })
    times = [f["time_ps"] for f in frames]
    all_mask = [True]*len(frames)
    nac_mask = [f["nac"] for f in frames]
    summary = {
        "schema_version": 1,
        "case": case,
        "technical_status": "PASS",
        "scientific_scope": "geometry_screen_only_no_classical_proton_transfer",
        "input_sha256": {"tpr": sha256(tpr), "xtc": sha256(xtc)},
        "atom_mapping": {
            "thr267": {"global_resid": int(thr.resid), "segid": thr.segid,
                       "og1_index1": int(og.index+1), "n_index1": int(nalpha[0].index+1),
                       "hg1_index1": int(hg[0].index+1),
                       "n_hydrogen_names": [x.name for x in nh],
                       "n_hydrogen_indices1": [int(x.index+1) for x in nh]},
            "asp306": {"global_resid": int(asp306.resid), "segid": asp306.segid,
                       "oxygen_indices1": [int(x.index+1) for x in asp306_o]},
            "asp308": {"global_resid": int(asp308.resid), "segid": asp308.segid,
                       "oxygen_indices1": [int(x.index+1) for x in asp308_o]},
            "reactive_carbonyl": {"c_index1": int(carbon_c.index+1), "c_name": carbon_c.name,
                                  "o_index1": int(carbon_o.index+1), "o_name": carbon_o.name,
                                  "resname": carbon_c.resname, "global_resid": int(carbon_c.resid)},
            "validation": "Thr267 first residue; Asp306 is ASH at Thr+39 and Asp308 is ASP at Thr+41",
        },
        "topology_state": {
            "thr267_nalpha_hydrogen_count": int(len(nh)),
            "thr267_nalpha_hydrogen_names": [x.name for x in nh],
            "asp306_proton_names": [x.name for x in asp306_h],
            "microstate_label": "M1_NalphaH2_Asp306H" if len(nh)==2 and asp306.resname=="ASH" and len(asp306_h)==1 else ("M0_like_NalphaH3" if len(nh)==3 else "other"),
            "explicit_water_oxygen_count": int(len(water_o)),
            "fixed_topology_warning": "Classical MD cannot move a proton between molecules or residues.",
        },
        "trajectory": {
            "frame_count": len(frames),
            "time_start_ps": float(times[0]) if times else None,
            "time_end_ps": float(times[-1]) if times else None,
            "median_frame_spacing_ps": float(np.median(np.diff(times))) if len(times)>1 else None,
        },
        "nac": {
            "definition": "OG1-carbonylC <= 3.5 A and O-C-OG1 angle 95-115 deg",
            "frame_count": int(sum(nac_mask)),
            "occupancy": float(sum(nac_mask)/len(frames)) if frames else None,
            "attack_distance_A": finite_stats([f["attack_distance_A"] for f in frames]),
            "attack_angle_deg": finite_stats([f["attack_angle_deg"] for f in frames]),
            "longest_episode": longest_run(nac_mask, times),
        },
        "route_statistics": {
            name: {
                "all_frames": route_summary(frames, name, all_mask),
                "nac_frames": route_summary(frames, name, nac_mask),
            } for name in route_names
        },
        "interpretation_gate": {
            "status": "GEOMETRY_SCREEN_COMPLETE",
            "can_prove_proton_transfer": False,
            "can_rank_qmmm_microstates": True,
            "nac_conditional_status": "EVALUATED" if sum(nac_mask) else "NOT_EVALUATED_NO_NAC_FRAMES",
        },
    }
    cdir = outdir/case["id"]
    cdir.mkdir(parents=True, exist_ok=True)
    with open(cdir/"frames.jsonl", "w", encoding="utf-8") as fh:
        for f in frames:
            fh.write(json.dumps(sanitize(f), sort_keys=True)+"\n")
    summary = sanitize(summary)
    with open(cdir/"summary.json", "w", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)
        fh.write("\n")
    with open(cdir/"frames.tsv", "w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        hdr = ["frame","time_ps","nac","attack_distance_A","attack_angle_deg"]
        for name in route_names:
            hdr += [name+"_fav",name+"_strong",name+"_h_a_A",name+"_dha_deg"]
        w.writerow(hdr)
        for f in frames:
            row=[f["frame"],f["time_ps"],int(f["nac"]),f["attack_distance_A"],f["attack_angle_deg"]]
            for name in route_names:
                x=f["routes"][name]
                row += [int(bool(x and x.get("favorable"))),int(bool(x and x.get("strong"))),
                        x.get("h_a_A") if x else None,x.get("dha_deg") if x else None]
            w.writerow(row)
    return summary

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cases", required=True)
    ap.add_argument("--outdir", required=True)
    args = ap.parse_args()
    cfg = json.loads(Path(args.cases).read_text())
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    combined = {
        "schema_version": 1,
        "thresholds": cfg["thresholds_angstrom_degree"],
        "analysis_scope": cfg["analysis_scope"],
        "excluded_superseded": cfg.get("excluded_superseded", []),
        "results": [],
        "failures": [],
    }
    for case in cfg["cases"]:
        try:
            print("ANALYZE", case["id"], flush=True)
            s = analyze(case, outdir)
            combined["results"].append(s)
            print("PASS", case["id"], flush=True)
        except Exception as e:
            combined["failures"].append({"case_id": case["id"], "type": type(e).__name__, "message": str(e)})
            print("FAIL", case["id"], type(e).__name__, str(e), file=sys.stderr, flush=True)
    combined["technical_status"] = "PASS" if not combined["failures"] else "PARTIAL"
    with open(outdir/"combined_summary.json", "w", encoding="utf-8") as fh:
        json.dump(sanitize(combined), fh, indent=2, sort_keys=True)
        fh.write("\n")
    return 0 if not combined["failures"] else 2

if __name__ == "__main__":
    raise SystemExit(main())

