#!/usr/bin/env python3
import argparse
import concurrent.futures
import csv
import json
import math
import subprocess
from collections import defaultdict
from pathlib import Path

TRUE_THR267_OG1 = 8961
WRONG_THR262_OG1 = 8896


def gro_atom(path: Path, atom_number: int):
    for line in path.read_text().splitlines()[2:-1]:
        try:
            nr = int(line[15:20])
        except ValueError:
            continue
        if nr == atom_number:
            return {
                "atom_number": nr,
                "residue_number": int(line[:5]),
                "residue_name": line[5:10].strip(),
                "atom_name": line[10:15].strip(),
            }
    raise ValueError(f"atom {atom_number} not found in {path}")


def require_identity(path: Path, atom_number: int, expected):
    atom = gro_atom(path, atom_number)
    observed = (atom["residue_number"], atom["residue_name"], atom["atom_name"])
    if observed != expected:
        raise ValueError(f"identity mismatch for {atom_number}: {observed} != {expected}")
    return atom


def read_xvg(path: Path):
    rows = []
    for line in path.read_text().splitlines():
        if not line.strip() or line[0] in "#@":
            continue
        parts = line.split()
        rows.append((float(parts[0]), float(parts[1])))
    return rows


def run_segment(task):
    seg, rows, directory, gmx, c_atom, o_atom, work = task
    xtc = directory / "run.xtc"
    tpr = directory / "run.tpr"
    if not xtc.is_file() or not tpr.is_file():
        raise FileNotFoundError(f"missing trajectory input for {seg}: {directory}")
    dist = work / f"{seg}.distance.xvg"
    angle = work / f"{seg}.angle.xvg"
    ndx = work / f"{seg}.angle.ndx"
    ndx.write_text(f"[ true_thr267_angle ]\n{TRUE_THR267_OG1} {c_atom} {o_atom}\n")
    subprocess.run(
        [gmx, "distance", "-f", str(xtc), "-s", str(tpr),
         "-select", f"atomnr {TRUE_THR267_OG1} plus atomnr {c_atom}",
         "-oall", str(dist)],
        check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    subprocess.run(
        [gmx, "angle", "-f", str(xtc), "-n", str(ndx),
         "-ov", str(angle), "-type", "angle"],
        input="0\n", text=True, check=True,
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
    )
    dvals = read_xvg(dist)
    avals = read_xvg(angle)
    if len(dvals) != len(avals):
        raise ValueError(f"distance/angle frame mismatch for {seg}")
    if len(rows) == len(dvals):
        paired = list(zip(rows, dvals, avals))
    elif len(rows) == len(dvals) - 1:
        paired = list(zip(rows, dvals[1:], avals[1:]))
    else:
        raise ValueError(f"geometry/XTC frame mismatch for {seg}: {len(rows)} vs {len(dvals)}")
    out = []
    for source, (_, distance_nm), (_, angle_deg) in paired:
        record = dict(source)
        record["true_thr267_distance_nm"] = distance_nm
        record["true_thr267_angle_deg"] = angle_deg
        record["true_joint_pass"] = int(
            int(record["substrate_restrained"]) == 0
            and int(record["gate_restrained"]) == 0
            and distance_nm <= 0.35
            and 95.0 <= angle_deg <= 115.0
        )
        out.append(record)
    return out


def longest_continuous_nac(rows):
    runs, current = [], []
    for row in rows:
        if int(row["true_joint_pass"]) == 1:
            if current and float(row["simulation_time_ps"]) - float(current[-1]["simulation_time_ps"]) > 2.01:
                runs.append(current)
                current = []
            current.append(row)
        elif current:
            runs.append(current)
            current = []
    if current:
        runs.append(current)
    if not runs:
        return None
    run = max(runs, key=lambda x: (len(x), float(x[-1]["simulation_time_ps"]) - float(x[0]["simulation_time_ps"])))
    return {
        "start_ps": float(run[0]["simulation_time_ps"]),
        "end_ps": float(run[-1]["simulation_time_ps"]),
        "duration_ps": float(run[-1]["simulation_time_ps"]) - float(run[0]["simulation_time_ps"]),
        "frame_count": len(run),
        "segments": sorted({x["segment"] for x in run}),
        "mean_distance_nm": sum(float(x["true_thr267_distance_nm"]) for x in run) / len(run),
        "mean_angle_deg": sum(float(x["true_thr267_angle_deg"]) for x in run) / len(run),
    }


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--branch-root", required=True, type=Path)
    p.add_argument("--carbonyl-c", required=True, type=int)
    p.add_argument("--carbonyl-o", required=True, type=int)
    p.add_argument("--gmx", required=True)
    p.add_argument("--output-dir", required=True, type=Path)
    p.add_argument("--workers", type=int, default=4)
    args = p.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    work = args.output_dir / "segment_xvg"
    work.mkdir(exist_ok=True)

    true_atom = require_identity(args.branch_root / "source.gro", TRUE_THR267_OG1, expected=(267, "THR", "OG1"))
    wrong_atom = require_identity(args.branch_root / "source.gro", WRONG_THR262_OG1, expected=(262, "THR", "OG1"))

    grouped = defaultdict(list)
    with (args.branch_root / "geometry.tsv").open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            grouped[row["segment"]].append(row)
    tasks = []
    for seg, rows in sorted(grouped.items()):
        matches = sorted((args.branch_root / "segments").glob(f"{seg}_*"))
        if len(matches) != 1:
            raise ValueError(f"expected one directory for {seg}, found {len(matches)}")
        tasks.append((seg, rows, matches[0], args.gmx, args.carbonyl_c, args.carbonyl_o, work))

    corrected = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
        for chunk in pool.map(run_segment, tasks):
            corrected.extend(chunk)
    corrected.sort(key=lambda x: float(x["simulation_time_ps"]))

    fields = list(corrected[0]) if corrected else []
    with (args.output_dir / "true_thr267_geometry.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(corrected)

    free = [x for x in corrected if int(x["substrate_restrained"]) == 0 and int(x["gate_restrained"]) == 0]
    nac = [x for x in free if int(x["true_joint_pass"]) == 1]
    closest = min(free, key=lambda x: float(x["true_thr267_distance_nm"])) if free else None
    angle_compatible = [x for x in free if 95.0 <= float(x["true_thr267_angle_deg"]) <= 115.0]
    closest_angle_compatible = min(angle_compatible, key=lambda x: float(x["true_thr267_distance_nm"])) if angle_compatible else None
    result = {
        "schema_version": 1,
        "status": "CORRECTED_TRUE_THR267_RESIDENCE_AUDIT",
        "branch": args.branch_root.name,
        "identity": {
            "true_catalytic_atom": true_atom,
            "superseded_wrong_atom": wrong_atom,
            "legacy_results_status": "SUPERSEDED_WRONG_NUCLEOPHILE",
        },
        "reaction_atoms": {
            "thr267_og1": TRUE_THR267_OG1,
            "carbonyl_c": args.carbonyl_c,
            "carbonyl_o": args.carbonyl_o,
        },
        "frame_counts": {
            "all": len(corrected),
            "fully_unrestrained": len(free),
            "true_nac": len(nac),
        },
        "true_nac_occupancy_unrestrained": (len(nac) / len(free)) if free else None,
        "longest_continuous_nac": longest_continuous_nac(corrected),
        "closest_unrestrained_frame": closest,
        "closest_angle_compatible_unrestrained_frame": closest_angle_compatible,
        "scientific_status": "TRUE_NAC_FOUND" if nac else "FAIL_NO_TRUE_THR267_NAC_IN_EXISTING_BRANCH",
    }
    (args.output_dir / "true_thr267_residence_audit.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
