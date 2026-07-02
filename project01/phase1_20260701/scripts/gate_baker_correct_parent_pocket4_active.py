#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

CATALYTIC_AND_OXY = {15, 16, 55, 73, 149}
CATALYTIC_ATOM_PAIRS = [
    ("ser15_og_to_his55_ne2", (15, "OG"), (55, "NE2")),
    ("his55_nd1_to_asp73_od1", (55, "ND1"), (73, "OD1")),
    ("his55_nd1_to_asp73_od2", (55, "ND1"), (73, "OD2")),
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open(), delimiter="\t")) if path.exists() else []


def heavy_xyz(line: str) -> tuple[float, float, float]:
    return float(line[30:38]), float(line[38:46]), float(line[46:54])


def ligand_4a_fixed_residues(path: Path, cutoff: float) -> set[int]:
    protein = []
    ligand = []
    for line in path.read_text(errors="ignore").splitlines():
        rec = line[:6].strip()
        atom = line[12:16].strip()
        if atom.startswith("H"):
            continue
        if rec == "ATOM" and (line[21].strip() or "A") == "A":
            protein.append((int(line[22:26]), heavy_xyz(line)))
        elif rec == "HETATM" and line[17:20].strip() == "bn1":
            ligand.append(heavy_xyz(line))
    fixed = set(CATALYTIC_AND_OXY)
    for resid, xyz in protein:
        if any(math.dist(xyz, lig_xyz) <= cutoff for lig_xyz in ligand):
            fixed.add(resid)
    return fixed


def parse_atoms(path: Path):
    ca = {}
    names = {}
    atoms = {}
    if not path.exists():
        return ca, names, atoms
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM") or (line[21].strip() or "A") != "A":
            continue
        resid = int(line[22:26])
        atom = line[12:16].strip()
        resn = line[17:20].strip()
        xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
        names[resid] = resn
        atoms[(resid, atom)] = xyz
        if atom == "CA":
            ca[resid] = xyz
    return ca, names, atoms


def kabsch_rmsd(ref_points: list[np.ndarray], cand_points: list[np.ndarray]) -> float:
    ref = np.array(ref_points)
    cand = np.array(cand_points)
    ref_centered = ref - ref.mean(0)
    cand_centered = cand - cand.mean(0)
    u, _s, vt = np.linalg.svd(cand_centered.T @ ref_centered)
    rot = u @ vt
    if np.linalg.det(rot) < 0:
        u[:, -1] *= -1
        rot = u @ vt
    aligned = cand_centered @ rot
    return float(np.sqrt(((ref_centered - aligned) ** 2).sum() / len(ref)))


def pair_delta(ref: dict[int, np.ndarray], cand: dict[int, np.ndarray], fixed: list[int]) -> tuple[float, float]:
    values = []
    for i, j in itertools.combinations(fixed, 2):
        values.append(abs(np.linalg.norm(ref[i] - ref[j]) - np.linalg.norm(cand[i] - cand[j])))
    return max(values), float(np.mean(values))


def catalytic_distance_deltas(ref_atoms, cand_atoms) -> tuple[float | str, str]:
    deltas = []
    missing = []
    for label, a, b in CATALYTIC_ATOM_PAIRS:
        if a not in ref_atoms or b not in ref_atoms or a not in cand_atoms or b not in cand_atoms:
            missing.append(label)
            continue
        ref_d = np.linalg.norm(ref_atoms[a] - ref_atoms[b])
        cand_d = np.linalg.norm(cand_atoms[a] - cand_atoms[b])
        deltas.append(abs(ref_d - cand_d))
    if missing:
        return "", "missing_catalytic_atoms:" + ",".join(missing)
    return max(deltas), ""


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref-pdb", required=True)
    parser.add_argument("--ligand-ref-pdb", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--status-tsv", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--pocket-cutoff", type=float, default=4.0)
    parser.add_argument("--ca-rmsd-cutoff", type=float, default=1.0)
    parser.add_argument("--pair-cutoff", type=float, default=1.0)
    parser.add_argument("--cat-distance-delta-cutoff", type=float, default=1.0)
    args = parser.parse_args()

    fixed = sorted(ligand_4a_fixed_residues(Path(args.ligand_ref_pdb), args.pocket_cutoff))
    ref_ca, ref_names, ref_atoms = parse_atoms(Path(args.ref_pdb))
    manifest_rows = read_tsv(Path(args.manifest))
    status_by_id = {row["sample_id"]: row for row in read_tsv(Path(args.status_tsv))}
    out_rows = []

    for row in manifest_rows:
        sid = row["sample_id"]
        status = status_by_id.get(sid, {})
        pdb = Path(row["postseq_predicted_pdb"])
        reasons = []
        rmsd = ""
        pmax = ""
        pmean = ""
        cat_delta = ""
        if status.get("status") not in {"OK", "SKIP_EXISTS"} or not pdb.exists():
            gate = "NOT_EVALUATED"
            reasons.append("missing_or_pending_esmfold")
        else:
            cand_ca, cand_names, cand_atoms = parse_atoms(pdb)
            missing = [resid for resid in fixed if resid not in cand_ca or resid not in ref_ca]
            fixed_mut = [
                f"A{resid}:{ref_names.get(resid)}>{cand_names.get(resid)}"
                for resid in fixed
                if cand_names.get(resid) != ref_names.get(resid)
            ]
            if missing:
                reasons.append("missing_fixed_ca")
            if fixed_mut:
                reasons.append("fixed_pocket_identity_mutation")
            if not missing:
                rmsd = kabsch_rmsd([ref_ca[resid] for resid in fixed], [cand_ca[resid] for resid in fixed])
                pmax, pmean = pair_delta(ref_ca, cand_ca, fixed)
                if rmsd > args.ca_rmsd_cutoff:
                    reasons.append("pocket_ca_rmsd_gt_cutoff")
                if pmax > args.pair_cutoff:
                    reasons.append("pocket_pair_delta_gt_cutoff")
            cat_delta, cat_missing = catalytic_distance_deltas(ref_atoms, cand_atoms)
            if cat_missing:
                reasons.append(cat_missing)
            elif cat_delta > args.cat_distance_delta_cutoff:
                reasons.append("catalytic_sidechain_distance_delta_gt_cutoff")
            gate = "FINAL_QUALIFIED_ACTIVE" if not reasons else "FAIL"
        out_rows.append({
            "sample_id": sid,
            "bin": row.get("bin", ""),
            "identity": row.get("identity", ""),
            "mutation_count": row.get("mutation_count", ""),
            "esmfold_status": status.get("status", "NOT_EVALUATED"),
            "gate": gate,
            "fail_reasons": ";".join(reasons),
            "pocket_ca_rmsd_A": "" if rmsd == "" else f"{rmsd:.4f}",
            "pocket_pair_max_delta_A": "" if pmax == "" else f"{pmax:.4f}",
            "pocket_pair_mean_delta_A": "" if pmean == "" else f"{pmean:.4f}",
            "catalytic_sidechain_max_delta_A": "" if cat_delta == "" else f"{cat_delta:.4f}",
            "pdb": str(pdb),
            "sequence": row["sequence"],
        })

    out_path = Path(args.out_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(out_rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(out_rows)

    by_bin = defaultdict(Counter)
    best_by_bin = defaultdict(list)
    for row in out_rows:
        by_bin[row["bin"]][row["gate"]] += 1
        if row["gate"] == "FINAL_QUALIFIED_ACTIVE":
            best_by_bin[row["bin"]].append(row)
    for rows in best_by_bin.values():
        rows.sort(key=lambda item: (float(item["pocket_ca_rmsd_A"]), float(item["catalytic_sidechain_max_delta_A"])))

    summary = {
        "status": "DONE",
        "evaluated_universe": {
            "manifest_rows": len(manifest_rows),
            "fixed_residue_count": len(fixed),
            "fixed_residues": fixed,
        },
        "gate_definition": (
            "FINAL_QUALIFIED_ACTIVE if ESMFold OK/SKIP_EXISTS, 4A ligand-pocket fixed identities retained, "
            f"pocket CA RMSD <= {args.ca_rmsd_cutoff} A, pocket pair max delta <= {args.pair_cutoff} A, "
            f"and catalytic sidechain distance max delta <= {args.cat_distance_delta_cutoff} A. "
            "This is the sequence-phase active-pocket gate, not a PLACER/QMMM success label."
        ),
        "counts_by_bin": {bin_id: dict(counts) for bin_id, counts in sorted(by_bin.items())},
        "best_passes_by_bin": {bin_id: rows[:10] for bin_id, rows in sorted(best_by_bin.items())},
        "out_tsv": str(out_path),
    }
    Path(args.summary_json).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
