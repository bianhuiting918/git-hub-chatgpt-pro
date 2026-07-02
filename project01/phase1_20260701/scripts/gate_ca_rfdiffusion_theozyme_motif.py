#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import itertools
import json
from pathlib import Path

import numpy as np

# Baker design_pipeline contig used for Project 01 formal new-backbone route:
# 12,A56-60,36,A83-85,15,A113-115,73,B145-147,10
# Output positions are sequential after each generated segment.
MOTIF_MAP = [
    (("A", 56), 13), (("A", 57), 14), (("A", 58), 15), (("A", 59), 16), (("A", 60), 17),
    (("A", 83), 54), (("A", 84), 55), (("A", 85), 56),
    (("A", 113), 72), (("A", 114), 73), (("A", 115), 74),
    (("B", 145), 148), (("B", 146), 149), (("B", 147), 150),
]


def ca_atoms(path: Path) -> dict[tuple[str, int], np.ndarray]:
    coords: dict[tuple[str, int], np.ndarray] = {}
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        chain = line[21].strip() or "A"
        try:
            resid = int(line[22:26])
            coords[(chain, resid)] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float)
        except ValueError:
            continue
    return coords


def ligand_count(path: Path, ligand: str) -> int:
    ligand_lower = ligand.lower()
    count = 0
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("HETATM") and line[17:20].strip().lower() == ligand_lower:
            count += 1
    return count


def kabsch_rmsd(ref_points: list[np.ndarray], cand_points: list[np.ndarray]) -> float:
    ref = np.array(ref_points, dtype=float)
    cand = np.array(cand_points, dtype=float)
    ref0 = ref - ref.mean(axis=0)
    cand0 = cand - cand.mean(axis=0)
    u, _s, vt = np.linalg.svd(cand0.T @ ref0)
    rot = u @ vt
    if np.linalg.det(rot) < 0:
        u[:, -1] *= -1
        rot = u @ vt
    aligned = cand0 @ rot
    return float(np.sqrt(((ref0 - aligned) ** 2).sum() / len(ref)))


def pair_delta(ref_points: list[np.ndarray], cand_points: list[np.ndarray]) -> tuple[float, float]:
    values = []
    for i, j in itertools.combinations(range(len(ref_points)), 2):
        ref_dist = np.linalg.norm(ref_points[i] - ref_points[j])
        cand_dist = np.linalg.norm(cand_points[i] - cand_points[j])
        values.append(abs(ref_dist - cand_dist))
    return float(max(values)), float(np.mean(values))


def evaluate_pdb(ref_coords: dict[tuple[str, int], np.ndarray], pdb: Path, ligand: str, rmsd_cutoff: float, pair_cutoff: float) -> dict[str, str]:
    cand_coords = ca_atoms(pdb)
    missing_ref = []
    missing_out = []
    ref_points = []
    cand_points = []
    for ref_key, out_resid in MOTIF_MAP:
        out_key = ("A", out_resid)
        if ref_key not in ref_coords:
            missing_ref.append(f"{ref_key[0]}{ref_key[1]}")
            continue
        if out_key not in cand_coords:
            missing_out.append(f"A{out_resid}")
            continue
        ref_points.append(ref_coords[ref_key])
        cand_points.append(cand_coords[out_key])

    lig_count = ligand_count(pdb, ligand)
    reasons = []
    if missing_ref:
        reasons.append("missing_reference_motif_ca")
    if missing_out:
        reasons.append("missing_output_motif_ca")
    if lig_count <= 0:
        reasons.append("missing_ligand_records")

    if reasons:
        return {
            "gate": "FAIL",
            "fail_reasons": ";".join(reasons),
            "motif_ca_rmsd_A": "",
            "motif_pair_max_delta_A": "",
            "motif_pair_mean_delta_A": "",
            "ligand_atom_records": str(lig_count),
            "missing_reference_motif_ca": ",".join(missing_ref),
            "missing_output_motif_ca": ",".join(missing_out),
        }

    rmsd = kabsch_rmsd(ref_points, cand_points)
    pair_max, pair_mean = pair_delta(ref_points, cand_points)
    if rmsd > rmsd_cutoff:
        reasons.append("motif_ca_rmsd_gt_cutoff")
    if pair_max > pair_cutoff:
        reasons.append("motif_pair_delta_gt_cutoff")
    return {
        "gate": "PASS" if not reasons else "FAIL",
        "fail_reasons": ";".join(reasons),
        "motif_ca_rmsd_A": f"{rmsd:.4f}",
        "motif_pair_max_delta_A": f"{pair_max:.4f}",
        "motif_pair_mean_delta_A": f"{pair_mean:.4f}",
        "ligand_atom_records": str(lig_count),
        "missing_reference_motif_ca": "",
        "missing_output_motif_ca": "",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-pdb", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--ligand", default="bn1")
    parser.add_argument("--rmsd-cutoff", type=float, default=1.0)
    parser.add_argument("--pair-cutoff", type=float, default=1.0)
    args = parser.parse_args()

    reference = Path(args.reference_pdb)
    output_dir = Path(args.output_dir)
    out_tsv = Path(args.out_tsv)
    summary_json = Path(args.summary_json)

    ref_coords = ca_atoms(reference)
    pdbs = sorted(p for p in output_dir.glob("*.pdb") if p.is_file() and "traj" not in p.name.lower())
    rows = []
    for pdb in pdbs:
        row = {
            "sample_id": pdb.stem,
            "pdb": str(pdb),
            "gate_definition": f"14 motif CA Kabsch RMSD <= {args.rmsd_cutoff} A and max pair-distance delta <= {args.pair_cutoff} A; ligand {args.ligand} present",
        }
        row.update(evaluate_pdb(ref_coords, pdb, args.ligand, args.rmsd_cutoff, args.pair_cutoff))
        rows.append(row)

    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    fields = [
        "sample_id", "gate", "fail_reasons", "motif_ca_rmsd_A", "motif_pair_max_delta_A",
        "motif_pair_mean_delta_A", "ligand_atom_records", "missing_reference_motif_ca",
        "missing_output_motif_ca", "pdb", "gate_definition",
    ]
    with out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    counts: dict[str, int] = {}
    for row in rows:
        counts[row["gate"]] = counts.get(row["gate"], 0) + 1
    passes = sorted((row for row in rows if row["gate"] == "PASS"), key=lambda row: float(row["motif_ca_rmsd_A"]))
    fails = sorted((row for row in rows if row["gate"] == "FAIL" and row["motif_ca_rmsd_A"]), key=lambda row: (float(row["motif_ca_rmsd_A"]), float(row["motif_pair_max_delta_A"])))
    summary = {
        "status": "DONE",
        "evaluated_universe": "CA_RFDiffusion output PDB files present in output_dir; absent future samples are NOT_EVALUATED, not FAIL",
        "reference_pdb": str(reference),
        "output_dir": str(output_dir),
        "out_tsv": str(out_tsv),
        "evaluated_pdb_count": len(rows),
        "counts": counts,
        "best_passes": passes[:20],
        "best_fails": fails[:20],
        "motif_map": [{"reference": f"{chain}{resid}", "output": f"A{out_resid}"} for (chain, resid), out_resid in MOTIF_MAP],
        "thresholds": {"motif_ca_rmsd_max_A": args.rmsd_cutoff, "motif_pair_max_delta_A": args.pair_cutoff},
    }
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
