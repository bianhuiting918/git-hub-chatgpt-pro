#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from pathlib import Path

import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest-json", required=True)
    parser.add_argument("--contig-set", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--rmsd-cutoff", type=float, default=1.0)
    parser.add_argument("--pair-cutoff", type=float, default=1.0)
    parser.add_argument("--ligand", default="bn1")
    return parser.parse_args()


def ca_atoms(path: Path) -> dict[tuple[str, int], np.ndarray]:
    coords: dict[tuple[str, int], np.ndarray] = {}
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        try:
            chain = line[21].strip() or "A"
            resid = int(line[22:26])
            coords[(chain, resid)] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float)
        except ValueError:
            continue
    return coords


def ligand_atom_count(path: Path, ligand: str) -> int:
    ligand_l = ligand.lower()
    count = 0
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("HETATM") and line[17:20].strip().lower() == ligand_l and not line[12:16].strip().startswith("H"):
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
        values.append(abs(float(np.linalg.norm(ref_points[i] - ref_points[j]) - np.linalg.norm(cand_points[i] - cand_points[j]))))
    return max(values), float(np.mean(values))


def mapping_for_contig(manifest: dict[str, object], contig_set: str) -> list[dict[str, object]]:
    mapping_by_contig = manifest.get("mapping_by_contig", {})
    if not isinstance(mapping_by_contig, dict) or contig_set not in mapping_by_contig:
        raise ValueError(f"Contig set not found in manifest mapping: {contig_set}")
    mapping = mapping_by_contig[contig_set]
    if not isinstance(mapping, list):
        raise ValueError(f"Invalid mapping for contig set: {contig_set}")
    return mapping


def evaluate(
    pdb: Path,
    ref_ca: dict[tuple[str, int], np.ndarray],
    mapping: list[dict[str, object]],
    ligand: str,
    rmsd_cutoff: float,
    pair_cutoff: float,
) -> dict[str, object]:
    cand_ca = ca_atoms(pdb)
    ref_points: list[np.ndarray] = []
    cand_points: list[np.ndarray] = []
    missing_ref: list[str] = []
    missing_out: list[str] = []
    for item in mapping:
        ref_key = (str(item["reference_chain"]), int(item["reference_resid"]))
        out_key = (str(item["output_chain"]), int(item["output_resid"]))
        if ref_key not in ref_ca:
            missing_ref.append(f"{ref_key[0]}{ref_key[1]}")
            continue
        if out_key not in cand_ca:
            missing_out.append(f"{out_key[0]}{out_key[1]}")
            continue
        ref_points.append(ref_ca[ref_key])
        cand_points.append(cand_ca[out_key])

    reasons: list[str] = []
    if missing_ref:
        reasons.append("missing_reference_pocket_ca")
    if missing_out:
        reasons.append("missing_output_pocket_ca")
    lig_count = ligand_atom_count(pdb, ligand)
    if lig_count <= 0:
        reasons.append("missing_ligand_records")

    rmsd = pair_max = pair_mean = None
    if not missing_ref and not missing_out and ref_points:
        rmsd = kabsch_rmsd(ref_points, cand_points)
        pair_max, pair_mean = pair_delta(ref_points, cand_points)
        if rmsd > rmsd_cutoff:
            reasons.append("pocket_ca_rmsd_gt_cutoff")
        if pair_max > pair_cutoff:
            reasons.append("pocket_pair_delta_gt_cutoff")

    return {
        "sample_id": pdb.stem,
        "gate": "PASS" if not reasons else "FAIL",
        "fail_reasons": ";".join(reasons),
        "pocket_ca_count": len(mapping),
        "pocket_ca_rmsd_A": "" if rmsd is None else f"{rmsd:.4f}",
        "pocket_pair_max_delta_A": "" if pair_max is None else f"{pair_max:.4f}",
        "pocket_pair_mean_delta_A": "" if pair_mean is None else f"{pair_mean:.4f}",
        "ligand_atom_records": lig_count,
        "missing_reference_pocket_ca": ",".join(missing_ref),
        "missing_output_pocket_ca": ",".join(missing_out),
        "pdb": str(pdb),
    }


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    manifest = json.loads(Path(args.manifest_json).read_text())
    reference_pdb = Path(str(manifest["input_pdb"]))
    ref_ca = ca_atoms(reference_pdb)
    mapping = mapping_for_contig(manifest, args.contig_set)
    pdbs = sorted(p for p in Path(args.output_dir).glob("*.pdb") if p.is_file() and "traj" not in p.name.lower())
    rows = [evaluate(pdb, ref_ca, mapping, args.ligand, args.rmsd_cutoff, args.pair_cutoff) for pdb in pdbs]

    fields = [
        "sample_id",
        "gate",
        "fail_reasons",
        "pocket_ca_count",
        "pocket_ca_rmsd_A",
        "pocket_pair_max_delta_A",
        "pocket_pair_mean_delta_A",
        "ligand_atom_records",
        "missing_reference_pocket_ca",
        "missing_output_pocket_ca",
        "pdb",
    ]
    write_tsv(Path(args.out_tsv), rows, fields)

    counts: dict[str, int] = {}
    for row in rows:
        counts[str(row["gate"])] = counts.get(str(row["gate"]), 0) + 1
    passes = sorted((row for row in rows if row["gate"] == "PASS"), key=lambda row: float(row["pocket_ca_rmsd_A"]))
    fails = sorted(
        (row for row in rows if row["gate"] == "FAIL" and row["pocket_ca_rmsd_A"]),
        key=lambda row: (float(row["pocket_ca_rmsd_A"]), float(row["pocket_pair_max_delta_A"])),
    )
    summary = {
        "status": "DONE",
        "route": "pocket4_first_layered_ca_rfdiffusion_gate",
        "evaluated_universe": "CA_RFDiffusion output PDB files present in output_dir for this pocket4 route; absent future samples are NOT_EVALUATED, not FAIL",
        "manifest_json": args.manifest_json,
        "reference_pdb": str(reference_pdb),
        "output_dir": args.output_dir,
        "contig_set": args.contig_set,
        "evaluated_pdb_count": len(rows),
        "counts": counts,
        "best_passes": passes[:20],
        "best_fails": fails[:20],
        "thresholds": {
            "pocket_ca_rmsd_max_A": args.rmsd_cutoff,
            "pocket_pair_max_delta_A": args.pair_cutoff,
            "ligand_required": args.ligand,
        },
        "pocket_mapping": mapping,
    }
    Path(args.summary_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary_json).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
