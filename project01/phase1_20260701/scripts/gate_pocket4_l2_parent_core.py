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
    parser.add_argument("--parent-pdb", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--accepted-tsv", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--n-flank-min", type=int, default=24)
    parser.add_argument("--n-flank-max", type=int, default=36)
    parser.add_argument("--parent-chain", default="A")
    parser.add_argument("--output-chain", default="A")
    parser.add_argument("--ligand", default="bn1")
    parser.add_argument("--rmsd-cutoff", type=float, default=1.0)
    parser.add_argument("--pair-cutoff", type=float, default=1.0)
    return parser.parse_args()


def ca_atoms(path: Path) -> dict[tuple[str, int], np.ndarray]:
    coords: dict[tuple[str, int], np.ndarray] = {}
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        try:
            chain = line[21].strip() or "A"
            resid = int(line[22:26])
            coords[(chain, resid)] = np.array(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                dtype=float,
            )
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


def contiguous_resids(coords: dict[tuple[str, int], np.ndarray], chain: str) -> list[int]:
    return sorted(resid for c, resid in coords if c == chain)


def evaluate_offset(
    ref_coords: dict[tuple[str, int], np.ndarray],
    cand_coords: dict[tuple[str, int], np.ndarray],
    parent_resids: list[int],
    parent_chain: str,
    output_chain: str,
    offset: int,
) -> tuple[float, float, float, int]:
    ref_points: list[np.ndarray] = []
    cand_points: list[np.ndarray] = []
    missing = 0
    for parent_resid in parent_resids:
        ref_key = (parent_chain, parent_resid)
        out_key = (output_chain, parent_resid + offset)
        if ref_key not in ref_coords or out_key not in cand_coords:
            missing += 1
            continue
        ref_points.append(ref_coords[ref_key])
        cand_points.append(cand_coords[out_key])
    if missing or len(ref_points) < 3:
        return math.inf, math.inf, math.inf, missing
    rmsd = kabsch_rmsd(ref_points, cand_points)
    max_delta, mean_delta = pair_delta(ref_points, cand_points)
    return rmsd, max_delta, mean_delta, missing


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    parent_pdb = Path(args.parent_pdb)
    output_dir = Path(args.output_dir)
    parent_ca = ca_atoms(parent_pdb)
    parent_resids = contiguous_resids(parent_ca, args.parent_chain)
    if not parent_resids:
        raise SystemExit(f"No parent CA atoms found for chain {args.parent_chain}: {parent_pdb}")

    rows: list[dict[str, object]] = []
    for pdb in sorted(p for p in output_dir.glob("*.pdb") if p.is_file() and "traj" not in p.name.lower()):
        cand_ca = ca_atoms(pdb)
        best = None
        for n_flank in range(args.n_flank_min, args.n_flank_max + 1):
            offset = n_flank
            rmsd, max_delta, mean_delta, missing = evaluate_offset(
                parent_ca,
                cand_ca,
                parent_resids,
                args.parent_chain,
                args.output_chain,
                offset,
            )
            candidate = (rmsd, max_delta, missing, offset, mean_delta)
            if best is None or candidate < best:
                best = candidate
        if best is None:
            best = (math.inf, math.inf, len(parent_resids), -1, math.inf)
        rmsd, max_delta, missing, offset, mean_delta = best
        ligand_count = ligand_atom_count(pdb, args.ligand)
        reasons: list[str] = []
        if missing:
            reasons.append("missing_parent_core_ca")
        if ligand_count <= 0:
            reasons.append("missing_ligand_records")
        if not math.isfinite(rmsd) or rmsd > args.rmsd_cutoff:
            reasons.append("parent_core_ca_rmsd_gt_cutoff")
        if not math.isfinite(max_delta) or max_delta > args.pair_cutoff:
            reasons.append("parent_core_pair_delta_gt_cutoff")
        rows.append(
            {
                "sample_id": pdb.stem,
                "gate": "PASS" if not reasons else "FAIL",
                "fail_reasons": ";".join(reasons),
                "parent_pdb": str(parent_pdb),
                "pdb": str(pdb),
                "parent_ca_count": len(parent_resids),
                "best_n_flank": offset,
                "parent_core_output_start": parent_resids[0] + offset,
                "parent_core_output_end": parent_resids[-1] + offset,
                "parent_core_ca_rmsd_A": "" if not math.isfinite(rmsd) else f"{rmsd:.4f}",
                "parent_core_pair_max_delta_A": "" if not math.isfinite(max_delta) else f"{max_delta:.4f}",
                "parent_core_pair_mean_delta_A": "" if not math.isfinite(mean_delta) else f"{mean_delta:.4f}",
                "missing_parent_core_ca": missing,
                "ligand_atom_records": ligand_count,
            }
        )

    fields = [
        "sample_id",
        "gate",
        "fail_reasons",
        "parent_pdb",
        "pdb",
        "parent_ca_count",
        "best_n_flank",
        "parent_core_output_start",
        "parent_core_output_end",
        "parent_core_ca_rmsd_A",
        "parent_core_pair_max_delta_A",
        "parent_core_pair_mean_delta_A",
        "missing_parent_core_ca",
        "ligand_atom_records",
    ]
    write_tsv(Path(args.out_tsv), fields, rows)
    accepted = [row for row in rows if row["gate"] == "PASS"]
    write_tsv(Path(args.accepted_tsv), fields, accepted)

    counts: dict[str, int] = {}
    for row in rows:
        counts[str(row["gate"])] = counts.get(str(row["gate"]), 0) + 1
    summary = {
        "status": "DONE",
        "route": "pocket4_l2_parent_core_ca_gate",
        "evaluated_universe": "L2 CA_RFDiffusion output PDB files present in output_dir; absent future samples are NOT_EVALUATED, not FAIL",
        "parent_pdb": str(parent_pdb),
        "output_dir": str(output_dir),
        "evaluated_pdb_count": len(rows),
        "counts": counts,
        "thresholds": {
            "parent_core_ca_rmsd_max_A": args.rmsd_cutoff,
            "parent_core_pair_max_delta_A": args.pair_cutoff,
            "ligand_required": args.ligand,
            "n_flank_scan": [args.n_flank_min, args.n_flank_max],
        },
        "best_passes": accepted[:20],
        "best_fails": [row for row in rows if row["gate"] == "FAIL"][:20],
    }
    Path(args.summary_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary_json).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
