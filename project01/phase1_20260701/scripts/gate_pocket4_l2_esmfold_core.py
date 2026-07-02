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
    parser.add_argument("--selected-tsv", required=True)
    parser.add_argument("--esmfold-status-tsv", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--accepted-tsv", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--rmsd-cutoff", type=float, default=1.0)
    parser.add_argument("--pair-cutoff", type=float, default=1.0)
    parser.add_argument("--plddt-cutoff", type=float, default=70.0)
    parser.add_argument("--target-per-bin", type=int, default=10)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    return list(csv.DictReader(path.open(), delimiter="\t"))


def ca_atoms(path: Path) -> dict[tuple[str, int], tuple[np.ndarray, float]]:
    coords: dict[tuple[str, int], tuple[np.ndarray, float]] = {}
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM") or line[12:16].strip() != "CA":
            continue
        try:
            chain = line[21].strip() or "A"
            resid = int(line[22:26])
            xyz = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float)
            bfac = float(line[60:66])
            coords[(chain, resid)] = (xyz, bfac)
        except ValueError:
            continue
    return coords


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


def evaluate(row: dict[str, str], status_by_sample: dict[str, dict[str, str]], args: argparse.Namespace) -> dict[str, object]:
    sample_id = row["sample_id"]
    status = status_by_sample.get(sample_id)
    base = {
        "sample_id": sample_id,
        "parent_sample_id": row.get("parent_sample_id", ""),
        "bin": row.get("bin", ""),
        "identity": row.get("identity", ""),
        "mutation_count": row.get("mutation_count", ""),
        "gate": "NOT_EVALUATED",
        "fail_reasons": "missing_or_failed_esmfold",
        "core_ca_count": "",
        "core_ca_rmsd_A": "",
        "core_pair_max_delta_A": "",
        "core_pair_mean_delta_A": "",
        "core_mean_plddt": "",
        "predicted_pdb": "",
        "parent_pdb": row.get("parent_pdb", ""),
        "fixed_core_start": row.get("fixed_core_start", ""),
        "fixed_core_end": row.get("fixed_core_end", ""),
    }
    if not status or status.get("status") not in {"OK", "SKIP_EXISTS"}:
        return base

    parent_pdb = Path(row["parent_pdb"])
    pred_pdb = Path(status["pdb"])
    if not parent_pdb.exists() or not pred_pdb.exists():
        return base
    start = int(row["fixed_core_start"])
    end = int(row["fixed_core_end"])
    parent_ca = ca_atoms(parent_pdb)
    pred_ca = ca_atoms(pred_pdb)
    ref_points: list[np.ndarray] = []
    cand_points: list[np.ndarray] = []
    plddts: list[float] = []
    missing: list[str] = []
    for resid in range(start, end + 1):
        key = ("A", resid)
        if key not in parent_ca or key not in pred_ca:
            missing.append(f"A{resid}")
            continue
        ref_points.append(parent_ca[key][0])
        cand_points.append(pred_ca[key][0])
        plddts.append(pred_ca[key][1])

    reasons: list[str] = []
    if missing:
        reasons.append("missing_core_ca")
    rmsd = pair_max = pair_mean = mean_plddt = None
    if not missing and ref_points:
        rmsd = kabsch_rmsd(ref_points, cand_points)
        pair_max, pair_mean = pair_delta(ref_points, cand_points)
        mean_plddt = float(np.mean(plddts))
        if rmsd > args.rmsd_cutoff:
            reasons.append("core_ca_rmsd_gt_cutoff")
        if pair_max > args.pair_cutoff:
            reasons.append("core_pair_delta_gt_cutoff")
        if mean_plddt < args.plddt_cutoff:
            reasons.append("core_plddt_lt_cutoff")

    base.update(
        {
            "gate": "PASS" if not reasons else "FAIL",
            "fail_reasons": ";".join(reasons),
            "core_ca_count": end - start + 1,
            "core_ca_rmsd_A": "" if rmsd is None else f"{rmsd:.4f}",
            "core_pair_max_delta_A": "" if pair_max is None else f"{pair_max:.4f}",
            "core_pair_mean_delta_A": "" if pair_mean is None else f"{pair_mean:.4f}",
            "core_mean_plddt": "" if mean_plddt is None else f"{mean_plddt:.2f}",
            "predicted_pdb": str(pred_pdb),
        }
    )
    return base


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    selected_rows = read_tsv(Path(args.selected_tsv))
    status_rows = read_tsv(Path(args.esmfold_status_tsv))
    status_by_sample = {row["sample_id"]: row for row in status_rows}
    rows = [evaluate(row, status_by_sample, args) for row in selected_rows]
    fields = [
        "sample_id", "parent_sample_id", "bin", "identity", "mutation_count", "gate",
        "fail_reasons", "core_ca_count", "core_ca_rmsd_A", "core_pair_max_delta_A",
        "core_pair_mean_delta_A", "core_mean_plddt", "predicted_pdb", "parent_pdb",
        "fixed_core_start", "fixed_core_end",
    ]
    write_tsv(Path(args.out_tsv), fields, rows)
    accepted = [row for row in rows if row["gate"] == "PASS"]
    write_tsv(Path(args.accepted_tsv), fields, accepted)

    counts_by_bin: dict[str, dict[str, int]] = {}
    for row in rows:
        bin_id = str(row["bin"])
        gate = str(row["gate"])
        counts_by_bin.setdefault(bin_id, {"PASS": 0, "FAIL": 0, "NOT_EVALUATED": 0})
        counts_by_bin[bin_id][gate] = counts_by_bin[bin_id].get(gate, 0) + 1
    target_complete = all(counts_by_bin.get(bin_id, {}).get("PASS", 0) >= args.target_per_bin for bin_id in ["90", "80", "70", "60", "50"])
    summary = {
        "status": "DONE",
        "route": "Route B L2 ESMFold parent-core gate",
        "selected_tsv": args.selected_tsv,
        "esmfold_status_tsv": args.esmfold_status_tsv,
        "out_tsv": args.out_tsv,
        "accepted_tsv": args.accepted_tsv,
        "rows_total": len(rows),
        "evaluated_universe": "Rows in selected-tsv with OK/SKIP_EXISTS ESMFold output; missing or failed predictions are NOT_EVALUATED, not structural FAIL.",
        "thresholds": {
            "core_ca_rmsd_max_A": args.rmsd_cutoff,
            "core_pair_max_delta_A": args.pair_cutoff,
            "core_mean_plddt_min": args.plddt_cutoff,
        },
        "counts_by_bin": counts_by_bin,
        "target_per_bin": args.target_per_bin,
        "target_complete": target_complete,
    }
    Path(args.summary_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary_json).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
