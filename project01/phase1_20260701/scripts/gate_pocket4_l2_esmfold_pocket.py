#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import itertools
import json
import math
from pathlib import Path

import numpy as np


DEFAULT_POCKET = "38,39,40,43,61,62,65,76,79,80,87,91,102,103,104,105"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected-tsv", required=True)
    parser.add_argument("--esmfold-status-tsv", required=True)
    parser.add_argument("--parent-pdb", default="outputs/ca_rfd_baker_pocket4_l2_extend_expanded_n1_20260703/sample_7500.pdb")
    parser.add_argument("--pocket-resids", default=DEFAULT_POCKET)
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
            coords[(chain, resid)] = (xyz, float(line[60:66]))
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


def evaluate(
    row: dict[str, str],
    status_by_sample: dict[str, dict[str, str]],
    parent_ca: dict[tuple[str, int], tuple[np.ndarray, float]],
    pocket_resids: list[int],
    args: argparse.Namespace,
) -> dict[str, object]:
    sample_id = row["sample_id"]
    out = {
        "sample_id": sample_id,
        "bin": row.get("bin", ""),
        "identity": row.get("identity", ""),
        "gate": "NOT_EVALUATED",
        "fail_reasons": "missing_or_failed_esmfold",
        "pocket_ca_count": "",
        "pocket_ca_rmsd_A": "",
        "pocket_pair_max_delta_A": "",
        "pocket_pair_mean_delta_A": "",
        "pocket_mean_plddt": "",
        "predicted_pdb": "",
        "parent_pdb": args.parent_pdb,
        "pocket_resids": ",".join(map(str, pocket_resids)),
    }
    status = status_by_sample.get(sample_id)
    if not status or status.get("status") not in {"OK", "SKIP_EXISTS"}:
        return out
    pred_pdb = Path(status["pdb"])
    if not pred_pdb.exists():
        return out
    pred_ca = ca_atoms(pred_pdb)
    ref_points: list[np.ndarray] = []
    cand_points: list[np.ndarray] = []
    plddts: list[float] = []
    missing: list[str] = []
    for resid in pocket_resids:
        key = ("A", resid)
        if key not in parent_ca or key not in pred_ca:
            missing.append(f"A{resid}")
            continue
        ref_points.append(parent_ca[key][0])
        cand_points.append(pred_ca[key][0])
        plddts.append(pred_ca[key][1])
    reasons: list[str] = []
    if missing:
        reasons.append("missing_pocket_ca")
    rmsd = pair_max = pair_mean = mean_plddt = None
    if not missing and ref_points:
        rmsd = kabsch_rmsd(ref_points, cand_points)
        pair_max, pair_mean = pair_delta(ref_points, cand_points)
        mean_plddt = float(np.mean(plddts))
        if rmsd > args.rmsd_cutoff:
            reasons.append("pocket_ca_rmsd_gt_cutoff")
        if pair_max > args.pair_cutoff:
            reasons.append("pocket_pair_delta_gt_cutoff")
        if mean_plddt < args.plddt_cutoff:
            reasons.append("pocket_plddt_lt_cutoff")
    out.update({
        "gate": "PASS" if not reasons else "FAIL",
        "fail_reasons": ";".join(reasons),
        "pocket_ca_count": len(pocket_resids),
        "pocket_ca_rmsd_A": "" if rmsd is None else f"{rmsd:.4f}",
        "pocket_pair_max_delta_A": "" if pair_max is None else f"{pair_max:.4f}",
        "pocket_pair_mean_delta_A": "" if pair_mean is None else f"{pair_mean:.4f}",
        "pocket_mean_plddt": "" if mean_plddt is None else f"{mean_plddt:.2f}",
        "predicted_pdb": str(pred_pdb),
    })
    return out


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
    pocket_resids = [int(item.strip()) for item in args.pocket_resids.split(",") if item.strip()]
    parent_ca = ca_atoms(Path(args.parent_pdb))
    rows = [evaluate(row, status_by_sample, parent_ca, pocket_resids, args) for row in selected_rows]
    fields = [
        "sample_id", "bin", "identity", "gate", "fail_reasons", "pocket_ca_count",
        "pocket_ca_rmsd_A", "pocket_pair_max_delta_A", "pocket_pair_mean_delta_A",
        "pocket_mean_plddt", "predicted_pdb", "parent_pdb", "pocket_resids",
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
        "route": "Route B2 L2 ESMFold pocket gate",
        "selected_tsv": args.selected_tsv,
        "esmfold_status_tsv": args.esmfold_status_tsv,
        "rows_total": len(rows),
        "evaluated_universe": "Rows in selected-tsv with OK/SKIP_EXISTS ESMFold output; missing or failed predictions are NOT_EVALUATED, not structural FAIL.",
        "thresholds": {
            "pocket_ca_rmsd_max_A": args.rmsd_cutoff,
            "pocket_pair_max_delta_A": args.pair_cutoff,
            "pocket_mean_plddt_min": args.plddt_cutoff,
        },
        "pocket_resids": pocket_resids,
        "counts_by_bin": counts_by_bin,
        "target_per_bin": args.target_per_bin,
        "target_complete": target_complete,
        "out_tsv": args.out_tsv,
        "accepted_tsv": args.accepted_tsv,
    }
    Path(args.summary_json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.summary_json).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
