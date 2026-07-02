#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


THRESHOLDS = {
    "motif_ca_rmsd_max_A": 1.0,
    "motif_pair_max_delta_max_A": 1.0,
    "motif_mean_plddt_min": 70.0,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--selected-tsv", required=True)
    parser.add_argument("--esmfold-status-tsv", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--accepted-tsv", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--target-per-bin", type=int, default=10)
    parser.add_argument("--motif-rmsd-cutoff", type=float, default=1.0)
    parser.add_argument("--pair-delta-cutoff", type=float, default=1.0)
    parser.add_argument("--motif-plddt-min", type=float, default=70.0)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def motif_positions_from_contig(contig: str) -> list[int]:
    out_pos = 1
    fixed: list[int] = []
    for token in [item.strip() for item in contig.split(",") if item.strip()]:
        if token.isdigit():
            out_pos += int(token)
            continue
        chain = token[0]
        if chain not in {"A", "B"}:
            raise ValueError(f"Unsupported contig token: {token}")
        bounds = token[1:]
        if "-" not in bounds:
            raise ValueError(f"Unsupported contig token: {token}")
        start_s, end_s = bounds.split("-", 1)
        start, end = int(start_s), int(end_s)
        step = 1 if end >= start else -1
        for _ in range(start, end + step, step):
            fixed.append(out_pos)
            out_pos += 1
    return fixed


def parse_pdb(path: Path) -> list[dict[str, object]]:
    atoms: list[dict[str, object]] = []
    for line in path.read_text(errors="ignore").splitlines():
        rec = line[:6].strip()
        if rec not in {"ATOM", "HETATM"} or len(line) < 54:
            continue
        try:
            atoms.append(
                {
                    "rec": rec,
                    "atom": line[12:16].strip(),
                    "resname": line[17:20].strip(),
                    "chain": line[21].strip() or "A",
                    "resseq": int(line[22:26]),
                    "coord": np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])], dtype=float),
                    "bfactor": float(line[60:66]) if len(line) >= 66 else math.nan,
                }
            )
        except ValueError:
            continue
    return atoms


def ca_map(atoms: list[dict[str, object]]) -> dict[int, dict[str, object]]:
    return {
        int(atom["resseq"]): atom
        for atom in atoms
        if atom["rec"] == "ATOM" and atom["chain"] == "A" and atom["atom"] == "CA"
    }


def kabsch_rmsd(ref_xyz: np.ndarray, cand_xyz: np.ndarray) -> float:
    ref_center = ref_xyz.mean(axis=0)
    cand_center = cand_xyz.mean(axis=0)
    ref0 = ref_xyz - ref_center
    cand0 = cand_xyz - cand_center
    u, _s, vt = np.linalg.svd(cand0.T @ ref0)
    rot = u @ vt
    if np.linalg.det(rot) < 0:
        u[:, -1] *= -1
        rot = u @ vt
    aligned = cand0 @ rot + ref_center
    return float(np.sqrt(((aligned - ref_xyz) ** 2).sum() / len(ref_xyz)))


def pair_deltas(ref_xyz: np.ndarray, cand_xyz: np.ndarray) -> tuple[float, float]:
    values: list[float] = []
    for i in range(len(ref_xyz)):
        for j in range(i + 1, len(ref_xyz)):
            values.append(abs(float(np.linalg.norm(ref_xyz[i] - ref_xyz[j]) - np.linalg.norm(cand_xyz[i] - cand_xyz[j]))))
    return max(values), float(np.mean(values))


def evaluate(row: dict[str, str], status: dict[str, str], parent_cache: dict[str, tuple[list[int], np.ndarray]]) -> dict[str, object]:
    sample_id = row["sample_id"]
    if not status or status.get("status") not in {"OK", "SKIP_EXISTS"}:
        return {
            "sample_id": sample_id,
            "parent_sample_id": row.get("parent_sample_id", ""),
            "bin": row.get("bin", ""),
            "identity": row.get("identity", ""),
            "mutation_count": row.get("mutation_count", ""),
            "gate": "NOT_EVALUATED",
            "fail_reasons": "missing_or_failed_esmfold",
            "motif_ca_rmsd_A": "",
            "motif_pair_max_delta_A": "",
            "motif_pair_mean_delta_A": "",
            "motif_mean_plddt": "",
            "predicted_pdb": status.get("pdb", ""),
            "parent_pdb": row.get("parent_pdb", ""),
        }
    parent_pdb = row["parent_pdb"]
    if parent_pdb not in parent_cache:
        positions = motif_positions_from_contig(row["contig"])
        parent_ca = ca_map(parse_pdb(Path(parent_pdb)))
        missing_parent = [pos for pos in positions if pos not in parent_ca]
        if missing_parent:
            raise ValueError(f"Missing parent motif CA positions in {parent_pdb}: {missing_parent}")
        parent_cache[parent_pdb] = (positions, np.array([parent_ca[pos]["coord"] for pos in positions], dtype=float))
    positions, parent_xyz = parent_cache[parent_pdb]
    pred = Path(status["pdb"])
    cand_ca = ca_map(parse_pdb(pred))
    missing = [pos for pos in positions if pos not in cand_ca]
    reasons: list[str] = []
    if missing:
        reasons.append("missing_motif_ca")
        rmsd = None
        pair_max = None
        pair_mean = None
        mean_plddt = None
    else:
        cand_xyz = np.array([cand_ca[pos]["coord"] for pos in positions], dtype=float)
        rmsd = kabsch_rmsd(parent_xyz, cand_xyz)
        pair_max, pair_mean = pair_deltas(parent_xyz, cand_xyz)
        plddts = [float(cand_ca[pos]["bfactor"]) for pos in positions if not math.isnan(float(cand_ca[pos]["bfactor"]))]
        mean_plddt = float(np.mean(plddts)) if plddts else None
        if rmsd > THRESHOLDS["motif_ca_rmsd_max_A"]:
            reasons.append("motif_ca_rmsd_gt_cutoff")
        if pair_max > THRESHOLDS["motif_pair_max_delta_max_A"]:
            reasons.append("motif_pair_delta_gt_cutoff")
        if mean_plddt is None or mean_plddt < THRESHOLDS["motif_mean_plddt_min"]:
            reasons.append("motif_plddt_lt_cutoff")
    return {
        "sample_id": sample_id,
        "parent_sample_id": row.get("parent_sample_id", ""),
        "bin": row.get("bin", ""),
        "identity": row.get("identity", ""),
        "mutation_count": row.get("mutation_count", ""),
        "gate": "PASS" if not reasons else "FAIL",
        "fail_reasons": ";".join(reasons),
        "motif_ca_rmsd_A": "" if rmsd is None else f"{rmsd:.4f}",
        "motif_pair_max_delta_A": "" if pair_max is None else f"{pair_max:.4f}",
        "motif_pair_mean_delta_A": "" if pair_mean is None else f"{pair_mean:.4f}",
        "motif_mean_plddt": "" if mean_plddt is None else f"{mean_plddt:.2f}",
        "predicted_pdb": str(pred),
        "parent_pdb": parent_pdb,
    }


def main() -> int:
    args = parse_args()
    THRESHOLDS["motif_ca_rmsd_max_A"] = args.motif_rmsd_cutoff
    THRESHOLDS["motif_pair_max_delta_max_A"] = args.pair_delta_cutoff
    THRESHOLDS["motif_mean_plddt_min"] = args.motif_plddt_min

    selected_rows = read_tsv(Path(args.selected_tsv))
    status_by_id = {row["sample_id"]: row for row in read_tsv(Path(args.esmfold_status_tsv))}
    parent_cache: dict[str, tuple[list[int], np.ndarray]] = {}
    rows = [evaluate(row, status_by_id.get(row["sample_id"], {}), parent_cache) for row in selected_rows]
    fields = list(rows[0].keys()) if rows else ["sample_id"]
    write_tsv(Path(args.out_tsv), fields, rows)

    by_bin: dict[str, dict[str, int]] = defaultdict(lambda: {"evaluated": 0, "not_evaluated": 0, "pass": 0, "fail": 0})
    for row in rows:
        bucket = by_bin[str(row.get("bin", ""))]
        if row["gate"] == "NOT_EVALUATED":
            bucket["not_evaluated"] += 1
        else:
            bucket["evaluated"] += 1
        if row["gate"] == "PASS":
            bucket["pass"] += 1
        elif row["gate"] == "FAIL":
            bucket["fail"] += 1

    pass_rows = [row for row in rows if row["gate"] == "PASS"]
    pass_rows.sort(key=lambda row: (int(row["bin"]), float(row["motif_ca_rmsd_A"] or 999.0), float(row["motif_pair_max_delta_A"] or 999.0)))
    accepted: list[dict[str, object]] = []
    selected_by_bin: dict[str, int] = defaultdict(int)
    for row in pass_rows:
        if selected_by_bin[str(row["bin"])] < args.target_per_bin:
            accepted.append(row)
            selected_by_bin[str(row["bin"])] += 1
    write_tsv(Path(args.accepted_tsv), fields, accepted)
    for bin_id, count in selected_by_bin.items():
        by_bin[bin_id]["accepted"] = count

    summary = {
        "status": "DONE",
        "selected_tsv": args.selected_tsv,
        "esmfold_status_tsv": args.esmfold_status_tsv,
        "out_tsv": args.out_tsv,
        "accepted_tsv": args.accepted_tsv,
        "rows_total": len(rows),
        "counts_by_bin": {key: dict(value) for key, value in sorted(by_bin.items())},
        "target_per_bin": args.target_per_bin,
        "target_complete": all(by_bin[str(bin_id)].get("accepted", 0) >= args.target_per_bin for bin_id in [90, 80, 70, 60, 50]),
        "thresholds": THRESHOLDS,
        "evaluated_universe": "Rows in selected-tsv with OK/SKIP_EXISTS ESMFold output; missing or failed predictions are NOT_EVALUATED, not structural FAIL.",
        "gate_definition": "Predicted apo protein motif CA geometry is compared to the parent CA_RFDiffusion scaffold motif positions from the layered contig.",
    }
    Path(args.summary_json).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
