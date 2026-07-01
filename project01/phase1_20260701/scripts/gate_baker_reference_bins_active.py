#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import itertools
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

FIXED = [13, 14, 15, 16, 17, 54, 55, 56, 72, 73, 74, 148, 149, 150]


def read_tsv(path: Path) -> list[dict[str, str]]:
    return list(csv.DictReader(path.open(), delimiter="\t")) if path.exists() else []


def ca_atoms(path: Path) -> tuple[dict[int, np.ndarray], dict[int, str]]:
    coords = {}
    names = {}
    if not path.exists():
        return coords, names
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and (line[21].strip() or "A") == "A":
            resid = int(line[22:26])
            coords[resid] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
            names[resid] = line[17:20].strip()
    return coords, names


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


def pair_delta(ref: dict[int, np.ndarray], cand: dict[int, np.ndarray]) -> tuple[float, float]:
    values = []
    for i, j in itertools.combinations(FIXED, 2):
        values.append(abs(np.linalg.norm(ref[i] - ref[j]) - np.linalg.norm(cand[i] - cand[j])))
    return max(values), float(np.mean(values))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ref-pdb", required=True)
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--status-tsv", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--rmsd-cutoff", type=float, default=1.0)
    parser.add_argument("--pair-cutoff", type=float, default=1.0)
    args = parser.parse_args()

    ref_coords, ref_names = ca_atoms(Path(args.ref_pdb))
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
        if status.get("status") not in {"OK", "SKIP_EXISTS"} or not pdb.exists():
            gate = "NOT_EVALUATED"
            reasons.append("missing_or_pending_esmfold")
        else:
            cand_coords, cand_names = ca_atoms(pdb)
            missing = [resid for resid in FIXED if resid not in cand_coords or resid not in ref_coords]
            fixed_mut = [
                f"A{resid}:{ref_names.get(resid)}>{cand_names.get(resid)}"
                for resid in FIXED
                if cand_names.get(resid) != ref_names.get(resid)
            ]
            if missing:
                reasons.append("missing_fixed_ca")
            if fixed_mut:
                reasons.append("fixed_identity_mutation")
            if not missing:
                rmsd = kabsch_rmsd([ref_coords[resid] for resid in FIXED], [cand_coords[resid] for resid in FIXED])
                pmax, pmean = pair_delta(ref_coords, cand_coords)
                if rmsd > args.rmsd_cutoff:
                    reasons.append("motif_ca_rmsd_gt_cutoff")
                if pmax > args.pair_cutoff:
                    reasons.append("motif_pair_delta_gt_cutoff")
            gate = "FINAL_QUALIFIED_ACTIVE" if not reasons else "FAIL"
        out_rows.append({
            "sample_id": sid,
            "bin": row.get("bin", ""),
            "identity": row.get("identity", ""),
            "mutation_count": row.get("mutation_count", ""),
            "esmfold_status": status.get("status", "NOT_EVALUATED"),
            "gate": gate,
            "fail_reasons": ";".join(reasons),
            "motif_ca_rmsd_A": "" if rmsd == "" else f"{rmsd:.4f}",
            "motif_pair_max_delta_A": "" if pmax == "" else f"{pmax:.4f}",
            "motif_pair_mean_delta_A": "" if pmean == "" else f"{pmean:.4f}",
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
        rows.sort(key=lambda item: float(item["motif_ca_rmsd_A"]))

    summary = {
        "status": "DONE",
        "gate_definition": (
            "FINAL_QUALIFIED_ACTIVE if ESMFold OK, fixed identities retained, "
            f"motif CA RMSD <= {args.rmsd_cutoff} A and motif pair max delta <= {args.pair_cutoff} A; "
            "sidechain/ligand compatibility remains a downstream stricter gate"
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
