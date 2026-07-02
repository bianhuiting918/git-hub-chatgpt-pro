#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path("/data/bht/project01_baker_serhyd_routeB_20260701")
SCRIPT_DIR = ROOT / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import generate_baker_reference_ligandmpnn_bins as base  # noqa: E402

TEMPERATURES = {"90": 0.10, "80": 0.12, "70": 0.15, "60": 0.18, "50": 0.22}
DEFAULT_SELECT = {"90": 100, "80": 200, "70": 200, "60": 200, "50": 200}


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    return list(csv.DictReader(path.open(), delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def fixed_core_positions(row: dict[str, str]) -> set[tuple[str, int]]:
    start = int(row["parent_core_output_start"])
    end = int(row["parent_core_output_end"])
    return {("A", resid) for resid in range(start, end + 1)}


def select_mutable_positions(residues, fixed: set[tuple[str, int]], target_mutations: int) -> set[tuple[str, int]]:
    mutable = [(chain, resid) for chain, resid, aa in residues if (chain, resid) not in fixed and aa != "X"]
    target = min(target_mutations, len(mutable))
    return set(base.evenly_select(mutable, target))


def ligandmpnn_counts(target_records: int) -> tuple[int, int]:
    batch_size = min(250, max(10, target_records))
    batches = max(1, math.ceil(target_records / batch_size))
    return batch_size, batches


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=str(ROOT))
    parser.add_argument("--gate-accepted-tsv", required=True, help="Accepted TSV from gate_pocket4_l2_parent_core.py")
    parser.add_argument("--out-dir", default="sequence_design/pocket4_l2_ligandmpnn")
    parser.add_argument("--manifest-tsv", default="manifests/pocket4_l2_ligandmpnn_selected.tsv")
    parser.add_argument("--summary-json", default="manifests/pocket4_l2_ligandmpnn_summary.json")
    parser.add_argument("--parent-limit", type=int, default=5)
    parser.add_argument("--rounds", type=int, default=2)
    parser.add_argument("--seed-base", type=int, default=2703000)
    parser.add_argument("--bins", nargs="*", default=["90", "80", "70", "60", "50"], choices=["90", "80", "70", "60", "50"])
    parser.add_argument("--select-90", type=int, default=DEFAULT_SELECT["90"])
    parser.add_argument("--select-other", type=int, default=200)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    root = Path(args.root)
    gate_tsv = root / args.gate_accepted_tsv if not Path(args.gate_accepted_tsv).is_absolute() else Path(args.gate_accepted_tsv)
    out_dir = root / args.out_dir
    manifest_tsv = root / args.manifest_tsv
    summary_json = root / args.summary_json
    pred_dir = root / "structure_prediction/esmfold_pocket4_l2_ligandmpnn_bins"
    out_dir.mkdir(parents=True, exist_ok=True)
    pred_dir.mkdir(parents=True, exist_ok=True)

    parents = [row for row in read_tsv(gate_tsv) if row.get("gate") == "PASS" and row.get("pdb")]
    parents = parents[: args.parent_limit]
    if not parents:
        raise SystemExit(f"No PASS L2 parent PDBs found in {gate_tsv}")

    ligand_py = Path("/data/bht/design_tools/envs/ligandmpnn_venv/bin/python")
    ligand_run = Path("/data/bht/design_tools/src/LigandMPNN/run.py")
    selected_rows: list[dict[str, object]] = []
    summary = {
        "status": "RUNNING",
        "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "route": "Route B L2 LigandMPNN sequence regeneration from parent-core-gated L2 backbones",
        "gate_accepted_tsv": str(gate_tsv),
        "parent_limit": args.parent_limit,
        "rounds": args.rounds,
        "bins": args.bins,
        "parents": {},
    }

    filter_fields = [
        "bin", "record_index", "status", "identity", "mutation_count",
        "fixed_mutation_count", "outside_redesignable_mutation_count",
        "header", "fail_reasons", "sequence", "mutations",
    ]
    selected_fields = [
        "sample_id", "parent_sample_id", "parent_pdb", "fixed_core_start", "fixed_core_end",
        "round", "bin", "identity", "mutation_count", "sequence_length", "sequence",
        "postseq_predicted_pdb", "record_index", "fixed_mutation_count",
        "mutations", "header", "source_filter_tsv",
    ]

    for parent_index, parent in enumerate(parents, 1):
        parent_sample = parent.get("sample_id") or f"l2_parent{parent_index}"
        parent_pdb = Path(parent["pdb"])
        if not parent_pdb.is_absolute():
            parent_pdb = root / parent_pdb
        residues, ref_seq = base.parse_pdb_sequence(parent_pdb)
        fixed = fixed_core_positions(parent)
        base.FIXED = fixed
        parent_key = parent_sample.replace("/", "_")
        summary["parents"][parent_sample] = {
            "pdb": str(parent_pdb),
            "sequence_length": len(ref_seq),
            "fixed_core_start": parent["parent_core_output_start"],
            "fixed_core_end": parent["parent_core_output_end"],
            "fixed_residue_count": len(fixed),
            "bins": {},
        }
        for round_id in range(1, args.rounds + 1):
            for bin_id in args.bins:
                target_identity = int(bin_id) / 100.0
                requested_mutations = round(len(ref_seq) * (1.0 - target_identity))
                selected = select_mutable_positions(residues, fixed, requested_mutations)
                target_mutations = len(selected)
                base.BINS[bin_id]["target_mutations"] = target_mutations
                select_target = args.select_90 if bin_id == "90" else args.select_other
                batch_size, batches = ligandmpnn_counts(select_target)
                seed = args.seed_base + parent_index * 10000 + round_id * 1000 + int(bin_id)
                run_dir = out_dir / "runs" / parent_key / f"round{round_id:02d}_bin{bin_id}_seed{seed}"
                run_dir.mkdir(parents=True, exist_ok=True)
                filter_tsv = run_dir / "filter.tsv"
                log = run_dir / "ligandmpnn.log"
                fasta = run_dir / "seqs" / f"{parent_pdb.stem}.fa"
                elapsed_s = 0.0
                if filter_tsv.exists() and not args.force:
                    rows = read_tsv(filter_tsv)
                    proc_returncode: int | str = "SKIP_EXISTS"
                else:
                    omit = run_dir / "omit_native.json"
                    base.write_omit_native(omit, residues, selected)
                    redesigned = " ".join(f"{chain}{resid}" for chain, resid in sorted(selected, key=lambda item: item[1]))
                    cmd = [
                        str(ligand_py), str(ligand_run),
                        "--model_type", "ligand_mpnn",
                        "--seed", str(seed),
                        "--pdb_path", str(parent_pdb),
                        "--out_folder", str(run_dir),
                        "--chains_to_design", "A",
                        "--parse_these_chains_only", "A",
                        "--redesigned_residues", redesigned,
                        "--temperature", str(TEMPERATURES[bin_id]),
                        "--batch_size", str(batch_size),
                        "--number_of_batches", str(batches),
                        "--omit_AA_per_residue", str(omit),
                        "--save_stats", "0",
                        "--verbose", "0",
                    ]
                    start = time.time()
                    with log.open("w") as handle:
                        proc = subprocess.run(cmd, cwd=str(ligand_run.parent), stdout=handle, stderr=subprocess.STDOUT)
                    elapsed_s = round(time.time() - start, 2)
                    proc_returncode = proc.returncode
                    rows = base.filter_records(bin_id, base.parse_fasta(fasta), residues, ref_seq, selected)
                    write_tsv(filter_tsv, filter_fields, rows)
                passed = [row for row in rows if row["status"] == "pass"][:select_target]
                for rank, row in enumerate(passed, 1):
                    sample_id = f"pocket4_l2_{parent_key}_r{round_id:02d}_bin{bin_id}_rank{rank}_rec{row['record_index']}"
                    selected_rows.append({
                        "sample_id": sample_id,
                        "parent_sample_id": parent_sample,
                        "parent_pdb": str(parent_pdb),
                        "fixed_core_start": parent["parent_core_output_start"],
                        "fixed_core_end": parent["parent_core_output_end"],
                        "round": round_id,
                        "bin": bin_id,
                        "identity": row["identity"],
                        "mutation_count": row["mutation_count"],
                        "sequence_length": len(row["sequence"]),
                        "sequence": row["sequence"],
                        "postseq_predicted_pdb": str(pred_dir / f"{sample_id}.esmfold.pdb"),
                        "record_index": row["record_index"],
                        "fixed_mutation_count": row["fixed_mutation_count"],
                        "mutations": row["mutations"],
                        "header": row["header"],
                        "source_filter_tsv": str(filter_tsv),
                    })
                summary["parents"][parent_sample]["bins"][f"round{round_id:02d}_bin{bin_id}"] = {
                    "exit_code": proc_returncode,
                    "elapsed_s": elapsed_s,
                    "target_identity": bin_id,
                    "requested_mutations": requested_mutations,
                    "target_mutations": target_mutations,
                    "redesigned_residues": len(selected),
                    "records": len(rows),
                    "pass_filter_count": len([row for row in rows if row["status"] == "pass"]),
                    "selected_count": len(passed),
                    "selected_target": select_target,
                    "run_dir": str(run_dir),
                    "filter_tsv": str(filter_tsv),
                    "log": str(log),
                }
                summary_json.parent.mkdir(parents=True, exist_ok=True)
                summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    selected_rows.sort(key=lambda row: (str(row["parent_sample_id"]), int(row["round"]), int(row["bin"]), int(row["record_index"])))
    write_tsv(manifest_tsv, selected_fields, selected_rows)
    counts: dict[str, int] = {}
    for row in selected_rows:
        counts[str(row["bin"])] = counts.get(str(row["bin"]), 0) + 1
    summary.update({
        "status": "DONE",
        "finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "selected_tsv": str(manifest_tsv),
        "selected_total": len(selected_rows),
        "selected_counts_by_bin": counts,
    })
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
