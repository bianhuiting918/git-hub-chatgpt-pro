#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path("/data/bht/project01_baker_serhyd_routeB_20260701")
SCRIPT_DIR = ROOT / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import generate_baker_reference_ligandmpnn_bins as base  # noqa: E402
import generate_baker_theozyme_sample1000_ligandmpnn_bins as primary  # noqa: E402


def read_selected(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    return list(csv.DictReader(path.open(), delimiter="\t"))


def write_selected(path: Path, rows: list[dict[str, str]]) -> None:
    fields = [
        "sample_id", "bin", "identity", "mutation_count", "sequence_length",
        "sequence", "postseq_predicted_pdb", "record_index", "fixed_mutation_count",
        "mutations", "header", "source_backbone_pdb",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run_seed(bin_id: str, seed: int, batch_size: int, batches: int, temperature: float, force: bool):
    residues, ref_seq = base.parse_pdb_sequence(primary.PARENT_PDB)
    selected = primary.evenly_select_nonfixed(residues, primary.BINS[bin_id]["target_mutations"])
    run_dir = primary.OUT / "runs" / f"bin{bin_id}_sample1000_refined_seed{seed}"
    run_dir.mkdir(parents=True, exist_ok=True)
    filter_tsv = run_dir / f"filter_bin{bin_id}.tsv"
    fasta = run_dir / "seqs" / f"{primary.PARENT_PDB.stem}.fa"
    log = run_dir / "ligandmpnn.log"
    if filter_tsv.exists() and not force:
        return list(csv.DictReader(filter_tsv.open(), delimiter="\t")), "SKIP_EXISTS", 0.0, str(run_dir), str(log), str(fasta), str(filter_tsv)

    omit = run_dir / f"bin{bin_id}_omit_native.json"
    base.write_omit_native(omit, residues, selected)
    redesigned = " ".join(f"{chain}{resid}" for chain, resid in sorted(selected, key=lambda item: item[1]))
    ligand_py = Path("/data/bht/design_tools/envs/ligandmpnn_venv/bin/python")
    ligand_run = Path("/data/bht/design_tools/src/LigandMPNN/run.py")
    cmd = [
        str(ligand_py), str(ligand_run),
        "--model_type", "ligand_mpnn",
        "--seed", str(seed),
        "--pdb_path", str(primary.PARENT_PDB),
        "--out_folder", str(run_dir),
        "--chains_to_design", "A",
        "--parse_these_chains_only", "A",
        "--redesigned_residues", redesigned,
        "--temperature", str(temperature),
        "--batch_size", str(batch_size),
        "--number_of_batches", str(batches),
        "--omit_AA_per_residue", str(omit),
        "--save_stats", "0",
        "--verbose", "0",
    ]
    start = time.time()
    with log.open("w") as handle:
        proc = subprocess.run(cmd, cwd=str(ligand_run.parent), stdout=handle, stderr=subprocess.STDOUT)
    elapsed = round(time.time() - start, 2)
    rows = base.filter_records(bin_id, base.parse_fasta(fasta), residues, ref_seq, selected)
    fields = [
        "bin", "record_index", "status", "identity", "mutation_count",
        "fixed_mutation_count", "outside_redesignable_mutation_count",
        "header", "fail_reasons", "sequence", "mutations",
    ]
    with filter_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows, proc.returncode, elapsed, str(run_dir), str(log), str(fasta), str(filter_tsv)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bin", required=True, choices=list(primary.BINS))
    parser.add_argument("--target", type=int, default=None)
    parser.add_argument("--seeds", nargs="+", type=int, required=True)
    parser.add_argument("--batch-size", type=int, default=None)
    parser.add_argument("--batches", type=int, default=None)
    parser.add_argument("--temperature", type=float, default=None)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    bin_id = args.bin
    spec = primary.BINS[bin_id]
    target = args.target or spec["select"]
    batch_size = args.batch_size or spec["batch_size"]
    batches = args.batches or spec["batches"]
    temperature = args.temperature if args.temperature is not None else spec["temperature"]

    primary.MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    primary.OUT.mkdir(parents=True, exist_ok=True)
    primary.PRED_DIR.mkdir(parents=True, exist_ok=True)
    base.FIXED = primary.FIXED
    base.BINS = primary.BINS

    selected = read_selected(primary.MANIFEST)
    kept_other_bins = [row for row in selected if row.get("bin") != bin_id]
    existing_bin_rows = [row for row in selected if row.get("bin") == bin_id]
    by_sequence = {row["sequence"]: row for row in existing_bin_rows}
    seed_summaries = []

    for seed in args.seeds:
        if len(by_sequence) >= target:
            break
        rows, exit_code, elapsed, run_dir, log, fasta, filter_tsv = run_seed(bin_id, seed, batch_size, batches, temperature, args.force)
        pass_rows = [row for row in rows if row["status"] == "pass"]
        added = 0
        for row in pass_rows:
            if len(by_sequence) >= target:
                break
            if row["sequence"] in by_sequence:
                continue
            rank = len(by_sequence) + 1
            sample_id = f"baker_theozyme_sample1000_refined_bin{bin_id}_rank{rank}_seed{seed}_rec{row['record_index']}"
            by_sequence[row["sequence"]] = {
                "sample_id": sample_id,
                "bin": row["bin"],
                "identity": row["identity"],
                "mutation_count": row["mutation_count"],
                "sequence_length": len(row["sequence"]),
                "sequence": row["sequence"],
                "postseq_predicted_pdb": str(primary.PRED_DIR / f"{sample_id}.esmfold.pdb"),
                "record_index": row["record_index"],
                "fixed_mutation_count": row["fixed_mutation_count"],
                "mutations": row["mutations"],
                "header": row["header"],
                "source_backbone_pdb": str(primary.PARENT_PDB),
            }
            added += 1
        seed_summaries.append({
            "seed": seed,
            "exit_code": exit_code,
            "elapsed_s": elapsed,
            "records": len(rows),
            "pass_filter_count": len(pass_rows),
            "added_unique": added,
            "merged_unique_count": len(by_sequence),
            "run_dir": run_dir,
            "log": log,
            "fasta": fasta,
            "filter_tsv": filter_tsv,
        })

    merged_bin_rows = list(by_sequence.values())[:target]
    all_rows = kept_other_bins + merged_bin_rows
    all_rows.sort(key=lambda row: (int(row["bin"]), row["sample_id"]))
    write_selected(primary.MANIFEST, all_rows)

    summary = json.loads(primary.SUMMARY.read_text()) if primary.SUMMARY.exists() else {}
    counts = {}
    for row in all_rows:
        counts[row["bin"]] = counts.get(row["bin"], 0) + 1
    summary.setdefault("bins", {})
    summary["bins"][bin_id] = {
        **summary.get("bins", {}).get(bin_id, {}),
        "supplement_status": "TARGET_MET" if len(merged_bin_rows) >= target else "TARGET_NOT_MET",
        "selected_count": len(merged_bin_rows),
        "selected_target": target,
        "supplement_seeds": seed_summaries,
        "supplement_finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
    }
    summary.update({
        "status": "DONE",
        "selected_total": len(all_rows),
        "selected_counts_by_bin": counts,
        "selected_tsv": str(primary.MANIFEST),
        "finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
    })
    primary.SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
