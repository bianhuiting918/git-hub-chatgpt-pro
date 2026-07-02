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

PARENT_PDB = ROOT / "outputs/ca_rfd_baker_theozyme_formal_constraints_batch50_20260702/sample_1000_refined_0.pdb"
OUT = ROOT / "sequence_design/baker_theozyme_sample1000_refined_ligandmpnn_bins"
MANIFEST = ROOT / "manifests/baker_theozyme_sample1000_refined_ligandmpnn_bins_selected.tsv"
SUMMARY = ROOT / "manifests/baker_theozyme_sample1000_refined_ligandmpnn_bins_summary.json"
PRED_DIR = ROOT / "structure_prediction/esmfold_baker_theozyme_sample1000_refined_bins"

FIXED = {("A", r) for r in (13, 14, 15, 16, 17, 54, 55, 56, 72, 73, 74, 148, 149, 150)}
BINS = {
    "90": {"target_mutations": 16, "seed": 991001, "temperature": 0.10, "batch_size": 100, "batches": 4, "select": 200},
    "80": {"target_mutations": 32, "seed": 891001, "temperature": 0.12, "batch_size": 250, "batches": 5, "select": 1000},
    "70": {"target_mutations": 48, "seed": 791001, "temperature": 0.15, "batch_size": 250, "batches": 5, "select": 1000},
    "60": {"target_mutations": 64, "seed": 691001, "temperature": 0.18, "batch_size": 250, "batches": 5, "select": 1000},
    "50": {"target_mutations": 80, "seed": 591001, "temperature": 0.22, "batch_size": 250, "batches": 5, "select": 1000},
}


def load_existing_selected(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    return list(csv.DictReader(path.open(), delimiter="\t"))


def evenly_select_nonfixed(residues, count: int) -> set[tuple[str, int]]:
    mutable = [(chain, resid) for chain, resid, aa in residues if (chain, resid) not in FIXED and aa != "X"]
    return set(base.evenly_select(mutable, count))


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bins", nargs="*", default=list(BINS), choices=list(BINS))
    parser.add_argument("--force", action="store_true", help="rerun bins even when their filter TSV already exists")
    args = parser.parse_args()

    if not PARENT_PDB.exists():
        raise SystemExit(f"missing parent PDB: {PARENT_PDB}")

    OUT.mkdir(parents=True, exist_ok=True)
    PRED_DIR.mkdir(parents=True, exist_ok=True)
    runs = OUT / "runs"
    runs.mkdir(parents=True, exist_ok=True)
    MANIFEST.parent.mkdir(parents=True, exist_ok=True)

    base.FIXED = FIXED
    base.BINS = BINS
    residues, ref_seq = base.parse_pdb_sequence(PARENT_PDB)
    ligand_py = Path("/data/bht/design_tools/envs/ligandmpnn_venv/bin/python")
    ligand_run = Path("/data/bht/design_tools/src/LigandMPNN/run.py")

    selected_all = [row for row in load_existing_selected(MANIFEST) if row.get("bin") not in set(args.bins) or not args.force]
    summary = json.loads(SUMMARY.read_text()) if SUMMARY.exists() else {
        "status": "RUNNING",
        "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "reference_parent_pdb": str(PARENT_PDB),
        "route": "baker_theozyme_new_backbone_sample1000_refined",
        "sequence_length": len(ref_seq),
        "fixed_policy": "Baker theozyme contig motif residues only; no ligand-4A shell lock",
        "fixed_residue_count": len(FIXED),
        "fixed_residues": [f"{chain}{resid}" for chain, resid in sorted(FIXED, key=lambda item: item[1])],
        "selected_per_bin_target": {bin_id: spec["select"] for bin_id, spec in BINS.items()},
        "bins": {},
    }
    summary["status"] = "RUNNING"
    summary["last_started"] = time.strftime("%Y-%m-%dT%H:%M:%S")

    filter_fields = [
        "bin", "record_index", "status", "identity", "mutation_count",
        "fixed_mutation_count", "outside_redesignable_mutation_count",
        "header", "fail_reasons", "sequence", "mutations",
    ]
    selected_fields = [
        "sample_id", "bin", "identity", "mutation_count", "sequence_length",
        "sequence", "postseq_predicted_pdb", "record_index", "fixed_mutation_count",
        "mutations", "header", "source_backbone_pdb",
    ]

    for bin_id in args.bins:
        spec = BINS[bin_id]
        selected = evenly_select_nonfixed(residues, spec["target_mutations"])
        run_dir = runs / f"bin{bin_id}_sample1000_refined_seed{spec['seed']}"
        run_dir.mkdir(parents=True, exist_ok=True)
        filter_tsv = run_dir / f"filter_bin{bin_id}.tsv"
        if filter_tsv.exists() and not args.force:
            rows = list(csv.DictReader(filter_tsv.open(), delimiter="\t"))
            proc_returncode = "SKIP_EXISTS"
            elapsed_s = 0.0
            log = run_dir / "ligandmpnn.log"
            fasta = run_dir / "seqs" / f"{PARENT_PDB.stem}.fa"
        else:
            omit = run_dir / f"bin{bin_id}_omit_native.json"
            base.write_omit_native(omit, residues, selected)
            redesigned = " ".join(f"{chain}{resid}" for chain, resid in sorted(selected, key=lambda item: item[1]))
            cmd = [
                str(ligand_py), str(ligand_run),
                "--model_type", "ligand_mpnn",
                "--seed", str(spec["seed"]),
                "--pdb_path", str(PARENT_PDB),
                "--out_folder", str(run_dir),
                "--chains_to_design", "A",
                "--parse_these_chains_only", "A",
                "--redesigned_residues", redesigned,
                "--temperature", str(spec["temperature"]),
                "--batch_size", str(spec["batch_size"]),
                "--number_of_batches", str(spec["batches"]),
                "--omit_AA_per_residue", str(omit),
                "--save_stats", "0",
                "--verbose", "0",
            ]
            log = run_dir / "ligandmpnn.log"
            start = time.time()
            with log.open("w") as handle:
                proc = subprocess.run(cmd, cwd=str(ligand_run.parent), stdout=handle, stderr=subprocess.STDOUT)
            elapsed_s = round(time.time() - start, 2)
            proc_returncode = proc.returncode
            fasta = run_dir / "seqs" / f"{PARENT_PDB.stem}.fa"
            rows = base.filter_records(bin_id, base.parse_fasta(fasta), residues, ref_seq, selected)
            with filter_tsv.open("w", newline="") as handle:
                writer = csv.DictWriter(handle, filter_fields, delimiter="\t", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)

        passed = [row for row in rows if row["status"] == "pass"][: spec["select"]]
        selected_all = [row for row in selected_all if row.get("bin") != bin_id]
        for rank, row in enumerate(passed, 1):
            sample_id = f"baker_theozyme_sample1000_refined_bin{bin_id}_rank{rank}_rec{row['record_index']}"
            selected_all.append({
                "sample_id": sample_id,
                "bin": row["bin"],
                "identity": row["identity"],
                "mutation_count": row["mutation_count"],
                "sequence_length": len(row["sequence"]),
                "sequence": row["sequence"],
                "postseq_predicted_pdb": str(PRED_DIR / f"{sample_id}.esmfold.pdb"),
                "record_index": row["record_index"],
                "fixed_mutation_count": row["fixed_mutation_count"],
                "mutations": row["mutations"],
                "header": row["header"],
                "source_backbone_pdb": str(PARENT_PDB),
            })
        summary["bins"][bin_id] = {
            "exit_code": proc_returncode,
            "elapsed_s": elapsed_s,
            "target_mutations": spec["target_mutations"],
            "target_identity": bin_id,
            "redesigned_residues": len(selected),
            "records": len(rows),
            "pass_filter_count": len([row for row in rows if row["status"] == "pass"]),
            "selected_count": len(passed),
            "selected_target": spec["select"],
            "run_dir": str(run_dir),
            "fasta": str(fasta),
            "filter_tsv": str(filter_tsv),
            "log": str(log),
        }
        SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    selected_all.sort(key=lambda row: (int(row["bin"]), int(row["record_index"])))
    with MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, selected_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected_all)
    counts_by_bin = {}
    for row in selected_all:
        counts_by_bin[row["bin"]] = counts_by_bin.get(row["bin"], 0) + 1
    summary.update({
        "status": "DONE",
        "selected_total": len(selected_all),
        "selected_counts_by_bin": counts_by_bin,
        "selected_tsv": str(MANIFEST),
        "finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
    })
    SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
