#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("/data/bht/project01_baker_serhyd_routeB_20260701")
SCRIPT_DIR = ROOT / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import generate_baker_reference_ligandmpnn_bins as base  # noqa: E402

PDB = ROOT / "external/serine-hydrolase-design/design_pipeline/04_af2/outputs/super_design_model_4_ptm_seed_0_unrelaxed.pdb"
OUT = ROOT / "sequence_design/baker_reference_shell8_ligandmpnn_bins"
MANIFEST = ROOT / "manifests/baker_reference_shell8_ligandmpnn_bins_selected.tsv"
SUMMARY = ROOT / "manifests/baker_reference_shell8_ligandmpnn_bins_summary.json"
PRED_DIR = ROOT / "structure_prediction/esmfold_baker_reference_shell8_bins"
MOTIF = {13, 14, 15, 16, 17, 54, 55, 56, 72, 73, 74, 148, 149, 150}
SHELL_CUTOFF = 8.0
SELECT_PER_BIN = 20
BINS = {
    "90": {"target_mutations": 16, "seed": 94201, "temperature": 0.10, "batch_size": 20, "batches": 3},
    "80": {"target_mutations": 32, "seed": 84201, "temperature": 0.12, "batch_size": 20, "batches": 3},
    "70": {"target_mutations": 48, "seed": 74201, "temperature": 0.15, "batch_size": 20, "batches": 3},
    "60": {"target_mutations": 64, "seed": 64201, "temperature": 0.18, "batch_size": 20, "batches": 3},
    "50": {"target_mutations": 80, "seed": 54201, "temperature": 0.22, "batch_size": 20, "batches": 3},
}


def ca_shell_fixed(path: Path) -> set[tuple[str, int]]:
    coords = {}
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and (line[21].strip() or "A") == "A":
            coords[int(line[22:26])] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    fixed = set(MOTIF)
    for resid, coord in coords.items():
        if min(np.linalg.norm(coord - coords[m]) for m in MOTIF if m in coords) <= SHELL_CUTOFF:
            fixed.add(resid)
    return {("A", resid) for resid in fixed}


def main() -> int:
    fixed = ca_shell_fixed(PDB)
    base.FIXED = fixed
    base.BINS = BINS
    OUT.mkdir(parents=True, exist_ok=True)
    PRED_DIR.mkdir(parents=True, exist_ok=True)
    runs = OUT / "runs"
    runs.mkdir(parents=True, exist_ok=True)
    residues, ref_seq = base.parse_pdb_sequence(PDB)
    mutable = [(chain, resid) for chain, resid, aa in residues if (chain, resid) not in fixed and aa != "X"]
    selected_all = []
    ligand_py = Path("/data/bht/design_tools/envs/ligandmpnn_venv/bin/python")
    ligand_run = Path("/data/bht/design_tools/src/LigandMPNN/run.py")
    summary = {
        "status": "RUNNING",
        "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "reference_pdb": str(PDB),
        "sequence_length": len(ref_seq),
        "shell_cutoff_ca_A": SHELL_CUTOFF,
        "fixed_residue_count": len(fixed),
        "mutable_positions": len(mutable),
        "selected_per_bin": SELECT_PER_BIN,
        "bins": {},
    }
    filter_fields = [
        "bin", "record_index", "status", "identity", "mutation_count",
        "fixed_mutation_count", "outside_redesignable_mutation_count",
        "header", "fail_reasons", "sequence", "mutations",
    ]
    for bin_id, spec in BINS.items():
        selected = set(base.evenly_select(mutable, spec["target_mutations"]))
        run_dir = runs / f"bin{bin_id}_shell8_seed{spec['seed']}"
        run_dir.mkdir(parents=True, exist_ok=True)
        omit = run_dir / f"bin{bin_id}_omit_native.json"
        base.write_omit_native(omit, residues, selected)
        redesigned = " ".join(f"{chain}{resid}" for chain, resid in sorted(selected, key=lambda item: item[1]))
        cmd = [
            str(ligand_py), str(ligand_run),
            "--model_type", "ligand_mpnn",
            "--seed", str(spec["seed"]),
            "--pdb_path", str(PDB),
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
        fasta = run_dir / "seqs" / f"{PDB.stem}.fa"
        rows = base.filter_records(bin_id, base.parse_fasta(fasta), residues, ref_seq, selected)
        filter_tsv = run_dir / f"filter_bin{bin_id}.tsv"
        with filter_tsv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, filter_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        passed = [row for row in rows if row["status"] == "pass"][:SELECT_PER_BIN]
        for rank, row in enumerate(passed, 1):
            sample_id = f"baker_reference_shell8_bin{bin_id}_rank{rank}_rec{row['record_index']}"
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
            })
        summary["bins"][bin_id] = {
            "exit_code": proc.returncode,
            "elapsed_s": round(time.time() - start, 2),
            "target_mutations": spec["target_mutations"],
            "redesigned_residues": len(selected),
            "records": len(rows),
            "pass_filter_count": len([row for row in rows if row["status"] == "pass"]),
            "selected_count": len(passed),
            "run_dir": str(run_dir),
        }
        SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    fields = [
        "sample_id", "bin", "identity", "mutation_count", "sequence_length",
        "sequence", "postseq_predicted_pdb", "record_index", "fixed_mutation_count",
        "mutations", "header",
    ]
    with MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected_all)
    summary.update({
        "status": "DONE",
        "selected_total": len(selected_all),
        "selected_tsv": str(MANIFEST),
        "finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
    })
    SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
