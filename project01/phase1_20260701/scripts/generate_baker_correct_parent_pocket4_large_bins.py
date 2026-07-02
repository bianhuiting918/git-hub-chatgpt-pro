#!/usr/bin/env python3
from __future__ import annotations

import csv
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path("/data/bht/project01_baker_serhyd_routeB_20260701")
SCRIPT_DIR = ROOT / "scripts"
sys.path.insert(0, str(SCRIPT_DIR))

import generate_baker_reference_ligandmpnn_bins as base  # noqa: E402

PARENT_PDB = ROOT / "external/serine-hydrolase-design/design_pipeline/04_af2/outputs/refined_out_1_bn1_1_20.pdb"
LIGAND_REF_PDB = ROOT / "external/serine-hydrolase-design/design_pipeline/03_design/outputs/refined_out_1_bn1_1.pdb"
OUT = ROOT / "sequence_design/baker_correct_parent_pocket4_large_ligandmpnn_bins"
MANIFEST = ROOT / "manifests/baker_correct_parent_pocket4_large_ligandmpnn_bins_selected.tsv"
SUMMARY = ROOT / "manifests/baker_correct_parent_pocket4_large_ligandmpnn_bins_summary.json"
PRED_DIR = ROOT / "structure_prediction/esmfold_baker_correct_parent_pocket4_large_bins"

CATALYTIC_AND_OXY = {15, 16, 55, 73, 149}
POCKET_HEAVY_CUTOFF_A = 4.0
SELECT_PER_BIN = None
BINS = {
    "90": {"target_mutations": 16, "seed": 98301, "temperature": 0.10, "batch_size": 100, "batches": 4, "select": 200},
    "80": {"target_mutations": 32, "seed": 88301, "temperature": 0.12, "batch_size": 250, "batches": 5, "select": 1000},
    "70": {"target_mutations": 48, "seed": 78301, "temperature": 0.15, "batch_size": 250, "batches": 5, "select": 1000},
    "60": {"target_mutations": 64, "seed": 68301, "temperature": 0.18, "batch_size": 250, "batches": 5, "select": 1000},
    "50": {"target_mutations": 80, "seed": 58301, "temperature": 0.22, "batch_size": 250, "batches": 5, "select": 1000},
}


def heavy_xyz(line: str) -> tuple[float, float, float]:
    return float(line[30:38]), float(line[38:46]), float(line[46:54])


def ligand_4a_fixed_residues(path: Path) -> set[int]:
    protein = []
    ligand = []
    for line in path.read_text(errors="ignore").splitlines():
        rec = line[:6].strip()
        atom = line[12:16].strip()
        if atom.startswith("H"):
            continue
        if rec == "ATOM" and (line[21].strip() or "A") == "A":
            protein.append((int(line[22:26]), heavy_xyz(line)))
        elif rec == "HETATM" and line[17:20].strip() == "bn1":
            ligand.append(heavy_xyz(line))
    fixed = set(CATALYTIC_AND_OXY)
    for resid, xyz in protein:
        if any(math.dist(xyz, lig_xyz) <= POCKET_HEAVY_CUTOFF_A for lig_xyz in ligand):
            fixed.add(resid)
    return fixed


def ca_coords(path: Path) -> dict[int, np.ndarray]:
    coords = {}
    for line in path.read_text(errors="ignore").splitlines():
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and (line[21].strip() or "A") == "A":
            coords[int(line[22:26])] = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    return coords


def farthest_mutable_positions(residues, fixed: set[tuple[str, int]], coords: dict[int, np.ndarray], count: int):
    fixed_resids = {resid for _chain, resid in fixed}
    scored = []
    for chain, resid, aa in residues:
        if (chain, resid) in fixed or aa == "X" or resid not in coords:
            continue
        distance = min(np.linalg.norm(coords[resid] - coords[pocket]) for pocket in fixed_resids if pocket in coords)
        scored.append((distance, chain, resid))
    scored.sort(reverse=True)
    return {(chain, resid) for _distance, chain, resid in scored[:count]}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    PRED_DIR.mkdir(parents=True, exist_ok=True)
    runs = OUT / "runs"
    runs.mkdir(parents=True, exist_ok=True)

    fixed_residues = ligand_4a_fixed_residues(LIGAND_REF_PDB)
    fixed = {("A", resid) for resid in fixed_residues}
    base.FIXED = fixed
    base.BINS = BINS
    residues, ref_seq = base.parse_pdb_sequence(PARENT_PDB)
    coords = ca_coords(PARENT_PDB)

    ligand_py = Path("/data/bht/design_tools/envs/ligandmpnn_venv/bin/python")
    ligand_run = Path("/data/bht/design_tools/src/LigandMPNN/run.py")
    selected_all = []
    summary = {
        "status": "RUNNING",
        "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "reference_parent_pdb": str(PARENT_PDB),
        "ligand_reference_pdb": str(LIGAND_REF_PDB),
        "ligand_resname": "bn1",
        "sequence_length": len(ref_seq),
        "fixed_policy": "bn1 heavy-atom 4 A protein residues plus catalytic/oxyanion residues",
        "fixed_residue_count": len(fixed),
        "fixed_residues": sorted(fixed_residues),
        "selected_per_bin": {bin_id: spec["select"] for bin_id, spec in BINS.items()},
        "bins": {},
    }
    filter_fields = [
        "bin", "record_index", "status", "identity", "mutation_count",
        "fixed_mutation_count", "outside_redesignable_mutation_count",
        "header", "fail_reasons", "sequence", "mutations",
    ]
    for bin_id, spec in BINS.items():
        selected = farthest_mutable_positions(residues, fixed, coords, spec["target_mutations"])
        run_dir = runs / f"bin{bin_id}_pocket4_seed{spec['seed']}"
        run_dir.mkdir(parents=True, exist_ok=True)
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
        fasta = run_dir / "seqs" / f"{PARENT_PDB.stem}.fa"
        rows = base.filter_records(bin_id, base.parse_fasta(fasta), residues, ref_seq, selected)
        filter_tsv = run_dir / f"filter_bin{bin_id}.tsv"
        with filter_tsv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, filter_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        passed = [row for row in rows if row["status"] == "pass"][: spec["select"]]
        for rank, row in enumerate(passed, 1):
            sample_id = f"baker_correct_parent_pocket4_large_bin{bin_id}_rank{rank}_rec{row['record_index']}"
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
            "filter_tsv": str(filter_tsv),
            "log": str(log),
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
