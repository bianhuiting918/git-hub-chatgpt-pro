#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import time
from pathlib import Path

AA3 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C",
    "GLN": "Q", "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I",
    "LEU": "L", "LYS": "K", "MET": "M", "PHE": "F", "PRO": "P",
    "SER": "S", "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}
BINS = {
    "90": {"target_mutations": 16, "seed": 93201, "temperature": 0.10, "batch_size": 10, "batches": 3},
    "80": {"target_mutations": 32, "seed": 83201, "temperature": 0.12, "batch_size": 10, "batches": 3},
    "70": {"target_mutations": 48, "seed": 73201, "temperature": 0.15, "batch_size": 10, "batches": 3},
    "60": {"target_mutations": 64, "seed": 63201, "temperature": 0.18, "batch_size": 10, "batches": 3},
    "50": {"target_mutations": 80, "seed": 53201, "temperature": 0.22, "batch_size": 10, "batches": 3},
}
FIXED = {("A", r) for r in (13, 14, 15, 16, 17, 54, 55, 56, 72, 73, 74, 148, 149, 150)}


def parse_pdb_sequence(path: Path) -> tuple[list[tuple[str, int, str]], str]:
    residues = []
    seen = set()
    for line in path.read_text(errors="ignore").splitlines():
        if not line.startswith("ATOM") or (line[21].strip() or "A") != "A":
            continue
        key = (line[21].strip() or "A", int(line[22:26]), line[26].strip())
        if key in seen:
            continue
        seen.add(key)
        residues.append((key[0], key[1], AA3.get(line[17:20].strip(), "X")))
    return residues, "".join(item[2] for item in residues)


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    header = None
    seq_parts = []
    if not path.exists():
        return records
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_parts)))
            header = line[1:]
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        records.append((header, "".join(seq_parts)))
    return records


def evenly_select(items: list[tuple[str, int]], n: int) -> list[tuple[str, int]]:
    if n >= len(items):
        return list(items)
    if n <= 1:
        return items[:n]
    return [items[round(i * (len(items) - 1) / (n - 1))] for i in range(n)]


def write_omit_native(path: Path, residues: list[tuple[str, int, str]], selected: set[tuple[str, int]]) -> None:
    data = {f"{chain}{resid}": aa for chain, resid, aa in residues if (chain, resid) in selected}
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def filter_records(bin_id: str, records: list[tuple[str, str]], residues, ref_seq: str, selected: set[tuple[str, int]]):
    rows = []
    seen = set()
    target = BINS[bin_id]["target_mutations"]
    for record_index, (header, seq) in enumerate(records):
        reasons = []
        mutations = []
        fixed_muts = []
        outside = []
        if len(seq) != len(ref_seq):
            reasons.append("length_mismatch")
        else:
            for (chain, resid, ref_aa), aa in zip(residues, seq):
                if aa == ref_aa:
                    continue
                label = f"{chain}{resid}:{ref_aa}>{aa}"
                mutations.append(label)
                if (chain, resid) in FIXED:
                    fixed_muts.append(label)
                if (chain, resid) not in selected:
                    outside.append(label)
        if len(mutations) != target:
            reasons.append("mutation_count_not_target")
        if fixed_muts:
            reasons.append("fixed_motif_mutation")
        if outside:
            reasons.append("outside_redesigned_mutation")
        if seq in seen:
            reasons.append("duplicate_sequence")
        seen.add(seq)
        rows.append({
            "bin": bin_id,
            "record_index": record_index,
            "status": "pass" if not reasons else "fail",
            "identity": round((len(ref_seq) - len(mutations)) / len(ref_seq), 6),
            "mutation_count": len(mutations),
            "fixed_mutation_count": len(fixed_muts),
            "outside_redesignable_mutation_count": len(outside),
            "header": header,
            "fail_reasons": ";".join(reasons),
            "sequence": seq,
            "mutations": ";".join(mutations),
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default="/data/bht/project01_baker_serhyd_routeB_20260701")
    parser.add_argument("--pdb", default="external/serine-hydrolase-design/design_pipeline/04_af2/outputs/super_design_model_4_ptm_seed_0_unrelaxed.pdb")
    parser.add_argument("--out-dir", default="sequence_design/baker_reference_ligandmpnn_bins")
    args = parser.parse_args()

    root = Path(args.root)
    pdb = root / args.pdb
    out = root / args.out_dir
    runs = out / "runs"
    manifests = root / "manifests"
    pred_dir = root / "structure_prediction/esmfold_baker_reference_bins"
    runs.mkdir(parents=True, exist_ok=True)
    manifests.mkdir(parents=True, exist_ok=True)
    pred_dir.mkdir(parents=True, exist_ok=True)
    ligand_py = Path("/data/bht/design_tools/envs/ligandmpnn_venv/bin/python")
    ligand_run = Path("/data/bht/design_tools/src/LigandMPNN/run.py")
    residues, ref_seq = parse_pdb_sequence(pdb)
    mutable = [(chain, resid) for chain, resid, aa in residues if (chain, resid) not in FIXED and aa != "X"]
    selected_all = []
    summary = {
        "status": "RUNNING",
        "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "reference_pdb": str(pdb),
        "sequence_length": len(ref_seq),
        "fixed_motif_residues": len(FIXED),
        "mutable_positions": len(mutable),
        "bins": {},
    }
    filter_fields = [
        "bin", "record_index", "status", "identity", "mutation_count",
        "fixed_mutation_count", "outside_redesignable_mutation_count",
        "header", "fail_reasons", "sequence", "mutations",
    ]
    for bin_id, spec in BINS.items():
        selected = set(evenly_select(mutable, spec["target_mutations"]))
        run_dir = runs / f"bin{bin_id}_baker_reference_seed{spec['seed']}"
        run_dir.mkdir(parents=True, exist_ok=True)
        omit = run_dir / f"bin{bin_id}_omit_native.json"
        write_omit_native(omit, residues, selected)
        redesigned = " ".join(f"{chain}{resid}" for chain, resid in sorted(selected, key=lambda item: item[1]))
        cmd = [
            str(ligand_py), str(ligand_run),
            "--model_type", "ligand_mpnn",
            "--seed", str(spec["seed"]),
            "--pdb_path", str(pdb),
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
        fasta = run_dir / "seqs" / f"{pdb.stem}.fa"
        rows = filter_records(bin_id, parse_fasta(fasta), residues, ref_seq, selected)
        filter_tsv = run_dir / f"filter_bin{bin_id}.tsv"
        with filter_tsv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, filter_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        passed = [row for row in rows if row["status"] == "pass"][:10]
        for rank, row in enumerate(passed, 1):
            sample_id = f"baker_reference_bin{bin_id}_rank{rank}_rec{row['record_index']}"
            selected_all.append({
                "sample_id": sample_id,
                "bin": row["bin"],
                "identity": row["identity"],
                "mutation_count": row["mutation_count"],
                "sequence_length": len(row["sequence"]),
                "sequence": row["sequence"],
                "postseq_predicted_pdb": str(pred_dir / f"{sample_id}.esmfold.pdb"),
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
            "fasta": str(fasta),
            "filter_tsv": str(filter_tsv),
            "log": str(log),
        }
        (manifests / "baker_reference_ligandmpnn_bins_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    selected_tsv = manifests / "baker_reference_ligandmpnn_bins_selected.tsv"
    selected_fields = [
        "sample_id", "bin", "identity", "mutation_count", "sequence_length",
        "sequence", "postseq_predicted_pdb", "record_index", "fixed_mutation_count",
        "mutations", "header",
    ]
    with selected_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, selected_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(selected_all)
    summary.update({
        "status": "DONE",
        "selected_total": len(selected_all),
        "selected_tsv": str(selected_tsv),
        "finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
    })
    (manifests / "baker_reference_ligandmpnn_bins_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
