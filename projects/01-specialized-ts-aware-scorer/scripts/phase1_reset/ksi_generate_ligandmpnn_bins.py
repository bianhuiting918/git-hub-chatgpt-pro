#!/usr/bin/env python3
"""Generate KSI natural-scaffold sequence bins with LigandMPNN.

Inputs are the PyFAMSA-derived fixed/mutable masks. PLACER is intentionally not
used here; this is a sequence-generation stage only.
"""
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
import time
from pathlib import Path

AA3 = {
    "ALA":"A","ARG":"R","ASN":"N","ASP":"D","CYS":"C","GLN":"Q","GLU":"E","GLY":"G","HIS":"H","ILE":"I",
    "LEU":"L","LYS":"K","MET":"M","PHE":"F","PRO":"P","SER":"S","THR":"T","TRP":"W","TYR":"Y","VAL":"V",
}

BIN_SPECS = {
    "90": {"target_mutations": 13, "mutation_min": 10, "mutation_max": 16, "seed": 9101, "temperature": 0.10, "batch_size": 10, "batches": 4},
    "80": {"target_mutations": 26, "mutation_min": 23, "mutation_max": 30, "seed": 8101, "temperature": 0.10, "batch_size": 10, "batches": 3},
    "70": {"target_mutations": 39, "mutation_min": 36, "mutation_max": 43, "seed": 7101, "temperature": 0.12, "batch_size": 10, "batches": 4},
    "60": {"target_mutations": 52, "mutation_min": 49, "mutation_max": 56, "seed": 6101, "temperature": 0.15, "batch_size": 10, "batches": 5},
    "50": {"target_mutations": 65, "mutation_min": 62, "mutation_max": 69, "seed": 5101, "temperature": 0.18, "batch_size": 10, "batches": 6},
}


def parse_mask_tsv(path: Path) -> list[tuple[str, int]]:
    rows: list[tuple[str, int]] = []
    with path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(("A", int(row["ref_pos"])))
    return rows


def parse_pdb_sequence(path: Path, chain: str) -> tuple[list[tuple[str, int, str]], str]:
    residues = []
    seen = set()
    with path.open() as handle:
        for line in handle:
            if not line.startswith("ATOM") or line[21].strip() != chain:
                continue
            key = (line[21].strip(), int(line[22:26]), line[26].strip())
            if key in seen:
                continue
            seen.add(key)
            residues.append((key[0], key[1], AA3.get(line[17:20].strip(), "X")))
    return residues, "".join(aa for _, _, aa in residues)


def parse_fasta(path: Path) -> list[tuple[str, str]]:
    records = []
    header = None
    seq_parts = []
    with path.open() as handle:
        for raw in handle:
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


def write_omit_native(path: Path, residues: list[tuple[str, int, str]], selected: set[tuple[str, int]]) -> None:
    data = {f"{chain}{resid}": aa for chain, resid, aa in residues if (chain, resid) in selected}
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n")


def evenly_select(items: list[tuple[str, int]], n: int) -> list[tuple[str, int]]:
    if n >= len(items):
        return list(items)
    if n <= 1:
        return items[:n]
    return [items[round(i * (len(items) - 1) / (n - 1))] for i in range(n)]


def filter_records(bin_id: str, fasta: Path, residues: list[tuple[str, int, str]], ref_seq: str, fixed: set[tuple[str, int]], redesign: set[tuple[str, int]]) -> list[dict[str, object]]:
    spec = BIN_SPECS[bin_id]
    rows: list[dict[str, object]] = []
    seen = set()
    for idx, (header, seq) in enumerate(parse_fasta(fasta)):
        reasons = []
        mutations = []
        fixed_mut = []
        outside_mut = []
        if len(seq) != len(ref_seq):
            rows.append({"bin": bin_id, "record_index": idx, "status": "fail", "identity": "", "mutation_count": "", "fixed_mutation_count": "", "outside_redesignable_mutation_count": "", "header": header, "fail_reasons": "length_mismatch", "sequence": seq, "mutations": ""})
            continue
        for (chain, resid, ref_aa), aa in zip(residues, seq):
            if aa == ref_aa:
                continue
            label = f"{chain}{resid}:{ref_aa}>{aa}"
            mutations.append(label)
            if (chain, resid) in fixed:
                fixed_mut.append(label)
            if (chain, resid) not in redesign:
                outside_mut.append(label)
        mut_count = len(mutations)
        if not (spec["mutation_min"] <= mut_count <= spec["mutation_max"]):
            reasons.append("mutation_count_out_of_bin")
        if fixed_mut:
            reasons.append("fixed_residue_mutation")
        if outside_mut:
            reasons.append("outside_redesignable_mutation")
        if seq in seen:
            reasons.append("duplicate_sequence")
        seen.add(seq)
        rows.append({
            "bin": bin_id,
            "record_index": idx,
            "status": "pass" if not reasons else "fail",
            "identity": round((len(ref_seq) - mut_count) / len(ref_seq), 6),
            "mutation_count": mut_count,
            "fixed_mutation_count": len(fixed_mut),
            "outside_redesignable_mutation_count": len(outside_mut),
            "header": header,
            "fail_reasons": ";".join(reasons),
            "sequence": seq,
            "mutations": ";".join(mutations),
        })
    return rows


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default="/data/bht/project01_phase1_reset_gpu")
    parser.add_argument("--pdb", default="natural_scaffold/KSI/refs/6UBQ.pdb")
    parser.add_argument("--mutable-tsv", default="natural_scaffold/KSI/manifests/natural_scaffold_msa_pyfamsa_mutable_positions.tsv")
    parser.add_argument("--fixed-tsv", default="natural_scaffold/KSI/manifests/natural_scaffold_msa_pyfamsa_fixed_positions.tsv")
    parser.add_argument("--out-dir", default="natural_scaffold/KSI/sequence_generation_chainA")
    parser.add_argument("--chain", default="A")
    args = parser.parse_args()

    root = Path(args.root)
    pdb = root / args.pdb
    out_dir = root / args.out_dir
    runs_dir = out_dir / "ligandmpnn_runs"
    manifest_dir = root / "natural_scaffold/KSI/manifests"
    runs_dir.mkdir(parents=True, exist_ok=True)
    manifest_dir.mkdir(parents=True, exist_ok=True)

    ligand_py = Path("/data/bht/design_tools/envs/ligandmpnn_venv/bin/python")
    ligand_run = Path("/data/bht/design_tools/src/LigandMPNN/run.py")
    ligand_cwd = ligand_run.parent

    residues, ref_seq = parse_pdb_sequence(pdb, args.chain)
    mutable = parse_mask_tsv(root / args.mutable_tsv)
    fixed = set(parse_mask_tsv(root / args.fixed_tsv))
    mutable_set = set(mutable)
    summary = {"started": time.strftime("%Y-%m-%dT%H:%M:%S"), "reference_pdb": str(pdb), "sequence_length": len(ref_seq), "mutable_positions": len(mutable), "fixed_positions": len(fixed), "bins": {}}
    all_pass_rows = []

    fields = ["bin", "record_index", "status", "identity", "mutation_count", "fixed_mutation_count", "outside_redesignable_mutation_count", "header", "fail_reasons", "sequence", "mutations"]
    for bin_id, spec in BIN_SPECS.items():
        selected = set(evenly_select(mutable, spec["target_mutations"]))
        run_dir = runs_dir / f"bin{bin_id}_ksi_seed{spec['seed']}"
        run_dir.mkdir(parents=True, exist_ok=True)
        omit_json = run_dir / f"bin{bin_id}_omit_native.json"
        write_omit_native(omit_json, residues, selected)
        redesigned = " ".join(f"{chain}{resid}" for chain, resid in sorted(selected, key=lambda x: x[1]))
        cmd = [str(ligand_py), str(ligand_run), "--model_type", "ligand_mpnn", "--seed", str(spec["seed"]), "--pdb_path", str(pdb), "--out_folder", str(run_dir), "--chains_to_design", args.chain, "--parse_these_chains_only", args.chain, "--redesigned_residues", redesigned, "--temperature", str(spec["temperature"]), "--batch_size", str(spec["batch_size"]), "--number_of_batches", str(spec["batches"]), "--omit_AA_per_residue", str(omit_json), "--save_stats", "0", "--verbose", "0"]
        log = run_dir / "ligandmpnn.log"
        t0 = time.time()
        with log.open("w") as handle:
            proc = subprocess.run(cmd, cwd=str(ligand_cwd), stdout=handle, stderr=subprocess.STDOUT)
        fasta = run_dir / "seqs" / f"{pdb.stem}.fa"
        rows = filter_records(bin_id, fasta, residues, ref_seq, fixed, selected) if fasta.exists() else []
        filter_tsv = run_dir / f"filter{bin_id}_ksi.tsv"
        with filter_tsv.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        pass_rows = [row for row in rows if row["status"] == "pass"]
        for rank, row in enumerate(pass_rows, start=1):
            all_pass_rows.append({"sample_id": f"ksi_bin{bin_id}_rank{rank}_rec{row['record_index']}", **row})
        summary["bins"][bin_id] = {"exit_code": proc.returncode, "elapsed_s": round(time.time() - t0, 2), "target_mutations": spec["target_mutations"], "redesigned_residues": len(selected), "records": len(rows), "pass_filter_count": len(pass_rows), "run_dir": str(run_dir), "filter_tsv": str(filter_tsv), "fasta": str(fasta)}

    selected_tsv = manifest_dir / "ksi_ligandmpnn_sequence_bins_selected.tsv"
    selected_fields = ["sample_id"] + fields
    with selected_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, selected_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(all_pass_rows)
    selected_fasta = manifest_dir / "ksi_ligandmpnn_sequence_bins_selected.fasta"
    with selected_fasta.open("w") as handle:
        for row in all_pass_rows:
            handle.write(f">{row['sample_id']}|bin={row['bin']}|identity={row['identity']}|mutations={row['mutation_count']}\n{row['sequence']}\n")
    summary["selected_total"] = len(all_pass_rows)
    summary["selected_tsv"] = str(selected_tsv)
    summary["selected_fasta"] = str(selected_fasta)
    summary["finished"] = time.strftime("%Y-%m-%dT%H:%M:%S")
    summary_path = manifest_dir / "ksi_ligandmpnn_sequence_bins_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


