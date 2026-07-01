#!/usr/bin/env python3
from __future__ import annotations

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

ROOT = Path("/data/bht/project01_baker_serhyd_routeB_20260701")
PDB = ROOT / "outputs/ca_rfd_baker_theozyme_diffusion_shared_n1_publicckpt_v2/sample_21_refined_0.pdb"
OUT = ROOT / "sequence_design/sample21_parent_round02_diverse"
MAN = ROOT / "manifests/baker_sample21_parent_round02_candidates.tsv"
SUMMARY = ROOT / "manifests/baker_sample21_parent_round02_summary.json"
FIXED = {13, 14, 15, 16, 17, 54, 55, 56, 72, 73, 74, 148, 149, 150}

RUN = Path("/data/bht/design_tools/src/LigandMPNN/run.py")
PY = Path("/data/bht/design_tools/envs/ligandmpnn_venv/bin/python")

STRATEGIES = [
    {
        "name": "soluble_bias_t020",
        "model_type": "soluble_mpnn",
        "seed": "43001",
        "temperature": "0.20",
        "batch_size": "10",
        "number_of_batches": "3",
        "bias_AA": "A:-1.2,G:-0.3,P:-0.8,C:-4.0",
    },
    {
        "name": "protein_bias_t025",
        "model_type": "protein_mpnn",
        "seed": "43002",
        "temperature": "0.25",
        "batch_size": "10",
        "number_of_batches": "3",
        "bias_AA": "A:-1.0,G:-0.2,P:-0.6,C:-4.0",
    },
    {
        "name": "ligand_bias_t030",
        "model_type": "ligand_mpnn",
        "seed": "43003",
        "temperature": "0.30",
        "batch_size": "10",
        "number_of_batches": "4",
        "bias_AA": "A:-1.0,G:-0.2,P:-0.6,C:-4.0",
    },
]


def residues(path: Path) -> tuple[list[tuple[str, int, str]], str]:
    seen = set()
    out = []
    for line in path.read_text().splitlines():
        if not line.startswith("ATOM") or (line[21].strip() or "A") != "A":
            continue
        key = (int(line[22:26]), line[26].strip())
        if key in seen:
            continue
        seen.add(key)
        out.append(("A", key[0], AA3.get(line[17:20].strip(), "X")))
    return out, "".join(item[2] for item in out)


def fasta_records(path: Path) -> list[tuple[str, str]]:
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


def main() -> int:
    res, ref = residues(PDB)
    fixed = " ".join(f"A{r}" for r in sorted(FIXED))
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    run_summaries = []
    t0 = time.time()

    for strategy in STRATEGIES:
        run_out = OUT / strategy["name"]
        run_out.mkdir(parents=True, exist_ok=True)
        cmd = [
            str(PY), str(RUN),
            "--model_type", strategy["model_type"],
            "--seed", strategy["seed"],
            "--pdb_path", str(PDB),
            "--out_folder", str(run_out),
            "--chains_to_design", "A",
            "--parse_these_chains_only", "A",
            "--fixed_residues", fixed,
            "--temperature", strategy["temperature"],
            "--batch_size", strategy["batch_size"],
            "--number_of_batches", strategy["number_of_batches"],
            "--bias_AA", strategy["bias_AA"],
            "--save_stats", "0",
            "--verbose", "0",
        ]
        log = run_out / "ligandmpnn_parent_round02.log"
        start = time.time()
        with log.open("w") as handle:
            proc = subprocess.run(cmd, cwd=str(RUN.parent), stdout=handle, stderr=subprocess.STDOUT)
        fa = run_out / "seqs" / f"{PDB.stem}.fa"
        records = fasta_records(fa)
        produced = 0
        for record_index, (header, seq) in enumerate(records):
            if record_index == 0:
                continue
            if len(seq) != len(ref):
                continue
            mutations = []
            fixed_mutations = []
            for (chain, resid, aa), new in zip(res, seq):
                if new != aa:
                    label = f"{chain}{resid}:{aa}>{new}"
                    mutations.append(label)
                    if resid in FIXED:
                        fixed_mutations.append(label)
            sample_id = f"baker_sample21_parent_round02_{strategy['name']}_rank{produced + 1}_global{len(rows) + 1}"
            rows.append({
                "sample_id": sample_id,
                "strategy": strategy["name"],
                "model_type": strategy["model_type"],
                "temperature": strategy["temperature"],
                "identity_to_allA_ref": round((len(ref) - len(mutations)) / len(ref), 6),
                "mutation_count_vs_allA_ref": len(mutations),
                "fixed_mutation_count": len(fixed_mutations),
                "fixed_mutations": ";".join(fixed_mutations),
                "sequence_length": len(seq),
                "sequence": seq,
                "ligandmpnn_header": header,
                "mutations_vs_allA_ref": ";".join(mutations),
                "postseq_predicted_pdb": str(ROOT / "structure_prediction/esmfold_sample21_parent_round02" / f"{sample_id}.esmfold.pdb"),
            })
            produced += 1
        run_summaries.append({
            "strategy": strategy["name"],
            "exit_code": proc.returncode,
            "elapsed_s": round(time.time() - start, 2),
            "records_after_reference": produced,
            "log": str(log),
            "fasta": str(fa),
        })

    fields = list(rows[0].keys()) if rows else ["sample_id"]
    with MAN.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "status": "DONE",
        "candidate_count": len(rows),
        "elapsed_s": round(time.time() - t0, 2),
        "fixed_residues": sorted(FIXED),
        "gate_target": "downstream parent active-pocket PASS: motif CA RMSD <= 1.0 A and max motif pair-distance delta <= 1.0 A",
        "manifest": str(MAN),
        "runs": run_summaries,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
