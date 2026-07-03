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


TARGETS = {"90": 0.90, "80": 0.80, "70": 0.70, "60": 0.60, "50": 0.50}


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


def identity(a: str, b: str) -> float:
    if len(a) != len(b) or not a:
        return 0.0
    return sum(x == y for x, y in zip(a, b)) / len(a)


def bin_for_identity(value: float, tolerance: float) -> str | None:
    best = min(TARGETS, key=lambda item: abs(value - TARGETS[item]))
    if abs(value - TARGETS[best]) <= tolerance:
        return best
    return None


def run_ligandmpnn(
    ligand_py: Path,
    ligand_run: Path,
    parent_pdb: Path,
    out_dir: Path,
    redesigned: str,
    seed: int,
    temperature: float,
    batch_size: int,
    batches: int,
) -> tuple[int, Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    log = out_dir / "ligandmpnn.log"
    cmd = [
        str(ligand_py), str(ligand_run),
        "--model_type", "ligand_mpnn",
        "--seed", str(seed),
        "--pdb_path", str(parent_pdb),
        "--out_folder", str(out_dir),
        "--chains_to_design", "A",
        "--parse_these_chains_only", "A",
        "--redesigned_residues", redesigned,
        "--temperature", str(temperature),
        "--batch_size", str(batch_size),
        "--number_of_batches", str(batches),
        "--save_stats", "0",
        "--verbose", "0",
    ]
    with log.open("w") as handle:
        proc = subprocess.run(cmd, cwd=str(ligand_run.parent), stdout=handle, stderr=subprocess.STDOUT)
    fasta = out_dir / "seqs" / f"{parent_pdb.stem}.fa"
    return proc.returncode, fasta, log


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", default=str(ROOT))
    parser.add_argument("--parent-pdb", default="outputs/ca_rfd_baker_pocket4_l2_extend_expanded_n1_20260703/sample_7500.pdb")
    parser.add_argument("--out-dir", default="sequence_design/pocket4_l2_refseq_ligandmpnn")
    parser.add_argument("--manifest-tsv", default="manifests/pocket4_l2_refseq_ligandmpnn_selected.tsv")
    parser.add_argument("--summary-json", default="manifests/pocket4_l2_refseq_ligandmpnn_summary.json")
    parser.add_argument("--fixed-residues", default="61,87", help="Comma-separated output residue ids fixed during design; default catalytic Ser/His")
    parser.add_argument("--temps", default="0.05,0.08,0.10,0.15,0.20,0.30,0.45")
    parser.add_argument("--batch-size", type=int, default=160)
    parser.add_argument("--batches", type=int, default=3)
    parser.add_argument("--seed-base", type=int, default=2704000)
    parser.add_argument("--identity-tolerance", type=float, default=0.055)
    parser.add_argument("--select-per-bin", type=int, default=120)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()

    root = Path(args.root)
    parent_pdb = root / args.parent_pdb if not Path(args.parent_pdb).is_absolute() else Path(args.parent_pdb)
    out_dir = root / args.out_dir
    manifest_tsv = root / args.manifest_tsv
    summary_json = root / args.summary_json
    pred_dir = root / "structure_prediction/esmfold_pocket4_l2_refseq_bins"
    out_dir.mkdir(parents=True, exist_ok=True)
    pred_dir.mkdir(parents=True, exist_ok=True)

    fixed = {("A", int(item.strip())) for item in args.fixed_residues.split(",") if item.strip()}
    residues, parent_seq = base.parse_pdb_sequence(parent_pdb)
    designable = [(chain, resid) for chain, resid, aa in residues if (chain, resid) not in fixed and aa != "X"]
    redesigned = " ".join(f"{chain}{resid}" for chain, resid in designable)
    ligand_py = Path("/data/bht/design_tools/envs/ligandmpnn_venv/bin/python")
    ligand_run = Path("/data/bht/design_tools/src/LigandMPNN/run.py")
    temps = [float(item) for item in args.temps.split(",") if item.strip()]

    summary = {
        "status": "RUNNING",
        "started": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "route": "Route B2 L2 reference-sequence LigandMPNN design",
        "parent_pdb": str(parent_pdb),
        "sequence_length": len(parent_seq),
        "fixed_residues": sorted(f"A{resid}" for _chain, resid in fixed),
        "designable_positions": len(designable),
        "identity_reference": "PENDING",
        "identity_tolerance": args.identity_tolerance,
        "temps": temps,
        "runs": {},
    }
    summary_json.parent.mkdir(parents=True, exist_ok=True)
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    all_records: list[dict[str, object]] = []
    seen: set[str] = set()
    reference_seq = ""
    for idx, temp in enumerate(temps):
        run_dir = out_dir / "runs" / f"temp{str(temp).replace('.', 'p')}_seed{args.seed_base + idx}"
        fasta = run_dir / "seqs" / f"{parent_pdb.stem}.fa"
        log = run_dir / "ligandmpnn.log"
        if fasta.exists() and not args.force:
            returncode = "SKIP_EXISTS"
        else:
            start = time.time()
            returncode, fasta, log = run_ligandmpnn(
                ligand_py,
                ligand_run,
                parent_pdb,
                run_dir,
                redesigned,
                args.seed_base + idx,
                temp,
                args.batch_size,
                args.batches,
            )
            elapsed_s = round(time.time() - start, 2)
        records = base.parse_fasta(fasta)
        valid = []
        for rec_index, (header, seq) in enumerate(records):
            if len(seq) != len(parent_seq) or seq == parent_seq or seq in seen:
                continue
            fixed_ok = all(seq[resid - 1] == parent_seq[resid - 1] for _chain, resid in fixed)
            if not fixed_ok:
                continue
            seen.add(seq)
            valid.append((rec_index, header, seq))
            if not reference_seq:
                reference_seq = seq
        for rec_index, header, seq in valid:
            if not reference_seq:
                continue
            ident = identity(seq, reference_seq)
            bin_id = bin_for_identity(ident, args.identity_tolerance)
            all_records.append({
                "sample_id": "",
                "bin": "" if bin_id is None else bin_id,
                "identity": round(ident, 6),
                "sequence_length": len(seq),
                "sequence": seq,
                "postseq_predicted_pdb": "",
                "record_index": rec_index,
                "temperature": temp,
                "header": header,
                "source_fasta": str(fasta),
            })
        summary["identity_reference"] = reference_seq
        summary["runs"][str(temp)] = {
            "returncode": returncode,
            "fasta": str(fasta),
            "log": str(log),
            "records_total": len(records),
            "records_valid_unique": len(valid),
            "elapsed_s": locals().get("elapsed_s", 0),
        }
        summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    selected: list[dict[str, object]] = []
    selected_counts = {bin_id: 0 for bin_id in TARGETS}
    for row in sorted(all_records, key=lambda item: (str(item["bin"]) == "", str(item["bin"]), abs(float(item["identity"]) - TARGETS.get(str(item["bin"]), 0.0)))):
        bin_id = str(row["bin"])
        if bin_id not in selected_counts or selected_counts[bin_id] >= args.select_per_bin:
            continue
        selected_counts[bin_id] += 1
        sample_id = f"pocket4_l2_refseq_bin{bin_id}_rank{selected_counts[bin_id]}_t{str(row['temperature']).replace('.', 'p')}_rec{row['record_index']}"
        row["sample_id"] = sample_id
        row["postseq_predicted_pdb"] = str(pred_dir / f"{sample_id}.esmfold.pdb")
        selected.append(row)

    fields = [
        "sample_id", "bin", "identity", "sequence_length", "sequence",
        "postseq_predicted_pdb", "record_index", "temperature", "header", "source_fasta",
    ]
    write_tsv(manifest_tsv, fields, selected)
    all_counts = {bin_id: 0 for bin_id in TARGETS}
    for row in all_records:
        if row["bin"] in all_counts:
            all_counts[str(row["bin"])] += 1
    summary.update({
        "status": "DONE",
        "finished": time.strftime("%Y-%m-%dT%H:%M:%S"),
        "selected_tsv": str(manifest_tsv),
        "selected_total": len(selected),
        "selected_counts_by_bin": selected_counts,
        "candidate_counts_by_bin": all_counts,
        "candidate_total": len(all_records),
    })
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
