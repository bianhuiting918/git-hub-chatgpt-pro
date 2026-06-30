#!/usr/bin/env python3
"""Reference-anchored KSI conservation screen for Project 01 natural scaffolds.

This is a first-pass fallback when a full MSA toolchain is unavailable. It
fetches UniProt FASTA records, filters to KSI-like sequence lengths, aligns each
candidate to the P07445 reference, and reports fixed/mutable positions.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path


CATALYTIC_POSITIONS = {16, 40, 103}


def fetch_uniprot_fasta(query: str, timeout: int = 90) -> str:
    url = (
        "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query="
        + urllib.parse.quote(query)
    )
    with urllib.request.urlopen(url, timeout=timeout) as handle:
        return handle.read().decode("utf-8", "replace")


def parse_fasta(text: str) -> list[dict[str, str]]:
    records: list[dict[str, str]] = []
    header = None
    seq_parts: list[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                records.append({"header": header, "seq": "".join(seq_parts)})
            header = line[1:]
            seq_parts = []
        else:
            seq_parts.append(line)
    if header is not None:
        records.append({"header": header, "seq": "".join(seq_parts)})
    return records


def accession_from_header(header: str) -> str:
    if "|" in header:
        parts = header.split("|")
        if len(parts) >= 2:
            return parts[1]
    return header.split()[0]


def global_align_to_ref(ref: str, query: str) -> tuple[str, str]:
    match_score = 2
    mismatch_score = -1
    gap_score = -4
    n = len(ref)
    m = len(query)
    scores = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[""] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        scores[i][0] = scores[i - 1][0] + gap_score
        trace[i][0] = "U"
    for j in range(1, m + 1):
        scores[0][j] = scores[0][j - 1] + gap_score
        trace[0][j] = "L"
    for i in range(1, n + 1):
        ri = ref[i - 1]
        for j in range(1, m + 1):
            qj = query[j - 1]
            diag = scores[i - 1][j - 1] + (match_score if ri == qj else mismatch_score)
            up = scores[i - 1][j] + gap_score
            left = scores[i][j - 1] + gap_score
            best = max(diag, up, left)
            scores[i][j] = best
            trace[i][j] = "D" if best == diag else ("U" if best == up else "L")
    i, j = n, m
    a_ref: list[str] = []
    a_query: list[str] = []
    while i > 0 or j > 0:
        step = trace[i][j]
        if step == "D":
            a_ref.append(ref[i - 1])
            a_query.append(query[j - 1])
            i -= 1
            j -= 1
        elif step == "U":
            a_ref.append(ref[i - 1])
            a_query.append("-")
            i -= 1
        else:
            a_ref.append("-")
            a_query.append(query[j - 1])
            j -= 1
    return "".join(reversed(a_ref)), "".join(reversed(a_query))


def ref_position_residues(ref: str, query: str) -> tuple[dict[int, str], float, float]:
    a_ref, a_query = global_align_to_ref(ref, query)
    ref_pos = 0
    residues: dict[int, str] = {}
    matches = 0
    covered = 0
    for r, q in zip(a_ref, a_query):
        if r == "-":
            continue
        ref_pos += 1
        if q != "-":
            residues[ref_pos] = q
            covered += 1
            if q == r:
                matches += 1
    identity = matches / len(ref) if ref else math.nan
    coverage = covered / len(ref) if ref else math.nan
    return residues, identity, coverage


def read_ligand_contacts(pdb_path: Path, chain: str, ligand: str, cutoff: float) -> set[int]:
    if not pdb_path.exists():
        return set()
    protein_atoms: list[tuple[int, float, float, float]] = []
    ligand_atoms: list[tuple[float, float, float]] = []
    for line in pdb_path.read_text(errors="replace").splitlines():
        rec = line[:6].strip()
        if rec not in {"ATOM", "HETATM"}:
            continue
        atom_chain = line[21].strip()
        if atom_chain != chain:
            continue
        resname = line[17:20].strip()
        element = line[76:78].strip() or line[12:16].strip()[0]
        if element.upper() == "H":
            continue
        try:
            resid = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
        except ValueError:
            continue
        if rec == "ATOM":
            protein_atoms.append((resid, x, y, z))
        elif resname == ligand:
            ligand_atoms.append((x, y, z))
    contacts: set[int] = set()
    cutoff2 = cutoff * cutoff
    for resid, px, py, pz in protein_atoms:
        for lx, ly, lz in ligand_atoms:
            dx = px - lx
            dy = py - ly
            dz = pz - lz
            if dx * dx + dy * dy + dz * dz <= cutoff2:
                contacts.add(resid)
                break
    return contacts


def write_tsv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write("\t".join(columns) + "\n")
        for row in rows:
            handle.write("\t".join(str(row.get(col, "")) for col in columns) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--reference-query", default="accession:P07445")
    parser.add_argument("--homolog-query", default="(ec:5.3.3.1)")
    parser.add_argument("--pdb", default="")
    parser.add_argument("--chain", default="A")
    parser.add_argument("--ligand", default="ASD")
    parser.add_argument("--contact-cutoff", type=float, default=4.0)
    parser.add_argument("--min-length-ratio", type=float, default=0.55)
    parser.add_argument("--max-length-ratio", type=float, default=1.80)
    parser.add_argument("--min-identity", type=float, default=0.20)
    parser.add_argument("--max-records", type=int, default=200)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    raw_dir = out_dir / "raw"
    manifest_dir = out_dir / "manifests"
    raw_dir.mkdir(parents=True, exist_ok=True)
    manifest_dir.mkdir(parents=True, exist_ok=True)

    ref_text = fetch_uniprot_fasta(args.reference_query)
    homolog_text = fetch_uniprot_fasta(args.homolog_query)
    (raw_dir / "reference_P07445.fasta").write_text(ref_text, encoding="utf-8")
    (raw_dir / "uniprot_ec_5_3_3_1.fasta").write_text(homolog_text, encoding="utf-8")

    ref_records = parse_fasta(ref_text)
    if len(ref_records) != 1:
        raise RuntimeError(f"Expected one reference record, found {len(ref_records)}")
    ref_seq = ref_records[0]["seq"]
    ref_accession = accession_from_header(ref_records[0]["header"])
    raw_homologs = parse_fasta(homolog_text)

    min_len = int(len(ref_seq) * args.min_length_ratio)
    max_len = int(len(ref_seq) * args.max_length_ratio)
    seen: set[str] = set()
    selection_rows: list[dict[str, object]] = []
    selected: list[dict[str, object]] = []
    for record in raw_homologs:
        accession = accession_from_header(record["header"])
        if accession in seen:
            continue
        seen.add(accession)
        seq = re.sub("[^A-Z]", "", record["seq"].upper())
        length_ok = min_len <= len(seq) <= max_len
        residues: dict[int, str] = {}
        identity = math.nan
        coverage = math.nan
        if length_ok:
            residues, identity, coverage = ref_position_residues(ref_seq, seq)
        selected_flag = bool(length_ok and identity >= args.min_identity and coverage >= 0.75)
        reason = "SELECTED" if selected_flag else "NOT_SELECTED"
        if not length_ok:
            reason = "NOT_SELECTED_LENGTH"
        elif identity < args.min_identity:
            reason = "NOT_SELECTED_LOW_IDENTITY"
        elif coverage < 0.75:
            reason = "NOT_SELECTED_LOW_COVERAGE"
        row = {
            "accession": accession,
            "length": len(seq),
            "identity_to_ref": f"{identity:.4f}" if not math.isnan(identity) else "NA",
            "coverage_to_ref": f"{coverage:.4f}" if not math.isnan(coverage) else "NA",
            "selected": "true" if selected_flag else "false",
            "reason": reason,
            "header": record["header"],
        }
        selection_rows.append(row)
        if selected_flag:
            selected.append({"accession": accession, "seq": seq, "residues": residues, "row": row})

    selected.sort(key=lambda x: float(x["row"]["identity_to_ref"]), reverse=True)
    selected = selected[: args.max_records]
    selected_accessions = {x["accession"] for x in selected}
    for row in selection_rows:
        if row["selected"] == "true" and row["accession"] not in selected_accessions:
            row["selected"] = "false"
            row["reason"] = "NOT_SELECTED_MAX_RECORDS"

    ligand_contacts = read_ligand_contacts(
        Path(args.pdb) if args.pdb else Path(),
        args.chain,
        args.ligand,
        args.contact_cutoff,
    )

    conservation_rows: list[dict[str, object]] = []
    fixed_rows: list[dict[str, object]] = []
    mutable_rows: list[dict[str, object]] = []
    for pos, ref_aa in enumerate(ref_seq, start=1):
        counts = Counter()
        for item in selected:
            qaa = item["residues"].get(pos)
            if qaa and qaa != "-":
                counts[qaa] += 1
        total = sum(counts.values())
        major_aa = ""
        major_freq = 0.0
        if total:
            major_aa, major_count = counts.most_common(1)[0]
            major_freq = major_count / total
        is_catalytic = pos in CATALYTIC_POSITIONS
        is_ligand_contact = pos in ligand_contacts
        if is_catalytic:
            cls = "FIXED_CATALYTIC"
        elif is_ligand_contact:
            cls = "FIXED_LIGAND_CONTACT"
        elif major_freq >= 0.80:
            cls = "FIXED_CONSERVED"
        elif major_freq >= 0.60:
            cls = "REVIEW_CONSERVATION"
        else:
            cls = "MUTABLE_BACKGROUND"
        row = {
            "ref_pos": pos,
            "ref_aa": ref_aa,
            "coverage_count": total,
            "major_aa": major_aa,
            "major_freq": f"{major_freq:.4f}",
            "is_catalytic": "true" if is_catalytic else "false",
            "is_ligand_contact": "true" if is_ligand_contact else "false",
            "class": cls,
        }
        conservation_rows.append(row)
        if cls.startswith("FIXED"):
            fixed_rows.append(row)
        elif cls == "MUTABLE_BACKGROUND":
            mutable_rows.append(row)

    columns_selection = [
        "accession",
        "length",
        "identity_to_ref",
        "coverage_to_ref",
        "selected",
        "reason",
        "header",
    ]
    columns_pos = [
        "ref_pos",
        "ref_aa",
        "coverage_count",
        "major_aa",
        "major_freq",
        "is_catalytic",
        "is_ligand_contact",
        "class",
    ]
    write_tsv(manifest_dir / "natural_scaffold_homolog_selection.tsv", selection_rows, columns_selection)
    write_tsv(manifest_dir / "natural_scaffold_conservation_by_position.tsv", conservation_rows, columns_pos)
    write_tsv(manifest_dir / "natural_scaffold_fixed_positions.tsv", fixed_rows, columns_pos)
    write_tsv(manifest_dir / "natural_scaffold_mutable_positions.tsv", mutable_rows, columns_pos)

    summary = {
        "method": "reference_anchored_pairwise_alignment_fallback",
        "full_msa_status": "NOT_RUN_TOOLCHAIN_UNAVAILABLE",
        "reference_accession": ref_accession,
        "reference_length": len(ref_seq),
        "homolog_query": args.homolog_query,
        "raw_homolog_records": len(raw_homologs),
        "selected_homolog_records": len(selected),
        "length_filter": [min_len, max_len],
        "min_identity": args.min_identity,
        "max_records": args.max_records,
        "catalytic_positions": sorted(CATALYTIC_POSITIONS),
        "ligand_contact_positions": sorted(ligand_contacts),
        "fixed_positions_count": len(fixed_rows),
        "mutable_positions_count": len(mutable_rows),
        "review_positions_count": sum(1 for r in conservation_rows if r["class"] == "REVIEW_CONSERVATION"),
        "outputs": {
            "homolog_selection": str(manifest_dir / "natural_scaffold_homolog_selection.tsv"),
            "conservation_by_position": str(manifest_dir / "natural_scaffold_conservation_by_position.tsv"),
            "fixed_positions": str(manifest_dir / "natural_scaffold_fixed_positions.tsv"),
            "mutable_positions": str(manifest_dir / "natural_scaffold_mutable_positions.tsv"),
        },
    }
    (manifest_dir / "natural_scaffold_msa_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
