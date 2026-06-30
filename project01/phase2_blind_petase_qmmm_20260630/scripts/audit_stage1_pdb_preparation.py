from __future__ import annotations

import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path

TRIAD_BY_PDB = {
    "5XJH": ("A", 160, "SER", 237, "HIS", 206, "ASP"),
    "5YFE": ("A", 134, "SER", 211, "HIS", 180, "ASP"),
    "6EQD": ("A", 160, "SER", 237, "HIS", 206, "ASP"),
    "6EQE": ("A", 160, "SER", 237, "HIS", 206, "ASP"),
    "6EQF": ("A", 160, "SER", 237, "HIS", 206, "ASP"),
    "6EQG": ("A", 160, "SER", 237, "HIS", 206, "ASP"),
    "6EQH": ("A", 160, "SER", 237, "HIS", 206, "ASP"),
    "6ILW": ("A", 160, "SER", 237, "HIS", 206, "ASP"),
    "6QGC": ("A", 160, "SER", 237, "HIS", 206, "ASP"),
}

AA3 = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL",
    "SEC", "PYL", "MSE",
}


def distance(a, b) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def parse_atom_line(line: str):
    return {
        "record": line[0:6].strip(),
        "atom": line[12:16].strip(),
        "altloc": line[16].strip(),
        "resname": line[17:20].strip(),
        "chain": line[21].strip() or "_",
        "resseq": int(line[22:26]),
        "icode": line[26].strip(),
        "x": float(line[30:38]),
        "y": float(line[38:46]),
        "z": float(line[46:54]),
        "element": line[76:78].strip() if len(line) >= 78 else "",
    }


def parse_missing_residue(line: str):
    if len(line) < 27:
        return None
    resname = line[15:18].strip()
    chain = line[19].strip() or "_"
    seq_text = line[21:26].strip()
    if not resname or not seq_text or not re.match(r"^-?\d+$", seq_text):
        return None
    return f"{chain}:{resname}{seq_text}{line[26].strip()}"


def parse_missing_atom(line: str):
    if len(line) < 32:
        return None
    resname = line[15:18].strip()
    chain = line[19].strip() or "_"
    seq_text = line[20:26].strip()
    atoms = line[30:].strip()
    if not resname or not seq_text or not re.match(r"^-?\d+$", seq_text) or not atoms:
        return None
    return f"{chain}:{resname}{seq_text}{line[26].strip()}:{atoms}"


def residue_label(chain: str, resname: str, resseq: int, icode: str = "") -> str:
    return f"{chain}:{resname}{resseq}{icode}"


def load_pdb(path: Path):
    atoms = []
    ssbond_records = []
    missing_residues = []
    missing_atoms = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")):
                try:
                    atoms.append(parse_atom_line(line))
                except ValueError:
                    continue
            elif line.startswith("SSBOND"):
                ssbond_records.append(line.rstrip())
            elif line.startswith("REMARK 465"):
                item = parse_missing_residue(line)
                if item:
                    missing_residues.append(item)
            elif line.startswith("REMARK 470"):
                item = parse_missing_atom(line)
                if item:
                    missing_atoms.append(item)
    return atoms, ssbond_records, missing_residues, missing_atoms


def residue_key(atom):
    return (atom["chain"], atom["resseq"], atom["icode"], atom["resname"])


def point(atom):
    return (atom["x"], atom["y"], atom["z"])


def audit_pdb(path: Path):
    pdb = path.stem.upper()
    atoms, ssbond_records, missing_residues, missing_atoms = load_pdb(path)
    protein_residues_by_chain = defaultdict(set)
    residue_altlocs = defaultdict(set)
    nonwater_hets = defaultdict(set)
    waters = []
    cys_sg = []
    residues = defaultdict(dict)

    for atom in atoms:
        key = residue_key(atom)
        residues[key][atom["atom"]] = atom
        if atom["altloc"]:
            residue_altlocs[key].add(atom["altloc"])
        if atom["resname"] in AA3 and atom["record"] == "ATOM":
            protein_residues_by_chain[atom["chain"]].add((atom["resseq"], atom["icode"], atom["resname"]))
        if atom["record"] == "HETATM" and atom["resname"] not in {"HOH", "WAT", "DOD"}:
            nonwater_hets[atom["resname"]].add(atom["chain"])
        if atom["record"] == "HETATM" and atom["resname"] in {"HOH", "WAT", "DOD"} and atom["atom"] in {"O", "OW"}:
            waters.append(atom)
        if atom["resname"] == "CYS" and atom["atom"] == "SG" and atom["altloc"] in {"", "A"}:
            cys_sg.append(atom)

    triad = TRIAD_BY_PDB.get(pdb)
    triad_atoms = []
    triad_label = "not_assigned"
    triad_complete = "no"
    catalytic_waters = []

    if triad:
        chain, ser_i, ser_name, his_i, his_name, asp_i, asp_name = triad
        triad_label = f"{chain}:{ser_name}{ser_i}-{his_name}{his_i}-{asp_name}{asp_i}"
        triad_keys = [(chain, ser_i, "", ser_name), (chain, his_i, "", his_name), (chain, asp_i, "", asp_name)]
        wanted_atoms = {"SER": {"OG"}, "HIS": {"ND1", "NE2"}, "ASP": {"OD1", "OD2"}}
        for key in triad_keys:
            res = residues.get(key, {})
            for atom_name in wanted_atoms.get(key[3], set()):
                atom = res.get(atom_name)
                if atom:
                    triad_atoms.append(atom)
        triad_complete = "yes" if len(triad_atoms) >= 5 else "no"
        for wat in waters:
            best = min((distance(point(wat), point(atom)) for atom in triad_atoms), default=999.0)
            if best <= 4.0:
                catalytic_waters.append({
                    "pdb": pdb,
                    "water": residue_label(wat["chain"], wat["resname"], wat["resseq"], wat["icode"]),
                    "nearest_triad_distance_A": f"{best:.2f}",
                    "retain_candidate": "yes",
                    "reason": "within_4A_of_SerHisAsp_hetero_atoms",
                })

    cys_pairs = []
    for i, atom_a in enumerate(cys_sg):
        for atom_b in cys_sg[i + 1:]:
            if atom_a["chain"] != atom_b["chain"]:
                continue
            d = distance(point(atom_a), point(atom_b))
            if d <= 2.35:
                cys_pairs.append(
                    f"{residue_label(atom_a['chain'], 'CYS', atom_a['resseq'], atom_a['icode'])}-"
                    f"{residue_label(atom_b['chain'], 'CYS', atom_b['resseq'], atom_b['icode'])}:{d:.2f}A"
                )

    chains = sorted(protein_residues_by_chain)
    chain_counts = ";".join(f"{chain}:{len(protein_residues_by_chain[chain])}" for chain in chains)
    selected_chain = triad[0] if triad else (chains[0] if chains else "")
    selected_chain_reason = "contains_consistent_SerHisAsp_triad" if triad else "first_protein_chain_fallback"
    altloc_residues = [
        residue_label(chain, resname, resseq, icode) + ":" + ",".join(sorted(values))
        for (chain, resseq, icode, resname), values in sorted(residue_altlocs.items())
    ]
    audit_row = {
        "pdb": pdb,
        "selected_chain": selected_chain,
        "selected_chain_reason": selected_chain_reason,
        "protein_chains": ";".join(chains),
        "protein_residue_counts": chain_counts,
        "triad": triad_label,
        "triad_complete": triad_complete,
        "missing_residue_count": str(len(missing_residues)),
        "missing_residues": ";".join(missing_residues[:20]),
        "missing_atom_count": str(len(missing_atoms)),
        "missing_atoms": ";".join(missing_atoms[:20]),
        "altloc_residue_count": str(len(altloc_residues)),
        "altloc_residues": ";".join(altloc_residues[:20]),
        "nonwater_heterogens": ";".join(f"{name}:{','.join(sorted(chains))}" for name, chains in sorted(nonwater_hets.items())),
        "ssbond_record_count": str(len(ssbond_records)),
        "geometric_cys_sg_pairs": ";".join(cys_pairs),
        "water_count": str(len(waters)),
        "catalytic_water_candidate_count": str(len(catalytic_waters)),
        "prep_status": "audit_only_needs_manual_preparation",
    }
    chain_row = {
        "pdb": pdb,
        "selected_chain": selected_chain,
        "chains_available": ";".join(chains),
        "decision": "use_selected_chain_for_initial_blind_setup",
        "reason": selected_chain_reason,
        "sensitivity_check": "repeat_key_setup_on_secondary_templates",
    }
    return audit_row, chain_row, catalytic_waters


def write_tsv(path: Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    if len(sys.argv) != 3:
        print("usage: audit_stage1_pdb_preparation.py <pdb-dir> <out-dir>", file=sys.stderr)
        return 2
    pdb_dir = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])
    audit_rows = []
    chain_rows = []
    water_rows = []
    for path in sorted(pdb_dir.glob("*.pdb")):
        audit_row, chain_row, catalytic_waters = audit_pdb(path)
        audit_rows.append(audit_row)
        chain_rows.append(chain_row)
        water_rows.extend(catalytic_waters)

    write_tsv(out_dir / "pdb_preparation_audit.tsv", audit_rows, [
        "pdb", "selected_chain", "selected_chain_reason", "protein_chains", "protein_residue_counts",
        "triad", "triad_complete", "missing_residue_count", "missing_residues", "missing_atom_count",
        "missing_atoms", "altloc_residue_count", "altloc_residues", "nonwater_heterogens",
        "ssbond_record_count", "geometric_cys_sg_pairs", "water_count", "catalytic_water_candidate_count",
        "prep_status",
    ])
    write_tsv(out_dir / "template_chain_decisions.tsv", chain_rows, [
        "pdb", "selected_chain", "chains_available", "decision", "reason", "sensitivity_check",
    ])
    write_tsv(out_dir / "retained_water_candidates.tsv", water_rows, [
        "pdb", "water", "nearest_triad_distance_A", "retain_candidate", "reason",
    ])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
