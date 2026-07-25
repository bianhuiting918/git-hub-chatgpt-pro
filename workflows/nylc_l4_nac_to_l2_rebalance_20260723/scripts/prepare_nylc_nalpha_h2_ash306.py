#!/usr/bin/env python3
"""Build an audited NylC chain-H NalphaH2/Asp306H topology proxy."""
from __future__ import annotations

import argparse
import json
import pathlib
import re
from typing import Dict, Tuple

INDEX_FIELD_COUNTS = {
    "bonds": 2, "pairs": 2, "pairs_nb": 2, "angles": 3,
    "dihedrals": 4, "constraints": 2, "cmap": 5,
    "position_restraints": 1, "distance_restraints": 2,
    "orientation_restraints": 2, "angle_restraints": 4,
    "dihedral_restraints": 4,
}


def _section_name(line: str) -> str | None:
    match = re.match(r"\s*\[\s*([^]]+?)\s*\]", line)
    return match.group(1).strip().lower() if match else None


def _split_comment(line: str) -> Tuple[str, str]:
    if ";" not in line:
        return line, ""
    body, comment = line.split(";", 1)
    return body, ";" + comment


def rewrite_chain_itp(text: str, remove_index1: int = 4, nalpha_charge: float = -0.62540):
    """Remove N-terminal H3, neutralize Nalpha, and renumber topology references."""
    section = ""
    output = []
    old_atom_count = new_atom_count = 0
    old_charge = new_charge = 0.0
    removed_atom = None
    dropped: Dict[str, int] = {"bonds": 0, "angles": 0, "dihedrals": 0}

    for original in text.splitlines():
        detected = _section_name(original)
        if detected is not None:
            section = detected
            output.append(original)
            continue
        stripped = original.strip()
        if not stripped or stripped.startswith((";", "#")):
            if section == "atoms" and "rtp NTHR q +1.0" in original:
                original = original.replace("rtp NTHR q +1.0", "rtp NTHR_NH2_PROXY q 0.0")
            output.append(original)
            continue

        body, comment = _split_comment(original)
        tokens = body.split()
        if section == "atoms":
            if len(tokens) < 8 or not tokens[0].isdigit():
                output.append(original)
                continue
            old_atom_count += 1
            atom_index = int(tokens[0])
            atom_name = tokens[4]
            charge = float(tokens[6])
            old_charge += charge
            if atom_index == remove_index1:
                if atom_name != "H3":
                    raise ValueError(f"atom {remove_index1} must be H3, found {atom_name}")
                removed_atom = {"index1": atom_index, "name": atom_name, "charge_e": charge}
                continue
            if atom_index == 1:
                if atom_name != "N":
                    raise ValueError("atom 1 must be Nalpha N")
                tokens[6] = f"{nalpha_charge:.5f}"
                charge = nalpha_charge
            if atom_index > remove_index1:
                tokens[0] = str(atom_index - 1)
            if tokens[5].isdigit() and int(tokens[5]) > remove_index1:
                tokens[5] = str(int(tokens[5]) - 1)
            new_atom_count += 1
            new_charge += charge
            if "qtot" in comment:
                match = re.search(r"qtot\s+([+-]?[0-9]+(?:\.[0-9]+)?)", comment)
                if match:
                    corrected = float(match.group(1)) - 1.0
                    comment = comment[:match.start(1)] + f"{corrected:g}" + comment[match.end(1):]
            output.append(" ".join(tokens) + ((" " + comment) if comment else ""))
            continue

        index_count = INDEX_FIELD_COUNTS.get(section)
        if section == "exclusions":
            index_count = len(tokens)
        if index_count is None or len(tokens) < index_count:
            output.append(original)
            continue
        try:
            indices = [int(token) for token in tokens[:index_count]]
        except ValueError:
            output.append(original)
            continue
        if remove_index1 in indices:
            if section in dropped:
                dropped[section] += 1
            continue
        for index, value in enumerate(indices):
            if value > remove_index1:
                tokens[index] = str(value - 1)
        output.append(" ".join(tokens) + ((" " + comment) if comment else ""))

    if removed_atom is None:
        raise ValueError(f"H3 atom {remove_index1} was not found")
    if new_atom_count != old_atom_count - 1:
        raise ValueError("atom-count delta is not exactly one")
    charge_delta = new_charge - old_charge
    if abs(charge_delta + 1.0) > 1e-8:
        raise ValueError(f"expected topology charge delta -1, observed {charge_delta:.12f}")
    audit = {
        "old_atom_count": old_atom_count, "new_atom_count": new_atom_count,
        "old_charge_e": round(old_charge, 8), "new_charge_e": round(new_charge, 8),
        "charge_delta_e": round(charge_delta, 8), "nalpha_charge_e": nalpha_charge,
        "removed_atom": removed_atom, "dropped_interaction_counts": dropped,
    }
    return "\n".join(output) + "\n", audit


def drop_gro_atom(text: str, remove_index1: int = 4, expected_name: str = "H3"):
    lines = text.splitlines()
    if len(lines) < 4:
        raise ValueError("invalid GRO")
    count = int(lines[1])
    atom_lines = lines[2:2 + count]
    if len(atom_lines) != count or len(lines) != count + 3:
        raise ValueError("GRO atom-count mismatch")
    selected = atom_lines[remove_index1 - 1]
    found_name = selected[10:15].strip()
    if found_name != expected_name:
        raise ValueError(f"GRO atom {remove_index1} must be {expected_name}, found {found_name}")
    kept = atom_lines[:remove_index1 - 1] + atom_lines[remove_index1:]
    renumbered = [line[:15] + f"{index % 100000:5d}" + line[20:] for index, line in enumerate(kept, 1)]
    result = [lines[0], str(count - 1), *renumbered, lines[-1]]
    audit = {"old_atom_count": count, "new_atom_count": count - 1,
             "removed_atom": {"index1": remove_index1, "name": found_name},
             "box_line": lines[-1]}
    return "\n".join(result) + "\n", audit



def extract_molecule_itp(standalone_topology: str) -> str:
    """Extract one protein moleculetype block from a standalone pdb2gmx topology."""
    start = re.search(r"(?m)^[ \t]*\[\s*moleculetype\s*\][ \t]*$", standalone_topology)
    if not start:
        raise ValueError("standalone topology has no moleculetype block")
    tail = standalone_topology[start.start():]
    stop = re.search(r"(?m)^;\s*Include water topology\s*$", tail)
    if not stop:
        stop = re.search(r"(?m)^\s*\[\s*system\s*\]\s*$", tail)
    if not stop:
        raise ValueError("cannot locate end of protein moleculetype block")
    result = tail[:stop.start()].rstrip() + "\n"
    if "[ atoms ]" not in result or "Protein_chain_H" not in result:
        raise ValueError("extracted block is not Protein_chain_H")
    return result


def replace_gro_slice(full_text: str, chain_text: str, start_index1: int, end_index1: int):
    """Replace an equal-length contiguous atom slice and renumber the full GRO."""
    full_lines = full_text.splitlines()
    chain_lines = chain_text.splitlines()
    full_count = int(full_lines[1])
    chain_count = int(chain_lines[1])
    if len(full_lines) != full_count + 3 or len(chain_lines) != chain_count + 3:
        raise ValueError("GRO atom-count mismatch")
    expected = end_index1 - start_index1 + 1
    if chain_count != expected:
        raise ValueError(f"replacement chain has {chain_count} atoms, expected {expected}")
    if start_index1 < 1 or end_index1 > full_count:
        raise ValueError("replacement slice outside full system")
    atoms = full_lines[2:2 + full_count]
    replacement = chain_lines[2:2 + chain_count]
    merged = atoms[:start_index1 - 1] + replacement + atoms[end_index1:]
    if len(merged) != full_count:
        raise ValueError("full-system atom count changed")
    renumbered = [line[:15] + f"{index % 100000:5d}" + line[20:] for index, line in enumerate(merged, 1)]
    result = [full_lines[0] + " | NalphaH2_Asp306H", str(full_count), *renumbered, full_lines[-1]]
    audit = {"full_atom_count": full_count, "replacement_start_index1": start_index1,
             "replacement_end_index1": end_index1, "replaced_atom_count": chain_count,
             "box_line": full_lines[-1]}
    return "\n".join(result) + "\n", audit


def export_chain_pdb(full_text: str, start_index1: int, end_index1: int, ash_resid: int = 306, chain_id: str = "H"):
    """Export one GRO atom slice as PDB, renaming the target ASP to ASH."""
    lines = full_text.splitlines()
    count = int(lines[1])
    if len(lines) != count + 3:
        raise ValueError("GRO atom-count mismatch")
    if start_index1 < 1 or end_index1 > count or start_index1 > end_index1:
        raise ValueError("chain slice outside full system")
    selected = lines[2 + start_index1 - 1:2 + end_index1]
    pdb_lines = []
    renamed = 0
    heavy = 0
    for serial, line in enumerate(selected, 1):
        resid = int(line[:5])
        resname = line[5:10].strip()
        atomname = line[10:15].strip()
        if resid == ash_resid and resname == "ASP":
            resname = "ASH"
            renamed += 1
        letters = "".join(character for character in atomname if character.isalpha()).upper()
        if not letters:
            raise ValueError(f"cannot infer element for {atomname}")
        element = letters[0]
        if element != "H":
            heavy += 1
        x, y, z = (float(line[a:b]) * 10.0 for a, b in ((20, 28), (28, 36), (36, 44)))
        pdb_lines.append(
            f"ATOM  {serial:5d} {atomname:>4s} {resname:>3s} {chain_id}{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2s}"
        )
    if renamed == 0:
        raise ValueError(f"ASP residue {ash_resid} was not found")
    audit = {"source_start_index1": start_index1, "source_end_index1": end_index1,
             "chain_atom_count": len(selected), "chain_heavy_atom_count": heavy,
             "ash_resid": ash_resid, "renamed_ash306_atom_count": renamed}
    return "\n".join(pdb_lines) + "\nTER\nEND\n", audit

def main() -> int:
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest="command", required=True)
    itp = sub.add_parser("rewrite-itp")
    itp.add_argument("--input", type=pathlib.Path, required=True)
    itp.add_argument("--output", type=pathlib.Path, required=True)
    itp.add_argument("--audit", type=pathlib.Path, required=True)
    gro = sub.add_parser("drop-gro-h3")
    gro.add_argument("--input", type=pathlib.Path, required=True)
    gro.add_argument("--output", type=pathlib.Path, required=True)
    gro.add_argument("--audit", type=pathlib.Path, required=True)
    args = parser.parse_args()
    if args.command == "rewrite-itp":
        transformed, audit = rewrite_chain_itp(args.input.read_text())
    else:
        transformed, audit = drop_gro_atom(args.input.read_text())
    args.output.write_text(transformed)
    args.audit.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())