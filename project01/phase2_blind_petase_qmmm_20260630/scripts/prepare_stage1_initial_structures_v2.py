from __future__ import annotations

import csv
import hashlib
import sys
from collections import defaultdict
from pathlib import Path

SELECTED_CHAINS = {
    "5XJH": "A",
    "5YFE": "A",
    "6EQD": "A",
    "6EQE": "A",
    "6EQF": "A",
    "6EQG": "A",
    "6EQH": "A",
    "6ILW": "A",
    "6QGC": "A",
}

DEFAULT_DISULFIDES = [("CYS203", "CYS239"), ("CYS273", "CYS289")]
DISULFIDES = {"5YFE": [("CYS177", "CYS213"), ("CYS247", "CYS263")]}


def read_retained_waters(path: Path) -> dict[str, set[tuple[str, int, str]]]:
    waters: dict[str, set[tuple[str, int, str]]] = defaultdict(set)
    with path.open("r", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            chain, water = row["water"].split(":", 1)
            resname = "".join(ch for ch in water if ch.isalpha())
            resseq = int("".join(ch for ch in water if ch.isdigit()))
            waters[row["pdb"].upper()].add((chain, resseq, resname))
    return waters


def atom_identity(line: str):
    return {
        "record": line[0:6],
        "atom": line[12:16],
        "altloc": line[16],
        "resname": line[17:20].strip(),
        "chain": line[21],
        "resseq": int(line[22:26]),
        "icode": line[26].strip(),
        "occupancy": float(line[54:60]) if line[54:60].strip() else 0.0,
    }


def normalize_atom_line(line: str, serial: int) -> str:
    return f"{line[:6]}{serial:5d}{line[11:16]} {line[17:]}"


def residue_altloc_choice(lines: list[str]) -> tuple[str, str]:
    occ_by_altloc: dict[str, float] = defaultdict(float)
    for line in lines:
        atom = atom_identity(line)
        altloc = atom["altloc"].strip()
        if altloc:
            occ_by_altloc[altloc] += atom["occupancy"]
    if not occ_by_altloc:
        return "blank", "no_altloc"
    ranked = sorted(occ_by_altloc.items(), key=lambda item: (-item[1], item[0]))
    return ranked[0][0], "highest_altloc_atom_occupancy"


def format_ssbond(serial: int, chain: str, pair: tuple[str, str]) -> str:
    cys1 = int(pair[0].replace("CYS", ""))
    cys2 = int(pair[1].replace("CYS", ""))
    return f"SSBOND {serial:3d} CYS {chain:1s} {cys1:4d}    CYS {chain:1s} {cys2:4d}                          2.05"


def prepare_one(pdb_path: Path, retained_waters: dict[str, set[tuple[str, int, str]]], out_dir: Path):
    pdb = pdb_path.stem.upper()
    chain = SELECTED_CHAINS[pdb]
    water_keys = retained_waters.get(pdb, set())
    protein_lines_by_residue: dict[tuple[str, int, str, str], list[str]] = defaultdict(list)
    kept_water_lines = []
    input_atom_count = 0

    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom = atom_identity(line)
            if atom["record"].startswith("ATOM") and atom["chain"] == chain:
                input_atom_count += 1
                key = (atom["chain"], atom["resseq"], atom["icode"], atom["resname"])
                protein_lines_by_residue[key].append(line.rstrip("\n"))
            elif atom["record"].startswith("HETATM") and (atom["chain"], atom["resseq"], atom["resname"]) in water_keys:
                kept_water_lines.append(line.rstrip("\n"))

    output_lines = [
        f"REMARK BLIND_PETASE_STAGE1_INITIAL_PREP {pdb}",
        "REMARK SOURCE RCSB_PDB_CHAIN_SELECTION_ONLY",
        "REMARK SCRIPT prepare_stage1_initial_structures_v2.py",
        "REMARK POLICY KEEP_SELECTED_CHAIN_A_KEEP_TRIAD_WATERS_DROP_OTHER_HETATM",
        "REMARK POLICY ALTLOC_HIGHEST_ALTLOC_ATOM_OCCUPANCY_KEEP_BLANK_ATOMS",
        "REMARK PROTONATION NOT_ASSIGNED",
    ]
    for serial, pair in enumerate(DISULFIDES.get(pdb, DEFAULT_DISULFIDES), start=1):
        output_lines.append(format_ssbond(serial, chain, pair))

    altloc_decisions = []
    serial = 1
    output_atom_count = 0
    for key in sorted(protein_lines_by_residue, key=lambda item: (item[0], item[1], item[2], item[3])):
        lines = protein_lines_by_residue[key]
        altlocs = sorted({atom_identity(line)["altloc"].strip() for line in lines if atom_identity(line)["altloc"].strip()})
        choice, reason = residue_altloc_choice(lines)
        if altlocs:
            altloc_decisions.append(
                {
                    "pdb": pdb,
                    "residue": f"{key[0]}:{key[3]}{key[1]}{key[2]}",
                    "available_altlocs": ",".join(altlocs),
                    "selected_altloc": choice,
                    "reason": reason,
                }
            )
        for line in lines:
            atom = atom_identity(line)
            altloc = atom["altloc"].strip()
            if not altloc or altloc == choice:
                output_lines.append(normalize_atom_line(line, serial))
                serial += 1
                output_atom_count += 1

    if kept_water_lines:
        output_lines.append("TER")
    for line in kept_water_lines:
        output_lines.append(normalize_atom_line(line, serial))
        serial += 1
        output_atom_count += 1
    output_lines.append("END")

    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{pdb}_chainA_initial_clean_v2.pdb"
    text = "\n".join(output_lines) + "\n"
    out_path.write_text(text, encoding="utf-8")
    sha256 = hashlib.sha256(text.encode("utf-8")).hexdigest()
    return {
        "pdb": pdb,
        "input_path": str(pdb_path),
        "prepared_path": str(out_path),
        "selected_chain": chain,
        "input_chain_atom_count": str(input_atom_count),
        "output_atom_count": str(output_atom_count),
        "retained_water_count": str(len(kept_water_lines)),
        "altloc_decision_count": str(len(altloc_decisions)),
        "protonation_status": "not_assigned",
        "heterogen_policy": "drop_nonwater_keep_listed_triad_waters",
        "disulfide_policy": "write_geometric_cys_pairs_from_audit",
        "sha256": sha256,
        "status": "initial_cleaned_v2_not_protonated",
    }, altloc_decisions


def write_tsv(path: Path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    if len(sys.argv) != 4:
        print("usage: prepare_stage1_initial_structures_v2.py <pdb-dir> <retained-water-tsv> <out-dir>", file=sys.stderr)
        return 2
    pdb_dir = Path(sys.argv[1])
    retained_water_tsv = Path(sys.argv[2])
    out_dir = Path(sys.argv[3])
    retained_waters = read_retained_waters(retained_water_tsv)
    rows = []
    altloc_rows = []
    for pdb_path in sorted(pdb_dir.glob("*.pdb")):
        pdb = pdb_path.stem.upper()
        if pdb in SELECTED_CHAINS:
            row, decisions = prepare_one(pdb_path, retained_waters, out_dir / "prepared_initial_pdbs")
            rows.append(row)
            altloc_rows.extend(decisions)
    write_tsv(
        out_dir / "prepared_structure_manifest.tsv",
        rows,
        [
            "pdb",
            "input_path",
            "prepared_path",
            "selected_chain",
            "input_chain_atom_count",
            "output_atom_count",
            "retained_water_count",
            "altloc_decision_count",
            "protonation_status",
            "heterogen_policy",
            "disulfide_policy",
            "sha256",
            "status",
        ],
    )
    write_tsv(
        out_dir / "altloc_resolution_decisions.tsv",
        altloc_rows,
        ["pdb", "residue", "available_altlocs", "selected_altloc", "reason"],
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
