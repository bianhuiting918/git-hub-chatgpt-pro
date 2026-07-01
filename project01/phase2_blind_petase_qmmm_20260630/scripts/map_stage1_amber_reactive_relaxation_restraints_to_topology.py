#!/usr/bin/env python3
"""Map Stage 1 reactive-relaxation restraints onto Amber topology indices.

Amber `tleap` can reorder or renumber atoms when it builds `prmtop`/`inpcrd`.
The pre-topology relaxation stage therefore cannot use its renumbered PDB
serials directly as Sander NMR-restraint `iat` values. This mapper rewrites
those `iat` values using the `complex_leap.pdb` serials, which match Amber atom
indices after topology preparation.
"""

from __future__ import annotations

import argparse
import csv
import math
import re
import shutil
from dataclasses import dataclass
from pathlib import Path


FIELDS = [
    "mapped_relaxation_job_id",
    "relaxation_job_id",
    "topology_prep_id",
    "pose_id",
    "mapped_iats",
    "mapped_relaxation_job_dir",
    "mapped_restraint_path",
    "mapped_mm_min_mdin_path",
    "runner_path",
    "prmtop_path",
    "inpcrd_path",
    "status",
]

IAT_RE = re.compile(r"iat\s*=\s*([0-9,\s-]+),")
DISANG_RE = re.compile(r"DISANG\s*=\s*\S+")


@dataclass(frozen=True)
class Atom:
    serial: int
    name: str
    resname: str
    chain: str
    resseq: int
    x: float
    y: float
    z: float
    element: str


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_pdb(path: Path) -> list[Atom]:
    atoms = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line[:6].strip() not in {"ATOM", "HETATM"}:
                continue
            try:
                atoms.append(
                    Atom(
                        serial=int(line[6:11]),
                        name=line[12:16].strip(),
                        resname=line[17:20].strip().upper(),
                        chain=line[21].strip(),
                        resseq=int(line[22:26]),
                        x=float(line[30:38]),
                        y=float(line[38:46]),
                        z=float(line[46:54]),
                        element=(line[76:78].strip().upper() if len(line) >= 78 else ""),
                    )
                )
            except ValueError:
                continue
    return atoms


def atom_by_serial(atoms: list[Atom]) -> dict[int, Atom]:
    return {atom.serial: atom for atom in atoms}


def element_or_name(atom: Atom) -> str:
    if atom.element:
        return atom.element
    for char in atom.name:
        if char.isalpha():
            return char.upper()
    return ""


HISTIDINE_NAMES = {"HIS", "HID", "HIE", "HIP"}


def canonical_resname(resname: str) -> str:
    return "HIS" if resname in HISTIDINE_NAMES else resname


def atom_family(atom: Atom, use_element: bool = False) -> tuple[str, str]:
    name = element_or_name(atom) if use_element else atom.name
    return canonical_resname(atom.resname), name


def occurrence_index(source_atom: Atom, atoms: list[Atom], use_element: bool = False) -> int:
    family = atom_family(source_atom, use_element)
    matches = [atom for atom in atoms if atom_family(atom, use_element) == family]
    for index, atom in enumerate(matches):
        if atom.serial == source_atom.serial:
            return index
    raise ValueError(f"Source atom serial {source_atom.serial} not found in its occurrence family")


def dist(a: Atom, b: Atom) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def map_atom(source_atom: Atom, source_atoms: list[Atom], topology_atoms: list[Atom], tolerance: float) -> Atom:
    family = atom_family(source_atom)
    occurrence = occurrence_index(source_atom, source_atoms)
    exact = [atom for atom in topology_atoms if atom_family(atom) == family]
    if occurrence < len(exact):
        return exact[occurrence]

    element_occurrence = occurrence_index(source_atom, source_atoms, use_element=True)
    element_family = atom_family(source_atom, use_element=True)
    element_matches = [atom for atom in topology_atoms if atom_family(atom, use_element=True) == element_family]
    if element_occurrence < len(element_matches):
        return element_matches[element_occurrence]

    candidates = exact or element_matches
    if not candidates:
        raise ValueError(f"No topology candidates for {source_atom.resname} {source_atom.name} serial {source_atom.serial}")
    mapped = min(candidates, key=lambda atom: dist(source_atom, atom))
    separation = dist(source_atom, mapped)
    if separation > tolerance:
        raise ValueError(
            f"Closest topology atom for {source_atom.resname} {source_atom.name} serial {source_atom.serial} "
            f"is {separation:.3f} A away, above tolerance {tolerance:.3f} A"
        )
    return mapped


def topology_rows_by_source(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    by_source = {}
    for row in rows:
        source = str(Path(row["source_structure_path"]).resolve())
        by_source[source] = row
    return by_source


def map_iat_text(
    text: str,
    source_by_serial: dict[int, Atom],
    source_atoms: list[Atom],
    topology_atoms: list[Atom],
    tolerance: float,
) -> tuple[str, list[str]]:
    mapped_iats = []

    def replace(match: re.Match[str]) -> str:
        source_serials = [int(value.strip()) for value in match.group(1).split(",") if value.strip()]
        mapped_serials = [
            map_atom(source_by_serial[source_serial], source_atoms, topology_atoms, tolerance).serial
            for source_serial in source_serials
        ]
        mapped_iats.append(",".join(str(serial) for serial in mapped_serials))
        return f"iat={','.join(str(serial) for serial in mapped_serials)},"

    return IAT_RE.sub(replace, text), mapped_iats


def mapped_mdin(text: str) -> str:
    mapped_disang = "DISANG=reactive_relaxation_restraints_mapped.RST"
    if DISANG_RE.search(text):
        return DISANG_RE.sub(mapped_disang, text, count=1)
    return text.rstrip() + "\n" + mapped_disang + "\n"


def runner(prmtop: Path, inpcrd: Path) -> str:
    return f"""#!/usr/bin/env bash
set -euo pipefail

SANDER=${{SANDER:-sander}}
PRMTOP=${{PRMTOP:-{prmtop}}}
INPCRD=${{INPCRD:-{inpcrd}}}

"${{SANDER}}" -O -i 00_restrained_mm_min_mapped.in -o 00_restrained_mm_min_mapped.out \\
  -p "${{PRMTOP}}" -c "${{INPCRD}}" -r 00_restrained_mm_min_mapped.rst7 -ref "${{INPCRD}}"
"""


def map_job(row: dict[str, str], topology_row: dict[str, str], out_dir: Path, tolerance: float) -> dict[str, str]:
    source_pdb = Path(row["complex_for_amber_pdb_path"])
    source_atoms = parse_pdb(source_pdb)
    source_by_serial = atom_by_serial(source_atoms)
    topology_atoms = parse_pdb(Path(topology_row["prep_dir"]) / "complex_leap.pdb")

    mapped_job_id = f"MAPPED_{row['relaxation_job_id']}"
    job_dir = out_dir / mapped_job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    source_rst = Path(row["restraint_path"])
    mapped_rst_text, mapped_iats = map_iat_text(
        source_rst.read_text(encoding="utf-8"),
        source_by_serial,
        source_atoms,
        topology_atoms,
        tolerance,
    )
    mapped_rst = job_dir / "reactive_relaxation_restraints_mapped.RST"
    mapped_rst.write_text(mapped_rst_text, encoding="utf-8")

    source_mdin = Path(row["mm_min_mdin_path"])
    mapped_mdin_path = job_dir / "00_restrained_mm_min_mapped.in"
    mapped_mdin_path.write_text(mapped_mdin(source_mdin.read_text(encoding="utf-8")), encoding="utf-8")

    copied_pdb = job_dir / "complex_for_amber.pdb"
    if source_pdb.exists():
        shutil.copyfile(source_pdb, copied_pdb)

    prmtop = Path(topology_row["prmtop_path"])
    inpcrd = Path(topology_row["inpcrd_path"])
    runner_path = job_dir / "run_reactive_mm_relaxation_mapped.sh"
    runner_path.write_text(runner(prmtop, inpcrd), encoding="utf-8")
    runner_path.chmod(0o755)

    return {
        "mapped_relaxation_job_id": mapped_job_id,
        "relaxation_job_id": row["relaxation_job_id"],
        "topology_prep_id": topology_row["topology_prep_id"],
        "pose_id": row.get("pose_id", ""),
        "mapped_iats": ";".join(mapped_iats),
        "mapped_relaxation_job_dir": str(job_dir),
        "mapped_restraint_path": str(mapped_rst),
        "mapped_mm_min_mdin_path": str(mapped_mdin_path),
        "runner_path": str(runner_path),
        "prmtop_path": str(prmtop),
        "inpcrd_path": str(inpcrd),
        "status": "reactive_relaxation_restraints_mapped_to_amber_topology",
    }


def generate(relaxation_manifest: Path, topology_prep_manifest: Path, out_dir: Path, tolerance: float) -> Path:
    topology_by_source = topology_rows_by_source(read_tsv(topology_prep_manifest))
    rows = []
    for row in read_tsv(relaxation_manifest):
        source = str(Path(row["complex_for_amber_pdb_path"]).resolve())
        if source in topology_by_source:
            rows.append(map_job(row, topology_by_source[source], out_dir, tolerance))
    manifest = out_dir / "amber_reactive_relaxation_topology_mapped_manifest.tsv"
    write_tsv(manifest, rows, FIELDS)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--relaxation-manifest", required=True)
    parser.add_argument("--topology-prep-manifest", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--coordinate-tolerance", type=float, default=0.08)
    args = parser.parse_args()
    manifest = generate(
        Path(args.relaxation_manifest),
        Path(args.topology_prep_manifest),
        Path(args.out_dir),
        args.coordinate_tolerance,
    )
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
