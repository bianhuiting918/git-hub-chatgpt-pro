#!/usr/bin/env python3
"""Map Stage 4 QM atom selections onto Amber topology atom indices after tleap."""

from __future__ import annotations

import argparse
import csv
import math
import re
from dataclasses import dataclass
from pathlib import Path


FIELDS = [
    "mapped_job_id",
    "qmmm_job_id",
    "topology_prep_id",
    "pose_id",
    "mapped_atom_count",
    "mapped_qmmask",
    "mapped_job_dir",
    "smoke_mdin_path",
    "min_mdin_path",
    "equil_mdin_path",
    "runner_path",
    "prmtop_path",
    "inpcrd_path",
    "status",
]


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
    with path.open("r", encoding="utf-8") as handle:
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


def dist(a: Atom, b: Atom) -> float:
    return math.sqrt((a.x - b.x) ** 2 + (a.y - b.y) ** 2 + (a.z - b.z) ** 2)


def atom_by_serial(atoms: list[Atom]) -> dict[int, Atom]:
    return {atom.serial: atom for atom in atoms}


def element_or_name(atom: Atom) -> str:
    if atom.element:
        return atom.element
    for char in atom.name:
        if char.isalpha():
            return char.upper()
    return ""


def map_atom(source_atom: Atom, leap_atoms: list[Atom], tolerance: float) -> Atom:
    strict = [atom for atom in leap_atoms if atom.name == source_atom.name and atom.resname == source_atom.resname]
    candidates = strict or [
        atom
        for atom in leap_atoms
        if atom.resname == source_atom.resname and element_or_name(atom) == element_or_name(source_atom)
    ]
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


def replace_qmmask(text: str, qmmask: str) -> str:
    return re.sub(r"qmmask='@[^']+'", f"qmmask='{qmmask}'", text, count=1)


def smoke_from_min(min_text: str, qmmask: str) -> str:
    text = replace_qmmask(min_text, qmmask)
    text = re.sub(r"maxcyc\s*=\s*\d+", "maxcyc=1", text, count=1)
    text = re.sub(r"ncyc\s*=\s*\d+", "ncyc=1", text, count=1)
    if "maxcyc=1" not in text:
        text = text.replace("imin=1,", "imin=1, maxcyc=1, ncyc=1,", 1)
    return text


def runner(prmtop: Path, inpcrd: Path) -> str:
    return f"""#!/usr/bin/env bash
set -euo pipefail

SANDER=${{SANDER:-sander}}
PRMTOP=${{PRMTOP:-{prmtop}}}
INPCRD=${{INPCRD:-{inpcrd}}}

"${{SANDER}}" -O \
  -i 00_qmmm_smoke_min_1step.in \
  -o 00_qmmm_smoke_min_1step.out \
  -p "${{PRMTOP}}" \
  -c "${{INPCRD}}" \
  -r 00_qmmm_smoke_min_1step.rst7 \
  -ref "${{INPCRD}}"
"""


def topology_rows_by_job(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row["qmmm_job_id"]: row for row in rows}


def map_job(qrow: dict[str, str], trow: dict[str, str], out_dir: Path, tolerance: float) -> dict[str, str]:
    source_atoms = atom_by_serial(parse_pdb(Path(qrow["source_structure_path"])))
    leap_pdb = Path(trow["prep_dir"]) / "complex_leap.pdb"
    leap_atoms = parse_pdb(leap_pdb)
    selection_rows = read_tsv(Path(qrow["qm_selection_path"]))

    mapped_serials = []
    for selection in selection_rows:
        source_serial = int(selection["serial"])
        source_atom = source_atoms[source_serial]
        mapped_serials.append(map_atom(source_atom, leap_atoms, tolerance).serial)
    qmmask = "@" + ",".join(str(serial) for serial in mapped_serials)

    mapped_job_id = f"MAPPED_{qrow['qmmm_job_id']}"
    job_dir = out_dir / mapped_job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    min_text = replace_qmmask(Path(qrow["min_mdin_path"]).read_text(encoding="utf-8"), qmmask)
    equil_text = replace_qmmask(Path(qrow["equil_mdin_path"]).read_text(encoding="utf-8"), qmmask)
    smoke_text = smoke_from_min(min_text, qmmask)

    smoke_path = job_dir / "00_qmmm_smoke_min_1step.in"
    min_path = job_dir / "01_qmmm_min_mapped.in"
    equil_path = job_dir / "02_qmmm_equil_200ps_mapped.in"
    runner_path = job_dir / "run_sander_qmmm_smoke.sh"
    smoke_path.write_text(smoke_text, encoding="utf-8")
    min_path.write_text(min_text, encoding="utf-8")
    equil_path.write_text(equil_text, encoding="utf-8")
    prmtop = Path(trow["prmtop_path"])
    inpcrd = Path(trow["inpcrd_path"])
    runner_path.write_text(runner(prmtop, inpcrd), encoding="utf-8")
    runner_path.chmod(0o755)

    return {
        "mapped_job_id": mapped_job_id,
        "qmmm_job_id": qrow["qmmm_job_id"],
        "topology_prep_id": trow["topology_prep_id"],
        "pose_id": qrow.get("pose_id", ""),
        "mapped_atom_count": str(len(mapped_serials)),
        "mapped_qmmask": qmmask,
        "mapped_job_dir": str(job_dir),
        "smoke_mdin_path": str(smoke_path),
        "min_mdin_path": str(min_path),
        "equil_mdin_path": str(equil_path),
        "runner_path": str(runner_path),
        "prmtop_path": str(prmtop),
        "inpcrd_path": str(inpcrd),
        "status": "qmmm_inputs_mapped_to_amber_topology",
    }


def generate(amber_qmmm_manifest: Path, topology_prep_manifest: Path, out_dir: Path, tolerance: float) -> Path:
    qrows = read_tsv(amber_qmmm_manifest)
    trows = topology_rows_by_job(read_tsv(topology_prep_manifest))
    rows = []
    for qrow in qrows:
        if qrow["qmmm_job_id"] in trows:
            rows.append(map_job(qrow, trows[qrow["qmmm_job_id"]], out_dir, tolerance))
    manifest = out_dir / "amber_qmmm_topology_mapped_manifest.tsv"
    write_tsv(manifest, rows, FIELDS)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amber-qmmm-manifest", required=True)
    parser.add_argument("--topology-prep-manifest", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--coordinate-tolerance", type=float, default=0.08)
    args = parser.parse_args()
    manifest = generate(
        Path(args.amber_qmmm_manifest),
        Path(args.topology_prep_manifest),
        Path(args.out_dir),
        args.coordinate_tolerance,
    )
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
