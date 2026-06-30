#!/usr/bin/env python3
"""Generate Amber/Sander DFTB3 QM/MM inputs for blind PETase Stage 4."""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path


DEFAULT_QM_RESIDUES = ["A:TRP159", "A:SER160", "A:ASP206", "A:HIS237", "A:SER238"]
BACKBONE_NAMES = {"N", "CA", "C", "O", "OXT", "H", "H1", "H2", "H3", "HA", "HA2", "HA3"}
FIELDS = [
    "qmmm_job_id",
    "pose_id",
    "ligand_model",
    "source_structure_path",
    "engine",
    "qm_theory",
    "qmcharge",
    "spin",
    "qm_atom_count",
    "qmmask",
    "job_dir",
    "qm_selection_path",
    "min_mdin_path",
    "equil_mdin_path",
    "runner_path",
    "status",
    "source",
]


@dataclass(frozen=True)
class Atom:
    serial: int
    name: str
    resname: str
    chain: str
    resseq: int
    record: str
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


def accepted_seed_rows(rows: list[dict[str, str]], max_jobs: int | None) -> list[dict[str, str]]:
    selected = []
    for row in rows:
        structure = row.get("relaxed_structure_path", "")
        if row.get("pass_fail") == "pass" and structure and structure != "pending":
            selected.append(row)
    if max_jobs is not None:
        selected = selected[:max_jobs]
    return selected


def parse_residue_selector(raw: str) -> tuple[str, str, int]:
    chain = ""
    residue = raw
    if ":" in raw:
        chain, residue = raw.split(":", 1)
    resname = residue[:3].upper()
    resseq = int(residue[3:])
    return chain, resname, resseq


def parse_pdb(path: Path) -> list[Atom]:
    atoms = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            record = line[:6].strip()
            if record not in {"ATOM", "HETATM"}:
                continue
            try:
                serial = int(line[6:11])
                name = line[12:16].strip()
                resname = line[17:20].strip().upper()
                chain = line[21].strip()
                resseq = int(line[22:26])
            except ValueError:
                continue
            element = line[76:78].strip() if len(line) >= 78 else ""
            atoms.append(Atom(serial, name, resname, chain, resseq, record, element))
    return atoms


def select_qm_atoms(atoms: list[Atom], selectors: list[str]) -> list[tuple[Atom, str]]:
    residue_selectors = {parse_residue_selector(raw) for raw in selectors}
    selected: list[tuple[Atom, str]] = []
    seen: set[int] = set()
    for atom in atoms:
        residue_key = (atom.chain, atom.resname, atom.resseq)
        chainless_key = ("", atom.resname, atom.resseq)
        if atom.record == "HETATM" and atom.resname == "LIG":
            role = "substrate_ligand"
        elif residue_key in residue_selectors or chainless_key in residue_selectors:
            if atom.name in BACKBONE_NAMES:
                continue
            role = "catalytic_or_neighbor_sidechain"
        else:
            continue
        if atom.serial not in seen:
            selected.append((atom, role))
            seen.add(atom.serial)
    return selected


def amber_mask(serials: list[int]) -> str:
    return "@" + ",".join(str(serial) for serial in serials)


def qmmm_block(qmmask: str, qmcharge: int, spin: int, qm_theory: str) -> str:
    return f"""&qmmm
  qmmask='{qmmask}',
  qmcharge={qmcharge},
  spin={spin},
  qm_theory='{qm_theory}',
  qmshake=0,
/
"""


def min_mdin(qmmask: str, qmcharge: int, spin: int, qm_theory: str) -> str:
    return f"""Blind PETase Amber/Sander DFTB3 QM/MM minimization
&cntrl
  imin=1, maxcyc=2000, ncyc=1000,
  ntb=1, cut=10.0,
  ntpr=50,
  ifqnt=1,
/
{qmmm_block(qmmask, qmcharge, spin, qm_theory)}"""


def equil_mdin(qmmask: str, qmcharge: int, spin: int, qm_theory: str) -> str:
    return f"""Blind PETase Amber/Sander DFTB3 QM/MM 200 ps equilibration
&cntrl
  imin=0, irest=0, ntx=1,
  nstlim=200000, dt=0.001,
  ntc=1, ntf=1,
  tempi=310.0, temp0=310.0,
  ntt=3, gamma_ln=1.0,
  ntb=1, cut=10.0,
  ntpr=500, ntwx=5000, ntwr=5000,
  ifqnt=1,
/
{qmmm_block(qmmask, qmcharge, spin, qm_theory)}"""


def runner() -> str:
    return """#!/usr/bin/env bash
set -euo pipefail

PRMTOP=${PRMTOP:?set PRMTOP to the Amber topology file}
INPCRD=${INPCRD:?set INPCRD to the Amber coordinate file}
SANDER=${SANDER:-sander}

"${SANDER}" -O \
  -i 01_qmmm_min.in \
  -o 01_qmmm_min.out \
  -p "${PRMTOP}" \
  -c "${INPCRD}" \
  -r 01_qmmm_min.rst7 \
  -ref "${INPCRD}"

"${SANDER}" -O \
  -i 02_qmmm_equil_200ps.in \
  -o 02_qmmm_equil_200ps.out \
  -p "${PRMTOP}" \
  -c 01_qmmm_min.rst7 \
  -r 02_qmmm_equil_200ps.rst7 \
  -x 02_qmmm_equil_200ps.nc \
  -ref 01_qmmm_min.rst7
"""


def readme(row: dict[str, str], qm_atom_count: int) -> str:
    return "\n".join(
        [
            f"# Amber/Sander QM/MM Job {row.get('pose_id', '')}",
            "",
            "Purpose: prepare the blind Stage 4 PETase QM/MM minimization and equilibration inputs from our accepted seed structures.",
            "",
            "Method route: Amber/Sander `ifqnt=1` with low-cost DFTB3 QM/MM, followed by 200 ps QM/MM equilibration before TS-search tooling is allowed to consume coordinates.",
            "",
            "Boundary: this input set is generated from structure, catalytic chemistry, and our own Stage 1/2 gates only. Literature coordinates, trajectories, CV formulas, windows, and numeric outputs must not be imported before final validation.",
            "",
            "Current status: Amber `prmtop` and `inpcrd` still need to be built or mapped for the selected seed before the Sander runner can execute.",
            "",
            f"- Source pose: `{row.get('pose_id', '')}`",
            f"- Ligand model: `{row.get('ligand_model', '')}`",
            f"- QM atom count: `{qm_atom_count}`",
        ]
    ) + "\n"


def resolve_structure_path(raw: str, manifest_path: Path) -> Path:
    path = Path(raw)
    if path.is_absolute():
        return path
    if path.exists():
        return path
    for base in [manifest_path.parent, *manifest_path.parents]:
        candidate = base / path
        if candidate.exists():
            return candidate
    return manifest_path.parent / path


def write_job(
    row: dict[str, str],
    structure_path: Path,
    out_dir: Path,
    index: int,
    selectors: list[str],
    qmcharge: int,
    spin: int,
    qm_theory: str,
) -> dict[str, str]:
    pose_id = row["pose_id"]
    job_id = f"AMBER_QMMM_{pose_id}_{index:03d}"
    job_dir = out_dir / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    selected = select_qm_atoms(parse_pdb(structure_path), selectors)
    if not selected:
        raise ValueError(f"No QM atoms selected from {structure_path}")
    serials = [atom.serial for atom, _role in selected]
    qmmask = amber_mask(serials)

    selection_rows = [
        {
            "serial": str(atom.serial),
            "atom_name": atom.name,
            "resname": atom.resname,
            "chain": atom.chain,
            "resseq": str(atom.resseq),
            "record": atom.record,
            "role": role,
        }
        for atom, role in selected
    ]
    selection_path = job_dir / "qm_atom_selection.tsv"
    write_tsv(selection_path, selection_rows, ["serial", "atom_name", "resname", "resseq", "chain", "record", "role"])

    files = {
        "00_README.md": readme(row, len(selected)),
        "01_qmmm_min.in": min_mdin(qmmask, qmcharge, spin, qm_theory),
        "02_qmmm_equil_200ps.in": equil_mdin(qmmask, qmcharge, spin, qm_theory),
        "run_sander_qmmm_smoke.sh": runner(),
    }
    for name, content in files.items():
        path = job_dir / name
        path.write_text(content, encoding="utf-8")
        if name.endswith(".sh"):
            path.chmod(0o755)

    return {
        "qmmm_job_id": job_id,
        "pose_id": pose_id,
        "ligand_model": row.get("ligand_model", ""),
        "source_structure_path": str(structure_path),
        "engine": "amber_sander",
        "qm_theory": qm_theory,
        "qmcharge": str(qmcharge),
        "spin": str(spin),
        "qm_atom_count": str(len(selected)),
        "qmmask": qmmask,
        "job_dir": str(job_dir),
        "qm_selection_path": str(selection_path),
        "min_mdin_path": str(job_dir / "01_qmmm_min.in"),
        "equil_mdin_path": str(job_dir / "02_qmmm_equil_200ps.in"),
        "runner_path": str(job_dir / "run_sander_qmmm_smoke.sh"),
        "status": "inputs_ready_needs_amber_prmtop_inpcrd",
        "source": "blind_stage1_accepted_seed",
    }


def generate(
    gs_pose_manifest: Path,
    out_dir: Path,
    max_jobs: int | None,
    selectors: list[str],
    qmcharge: int,
    spin: int,
    qm_theory: str,
) -> Path:
    rows = accepted_seed_rows(read_tsv(gs_pose_manifest), max_jobs)
    manifest_rows = []
    for index, row in enumerate(rows, start=1):
        structure_path = resolve_structure_path(row["relaxed_structure_path"], gs_pose_manifest)
        manifest_rows.append(write_job(row, structure_path, out_dir, index, selectors, qmcharge, spin, qm_theory))
    manifest = out_dir / "amber_qmmm_job_manifest.tsv"
    write_tsv(manifest, manifest_rows, FIELDS)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gs-pose-manifest", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--max-jobs", type=int)
    parser.add_argument("--qm-residue", action="append", default=[])
    parser.add_argument("--qmcharge", type=int, default=0)
    parser.add_argument("--spin", type=int, default=1)
    parser.add_argument("--qm-theory", default="DFTB3")
    args = parser.parse_args()
    selectors = args.qm_residue or DEFAULT_QM_RESIDUES
    manifest = generate(
        Path(args.gs_pose_manifest),
        Path(args.out_dir),
        args.max_jobs,
        selectors,
        args.qmcharge,
        args.spin,
        args.qm_theory,
    )
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


