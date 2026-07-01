#!/usr/bin/env python3
"""Prepare AmberTools topology-building inputs for Stage 4 Amber/Sander QM/MM jobs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "topology_prep_id",
    "qmmm_job_id",
    "pose_id",
    "ligand_model",
    "source_structure_path",
    "prep_dir",
    "complex_pdb_path",
    "ligand_pdb_path",
    "tleap_input_path",
    "runner_path",
    "prmtop_path",
    "inpcrd_path",
    "status",
    "source",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def selected_qmmm_rows(rows: list[dict[str, str]], max_jobs: int | None) -> list[dict[str, str]]:
    selected = []
    for row in rows:
        if row.get("engine") == "amber_sander" and row.get("source_structure_path"):
            selected.append(row)
    if max_jobs is not None:
        selected = selected[:max_jobs]
    return selected


def sanitize_id(raw: str) -> str:
    allowed = []
    for char in raw:
        if char.isalnum() or char in {"_", "-"}:
            allowed.append(char)
        else:
            allowed.append("_")
    return "".join(allowed)


def read_pdb_lines(path: Path) -> list[str]:
    return path.read_text(encoding="utf-8").splitlines()


def ligand_lines(lines: list[str]) -> list[str]:
    selected = []
    for line in lines:
        record = line[:6].strip()
        resname = line[17:20].strip().upper() if len(line) >= 20 else ""
        if record == "HETATM" and resname == "LIG":
            selected.append(line)
    if not selected:
        raise ValueError("No HETATM LIG records found for Amber ligand preparation")
    return selected


def write_pdb(path: Path, lines: list[str]) -> None:
    path.write_text("\n".join(lines + ["END"]) + "\n", encoding="utf-8")


def tleap_input(buffer_angstrom: float) -> str:
    return f"""source leaprc.protein.ff14SB
source leaprc.gaff2
source leaprc.water.tip3p

LIG = loadmol2 ligand.mol2
loadamberparams ligand.frcmod
complex = loadpdb complex_for_leap.pdb
check complex
solvatebox complex TIP3PBOX {buffer_angstrom:.1f}
addionsrand complex Na+ 0
addionsrand complex Cl- 0
saveamberparm complex complex.prmtop complex.inpcrd
savepdb complex complex_leap.pdb
quit
"""


def runner(ligand_charge: int) -> str:
    return f"""#!/usr/bin/env bash
set -euo pipefail

ANTECHAMBER=${{ANTECHAMBER:-antechamber}}
PARMCHK2=${{PARMCHK2:-parmchk2}}
TLEAP=${{TLEAP:-tleap}}
LIGAND_CHARGE=${{LIGAND_CHARGE:-{ligand_charge}}}

"${{ANTECHAMBER}}" \
  -i ligand.pdb -fi pdb \
  -o ligand.mol2 -fo mol2 \
  -at gaff2 -c bcc -s 2 \
  -rn LIG -nc "${{LIGAND_CHARGE}}"

"${{PARMCHK2}}" \
  -i ligand.mol2 -f mol2 \
  -o ligand.frcmod

"${{TLEAP}}" -f tleap.in > tleap.log

test -s complex.prmtop
test -s complex.inpcrd
"""


def readme(row: dict[str, str], ligand_charge: int, buffer_angstrom: float) -> str:
    return "\n".join(
        [
            f"# Amber Topology Preparation {row.get('qmmm_job_id', '')}",
            "",
            "Purpose: build Amber `prmtop` and `inpcrd` inputs needed before the Sander DFTB3/MM minimization and equilibration job can run.",
            "",
            "Method route: GAFF2/AM1-BCC ligand preparation with AmberTools, ff14SB protein parameters, TIP3P water, and a 15 A solvent buffer for the blind PETase QM/MM workflow.",
            "",
            "Boundary: this preparation uses our accepted seed structure and generic AmberTools parameterization only. Literature coordinates, trajectories, reaction coordinates, windows, barriers, and mechanism conclusions must not be imported before final validation.",
            "",
            f"- Source QM/MM job: `{row.get('qmmm_job_id', '')}`",
            f"- Source pose: `{row.get('pose_id', '')}`",
            f"- Ligand charge: `{ligand_charge}`",
            f"- Solvent buffer A: `{buffer_angstrom:.1f}`",
            "- Salt: neutralization is scripted now; 0.10 M NaCl ion-count audit is a required next gate after solvated water count is known.",
        ]
    ) + "\n"


def write_job(row: dict[str, str], out_dir: Path, ligand_charge: int, buffer_angstrom: float) -> dict[str, str]:
    qmmm_job_id = row["qmmm_job_id"]
    prep_id = f"TOPOLOGY_{sanitize_id(qmmm_job_id)}"
    prep_dir = out_dir / prep_id
    prep_dir.mkdir(parents=True, exist_ok=True)

    source_structure = Path(row["source_structure_path"])
    lines = read_pdb_lines(source_structure)
    complex_pdb = prep_dir / "complex_for_leap.pdb"
    ligand_pdb = prep_dir / "ligand.pdb"
    tleap_path = prep_dir / "tleap.in"
    runner_path = prep_dir / "run_amber_topology_prep.sh"

    write_pdb(complex_pdb, [line for line in lines if line[:6].strip() in {"ATOM", "HETATM", "TER"}])
    write_pdb(ligand_pdb, ligand_lines(lines))
    tleap_path.write_text(tleap_input(buffer_angstrom), encoding="utf-8")
    runner_path.write_text(runner(ligand_charge), encoding="utf-8")
    runner_path.chmod(0o755)
    (prep_dir / "00_README.md").write_text(readme(row, ligand_charge, buffer_angstrom), encoding="utf-8")

    return {
        "topology_prep_id": prep_id,
        "qmmm_job_id": qmmm_job_id,
        "pose_id": row.get("pose_id", ""),
        "ligand_model": row.get("ligand_model", ""),
        "source_structure_path": str(source_structure),
        "prep_dir": str(prep_dir),
        "complex_pdb_path": str(complex_pdb),
        "ligand_pdb_path": str(ligand_pdb),
        "tleap_input_path": str(tleap_path),
        "runner_path": str(runner_path),
        "prmtop_path": str(prep_dir / "complex.prmtop"),
        "inpcrd_path": str(prep_dir / "complex.inpcrd"),
        "status": "topology_prep_ready_requires_ambertools_execution",
        "source": "blind_stage4_amber_qmmm_manifest",
    }


def generate(
    amber_qmmm_manifest: Path,
    out_dir: Path,
    max_jobs: int | None,
    ligand_charge: int,
    buffer_angstrom: float,
) -> Path:
    rows = selected_qmmm_rows(read_tsv(amber_qmmm_manifest), max_jobs)
    manifest_rows = [write_job(row, out_dir, ligand_charge, buffer_angstrom) for row in rows]
    manifest = out_dir / "amber_topology_prep_manifest.tsv"
    write_tsv(manifest, manifest_rows, FIELDS)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amber-qmmm-manifest", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--max-jobs", type=int)
    parser.add_argument("--ligand-charge", type=int, default=0)
    parser.add_argument("--buffer-angstrom", type=float, default=15.0)
    args = parser.parse_args()
    manifest = generate(
        Path(args.amber_qmmm_manifest),
        Path(args.out_dir),
        args.max_jobs,
        args.ligand_charge,
        args.buffer_angstrom,
    )
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
