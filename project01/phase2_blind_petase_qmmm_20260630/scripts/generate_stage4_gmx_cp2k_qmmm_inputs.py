#!/usr/bin/env python3
"""Generate GMX-CP2K single-point QM/MM inputs for blind PETase Stage 4."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "qmmm_job_id",
    "conformer_id",
    "source_pose_id",
    "stage",
    "structure_path",
    "topology_path",
    "index_path",
    "qm_atom_group",
    "job_dir",
    "mdp_path",
    "cp2k_input_path",
    "runner_path",
    "extractor_path",
    "status",
    "energy_source",
    "source",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def productive_rows(rows: list[dict[str, str]], max_jobs: int | None) -> list[dict[str, str]]:
    selected = [
        row
        for row in rows
        if row.get("selection_status") == "productive"
        and row.get("representative_structure_path") not in {"", "pending", None}
    ]
    if max_jobs is not None:
        selected = selected[:max_jobs]
    return selected


def mdp(qm_group: str) -> str:
    return f"""integrator = steep
nsteps = 0
emtol = 1000.0
emstep = 0.01

cutoff-scheme = Verlet
coulombtype = Cut-off
rcoulomb = 1.2
rvdw = 1.2
pbc = xyz

qmmm = yes
qmmm-grps = {qm_group}
qmmm-qmgroup = {qm_group}
qmmm-qmmethod = cp2k
qmmm-qminput = topol-qmmm.inp
"""


def cp2k_input(conformer_id: str) -> str:
    return f"""&GLOBAL
  PROJECT PETASE_{conformer_id}_BLIND_QMMM_SP
  RUN_TYPE ENERGY
  PRINT_LEVEL LOW
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &DFT
    BASIS_SET_FILE_NAME BASIS_MOLOPT
    POTENTIAL_FILE_NAME GTH_POTENTIALS
    CHARGE 0
    MULTIPLICITY 1
    &MGRID
      CUTOFF 400
      REL_CUTOFF 50
    &END MGRID
    &SCF
      EPS_SCF 1.0E-6
      MAX_SCF 100
      SCF_GUESS ATOMIC
    &END SCF
    &XC
      &XC_FUNCTIONAL PBE
      &END XC_FUNCTIONAL
    &END XC
  &END DFT
&END FORCE_EVAL
"""


def readme(row: dict[str, str]) -> str:
    return "\n".join(
        [
            f"# GMX-CP2K QM/MM Job {row['conformer_id']}",
            "",
            "Purpose: blind Stage 4 QM/MM single-point or seed calculation for PETase mechanism exploration.",
            "",
            "Boundary: do not use article trajectories, article transition-state coordinates, article reaction coordinates, article selected CVs, or article barrier values.",
            "",
            "This generator follows the existing CPU-server PETase GMX-CP2K route:",
            "",
            "- use `/Dell/Dell14/bianht/gmx_cp2k_patched.sh` rather than PATH `gmx`;",
            "- run `grompp` with `-qmi topol-qmmm.inp` when the local build requires it;",
            "- run `mdrun -nsteps 1` for the single-point energy smoke calculation;",
            "- extract final comparison energy from GROMACS log `Potential Energy`, not CP2K `Total FORCE_EVAL`.",
            "",
            f"- Source pose: `{row.get('source_pose_id', '')}`",
            f"- Stage: `{row.get('stage', '')}`",
            f"- Structure: `{row.get('representative_structure_path', '')}`",
        ]
    ) + "\n"


def runner(wrapper: str) -> str:
    return f"""#!/usr/bin/env bash
set -euo pipefail

GMX_CP2K_WRAPPER=${{GMX_CP2K_WRAPPER:-{wrapper}}}
DEFFNM=${{DEFFNM:-qmmm_sp}}

"${{GMX_CP2K_WRAPPER}}" grompp \\
  -f qmmm_sp.mdp \\
  -c conf.gro \\
  -p topol.top \\
  -n index.ndx \\
  -qmi topol-qmmm.inp \\
  -o "${{DEFFNM}}.tpr" \\
  -maxwarn 8

OMP_NUM_THREADS=${{OMP_NUM_THREADS:-4}} "${{GMX_CP2K_WRAPPER}}" mdrun \\
  -s "${{DEFFNM}}.tpr" \\
  -deffnm "${{DEFFNM}}" \\
  -ntomp "${{OMP_NUM_THREADS}}" \\
  -nsteps 1
"""


def extractor() -> str:
    return """#!/usr/bin/env bash
set -euo pipefail

LOG=${1:-qmmm_sp.log}
grep 'Potential Energy' "$LOG" | tail -n 1
"""


def write_job(row: dict[str, str], out_dir: Path, wrapper: str) -> dict[str, str]:
    conformer_id = row["conformer_id"]
    job_id = f"QMMM_{conformer_id}_SP001"
    job_dir = out_dir / job_id
    job_dir.mkdir(parents=True, exist_ok=True)
    qm_group = row.get("qm_atom_group") or "QMatoms"

    files = {
        "00_README.md": readme(row),
        "qmmm_sp.mdp": mdp(qm_group),
        "topol-qmmm.inp": cp2k_input(conformer_id),
        "run_gmx_cp2k_qmmm.sh": runner(wrapper),
        "extract_potential_energy.sh": extractor(),
    }
    for name, content in files.items():
        path = job_dir / name
        path.write_text(content, encoding="utf-8")
        if name.endswith(".sh"):
            path.chmod(0o755)

    return {
        "qmmm_job_id": job_id,
        "conformer_id": conformer_id,
        "source_pose_id": row.get("source_pose_id", ""),
        "stage": row.get("stage", ""),
        "structure_path": row.get("representative_structure_path", ""),
        "topology_path": row.get("topology_path", "pending"),
        "index_path": row.get("index_path", "pending"),
        "qm_atom_group": qm_group,
        "job_dir": str(job_dir),
        "mdp_path": str(job_dir / "qmmm_sp.mdp"),
        "cp2k_input_path": str(job_dir / "topol-qmmm.inp"),
        "runner_path": str(job_dir / "run_gmx_cp2k_qmmm.sh"),
        "extractor_path": str(job_dir / "extract_potential_energy.sh"),
        "status": "inputs_ready_requires_grompp",
        "energy_source": "gromacs_log_potential_energy",
        "source": "blind_stage2_productive_conformer",
    }


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def generate(productive_conformers: Path, out_dir: Path, max_jobs: int | None, wrapper: str) -> Path:
    rows = productive_rows(read_tsv(productive_conformers), max_jobs)
    manifest_rows = [write_job(row, out_dir, wrapper) for row in rows]
    manifest = out_dir / "gmx_cp2k_qmmm_job_manifest.tsv"
    write_manifest(manifest, manifest_rows)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--productive-conformers", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--max-jobs", type=int)
    parser.add_argument("--wrapper", default="/Dell/Dell14/bianht/gmx_cp2k_patched.sh")
    args = parser.parse_args()
    manifest = generate(Path(args.productive_conformers), Path(args.out_dir), args.max_jobs, args.wrapper)
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
