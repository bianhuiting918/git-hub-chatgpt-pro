#!/usr/bin/env python3
"""Generate Amber MD input directories for blind PETase Stage 2 jobs."""

from __future__ import annotations

import argparse
import csv
import shlex
from pathlib import Path


FIELDS = [
    "md_job_id",
    "source_pose_id",
    "stage",
    "template_pdb",
    "substrate_model",
    "structure_path",
    "job_dir",
    "minimize_mdin",
    "heat_mdin",
    "density_mdin",
    "equilibrate_mdin",
    "production_mdin",
    "runner",
    "status",
    "source",
]


def read_queue(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def runnable_rows(rows: list[dict[str, str]], max_jobs: int | None) -> list[dict[str, str]]:
    selected = [
        row
        for row in rows
        if row.get("status") == "not_started"
        and row.get("source") == "blind_accepted_gs_pose"
        and row.get("structure_path") not in {"", "pending"}
    ]
    if max_jobs is not None:
        selected = selected[:max_jobs]
    return selected


def mdin_minimize() -> str:
    return """Stage 2 blind minimization
 &cntrl
  imin=1, maxcyc=5000, ncyc=2500,
  cut=10.0, ntb=1,
  ntpr=100,
 /
"""


def mdin_heat(temperature: str) -> str:
    return f"""Stage 2 blind NVT heating
 &cntrl
  imin=0, irest=0, ntx=1,
  nstlim=50000, dt=0.002,
  tempi=0.0, temp0={temperature},
  ntc=2, ntf=2,
  ntb=1, ntt=3, gamma_ln=2.0,
  ntpr=500, ntwx=500, ntwr=5000,
  cut=10.0,
 /
"""


def mdin_density(temperature: str, pressure: str) -> str:
    return f"""Stage 2 blind NPT density equilibration
 &cntrl
  imin=0, irest=1, ntx=5,
  nstlim=100000, dt=0.002,
  temp0={temperature},
  ntc=2, ntf=2,
  ntb=2, ntp=1, pres0={pressure},
  ntt=3, gamma_ln=2.0,
  ntpr=500, ntwx=500, ntwr=5000,
  cut=10.0,
 /
"""


def mdin_equilibrate(temperature: str, pressure: str) -> str:
    return f"""Stage 2 blind NPT equilibration
 &cntrl
  imin=0, irest=1, ntx=5,
  nstlim=250000, dt=0.002,
  temp0={temperature},
  ntc=2, ntf=2,
  ntb=2, ntp=1, pres0={pressure},
  ntt=3, gamma_ln=2.0,
  ntpr=1000, ntwx=1000, ntwr=10000,
  cut=10.0,
 /
"""


def mdin_production(temperature: str, pressure: str) -> str:
    return f"""Stage 2 blind short production replicate for active-site clustering
 &cntrl
  imin=0, irest=1, ntx=5,
  nstlim=250000, dt=0.002,
  temp0={temperature},
  ntc=2, ntf=2,
  ntb=2, ntp=1, pres0={pressure},
  ntt=3, gamma_ln=2.0,
  ntpr=1000, ntwx=1000, ntwr=10000,
  cut=10.0,
 /
"""


def readme(row: dict[str, str]) -> str:
    return "\n".join(
        [
            f"# Amber MD Job {row['md_job_id']}",
            "",
            "Purpose: run a blind Stage 2 PETase classical MD replicate before QM/MM exploration.",
            "",
            "Boundary: do not use article trajectories, article transition-state coordinates, article reaction coordinates, or article barrier values.",
            "",
            "Required before running:",
            "",
            "1. Build and inspect Amber topology/coordinates for the listed protein-ligand complex.",
            "2. Parameterize ligand atoms from the blind ligand model, not from paper structures.",
            "3. Solvate, neutralize, and manually inspect the active site.",
            "4. Run minimization, heating, density equilibration, equilibration, and short production.",
            "5. Cluster frames by blind catalytic geometry and update productive_conformer_manifest.tsv.",
            "",
            f"- Source pose: `{row['source_pose_id']}`",
            f"- Stage: `{row['stage']}`",
            f"- Structure path: `{row['structure_path']}`",
        ]
    ) + "\n"


def runner() -> str:
    return """#!/usr/bin/env bash
set -euo pipefail

ENGINE=${AMBER_ENGINE:-sander}
if command -v pmemd.cuda >/dev/null 2>&1 && [[ "${USE_CUDA:-0}" == "1" ]]; then
  ENGINE=pmemd.cuda
fi

PRMTOP=${PRMTOP:-complex.prmtop}
INPCRD=${INPCRD:-complex.inpcrd}

${ENGINE} -O -i 01_minimize.in -o 01_minimize.out -p ${PRMTOP} -c ${INPCRD} -r 01_minimize.rst7 -ref ${INPCRD}
${ENGINE} -O -i 02_heat.in -o 02_heat.out -p ${PRMTOP} -c 01_minimize.rst7 -r 02_heat.rst7 -x 02_heat.nc -ref 01_minimize.rst7
${ENGINE} -O -i 03_density.in -o 03_density.out -p ${PRMTOP} -c 02_heat.rst7 -r 03_density.rst7 -x 03_density.nc -ref 02_heat.rst7
${ENGINE} -O -i 04_equilibrate.in -o 04_equilibrate.out -p ${PRMTOP} -c 03_density.rst7 -r 04_equilibrate.rst7 -x 04_equilibrate.nc
${ENGINE} -O -i 05_production.in -o 05_production.out -p ${PRMTOP} -c 04_equilibrate.rst7 -r 05_production.rst7 -x 05_production.nc
"""


def write_job(row: dict[str, str], out_dir: Path) -> dict[str, str]:
    job_dir = out_dir / row["md_job_id"]
    job_dir.mkdir(parents=True, exist_ok=True)
    temperature = row.get("temperature_K") or "300"
    pressure = row.get("pressure_bar") or "1"

    files = {
        "00_README.md": readme(row),
        "01_minimize.in": mdin_minimize(),
        "02_heat.in": mdin_heat(temperature),
        "03_density.in": mdin_density(temperature, pressure),
        "04_equilibrate.in": mdin_equilibrate(temperature, pressure),
        "05_production.in": mdin_production(temperature, pressure),
        "run_amber_md.sh": runner(),
    }
    for name, content in files.items():
        path = job_dir / name
        path.write_text(content, encoding="utf-8")
        if name.endswith(".sh"):
            path.chmod(0o755)

    return {
        "md_job_id": row["md_job_id"],
        "source_pose_id": row["source_pose_id"],
        "stage": row["stage"],
        "template_pdb": row.get("template_pdb", ""),
        "substrate_model": row.get("substrate_model", ""),
        "structure_path": row["structure_path"],
        "job_dir": str(job_dir),
        "minimize_mdin": str(job_dir / "01_minimize.in"),
        "heat_mdin": str(job_dir / "02_heat.in"),
        "density_mdin": str(job_dir / "03_density.in"),
        "equilibrate_mdin": str(job_dir / "04_equilibrate.in"),
        "production_mdin": str(job_dir / "05_production.in"),
        "runner": str(job_dir / "run_amber_md.sh"),
        "status": "inputs_ready_needs_topology_and_md_engine",
        "source": "blind_stage2_md_queue",
    }


def write_manifest(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def generate(md_queue: Path, out_dir: Path, max_jobs: int | None) -> Path:
    rows = runnable_rows(read_queue(md_queue), max_jobs)
    manifest_rows = [write_job(row, out_dir) for row in rows]
    manifest = out_dir / "amber_md_job_manifest.tsv"
    write_manifest(manifest, manifest_rows)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--md-queue", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--max-jobs", type=int)
    args = parser.parse_args()
    manifest = generate(Path(args.md_queue), Path(args.out_dir), args.max_jobs)
    print(shlex.quote(str(manifest)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
