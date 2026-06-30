#!/usr/bin/env python3
"""Prepare staged Amber MM/QM-MM cleanup minimization inputs for Stage 4 jobs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "cleanup_job_id",
    "mapped_job_id",
    "pose_id",
    "mapped_qmmask",
    "cleanup_job_dir",
    "prmtop_path",
    "inpcrd_path",
    "mm_solvent_mdin_path",
    "mm_all_mdin_path",
    "qmmm_short_mdin_path",
    "qmmm_full_mdin_path",
    "runner_path",
    "status",
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


def sanitize_id(raw: str) -> str:
    return "".join(char if char.isalnum() or char in {"_", "-"} else "_" for char in raw)


def mm_solvent_min(maxcyc: int, restraint_wt: float) -> str:
    return f"""Stage4 cleanup 00: MM solvent/ion minimization with restrained solute
&cntrl
  imin=1, maxcyc={maxcyc}, ncyc={max(1, maxcyc // 2)},
  ntb=1, cut=10.0,
  ntpr=10,
  ntr=1, restraint_wt={restraint_wt:.1f}, restraintmask='!:WAT,Na+,Cl-',
/
"""


def mm_all_min(maxcyc: int) -> str:
    return f"""Stage4 cleanup 01: all-atom MM minimization
&cntrl
  imin=1, maxcyc={maxcyc}, ncyc={max(1, maxcyc // 2)},
  ntb=1, cut=10.0,
  ntpr=10,
/
"""


def qmmm_min(maxcyc: int, qmmask: str, qmcharge: int, spin: int, title: str) -> str:
    return f"""{title}
&cntrl
  imin=1, maxcyc={maxcyc}, ncyc={max(1, maxcyc // 2)},
  ntb=1, cut=10.0,
  ntpr=1,
  ifqnt=1,
/
&qmmm
  qmmask='{qmmask}',
  qmcharge={qmcharge},
  spin={spin},
  qm_theory='DFTB3',
  qmshake=0,
/
"""


def runner(prmtop: Path, inpcrd: Path, run_full_qmmm_default: int) -> str:
    return f"""#!/usr/bin/env bash
set -euo pipefail

SANDER=${{SANDER:-sander}}
PRMTOP=${{PRMTOP:-{prmtop}}}
INPCRD=${{INPCRD:-{inpcrd}}}
RUN_FULL_QMMM=${{RUN_FULL_QMMM:-{run_full_qmmm_default}}}

if [ -z "${{AMBERHOME:-}}" ]; then
  SANDER_PATH=$(command -v "${{SANDER}}")
  AMBERHOME=$(cd "$(dirname "${{SANDER_PATH}}")/.." && pwd)
  export AMBERHOME
fi

"${{SANDER}}" -O -i 00_mm_solvent_ion_min.in -o 00_mm_solvent_ion_min.out \
  -p "${{PRMTOP}}" -c "${{INPCRD}}" -r 00_mm_solvent_ion_min.rst7 -ref "${{INPCRD}}"

"${{SANDER}}" -O -i 01_mm_all_min.in -o 01_mm_all_min.out \
  -p "${{PRMTOP}}" -c 00_mm_solvent_ion_min.rst7 -r 01_mm_all_min.rst7 -ref 00_mm_solvent_ion_min.rst7

"${{SANDER}}" -O -i 02_qmmm_short_min.in -o 02_qmmm_short_min.out \
  -p "${{PRMTOP}}" -c 01_mm_all_min.rst7 -r 02_qmmm_short_min.rst7 -ref 01_mm_all_min.rst7

if [ "${{RUN_FULL_QMMM}}" = "1" ]; then
  "${{SANDER}}" -O -i 03_qmmm_min.in -o 03_qmmm_min.out \
    -p "${{PRMTOP}}" -c 02_qmmm_short_min.rst7 -r 03_qmmm_min.rst7 -ref 02_qmmm_short_min.rst7
fi
"""


def readme(row: dict[str, str]) -> str:
    return "\n".join(
        [
            f"# Stage4 Amber Cleanup {row.get('mapped_job_id', '')}",
            "",
            "Purpose: reduce obvious close contacts before treating any DFTB3/MM minimization or 200 ps equilibration as chemically meaningful.",
            "",
            "Stages:",
            "1. MM minimization of solvent/ions with solute restrained.",
            "2. Short all-atom MM minimization.",
            "3. Short DFTB3/MM minimization using the mapped QM mask.",
            "4. Optional longer DFTB3/MM minimization when RUN_FULL_QMMM=1.",
            "",
            "Boundary: no paper coordinates, trajectories, reaction coordinates, windows, barriers, or mechanism conclusions are used.",
            "This cleanup is not a TS search and does not produce a barrier.",
        ]
    ) + "\n"


def write_job(
    row: dict[str, str],
    out_dir: Path,
    mm_solvent_maxcyc: int,
    mm_all_maxcyc: int,
    qmmm_short_maxcyc: int,
    qmmm_full_maxcyc: int,
    restraint_wt: float,
    run_full_qmmm_default: int,
) -> dict[str, str]:
    mapped_job_id = row["mapped_job_id"]
    cleanup_job_id = f"CLEANUP_{sanitize_id(mapped_job_id)}"
    job_dir = out_dir / cleanup_job_id
    job_dir.mkdir(parents=True, exist_ok=True)
    qmmask = row["mapped_qmmask"]
    prmtop = Path(row["prmtop_path"])
    inpcrd = Path(row["inpcrd_path"])

    files = {
        "00_mm_solvent_ion_min.in": mm_solvent_min(mm_solvent_maxcyc, restraint_wt),
        "01_mm_all_min.in": mm_all_min(mm_all_maxcyc),
        "02_qmmm_short_min.in": qmmm_min(qmmm_short_maxcyc, qmmask, 0, 1, "Stage4 cleanup 02: short DFTB3/MM minimization"),
        "03_qmmm_min.in": qmmm_min(qmmm_full_maxcyc, qmmask, 0, 1, "Stage4 cleanup 03: longer DFTB3/MM minimization"),
        "run_amber_cleanup.sh": runner(prmtop, inpcrd, run_full_qmmm_default),
        "00_README.md": readme(row),
    }
    for name, content in files.items():
        path = job_dir / name
        path.write_text(content, encoding="utf-8")
        if name.endswith(".sh"):
            path.chmod(0o755)

    return {
        "cleanup_job_id": cleanup_job_id,
        "mapped_job_id": mapped_job_id,
        "pose_id": row.get("pose_id", ""),
        "mapped_qmmask": qmmask,
        "cleanup_job_dir": str(job_dir),
        "prmtop_path": str(prmtop),
        "inpcrd_path": str(inpcrd),
        "mm_solvent_mdin_path": str(job_dir / "00_mm_solvent_ion_min.in"),
        "mm_all_mdin_path": str(job_dir / "01_mm_all_min.in"),
        "qmmm_short_mdin_path": str(job_dir / "02_qmmm_short_min.in"),
        "qmmm_full_mdin_path": str(job_dir / "03_qmmm_min.in"),
        "runner_path": str(job_dir / "run_amber_cleanup.sh"),
        "status": "cleanup_inputs_ready_requires_sander_execution",
    }


def generate(
    mapped_qmmm_manifest: Path,
    out_dir: Path,
    mm_solvent_maxcyc: int,
    mm_all_maxcyc: int,
    qmmm_short_maxcyc: int,
    qmmm_full_maxcyc: int,
    restraint_wt: float,
    run_full_qmmm_default: int,
) -> Path:
    rows = [
        write_job(
            row,
            out_dir,
            mm_solvent_maxcyc,
            mm_all_maxcyc,
            qmmm_short_maxcyc,
            qmmm_full_maxcyc,
            restraint_wt,
            run_full_qmmm_default,
        )
        for row in read_tsv(mapped_qmmm_manifest)
        if row.get("status") == "qmmm_inputs_mapped_to_amber_topology"
    ]
    manifest = out_dir / "amber_cleanup_manifest.tsv"
    write_tsv(manifest, rows, FIELDS)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mapped-qmmm-manifest", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--mm-solvent-maxcyc", type=int, default=200)
    parser.add_argument("--mm-all-maxcyc", type=int, default=100)
    parser.add_argument("--qmmm-short-maxcyc", type=int, default=5)
    parser.add_argument("--qmmm-full-maxcyc", type=int, default=2000)
    parser.add_argument("--restraint-wt", type=float, default=50.0)
    parser.add_argument("--run-full-qmmm-default", type=int, choices=[0, 1], default=0)
    args = parser.parse_args()
    manifest = generate(
        Path(args.mapped_qmmm_manifest),
        Path(args.out_dir),
        args.mm_solvent_maxcyc,
        args.mm_all_maxcyc,
        args.qmmm_short_maxcyc,
        args.qmmm_full_maxcyc,
        args.restraint_wt,
        args.run_full_qmmm_default,
    )
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
