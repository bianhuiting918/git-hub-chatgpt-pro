#!/usr/bin/env python3
"""Launch the blind PETase Stage 1/2 compute gates and preserve run evidence."""

from __future__ import annotations

import argparse
import shlex
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path


BOUNDARY = (
    "Blind boundary: do not use article trajectories, article transition-state "
    "coordinates, article reaction coordinates, article selected CVs, article "
    "barriers, or article mechanism conclusions before final validation."
)


def quote_command(parts: list[str]) -> str:
    return " ".join(shlex.quote(part) for part in parts)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def build_runner_command(args: argparse.Namespace, out_root: Path) -> list[str]:
    command = [
        sys.executable,
        str(Path(args.runner)),
        "--project-root",
        str(Path(args.project_root)),
        "--phase-subdir",
        args.phase_subdir,
        "--out-root",
        str(out_root),
        "--conformers",
        str(args.conformers),
        "--replicates",
        str(args.replicates),
    ]
    if args.skip_shell_probes:
        command.append("--skip-shell-probes")
    return command


def launch(args: argparse.Namespace) -> int:
    project_root = Path(args.project_root)
    phase_root = project_root / args.phase_subdir
    out_root = Path(args.out_root) if args.out_root else phase_root / "blind_work"
    status_dir = out_root / "00_run_status"
    stdout_log = status_dir / "stage1_stage2_runner.stdout.log"
    stderr_log = status_dir / "stage1_stage2_runner.stderr.log"
    command_log = status_dir / "stage1_stage2_runner_command.txt"
    summary_path = status_dir / "compute_launch_summary.md"

    command = build_runner_command(args, out_root)
    started = datetime.now(timezone.utc)
    result = subprocess.run(command, cwd=str(project_root), text=True, capture_output=True)
    finished = datetime.now(timezone.utc)

    write_text(stdout_log, result.stdout)
    write_text(stderr_log, result.stderr)
    write_text(command_log, quote_command(command) + "\n")
    summary = "\n".join(
        [
            "# Blind Stage 1/2 Compute Launch Summary",
            "",
            BOUNDARY,
            "",
            f"- Started UTC: {started.isoformat()}",
            f"- Finished UTC: {finished.isoformat()}",
            f"- Exit code: {result.returncode}",
            f"- Project root: {project_root}",
            f"- Phase root: {phase_root}",
            f"- Output root: {out_root}",
            f"- Command log: {command_log}",
            f"- Stdout log: {stdout_log}",
            f"- Stderr log: {stderr_log}",
            "",
            "Interpretation:",
            "",
            "- Exit code 0 means all currently runnable gates completed.",
            "- Exit code 2 means the blind workflow is blocked at a recorded upstream gate.",
            "- Any other exit code requires inspecting stderr and the gate-status TSV.",
        ]
    )
    write_text(summary_path, summary + "\n")
    print(summary_path)
    return result.returncode


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", default=".")
    parser.add_argument("--phase-subdir", default="project01/phase2_blind_petase_qmmm_20260630")
    parser.add_argument("--out-root", default="")
    parser.add_argument("--runner", default="project01/phase2_blind_petase_qmmm_20260630/scripts/run_blind_stage1_stage2_gates.py")
    parser.add_argument("--conformers", type=int, default=20)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--skip-shell-probes", action="store_true")
    return launch(parser.parse_args())


if __name__ == "__main__":
    raise SystemExit(main())