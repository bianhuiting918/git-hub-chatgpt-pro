#!/usr/bin/env python3
"""Validate an explicit-probe smoke manifest before any Slurm submission."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shlex
import subprocess
from pathlib import Path


REQUIRED_COLUMNS = (
    "material_family",
    "record_id",
    "probe_ids",
    "shell_dir",
    "patch_dir",
    "maps_root",
    "receptor_pdbqt",
    "output_dir",
)


def _inside(path: Path, root: Path) -> bool:
    try:
        path.resolve().relative_to(root.resolve())
        return True
    except ValueError:
        return False


def validate_row(row: dict, root: Path) -> dict:
    missing_columns = [column for column in REQUIRED_COLUMNS if not row.get(column)]
    if missing_columns:
        raise ValueError("missing columns: " + ",".join(missing_columns))
    family = row["material_family"].upper()
    if family not in {"PET", "NYLON"}:
        raise ValueError("invalid material family")
    probe_ids = [value for value in row["probe_ids"].split(",") if value]
    if not probe_ids:
        raise ValueError("missing probe IDs")
    expected_prefix = "PET_" if family == "PET" else ("NYL_", "PA6_", "PA66_")
    if family == "PET":
        valid_probe_family = all(value.startswith(expected_prefix) for value in probe_ids)
    else:
        valid_probe_family = all(value.startswith(expected_prefix) for value in probe_ids)
    if not valid_probe_family:
        raise ValueError("probe family does not match material family")

    shell_dir = Path(row["shell_dir"]).resolve()
    patch_dir = Path(row["patch_dir"]).resolve()
    maps_root = Path(row["maps_root"]).resolve()
    receptor = Path(row["receptor_pdbqt"]).resolve()
    output = Path(row["output_dir"]).resolve()
    for path in (shell_dir, patch_dir, maps_root, receptor):
        if not _inside(path, root):
            raise ValueError(f"input path outside project root: {path}")
    if not _inside(output, root):
        raise ValueError(f"output path outside project root: {output}")
    required_inputs = (
        shell_dir / "SHELL_PASS.json",
        patch_dir / "PATCH_PASS.json",
        maps_root / "MAPS_PASS.json",
        receptor,
    )
    missing_inputs = [str(path) for path in required_inputs if not path.is_file()]
    if missing_inputs:
        raise ValueError("missing validated input: " + ",".join(missing_inputs))
    allowed_results = (root / "results").resolve()
    if not _inside(output, allowed_results) or not output.parent.name.startswith(
        "probe_match_smoke_"
    ):
        raise ValueError("output is outside an allowed probe_match_smoke root")
    if output.exists():
        raise ValueError(f"output already exists: {output}")

    validated = dict(row)
    validated.update(
        {
            "material_family": family,
            "probe_ids": ",".join(probe_ids),
            "shell_dir": str(shell_dir),
            "patch_dir": str(patch_dir),
            "maps_root": str(maps_root),
            "receptor_pdbqt": str(receptor),
            "output_dir": str(output),
        }
    )
    return validated


def build_sbatch_arguments(
    row: dict,
    *,
    root: Path,
    commit: str,
    wrapper: Path,
    log_dir: Path,
) -> list[str]:
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("invalid commit")
    job_slug = re.sub(r"[^A-Za-z0-9_]+", "_", row["record_id"]).lower()
    command = [
        "bash",
        str(wrapper),
        "--commit",
        commit,
        "--material-family",
        row["material_family"],
        "--record-id",
        row["record_id"],
        "--probe-ids",
        row["probe_ids"],
        "--shell-dir",
        row["shell_dir"],
        "--patch-dir",
        row["patch_dir"],
        "--maps-root",
        row["maps_root"],
        "--receptor-pdbqt",
        row["receptor_pdbqt"],
        "--output-dir",
        row["output_dir"],
    ]
    return [
        "sbatch",
        "--parsable",
        "-p",
        "xahcnormal",
        "--job-name",
        f"probe_{job_slug}",
        "--cpus-per-task",
        "1",
        "--mem",
        "3G",
        "--time",
        "02:00:00",
        "--output",
        str(log_dir / f"{job_slug}_%j.out"),
        "--error",
        str(log_dir / f"{job_slug}_%j.err"),
        "--wrap",
        shlex.join(command),
    ]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--commit", required=True)
    parser.add_argument("--wrapper", type=Path, required=True)
    parser.add_argument("--log-dir", type=Path, required=True)
    parser.add_argument("--audit-json", type=Path, required=True)
    parser.add_argument("--dry-run", action="store_true")
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.audit_json.exists():
        raise SystemExit(f"refusing existing audit: {args.audit_json}")
    with args.manifest.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = list(reader)
    if not rows:
        raise SystemExit("empty smoke manifest")
    validated = [validate_row(row, args.root) for row in rows]
    outputs = [row["output_dir"] for row in validated]
    if len(outputs) != len(set(outputs)):
        raise SystemExit("duplicate output directories in manifest")
    if not args.wrapper.is_file():
        raise SystemExit("missing wrapper")
    args.log_dir.mkdir(parents=True, exist_ok=True)

    submissions = []
    for row in validated:
        command = build_sbatch_arguments(
            row,
            root=args.root,
            commit=args.commit,
            wrapper=args.wrapper,
            log_dir=args.log_dir,
        )
        if args.dry_run:
            submissions.append(
                {"record_id": row["record_id"], "status": "DRY_RUN", "command": command}
            )
            continue
        completed = subprocess.run(command, text=True, capture_output=True, check=False)
        if completed.returncode != 0:
            raise SystemExit(
                f"submission failed for {row['record_id']}: {completed.stderr.strip()}"
            )
        submissions.append(
            {
                "record_id": row["record_id"],
                "status": "SUBMITTED",
                "job_id": completed.stdout.strip(),
                "command": command,
            }
        )
    payload = {
        "status": "SUBMISSION_DRY_RUN" if args.dry_run else "SUBMISSION_PASS",
        "manifest": str(args.manifest.resolve()),
        "commit": args.commit,
        "n_rows": len(validated),
        "submissions": submissions,
    }
    args.audit_json.parent.mkdir(parents=True, exist_ok=True)
    with args.audit_json.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")
    for item in submissions:
        print(item.get("job_id", "DRY_RUN"), item["record_id"])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
