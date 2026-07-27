#!/usr/bin/env python3
"""Run independent surface-field rows concurrently inside one Slurm job."""
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--family", choices=("PET", "NYLON"), required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--existing-pass-root", type=Path, required=True)
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--tool-root", type=Path, required=True)
    parser.add_argument("--adfr-root", type=Path, required=True)
    parser.add_argument("--row-worker", type=Path, required=True)
    parser.add_argument("--start-row", type=int, required=True)
    parser.add_argument("--end-row", type=int, required=True)
    parser.add_argument("--workers", type=int, required=True)
    return parser.parse_args()


def require_under_root(path: Path, root: Path) -> Path:
    resolved = path.resolve()
    allowed = root.resolve()
    if resolved != allowed and allowed not in resolved.parents:
        raise ValueError(f"path outside project root: {resolved}")
    return resolved


def load_rows(path: Path, start: int, end: int) -> list[tuple[int, dict[str, str]]]:
    if start < 0 or end <= start:
        raise ValueError(f"invalid row interval: [{start}, {end})")
    selected: list[tuple[int, dict[str, str]]] = []
    with path.open(newline="", encoding="utf-8") as handle:
        for index, row in enumerate(csv.DictReader(handle, delimiter="\t")):
            if index >= end:
                break
            if index >= start:
                selected.append((index, {key: (value or "").strip() for key, value in row.items()}))
    return selected


def valid_pass(root: Path, sequence_md5: str) -> bool:
    gate = root / sequence_md5 / "FIELD_PASS.json"
    if not gate.is_file():
        return False
    try:
        return json.loads(gate.read_text(encoding="utf-8")).get("status") == "FIELD_PASS"
    except (OSError, json.JSONDecodeError):
        return False


def run_row(args: argparse.Namespace, row_index: int) -> tuple[int, int, str, str]:
    command = [
        sys.executable,
        str(args.row_worker),
        "--manifest",
        str(args.manifest),
        "--row-index",
        str(row_index),
        "--family",
        args.family,
        "--output-root",
        str(args.output_root),
        "--project-root",
        str(args.project_root),
        "--tool-root",
        str(args.tool_root),
        "--adfr-root",
        str(args.adfr_root),
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    return row_index, completed.returncode, completed.stdout, completed.stderr


def main() -> int:
    args = parse_args()
    if not 1 <= args.workers <= 64:
        raise ValueError("--workers must be between 1 and 64")
    args.output_root = require_under_root(args.output_root, args.project_root)
    args.existing_pass_root = require_under_root(args.existing_pass_root, args.project_root)
    if args.output_root == args.existing_pass_root:
        raise ValueError("recovery output root must differ from existing pass root")
    if not args.manifest.is_file() or not args.row_worker.is_file():
        raise FileNotFoundError("manifest or row worker missing")

    rows = load_rows(args.manifest, args.start_row, args.end_row)
    runnable: list[int] = []
    skipped_primary_pass = 0
    for row_index, row in rows:
        sequence_md5 = row.get("sequence_md5", "")
        if len(sequence_md5) == 32 and valid_pass(args.existing_pass_root, sequence_md5):
            skipped_primary_pass += 1
        else:
            runnable.append(row_index)

    failed = 0
    completed_count = 0
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        futures = {pool.submit(run_row, args, row_index): row_index for row_index in runnable}
        for future in as_completed(futures):
            row_index, returncode, stdout, stderr = future.result()
            if stdout:
                print(stdout, end="" if stdout.endswith("\n") else "\n")
            if stderr:
                print(stderr, end="" if stderr.endswith("\n") else "\n", file=sys.stderr)
            completed_count += 1
            if returncode:
                failed += 1
                print(f"ROW_NOT_EVALUATED row={row_index} rc={returncode}", file=sys.stderr)

    print(
        "PARALLEL_CHUNK_COMPLETE "
        f"family={args.family} start={args.start_row} end={args.end_row} "
        f"workers={args.workers} manifest_rows={len(rows)} "
        f"skipped_primary_pass={skipped_primary_pass} "
        f"processed={completed_count} failed={failed}"
    )
    return 10 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
