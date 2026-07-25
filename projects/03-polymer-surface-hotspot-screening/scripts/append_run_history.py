#!/usr/bin/env python3
"""Append one audited run record to a TSV ledger without duplicating run IDs."""

from __future__ import annotations

import argparse
import csv
import os
import sys
from pathlib import Path

try:
    import fcntl
except ImportError:  # pragma: no cover - Windows development fallback
    fcntl = None


COLUMNS = [
    "timestamp_local",
    "timestamp_utc",
    "git_commit",
    "run_id",
    "script",
    "parameters",
    "input_manifest",
    "input_sha256",
    "output_path",
    "scheduler_job_id",
    "exit_status",
    "n_input",
    "n_evaluated",
    "n_failed",
    "n_not_evaluated",
    "summary",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--history", type=Path, required=True)
    for column in COLUMNS:
        parser.add_argument("--" + column.replace("_", "-"), required=True)
    return parser.parse_args()


def clean(value: str, field: str) -> str:
    if "\t" in value or "\r" in value or "\n" in value:
        raise ValueError(f"{field} contains a forbidden tab or newline")
    return value


def main() -> int:
    args = parse_args()
    row = {column: clean(str(getattr(args, column)), column) for column in COLUMNS}
    history = args.history
    history.parent.mkdir(parents=True, exist_ok=True)

    with history.open("a+", encoding="utf-8", newline="") as handle:
        if fcntl is not None:
            fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            handle.seek(0)
            existing_text = handle.read()
            if existing_text:
                handle.seek(0)
                reader = csv.DictReader(handle, delimiter="\t")
                if reader.fieldnames != COLUMNS:
                    print("existing run history header does not match schema", file=sys.stderr)
                    return 3
                for existing in reader:
                    if existing["run_id"] == row["run_id"]:
                        print(f"duplicate run_id: {row['run_id']}", file=sys.stderr)
                        return 4
            handle.seek(0, os.SEEK_END)
            writer = csv.DictWriter(
                handle,
                fieldnames=COLUMNS,
                delimiter="\t",
                lineterminator="\n",
                extrasaction="raise",
            )
            if not existing_text:
                writer.writeheader()
            writer.writerow(row)
            handle.flush()
            os.fsync(handle.fileno())
        finally:
            if fcntl is not None:
                fcntl.flock(handle.fileno(), fcntl.LOCK_UN)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
