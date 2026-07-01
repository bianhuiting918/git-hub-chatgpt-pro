#!/usr/bin/env python3
"""Analyze Amber/Sander QM/MM output health for Stage 4 cleanup gates."""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path


FIELDS = [
    "case_id",
    "out_path",
    "stderr_path",
    "status",
    "run_done",
    "final_results",
    "scc_convergence_warnings",
    "last_step",
    "last_vdwaals",
    "last_dftbescf",
    "stderr_nonempty",
]

STEP_RE = re.compile(r"^\s*(\d+)\s+[-+]?\d+\.\d+E[+-]\d+\s+", re.MULTILINE)
VDWAALS_RE = re.compile(r"VDWAALS\s*=\s*([-+*0-9.Ee]+)")
DFTBESCF_RE = re.compile(r"DFTBESCF\s*=\s*([-+*0-9.Ee]+)")


def read_text(path: Path) -> str:
    if not path.exists():
        return ""
    return path.read_text(encoding="utf-8", errors="replace")


def last_match(pattern: re.Pattern[str], text: str) -> str:
    value = ""
    for match in pattern.finditer(text):
        value = match.group(1)
    return value


def analyze_case(case_id: str, out_path: Path, stderr_path: Path) -> dict[str, str]:
    out_text = read_text(out_path)
    stderr_text = read_text(stderr_path)
    run_done = "Run   done" in out_text
    final_results = "FINAL RESULTS" in out_text
    warning_count = out_text.count("Convergence could not be achieved in this step")
    stderr_nonempty = bool(stderr_text.strip())

    if run_done and final_results and warning_count:
        status = "completed_with_scc_warnings"
    elif run_done and final_results:
        status = "completed"
    elif stderr_nonempty or "ERROR" in out_text or "Error" in out_text:
        status = "failed_or_error_output"
    else:
        status = "running_or_incomplete"

    return {
        "case_id": case_id,
        "out_path": str(out_path),
        "stderr_path": str(stderr_path),
        "status": status,
        "run_done": "1" if run_done else "0",
        "final_results": "1" if final_results else "0",
        "scc_convergence_warnings": str(warning_count),
        "last_step": last_match(STEP_RE, out_text),
        "last_vdwaals": last_match(VDWAALS_RE, out_text),
        "last_dftbescf": last_match(DFTBESCF_RE, out_text),
        "stderr_nonempty": "1" if stderr_nonempty else "0",
    }


def parse_case(raw: str) -> tuple[str, Path, Path]:
    try:
        case_id, paths = raw.split("=", 1)
        parts = re.split(r":(?=[A-Za-z]:[\\/])", paths, maxsplit=1)
        if len(parts) == 1:
            parts = paths.split(":", 1)
        out_raw, stderr_raw = parts
    except ValueError as exc:
        raise argparse.ArgumentTypeError("case must be CASE_ID=OUT_PATH:STDERR_PATH") from exc
    return case_id, Path(out_raw), Path(stderr_raw)


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-manifest", required=True)
    parser.add_argument("--case", action="append", default=[], type=parse_case, help="CASE_ID=OUT_PATH:STDERR_PATH")
    args = parser.parse_args()
    rows = [analyze_case(case_id, out_path, stderr_path) for case_id, out_path, stderr_path in args.case]
    write_tsv(Path(args.output_manifest), rows)
    print(args.output_manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())



