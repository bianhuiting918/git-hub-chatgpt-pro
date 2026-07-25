#!/usr/bin/env python3
"""Extract explicit chain contexts from a PDB file without changing coordinates."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path


COORDINATE_RECORDS = ("ATOM  ", "HETATM")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--chains", required=True, help="Comma-separated PDB chain IDs")
    parser.add_argument(
        "--protein-only",
        action="store_true",
        help="Retain ATOM records only; exclude all HETATM records explicitly.",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    requested = tuple(dict.fromkeys(part.strip() for part in args.chains.split(",") if part.strip()))
    if not requested:
        print("no requested chains", file=sys.stderr)
        return 2
    if any(len(chain) != 1 for chain in requested):
        print("PDB chain IDs must be one character", file=sys.stderr)
        return 2
    if args.output.exists():
        print(f"refusing to overwrite existing output: {args.output}", file=sys.stderr)
        return 3

    selected = []
    found = set()
    with args.input.open(encoding="ascii", errors="strict") as handle:
        for line in handle:
            record = line[:6]
            if record in COORDINATE_RECORDS:
                if args.protein_only and record != "ATOM  ":
                    continue
                chain = line[21] if len(line) > 21 else ""
                if chain in requested:
                    selected.append(line.rstrip("\r\n") + "\n")
                    found.add(chain)
            elif record == "TER   ":
                chain = line[21] if len(line) > 21 else ""
                if chain in requested:
                    selected.append(line.rstrip("\r\n") + "\n")

    missing = [chain for chain in requested if chain not in found]
    if missing:
        print(f"missing requested chains: {','.join(missing)}", file=sys.stderr)
        return 4

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="ascii", newline="\n") as handle:
        handle.write(f"REMARK 999 SOURCE {args.input.resolve()}\n")
        handle.write(f"REMARK 999 SELECTED CHAINS {','.join(requested)}\n")
        record_policy = "PROTEIN_ONLY" if args.protein_only else "ATOM_AND_HETATM"
        handle.write(f"REMARK 999 RECORD POLICY {record_policy}\n")
        handle.writelines(selected)
        handle.write("END\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
