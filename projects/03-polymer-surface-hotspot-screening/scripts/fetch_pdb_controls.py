#!/usr/bin/env python3
"""Fetch a frozen set of official RCSB PDB controls with auditable provenance."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import sys
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import NamedTuple


class Source(NamedTuple):
    pdb_id: str
    url: str
    filename: str
    compressed: bool = False
    allow_legacy_assembly_placeholder: bool = False


FROZEN_SOURCES = {
    "1WYB": Source(
        "1WYB",
        "https://files.rcsb.org/download/1WYB.pdb",
        "1WYB.pdb",
    ),
    "9DYS": Source(
        "9DYS",
        "https://files.rcsb.org/download/9DYS.pdb",
        "9DYS.pdb",
    ),
    "3AXG_BIOASSEMBLY1": Source(
        "3AXG",
        "https://files.rcsb.org/download/3AXG.pdb1.gz",
        "3AXG_BIOASSEMBLY1.pdb",
        True,
        True,
    ),
}


def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def decode_download(data: bytes, compressed: bool) -> bytes:
    return gzip.decompress(data) if compressed else data


def validate_pdb_bytes(
    data: bytes,
    expected_pdb_id: str,
    allow_legacy_assembly_placeholder: bool = False,
) -> dict:
    try:
        text = data.decode("ascii")
    except UnicodeDecodeError as error:
        raise ValueError("PDB payload is not ASCII text") from error

    prefix = text[:20000].upper()
    if expected_pdb_id.upper() in prefix:
        identity_validation = "expected_pdb_id_in_prefix"
    else:
        header = next(
            (line for line in text.splitlines() if line.startswith("HEADER")),
            "",
        )
        legacy_placeholder = header.rstrip().upper().endswith("XXXX")
        if not (allow_legacy_assembly_placeholder and legacy_placeholder):
            raise ValueError(f"expected PDB id {expected_pdb_id} not found")
        identity_validation = "rcsb_legacy_bioassembly_header_xxxx"

    n_coordinates = sum(
        line.startswith(("ATOM  ", "HETATM")) for line in text.splitlines()
    )
    if n_coordinates == 0:
        raise ValueError("PDB payload has no coordinate records")
    if not text.rstrip().endswith("END"):
        raise ValueError("PDB payload lacks terminal END record")
    return {
        "n_coordinate_records": int(n_coordinates),
        "identity_validation": identity_validation,
    }


def prepare_output_directory(path: Path) -> None:
    if path.exists():
        raise FileExistsError(f"refusing existing output directory: {path}")
    path.mkdir(parents=True)


def fetch(url: str, timeout_seconds: int = 60) -> bytes:
    request = urllib.request.Request(
        url,
        headers={"User-Agent": "polymer-surface-hotspot-screen/phase1"},
    )
    with urllib.request.urlopen(request, timeout=timeout_seconds) as response:
        if getattr(response, "status", 200) != 200:
            raise RuntimeError(f"HTTP status {response.status} for {url}")
        return response.read()


def write_json_exclusive(path: Path, payload: dict) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725"),
    )
    parser.add_argument("--timeout-seconds", type=int, default=60)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    root = args.project_root.resolve()
    output = args.output_dir.resolve()
    if root != output and root not in output.parents:
        print("refusing output outside project root", file=sys.stderr)
        return 2
    if args.timeout_seconds <= 0:
        print("timeout must be positive", file=sys.stderr)
        return 2
    try:
        prepare_output_directory(args.output_dir)
    except FileExistsError as error:
        print(str(error), file=sys.stderr)
        return 3

    downloaded_at = datetime.now(timezone.utc).isoformat()
    rows = []
    try:
        for source_id, source in FROZEN_SOURCES.items():
            wire = fetch(source.url, args.timeout_seconds)
            decoded = decode_download(wire, source.compressed)
            validation = validate_pdb_bytes(
                decoded,
                source.pdb_id,
                source.allow_legacy_assembly_placeholder,
            )
            destination = args.output_dir / source.filename
            destination.write_bytes(decoded)
            rows.append(
                {
                    "source_id": source_id,
                    "pdb_id": source.pdb_id,
                    "url": source.url,
                    "downloaded_at_utc": downloaded_at,
                    "wire_bytes": len(wire),
                    "decoded_bytes": len(decoded),
                    "sha256": sha256_bytes(decoded),
                    "n_coordinate_records": validation["n_coordinate_records"],
                    "identity_validation": validation["identity_validation"],
                }
            )
        provenance = args.output_dir / "provenance.tsv"
        with provenance.open("x", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        write_json_exclusive(
            args.output_dir / "PDB_CONTROLS_PASS.json",
            {
                "status": "PDB_CONTROLS_PASS",
                "scientific_scope": "official_structure_retrieval_only",
                "n_sources": len(rows),
                "source_ids": list(FROZEN_SOURCES),
                "downloaded_at_utc": downloaded_at,
            },
        )
        checksum_lines = [
            f"{sha256_file(path)}  {path.name}"
            for path in sorted(args.output_dir.iterdir())
            if path.is_file() and path.name != "SHA256SUMS"
        ]
        (args.output_dir / "SHA256SUMS").write_text(
            "\n".join(checksum_lines) + "\n", encoding="ascii"
        )
    except Exception as error:
        write_json_exclusive(
            args.output_dir / "PDB_CONTROLS_FAIL.json",
            {
                "status": "NOT_EVALUATED_PDB_CONTROL_DOWNLOAD",
                "error": str(error),
                "downloaded_at_utc": downloaded_at,
            },
        )
        print(str(error), file=sys.stderr)
        return 4
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
