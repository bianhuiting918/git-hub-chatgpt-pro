#!/usr/bin/env python3
"""Freeze auditable PET and nylon structure manifests for surface screening."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path
from typing import Dict, Iterable, List, Tuple


BASE_FIELDS = [
    "material_family",
    "manifest_role",
    "source_manifest",
    "source_row_number",
    "source_priority",
    "candidate_id",
    "sequence_md5",
    "selected_chain",
    "receptor_path",
    "receptor_sha256",
    "receptor_observed_sha256",
    "input_status",
    "structure_provenance",
]
REJECT_FIELDS = BASE_FIELDS + ["exclusion_reason"]


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_rows(path: Path) -> Iterable[Tuple[int, Dict[str, str]]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"candidate_id", "sequence_md5", "receptor_path", "receptor_sha256"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path}: missing required columns: {sorted(missing)}")
        for row_number, row in enumerate(reader, start=2):
            yield row_number, {key: (value or "").strip() for key, value in row.items()}


def normalized_row(
    row: Dict[str, str],
    *,
    family: str,
    role: str,
    manifest: Path,
    row_number: int,
    priority: int,
) -> Dict[str, str]:
    return {
        "material_family": family,
        "manifest_role": role,
        "source_manifest": str(manifest.resolve()),
        "source_row_number": str(row_number),
        "source_priority": str(priority),
        "candidate_id": row.get("candidate_id", ""),
        "sequence_md5": row.get("sequence_md5", ""),
        "selected_chain": row.get("selected_chain", ""),
        "receptor_path": row.get("receptor_path", ""),
        "receptor_sha256": row.get("receptor_sha256", "").lower(),
        "receptor_observed_sha256": row.get("receptor_observed_sha256", "").lower(),
        "input_status": row.get("input_status", ""),
        "structure_provenance": row.get("structure_provenance", ""),
    }


def validate(candidate: Dict[str, str]) -> str:
    if candidate["input_status"] != "READY_FOR_STRUCTURE_EVALUATION":
        return "NOT_EVALUATED_INPUT_STATUS"
    if not candidate["candidate_id"] or not candidate["sequence_md5"]:
        return "NOT_EVALUATED_MISSING_IDENTITY"
    receptor_text = candidate["receptor_path"]
    if not receptor_text:
        return "NOT_EVALUATED_MISSING_STRUCTURE"
    receptor = Path(receptor_text)
    try:
        if not receptor.is_file():
            return "NOT_EVALUATED_MISSING_STRUCTURE"
        with receptor.open("rb") as handle:
            handle.read(1)
    except OSError:
        return "NOT_EVALUATED_UNREADABLE_STRUCTURE"
    expected = candidate["receptor_sha256"]
    if len(expected) != 64:
        return "NOT_EVALUATED_MISSING_CHECKSUM"
    try:
        int(expected, 16)
    except ValueError:
        return "NOT_EVALUATED_INVALID_CHECKSUM"
    try:
        observed = file_sha256(receptor)
    except OSError:
        return "NOT_EVALUATED_UNREADABLE_STRUCTURE"
    candidate["receptor_observed_sha256"] = observed
    if observed != expected:
        return "NOT_EVALUATED_CHECKSUM_MISMATCH"
    return ""


def write_tsv(path: Path, rows: List[Dict[str, str]], fields: List[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pet-manifest", type=Path, required=True)
    parser.add_argument("--nylon-manifest", type=Path, required=True)
    parser.add_argument("--nylon-recovered-manifest", type=Path, required=True)
    parser.add_argument("--nylon-extra-manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    sources = [
        ("PET", "pet_primary", args.pet_manifest, 0),
        ("NYLON", "nylon_primary", args.nylon_manifest, 0),
        ("NYLON", "nylon_recovered", args.nylon_recovered_manifest, 1),
        ("NYLON", "nylon_extra_physical", args.nylon_extra_manifest, 2),
    ]
    included = {"PET": [], "NYLON": []}
    rejected: List[Dict[str, str]] = []
    seen = {"PET": set(), "NYLON": set()}
    input_counts = Counter()

    for family, role, manifest, priority in sources:
        for row_number, row in read_rows(manifest):
            input_counts[role] += 1
            candidate = normalized_row(
                row,
                family=family,
                role=role,
                manifest=manifest,
                row_number=row_number,
                priority=priority,
            )
            reason = validate(candidate)
            key = candidate["sequence_md5"]
            if not reason and key in seen[family]:
                reason = "DUPLICATE_LOWER_PRIORITY"
            if reason:
                candidate["exclusion_reason"] = reason
                rejected.append(candidate)
                continue
            seen[family].add(key)
            included[family].append(candidate)

    for family in included:
        included[family].sort(key=lambda row: (row["sequence_md5"], row["candidate_id"]))
    rejected.sort(
        key=lambda row: (
            row["material_family"],
            int(row["source_priority"]),
            int(row["source_row_number"]),
        )
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    write_tsv(args.output_dir / "pet_structures.tsv", included["PET"], BASE_FIELDS)
    write_tsv(args.output_dir / "nylon_structures.tsv", included["NYLON"], BASE_FIELDS)
    write_tsv(args.output_dir / "not_evaluated.tsv", rejected, REJECT_FIELDS)

    reason_counts = Counter(row["exclusion_reason"] for row in rejected)
    summary = {
        "schema_version": "structure_manifest_v1",
        "input_rows": dict(sorted(input_counts.items())),
        "included": {family: len(rows) for family, rows in included.items()},
        "not_evaluated_or_duplicate": len(rejected),
        "exclusion_reasons": dict(sorted(reason_counts.items())),
        "deduplication_key": "material_family + sequence_md5",
        "source_priority": {
            "pet_primary": 0,
            "nylon_primary": 0,
            "nylon_recovered": 1,
            "nylon_extra_physical": 2,
        },
    }
    (args.output_dir / "structure_manifest_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
