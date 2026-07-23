#!/usr/bin/env python3
"""Resolve and extract the five immutable FREE_MONITOR source frames."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
from pathlib import Path
from typing import Any, Dict


def load_manifest(path: Path) -> Dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    candidates = data.get("candidates", [])
    if data.get("formal_candidate_count") != 5 or len(candidates) != 5:
        raise ValueError("Formal candidate universe must contain exactly five rows")
    ids = [row["id"] for row in candidates]
    if len(set(ids)) != len(ids):
        raise ValueError("Candidate IDs are not unique")
    return data


def _load_geometry_rows(path: Path) -> list[Dict[str, Any]]:
    rows = []
    with path.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle, delimiter="\t"):
            rows.append(
                {
                    "segment": raw["segment"],
                    "action": raw["action"],
                    "simulation_time_ps": float(raw["simulation_time_ps"]),
                    "gate_opening_nm": float(raw["gate_opening_nm"]),
                    "distance_nm": float(raw["distance_nm"]),
                    "angle_deg": float(raw["angle_deg"]),
                    "joint_pass": int(raw["joint_pass"]),
                    "substrate_restrained": int(raw["substrate_restrained"]),
                    "gate_restrained": int(raw["gate_restrained"]),
                }
            )
    return rows


def resolve_candidate(candidate: Dict[str, Any]) -> Dict[str, Any]:
    branch = Path(candidate["branch_dir"])
    segment = Path(candidate["segment_dir"])
    geometry_path = branch / "geometry.tsv"
    rows = _load_geometry_rows(geometry_path)
    segment_rows = [
        row for row in rows if row["segment"] == candidate["segment"]
    ]
    if not segment_rows:
        raise ValueError(f"No geometry rows for {candidate['id']}")
    selected = [
        row
        for row in segment_rows
        if abs(row["simulation_time_ps"] - candidate["cumulative_time_ps"]) < 1e-6
    ]
    if len(selected) != 1:
        raise ValueError(
            f"Expected one selected geometry row for {candidate['id']}, got {len(selected)}"
        )
    selected_row = selected[0]
    if (
        selected_row["action"] != "FREE_MONITOR"
        or selected_row["joint_pass"] != 1
        or selected_row["substrate_restrained"] != 0
        or selected_row["gate_restrained"] != 0
    ):
        raise ValueError(f"Selected frame is not an unrestrained NAC: {candidate['id']}")
    geometry_times = sorted({row["simulation_time_ps"] for row in segment_rows})
    positive_steps = [
        right - left
        for left, right in zip(geometry_times, geometry_times[1:])
        if right > left
    ]
    if not positive_steps:
        raise ValueError(f"Cannot derive geometry sample interval for {candidate['id']}")
    sample_interval = min(positive_steps)
    if any(abs(step - sample_interval) > 1e-6 for step in positive_steps):
        raise ValueError(f"Nonuniform geometry sampling for {candidate['id']}")
    first_time = geometry_times[0]
    segment_start_time = first_time - sample_interval
    local_time = selected_row["simulation_time_ps"] - segment_start_time
    if abs(local_time - candidate["local_xtc_time_ps"]) > 1e-6:
        raise ValueError(
            f"Local time mismatch for {candidate['id']}: {local_time} vs "
            f"{candidate['local_xtc_time_ps']}"
        )
    resolved = {
        "candidate_id": candidate["id"],
        "geometry_tsv": str(geometry_path),
        "xtc": str(segment / "run.xtc"),
        "tpr": str(segment / "run.tpr"),
        "source_itp": str(branch / candidate["source_itp"]),
        "source_topology": str(branch / "topol.top"),
        "source_index": str(branch / "cycle.ndx"),
        "selected_row": selected_row,
        "segment_first_geometry_time_ps": first_time,
        "geometry_sample_interval_ps": sample_interval,
        "segment_start_time_ps": segment_start_time,
        "derived_local_xtc_time_ps": local_time,
    }
    missing = [
        value
        for key, value in resolved.items()
        if key in {"xtc", "tpr", "source_itp", "source_topology", "source_index"}
        and not Path(value).is_file()
    ]
    if missing:
        raise FileNotFoundError(f"Missing source inputs for {candidate['id']}: {missing}")
    return resolved


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def gro_atom_count(path: Path) -> int:
    with path.open(encoding="utf-8") as handle:
        handle.readline()
        return int(handle.readline().strip())


def temporary_gro_path(output_gro: Path) -> Path:
    return output_gro.with_name(output_gro.stem + ".tmp.gro")


def extract_all(
    manifest_path: Path, task_root: Path, gmx: str, force: bool = False
) -> Dict[str, Any]:
    manifest = load_manifest(manifest_path)
    inputs_root = task_root / "inputs"
    inputs_root.mkdir(parents=True, exist_ok=True)
    resolved_rows = []
    hash_rows = []
    for candidate in manifest["candidates"]:
        resolved = resolve_candidate(candidate)
        candidate_dir = inputs_root / candidate["id"]
        candidate_dir.mkdir(parents=True, exist_ok=True)
        output_gro = candidate_dir / "source.gro"
        if force or not output_gro.exists():
            temporary_gro = temporary_gro_path(output_gro)
            subprocess.run(
                [
                    gmx,
                    "trjconv",
                    "-f",
                    resolved["xtc"],
                    "-s",
                    resolved["tpr"],
                    "-dump",
                    str(resolved["derived_local_xtc_time_ps"]),
                    "-o",
                    str(temporary_gro),
                ],
                input="0\n",
                text=True,
                check=True,
            )
            temporary_gro.replace(output_gro)
        expected_atoms = gro_atom_count(Path(resolved["tpr"]).with_name("run.gro"))
        if gro_atom_count(output_gro) != expected_atoms:
            raise ValueError(f"Extracted GRO atom count mismatch for {candidate['id']}")
        resolved["extracted_gro"] = str(output_gro)
        resolved["extracted_atom_count"] = gro_atom_count(output_gro)
        resolved["hashes"] = {}
        for key in ("xtc", "tpr", "source_itp", "source_topology", "source_index"):
            value = Path(resolved[key])
            resolved["hashes"][key] = sha256(value)
            hash_rows.append((candidate["id"], key, str(value), resolved["hashes"][key]))
        resolved["hashes"]["extracted_gro"] = sha256(output_gro)
        hash_rows.append(
            (
                candidate["id"],
                "extracted_gro",
                str(output_gro),
                resolved["hashes"]["extracted_gro"],
            )
        )
        resolved_rows.append(resolved)

    result = {
        "schema_version": 1,
        "formal_candidate_count": 5,
        "source_manifest": str(manifest_path),
        "source_manifest_sha256": sha256(manifest_path),
        "candidates": resolved_rows,
    }
    resolved_path = inputs_root / "source_manifest.resolved.json"
    temp = resolved_path.with_suffix(".json.tmp")
    temp.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    temp.replace(resolved_path)
    hash_path = inputs_root / "source_sha256.tsv"
    lines = ["candidate_id\tinput_kind\tpath\tsha256"]
    lines.extend("\t".join(row) for row in hash_rows)
    hash_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--task-root", type=Path, required=True)
    parser.add_argument("--gmx", required=True)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    result = extract_all(args.manifest, args.task_root, args.gmx, force=args.force)
    print(json.dumps({"candidate_count": len(result["candidates"])}, indent=2))


if __name__ == "__main__":
    main()
