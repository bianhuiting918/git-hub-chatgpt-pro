#!/usr/bin/env python3
"""Assemble compact raw primitives for the independent NylC M1 final audit."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import pathlib
from collections import Counter
from typing import Callable

import numpy as np

EXPECTED_SEEDS = (26711, 26723, 26737)


def _read_table(path: pathlib.Path, columns: int) -> list[list[float]]:
    rows = []
    for raw in pathlib.Path(path).read_text().splitlines():
        line = raw.strip()
        if not line or line.startswith(("#", "@")):
            continue
        values = [float(value) for value in line.split()]
        if len(values) < columns:
            raise ValueError(f"{path}: expected at least {columns} columns")
        row = values[:columns]
        if not all(math.isfinite(value) for value in row):
            raise ValueError(f"{path}: non-finite value")
        rows.append(row)
    if not rows:
        raise ValueError(f"{path}: no data rows")
    return rows


def _same_times(reference: list[list[float]], other: list[list[float]]) -> bool:
    return len(reference) == len(other) and all(
        abs(left[0] - right[0]) <= 1.0e-6
        for left, right in zip(reference, other)
    )


def load_raw_frames(root: pathlib.Path) -> list[dict]:
    root = pathlib.Path(root)
    distance = _read_table(root / "nac_distance.xvg", 2)
    angle = _read_table(root / "nac_angle.xvg", 2)
    gate = _read_table(root / "gate_opening.xvg", 2)
    pocket = _read_table(root / "pocket_state.xvg", 3)
    for label, rows in (("angle", angle), ("gate", gate), ("pocket", pocket)):
        if not _same_times(distance, rows):
            raise ValueError(f"primitive time axes differ for {label}")
    if any(right[0] <= left[0] for left, right in zip(distance, distance[1:])):
        raise ValueError("primitive time axes are duplicated or unordered")
    return [
        {
            "time_ps": float(drow[0]),
            "distance_nm": float(drow[1]),
            "angle_deg": float(arow[1]),
            "gate_opening_nm": float(grow[1]),
            "pocket_contact_count": int(round(prow[1])),
            "ligand_pocket_com_nm": float(prow[2]),
        }
        for drow, arow, grow, prow in zip(distance, angle, gate, pocket)
    ]


def _add_proton_annotations(root: pathlib.Path, frames: list[dict]) -> None:
    path = pathlib.Path(root) / "proton_network.jsonl"
    if not path.is_file():
        raise ValueError(f"{path}: missing proton-network primitives")
    rows = {}
    for raw in path.read_text().splitlines():
        if not raw.strip():
            continue
        row = json.loads(raw)
        key = round(float(row["time_ps"]), 6)
        if key in rows:
            raise ValueError("duplicate proton-network time")
        rows[key] = row
    for frame in frames:
        key = round(float(frame["time_ps"]), 6)
        if key not in rows:
            raise ValueError("proton-network and geometry time axes differ")
        channels = rows[key].get("channels") or {}
        favorable = sorted(
            name
            for name, metrics in channels.items()
            if isinstance(metrics, dict) and bool(metrics.get("favorable"))
        )
        frame["proton_route_annotation"] = (
            "|".join(favorable)
            if favorable
            else "NO_FAVORABLE_FIXED_TOPOLOGY_CHANNEL"
        )


def _align_local(reference_alignment, alignment, local):
    reference_alignment = np.asarray(reference_alignment, dtype=float)
    alignment = np.asarray(alignment, dtype=float)
    local = np.asarray(local, dtype=float)
    if alignment.shape != reference_alignment.shape or alignment.ndim != 2:
        raise ValueError("alignment fingerprints have inconsistent shapes")
    moving_center = alignment.mean(axis=0)
    reference_center = reference_alignment.mean(axis=0)
    moving = alignment - moving_center
    reference = reference_alignment - reference_center
    left, _, right_t = np.linalg.svd(moving.T @ reference)
    rotation = left @ right_t
    if np.linalg.det(rotation) < 0:
        left[:, -1] *= -1
        rotation = left @ right_t
    return (local - moving_center) @ rotation + reference_center


def _components(adjacency: list[list[int]]) -> list[list[int]]:
    seen = set()
    result = []
    for start in range(len(adjacency)):
        if start in seen:
            continue
        stack = [start]
        seen.add(start)
        component = []
        while stack:
            node = stack.pop()
            component.append(node)
            for neighbor in adjacency[node]:
                if neighbor not in seen:
                    seen.add(neighbor)
                    stack.append(neighbor)
        result.append(sorted(component))
    return result


def cluster_candidate(
    candidate_id: str,
    fingerprint_paths: dict[int, pathlib.Path],
    frame_maps: dict[int, dict[float, dict]],
    coordinate_resolver: Callable[[int, float], dict],
    cutoff_nm: float = 0.20,
) -> dict:
    if cutoff_nm <= 0:
        raise ValueError("RMSD cutoff must be positive")
    frames = []
    for seed in sorted(fingerprint_paths):
        path = pathlib.Path(fingerprint_paths[seed])
        data = np.load(path)
        times = np.asarray(data["time_ps"], dtype=float)
        alignment = np.asarray(data["alignment_A"], dtype=float)
        local = np.asarray(data["local_A"], dtype=float)
        if (
            times.ndim != 1
            or alignment.ndim != 3
            or local.ndim != 3
            or len(times) != len(alignment)
            or len(times) != len(local)
        ):
            raise ValueError(f"{path}: invalid fingerprint arrays")
        for time_ps, frame_alignment, frame_local in zip(times, alignment, local):
            frames.append(
                {
                    "seed": int(seed),
                    "time_ps": float(time_ps),
                    "alignment": frame_alignment,
                    "local": frame_local,
                }
            )
    if not frames:
        return {"reproduced": False, "cluster_count": 0, "medoid": None}
    reference = frames[0]["alignment"]
    aligned = [
        _align_local(reference, frame["alignment"], frame["local"])
        for frame in frames
    ]
    if any(item.shape != aligned[0].shape for item in aligned):
        raise ValueError("local fingerprints have inconsistent shapes")
    distances = np.zeros((len(aligned), len(aligned)), dtype=float)
    adjacency = [[] for _ in aligned]
    cutoff_A = float(cutoff_nm) * 10.0
    for left in range(len(aligned)):
        for right in range(left + 1, len(aligned)):
            rmsd = float(np.sqrt(np.mean((aligned[left] - aligned[right]) ** 2)))
            distances[left, right] = distances[right, left] = rmsd
            if rmsd <= cutoff_A:
                adjacency[left].append(right)
                adjacency[right].append(left)
    qualifying = []
    for component in _components(adjacency):
        membership = Counter(frames[index]["seed"] for index in component)
        contributing = sum(count >= 2 for count in membership.values())
        if len(component) >= 4 and contributing >= 2:
            pairwise = distances[np.ix_(component, component)]
            mean_pairwise = float(pairwise.mean())
            qualifying.append((component, membership, contributing, mean_pairwise))
    if not qualifying:
        return {"reproduced": False, "cluster_count": 0, "medoid": None}
    qualifying.sort(
        key=lambda item: (
            -len(item[0]),
            -item[2],
            item[3],
            [
                (frames[index]["seed"], frames[index]["time_ps"])
                for index in item[0]
            ],
        )
    )
    component, membership, _, _ = qualifying[0]
    medoid_index = min(
        component,
        key=lambda index: (
            float(distances[index, component].mean()),
            frames[index]["seed"],
            frames[index]["time_ps"],
        ),
    )
    source = frames[medoid_index]
    seed = int(source["seed"])
    time_ps = float(source["time_ps"])
    frame = frame_maps.get(seed, {}).get(round(time_ps, 6))
    if frame is None:
        raise ValueError(f"{candidate_id}: medoid lacks matching raw primitive frame")
    coordinate = coordinate_resolver(seed, time_ps)
    sha = str(coordinate.get("source_coordinate_sha256", ""))
    if not len(sha) == 64 or any(char not in "0123456789abcdef" for char in sha):
        raise ValueError(f"{candidate_id}: medoid coordinate checksum is invalid")
    return {
        "reproduced": True,
        "cluster_count": len(qualifying),
        "medoid": {
            "source_replica_seed": seed,
            "time_ps": time_ps,
            **coordinate,
            "reaction_geometry": {
                "distance_nm": float(frame["distance_nm"]),
                "angle_deg": float(frame["angle_deg"]),
            },
            "gate_opening_nm": float(frame["gate_opening_nm"]),
            "proton_route_annotation": str(
                frame.get(
                    "proton_route_annotation",
                    "NOT_EVALUATED_PROTON_ROUTE_ANNOTATION",
                )
            ),
            "membership_by_replica": {
                str(key): int(value)
                for key, value in sorted(membership.items())
                if value >= 2
            },
        },
    }


def _sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _extract_coordinate(
    record: dict,
    time_ps: float,
    output_path: pathlib.Path,
) -> dict:
    import MDAnalysis as mda

    universe = mda.Universe(str(record["tpr"]), str(record["xtc"]))
    nearest = None
    for ts in universe.trajectory:
        delta = abs(float(ts.time) - float(time_ps))
        choice = (delta, int(ts.frame), float(ts.time))
        if nearest is None or choice < nearest:
            nearest = choice
    if nearest is None or nearest[0] > 1.0e-3:
        raise ValueError(
            f"cannot extract exact medoid time {time_ps} ps from {record['xtc']}"
        )
    universe.trajectory[nearest[1]]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        raise ValueError(f"refuse to overwrite medoid coordinate: {output_path}")
    universe.atoms.write(str(output_path))
    return {
        "source_coordinate_sha256": _sha256(output_path),
        "source_coordinate_path": str(output_path),
        "source_tpr": str(record["tpr"]),
        "source_xtc": str(record["xtc"]),
        "source_frame_time_ps": nearest[2],
    }


def _load_json(path) -> dict:
    return json.loads(pathlib.Path(path).read_text())


def assemble_payload(
    authority_path: pathlib.Path,
    stage_a_audits_path: pathlib.Path,
    stage_a_top6_path: pathlib.Path,
    stage_b_records_path: pathlib.Path,
    output_dir: pathlib.Path,
    forbidden_job_ids: list[str],
) -> dict:
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=False)
    authority = _load_json(authority_path)
    candidate_ids = [
        str(item["candidate_id"]) for item in authority["selected"]
    ]
    stage_a_audits = _load_json(stage_a_audits_path)
    if isinstance(stage_a_audits, dict):
        stage_a_audits = stage_a_audits["replica_audits"]
    top6 = _load_json(stage_a_top6_path)
    selected = [str(value) for value in top6["selected_candidate_ids"]]
    record_manifest = _load_json(stage_b_records_path)
    manifest_records = (
        record_manifest["records"]
        if isinstance(record_manifest, dict)
        else record_manifest
    )
    expected_pairs = {
        (candidate, seed)
        for candidate in selected
        for seed in EXPECTED_SEEDS
    }
    actual_pairs = {
        (str(record["candidate_id"]), int(record["velocity_seed"]))
        for record in manifest_records
    }
    if len(actual_pairs) != len(manifest_records) or actual_pairs != expected_pairs:
        raise ValueError("Stage B record manifest is incomplete or duplicated")

    audit_records = []
    fingerprint_paths = {}
    frame_maps = {}
    production_records = {}
    for record in manifest_records:
        candidate = str(record["candidate_id"])
        seed = int(record["velocity_seed"])
        production_records[(candidate, seed)] = record
        completion_status = str(record["completion_status"])
        base = {
            "candidate_id": candidate,
            "velocity_seed": seed,
            "path": str(record["path"]),
            "complete_json": {"status": completion_status},
            "advertised_scientific_status": (
                "NOT_EVALUATED_PENDING_ENSEMBLE_AUDIT"
            ),
        }
        if completion_status != "PASS_TECHNICAL_STAGE_B":
            audit_records.append(base)
            continue
        complete = _load_json(record["complete_json_path"])
        if complete.get("status") != "PASS_TECHNICAL_STAGE_B":
            raise ValueError("Stage B manifest and COMPLETE status disagree")
        primitives = pathlib.Path(record["primitives_dir"])
        frames = load_raw_frames(primitives)
        _add_proton_annotations(primitives, frames)
        frame_maps.setdefault(candidate, {})[seed] = {
            round(float(frame["time_ps"]), 6): frame for frame in frames
        }
        fingerprint_paths.setdefault(candidate, {})[seed] = (
            primitives / "nac_fingerprints.npz"
        )
        replica_audit = _load_json(record["replica_audit_path"])
        audit_records.append(
            {
                **base,
                "complete_json": complete,
                "free_tpr_contract": _load_json(
                    record["free_tpr_contract_path"]
                ),
                "numerical": _load_json(record["numerical_audit_path"]),
                "thermodynamics": replica_audit["thermodynamics"],
                "minimum_contact": _load_json(
                    primitives / "minimum_contact.json"
                ),
                "primitive_frames": frames,
            }
        )

    clusters = {}
    medoid_root = output_dir / "medoids"
    for candidate in selected:
        paths = fingerprint_paths.get(candidate, {})

        def resolve(seed, time_ps, candidate_id=candidate):
            record = production_records[(candidate_id, int(seed))]
            safe = candidate_id.replace("/", "_")
            coordinate_path = (
                medoid_root
                / safe
                / f"seed{int(seed)}_time{float(time_ps):.3f}ps.gro"
            )
            return _extract_coordinate(record, time_ps, coordinate_path)

        clusters[candidate] = cluster_candidate(
            candidate,
            paths,
            frame_maps.get(candidate, {}),
            coordinate_resolver=resolve,
            cutoff_nm=0.20,
        )

    payload = {
        "schema_version": 1,
        "authority": {
            "selected_candidate_ids": candidate_ids,
            "gate_definition": "NylC residues 261-266; Thr267 excluded",
        },
        "stageA": {
            "expected_replica_count": 36,
            "replica_records": stage_a_audits,
            "selected_candidate_ids": selected,
        },
        "stageB": {
            "replica_records": audit_records,
            "cluster_summaries": clusters,
        },
        "forbidden_job_ids": [str(value) for value in forbidden_job_ids],
    }
    (output_dir / "final_audit_input.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n"
    )
    (output_dir / "cluster_summaries.json").write_text(
        json.dumps(clusters, indent=2, sort_keys=True) + "\n"
    )
    return payload


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--authority", type=pathlib.Path, required=True)
    parser.add_argument("--stage-a-audits", type=pathlib.Path, required=True)
    parser.add_argument("--stage-a-top6", type=pathlib.Path, required=True)
    parser.add_argument("--stage-b-records", type=pathlib.Path, required=True)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    parser.add_argument("--forbidden-job-id", action="append", default=[])
    args = parser.parse_args()
    payload = assemble_payload(
        args.authority,
        args.stage_a_audits,
        args.stage_a_top6,
        args.stage_b_records,
        args.output_dir,
        args.forbidden_job_id,
    )
    print(
        json.dumps(
            {
                "candidate_universe": len(
                    payload["authority"]["selected_candidate_ids"]
                ),
                "stageB_replica_records": len(
                    payload["stageB"]["replica_records"]
                ),
                "cluster_candidate_count": len(
                    payload["stageB"]["cluster_summaries"]
                ),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
