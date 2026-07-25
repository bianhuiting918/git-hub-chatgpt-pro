#!/usr/bin/env python3
"""Select twelve deterministic, event-distinct NylC real-NAC source frames."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

from analyze_nac_series import read_xvg


TIME_TOLERANCE_PS = 1e-6
HYDROGEN_RE = re.compile(r"^[0-9]*H", re.IGNORECASE)


def _finish_event(events, members):
    if not members:
        return
    event = {
        "event_id": f"event_{len(events) + 1:04d}",
        "start_ps": members[0]["time_ps"],
        "end_ps": members[-1]["time_ps"],
        "duration_ps": members[-1]["time_ps"] - members[0]["time_ps"],
        "frame_count": len(members),
        "member_times_ps": [member["time_ps"] for member in members],
        "members": [dict(member) for member in members],
        "mean_distance_nm": sum(member["distance_nm"] for member in members) / len(members),
        "mean_angle_deg": sum(member["angle_deg"] for member in members) / len(members),
    }
    potentials = [
        member["potential_energy_kj_mol"]
        for member in members
        if member.get("potential_energy_kj_mol") is not None
    ]
    if potentials:
        event["minimum_potential_energy_kj_mol"] = min(potentials)
    events.append(event)


def group_nac_events(rows: list[dict], sample_interval_ps: float = 2.0) -> list[dict]:
    """Group adjacent NAC samples; missing time points always start a new event."""
    if not math.isfinite(sample_interval_ps) or sample_interval_ps <= 0:
        raise ValueError("sample_interval_ps must be positive and finite")

    events = []
    current = []
    previous_time = None
    for source_index, source_row in enumerate(rows):
        row = dict(source_row)
        time_ps = float(row["time_ps"])
        if previous_time is not None and time_ps <= previous_time + TIME_TOLERANCE_PS:
            raise ValueError("NAC input times must be strictly increasing")
        previous_time = time_ps
        row["time_ps"] = time_ps
        row.setdefault("frame_index", source_index)

        if not row.get("nac", True):
            _finish_event(events, current)
            current = []
            continue

        if (
            current
            and time_ps - current[-1]["time_ps"]
            > sample_interval_ps + TIME_TOLERANCE_PS
        ):
            _finish_event(events, current)
            current = []
        current.append(row)

    _finish_event(events, current)
    return events


def read_gro(path: Path) -> tuple[list[dict], tuple[float, ...]]:
    """Parse fixed-column GRO atoms while preserving sequential global identity."""
    lines = path.read_text().splitlines()
    if len(lines) < 3:
        raise ValueError(f"{path}: incomplete GRO")
    try:
        atom_count = int(lines[1].strip())
    except ValueError as error:
        raise ValueError(f"{path}: invalid GRO atom count") from error
    if len(lines) != atom_count + 3:
        raise ValueError(
            f"{path}: expected {atom_count + 3} lines, observed {len(lines)}"
        )

    atoms = []
    residue_instance = 0
    previous_residue = None
    for global_index1, raw in enumerate(lines[2 : 2 + atom_count], 1):
        if len(raw) < 44:
            raise ValueError(f"{path}: atom line {global_index1} is too short")
        try:
            resid = int(raw[0:5])
            atom_number = int(raw[15:20])
            coord = (float(raw[20:28]), float(raw[28:36]), float(raw[36:44]))
        except ValueError as error:
            raise ValueError(
                f"{path}: invalid fixed-column atom line {global_index1}"
            ) from error
        resname = raw[5:10].strip()
        atom_name = raw[10:15].strip()
        residue_key = (resid, resname)
        if residue_key != previous_residue:
            residue_instance += 1
            previous_residue = residue_key
        atoms.append(
            {
                "global_index1": global_index1,
                "resid": resid,
                "resname": resname,
                "atom_name": atom_name,
                "atom_number": atom_number,
                "residue_instance": residue_instance,
                "coord_nm": coord,
            }
        )
    try:
        box = tuple(float(value) for value in lines[-1].split())
    except ValueError as error:
        raise ValueError(f"{path}: invalid GRO box") from error
    if len(box) not in (3, 9):
        raise ValueError(f"{path}: GRO box must contain 3 or 9 values")
    return atoms, box


def _kabsch_transform(reference: np.ndarray, mobile: np.ndarray):
    reference = np.asarray(reference, dtype=float)
    mobile = np.asarray(mobile, dtype=float)
    if reference.shape != mobile.shape or reference.ndim != 2 or reference.shape[1] != 3:
        raise ValueError("Kabsch inputs must have identical N x 3 shapes")
    if len(reference) < 3:
        raise ValueError("Kabsch alignment requires at least three atoms")
    reference_center = reference.mean(axis=0)
    mobile_center = mobile.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (reference - reference_center)
    left, _, right_t = np.linalg.svd(covariance)
    correction = np.eye(3)
    correction[-1, -1] = np.sign(np.linalg.det(left @ right_t))
    rotation = left @ correction @ right_t
    return rotation, mobile_center, reference_center


def kabsch_rmsd(reference: np.ndarray, mobile: np.ndarray) -> float:
    rotation, mobile_center, reference_center = _kabsch_transform(reference, mobile)
    aligned = (np.asarray(mobile, dtype=float) - mobile_center) @ rotation + reference_center
    delta = aligned - np.asarray(reference, dtype=float)
    return float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))


def _rmsd_index(member):
    return int(member.get("rmsd_index", member["frame_index"]))


def _representative(event, rmsd_nm, predicate=None, exact_time_ps=None):
    members = event["members"]
    event_indices = [_rmsd_index(member) for member in members]
    eligible = []
    for member in members:
        if exact_time_ps is not None and abs(member["time_ps"] - exact_time_ps) > TIME_TOLERANCE_PS:
            continue
        if predicate is not None and not predicate(member):
            continue
        index = _rmsd_index(member)
        mean_local = float(np.mean(rmsd_nm[index, event_indices]))
        potential = member.get("potential_energy_kj_mol")
        potential_key = math.inf if potential is None else float(potential)
        geometry_key = (
            abs(float(member["distance_nm"]) - 0.30),
            abs(float(member["angle_deg"]) - 105.0),
        )
        eligible.append(
            (
                (
                    mean_local,
                    potential_key,
                    geometry_key[0],
                    geometry_key[1],
                    float(member["time_ps"]),
                ),
                member,
                mean_local,
            )
        )
    if not eligible:
        return None
    _, member, mean_local = min(eligible, key=lambda item: item[0])
    return {
        "event_id": event["event_id"],
        "event_start_ps": event["start_ps"],
        "event_end_ps": event["end_ps"],
        "event_frame_count": event["frame_count"],
        "time_ps": float(member["time_ps"]),
        "distance_nm": float(member["distance_nm"]),
        "angle_deg": float(member["angle_deg"]),
        "potential_energy_kj_mol": (
            None
            if member.get("potential_energy_kj_mol") is None
            else float(member["potential_energy_kj_mol"])
        ),
        "source_frame_index0": int(member["frame_index"]),
        "rmsd_index0": _rmsd_index(member),
        "mean_within_event_rmsd_nm": mean_local,
    }


def _candidate_rank(candidate):
    potential = candidate["potential_energy_kj_mol"]
    return (
        candidate["mean_within_event_rmsd_nm"],
        math.inf if potential is None else potential,
        abs(candidate["distance_nm"] - 0.30),
        abs(candidate["angle_deg"] - 105.0),
        candidate["time_ps"],
        candidate["event_id"],
    )


def _deduplicated_pick(candidates, count, selected, rmsd_nm, cutoff_nm):
    ordered = sorted(candidates, key=_candidate_rank)
    picked = []
    reference_indices = [item["rmsd_index0"] for item in selected]
    for candidate in ordered:
        index = candidate["rmsd_index0"]
        if all(rmsd_nm[index, other] >= cutoff_nm - 1e-12 for other in reference_indices):
            picked.append(candidate)
            reference_indices.append(index)
            if len(picked) == count:
                return picked
    # The approved design requires twelve distinct events even if a strict
    # structural de-duplication pass yields fewer than the category allocation.
    for candidate in ordered:
        if candidate in picked:
            continue
        picked.append(candidate)
        if len(picked) == count:
            return picked
    raise ValueError(f"only {len(picked)} candidates available for requested count {count}")


def _candidate_id(event_id, time_ps):
    event_number = int(event_id.rsplit("_", 1)[1])
    if abs(time_ps - round(time_ps)) <= TIME_TOLERANCE_PS:
        time_label = str(int(round(time_ps)))
    else:
        time_label = f"{time_ps:.3f}".rstrip("0").rstrip(".").replace(".", "p")
    return f"nac_evt{event_number:02d}_time{time_label}ps"


def select_twelve(
    events: list[dict],
    rmsd_nm: np.ndarray,
    *,
    old_control_time_ps: float = 1462.0,
    central_distance_max_nm: float = 0.32,
    central_angle_min_deg: float = 100.0,
    central_angle_max_deg: float = 110.0,
    local_rmsd_dedup_nm: float = 0.10,
) -> list[dict]:
    """Apply the frozen 4 + 4 + 3 + 1 event-distinct allocation."""
    rmsd_nm = np.asarray(rmsd_nm, dtype=float)
    if rmsd_nm.ndim != 2 or rmsd_nm.shape[0] != rmsd_nm.shape[1]:
        raise ValueError("rmsd_nm must be a square matrix")
    if not np.all(np.isfinite(rmsd_nm)):
        raise ValueError("rmsd_nm contains non-finite values")
    if not np.allclose(rmsd_nm, rmsd_nm.T, atol=1e-10):
        raise ValueError("rmsd_nm must be symmetric")

    legacy_event = None
    for event in events:
        if any(
            abs(member["time_ps"] - old_control_time_ps) <= TIME_TOLERANCE_PS
            for member in event["members"]
        ):
            legacy_event = event
            break
    if legacy_event is None:
        raise ValueError(f"legacy control time {old_control_time_ps} ps is absent")
    legacy = _representative(
        legacy_event, rmsd_nm, exact_time_ps=old_control_time_ps
    )
    legacy["category"] = "legacy_control"

    reserved_event_ids = {legacy_event["event_id"]}
    selected = []

    three_pool = []
    for event in events:
        if event["event_id"] in reserved_event_ids or event["frame_count"] < 3:
            continue
        candidate = _representative(event, rmsd_nm)
        candidate["category"] = "three_frame_medoid"
        three_pool.append(candidate)
    chosen = _deduplicated_pick(
        three_pool, 4, selected, rmsd_nm, local_rmsd_dedup_nm
    )
    selected.extend(chosen)
    reserved_event_ids.update(item["event_id"] for item in chosen)

    two_pool = []
    for event in events:
        if event["event_id"] in reserved_event_ids or event["frame_count"] != 2:
            continue
        candidate = _representative(event, rmsd_nm)
        candidate["category"] = "two_frame_medoid"
        two_pool.append(candidate)
    chosen = _deduplicated_pick(
        two_pool, 4, selected, rmsd_nm, local_rmsd_dedup_nm
    )
    selected.extend(chosen)
    reserved_event_ids.update(item["event_id"] for item in chosen)

    def is_central(member):
        return (
            member["distance_nm"] <= central_distance_max_nm
            and central_angle_min_deg <= member["angle_deg"] <= central_angle_max_deg
        )

    central_pool = []
    for event in events:
        if event["event_id"] in reserved_event_ids:
            continue
        candidate = _representative(event, rmsd_nm, predicate=is_central)
        if candidate is None:
            continue
        candidate["category"] = "central_low_energy"
        central_pool.append(candidate)
    chosen = _deduplicated_pick(
        central_pool, 3, selected, rmsd_nm, local_rmsd_dedup_nm
    )
    selected.extend(chosen)
    reserved_event_ids.update(item["event_id"] for item in chosen)

    selected.append(legacy)
    if len(selected) != 12 or len({item["event_id"] for item in selected}) != 12:
        raise AssertionError("selection did not produce twelve event-distinct candidates")
    for candidate in selected:
        candidate["candidate_id"] = _candidate_id(
            candidate["event_id"], candidate["time_ps"]
        )
    return selected


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _source_paths(authority):
    source_root = Path(authority["source_run_root"])
    selection_root = Path(authority["selection_root"])
    return {
        "tpr": source_root / authority["source_tpr"],
        "xtc": source_root / authority["source_xtc"],
        "edr": source_root / authority["source_edr"],
        "distance_xvg": selection_root / "nac_distance.xvg",
        "angle_xvg": selection_root / "nac_angle.xvg",
        "potential_xvg": selection_root / "potential_energy.xvg",
        "series_audit": selection_root / "nac_energy_series.json",
        "identity_gro": selection_root / "source_selected_NAC.gro",
    }


def _load_source_series(authority, paths):
    distance_rows = read_xvg(paths["distance_xvg"])
    angle_rows = read_xvg(paths["angle_xvg"])
    potential_rows = read_xvg(paths["potential_xvg"])
    if not (
        len(distance_rows) == len(angle_rows) == len(potential_rows)
        == authority["expected_frame_count"]
    ):
        raise ValueError("source series frame count disagrees with authority")

    rows = []
    gate = authority["nac_gate"]
    for frame_index, (distance_row, angle_row, potential_row) in enumerate(
        zip(distance_rows, angle_rows, potential_rows)
    ):
        distance_time, distance = distance_row
        angle_time, angle = angle_row
        potential_time, potential = potential_row
        if (
            abs(distance_time - angle_time) > TIME_TOLERANCE_PS
            or abs(distance_time - potential_time) > TIME_TOLERANCE_PS
        ):
            raise ValueError("source distance/angle/potential times do not align")
        rows.append(
            {
                "frame_index": frame_index,
                "time_ps": distance_time,
                "distance_nm": distance,
                "angle_deg": angle,
                "potential_energy_kj_mol": potential,
                "nac": (
                    distance <= gate["distance_max_nm"]
                    and gate["angle_min_deg"] <= angle <= gate["angle_max_deg"]
                ),
            }
        )
    expected_start, expected_end = authority["analysis_window_ps"]
    if (
        abs(rows[0]["time_ps"] - expected_start) > TIME_TOLERANCE_PS
        or abs(rows[-1]["time_ps"] - expected_end) > TIME_TOLERANCE_PS
    ):
        raise ValueError("source series analysis window disagrees with authority")
    return rows


def _validate_reactive_identities(authority, identity_gro):
    atoms, _ = read_gro(identity_gro)
    expected = {
        "thr267_og1": (267, "THR", "OG1"),
        "l2_c": (356, "L2", "C12"),
        "l2_o": (356, "L2", "O2"),
        "l2_n": (356, "L2", "N3"),
    }
    observed = {}
    for label, global_index1 in authority["reactive_global_index1_m0"].items():
        atom = atoms[int(global_index1) - 1]
        identity = (atom["resid"], atom["resname"], atom["atom_name"])
        if identity != expected[label]:
            raise ValueError(
                f"reactive identity mismatch for {label}: expected {expected[label]}, observed {identity}"
            )
        observed[label] = {
            "global_index1": int(global_index1),
            "resid": atom["resid"],
            "resname": atom["resname"],
            "atom_name": atom["atom_name"],
        }
    return observed


def validate_authority(authority_path: Path):
    authority = json.loads(authority_path.read_text())
    paths = _source_paths(authority)
    for label, path in paths.items():
        if not path.is_file() or path.stat().st_size == 0:
            raise FileNotFoundError(f"missing authority input {label}: {path}")

    if authority["gate_residues"] != [261, 262, 263, 264, 265, 266]:
        raise ValueError("NylC gate residues must be exactly 261-266")
    if authority["excluded_gate_residues"] != [267]:
        raise ValueError("Thr267 must be excluded from the gate group")

    source_audit = json.loads(paths["series_audit"].read_text())
    for key, expected in (
        ("frame_count", authority["expected_frame_count"]),
        ("nac_frame_count", authority["expected_nac_frame_count"]),
        ("nac_run_count", authority["expected_event_count"]),
    ):
        if source_audit.get(key) != expected:
            raise ValueError(f"source audit {key} disagrees with authority")

    rows = _load_source_series(authority, paths)
    events = group_nac_events(rows, authority["sample_interval_ps"])
    nac_frame_count = sum(event["frame_count"] for event in events)
    if nac_frame_count != authority["expected_nac_frame_count"]:
        raise ValueError("recomputed NAC frame count disagrees with authority")
    if len(events) != authority["expected_event_count"]:
        raise ValueError("recomputed NAC event count disagrees with authority")
    observed_atoms = _validate_reactive_identities(authority, paths["identity_gro"])
    return authority, paths, rows, events, observed_atoms


def _gro_time(path):
    title = path.read_text().splitlines()[0]
    matches = re.findall(r"(?:t\s*=\s*|time\s+)([-+0-9.eE]+)", title)
    return None if not matches else float(matches[-1])


def _write_frame_index(path, frame_indices):
    lines = ["[ frames ]"]
    for start in range(0, len(frame_indices), 15):
        lines.append(" ".join(str(index + 1) for index in frame_indices[start : start + 15]))
    path.write_text("\n".join(lines) + "\n")


def _extract_nac_frames(gmx, tpr, xtc, events, frames_dir):
    members = [member for event in events for member in event["members"]]
    for rmsd_index, member in enumerate(members):
        member["rmsd_index"] = rmsd_index
    frame_indices = [member["frame_index"] for member in members]
    frames_dir.mkdir()
    frame_ndx = frames_dir / "source_frame_indices.ndx"
    _write_frame_index(frame_ndx, frame_indices)
    command = [
        gmx,
        "trjconv",
        "-f",
        str(xtc),
        "-s",
        str(tpr),
        "-fr",
        str(frame_ndx),
        "-sep",
        "-o",
        str(frames_dir / "nac.gro"),
    ]
    result = subprocess.run(command, input="0\n", text=True, capture_output=True)
    (frames_dir / "trjconv.stdout").write_text(result.stdout)
    (frames_dir / "trjconv.stderr").write_text(result.stderr)
    if result.returncode != 0:
        raise RuntimeError(f"gmx trjconv failed with exit code {result.returncode}")

    gro_paths = list(frames_dir.glob("nac*.gro"))
    if len(gro_paths) != len(members):
        raise ValueError(
            f"expected {len(members)} extracted NAC GROs, observed {len(gro_paths)}"
        )
    by_time = {}
    for gro_path in gro_paths:
        time_ps = _gro_time(gro_path)
        if time_ps is None:
            raise ValueError(f"cannot recover trajectory time from {gro_path}")
        key = round(time_ps, 6)
        if key in by_time:
            raise ValueError(f"duplicate extracted time {time_ps} ps")
        by_time[key] = gro_path
    for member in members:
        key = round(member["time_ps"], 6)
        if key not in by_time:
            raise ValueError(f"missing extracted time {member['time_ps']} ps")
        member["extracted_gro"] = str(by_time[key])
    return members


def _is_hydrogen(atom_name):
    return bool(HYDROGEN_RE.match(atom_name.strip()))


def _box_matrix(box):
    if len(box) == 3:
        return np.diag(np.asarray(box, dtype=float))
    return np.asarray(
        [
            [box[0], box[3], box[4]],
            [box[5], box[1], box[6]],
            [box[7], box[8], box[2]],
        ],
        dtype=float,
    )


def _minimum_distances(points, targets, box):
    delta = points[:, None, :] - targets[None, :, :]
    matrix = _box_matrix(box)
    fractional = delta @ np.linalg.inv(matrix)
    fractional -= np.round(fractional)
    delta = fractional @ matrix
    return np.sqrt(np.sum(delta * delta, axis=2)).min(axis=1)


def _local_rmsd_matrix(events, authority):
    members = [member for event in events for member in event["members"]]
    first_atoms, _ = read_gro(Path(members[0]["extracted_gro"]))
    protein_start, protein_end = authority["protein_global_index1_range_m0"]
    active_start, active_end = authority["active_chain_global_index1_range_m0"]

    protein_heavy_idx0 = np.asarray(
        [
            index
            for index, atom in enumerate(first_atoms)
            if protein_start <= atom["global_index1"] <= protein_end
            and not _is_hydrogen(atom["atom_name"])
        ],
        dtype=int,
    )
    ligand_heavy_idx0 = np.asarray(
        [
            index
            for index, atom in enumerate(first_atoms)
            if atom["resname"] == authority["ligand_resname"]
            and not _is_hydrogen(atom["atom_name"])
        ],
        dtype=int,
    )
    backbone_idx0 = np.asarray(
        [
            index
            for index, atom in enumerate(first_atoms)
            if active_start <= atom["global_index1"] <= active_end
            and atom["atom_name"] in {"N", "CA", "C"}
        ],
        dtype=int,
    )
    if len(ligand_heavy_idx0) == 0 or len(backbone_idx0) < 3:
        raise ValueError("failed to identify ligand heavy atoms or active-chain backbone")

    protein_residue_instances = np.asarray(
        [first_atoms[index]["residue_instance"] for index in protein_heavy_idx0],
        dtype=int,
    )
    special_instances = {
        atom["residue_instance"]
        for atom in first_atoms
        if active_start <= atom["global_index1"] <= active_end
        and atom["resid"] in {267, 306, 308}
    }

    frames = []
    union_nearby_instances = set()
    reference_identity = [
        (atom["resid"], atom["resname"], atom["atom_name"]) for atom in first_atoms
    ]
    for member in members:
        gro_path = Path(member["extracted_gro"])
        atoms, box = read_gro(gro_path)
        identity = [(atom["resid"], atom["resname"], atom["atom_name"]) for atom in atoms]
        if identity != reference_identity:
            raise ValueError(f"atom identity/order differs in {gro_path}")
        coords = np.asarray([atom["coord_nm"] for atom in atoms], dtype=float)
        protein_coords = coords[protein_heavy_idx0]
        ligand_coords = coords[ligand_heavy_idx0]
        near = _minimum_distances(protein_coords, ligand_coords, box)
        union_nearby_instances.update(
            int(instance)
            for instance, distance in zip(protein_residue_instances, near)
            if distance <= authority["local_environment_cutoff_nm"] + 1e-12
        )
        frames.append(
            {
                "protein": protein_coords,
                "ligand": ligand_coords,
                "backbone": coords[backbone_idx0],
            }
        )

    local_instances = union_nearby_instances | special_instances
    local_protein_mask = np.asarray(
        [instance in local_instances for instance in protein_residue_instances],
        dtype=bool,
    )
    reference_backbone = frames[0]["backbone"]
    aligned_local = []
    for frame in frames:
        rotation, mobile_center, reference_center = _kabsch_transform(
            reference_backbone, frame["backbone"]
        )
        local = np.vstack((frame["protein"][local_protein_mask], frame["ligand"]))
        aligned = (local - mobile_center) @ rotation + reference_center
        aligned_local.append(aligned)

    count = len(aligned_local)
    matrix = np.zeros((count, count), dtype=float)
    for left in range(count):
        for right in range(left + 1, count):
            delta = aligned_local[left] - aligned_local[right]
            value = float(np.sqrt(np.mean(np.sum(delta * delta, axis=1))))
            matrix[left, right] = value
            matrix[right, left] = value

    residue_manifest = []
    for instance in sorted(local_instances):
        atom = next(
            atom
            for atom in first_atoms
            if atom["residue_instance"] == instance
            and protein_start <= atom["global_index1"] <= protein_end
        )
        residue_manifest.append(
            {
                "residue_instance": instance,
                "resid": atom["resid"],
                "resname": atom["resname"],
            }
        )
    metadata = {
        "alignment_atom_count": int(len(backbone_idx0)),
        "local_protein_heavy_atom_count": int(local_protein_mask.sum()),
        "ligand_heavy_atom_count": int(len(ligand_heavy_idx0)),
        "local_total_heavy_atom_count": int(
            local_protein_mask.sum() + len(ligand_heavy_idx0)
        ),
        "local_environment_cutoff_nm": authority["local_environment_cutoff_nm"],
        "local_residues": residue_manifest,
    }
    return matrix, metadata


def run_selection(authority_path: Path, output_root: Path, gmx: str):
    authority, paths, rows, events, observed_atoms = validate_authority(authority_path)
    if output_root.exists():
        raise FileExistsError(f"refusing to overwrite existing output root: {output_root}")
    output_root.mkdir(parents=True)

    members = _extract_nac_frames(
        gmx, paths["tpr"], paths["xtc"], events, output_root / "frames_all_nac"
    )
    rmsd_nm, local_metadata = _local_rmsd_matrix(events, authority)
    np.save(output_root / "local_rmsd_matrix_nm.npy", rmsd_nm)
    (output_root / "local_rmsd_matrix_nm.json").write_text(
        json.dumps(rmsd_nm.tolist(), separators=(",", ":")) + "\n"
    )

    central = authority["central_gate"]
    selected = select_twelve(
        events,
        rmsd_nm,
        old_control_time_ps=authority["old_control_time_ps"],
        central_distance_max_nm=central["distance_max_nm"],
        central_angle_min_deg=central["angle_min_deg"],
        central_angle_max_deg=central["angle_max_deg"],
        local_rmsd_dedup_nm=authority["local_rmsd_dedup_nm"],
    )

    frame_by_rmsd_index = {member["rmsd_index"]: member for member in members}
    frame_records = []
    for member in members:
        gro_path = Path(member["extracted_gro"])
        frame_records.append(
            {
                "rmsd_index0": member["rmsd_index"],
                "source_frame_index0": member["frame_index"],
                "time_ps": member["time_ps"],
                "distance_nm": member["distance_nm"],
                "angle_deg": member["angle_deg"],
                "potential_energy_kj_mol": member["potential_energy_kj_mol"],
                "gro_path": str(gro_path),
                "gro_sha256": sha256_file(gro_path),
            }
        )
    (output_root / "all_nac_frames.json").write_text(
        json.dumps({"schema_version": 1, "frames": frame_records}, indent=2) + "\n"
    )

    candidate_dir = output_root / "candidates"
    candidate_dir.mkdir()
    for candidate in selected:
        source_member = frame_by_rmsd_index[candidate["rmsd_index0"]]
        source_gro = Path(source_member["extracted_gro"])
        destination = candidate_dir / f"{candidate['candidate_id']}.gro"
        shutil.copy2(source_gro, destination)
        candidate["source_gro"] = str(destination)
        candidate["source_gro_sha256"] = sha256_file(destination)

    source_hashes = {
        label: {"path": str(paths[label]), "sha256": sha256_file(paths[label])}
        for label in (
            "tpr",
            "xtc",
            "edr",
            "distance_xvg",
            "angle_xvg",
            "potential_xvg",
            "series_audit",
        )
    }
    output = {
        "schema_version": 1,
        "authority_manifest": str(authority_path),
        "authority_manifest_sha256": sha256_file(authority_path),
        "candidate_id": authority["candidate_id"],
        "source_job_id": authority["source_job_id"],
        "selection_audit_job_id": authority["selection_audit_job_id"],
        "source_hashes": source_hashes,
        "reactive_atom_identities": observed_atoms,
        "gate_definition": "NylC residues 261-266; Thr267 excluded",
        "analysis_window_ps": authority["analysis_window_ps"],
        "sample_interval_ps": authority["sample_interval_ps"],
        "frame_count": len(rows),
        "nac_frame_count": len(members),
        "event_count": len(events),
        "local_rmsd": local_metadata,
        "local_rmsd_dedup_nm": authority["local_rmsd_dedup_nm"],
        "selection_allocation": authority["selection_allocation"],
        "selected": selected,
        "scientific_status": "SOURCE_SELECTION_ONLY_NOT_MD_VALIDATED",
    }
    (output_root / "ensemble_selection.json").write_text(
        json.dumps(output, indent=2, sort_keys=True) + "\n"
    )
    return output


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--authority", type=Path, required=True)
    parser.add_argument("--output-root", type=Path)
    parser.add_argument("--gmx", default="gmx")
    parser.add_argument("--validate-only", action="store_true")
    args = parser.parse_args()

    try:
        authority, paths, rows, events, observed_atoms = validate_authority(args.authority)
        if args.validate_only:
            print(
                json.dumps(
                    {
                        "status": "PASS_VALIDATE_ONLY",
                        "source_paths": {
                            label: str(paths[label]) for label in ("tpr", "xtc", "edr")
                        },
                        "frame_count": len(rows),
                        "nac_frame_count": sum(
                            event["frame_count"] for event in events
                        ),
                        "event_count": len(events),
                        "reactive_atom_identities": observed_atoms,
                    },
                    indent=2,
                    sort_keys=True,
                )
            )
            return 0
        if args.output_root is None:
            raise ValueError("--output-root is required unless --validate-only is used")
        result = run_selection(args.authority, args.output_root, args.gmx)
    except (FileExistsError, FileNotFoundError, RuntimeError, ValueError) as error:
        print(str(error), file=sys.stderr)
        return 2
    print(
        json.dumps(
            {
                "status": "PASS_SELECTION",
                "output_root": str(args.output_root),
                "selected_count": len(result["selected"]),
            },
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
