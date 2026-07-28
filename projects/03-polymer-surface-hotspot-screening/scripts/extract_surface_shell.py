#!/usr/bin/env python3
"""Extract a compact exterior-connected A/C/OA/HD affinity shell."""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
from scipy import ndimage
from scipy.spatial import cKDTree


CHANNELS = ("A", "C", "OA", "HD")
VDW_RADII = {
    "H": 1.20,
    "HD": 1.20,
    "HS": 1.20,
    "C": 1.70,
    "A": 1.70,
    "N": 1.55,
    "NA": 1.55,
    "NS": 1.55,
    "O": 1.52,
    "OA": 1.52,
    "OS": 1.52,
    "F": 1.47,
    "P": 1.80,
    "S": 1.80,
    "SA": 1.80,
    "CL": 1.75,
    "BR": 1.85,
    "I": 1.98,
    "MG": 1.73,
    "MN": 1.73,
    "ZN": 1.39,
    "CA": 1.94,
    "FE": 1.56,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def connectivity_structure(connectivity: int) -> np.ndarray:
    if connectivity == 6:
        return ndimage.generate_binary_structure(3, 1)
    if connectivity == 18:
        return ndimage.generate_binary_structure(3, 2)
    if connectivity == 26:
        return ndimage.generate_binary_structure(3, 3)
    raise ValueError("connectivity must be 6, 18, or 26")


def exterior_connected_solvent(
    occupied: np.ndarray, connectivity: int = 26
) -> np.ndarray:
    occupied = np.asarray(occupied, dtype=bool)
    if occupied.ndim != 3 or min(occupied.shape) < 2:
        raise ValueError("occupied must be a three-dimensional array")
    solvent = ~occupied
    seed = np.zeros_like(solvent)
    seed[0, :, :] = solvent[0, :, :]
    seed[-1, :, :] = solvent[-1, :, :]
    seed[:, 0, :] = solvent[:, 0, :]
    seed[:, -1, :] = solvent[:, -1, :]
    seed[:, :, 0] = solvent[:, :, 0]
    seed[:, :, -1] = solvent[:, :, -1]
    return ndimage.binary_propagation(
        seed, structure=connectivity_structure(connectivity), mask=solvent
    )


def exterior_shell_mask(
    occupied: np.ndarray,
    spacing: float,
    inner: float,
    outer: float,
    connectivity: int = 26,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    if spacing <= 0 or inner < 0 or outer <= inner:
        raise ValueError("invalid shell geometry")
    occupied = np.asarray(occupied, dtype=bool)
    exterior = exterior_connected_solvent(occupied, connectivity)
    distances = ndimage.distance_transform_edt(~occupied, sampling=float(spacing))
    mask = exterior & (distances > float(inner)) & (distances <= float(outer))
    return mask, distances, exterior


def shell_area_proxy(
    n_shell_points: int, spacing: float, inner: float, outer: float
) -> float:
    thickness = float(outer) - float(inner)
    if thickness <= 0:
        raise ValueError("shell thickness must be positive")
    return float(n_shell_points) * float(spacing) ** 3 / thickness


def steiner_surface_area_from_counts(
    cumulative_counts: list[int] | tuple[int, int, int] | np.ndarray,
    spacing: float,
) -> float:
    """Estimate the zero-offset area from cumulative h, 2h, and 3h shell volumes."""
    counts = np.asarray(cumulative_counts, dtype=float)
    if counts.shape != (3,) or np.any(counts < 0):
        raise ValueError("three nonnegative cumulative shell counts are required")
    if spacing <= 0:
        raise ValueError("spacing must be positive")
    if np.any(np.diff(counts) < 0):
        raise ValueError("cumulative shell counts must be nondecreasing")
    volumes = counts * float(spacing) ** 3
    estimate = (
        18.0 * volumes[0] - 9.0 * volumes[1] + 2.0 * volumes[2]
    ) / (6.0 * float(spacing))
    if estimate <= 0:
        raise ValueError("nonpositive Steiner surface-area estimate")
    return float(estimate)


def parse_autogrid_map(path: Path) -> dict:
    spacing = None
    n_elements = None
    center = None
    with path.open(encoding="ascii", errors="strict") as handle:
        header = [handle.readline() for _ in range(6)]
        for line in header:
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "SPACING":
                spacing = float(fields[1])
            elif fields[0] == "NELEMENTS":
                n_elements = np.array([int(value) for value in fields[1:4]], dtype=int)
            elif fields[0] == "CENTER":
                center = np.array([float(value) for value in fields[1:4]], dtype=float)
        values = np.fromiter(
            (float(line.strip()) for line in handle if line.strip()),
            dtype=np.float32,
        )
    if spacing is None or n_elements is None or center is None:
        raise ValueError(f"incomplete AutoGrid header: {path}")
    shape = n_elements + 1
    expected = int(np.prod(shape))
    if values.size != expected:
        raise ValueError(
            f"AutoGrid value count mismatch for {path}: {values.size} != {expected}"
        )
    values = values.reshape((shape[2], shape[1], shape[0])).transpose(2, 1, 0)
    origin = center - n_elements.astype(float) * spacing / 2.0
    if not np.isfinite(values).all():
        raise ValueError(f"nonfinite AutoGrid values: {path}")
    return {
        "spacing": float(spacing),
        "n_elements": n_elements,
        "center": center,
        "origin": origin,
        "values": values,
    }


def atom_radius(atom_type: str) -> float:
    key = atom_type.strip().upper()
    if key in VDW_RADII:
        return VDW_RADII[key]
    letters = "".join(character for character in key if character.isalpha())
    if letters[:2] in VDW_RADII:
        return VDW_RADII[letters[:2]]
    if letters[:1] in VDW_RADII:
        return VDW_RADII[letters[:1]]
    raise ValueError(f"unsupported AutoDock atom type: {atom_type}")


def parse_pdbqt_atoms(path: Path) -> dict:
    coordinates = []
    radii = []
    labels = []
    with path.open(encoding="ascii", errors="strict") as handle:
        for line in handle:
            if line[:6] not in {"ATOM  ", "HETATM"}:
                continue
            try:
                coordinate = [float(line[30:38]), float(line[38:46]), float(line[46:54])]
                atom_name = line[12:16].strip()
                residue_name = line[17:20].strip() or "UNK"
                chain = line[21:22].strip() or "_"
                residue_number = line[22:27].strip()
                atom_type = line.split()[-1]
            except (ValueError, IndexError) as error:
                raise ValueError(f"invalid PDBQT atom line: {line.rstrip()}") from error
            coordinates.append(coordinate)
            radii.append(atom_radius(atom_type))
            labels.append(f"{chain}:{residue_name}:{residue_number}:{atom_name}")
    if not coordinates:
        raise ValueError(f"no receptor atoms in {path}")
    return {
        "coordinates": np.asarray(coordinates, dtype=float),
        "radii": np.asarray(radii, dtype=float),
        "labels": np.asarray(labels),
    }


def build_occupancy(
    shape: tuple[int, int, int],
    origin: np.ndarray,
    spacing: float,
    atom_coords: np.ndarray,
    atom_radii: np.ndarray,
    probe_radius: float = 1.4,
) -> np.ndarray:
    shape_array = np.asarray(shape, dtype=int)
    if shape_array.shape != (3,) or np.any(shape_array <= 0):
        raise ValueError("invalid grid shape")
    origin = np.asarray(origin, dtype=float)
    atom_coords = np.asarray(atom_coords, dtype=float)
    atom_radii = np.asarray(atom_radii, dtype=float)
    if atom_coords.ndim != 2 or atom_coords.shape[1] != 3:
        raise ValueError("atom coordinates must be N by 3")
    if atom_coords.shape[0] != atom_radii.shape[0]:
        raise ValueError("atom coordinate and radius counts differ")
    if probe_radius < 0:
        raise ValueError("probe radius must be nonnegative")

    occupied = np.zeros(tuple(shape_array), dtype=bool)
    for coordinate, radius in zip(atom_coords, atom_radii):
        inflated = float(radius) + float(probe_radius)
        lower = np.maximum(
            0, np.floor((coordinate - inflated - origin) / spacing).astype(int)
        )
        upper = np.minimum(
            shape_array - 1,
            np.ceil((coordinate + inflated - origin) / spacing).astype(int),
        )
        if np.any(lower > upper):
            continue
        axes = [
            origin[axis]
            + np.arange(lower[axis], upper[axis] + 1, dtype=float) * spacing
            for axis in range(3)
        ]
        x, y, z = np.meshgrid(*axes, indexing="ij")
        local = (
            (x - coordinate[0]) ** 2
            + (y - coordinate[1]) ** 2
            + (z - coordinate[2]) ** 2
        ) <= inflated**2
        slices = tuple(
            slice(lower[axis], upper[axis] + 1) for axis in range(3)
        )
        occupied[slices] |= local
    return occupied


def _neighbor_offsets(connectivity: int = 26) -> list[tuple[int, int, int]]:
    structure = connectivity_structure(connectivity)
    center = np.array(structure.shape) // 2
    offsets = []
    for index in np.argwhere(structure):
        offset = tuple((index - center).tolist())
        if offset != (0, 0, 0):
            offsets.append(offset)
    return offsets


def coordinate_graph(
    coordinates: np.ndarray, spacing: float, connectivity: int = 26
) -> tuple[np.ndarray, np.ndarray]:
    """Build a symmetric distance graph without absolute-coordinate rounding."""
    coordinates = np.asarray(coordinates, dtype=float)
    if coordinates.size == 0:
        return np.array([0], dtype=np.int64), np.array([], dtype=np.int32)
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("coordinates must be N by 3")
    if spacing <= 0:
        raise ValueError("spacing must be positive")
    radius_factors = {6: 1.0, 18: np.sqrt(2.0), 26: np.sqrt(3.0)}
    if connectivity not in radius_factors:
        raise ValueError("connectivity must be 6, 18, or 26")
    tolerance = max(1e-6, abs(float(spacing)) * 1e-6)
    radius = radius_factors[connectivity] * float(spacing) + tolerance
    pairs = cKDTree(coordinates).query_pairs(radius, output_type="ndarray")
    if pairs.size == 0:
        return (
            np.zeros(coordinates.shape[0] + 1, dtype=np.int64),
            np.array([], dtype=np.int32),
        )
    pairs = np.asarray(pairs, dtype=np.int64).reshape(-1, 2)
    rows = np.concatenate((pairs[:, 0], pairs[:, 1]))
    columns = np.concatenate((pairs[:, 1], pairs[:, 0])).astype(np.int32)
    order = np.lexsort((columns, rows))
    rows = rows[order]
    columns = columns[order]
    counts = np.bincount(rows, minlength=coordinates.shape[0])
    indptr = np.concatenate(([0], np.cumsum(counts))).astype(np.int64)
    return indptr, columns


def _residue_from_atom_label(label: str) -> str:
    fields = str(label).split(":")
    return ":".join(fields[:3]) if len(fields) >= 3 else str(label)


def extract_shell_from_arrays(
    occupied: np.ndarray,
    channel_arrays: dict[str, np.ndarray],
    origin: np.ndarray,
    spacing: float,
    shell_inner: float,
    shell_outer: float,
    atom_coords: np.ndarray,
    atom_radii: np.ndarray,
    atom_labels: np.ndarray,
    tile_id: str,
) -> dict:
    occupied = np.asarray(occupied, dtype=bool)
    for channel in CHANNELS:
        if channel not in channel_arrays:
            raise ValueError(f"missing channel {channel}")
        if np.asarray(channel_arrays[channel]).shape != occupied.shape:
            raise ValueError(f"channel geometry mismatch: {channel}")
    atom_coords = np.asarray(atom_coords, dtype=float)
    atom_radii = np.asarray(atom_radii, dtype=float)
    atom_labels = np.asarray(atom_labels)
    if not (len(atom_coords) == len(atom_radii) == len(atom_labels)):
        raise ValueError("atom metadata lengths differ")

    mask, _, _ = exterior_shell_mask(
        occupied, spacing, shell_inner, shell_outer, connectivity=26
    )
    grid_indices = np.argwhere(mask)
    coordinates = np.asarray(origin, dtype=float) + grid_indices * float(spacing)
    if coordinates.size:
        nearest_atom = cKDTree(atom_coords).query(coordinates, k=1)[1].astype(np.int32)
        nearest_residue = np.asarray(
            [_residue_from_atom_label(atom_labels[index]) for index in nearest_atom]
        )
    else:
        nearest_atom = np.array([], dtype=np.int32)
        nearest_residue = np.array([], dtype="<U1")
    indptr, neighbor_indices = coordinate_graph(coordinates, spacing, connectivity=26)
    result = {
        "coordinates": coordinates.astype(np.float32),
        "nearest_atom_index": nearest_atom,
        "nearest_residue": nearest_residue,
        "tile_provenance": np.full(coordinates.shape[0], str(tile_id)),
        "neighbor_indptr": indptr,
        "neighbor_indices": neighbor_indices,
        "summary": {
            "n_shell_points": int(coordinates.shape[0]),
            "shell_area_proxy": shell_area_proxy(
                coordinates.shape[0], spacing, shell_inner, shell_outer
            ),
            "shell_inner": float(shell_inner),
            "shell_outer": float(shell_outer),
            "spacing": float(spacing),
        },
    }
    selection = tuple(grid_indices[:, axis] for axis in range(3))
    for channel in CHANNELS:
        result[channel] = np.asarray(channel_arrays[channel])[selection].astype(np.float32)
    return result


def _coordinate_key(coordinate: np.ndarray, spacing: float) -> tuple[int, int, int]:
    return tuple(np.rint(np.asarray(coordinate) / spacing * 1_000_000).astype(np.int64))


def merge_shell_records(
    records: list[dict],
    spacing: float,
    value_tolerance: float = 1e-5,
    phase_match_fraction: float = 0.5,
    phase_p999_tolerance: float = 0.10,
    phase_fraction_over_tolerance: float = 0.001,
    phase_max_tolerance: float = 0.25,
    atom_coords: np.ndarray | None = None,
    atom_labels: np.ndarray | None = None,
) -> dict:
    """Merge phase-shifted tile shells using cross-tile mutual nearest neighbors."""
    if not records:
        raise ValueError("no shell records to merge")
    if spacing <= 0 or not 0 < phase_match_fraction <= 1:
        raise ValueError("invalid tile merge geometry")

    lengths = []
    coordinates_by_record = []
    for record in records:
        coordinates = np.asarray(record["coordinates"], dtype=float)
        if coordinates.ndim != 2 or coordinates.shape[1] != 3:
            raise ValueError("tile coordinates must be N by 3")
        n_points = len(coordinates)
        for key in (
            "nearest_atom_index",
            "nearest_residue",
            "tile_provenance",
            *CHANNELS,
        ):
            if len(record[key]) != n_points:
                raise ValueError(f"tile shell length mismatch: {key}")
        lengths.append(n_points)
        coordinates_by_record.append(coordinates)

    offsets = np.concatenate(([0], np.cumsum(lengths))).astype(np.int64)
    total_points = int(offsets[-1])
    parents = np.arange(total_points, dtype=np.int64)
    ranks = np.zeros(total_points, dtype=np.int8)

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = int(parents[index])
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root == right_root:
            return
        if ranks[left_root] < ranks[right_root]:
            left_root, right_root = right_root, left_root
        parents[right_root] = left_root
        if ranks[left_root] == ranks[right_root]:
            ranks[left_root] += 1

    phase_differences = {channel: [] for channel in CHANNELS}
    n_exact_pairs = 0
    n_phase_pairs = 0
    maximum_pair_distance = 0.0
    match_tolerance = phase_match_fraction * float(spacing)
    exact_coordinate_tolerance = max(1e-7, abs(float(spacing)) * 1e-7)

    for left_record_index in range(len(records)):
        left_coordinates = coordinates_by_record[left_record_index]
        if len(left_coordinates) == 0:
            continue
        left_tree = cKDTree(left_coordinates)
        for right_record_index in range(left_record_index + 1, len(records)):
            right_coordinates = coordinates_by_record[right_record_index]
            if len(right_coordinates) == 0:
                continue
            right_tree = cKDTree(right_coordinates)
            left_distance, left_neighbor = right_tree.query(
                left_coordinates,
                k=1,
                distance_upper_bound=match_tolerance,
            )
            right_distance, right_neighbor = left_tree.query(
                right_coordinates,
                k=1,
                distance_upper_bound=match_tolerance,
            )
            for left_local_index, distance in enumerate(left_distance):
                if not np.isfinite(distance):
                    continue
                right_local_index = int(left_neighbor[left_local_index])
                if right_local_index >= len(right_coordinates):
                    continue
                if (
                    int(right_neighbor[right_local_index]) != left_local_index
                    or not np.isfinite(right_distance[right_local_index])
                ):
                    continue
                left_global = int(offsets[left_record_index] + left_local_index)
                right_global = int(offsets[right_record_index] + right_local_index)
                maximum_pair_distance = max(maximum_pair_distance, float(distance))
                if float(distance) <= exact_coordinate_tolerance:
                    n_exact_pairs += 1
                    for channel in CHANNELS:
                        left_value = float(
                            records[left_record_index][channel][left_local_index]
                        )
                        right_value = float(
                            records[right_record_index][channel][right_local_index]
                        )
                        if not np.isclose(
                            left_value,
                            right_value,
                            rtol=0.0,
                            atol=value_tolerance,
                        ):
                            raise ValueError(
                                "overlap value disagreement at exact tile coordinate "
                                f"channel {channel}"
                            )
                else:
                    n_phase_pairs += 1
                    for channel in CHANNELS:
                        difference = abs(
                            float(records[left_record_index][channel][left_local_index])
                            - float(
                                records[right_record_index][channel][right_local_index]
                            )
                        )
                        phase_differences[channel].append(difference)
                union(left_global, right_global)

    phase_statistics = {}
    for channel in CHANNELS:
        differences = np.asarray(phase_differences[channel], dtype=float)
        if differences.size:
            p999 = float(np.quantile(differences, 0.999))
            fraction_over = float(np.mean(differences > phase_p999_tolerance))
            maximum = float(np.max(differences))
        else:
            p999 = 0.0
            fraction_over = 0.0
            maximum = 0.0
        phase_statistics[channel] = {
            "p999_absolute_difference": p999,
            "fraction_over_0_1": fraction_over,
            "maximum_absolute_difference": maximum,
        }
        if (
            p999 > phase_p999_tolerance
            or fraction_over > phase_fraction_over_tolerance
            or maximum > phase_max_tolerance
        ):
            raise ValueError(
                "tile field disagreement after phase matching: "
                f"{channel} p99.9={p999:.6g}, "
                f"fraction_over_0.1={fraction_over:.6g}, max={maximum:.6g}"
            )

    all_coordinates = np.concatenate(coordinates_by_record, axis=0)
    all_channels = {
        channel: np.concatenate(
            [np.asarray(record[channel], dtype=float) for record in records]
        )
        for channel in CHANNELS
    }
    all_atom_indices = np.concatenate(
        [np.asarray(record["nearest_atom_index"], dtype=np.int32) for record in records]
    )
    all_residues = np.concatenate(
        [np.asarray(record["nearest_residue"]).astype(str) for record in records]
    )
    all_tiles = np.concatenate(
        [np.asarray(record["tile_provenance"]).astype(str) for record in records]
    )

    groups: dict[int, list[int]] = {}
    for index in range(total_points):
        groups.setdefault(find(index), []).append(index)

    representative_coordinates = []
    representative_channels = {channel: [] for channel in CHANNELS}
    representative_tiles = []
    fallback_atom_indices = []
    fallback_residues = []
    for members in groups.values():
        member_indices = np.asarray(members, dtype=np.int64)
        representative_coordinates.append(
            np.mean(all_coordinates[member_indices], axis=0)
        )
        for channel in CHANNELS:
            representative_channels[channel].append(
                float(np.mean(all_channels[channel][member_indices]))
            )
        tile_set = set()
        for tile_value in all_tiles[member_indices]:
            tile_set.update(str(tile_value).split("|"))
        representative_tiles.append("|".join(sorted(tile_set)))
        first_member = int(member_indices.min())
        fallback_atom_indices.append(int(all_atom_indices[first_member]))
        fallback_residues.append(str(all_residues[first_member]))

    coordinates = np.asarray(representative_coordinates, dtype=float)
    if len(coordinates):
        order = np.lexsort((coordinates[:, 2], coordinates[:, 1], coordinates[:, 0]))
        coordinates = coordinates[order]
    else:
        order = np.array([], dtype=np.int64)

    if (atom_coords is None) != (atom_labels is None):
        raise ValueError("atom_coords and atom_labels must be provided together")
    if atom_coords is not None:
        atom_coords = np.asarray(atom_coords, dtype=float)
        atom_labels = np.asarray(atom_labels)
        if atom_coords.ndim != 2 or atom_coords.shape[1] != 3:
            raise ValueError("atom coordinates must be N by 3")
        if len(atom_coords) != len(atom_labels):
            raise ValueError("atom coordinate and label counts differ")
        nearest_atom = cKDTree(atom_coords).query(coordinates, k=1)[1].astype(np.int32)
        nearest_residue = np.asarray(
            [_residue_from_atom_label(atom_labels[index]) for index in nearest_atom]
        )
    else:
        nearest_atom = np.asarray(fallback_atom_indices, dtype=np.int32)[order]
        nearest_residue = np.asarray(fallback_residues)[order]

    output = {
        "coordinates": coordinates.astype(np.float32),
        "nearest_atom_index": nearest_atom,
        "nearest_residue": nearest_residue,
        "tile_provenance": np.asarray(representative_tiles)[order],
    }
    for channel in CHANNELS:
        output[channel] = np.asarray(
            representative_channels[channel], dtype=np.float32
        )[order]
    output["neighbor_indptr"], output["neighbor_indices"] = coordinate_graph(
        coordinates, spacing, connectivity=26
    )
    output["merge_diagnostics"] = {
        "n_input_points": total_points,
        "n_output_points": int(len(coordinates)),
        "n_exact_cross_tile_pairs": int(n_exact_pairs),
        "n_phase_cross_tile_pairs": int(n_phase_pairs),
        "maximum_pair_distance": float(maximum_pair_distance),
        "match_tolerance": float(match_tolerance),
        "phase_field_statistics": phase_statistics,
    }
    return output

def classify_sasa_crosscheck(
    shell_area: float,
    freesasa_area: float,
    shell_residues: set[str],
    freesasa_residues: set[str],
    max_relative_area_difference: float = 0.50,
    min_residue_jaccard: float = 0.50,
) -> dict:
    if freesasa_area <= 0:
        raise ValueError("FreeSASA area must be positive")
    relative = abs(float(shell_area) - float(freesasa_area)) / float(freesasa_area)
    union = set(shell_residues) | set(freesasa_residues)
    intersection = set(shell_residues) & set(freesasa_residues)
    jaccard = len(intersection) / len(union) if union else 1.0
    status = (
        "SHELL_SASA_PASS"
        if relative <= max_relative_area_difference
        and jaccard >= min_residue_jaccard
        else "SHELL_SASA_DISAGREEMENT"
    )
    return {
        "status": status,
        "shell_area_proxy": float(shell_area),
        "freesasa_total_area": float(freesasa_area),
        "relative_area_difference": float(relative),
        "residue_jaccard": float(jaccard),
        "n_shell_exposed_residues": len(shell_residues),
        "n_freesasa_exposed_residues": len(freesasa_residues),
        "max_relative_area_difference": float(max_relative_area_difference),
        "min_residue_jaccard": float(min_residue_jaccard),
    }


def run_freesasa(pdb_path: Path, minimum_residue_area: float = 1.0) -> dict:
    try:
        import freesasa
    except ImportError as error:
        raise RuntimeError("Python FreeSASA module is unavailable") from error
    structure = freesasa.Structure(str(pdb_path))
    result = freesasa.calc(structure)
    residues = set()
    for chain, chain_areas in result.residueAreas().items():
        for residue_number, areas in chain_areas.items():
            if float(areas.total) >= minimum_residue_area:
                residues.add(f"{chain or '_'}:{str(residue_number).strip()}")
    return {"total_area": float(result.totalArea()), "residues": residues}


def _shell_residue_keys(labels: np.ndarray) -> set[str]:
    keys = set()
    for label in labels:
        fields = str(label).split(":")
        if len(fields) >= 3:
            keys.add(f"{fields[0]}:{fields[2]}")
    return keys


def load_tile(
    tile_dir: Path,
    atom_data: dict,
    shell_inner: float,
    shell_outer: float,
    probe_radius: float,
) -> dict:
    parsed = {
        channel: parse_autogrid_map(tile_dir / f"whole_receptor.{channel}.map")
        for channel in CHANNELS
    }
    first = parsed["A"]
    for channel in CHANNELS[1:]:
        other = parsed[channel]
        if (
            other["spacing"] != first["spacing"]
            or not np.array_equal(other["n_elements"], first["n_elements"])
            or not np.allclose(other["origin"], first["origin"], atol=1e-6, rtol=0)
        ):
            raise ValueError(f"nonidentical map geometry in {tile_dir}")
    arrays = {channel: parsed[channel]["values"] for channel in CHANNELS}
    occupied = build_occupancy(
        arrays["A"].shape,
        first["origin"],
        first["spacing"],
        atom_data["coordinates"],
        atom_data["radii"],
        probe_radius,
    )
    record = extract_shell_from_arrays(
        occupied,
        arrays,
        first["origin"],
        first["spacing"],
        shell_inner,
        shell_outer,
        atom_data["coordinates"],
        atom_data["radii"],
        atom_data["labels"],
        tile_dir.name,
    )
    return {
        "record": record,
        "occupied": occupied,
        "spacing": first["spacing"],
        "origin": first["origin"],
        "arrays": arrays,
    }


def write_json_exclusive(path: Path, payload: dict) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--maps-root", type=Path, required=True)
    parser.add_argument("--receptor-pdbqt", type=Path, required=True)
    parser.add_argument("--receptor-pdb", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--record-id", required=True)
    parser.add_argument("--shell-inner", type=float, default=0.0)
    parser.add_argument("--shell-outer", type=float, default=2.25)
    parser.add_argument(
        "--sensitivity-thicknesses", default="0.75,1.50,2.25,3.00"
    )
    parser.add_argument("--probe-radius", type=float, default=1.4)
    parser.add_argument("--max-relative-area-difference", type=float, default=0.50)
    parser.add_argument("--min-residue-jaccard", type=float, default=0.50)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    required = [
        args.maps_root / "MAPS_PASS.json",
        args.receptor_pdbqt,
        args.receptor_pdb,
    ]
    if any(not path.is_file() for path in required):
        print("one or more required shell inputs are missing", file=sys.stderr)
        return 2
    root = args.project_root.resolve()
    output = args.output_dir.resolve()
    if root != output and root not in output.parents:
        print("refusing output outside project root", file=sys.stderr)
        return 2
    if args.output_dir.exists():
        pass_gate = args.output_dir / "SHELL_PASS.json"
        if pass_gate.is_file():
            print(f"existing SHELL_PASS retained: {pass_gate}")
            return 0
        print("refusing partial or unverified shell output", file=sys.stderr)
        return 3

    args.output_dir.mkdir(parents=True)
    try:
        map_gate = json.loads(
            (args.maps_root / "MAPS_PASS.json").read_text(encoding="utf-8")
        )
        if map_gate.get("status") != "MAPS_PASS" or map_gate.get("channels") != list(CHANNELS):
            raise ValueError("incompatible MAPS_PASS gate")
        atom_data = parse_pdbqt_atoms(args.receptor_pdbqt)
        tile_dirs = sorted(
            path for path in args.maps_root.glob("tile_*") if path.is_dir()
        )
        if not tile_dirs:
            raise ValueError("no map tile directories")
        loaded = [
            load_tile(
                tile_dir,
                atom_data,
                args.shell_inner,
                args.shell_outer,
                args.probe_radius,
            )
            for tile_dir in tile_dirs
        ]
        spacings = {round(item["spacing"], 9) for item in loaded}
        if len(spacings) != 1:
            raise ValueError("tile spacings differ")
        spacing = loaded[0]["spacing"]
        merged = merge_shell_records(
            [item["record"] for item in loaded],
            spacing=spacing,
            atom_coords=atom_data["coordinates"],
            atom_labels=atom_data["labels"],
        )
        if len(merged["coordinates"]) == 0:
            raise ValueError("empty exterior shell")
        np.savez_compressed(
            args.output_dir / "surface_shell.npz",
            **{
                key: value
                for key, value in merged.items()
                if key != "merge_diagnostics"
            },
        )

        sensitivity = []
        for outer in [
            float(value) for value in args.sensitivity_thicknesses.split(",")
        ]:
            if outer <= args.shell_inner:
                raise ValueError("sensitivity thickness must exceed shell inner")
            records = []
            for tile_dir, item in zip(tile_dirs, loaded):
                records.append(
                    extract_shell_from_arrays(
                        item["occupied"],
                        item["arrays"],
                        item["origin"],
                        spacing,
                        args.shell_inner,
                        outer,
                        atom_data["coordinates"],
                        atom_data["radii"],
                        atom_data["labels"],
                        tile_dir.name,
                    )
                )
            merged_sensitivity = merge_shell_records(
                records,
                spacing=spacing,
                atom_coords=atom_data["coordinates"],
                atom_labels=atom_data["labels"],
            )
            sensitivity.append(
                {
                    "shell_outer": outer,
                    "n_shell_points": int(len(merged_sensitivity["coordinates"])),
                    "shell_area_proxy": shell_area_proxy(
                        len(merged_sensitivity["coordinates"]),
                        spacing,
                        args.shell_inner,
                        outer,
                    ),
                }
            )

        primary_area = shell_area_proxy(
            len(merged["coordinates"]),
            spacing,
            args.shell_inner,
            args.shell_outer,
        )
        sensitivity_counts = {
            round(item["shell_outer"], 9): item["n_shell_points"]
            for item in sensitivity
        }
        steiner_area = None
        if abs(args.shell_inner) <= 1e-12:
            steiner_keys = [round(spacing * multiple, 9) for multiple in (1, 2, 3)]
            if all(key in sensitivity_counts for key in steiner_keys):
                steiner_area = steiner_surface_area_from_counts(
                    [sensitivity_counts[key] for key in steiner_keys],
                    spacing,
                )
        freesasa_result = run_freesasa(args.receptor_pdb)
        crosscheck = classify_sasa_crosscheck(
            primary_area,
            freesasa_result["total_area"],
            _shell_residue_keys(merged["nearest_residue"]),
            freesasa_result["residues"],
            args.max_relative_area_difference,
            args.min_residue_jaccard,
        )
        payload = {
            "status": (
                "SHELL_PASS"
                if crosscheck["status"] == "SHELL_SASA_PASS"
                else "SHELL_SASA_DISAGREEMENT"
            ),
            "scientific_scope": "exterior_surface_shell_only",
            "record_id": args.record_id,
            "maps_gate_sha256": sha256(args.maps_root / "MAPS_PASS.json"),
            "receptor_pdbqt_sha256": sha256(args.receptor_pdbqt),
            "receptor_pdb_sha256": sha256(args.receptor_pdb),
            "channels": list(CHANNELS),
            "tile_count": len(tile_dirs),
            "spacing": spacing,
            "probe_radius": args.probe_radius,
            "shell_inner": args.shell_inner,
            "shell_outer": args.shell_outer,
            "connectivity": 26,
            "n_shell_points": int(len(merged["coordinates"])),
            "shell_area_proxy": primary_area,
            "steiner_surface_area_estimate": steiner_area,
            "tile_merge": merged["merge_diagnostics"],
            "sensitivity": sensitivity,
            "sasa_crosscheck": crosscheck,
        }
        gate_name = (
            "SHELL_PASS.json"
            if payload["status"] == "SHELL_PASS"
            else "SHELL_REVIEW.json"
        )
        write_json_exclusive(args.output_dir / gate_name, payload)
        checksums = []
        for path in sorted(args.output_dir.iterdir()):
            if path.is_file() and path.name != "SHA256SUMS":
                checksums.append(f"{sha256(path)}  {path.name}")
        (args.output_dir / "SHA256SUMS").write_text(
            "\n".join(checksums) + "\n", encoding="ascii"
        )
        print(json.dumps(payload, sort_keys=True))
        return 0 if payload["status"] == "SHELL_PASS" else 6
    except Exception as error:
        write_json_exclusive(
            args.output_dir / "SHELL_FAIL.json",
            {
                "status": "NOT_EVALUATED_SURFACE_SHELL",
                "record_id": args.record_id,
                "error": str(error),
            },
        )
        print(str(error), file=sys.stderr)
        return 5


if __name__ == "__main__":
    raise SystemExit(main())
