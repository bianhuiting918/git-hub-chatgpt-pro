#!/usr/bin/env python3
"""Core geometry utilities for explicit probe-to-surface-field matching.

The CLI and server I/O layer are added only after these scientific invariants pass
synthetic tests. Scores are screening proxies, never binding free energies.
"""

from __future__ import annotations

from typing import Iterable

import argparse
import csv
import hashlib
import json
import time
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree







def match_probe_conformer(
    *,
    complete_xyz: np.ndarray,
    typed_atom_indices: np.ndarray,
    typed_channels: Iterable[str],
    heavy_atom_indices: np.ndarray,
    heavy_atom_radii: np.ndarray,
    region_indices: np.ndarray,
    shell_xyz: np.ndarray,
    shell_raw_channels: dict[str, np.ndarray],
    maps: dict[str, dict],
    protein_xyz: np.ndarray,
    protein_radii: np.ndarray,
    relative_strain: float,
    anchor_minimum_separation: float,
    anchor_distance_tolerance: float,
    hard_clash_fraction: float,
    maximum_outside_shell_fraction: float,
    shell_cutoff: float,
    maximum_field_anchors: int,
    maximum_probe_triplets: int,
    maximum_matches: int,
) -> list[dict]:
    """Generate deterministic rigid placements for one explicit probe conformer."""
    complete_xyz = np.asarray(complete_xyz, dtype=float)
    typed_atom_indices = np.asarray(typed_atom_indices, dtype=np.int64)
    typed_channels = tuple(typed_channels)
    heavy_atom_indices = np.asarray(heavy_atom_indices, dtype=np.int64)
    heavy_atom_radii = np.asarray(heavy_atom_radii, dtype=float)
    region_indices = np.asarray(region_indices, dtype=np.int64)
    shell_xyz = np.asarray(shell_xyz, dtype=float)
    if len(typed_atom_indices) != len(typed_channels):
        raise ValueError("typed atom indices and channels differ")
    if len(heavy_atom_indices) != len(heavy_atom_radii):
        raise ValueError("heavy atom indices and radii differ")
    typed_xyz = complete_xyz[typed_atom_indices]
    field_coordinates: list[np.ndarray] = []
    field_types: list[str] = []
    for channel in sorted(set(typed_channels)):
        if channel not in shell_raw_channels:
            raise ValueError(f"missing shell channel {channel}")
        anchors = extract_field_anchors(
            shell_xyz,
            np.asarray(shell_raw_channels[channel], dtype=float),
            region_indices,
            minimum_separation=anchor_minimum_separation,
            maximum_anchors=maximum_field_anchors,
        )
        for anchor in anchors:
            field_coordinates.append(shell_xyz[int(anchor)])
            field_types.append(channel)
    if len(field_coordinates) < 3:
        return []
    field_xyz = np.vstack(field_coordinates)
    triplets = select_probe_anchor_triplets(
        typed_xyz,
        typed_channels,
        maximum_triplets=maximum_probe_triplets,
        minimum_triangle_area=0.1,
    )
    records: list[dict] = []
    seen = set()
    for triplet in triplets:
        triplet_array = np.asarray(triplet, dtype=np.int64)
        source = typed_xyz[triplet_array]
        source_types = tuple(typed_channels[index] for index in triplet)
        matches = enumerate_typed_anchor_matches(
            source,
            source_types,
            field_xyz,
            tuple(field_types),
            tolerance=anchor_distance_tolerance,
            maximum_matches=maximum_matches,
        )
        for assignment in matches:
            target = field_xyz[np.asarray(assignment, dtype=np.int64)]
            placed, anchor_rmsd = place_from_anchor_match(
                complete_xyz,
                typed_atom_indices[triplet_array],
                target,
            )
            heavy_xyz = placed[heavy_atom_indices]
            if has_hard_clash(
                heavy_xyz,
                heavy_atom_radii,
                protein_xyz,
                protein_radii,
                fraction=hard_clash_fraction,
            ):
                continue
            outside = outside_shell_fraction(
                heavy_xyz,
                shell_xyz[region_indices],
                cutoff=shell_cutoff,
            )
            if reject_outside_shell(outside, maximum_outside_shell_fraction):
                continue
            typed_placed = placed[typed_atom_indices]
            field_score, inside = score_pose_on_grids(
                typed_placed,
                typed_channels,
                maps,
            )
            if not bool(np.all(inside)):
                continue
            key = tuple(np.round(heavy_xyz.reshape(-1), 3).tolist())
            if key in seen:
                continue
            seen.add(key)
            records.append(
                {
                    "pose_id": f"pose_{len(records) + 1:06d}",
                    "score": total_pose_score(field_score, relative_strain),
                    "field_score": float(field_score),
                    "relative_strain": float(relative_strain),
                    "anchor_rmsd": float(anchor_rmsd),
                    "outside_shell_fraction": float(outside),
                    "heavy_xyz": heavy_xyz,
                    "complete_xyz": placed,
                    "probe_anchor_triplet": tuple(int(value) for value in triplet),
                    "field_anchor_assignment": tuple(int(value) for value in assignment),
                }
            )
    records.sort(
        key=lambda row: (
            -float(row["score"]),
            float(row["anchor_rmsd"]),
            str(row["pose_id"]),
        )
    )
    return records

def conformer_budget(rotatable_bonds: int) -> int:
    if rotatable_bonds < 0:
        raise ValueError("rotatable bond count must be nonnegative")
    if rotatable_bonds <= 2:
        return 8
    if rotatable_bonds <= 6:
        return 32
    return 64


def place_from_anchor_match(
    complete_xyz: np.ndarray,
    anchor_indices: np.ndarray,
    target_anchor_xyz: np.ndarray,
) -> tuple[np.ndarray, float]:
    complete_xyz = np.asarray(complete_xyz, dtype=float)
    anchor_indices = np.asarray(anchor_indices, dtype=np.int64)
    target_anchor_xyz = np.asarray(target_anchor_xyz, dtype=float)
    source_anchors = complete_xyz[anchor_indices]
    if source_anchors.shape != target_anchor_xyz.shape or len(source_anchors) < 3:
        raise ValueError("invalid anchor placement")
    source_center = source_anchors.mean(axis=0)
    target_center = target_anchor_xyz.mean(axis=0)
    source_zero = source_anchors - source_center
    target_zero = target_anchor_xyz - target_center
    left, _, right_t = np.linalg.svd(source_zero.T @ target_zero)
    rotation = right_t.T @ left.T
    if np.linalg.det(rotation) < 0:
        right_t[-1, :] *= -1
        rotation = right_t.T @ left.T
    placed = (complete_xyz - source_center) @ rotation.T + target_center
    placed_anchors = placed[anchor_indices]
    rmsd = float(
        np.sqrt(np.mean(np.sum((placed_anchors - target_anchor_xyz) ** 2, axis=1)))
    )
    return placed, rmsd


def score_pose_on_grids(
    atom_xyz: np.ndarray,
    atom_channels: Iterable[str],
    maps: dict[str, dict],
) -> tuple[float, np.ndarray]:
    atom_xyz = np.asarray(atom_xyz, dtype=float)
    atom_channels = tuple(atom_channels)
    if atom_xyz.shape != (len(atom_channels), 3):
        raise ValueError("atom channels and coordinates differ")
    inside = np.zeros(len(atom_xyz), dtype=bool)
    total = 0.0
    for channel in sorted(set(atom_channels)):
        if channel not in maps:
            raise ValueError(f"missing map channel {channel}")
        indices = np.asarray(
            [index for index, value in enumerate(atom_channels) if value == channel],
            dtype=np.int64,
        )
        grid = maps[channel]
        sampled, channel_inside = trilinear_interpolate(
            np.asarray(grid["values"], dtype=float),
            np.asarray(grid["origin"], dtype=float),
            float(grid["spacing"]),
            atom_xyz[indices],
        )
        inside[indices] = channel_inside
        total += float(np.sum(-sampled[channel_inside]))
    return total, inside

def trilinear_interpolate(
    values: np.ndarray,
    origin: np.ndarray,
    spacing: float,
    points: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Sample a regular xyz AutoGrid array and mark points inside its union."""
    values = np.asarray(values, dtype=float)
    origin = np.asarray(origin, dtype=float)
    points = np.asarray(points, dtype=float)
    if values.ndim != 3 or origin.shape != (3,) or points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("invalid grid geometry")
    if spacing <= 0:
        raise ValueError("grid spacing must be positive")
    fractional = (points - origin) / float(spacing)
    upper = np.asarray(values.shape, dtype=float) - 1.0
    inside = np.all((fractional >= 0.0) & (fractional <= upper), axis=1)
    sampled = np.full(len(points), np.nan, dtype=float)
    for row in np.flatnonzero(inside):
        coordinate = fractional[row]
        lower = np.floor(coordinate).astype(int)
        upper_index = np.minimum(lower + 1, np.asarray(values.shape) - 1)
        weight = coordinate - lower
        x0, y0, z0 = lower
        x1, y1, z1 = upper_index
        wx, wy, wz = weight
        sampled[row] = (
            values[x0, y0, z0] * (1 - wx) * (1 - wy) * (1 - wz)
            + values[x1, y0, z0] * wx * (1 - wy) * (1 - wz)
            + values[x0, y1, z0] * (1 - wx) * wy * (1 - wz)
            + values[x1, y1, z0] * wx * wy * (1 - wz)
            + values[x0, y0, z1] * (1 - wx) * (1 - wy) * wz
            + values[x1, y0, z1] * wx * (1 - wy) * wz
            + values[x0, y1, z1] * (1 - wx) * wy * wz
            + values[x1, y1, z1] * wx * wy * wz
        )
    return sampled, inside


def select_probe_anchor_triplets(
    coordinates: np.ndarray,
    channels: Iterable[str],
    maximum_triplets: int,
    minimum_triangle_area: float,
) -> list[tuple[int, int, int]]:
    """Return deterministic non-collinear typed atom triples."""
    from itertools import combinations

    coordinates = np.asarray(coordinates, dtype=float)
    channels = tuple(channels)
    if coordinates.shape != (len(channels), 3):
        raise ValueError("channel and coordinate counts differ")
    if maximum_triplets < 1 or minimum_triangle_area < 0:
        raise ValueError("invalid triplet parameters")
    output: list[tuple[int, int, int]] = []
    for triplet in combinations(range(len(channels)), 3):
        left, middle, right = coordinates[np.asarray(triplet, dtype=np.int64)]
        area = 0.5 * np.linalg.norm(np.cross(middle - left, right - left))
        if area < minimum_triangle_area - 1e-12:
            continue
        output.append(tuple(int(value) for value in triplet))
        if len(output) == maximum_triplets:
            break
    return output

def extract_field_anchors(
    coordinates: np.ndarray,
    raw_map_energy: np.ndarray,
    region_indices: np.ndarray,
    minimum_separation: float,
    maximum_anchors: int,
) -> np.ndarray:
    """Select deterministic local representatives from most to least favorable."""
    coordinates = np.asarray(coordinates, dtype=float)
    raw_map_energy = np.asarray(raw_map_energy, dtype=float)
    region_indices = np.unique(np.asarray(region_indices, dtype=np.int64))
    if coordinates.ndim != 2 or coordinates.shape[1] != 3:
        raise ValueError("coordinates must be N by 3")
    if len(raw_map_energy) != len(coordinates):
        raise ValueError("map values do not match coordinates")
    if minimum_separation < 0 or maximum_anchors < 1:
        raise ValueError("invalid anchor parameters")
    order = sorted(
        region_indices.tolist(),
        key=lambda index: (float(raw_map_energy[index]), int(index)),
    )
    selected: list[int] = []
    for index in order:
        if all(
            np.linalg.norm(coordinates[index] - coordinates[other])
            >= minimum_separation - 1e-12
            for other in selected
        ):
            selected.append(int(index))
            if len(selected) == maximum_anchors:
                break
    return np.asarray(selected, dtype=np.int64)


def enumerate_typed_anchor_matches(
    probe_xyz: np.ndarray,
    probe_types: Iterable[str],
    field_xyz: np.ndarray,
    field_types: Iterable[str],
    tolerance: float,
    maximum_matches: int,
) -> list[tuple[int, ...]]:
    """Enumerate deterministic typed field tuples satisfying all pair distances."""
    from itertools import product

    probe_xyz = np.asarray(probe_xyz, dtype=float)
    field_xyz = np.asarray(field_xyz, dtype=float)
    probe_types = tuple(probe_types)
    field_types = tuple(field_types)
    if len(probe_xyz) != len(probe_types) or len(field_xyz) != len(field_types):
        raise ValueError("type and coordinate counts differ")
    if maximum_matches < 1:
        raise ValueError("maximum matches must be positive")
    candidates = [
        [index for index, channel in enumerate(field_types) if channel == probe_type]
        for probe_type in probe_types
    ]
    if any(not group for group in candidates):
        return []
    matches: list[tuple[int, ...]] = []
    for assignment in product(*candidates):
        if len(set(assignment)) != len(assignment):
            continue
        canonical = True
        for left in range(len(probe_types)):
            for right in range(left + 1, len(probe_types)):
                if probe_types[left] == probe_types[right] and assignment[left] > assignment[right]:
                    canonical = False
                    break
            if not canonical:
                break
        if not canonical:
            continue
        target_xyz = field_xyz[np.asarray(assignment, dtype=np.int64)]
        if typed_geometry_compatible(
            probe_xyz,
            probe_types,
            target_xyz,
            probe_types,
            tolerance,
        ):
            matches.append(tuple(int(value) for value in assignment))
            if len(matches) == maximum_matches:
                break
    return matches


def probe_atom_channels(molecule) -> list[str]:
    """Map explicit probe atoms to AutoGrid A/C/OA/HD field channels."""
    channels: list[str] = []
    for atom in molecule.GetAtoms():
        atomic_number = atom.GetAtomicNum()
        if atomic_number == 6:
            channels.append("A" if atom.GetIsAromatic() else "C")
        elif atomic_number == 8 and atom.GetFormalCharge() <= 0:
            channels.append("OA")
        elif atomic_number == 1:
            neighbors = list(atom.GetNeighbors())
            if len(neighbors) == 1 and neighbors[0].GetAtomicNum() in {7, 8, 16}:
                channels.append("HD")
    return channels

def typed_geometry_compatible(
    probe_xyz: np.ndarray,
    probe_types: Iterable[str],
    field_xyz: np.ndarray,
    field_types: Iterable[str],
    tolerance: float,
) -> bool:
    probe_xyz = np.asarray(probe_xyz, dtype=float)
    field_xyz = np.asarray(field_xyz, dtype=float)
    probe_types = tuple(probe_types)
    field_types = tuple(field_types)
    if probe_xyz.shape != field_xyz.shape or probe_xyz.ndim != 2 or probe_xyz.shape[1] != 3:
        return False
    if len(probe_types) != len(field_types) or len(probe_types) != len(probe_xyz):
        return False
    if probe_types != field_types or tolerance < 0:
        return False
    probe_distances = np.linalg.norm(
        probe_xyz[:, None, :] - probe_xyz[None, :, :], axis=2
    )
    field_distances = np.linalg.norm(
        field_xyz[:, None, :] - field_xyz[None, :, :], axis=2
    )
    return bool(np.all(np.abs(probe_distances - field_distances) <= tolerance + 1e-12))


def rigid_fit(source_xyz: np.ndarray, target_xyz: np.ndarray) -> tuple[np.ndarray, float]:
    source_xyz = np.asarray(source_xyz, dtype=float)
    target_xyz = np.asarray(target_xyz, dtype=float)
    if source_xyz.shape != target_xyz.shape or source_xyz.ndim != 2 or source_xyz.shape[1] != 3:
        raise ValueError("source and target coordinates must be matching N by 3 arrays")
    if len(source_xyz) < 3:
        raise ValueError("at least three anchors are required")
    source_center = source_xyz.mean(axis=0)
    target_center = target_xyz.mean(axis=0)
    source_zero = source_xyz - source_center
    target_zero = target_xyz - target_center
    covariance = source_zero.T @ target_zero
    left, _, right_t = np.linalg.svd(covariance)
    rotation = right_t.T @ left.T
    if np.linalg.det(rotation) < 0:
        right_t[-1, :] *= -1
        rotation = right_t.T @ left.T
    placed = source_zero @ rotation.T + target_center
    rmsd = float(np.sqrt(np.mean(np.sum((placed - target_xyz) ** 2, axis=1))))
    return placed, rmsd


def has_hard_clash(
    ligand_xyz: np.ndarray,
    ligand_radii: np.ndarray,
    protein_xyz: np.ndarray,
    protein_radii: np.ndarray,
    fraction: float,
) -> bool:
    ligand_xyz = np.asarray(ligand_xyz, dtype=float)
    protein_xyz = np.asarray(protein_xyz, dtype=float)
    ligand_radii = np.asarray(ligand_radii, dtype=float)
    protein_radii = np.asarray(protein_radii, dtype=float)
    if not 0 < fraction <= 1:
        raise ValueError("clash fraction must be in (0, 1]")
    if ligand_xyz.shape != (len(ligand_radii), 3):
        raise ValueError("ligand radii do not match coordinates")
    if protein_xyz.shape != (len(protein_radii), 3):
        raise ValueError("protein radii do not match coordinates")
    if len(ligand_xyz) == 0 or len(protein_xyz) == 0:
        return False
    tree = cKDTree(protein_xyz)
    maximum_protein_radius = float(protein_radii.max())
    for coordinate, ligand_radius in zip(ligand_xyz, ligand_radii):
        cutoff = fraction * (float(ligand_radius) + maximum_protein_radius)
        for protein_index in tree.query_ball_point(coordinate, cutoff):
            distance = float(np.linalg.norm(coordinate - protein_xyz[protein_index]))
            threshold = fraction * (
                float(ligand_radius) + float(protein_radii[protein_index])
            )
            if distance < threshold:
                return True
    return False


def outside_shell_fraction(
    ligand_xyz: np.ndarray, shell_xyz: np.ndarray, cutoff: float
) -> float:
    ligand_xyz = np.asarray(ligand_xyz, dtype=float)
    shell_xyz = np.asarray(shell_xyz, dtype=float)
    if ligand_xyz.ndim != 2 or ligand_xyz.shape[1] != 3 or len(ligand_xyz) == 0:
        raise ValueError("ligand coordinates must be a nonempty N by 3 array")
    if shell_xyz.ndim != 2 or shell_xyz.shape[1] != 3 or len(shell_xyz) == 0:
        raise ValueError("shell coordinates must be a nonempty N by 3 array")
    if cutoff < 0:
        raise ValueError("shell cutoff must be nonnegative")
    distances, _ = cKDTree(shell_xyz).query(ligand_xyz, k=1)
    return float(np.mean(distances > cutoff))


def reject_outside_shell(fraction: float, maximum_fraction: float) -> bool:
    if not 0 <= fraction <= 1 or not 0 <= maximum_fraction <= 1:
        raise ValueError("fractions must be in [0, 1]")
    return bool(fraction > maximum_fraction + 1e-12)


def total_pose_score(field_score: float, relative_strain: float) -> float:
    if relative_strain < 0:
        raise ValueError("relative conformer strain must be nonnegative")
    return float(field_score - relative_strain)


def _pose_rmsd(left: np.ndarray, right: np.ndarray) -> float:
    """RMSD in the fixed receptor frame; translated sites are distinct poses."""
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if left.shape != right.shape:
        return float("inf")
    return float(np.sqrt(np.mean(np.sum((left - right) ** 2, axis=1))))


def cluster_pose_records(
    records: Iterable[dict], rmsd_cutoff: float, maximum_clusters: int
) -> list[dict]:
    if rmsd_cutoff < 0 or maximum_clusters < 1:
        raise ValueError("invalid clustering parameters")
    ordered = sorted(
        records,
        key=lambda row: (-float(row["score"]), str(row["pose_id"])),
    )
    representatives: list[dict] = []
    for record in ordered:
        coordinates = np.asarray(record["heavy_xyz"], dtype=float)
        if any(
            _pose_rmsd(coordinates, np.asarray(rep["heavy_xyz"], dtype=float))
            <= rmsd_cutoff
            for rep in representatives
        ):
            continue
        representatives.append(record)
        if len(representatives) == maximum_clusters:
            break
    return representatives


def select_active_probe_rows(rows: Iterable[dict], material_family: str) -> list[dict]:
    family = material_family.upper()
    return [
        row
        for row in rows
        if str(row.get("material_family", "")).upper() == family
        and str(row.get("status", "")).lower() == "active"
    ]


def composite_weight(role: str) -> float:
    role = str(role).lower()
    if role == "raw_monomer_control":
        return 0.0
    if role in {"primary_proxy", "diagnostic", "secondary_proxy"}:
        return 1.0
    return 0.0


AUTODOCK_VDW_RADII = {
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
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


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
        raise ValueError(f"AutoGrid value count mismatch: {values.size} != {expected}")
    values = values.reshape((shape[2], shape[1], shape[0])).transpose(2, 1, 0)
    origin = center - n_elements.astype(float) * spacing / 2.0
    return {
        "spacing": float(spacing),
        "origin": origin,
        "center": center,
        "n_elements": n_elements,
        "values": values,
    }


def read_pdbqt_atoms(path: Path) -> tuple[np.ndarray, np.ndarray]:
    coordinates = []
    radii = []
    with path.open(encoding="ascii", errors="strict") as handle:
        for line in handle:
            if line[:6] not in {"ATOM  ", "HETATM"}:
                continue
            coordinates.append(
                [float(line[30:38]), float(line[38:46]), float(line[46:54])]
            )
            atom_type = line.split()[-1].upper()
            radii.append(AUTODOCK_VDW_RADII.get(atom_type, 1.70))
    if not coordinates:
        raise ValueError(f"no receptor atoms in {path}")
    return np.asarray(coordinates, dtype=float), np.asarray(radii, dtype=float)


def typed_probe_atom_indices(molecule) -> tuple[np.ndarray, tuple[str, ...]]:
    indices = []
    channels = []
    for atom in molecule.GetAtoms():
        atomic_number = atom.GetAtomicNum()
        channel = None
        if atomic_number == 6:
            channel = "A" if atom.GetIsAromatic() else "C"
        elif atomic_number == 8 and atom.GetFormalCharge() <= 0:
            channel = "OA"
        elif atomic_number == 1:
            neighbors = list(atom.GetNeighbors())
            if len(neighbors) == 1 and neighbors[0].GetAtomicNum() in {7, 8, 16}:
                channel = "HD"
        if channel is not None:
            indices.append(atom.GetIdx())
            channels.append(channel)
    return np.asarray(indices, dtype=np.int64), tuple(channels)


def generate_probe_conformers(molecule, seed: int) -> list[dict]:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Lipinski

    work = Chem.Mol(molecule)
    work.RemoveAllConformers()
    rotors = int(Lipinski.NumRotatableBonds(work))
    budget = conformer_budget(rotors)
    parameters = AllChem.ETKDGv3()
    parameters.randomSeed = int(seed)
    parameters.pruneRmsThresh = 0.5
    parameters.numThreads = 1
    conformer_ids = list(AllChem.EmbedMultipleConfs(work, numConfs=budget, params=parameters))
    if not conformer_ids:
        return []
    energies = []
    if AllChem.MMFFHasAllMoleculeParams(work):
        properties = AllChem.MMFFGetMoleculeProperties(work)
        for conformer_id in conformer_ids:
            AllChem.MMFFOptimizeMolecule(
                work, mmffVariant="MMFF94s", confId=int(conformer_id), maxIters=500
            )
            force_field = AllChem.MMFFGetMoleculeForceField(
                work, properties, confId=int(conformer_id)
            )
            energies.append(float(force_field.CalcEnergy()))
        force_field_name = "MMFF94s"
    else:
        for conformer_id in conformer_ids:
            AllChem.UFFOptimizeMolecule(work, confId=int(conformer_id), maxIters=500)
            force_field = AllChem.UFFGetMoleculeForceField(work, confId=int(conformer_id))
            energies.append(float(force_field.CalcEnergy()))
        force_field_name = "UFF"
    minimum = min(energies)
    output = []
    for conformer_id, energy in zip(conformer_ids, energies):
        conformer = work.GetConformer(int(conformer_id))
        coordinates = np.asarray(conformer.GetPositions(), dtype=float)
        output.append(
            {
                "conformer_id": int(conformer_id),
                "coordinates": coordinates,
                "energy": float(energy),
                "relative_strain": float(max(0.0, energy - minimum)),
                "force_field": force_field_name,
            }
        )
    output.sort(key=lambda row: (row["relative_strain"], row["conformer_id"]))
    return output


def _write_tsv(path: Path, rows: list[dict], columns: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _load_probe_manifest(path: Path, material_family: str, probe_ids: set[str]) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    output = []
    for row in rows:
        if row["material_family"].upper() != material_family:
            continue
        if probe_ids and row["probe_id"] not in probe_ids:
            continue
        if row.get("validation_status") != "PASS":
            continue
        output.append(row)
    return output


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--shell-dir", type=Path, required=True)
    parser.add_argument("--patch-dir", type=Path, required=True)
    parser.add_argument("--maps-root", type=Path, required=True)
    parser.add_argument("--probe-sdf", type=Path, required=True)
    parser.add_argument("--probe-manifest", type=Path, required=True)
    parser.add_argument("--receptor-pdbqt", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--record-id", required=True)
    parser.add_argument("--material-family", choices=("PET", "NYLON"), required=True)
    parser.add_argument("--probe-ids", default="")
    parser.add_argument("--region", default="catalytic_neighborhood")
    parser.add_argument("--conformer-seed", type=int, default=12648430)
    parser.add_argument("--anchor-spacing", type=float, default=1.5)
    parser.add_argument("--anchor-tolerance", type=float, default=1.0)
    parser.add_argument("--hard-clash-fraction", type=float, default=0.75)
    parser.add_argument("--maximum-outside-shell-fraction", type=float, default=0.20)
    parser.add_argument("--pose-cluster-rmsd", type=float, default=2.0)
    parser.add_argument("--maximum-pose-clusters", type=int, default=20)
    parser.add_argument("--maximum-field-anchors", type=int, default=12)
    parser.add_argument("--maximum-probe-triplets", type=int, default=8)
    parser.add_argument("--maximum-matches", type=int, default=100)
    return parser.parse_args()


def main() -> int:
    from rdkit import Chem

    args = parse_args()
    started = time.perf_counter()
    if args.output_dir.exists():
        raise SystemExit(f"refusing existing output directory: {args.output_dir}")
    required = [
        args.shell_dir / "SHELL_PASS.json",
        args.shell_dir / "surface_shell.npz",
        args.patch_dir / "PATCH_PASS.json",
        args.patch_dir / "patch_membership.npz",
        args.probe_sdf,
        args.probe_manifest,
        args.receptor_pdbqt,
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise SystemExit("missing required inputs: " + ", ".join(missing))
    args.output_dir.mkdir(parents=True, exist_ok=False)

    shell = np.load(args.shell_dir / "surface_shell.npz", allow_pickle=False)
    shell_xyz = np.asarray(shell["coordinates"], dtype=float)
    shell_raw = {
        channel: np.asarray(shell[channel], dtype=float)
        for channel in ("A", "C", "OA", "HD")
    }
    memberships = np.load(args.patch_dir / "patch_membership.npz", allow_pickle=False)
    catalytic_key = f"catalytic__{args.region}"
    off_target_key = f"off_target__{args.region}__composite"
    if catalytic_key not in memberships or off_target_key not in memberships:
        raise SystemExit("requested equal-area patch membership is unavailable")
    regions = {
        "catalytic": np.asarray(memberships[catalytic_key], dtype=np.int64),
        "off_target": np.asarray(memberships[off_target_key], dtype=np.int64),
    }

    maps = {}
    tile = args.maps_root / "tile_0001"
    for channel in ("A", "C", "OA", "HD"):
        maps[channel] = parse_autogrid_map(tile / f"whole_receptor.{channel}.map")
    reference_shape = maps["A"]["values"].shape
    reference_origin = maps["A"]["origin"]
    reference_spacing = maps["A"]["spacing"]
    for channel, grid in maps.items():
        if (
            grid["values"].shape != reference_shape
            or not np.allclose(grid["origin"], reference_origin)
            or not np.isclose(grid["spacing"], reference_spacing)
        ):
            raise SystemExit(f"map geometry mismatch for {channel}")

    receptor_xyz, receptor_radii = read_pdbqt_atoms(args.receptor_pdbqt)
    requested_probe_ids = {value for value in args.probe_ids.split(",") if value}
    manifest_rows = _load_probe_manifest(
        args.probe_manifest, args.material_family, requested_probe_ids
    )
    if not manifest_rows:
        raise SystemExit("no validated probes selected")
    molecules = {
        molecule.GetProp("_Name"): molecule
        for molecule in Chem.SDMolSupplier(str(args.probe_sdf), removeHs=False)
        if molecule is not None and molecule.HasProp("_Name")
    }

    summaries = []
    failures = []
    sdf_path = args.output_dir / "probe_poses.sdf"
    writer = Chem.SDWriter(str(sdf_path))
    accepted_pose_count = 0
    for manifest_row in manifest_rows:
        probe_id = manifest_row["probe_id"]
        molecule = molecules.get(probe_id)
        if molecule is None:
            failures.append(
                {
                    "record_id": args.record_id,
                    "material_family": args.material_family,
                    "probe_id": probe_id,
                    "region": "all",
                    "reason": "NOT_EVALUATED_PROBE_SDF_MISSING",
                }
            )
            continue
        conformers = generate_probe_conformers(molecule, args.conformer_seed)
        if not conformers:
            failures.append(
                {
                    "record_id": args.record_id,
                    "material_family": args.material_family,
                    "probe_id": probe_id,
                    "region": "all",
                    "reason": "NOT_EVALUATED_PROBE_CONFORMER",
                }
            )
            continue
        typed_indices, typed_channels = typed_probe_atom_indices(molecule)
        heavy_indices = np.asarray(
            [atom.GetIdx() for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1],
            dtype=np.int64,
        )
        periodic = Chem.GetPeriodicTable()
        heavy_radii = np.asarray(
            [float(periodic.GetRvdw(molecule.GetAtomWithIdx(int(index)).GetAtomicNum())) for index in heavy_indices],
            dtype=float,
        )
        for region_name, region_indices in regions.items():
            pose_records = []
            for conformer in conformers:
                poses = match_probe_conformer(
                    complete_xyz=conformer["coordinates"],
                    typed_atom_indices=typed_indices,
                    typed_channels=typed_channels,
                    heavy_atom_indices=heavy_indices,
                    heavy_atom_radii=heavy_radii,
                    region_indices=region_indices,
                    shell_xyz=shell_xyz,
                    shell_raw_channels=shell_raw,
                    maps=maps,
                    protein_xyz=receptor_xyz,
                    protein_radii=receptor_radii,
                    relative_strain=conformer["relative_strain"],
                    anchor_minimum_separation=args.anchor_spacing,
                    anchor_distance_tolerance=args.anchor_tolerance,
                    hard_clash_fraction=args.hard_clash_fraction,
                    maximum_outside_shell_fraction=args.maximum_outside_shell_fraction,
                    shell_cutoff=1.125,
                    maximum_field_anchors=args.maximum_field_anchors,
                    maximum_probe_triplets=args.maximum_probe_triplets,
                    maximum_matches=args.maximum_matches,
                )
                for pose in poses:
                    pose["conformer_id"] = conformer["conformer_id"]
                    pose["force_field"] = conformer["force_field"]
                pose_records.extend(poses)
            clustered = cluster_pose_records(
                pose_records,
                rmsd_cutoff=args.pose_cluster_rmsd,
                maximum_clusters=args.maximum_pose_clusters,
            )
            if not clustered:
                failures.append(
                    {
                        "record_id": args.record_id,
                        "material_family": args.material_family,
                        "probe_id": probe_id,
                        "region": region_name,
                        "reason": "NOT_EVALUATED_PROBE_NO_ACCEPTED_POSE",
                    }
                )
                continue
            best = float(clustered[0]["score"])
            top5 = [float(row["score"]) for row in clustered[:5]]
            summaries.append(
                {
                    "record_id": args.record_id,
                    "material_family": args.material_family,
                    "probe_id": probe_id,
                    "probe_role": manifest_row["role"],
                    "primary_or_control": manifest_row["primary_or_control"],
                    "region": region_name,
                    "n_conformers": len(conformers),
                    "n_raw_poses": len(pose_records),
                    "n_pose_clusters": len(clustered),
                    "best_total_score": best,
                    "top5_mean_score": float(np.mean(top5)),
                    "count_within_2_of_best": sum(
                        float(row["score"]) >= best - 2.0 for row in clustered
                    ),
                    "best_field_score": float(clustered[0]["field_score"]),
                    "best_relative_strain": float(clustered[0]["relative_strain"]),
                    "best_anchor_rmsd": float(clustered[0]["anchor_rmsd"]),
                    "scientific_scope": "explicit_probe_surface_field_screening_proxy_only",
                }
            )
            for cluster_index, pose in enumerate(clustered, start=1):
                output_molecule = Chem.Mol(molecule)
                output_molecule.RemoveAllConformers()
                conformer = Chem.Conformer(output_molecule.GetNumAtoms())
                for atom_index, coordinate in enumerate(pose["complete_xyz"]):
                    conformer.SetAtomPosition(
                        atom_index,
                        (float(coordinate[0]), float(coordinate[1]), float(coordinate[2])),
                    )
                output_molecule.AddConformer(conformer, assignId=True)
                output_molecule.SetProp("_Name", f"{probe_id}_{region_name}_{cluster_index:03d}")
                output_molecule.SetProp("record_id", args.record_id)
                output_molecule.SetProp("material_family", args.material_family)
                output_molecule.SetProp("region", region_name)
                output_molecule.SetProp("score", f"{float(pose['score']):.8f}")
                output_molecule.SetProp("field_score", f"{float(pose['field_score']):.8f}")
                output_molecule.SetProp("relative_strain", f"{float(pose['relative_strain']):.8f}")
                writer.write(output_molecule)
                accepted_pose_count += 1
    writer.close()

    summary_columns = [
        "record_id", "material_family", "probe_id", "probe_role",
        "primary_or_control", "region", "n_conformers", "n_raw_poses",
        "n_pose_clusters", "best_total_score", "top5_mean_score",
        "count_within_2_of_best", "best_field_score",
        "best_relative_strain", "best_anchor_rmsd", "scientific_scope",
    ]
    failure_columns = [
        "record_id", "material_family", "probe_id", "region", "reason"
    ]
    _write_tsv(args.output_dir / "probe_pose_summary.tsv", summaries, summary_columns)
    _write_tsv(args.output_dir / "probe_failures.tsv", failures, failure_columns)
    elapsed = time.perf_counter() - started
    _write_tsv(
        args.output_dir / "runtime.tsv",
        [{
            "record_id": args.record_id,
            "material_family": args.material_family,
            "elapsed_seconds": f"{elapsed:.6f}",
            "n_selected_probes": len(manifest_rows),
            "n_summary_rows": len(summaries),
            "n_failures": len(failures),
            "n_written_pose_clusters": accepted_pose_count,
        }],
        [
            "record_id", "material_family", "elapsed_seconds", "n_selected_probes",
            "n_summary_rows", "n_failures", "n_written_pose_clusters",
        ],
    )
    files = [
        args.output_dir / "probe_pose_summary.tsv",
        args.output_dir / "probe_poses.sdf",
        args.output_dir / "probe_failures.tsv",
        args.output_dir / "runtime.tsv",
    ]
    with (args.output_dir / "SHA256SUMS").open("x", encoding="ascii") as handle:
        for output in files:
            handle.write(f"{sha256(output)}  {output.name}\n")
    gate = {
        "status": "PROBE_MATCH_TECHNICAL_PASS" if summaries else "NOT_EVALUATED_PROBE_MATCH",
        "record_id": args.record_id,
        "material_family": args.material_family,
        "region": args.region,
        "selected_probe_ids": [row["probe_id"] for row in manifest_rows],
        "n_summary_rows": len(summaries),
        "n_failure_rows": len(failures),
        "elapsed_seconds": elapsed,
        "scientific_scope": "technical_smoke_only_not_binding_energy_or_activity",
    }
    gate_name = "PROBE_MATCH_PASS.json" if summaries else "PROBE_MATCH_NOT_EVALUATED.json"
    with (args.output_dir / gate_name).open("x", encoding="utf-8") as handle:
        json.dump(gate, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0 if summaries else 4


if __name__ == "__main__":
    raise SystemExit(main())
