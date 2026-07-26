#!/usr/bin/env python3
"""Core geometry utilities for explicit probe-to-surface-field matching.

The CLI and server I/O layer are added only after these scientific invariants pass
synthetic tests. Scores are screening proxies, never binding free energies.
"""

from __future__ import annotations

from typing import Iterable

import numpy as np
from scipy.spatial import cKDTree


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
    distances = np.linalg.norm(
        ligand_xyz[:, None, :] - protein_xyz[None, :, :], axis=2
    )
    thresholds = fraction * (
        ligand_radii[:, None] + protein_radii[None, :]
    )
    return bool(np.any(distances < thresholds))


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
    left = np.asarray(left, dtype=float)
    right = np.asarray(right, dtype=float)
    if left.shape != right.shape:
        return float("inf")
    centered_left = left - left.mean(axis=0)
    centered_right = right - right.mean(axis=0)
    covariance = centered_left.T @ centered_right
    u, _, vt = np.linalg.svd(covariance)
    rotation = vt.T @ u.T
    if np.linalg.det(rotation) < 0:
        vt[-1, :] *= -1
        rotation = vt.T @ u.T
    aligned = centered_left @ rotation.T
    return float(np.sqrt(np.mean(np.sum((aligned - centered_right) ** 2, axis=1))))


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
