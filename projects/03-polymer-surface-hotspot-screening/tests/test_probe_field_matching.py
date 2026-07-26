from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "match_probe_fields.py"
SPEC = importlib.util.spec_from_file_location("match_probe_fields", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load {SCRIPT}")
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def test_dmt_like_pattern_requires_joint_aromatic_and_acceptor_geometry():
    probe_xyz = np.array(
        [[0.0, 0.0, 0.0], [2.8, 0.0, 0.0], [-2.8, 0.0, 0.0]],
        dtype=float,
    )
    probe_types = ("A", "OA", "OA")
    compatible_xyz = probe_xyz + np.array([5.0, -3.0, 2.0])
    incompatible_xyz = np.array(
        [[5.0, -3.0, 2.0], [11.0, -3.0, 2.0], [-1.0, -3.0, 2.0]],
        dtype=float,
    )
    assert mod.typed_geometry_compatible(
        probe_xyz, probe_types, compatible_xyz, probe_types, tolerance=1.0
    )
    assert not mod.typed_geometry_compatible(
        probe_xyz, probe_types, incompatible_xyz, probe_types, tolerance=1.0
    )


def test_rigid_fit_is_translation_and_rotation_invariant():
    source = np.array(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
        dtype=float,
    )
    rotation = np.array(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]],
        dtype=float,
    )
    target = source @ rotation.T + np.array([4.0, -2.0, 3.0])
    placed, rmsd = mod.rigid_fit(source, target)
    assert rmsd < 1e-8
    assert np.allclose(placed, target, atol=1e-8)


def test_hard_clash_uses_frozen_vdw_fraction():
    protein_xyz = np.array([[0.0, 0.0, 0.0]])
    protein_radii = np.array([1.70])
    ligand_radii = np.array([1.70])
    assert mod.has_hard_clash(
        np.array([[2.0, 0.0, 0.0]]),
        ligand_radii,
        protein_xyz,
        protein_radii,
        fraction=0.75,
    )
    assert not mod.has_hard_clash(
        np.array([[3.0, 0.0, 0.0]]),
        ligand_radii,
        protein_xyz,
        protein_radii,
        fraction=0.75,
    )


def test_outside_shell_fraction_rejects_more_than_twenty_percent():
    shell_xyz = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    ligand_xyz = np.array(
        [[0.1, 0.0, 0.0], [0.9, 0.0, 0.0], [5.0, 0.0, 0.0]]
    )
    fraction = mod.outside_shell_fraction(ligand_xyz, shell_xyz, cutoff=0.5)
    assert fraction == pytest.approx(1.0 / 3.0)
    assert mod.reject_outside_shell(fraction, maximum_fraction=0.20)


def test_conformer_strain_penalizes_equal_field_fit():
    assert mod.total_pose_score(field_score=12.0, relative_strain=1.5) == pytest.approx(10.5)
    assert mod.total_pose_score(field_score=12.0, relative_strain=0.2) > mod.total_pose_score(
        field_score=12.0, relative_strain=1.5
    )


def test_pose_clustering_is_deterministic_and_keeps_best():
    base = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    records = [
        {"pose_id": "b", "score": 8.0, "heavy_xyz": base + 0.10},
        {"pose_id": "a", "score": 9.0, "heavy_xyz": base},
        {"pose_id": "c", "score": 7.0, "heavy_xyz": base + 5.0},
    ]
    first = mod.cluster_pose_records(records, rmsd_cutoff=0.5, maximum_clusters=20)
    second = mod.cluster_pose_records(list(reversed(records)), rmsd_cutoff=0.5, maximum_clusters=20)
    assert [x["pose_id"] for x in first] == ["a", "c"]
    assert [x["pose_id"] for x in second] == ["a", "c"]


def test_material_probe_universes_remain_separate():
    rows = [
        {"probe_id": "PET_DMT", "material_family": "PET", "status": "active"},
        {"probe_id": "NYL_NMA", "material_family": "NYLON", "status": "active"},
        {"probe_id": "PA66_CAPPED_DIMER", "material_family": "NYLON", "status": "deferred"},
    ]
    assert [r["probe_id"] for r in mod.select_active_probe_rows(rows, "PET")] == ["PET_DMT"]
    assert [r["probe_id"] for r in mod.select_active_probe_rows(rows, "NYLON")] == ["NYL_NMA"]


def test_raw_monomer_controls_are_not_primary_composite_evidence():
    assert mod.composite_weight("primary_proxy") > 0.0
    assert mod.composite_weight("diagnostic") > 0.0
    assert mod.composite_weight("raw_monomer_control") == 0.0
