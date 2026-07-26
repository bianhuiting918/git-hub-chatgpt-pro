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



def test_field_anchor_extraction_prefers_favorable_values_and_enforces_spacing():
    coordinates = np.array(
        [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [3.0, 0.0, 0.0], [7.0, 0.0, 0.0]]
    )
    raw_map_energy = np.array([-5.0, -4.0, -3.0, 2.0])
    anchors = mod.extract_field_anchors(
        coordinates,
        raw_map_energy,
        np.arange(4),
        minimum_separation=1.5,
        maximum_anchors=3,
    )
    assert anchors.tolist() == [0, 2, 3]


def test_typed_anchor_enumeration_uses_joint_pairwise_geometry():
    probe_xyz = np.array(
        [[0.0, 0.0, 0.0], [2.8, 0.0, 0.0], [-2.8, 0.0, 0.0]]
    )
    probe_types = ("A", "OA", "OA")
    field_xyz = np.array(
        [
            [10.0, 0.0, 0.0],
            [12.8, 0.0, 0.0],
            [7.2, 0.0, 0.0],
            [16.0, 0.0, 0.0],
            [4.0, 0.0, 0.0],
        ]
    )
    field_types = ("A", "OA", "OA", "OA", "OA")
    matches = mod.enumerate_typed_anchor_matches(
        probe_xyz,
        probe_types,
        field_xyz,
        field_types,
        tolerance=1.0,
        maximum_matches=20,
    )
    assert len(matches) == 1
    assert matches[0] == (0, 1, 2)


def test_probe_atom_channels_recover_pet_aromatic_and_ester_pattern():
    from rdkit import Chem

    molecule = Chem.AddHs(Chem.MolFromSmiles("COC(=O)c1ccc(C(=O)OC)cc1"))
    channels = mod.probe_atom_channels(molecule)
    assert channels.count("A") == 6
    assert channels.count("OA") == 4
    assert channels.count("C") == 4
    assert channels.count("HD") == 0


def test_probe_atom_channels_recover_nylon_amide_donor_and_acceptor():
    from rdkit import Chem

    molecule = Chem.AddHs(Chem.MolFromSmiles("CNC(C)=O"))
    channels = mod.probe_atom_channels(molecule)
    assert channels.count("OA") == 1
    assert channels.count("HD") == 1
    assert channels.count("C") == 3



def test_trilinear_interpolation_recovers_linear_field():
    x, y, z = np.meshgrid(
        np.arange(3.0), np.arange(3.0), np.arange(3.0), indexing="ij"
    )
    values = x + 2.0 * y + 3.0 * z
    points = np.array([[0.2, 0.3, 0.4], [1.25, 1.5, 0.75]])
    sampled, inside = mod.trilinear_interpolate(
        values, origin=np.zeros(3), spacing=1.0, points=points
    )
    assert inside.tolist() == [True, True]
    assert sampled.tolist() == pytest.approx([2.0, 6.5])


def test_trilinear_interpolation_marks_points_outside_grid():
    values = np.zeros((2, 2, 2), dtype=float)
    sampled, inside = mod.trilinear_interpolate(
        values,
        origin=np.zeros(3),
        spacing=1.0,
        points=np.array([[0.5, 0.5, 0.5], [2.0, 0.0, 0.0]]),
    )
    assert inside.tolist() == [True, False]
    assert sampled[0] == pytest.approx(0.0)
    assert np.isnan(sampled[1])


def test_noncollinear_probe_triplets_are_deterministic():
    coordinates = np.array(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0], [4.0, 0.0, 0.0]]
    )
    channels = ("A", "OA", "OA", "C")
    triplets = mod.select_probe_anchor_triplets(
        coordinates, channels, maximum_triplets=8, minimum_triangle_area=0.1
    )
    assert triplets[0] == (0, 1, 2)
    assert all(len(set(item)) == 3 for item in triplets)
    assert (0, 1, 3) not in triplets



@pytest.mark.parametrize("rotatable_bonds,expected", [(0, 8), (2, 8), (3, 32), (6, 32), (7, 64)])
def test_conformer_budget_follows_frozen_rotatable_bond_tiers(rotatable_bonds, expected):
    assert mod.conformer_budget(rotatable_bonds) == expected


def test_rigid_transform_places_complete_probe_from_three_anchors():
    complete = np.array(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 2.0]]
    )
    anchor_indices = np.array([0, 1, 2])
    rotation = np.array(
        [[0.0, -1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0]]
    )
    target_anchors = complete[anchor_indices] @ rotation.T + np.array([5.0, 2.0, -1.0])
    placed, anchor_rmsd = mod.place_from_anchor_match(
        complete, anchor_indices, target_anchors
    )
    expected = complete @ rotation.T + np.array([5.0, 2.0, -1.0])
    assert anchor_rmsd < 1e-8
    assert np.allclose(placed, expected, atol=1e-8)


def test_pose_field_score_sums_typed_favorable_autogrid_values():
    grid_a = np.full((3, 3, 3), -2.0)
    grid_oa = np.full((3, 3, 3), -3.0)
    maps = {
        "A": {"values": grid_a, "origin": np.zeros(3), "spacing": 1.0},
        "OA": {"values": grid_oa, "origin": np.zeros(3), "spacing": 1.0},
    }
    score, inside = mod.score_pose_on_grids(
        np.array([[0.5, 0.5, 0.5], [1.5, 1.5, 1.5]]),
        ("A", "OA"),
        maps,
    )
    assert inside.tolist() == [True, True]
    assert score == pytest.approx(5.0)



def test_match_probe_conformer_accepts_joint_pattern_and_scores_complete_pose():
    probe_xyz = np.array(
        [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0]]
    )
    channels = ("A", "OA", "OA")
    shell_xyz = probe_xyz + np.array([4.0, 4.0, 4.0])
    shell_raw = {
        "A": np.array([-2.0, 0.0, 0.0]),
        "OA": np.array([0.0, -3.0, -3.0]),
    }
    values_a = np.full((10, 10, 10), -2.0)
    values_oa = np.full((10, 10, 10), -3.0)
    maps = {
        "A": {"values": values_a, "origin": np.zeros(3), "spacing": 1.0},
        "OA": {"values": values_oa, "origin": np.zeros(3), "spacing": 1.0},
    }
    poses = mod.match_probe_conformer(
        complete_xyz=probe_xyz,
        typed_atom_indices=np.array([0, 1, 2]),
        typed_channels=channels,
        heavy_atom_indices=np.array([0, 1, 2]),
        heavy_atom_radii=np.array([1.7, 1.52, 1.52]),
        region_indices=np.array([0, 1, 2]),
        shell_xyz=shell_xyz,
        shell_raw_channels=shell_raw,
        maps=maps,
        protein_xyz=np.empty((0, 3)),
        protein_radii=np.empty(0),
        relative_strain=0.5,
        anchor_minimum_separation=1.5,
        anchor_distance_tolerance=1.0,
        hard_clash_fraction=0.75,
        maximum_outside_shell_fraction=0.20,
        shell_cutoff=0.6,
        maximum_field_anchors=8,
        maximum_probe_triplets=8,
        maximum_matches=20,
    )
    assert poses
    assert poses[0]["field_score"] == pytest.approx(8.0)
    assert poses[0]["score"] == pytest.approx(7.5)
    assert poses[0]["outside_shell_fraction"] == pytest.approx(0.0)
