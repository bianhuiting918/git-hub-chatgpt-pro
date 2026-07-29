import importlib.util
import unittest
from pathlib import Path

import numpy as np


PROJECT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT / "scripts" / "score_surface_patches.py"

patches = None
if SCRIPT.is_file():
    spec = importlib.util.spec_from_file_location("score_surface_patches", SCRIPT)
    patches = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(patches)


def csr_graph(n_nodes, edges):
    adjacency = [set() for _ in range(n_nodes)]
    for left, right in edges:
        adjacency[left].add(right)
        adjacency[right].add(left)
    indices = []
    indptr = [0]
    for neighbors in adjacency:
        indices.extend(sorted(neighbors))
        indptr.append(len(indices))
    return np.asarray(indptr, dtype=np.int64), np.asarray(indices, dtype=np.int32)


def line_shell(raw_a, labels=None, disconnected=False):
    raw_a = np.asarray(raw_a, dtype=float)
    n_nodes = len(raw_a)
    edges = [(i, i + 1) for i in range(n_nodes - 1)]
    if disconnected:
        edges = [edge for edge in edges if edge != (4, 5)]
    indptr, indices = csr_graph(n_nodes, edges)
    if labels is None:
        labels = ["A:ALA:1"] + [f"A:GLY:{i + 2}" for i in range(n_nodes - 1)]
    return {
        "coordinates": np.column_stack(
            [np.arange(n_nodes, dtype=float), np.zeros(n_nodes), np.zeros(n_nodes)]
        ),
        "neighbor_indptr": indptr,
        "neighbor_indices": indices,
        "nearest_residue": np.asarray(labels),
        "channels": {
            "A": raw_a,
            "C": raw_a.copy(),
            "OA": raw_a.copy(),
            "HD": raw_a.copy(),
        },
    }


@unittest.skipIf(patches is None, "implementation intentionally absent during RED")
class PatchScoringUnitTests(unittest.TestCase):
    def test_top_fraction_distinguishes_equal_mean_patches(self):
        concentrated = patches.summarize_values(np.array([10.0, 0.0, 0.0, 0.0, 0.0]), 0.20)
        diffuse = patches.summarize_values(np.array([2.0, 2.0, 2.0, 2.0, 2.0]), 0.20)
        self.assertEqual(concentrated["mean"], diffuse["mean"])
        self.assertGreater(concentrated["top_fraction_mean"], diffuse["top_fraction_mean"])

    def test_farthest_point_seeds_are_deterministic(self):
        coordinates = np.array(
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.0, 2.0, 0.0],
             [2.0, 2.0, 0.0], [1.0, 1.0, 0.0]]
        )
        eligible = np.array([4, 3, 2, 1, 0])
        first = patches.deterministic_farthest_point_seeds(coordinates, eligible, 4)
        second = patches.deterministic_farthest_point_seeds(coordinates, eligible, 4)
        np.testing.assert_array_equal(first, second)
        self.assertEqual(len(set(first.tolist())), 4)
        self.assertEqual(first[0], 0)

    def test_catalytic_hotspot_scores_above_equal_area_off_target(self):
        shell = line_shell([-10, -10, -10, -10, -1, -1, -1, -1, -1, -1, -1, -1])
        result = patches.score_shell_data(
            **shell,
            catalytic_residues={"A:1"},
            material_family="PET",
            region_radii={"catalytic_core": 3.1},
            off_target_buffer=1.0,
            top_fraction=0.20,
            max_seeds=32,
        )
        self.assertEqual(result["status"], "PATCH_SCORE_PASS")
        channel = result["regions"]["catalytic_core"]["channels"]["A"]
        self.assertGreater(channel["delta_top_fraction_mean"], 0.0)
        self.assertEqual(
            result["regions"]["catalytic_core"]["n_catalytic_points"],
            result["regions"]["catalytic_core"]["n_off_target_points"],
        )

    def test_off_target_hotspot_can_beat_catalytic_patch(self):
        shell = line_shell([-5, -5, -5, -5, -1, -1, -12, -12, -12, -12, -1, -1])
        result = patches.score_shell_data(
            **shell,
            catalytic_residues={"A:1"},
            material_family="PET",
            region_radii={"catalytic_core": 3.1},
            off_target_buffer=1.0,
            top_fraction=0.20,
            max_seeds=32,
        )
        channel = result["regions"]["catalytic_core"]["channels"]["A"]
        self.assertLess(channel["delta_top_fraction_mean"], 0.0)
        self.assertGreater(channel["off_target"]["top_fraction_mean"], channel["catalytic"]["top_fraction_mean"])

    def test_uniformly_sticky_surface_is_not_mistaken_for_local_selectivity(self):
        shell = line_shell([-8] * 16)
        result = patches.score_shell_data(
            **shell,
            catalytic_residues={"A:1"},
            material_family="PET",
            region_radii={"catalytic_core": 3.1},
            off_target_buffer=1.0,
            top_fraction=0.20,
            max_seeds=32,
        )
        region = result["regions"]["catalytic_core"]
        self.assertAlmostEqual(region["global_sticky_fraction"], 1.0)
        self.assertAlmostEqual(
            region["channels"]["A"]["delta_top_fraction_mean"], 0.0
        )

    def test_disconnected_off_target_smaller_than_target_is_not_evaluated(self):
        shell = line_shell([-5] * 8, disconnected=True)
        result = patches.score_shell_data(
            **shell,
            catalytic_residues={"A:1"},
            material_family="PET",
            region_radii={"catalytic_core": 4.1},
            off_target_buffer=0.0,
            top_fraction=0.20,
            max_seeds=16,
        )
        self.assertEqual(result["status"], "NOT_EVALUATED_EQUAL_AREA_OFF_TARGET")
        self.assertNotIn("inactive", str(result).lower())

    def test_missing_catalytic_mapping_is_not_a_biological_failure(self):
        shell = line_shell([-5] * 10)
        result = patches.score_shell_data(
            **shell,
            catalytic_residues={"A:999"},
            material_family="NYLON",
            region_radii={"catalytic_core": 3.0},
            off_target_buffer=1.0,
            top_fraction=0.20,
            max_seeds=16,
        )
        self.assertEqual(result["status"], "NOT_EVALUATED_CATALYTIC_MAPPING")
        self.assertEqual(result["missing_catalytic_residues"], ["A:999"])
        self.assertNotIn("inactive", str(result).lower())

    def test_pet_and_nylon_use_separate_frozen_channel_sets(self):
        shell = line_shell([-5] * 12)
        pet = patches.score_shell_data(
            **shell,
            catalytic_residues={"A:1"},
            material_family="PET",
            region_radii={"catalytic_core": 3.1},
            off_target_buffer=1.0,
            top_fraction=0.20,
            max_seeds=16,
        )
        nylon = patches.score_shell_data(
            **shell,
            catalytic_residues={"A:1"},
            material_family="NYLON",
            region_radii={"catalytic_core": 3.1},
            off_target_buffer=1.0,
            top_fraction=0.20,
            max_seeds=16,
        )
        self.assertEqual(pet["material_channels"], ["A", "C", "OA"])
        self.assertEqual(nylon["material_channels"], ["C", "OA", "HD"])

    def test_connected_patch_growth_has_exact_requested_size(self):
        indptr, indices = csr_graph(10, [(i, i + 1) for i in range(9)])
        eligible = np.ones(10, dtype=bool)
        patch = patches.grow_connected_patch(5, indptr, indices, eligible, 6)
        self.assertEqual(len(patch), 6)
        self.assertIn(5, patch)


class PatchScoringRedGate(unittest.TestCase):
    def test_patch_scoring_implementation_exists(self):
        self.assertTrue(
            SCRIPT.is_file(),
            "score_surface_patches.py is intentionally missing: RED gate",
        )


if __name__ == "__main__":
    unittest.main()
