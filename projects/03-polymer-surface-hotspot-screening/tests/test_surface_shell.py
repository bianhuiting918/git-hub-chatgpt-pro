import importlib.util
import math
import tempfile
import unittest
from pathlib import Path

import numpy as np


PROJECT = Path(__file__).resolve().parents[1]
SCRIPT = PROJECT / "scripts" / "extract_surface_shell.py"

shell = None
if SCRIPT.is_file():
    spec = importlib.util.spec_from_file_location("extract_surface_shell", SCRIPT)
    shell = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(shell)


def sphere_occupancy(shape, center, radius, origin=(0.0, 0.0, 0.0), spacing=1.0):
    axes = [
        origin[i] + np.arange(shape[i], dtype=float) * spacing
        for i in range(3)
    ]
    x, y, z = np.meshgrid(*axes, indexing="ij")
    return (
        (x - center[0]) ** 2
        + (y - center[1]) ** 2
        + (z - center[2]) ** 2
    ) <= radius**2


@unittest.skipIf(shell is None, "implementation intentionally absent during RED")
class SurfaceShellUnitTests(unittest.TestCase):
    def test_boundary_flood_fill_excludes_enclosed_cavity(self):
        occupied = np.zeros((21, 21, 21), dtype=bool)
        occupied[4:17, 4:17, 4:17] = True
        occupied[8:13, 8:13, 8:13] = False

        exterior = shell.exterior_connected_solvent(occupied, connectivity=26)

        self.assertTrue(exterior[0, 0, 0])
        self.assertFalse(exterior[10, 10, 10])
        self.assertFalse(np.any(exterior & occupied))

    def test_shell_band_is_bounded_within_one_grid_cell(self):
        occupied = np.zeros((21, 21, 21), dtype=bool)
        occupied[10, 10, 10] = True

        mask, distances, exterior = shell.exterior_shell_mask(
            occupied, spacing=1.0, inner=0.0, outer=2.25, connectivity=26
        )

        selected = distances[mask]
        self.assertGreater(selected.size, 0)
        self.assertGreater(selected.min(), 0.0)
        self.assertLessEqual(selected.max(), 2.25 + 1e-12)
        self.assertFalse(mask[10, 10, 10])
        self.assertTrue(np.all(exterior[mask]))

    def test_rigid_translation_preserves_values_and_summary(self):
        occupied = sphere_occupancy((25, 25, 25), (12.0, 12.0, 12.0), 4.0)
        channels = {
            name: np.full(occupied.shape, value, dtype=np.float32)
            for name, value in {"A": 1.0, "C": 2.0, "OA": 3.0, "HD": 4.0}.items()
        }
        atom_coords = np.array([[12.0, 12.0, 12.0]])
        common = dict(
            occupied=occupied,
            channel_arrays=channels,
            spacing=1.0,
            shell_inner=0.0,
            shell_outer=2.25,
            atom_radii=np.array([4.0]),
            atom_labels=np.array(["A:SPH:1:C"]),
            tile_id="tile",
        )

        first = shell.extract_shell_from_arrays(
            origin=np.array([0.0, 0.0, 0.0]),
            atom_coords=atom_coords,
            **common,
        )
        shift = np.array([7.5, -3.0, 11.25])
        second = shell.extract_shell_from_arrays(
            origin=shift,
            atom_coords=atom_coords + shift,
            **common,
        )

        self.assertEqual(first["coordinates"].shape, second["coordinates"].shape)
        np.testing.assert_allclose(
            second["coordinates"] - first["coordinates"],
            np.broadcast_to(shift, first["coordinates"].shape),
        )
        for channel in ("A", "C", "OA", "HD"):
            np.testing.assert_array_equal(first[channel], second[channel])
        self.assertEqual(first["summary"]["n_shell_points"], second["summary"]["n_shell_points"])
        self.assertAlmostEqual(
            first["summary"]["shell_area_proxy"],
            second["summary"]["shell_area_proxy"],
        )

    def test_two_lobed_object_keeps_connected_exterior_shell(self):
        left = sphere_occupancy((33, 25, 25), (10.0, 12.0, 12.0), 5.0)
        right = sphere_occupancy((33, 25, 25), (22.0, 12.0, 12.0), 5.0)
        bridge = np.zeros_like(left)
        bridge[10:23, 11:14, 11:14] = True
        occupied = left | right | bridge

        mask, _, _ = shell.exterior_shell_mask(
            occupied, spacing=1.0, inner=0.0, outer=2.25, connectivity=26
        )

        self.assertTrue(np.any(mask[:16]))
        self.assertTrue(np.any(mask[17:]))
        self.assertTrue(np.any(mask[14:19, 9:16, 9:16]))

    def test_graph_index_is_symmetric_and_has_no_self_edges(self):
        occupied = np.zeros((9, 9, 9), dtype=bool)
        occupied[4, 4, 4] = True
        channels = {name: np.zeros(occupied.shape) for name in ("A", "C", "OA", "HD")}
        result = shell.extract_shell_from_arrays(
            occupied=occupied,
            channel_arrays=channels,
            origin=np.zeros(3),
            spacing=1.0,
            shell_inner=0.0,
            shell_outer=1.5,
            atom_coords=np.array([[4.0, 4.0, 4.0]]),
            atom_radii=np.array([1.0]),
            atom_labels=np.array(["A:ONE:1:C"]),
            tile_id="tile_0001",
        )

        indptr = result["neighbor_indptr"]
        indices = result["neighbor_indices"]
        adjacency = {
            i: set(indices[indptr[i] : indptr[i + 1]].tolist())
            for i in range(len(indptr) - 1)
        }
        for i, neighbors in adjacency.items():
            self.assertNotIn(i, neighbors)
            for j in neighbors:
                self.assertIn(i, adjacency[j])

    def test_tiled_overlap_merges_identically_to_untiled_coordinates(self):
        spacing = 1.0
        full_shape = (25, 17, 17)
        full_origin = np.array([0.0, 0.0, 0.0])
        occupied = sphere_occupancy(full_shape, (12.0, 8.0, 8.0), 4.0)
        base = np.indices(full_shape).sum(axis=0).astype(np.float32)
        channels = {name: base + i for i, name in enumerate(("A", "C", "OA", "HD"))}
        atom_coords = np.array([[12.0, 8.0, 8.0]])
        atom_radii = np.array([4.0])
        atom_labels = np.array(["A:SPH:1:C"])

        whole = shell.extract_shell_from_arrays(
            occupied, channels, full_origin, spacing, 0.0, 2.25,
            atom_coords, atom_radii, atom_labels, "whole"
        )
        records = []
        for tile_id, slc in (
            ("tile_left", slice(0, 17)),
            ("tile_right", slice(8, 25)),
        ):
            records.append(shell.extract_shell_from_arrays(
                occupied[slc, :, :],
                {name: array[slc, :, :] for name, array in channels.items()},
                np.array([float(slc.start), 0.0, 0.0]),
                spacing,
                0.0,
                2.25,
                atom_coords,
                atom_radii,
                atom_labels,
                tile_id,
            ))
        merged = shell.merge_shell_records(records, spacing=spacing)

        whole_lookup = {
            tuple(np.rint(coord / spacing).astype(int)): i
            for i, coord in enumerate(whole["coordinates"])
        }
        merged_lookup = {
            tuple(np.rint(coord / spacing).astype(int)): i
            for i, coord in enumerate(merged["coordinates"])
        }
        self.assertEqual(set(whole_lookup), set(merged_lookup))
        for key in whole_lookup:
            wi, mi = whole_lookup[key], merged_lookup[key]
            for channel in ("A", "C", "OA", "HD"):
                self.assertAlmostEqual(float(whole[channel][wi]), float(merged[channel][mi]))

    def test_phase_shifted_tile_overlap_is_deduplicated(self):
        def record(coordinates, values, residues, tile):
            coordinates = np.asarray(coordinates, dtype=float)
            values = np.asarray(values, dtype=float)
            return {
                "coordinates": coordinates,
                "A": values.copy(),
                "C": values.copy(),
                "OA": values.copy(),
                "HD": values.copy(),
                "nearest_atom_index": np.arange(len(coordinates), dtype=np.int32),
                "nearest_residue": np.asarray(residues),
                "tile_provenance": np.full(len(coordinates), tile),
            }

        left = record(
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            [0.0, 1.0, 2.0],
            ["A:RES:1", "A:RES:2", "A:RES:3"],
            "tile_left",
        )
        shifted = record(
            [[1.25, 0.0, 0.0], [2.25, 0.0, 0.0], [3.25, 0.0, 0.0]],
            [1.0, 2.0, 3.0],
            ["A:RES:2", "A:RES:3", "A:RES:4"],
            "tile_shifted",
        )

        merged = shell.merge_shell_records([left, shifted], spacing=1.0)

        self.assertEqual(len(merged["coordinates"]), 4)
        self.assertEqual(
            sum("|" in str(value) for value in merged["tile_provenance"]), 2
        )

    def test_grid_origin_phase_rejects_nearby_nonduplicate_pair(self):
        def one_point(coordinate, origin, tile):
            return {
                "coordinates": np.asarray([coordinate], dtype=float),
                "A": np.asarray([1.0]),
                "C": np.asarray([1.0]),
                "OA": np.asarray([1.0]),
                "HD": np.asarray([1.0]),
                "nearest_atom_index": np.asarray([0], dtype=np.int32),
                "nearest_residue": np.asarray(["A:RES:1"]),
                "tile_provenance": np.asarray([tile]),
                "grid_origin": np.asarray(origin, dtype=float),
            }

        left = one_point([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], "left")
        shifted_but_not_duplicate = one_point(
            [0.35, 0.0, 0.0], [0.25, 0.0, 0.0], "right"
        )

        merged = shell.merge_shell_records(
            [left, shifted_but_not_duplicate], spacing=1.0
        )

        self.assertEqual(len(merged["coordinates"]), 2)
        self.assertEqual(
            sum("|" in str(value) for value in merged["tile_provenance"]), 0
        )

    def test_geometry_only_sensitivity_merge_does_not_apply_field_gate(self):
        def one_point(coordinate, origin, value, tile):
            return {
                "coordinates": np.asarray([coordinate], dtype=float),
                "A": np.asarray([value]),
                "C": np.asarray([value]),
                "OA": np.asarray([value]),
                "HD": np.asarray([value]),
                "nearest_atom_index": np.asarray([0], dtype=np.int32),
                "nearest_residue": np.asarray(["A:RES:1"]),
                "tile_provenance": np.asarray([tile]),
                "grid_origin": np.asarray(origin, dtype=float),
            }

        left = one_point([0.0, 0.0, 0.0], [0.0, 0.0, 0.0], 0.0, "left")
        shifted = one_point([0.25, 0.0, 0.0], [0.25, 0.0, 0.0], 1.0, "right")

        with self.assertRaisesRegex(ValueError, "tile field disagreement"):
            shell.merge_shell_records([left, shifted], spacing=1.0)
        geometry_only = shell.merge_shell_records(
            [left, shifted],
            spacing=1.0,
            validate_phase_fields=False,
            build_graph=False,
        )

        self.assertEqual(len(geometry_only["coordinates"]), 1)
        self.assertEqual(len(geometry_only["neighbor_indices"]), 0)

    def test_coordinate_graph_connects_phase_shifted_tile_seam(self):
        coordinates = np.array(
            [[0.0, 0.0, 0.0], [1.10, 0.0, 0.0], [2.10, 0.0, 0.0]]
        )
        indptr, indices = shell.coordinate_graph(
            coordinates, spacing=1.0, connectivity=26
        )
        visited = {0}
        frontier = [0]
        while frontier:
            node = frontier.pop()
            for neighbor in indices[indptr[node] : indptr[node + 1]]:
                neighbor = int(neighbor)
                if neighbor not in visited:
                    visited.add(neighbor)
                    frontier.append(neighbor)

        self.assertEqual(visited, {0, 1, 2})

    def test_merge_rejects_disagreeing_overlap_values(self):
        record = {
            "coordinates": np.array([[0.0, 0.0, 0.0]]),
            "A": np.array([1.0]),
            "C": np.array([1.0]),
            "OA": np.array([1.0]),
            "HD": np.array([1.0]),
            "nearest_atom_index": np.array([0]),
            "nearest_residue": np.array(["A:ALA:1"]),
            "tile_provenance": np.array(["a"]),
        }
        conflict = {key: value.copy() if hasattr(value, "copy") else value for key, value in record.items()}
        conflict["A"][0] = 2.0
        with self.assertRaisesRegex(ValueError, "overlap"):
            shell.merge_shell_records([record, conflict], spacing=1.0)

    def test_autogrid_parser_uses_n_elements_plus_one_points(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "tiny.map"
            path.write_text(
                "GRID_PARAMETER_FILE tiny.gpf\n"
                "GRID_DATA_FILE tiny.fld\n"
                "MACROMOLECULE rec.pdbqt\n"
                "SPACING 1.0\n"
                "NELEMENTS 1 1 1\n"
                "CENTER 0.0 0.0 0.0\n"
                "0\n1\n2\n3\n4\n5\n6\n7\n",
                encoding="ascii",
            )
            parsed = shell.parse_autogrid_map(path)
            self.assertEqual(parsed["values"].shape, (2, 2, 2))
            np.testing.assert_allclose(parsed["origin"], [-0.5, -0.5, -0.5])
            self.assertEqual(float(parsed["values"][0, 0, 0]), 0.0)
            self.assertEqual(float(parsed["values"][1, 1, 1]), 7.0)

    def test_shell_area_proxy_is_reasonable_for_synthetic_sphere(self):
        spacing = 0.5
        shape = (41, 41, 41)
        center = (10.0, 10.0, 10.0)
        radius = 5.0
        occupied = sphere_occupancy(shape, center, radius, spacing=spacing)
        mask, _, _ = shell.exterior_shell_mask(
            occupied, spacing=spacing, inner=0.0, outer=1.0, connectivity=26
        )
        estimate = shell.shell_area_proxy(mask.sum(), spacing, 0.0, 1.0)
        expected = 4.0 * math.pi * radius * radius
        self.assertLess(abs(estimate - expected) / expected, 0.35)

    def test_steiner_estimator_removes_offset_curvature_bias(self):
        spacing = 0.75
        radius = 3.1
        margin = radius + 4.0 * spacing
        n = int(math.ceil(2.0 * margin / spacing)) + 1
        center = (margin, margin, margin)
        occupied = sphere_occupancy(
            (n, n, n), center, radius, spacing=spacing
        )
        cumulative_counts = []
        for multiple in (1.0, 2.0, 3.0):
            mask, _, _ = shell.exterior_shell_mask(
                occupied,
                spacing=spacing,
                inner=0.0,
                outer=multiple * spacing,
                connectivity=26,
            )
            cumulative_counts.append(int(mask.sum()))

        estimate = shell.steiner_surface_area_from_counts(
            cumulative_counts, spacing
        )
        old_offset_volume_proxy = shell.shell_area_proxy(
            cumulative_counts[-1], spacing, 0.0, 3.0 * spacing
        )
        expected = 4.0 * math.pi * radius * radius

        self.assertGreater(
            abs(old_offset_volume_proxy - expected) / expected, 0.50
        )
        self.assertLess(abs(estimate - expected) / expected, 0.20)

    def test_sasa_disagreement_is_a_review_status_not_biological_failure(self):
        good = shell.classify_sasa_crosscheck(
            shell_area=100.0,
            freesasa_area=110.0,
            shell_residues={"A:1", "A:2"},
            freesasa_residues={"A:1", "A:2"},
            max_relative_area_difference=0.35,
            min_residue_jaccard=0.6,
        )
        bad = shell.classify_sasa_crosscheck(
            shell_area=30.0,
            freesasa_area=110.0,
            shell_residues={"A:1"},
            freesasa_residues={"A:2"},
            max_relative_area_difference=0.35,
            min_residue_jaccard=0.6,
        )
        self.assertEqual(good["status"], "SHELL_SASA_PASS")
        self.assertEqual(bad["status"], "SHELL_SASA_DISAGREEMENT")
        self.assertNotIn("inactive", str(bad).lower())


class SurfaceShellRedGate(unittest.TestCase):
    def test_surface_shell_implementation_exists(self):
        self.assertTrue(
            SCRIPT.is_file(),
            "extract_surface_shell.py is intentionally missing: RED gate",
        )


if __name__ == "__main__":
    unittest.main()
