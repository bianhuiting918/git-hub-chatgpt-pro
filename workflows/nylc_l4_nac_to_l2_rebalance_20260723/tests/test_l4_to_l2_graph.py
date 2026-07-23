#!/usr/bin/env python3
import os
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "scripts"))

from build_l4_to_l2 import (  # noqa: E402
    build_atom_mapping,
    build_l2_coordinates,
    discover_heavy_mapping,
    parse_itp,
    validate_reactive_triplet,
)


class ItpGraphTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.l2 = parse_itp(Path(os.environ["TASK_L2_ITP"]))
        cls.nylc = parse_itp(Path(os.environ["NYLC_L4_ITP"]))
        cls.nyl12 = parse_itp(Path(os.environ["NYL12_L4_ITP"]))

    def test_audited_l2_counts_and_charge(self):
        self.assertEqual(len(self.l2.atoms), 79)
        heavy = [a for a in self.l2.atoms.values() if a.element != "H"]
        self.assertEqual(len(heavy), 33)
        self.assertAlmostEqual(sum(a.charge for a in self.l2.atoms.values()), 0.0, places=6)

    def test_reactive_triplets_are_bonded_graph_facts(self):
        for topology, triplet in [
            (self.nylc, (25, 26, 24)),
            (self.nylc, (31, 32, 33)),
            (self.nyl12, (18, 19, 20)),
            (self.nyl12, (34, 35, 36)),
        ]:
            validate_reactive_triplet(topology, *triplet)

    def assert_formal_mapping(self, topology, triplet):
        result = discover_heavy_mapping(topology, self.l2, triplet)
        self.assertEqual(len(result.target_to_source), 32)
        self.assertNotIn(33, result.target_to_source)
        self.assertEqual(len(result.cut_edges), 2)
        allowed = {1, 31}
        self.assertEqual({edge.target_endpoint for edge in result.cut_edges}, allowed)
        source_c, source_o, source_n = triplet
        reverse = {source: target for target, source in result.target_to_source.items()}
        target_c = reverse[source_c]
        target_o = reverse[source_o]
        target_n = reverse[source_n]
        self.assertIn((target_c, target_o, target_n), {(9, 10, 8), (15, 16, 17), (25, 26, 24)})
        return result

    def test_nylc_c23_matches_audited_regression_oracle(self):
        result = self.assert_formal_mapping(self.nylc, (31, 32, 33))
        expected = {
            1: 56,
            2: 55, 3: 54, 4: 53, 5: 52, 6: 51, 7: 50,
            8: 49, 9: 47, 10: 48,
            11: 46, 12: 45, 13: 44, 14: 43,
            15: 41, 16: 42, 17: 40,
            18: 39, 19: 38, 20: 37, 21: 36, 22: 35, 23: 34,
            24: 33, 25: 31, 26: 32,
            27: 30, 28: 29, 29: 28, 30: 27,
            31: 25, 32: 26,
        }
        self.assertEqual(result.target_to_source, expected)

    def test_nylc_c18_is_discovered_independently(self):
        c18 = self.assert_formal_mapping(self.nylc, (25, 26, 24))
        c23 = discover_heavy_mapping(self.nylc, self.l2, (31, 32, 33))
        self.assertNotEqual(c18.target_to_source, c23.target_to_source)
        self.assertEqual({c18.target_to_source[25], c18.target_to_source[15], c18.target_to_source[9]} & {25}, {25})

    def test_nyl12_j1_and_j2_use_nyl12_graph(self):
        j1 = self.assert_formal_mapping(self.nyl12, (18, 19, 20))
        j2 = self.assert_formal_mapping(self.nyl12, (34, 35, 36))
        self.assertNotEqual(j1.target_to_source, j2.target_to_source)
        self.assertIn(18, j1.target_to_source.values())
        self.assertIn(34, j2.target_to_source.values())

    def test_nonterminal_mapped_atoms_have_no_unmapped_heavy_neighbors(self):
        result = self.assert_formal_mapping(self.nylc, (31, 32, 33))
        mapped_sources = set(result.target_to_source.values())
        for target_id, source_id in result.target_to_source.items():
            extras = {
                neighbor
                for neighbor in self.nylc.neighbors[source_id]
                if self.nylc.atoms[neighbor].element != "H" and neighbor not in mapped_sources
            }
            if target_id not in {1, 31}:
                self.assertEqual(extras, set())

    def test_hydrogens_follow_bonded_parent_equivalence(self):
        result = discover_heavy_mapping(self.nylc, self.l2, (31, 32, 33))
        atom_mapping = build_atom_mapping(self.nylc, self.l2, result)
        self.assertEqual(set(atom_mapping), set(self.l2.atoms))
        rebuilt = {target for target, source in atom_mapping.items() if source is None}
        self.assertEqual(rebuilt, {33, 34, 35, 36})
        for target_id, source_id in atom_mapping.items():
            if target_id in rebuilt:
                continue
            if self.l2.atoms[target_id].element == "H":
                target_parent = next(iter(self.l2.neighbors[target_id]))
                source_parent = next(iter(self.nylc.neighbors[source_id]))
                self.assertEqual(atom_mapping[target_parent], source_parent)

    def test_coordinate_builder_copies_inherited_atoms_and_rebuilds_endpoints(self):
        result = discover_heavy_mapping(self.nylc, self.l2, (31, 32, 33))
        source_xyz = {
            atom_id: (0.01 * atom_id, 0.003 * atom_id, -0.002 * atom_id)
            for atom_id in self.nylc.atoms
        }
        # Give endpoint neighborhoods non-collinear geometry.
        source_xyz[56] = (0.0, 0.0, 0.0)
        source_xyz[55] = (0.145, 0.0, 0.0)
        source_xyz[57] = (-0.08, 0.10, 0.0)
        source_xyz[148] = (-0.03, -0.02, 0.095)
        source_xyz[25] = (1.0, 1.0, 1.0)
        source_xyz[27] = (1.15, 1.0, 1.0)
        source_xyz[24] = (0.94, 1.11, 1.0)
        built, atom_mapping = build_l2_coordinates(
            self.nylc, self.l2, result, source_xyz
        )
        for target_id, source_id in atom_mapping.items():
            if source_id is not None:
                self.assertEqual(built[target_id], source_xyz[source_id])
        self.assertEqual(built[25], source_xyz[31])
        self.assertEqual(built[26], source_xyz[32])
        self.assertEqual(built[24], source_xyz[33])

        def distance(a, b):
            return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5

        self.assertAlmostEqual(distance(built[31], built[33]), 0.12190, places=5)
        for hydrogen in (34, 35, 36):
            self.assertAlmostEqual(distance(built[1], built[hydrogen]), 0.10271, places=5)


if __name__ == "__main__":
    unittest.main(verbosity=2)
