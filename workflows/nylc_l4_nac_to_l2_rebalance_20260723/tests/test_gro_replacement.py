#!/usr/bin/env python3
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "scripts"))

from build_l4_to_l2 import (  # noqa: E402
    build_l2_coordinates,
    discover_heavy_mapping,
    extract_ligand_coordinates,
    minimum_ligand_environment_distance,
    parse_gro,
    parse_itp,
    replace_ligand,
    write_gro,
)


class GroReplacementTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.task_root = Path(os.environ["TASK_ROOT"])
        cls.manifest = json.loads(
            (HERE / "manifests" / "candidates.json").read_text(encoding="utf-8")
        )
        cls.candidate = next(
            row for row in cls.manifest["candidates"] if row["id"] == "nylc_c23_29684ps"
        )
        cls.source_itp = parse_itp(
            Path(cls.candidate["branch_dir"]) / cls.candidate["source_itp"]
        )
        cls.l2_itp = parse_itp(Path(os.environ["TASK_L2_ITP"]))
        cls.source_gro = parse_gro(
            cls.task_root / "inputs" / cls.candidate["id"] / "source.gro"
        )
        cls.heavy_mapping = discover_heavy_mapping(
            cls.source_itp,
            cls.l2_itp,
            (
                cls.candidate["reactive_local_atoms"]["carbonyl_c"],
                cls.candidate["reactive_local_atoms"]["carbonyl_o"],
                cls.candidate["reactive_local_atoms"]["amide_n"],
            ),
        )
        source_coordinates = extract_ligand_coordinates(
            cls.source_gro,
            cls.candidate["ligand_first_global_atom"],
            cls.candidate["source_atom_count"],
        )
        cls.l2_coordinates, cls.atom_mapping = build_l2_coordinates(
            cls.source_itp, cls.l2_itp, cls.heavy_mapping, source_coordinates
        )
        cls.replaced = replace_ligand(
            cls.source_gro,
            cls.candidate["ligand_first_global_atom"],
            cls.candidate["source_atom_count"],
            cls.l2_itp,
            cls.l2_coordinates,
        )

    def test_replacement_changes_only_ligand_block_and_atom_count(self):
        first = self.candidate["ligand_first_global_atom"]
        source_count = self.candidate["source_atom_count"]
        self.assertEqual(
            len(self.replaced.atoms),
            len(self.source_gro.atoms) - source_count + 79,
        )
        for before, after in zip(
            self.source_gro.atoms[: first - 1], self.replaced.atoms[: first - 1]
        ):
            self.assertEqual(before.coordinate, after.coordinate)
            self.assertEqual((before.resid, before.resname, before.atom_name),
                             (after.resid, after.resname, after.atom_name))
        source_tail = self.source_gro.atoms[first - 1 + source_count :]
        replaced_tail = self.replaced.atoms[first - 1 + 79 :]
        self.assertEqual(
            [atom.coordinate for atom in source_tail],
            [atom.coordinate for atom in replaced_tail],
        )

    def test_reactive_coordinates_are_exactly_preserved(self):
        first = self.candidate["ligand_first_global_atom"]
        for target_id in self.heavy_mapping.target_reactive_triplet:
            source_id = self.heavy_mapping.target_to_source[target_id]
            source_atom = self.source_gro.atoms[first + source_id - 2]
            target_atom = self.replaced.atoms[first + target_id - 2]
            self.assertEqual(source_atom.coordinate, target_atom.coordinate)

    def test_written_gro_round_trips(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "rebuilt.gro"
            write_gro(self.replaced, path)
            reread = parse_gro(path)
            self.assertEqual(len(reread.atoms), len(self.replaced.atoms))
            self.assertEqual(reread.box, self.replaced.box)
            for observed, expected in zip(reread.atoms, self.replaced.atoms):
                for observed_value, expected_value in zip(
                    observed.coordinate, expected.coordinate
                ):
                    self.assertAlmostEqual(observed_value, expected_value, places=3)

    def test_minimum_ligand_environment_distance_is_finite(self):
        first = self.candidate["ligand_first_global_atom"]
        distance = minimum_ligand_environment_distance(self.replaced, first, 79)
        self.assertGreater(distance, 0.0)
        self.assertLess(distance, 2.0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
