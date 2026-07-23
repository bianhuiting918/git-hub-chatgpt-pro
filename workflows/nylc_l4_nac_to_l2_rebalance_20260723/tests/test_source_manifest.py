#!/usr/bin/env python3
import os
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "scripts"))

from extract_sources import load_manifest, resolve_candidate  # noqa: E402


class SourceManifestTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.manifest_path = HERE / "manifests" / "candidates.json"
        cls.manifest = load_manifest(cls.manifest_path)

    def test_formal_universe_is_exactly_five_unique_candidates(self):
        candidates = self.manifest["candidates"]
        self.assertEqual(self.manifest["formal_candidate_count"], 5)
        self.assertEqual(len(candidates), 5)
        self.assertEqual(len({row["id"] for row in candidates}), 5)
        self.assertEqual(
            {row["id"] for row in candidates},
            {
                "nylc_c18_11854ps",
                "nylc_c23_29684ps",
                "nyl50_c18_70792ps",
                "nyl12_j1_37848ps",
                "nyl12_j2_87418ps",
            },
        )

    def test_each_candidate_resolves_to_unrestrained_selected_frame(self):
        for candidate in self.manifest["candidates"]:
            resolved = resolve_candidate(candidate)
            self.assertEqual(resolved["selected_row"]["segment"], candidate["segment"])
            self.assertEqual(resolved["selected_row"]["action"], "FREE_MONITOR")
            self.assertEqual(resolved["selected_row"]["joint_pass"], 1)
            self.assertEqual(resolved["selected_row"]["substrate_restrained"], 0)
            self.assertEqual(resolved["selected_row"]["gate_restrained"], 0)
            self.assertAlmostEqual(
                resolved["derived_local_xtc_time_ps"],
                candidate["local_xtc_time_ps"],
                places=6,
            )
            self.assertTrue(Path(resolved["xtc"]).is_file())
            self.assertTrue(Path(resolved["tpr"]).is_file())
            self.assertTrue(Path(resolved["source_itp"]).is_file())
            self.assertTrue(Path(resolved["source_topology"]).is_file())
            self.assertTrue(Path(resolved["source_index"]).is_file())

    def test_global_and_local_reactive_indices_share_one_ligand_offset(self):
        for candidate in self.manifest["candidates"]:
            first = candidate["ligand_first_global_atom"]
            local = candidate["reactive_local_atoms"]
            global_atoms = candidate["reactive_global_atoms"]
            for key in ("carbonyl_c", "carbonyl_o", "amide_n"):
                self.assertEqual(global_atoms[key], first + local[key] - 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
