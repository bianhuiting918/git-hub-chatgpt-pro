#!/usr/bin/env python3
"""Contract for the NylC A1 four-anchor bidirectional Step1 gap fill."""

from __future__ import annotations

import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "prepare_audit_nylc_a1_step1_four_anchor_bidirectional_gapfill.py"
RUNNER = ROOT / "scripts" / "run_nylc_a1_step1_four_anchor_bidirectional_gapfill.sh"
SBATCH = ROOT / "slurm" / "run_nylc_a1_step1_four_anchor_bidirectional_gapfill.sbatch"


def load_driver():
    spec = importlib.util.spec_from_file_location("_a1_four_anchor_gapfill", DRIVER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {DRIVER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class FourAnchorBidirectionalGapfillContract(unittest.TestCase):
    def test_four_anchors_two_directions_two_seeds_make_sixteen_unthrottled_tasks(self):
        self.assertTrue(DRIVER.is_file(), "RED: four-anchor bidirectional driver is missing")
        self.assertTrue(RUNNER.is_file(), "RED: four-anchor bidirectional runner is missing")
        self.assertTrue(SBATCH.is_file(), "RED: four-anchor bidirectional sbatch is missing")
        driver = load_driver()
        description = driver.describe()
        self.assertEqual(description["anchors"], ["R", "I1", "I2", "P"])
        self.assertEqual(description["directions"], ["REACTANT", "PRODUCT"])
        self.assertEqual(description["seed_count"], 2)
        self.assertEqual(description["array_task_count"], 16)
        mappings = [driver.task_spec(i) for i in range(16)]
        self.assertEqual(
            {(x["seed_index"], x["anchor"], x["direction"]) for x in mappings},
            {(s, a, d) for s in range(2) for a in ("R", "I1", "I2", "P")
             for d in ("REACTANT", "PRODUCT")},
        )
        self.assertEqual(
            driver.REACTANT_SOURCES[0]["restart_sha256"],
            "9d93b60aee9e8d6fa97493396d757f3bf4d0a47595591968e3debaebb2640e86",
        )
        self.assertEqual(
            driver.REACTANT_SOURCES[1]["restart_sha256"],
            "a12c018e195c32b34239004aea36ebc80de8275c9f9332b44084ad8d88d4db0e",
        )
        self.assertEqual(
            driver.PRODUCT_SOURCES[0]["restart_sha256"],
            "5d8f76d2c90e3e8c707b640c55a93938f53e18dc30d6f93d54adda467da25f41",
        )
        self.assertEqual(
            driver.PRODUCT_SOURCES[1]["restart_sha256"],
            "4cc60ad4d7be099b3f76040f1ee8c49b91deecebb20ed98570bc192b3d511132",
        )
        self.assertEqual(driver.ANCHOR_TARGETS["I1"]["proton_site"], "NALPHA")
        self.assertEqual(driver.ANCHOR_TARGETS["I2"]["proton_site"], "N3")
        self.assertTrue(driver.ANCHOR_TARGETS["I1"]["tetrahedral_required"])
        self.assertTrue(driver.ANCHOR_TARGETS["I2"]["tetrahedral_required"])
        self.assertTrue(description["stage_restart_persisted"])
        self.assertTrue(description["restrained_structures_are_not_ts"])
        runner = RUNNER.read_text(encoding="utf-8")
        self.assertIn("SLURM_ARRAY_TASK_ID", runner)
        self.assertIn("stage.rst7", runner)
        self.assertIn("SHA256.tsv", runner)
        sbatch = SBATCH.read_text(encoding="utf-8")
        self.assertIn("#SBATCH -n 8", sbatch)
        self.assertIn("#SBATCH --array=0-15", sbatch)
        self.assertNotIn("#SBATCH --array=0-15%", sbatch)


if __name__ == "__main__":
    unittest.main()
