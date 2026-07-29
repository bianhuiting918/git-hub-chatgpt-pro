#!/usr/bin/env python3
"""Contract test for the NylC A1 Step1 balanced four-distance bridge scout."""

from __future__ import annotations

import importlib.util
import pathlib
import unittest


WORKFLOW_ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPTS = WORKFLOW_ROOT / "scripts"
SLURM = WORKFLOW_ROOT / "slurm"
DRIVER = SCRIPTS / "prepare_audit_nylc_a1_step1_balanced_four_distance_bridge.py"
RUNNER = SCRIPTS / "run_nylc_a1_step1_balanced_four_distance_bridge.sh"
SBATCH = SLURM / "run_nylc_a1_step1_balanced_four_distance_bridge.sbatch"


def load_driver():
    spec = importlib.util.spec_from_file_location("_balanced_bridge", DRIVER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {DRIVER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class BalancedFourDistanceBridgeContractTest(unittest.TestCase):
    def test_six_independent_tasks_preserve_the_four_distance_bridge_contract(self):
        self.assertTrue(DRIVER.is_file(), "RED: balanced bridge driver is missing")
        self.assertTrue(RUNNER.is_file(), "RED: balanced bridge runner is missing")
        self.assertTrue(SBATCH.is_file(), "RED: balanced bridge sbatch is missing")

        driver = load_driver()
        self.assertEqual(driver.SOURCES[0]["restart_sha256"], "5d8f76d2c90e3e8c707b640c55a93938f53e18dc30d6f93d54adda467da25f41")
        self.assertEqual(driver.SOURCES[1]["restart_sha256"], "4cc60ad4d7be099b3f76040f1ee8c49b91deecebb20ed98570bc192b3d511132")
        self.assertEqual(
            [(item["attack_A"], item["cn_A"]) for item in driver.BRANCHES],
            [(1.675, 1.925), (1.800, 1.800), (1.925, 1.675)],
        )
        self.assertTrue(all(item["nalpha_hg1_A"] == item["hg1_n3_A"] == 1.30 for item in driver.BRANCHES))
        self.assertTrue(all(item["force_attack_cn"] == 15.0 and item["force_proton"] == 10.0 for item in driver.BRANCHES))
        self.assertEqual(driver.MAX_CANDIDATES_PER_SEED, 1)
        self.assertEqual(driver.MIN_BOUNDARY_RUN_FRAMES, 3)
        self.assertEqual(driver.describe()["bias_representation"], "four_independent_distance_restraints_not_a_linear_combination_cv")

        restraints = driver.bridge_restraints(driver.BRANCHES[1])
        self.assertIn("&rst iat=8960,10287", restraints)
        self.assertIn("&rst iat=10287,10289", restraints)
        self.assertIn("&rst iat=8949,8961", restraints)
        self.assertIn("&rst iat=8961,10289", restraints)

        runner = RUNNER.read_text(encoding="utf-8")
        for index in range(6):
            self.assertIn(f"    {index})", runner)
        self.assertNotIn("for WINDOW_INDEX", runner)
        self.assertIn("CANDIDATE_LIMIT_REACHED", runner)

        sbatch = SBATCH.read_text(encoding="utf-8")
        self.assertIn("#SBATCH -n 8", sbatch)
        self.assertIn("#SBATCH --array=0-5", sbatch)
        self.assertNotIn("%", sbatch)


if __name__ == "__main__":
    unittest.main()
