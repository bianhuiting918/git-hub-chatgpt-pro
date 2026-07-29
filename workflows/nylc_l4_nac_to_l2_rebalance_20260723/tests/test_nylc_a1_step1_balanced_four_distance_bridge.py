#!/usr/bin/env python3
"""Contract test for the NylC A1 Step1 balanced four-distance bridge scout."""

from __future__ import annotations
import importlib.util
import json
import pathlib
import unittest

WORKFLOW_ROOT = pathlib.Path(__file__).resolve().parents[1]
SCRIPTS, SLURM = WORKFLOW_ROOT / "scripts", WORKFLOW_ROOT / "slurm"
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
        source = driver.serializable_source(driver.source_from_index(0))
        json.dumps(source)
        self.assertIsInstance(source["restart"], str)
        self.assertEqual(driver.SOURCES[0]["restart_sha256"], "5d8f76d2c90e3e8c707b640c55a93938f53e18dc30d6f93d54adda467da25f41")
        self.assertEqual(driver.SOURCES[1]["restart_sha256"], "4cc60ad4d7be099b3f76040f1ee8c49b91deecebb20ed98570bc192b3d511132")
        self.assertEqual([(x["attack_A"], x["cn_A"]) for x in driver.BRANCHES], [(1.675,1.925),(1.800,1.800),(1.925,1.675)])
        self.assertTrue(all(x["nalpha_hg1_A"] == x["hg1_n3_A"] == 1.30 for x in driver.BRANCHES))
        self.assertTrue(all(x["force_attack_cn"] == 15.0 and x["force_proton"] == 10.0 for x in driver.BRANCHES))
        self.assertEqual((driver.MAX_CANDIDATES_PER_SEED, driver.MIN_BOUNDARY_RUN_FRAMES), (1,3))
        self.assertEqual(driver.describe()["bias_representation"], "four_independent_distance_restraints_not_a_linear_combination_cv")
        restraints = driver.bridge_restraints(driver.BRANCHES[1])
        for pair in ("8960,10287", "10287,10289", "8949,8961", "8961,10289"):
            self.assertIn(f"&rst iat={pair}", restraints)
        runner = RUNNER.read_text(encoding="utf-8")
        for index in range(6): self.assertIn(f"    {index})", runner)
        self.assertNotIn("for WINDOW_INDEX", runner)
        self.assertIn("CANDIDATE_LIMIT_REACHED", runner)
        self.assertIn(
            'SCRATCH=${SLURM_TMPDIR:-/tmp}/nylc_a1_balanced_bridge_',
            runner,
        )
        self.assertNotIn('SCRATCH=$SLURM_TMPDIR/', runner)
        sbatch = SBATCH.read_text(encoding="utf-8")
        self.assertIn("#SBATCH -n 8", sbatch)
        self.assertIn("#SBATCH --array=0-5", sbatch)
        self.assertNotIn("#SBATCH --array=0-5%", sbatch)

if __name__ == "__main__":
    unittest.main()
