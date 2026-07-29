#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts/prepare_audit_nylc_a1_step1_bidirectional_2cv.py"
RUNNER = ROOT / "scripts/run_nylc_a1_step1_bidirectional_2cv.sh"
SBATCH = ROOT / "slurm/run_nylc_a1_step1_bidirectional_2cv.sbatch"


def load_driver():
    spec = importlib.util.spec_from_file_location("_a1_bidirectional_2cv", DRIVER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {DRIVER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class BidirectionalTwoCVContract(unittest.TestCase):
    def test_eight_chain_mapping_and_frozen_contract(self):
        self.assertTrue(DRIVER.is_file(), "RED: two-CV driver missing")
        self.assertTrue(RUNNER.is_file(), "RED: two-CV runner missing")
        self.assertTrue(SBATCH.is_file(), "RED: two-CV sbatch missing")
        driver = load_driver()
        contract = driver.describe()
        self.assertEqual(contract["array_task_count"], 8)
        self.assertEqual(contract["mpi_ranks_per_task"], 8)
        self.assertEqual(contract["mechanisms"], ["ADDITION_FIRST", "CONCERTED"])
        self.assertEqual(contract["source_basins"], ["R", "P"])
        self.assertEqual(
            contract["frozen_step1_contract"],
            {
                "qm_atom_count": 146,
                "qmcharge": 0,
                "electron_count": 510,
                "link_atom_count": 6,
                "qm_water_count": 0,
            },
        )
        self.assertEqual(driver.task_spec(0)["seed"], "seed26723")
        self.assertEqual(driver.task_spec(0)["source_basin"], "R")
        self.assertEqual(driver.task_spec(0)["mechanism"], "ADDITION_FIRST")
        self.assertEqual(driver.task_spec(3)["source_basin"], "P")
        self.assertEqual(driver.task_spec(3)["mechanism"], "CONCERTED")
        self.assertEqual(driver.task_spec(7)["seed"], "seed26737")
        self.assertEqual(driver.task_spec(7)["destination_basin"], "R")

    def test_path_endpoints_and_addition_first_ordering(self):
        driver = load_driver()
        for mechanism in driver.MECHANISMS:
            self.assertEqual(driver.path_targets(mechanism, 0.0), driver.R_TARGET)
            self.assertEqual(driver.path_targets(mechanism, 1.0), driver.P_TARGET)
        early = driver.path_targets("ADDITION_FIRST", 0.30)
        concerted = driver.path_targets("CONCERTED", 0.30)
        attack_progress = (
            (early["attack_A"] - driver.R_TARGET["attack_A"])
            / (driver.P_TARGET["attack_A"] - driver.R_TARGET["attack_A"])
        )
        cn_progress = (
            (early["cn_A"] - driver.R_TARGET["cn_A"])
            / (driver.P_TARGET["cn_A"] - driver.R_TARGET["cn_A"])
        )
        pt_progress = (
            (early["hg1_n3_A"] - driver.R_TARGET["hg1_n3_A"])
            / (driver.P_TARGET["hg1_n3_A"] - driver.R_TARGET["hg1_n3_A"])
        )
        self.assertGreater(attack_progress, cn_progress)
        self.assertGreater(attack_progress, pt_progress)
        self.assertAlmostEqual(
            (concerted["attack_A"] - driver.R_TARGET["attack_A"])
            / (driver.P_TARGET["attack_A"] - driver.R_TARGET["attack_A"]),
            0.30,
        )

    def test_actual_response_gate_is_directional(self):
        driver = load_driver()
        previous = {
            "attack_A": 2.10,
            "c12_n3_A": 1.40,
            "nalpha_hg1_A": 1.04,
            "hg1_n3_A": 2.40,
            "c12_o2_A": 1.26,
            "hg1_nearest_qm_heavy_atom": 8949,
        }
        toward_product = dict(
            previous,
            attack_A=1.95,
            c12_n3_A=1.50,
            nalpha_hg1_A=1.12,
            hg1_n3_A=2.28,
        )
        away_from_product = dict(
            previous,
            attack_A=2.18,
            c12_n3_A=1.34,
            nalpha_hg1_A=1.00,
            hg1_n3_A=2.47,
        )
        accepted = driver.actual_response_gate(
            driver.task_spec(1), previous, toward_product, 0.0, 0.125
        )
        rejected = driver.actual_response_gate(
            driver.task_spec(1), previous, away_from_product, 0.0, 0.125
        )
        self.assertTrue(accepted["all"])
        self.assertFalse(rejected["all"])
        self.assertTrue(rejected["wrong_direction"])

    def test_runner_inherits_only_accepted_restart_and_array_has_no_throttle(self):
        runner = RUNNER.read_text(encoding="utf-8")
        sbatch = SBATCH.read_text(encoding="utf-8")
        self.assertIn("attempt_${ARRAY_JOB}_${INDEX}", runner)
        self.assertIn("accepted.rst7", runner)
        self.assertIn("for WINDOW in", runner)
        self.assertIn("mpirun --bind-to none -np 8 sander.MPI", runner)
        self.assertIn("sha256sum", runner)
        self.assertNotIn("62021985", runner)
        self.assertIn("#SBATCH -n 8", sbatch)
        self.assertIn("#SBATCH --array=0-7", sbatch)
        self.assertNotIn("#SBATCH --array=0-7%", sbatch)
        self.assertIn("PYTHONDONTWRITEBYTECODE=1", sbatch)
        self.assertNotIn("tests/test_", sbatch)


if __name__ == "__main__":
    unittest.main()
