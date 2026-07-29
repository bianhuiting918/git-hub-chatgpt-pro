#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts/prepare_audit_nylc_a1_step1_i1_i2_response_builder.py"
RUNNER = ROOT / "scripts/run_nylc_a1_step1_i1_i2_response_builder.sh"
SBATCH = ROOT / "slurm/run_nylc_a1_step1_i1_i2_response_builder.sbatch"


def load_driver():
    spec = importlib.util.spec_from_file_location("_i1_i2_response_builder", DRIVER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {DRIVER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class I1I2ResponseBuilderContract(unittest.TestCase):
    def test_actual_response_adaptive_contract(self):
        self.assertTrue(DRIVER.is_file(), "RED: I1/I2 response driver missing")
        driver = load_driver()
        self.assertEqual(len(driver.SOURCE_SPECS), 4)
        self.assertEqual([item["anchor"] for item in driver.SOURCE_SPECS], ["I1", "I2", "I1", "I2"])
        self.assertEqual(driver.step_toward(2.35, 1.60, 0.10), 2.25)
        self.assertEqual(driver.step_toward(1.40, 1.60, 0.10), 1.50)
        self.assertEqual(driver.step_toward(1.58, 1.60, 0.10), 1.60)
        contract = driver.describe()
        self.assertEqual(contract["array_task_count"], 4)
        self.assertEqual(contract["mpi_ranks_per_task"], 8)
        self.assertEqual(contract["maximum_windows"], 12)
        self.assertEqual(contract["minimum_step_scale"], 0.25)
        self.assertTrue(contract["actual_target_residual_is_required"])
        self.assertTrue(contract["failed_restart_inheritance"] is False)
        self.assertEqual(contract["automatic_downstream_action"], "NONE")

    def test_serial_inheritance_and_parallel_array(self):
        self.assertTrue(RUNNER.is_file(), "RED: I1/I2 response runner missing")
        self.assertTrue(SBATCH.is_file(), "RED: I1/I2 response sbatch missing")
        runner = RUNNER.read_text(encoding="utf-8")
        self.assertIn("accepted.rst7", runner)
        self.assertIn("for WINDOW in 0 1 2 3 4 5 6 7 8 9 10 11", runner)
        self.assertIn("mpirun --bind-to none -np 8 sander.MPI", runner)
        self.assertIn("accepted_for_inheritance", runner)
        self.assertNotIn("62021985", runner)
        sbatch = SBATCH.read_text(encoding="utf-8")
        self.assertIn("#SBATCH --array=0-3", sbatch)
        self.assertNotIn("#SBATCH --array=0-3%", sbatch)
        self.assertIn("#SBATCH -n 8", sbatch)
        marker = "PYTHONPYCACHEPREFIX"
        self.assertIn(marker, sbatch)
        self.assertLess(sbatch.index(marker), sbatch.index("test_nylc_a1_step1_i1_i2_response_builder.py"))


if __name__ == "__main__":
    unittest.main()
