#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts/prepare_audit_nylc_a1_step1_2cv_continuation.py"
RUNNER = ROOT / "scripts/run_nylc_a1_step1_2cv_continuation.sh"
SBATCH = ROOT / "slurm/run_nylc_a1_step1_2cv_continuation.sbatch"


def load_driver():
    spec = importlib.util.spec_from_file_location("_a1_2cv_continuation", DRIVER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {DRIVER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TwoCVContinuationContract(unittest.TestCase):
    def test_four_anchors_two_force_tiers_and_fixed_sources(self):
        self.assertTrue(DRIVER.is_file(), "RED: continuation driver missing")
        self.assertTrue(RUNNER.is_file(), "RED: continuation runner missing")
        self.assertTrue(SBATCH.is_file(), "RED: continuation sbatch missing")
        driver = load_driver()
        contract = driver.describe()
        self.assertEqual(contract["array_task_count"], 8)
        self.assertEqual(contract["anchor_task_indices"], [2, 3, 6, 7])
        self.assertEqual(contract["force_tiers"], ["MODERATE", "FIRM"])
        self.assertEqual(contract["minimum_lambda_step"], 1.0 / 256.0)
        self.assertEqual(contract["maximum_windows_per_chain"], 12)
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
        first = driver.task_spec(0)
        last = driver.task_spec(7)
        self.assertEqual((first["source_task_index"], first["force_tier"]), (2, "MODERATE"))
        self.assertEqual((last["source_task_index"], last["force_tier"]), (7, "FIRM"))
        self.assertEqual(first["source_fraction"], 0.75)
        self.assertEqual(last["source_fraction"], 0.625)
        self.assertEqual(
            first["source_restart_sha256"],
            "cda206065902568f1406fa6f227f2d924b581a0140b12a2c0ba1bc639a037823",
        )
        self.assertEqual(
            last["source_restart_sha256"],
            "a8037eb79eed6e3007f6ca14b51f5bcfe4d25459a306a77510562a22de838626",
        )

    def test_runner_keeps_strict_inheritance_and_parallel_array(self):
        runner = RUNNER.read_text(encoding="utf-8")
        sbatch = SBATCH.read_text(encoding="utf-8")
        self.assertIn("attempt_${ARRAY_JOB}_${INDEX}", runner)
        self.assertIn("accepted.rst7", runner)
        self.assertIn("mpirun --bind-to none -np 8 sander.MPI", runner)
        self.assertIn("PYTHONDONTWRITEBYTECODE=1", runner)
        self.assertNotIn("62021985", runner)
        self.assertIn("#SBATCH -n 8", sbatch)
        self.assertIn("#SBATCH --array=0-7", sbatch)
        self.assertNotIn("#SBATCH --array=0-7%", sbatch)


if __name__ == "__main__":
    unittest.main()
