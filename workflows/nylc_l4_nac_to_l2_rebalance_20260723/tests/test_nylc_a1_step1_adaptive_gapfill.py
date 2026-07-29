#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts/prepare_audit_nylc_a1_step1_adaptive_gapfill.py"
RUNNER = ROOT / "scripts/run_nylc_a1_step1_adaptive_gapfill.sh"
SBATCH = ROOT / "slurm/run_nylc_a1_step1_adaptive_gapfill.sbatch"


def load_driver():
    spec = importlib.util.spec_from_file_location("_a1_adaptive_gapfill", DRIVER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {DRIVER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class AdaptiveGapfillContract(unittest.TestCase):
    def test_mapping_and_midpoint_contract(self):
        self.assertTrue(DRIVER.is_file(), "RED: adaptive driver missing")
        self.assertTrue(RUNNER.is_file(), "RED: adaptive runner missing")
        self.assertTrue(SBATCH.is_file(), "RED: adaptive sbatch missing")
        driver = load_driver()
        self.assertEqual(driver.SOURCE_TASKS, (2, 3, 4, 5, 10, 11, 12, 13))
        contract = driver.describe()
        self.assertEqual(contract["array_task_count"], 8)
        self.assertEqual(contract["mpi_ranks_per_task"], 8)
        self.assertEqual(contract["maximum_new_windows"], 8)
        self.assertEqual(contract["minimum_bracket_fraction"], 1 / 16)
        self.assertEqual(contract["source_job"], "62267272")
        self.assertEqual(driver.task_spec(0)["source_task_index"], 2)
        self.assertEqual(driver.task_spec(7)["source_task_index"], 13)
        self.assertEqual(
            driver.midpoint({"attack": 2.0, "cn": 1.5}, {"attack": 1.8, "cn": 1.7}),
            {"attack": 1.9, "cn": 1.6},
        )
        self.assertEqual(driver.next_bracket(0.0, 1.0, True), (0.5, 1.0))
        self.assertEqual(driver.next_bracket(0.0, 1.0, False), (0.0, 0.5))

    def test_runner_is_serially_inherited_and_array_is_parallel(self):
        runner = RUNNER.read_text(encoding="utf-8")
        self.assertIn("attempt_${ARRAY_JOB}_${INDEX}", runner)
        self.assertIn("accepted.rst7", runner)
        self.assertIn('for WINDOW in 0 1 2 3 4 5 6 7', runner)
        self.assertIn("mpirun --bind-to none -np 8 sander.MPI", runner)
        self.assertIn("sha256sum", runner)
        self.assertIn("flock -x 9", runner)
        self.assertNotIn("62021985", runner)
        sbatch = SBATCH.read_text(encoding="utf-8")
        self.assertIn("#SBATCH -n 8", sbatch)
        self.assertIn("#SBATCH --array=0-7", sbatch)
        self.assertNotIn("#SBATCH --array=0-7%", sbatch)

    def test_python_checks_use_task_local_pycache(self):
        sbatch = SBATCH.read_text(encoding="utf-8")
        marker = 'export PYTHONPYCACHEPREFIX="${SLURM_TMPDIR:-/tmp}/nylc_a1_adaptive_pycache_${SLURM_ARRAY_JOB_ID:-manual}_${SLURM_ARRAY_TASK_ID:-manual}"'
        self.assertIn(marker, sbatch)
        self.assertLess(sbatch.index(marker), sbatch.index('"$PY" tests/test_nylc_a1_step1_adaptive_gapfill.py'))
        self.assertLess(sbatch.index(marker), sbatch.index('"$PY" -m py_compile'))


if __name__ == "__main__":
    unittest.main()
