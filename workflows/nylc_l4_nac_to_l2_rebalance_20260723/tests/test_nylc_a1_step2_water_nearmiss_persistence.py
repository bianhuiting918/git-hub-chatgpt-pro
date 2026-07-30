#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts/prepare_audit_nylc_a1_step2_water_nearmiss_persistence.py"
RUNNER = ROOT / "scripts/run_nylc_a1_step2_water_nearmiss_persistence.sh"
SBATCH = ROOT / "slurm/run_nylc_a1_step2_water_nearmiss_persistence.sbatch"


def load_driver():
    spec = importlib.util.spec_from_file_location("_a1_step2_water_persist", DRIVER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {DRIVER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Step2WaterNearMissPersistenceContract(unittest.TestCase):
    def test_persistence_contract_and_sources(self):
        self.assertTrue(DRIVER.is_file(), "RED: persistence driver missing")
        driver = load_driver()
        contract = driver.describe()
        self.assertEqual(contract["seed_count"], 2)
        self.assertEqual(contract["replicas_per_seed"], 4)
        self.assertEqual(contract["persistence_task_count"], 8)
        self.assertEqual(contract["combined_array_task_count"], 9)
        self.assertEqual(contract["mpi_ranks_per_task"], 8)
        self.assertEqual(contract["steps"], 4000)
        self.assertEqual(contract["duration_ps"], 2.0)
        self.assertFalse(contract["reactive_restraints"])
        self.assertFalse(contract["water_position_restraints"])
        self.assertTrue(contract["scan_all_complete_waters"])
        self.assertEqual(contract["minimum_consecutive_frames"], 3)
        self.assertEqual(
            contract["source_restart_sha256"]["seed26723"],
            "a65834ae6c0f82291fd0b8c38ece41983b501b2444deb98de1794b92d79622a7",
        )
        self.assertEqual(
            contract["source_restart_sha256"]["seed26737"],
            "34ad5af01f5a2cd7e7ce0342590ab66f201902a31fd49236546441afa82f5ea6",
        )
        self.assertEqual(driver.task_spec(0)["seed"], "seed26723")
        self.assertEqual(driver.task_spec(3)["replica"], 3)
        self.assertEqual(driver.task_spec(4)["seed"], "seed26737")
        self.assertEqual(driver.task_spec(7)["replica"], 3)
        self.assertEqual(len({driver.task_spec(i)["velocity_seed"] for i in range(8)}), 8)

    def test_combined_runner_is_unrestrained_and_task0_only_repair(self):
        self.assertTrue(RUNNER.is_file(), "RED: persistence runner missing")
        self.assertTrue(SBATCH.is_file(), "RED: persistence sbatch missing")
        driver = load_driver()
        mdin = driver.sampling_input("seed26723", 0, 26723101, "@1,2,3")
        self.assertIn("nstlim=4000", mdin)
        self.assertIn("ifqnt=1", mdin)
        self.assertNotIn("DISANG", mdin)
        self.assertNotIn("restraintmask", mdin)
        runner = RUNNER.read_text(encoding="utf-8")
        self.assertIn("ARRAY_INDEX == 0", runner)
        self.assertIn("run_nylc_a1_step2_water_reorganization_sampling.sh", runner)
        self.assertIn("PERSIST_INDEX=$((ARRAY_INDEX - 1))", runner)
        self.assertIn("mpirun --bind-to none -np 8 sander.MPI", runner)
        self.assertIn("attempt_${ARRAY_JOB}_${PERSIST_INDEX}", runner)
        self.assertNotIn('cp "$SCRATCH_ROOT/stage.nc"', runner)
        sbatch = SBATCH.read_text(encoding="utf-8")
        self.assertIn("#SBATCH --array=0-8", sbatch)
        self.assertNotIn("#SBATCH --array=0-8%", sbatch)
        self.assertIn("#SBATCH -n 8", sbatch)
        self.assertIn("PYTHONPYCACHEPREFIX", sbatch)
        self.assertNotIn("test_nylc_a1_step2_water_nearmiss_persistence.py", sbatch)


if __name__ == "__main__":
    unittest.main()
