#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts/prepare_audit_nylc_a1_step2_water_network_sampling.py"
RUNNER = ROOT / "scripts/run_nylc_a1_step2_water_network_sampling.sh"
SBATCH = ROOT / "slurm/run_nylc_a1_step2_water_network_sampling.sbatch"


def load_driver():
    spec = importlib.util.spec_from_file_location("_a1_step2_water_network", DRIVER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {DRIVER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Step2WaterNetworkSamplingContract(unittest.TestCase):
    def test_parallel_sampling_contract(self):
        self.assertTrue(DRIVER.is_file(), "RED: water-network driver missing")
        self.assertTrue(RUNNER.is_file(), "RED: water-network runner missing")
        self.assertTrue(SBATCH.is_file(), "RED: water-network sbatch missing")
        driver = load_driver()
        contract = driver.describe()
        self.assertEqual(contract["seed_count"], 2)
        self.assertEqual(contract["replicas_per_seed"], 4)
        self.assertEqual(contract["array_task_count"], 8)
        self.assertEqual(contract["mpi_ranks_per_task"], 8)
        self.assertEqual(contract["duration_ps"], 5.0)
        self.assertEqual(contract["expected_frames"], 500)
        self.assertEqual(contract["minimum_consecutive_frames"], 5)
        self.assertTrue(contract["scan_all_complete_waters"])
        self.assertTrue(contract["scan_direct_and_two_water_relay"])
        self.assertFalse(contract["water_identity_restraints"])
        self.assertFalse(contract["water_position_restraints"])
        self.assertFalse(contract["reactive_restraints"])
        self.assertEqual(contract["step1_qm_contract"], {
            "qm_atom_count": 146,
            "qmcharge": 0,
            "electron_count": 510,
            "link_atom_count": 6,
            "qm_water_count": 0,
        })

    def test_direct_and_two_water_relay_geometry_are_separate(self):
        driver = load_driver()
        attack = {
            "c12_ow_A": 3.05,
            "o2_c12_ow_deg": 110.0,
            "h_nalpha_A": 2.10,
            "ow_h_nalpha_deg": 155.0,
            "acyl_guard_pass": True,
            "proton_guard_pass": True,
        }
        self.assertTrue(driver.direct_water_candidate(attack))
        relay = {
            "c12_ow_A": 3.05,
            "o2_c12_ow_deg": 110.0,
            "attack_h_relay_o_A": 1.90,
            "attack_ow_attack_h_relay_o_deg": 160.0,
            "relay_h_nalpha_A": 2.10,
            "relay_ow_relay_h_nalpha_deg": 155.0,
            "acyl_guard_pass": True,
            "proton_guard_pass": True,
        }
        self.assertTrue(driver.relay_water_candidate(relay))
        relay["relay_h_nalpha_A"] = 3.20
        self.assertFalse(driver.relay_water_candidate(relay))
        self.assertNotEqual(
            driver.describe()["direct_filter"],
            driver.describe()["relay_filter"],
        )

    def test_runner_is_parallel_unrestrained_and_compact(self):
        driver = load_driver()
        mdin = driver.sampling_input("seed26723", 0, 26723201, "@1,2,3")
        self.assertIn("nstlim=10000", mdin)
        self.assertIn("dt=0.0005", mdin)
        self.assertIn("ntwx=20", mdin)
        self.assertIn("ifqnt=1", mdin)
        self.assertNotIn("DISANG", mdin)
        self.assertNotIn("nmropt=1", mdin)
        runner = RUNNER.read_text(encoding="utf-8")
        sbatch = SBATCH.read_text(encoding="utf-8")
        self.assertIn("attempt_${ARRAY_JOB}_${INDEX}", runner)
        self.assertIn("mpirun --bind-to none -np 8 sander.MPI", runner)
        self.assertNotIn('cp "$SCRATCH/stage.nc"', runner)
        self.assertIn("#SBATCH -n 8", sbatch)
        self.assertIn("#SBATCH --array=0-7", sbatch)
        self.assertNotIn("#SBATCH --array=0-7%", sbatch)


if __name__ == "__main__":
    unittest.main()
