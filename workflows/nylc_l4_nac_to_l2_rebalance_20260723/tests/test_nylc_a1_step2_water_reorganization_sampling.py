#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts/prepare_audit_nylc_a1_step2_water_reorganization_sampling.py"
RUNNER = ROOT / "scripts/run_nylc_a1_step2_water_reorganization_sampling.sh"
SBATCH = ROOT / "slurm/run_nylc_a1_step2_water_reorganization_sampling.sbatch"


def load_driver():
    spec = importlib.util.spec_from_file_location("_a1_step2_water_reorg", DRIVER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {DRIVER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Step2WaterReorganizationSamplingContract(unittest.TestCase):
    def test_sampling_contract_and_array_mapping(self):
        self.assertTrue(DRIVER.is_file(), "RED: sampling driver missing")
        self.assertTrue(RUNNER.is_file(), "RED: sampling runner missing")
        self.assertTrue(SBATCH.is_file(), "RED: sampling sbatch missing")
        driver = load_driver()
        contract = driver.describe()
        self.assertEqual(contract["seed_count"], 2)
        self.assertEqual(contract["replicas_per_seed"], 4)
        self.assertEqual(contract["array_task_count"], 8)
        self.assertEqual(contract["mpi_ranks_per_task"], 8)
        self.assertEqual(contract["steps"], 4000)
        self.assertEqual(contract["dt_ps"], 0.0005)
        self.assertEqual(contract["duration_ps"], 2.0)
        self.assertEqual(contract["step1_qm_contract"], {
            "qm_atom_count": 146,
            "qmcharge": 0,
            "electron_count": 510,
            "link_atom_count": 6,
            "qm_water_count": 0,
        })
        self.assertFalse(contract["reactive_restraints"])
        self.assertFalse(contract["water_position_restraints"])
        self.assertTrue(contract["scan_all_complete_waters"])
        self.assertTrue(contract["triclinic_minimum_image"])
        self.assertEqual(contract["minimum_consecutive_frames"], 3)
        self.assertEqual(contract["maximum_events_per_seed"], 2)
        self.assertEqual(driver.task_spec(0)["seed"], "seed26723")
        self.assertEqual(driver.task_spec(3)["replica"], 3)
        self.assertEqual(driver.task_spec(4)["seed"], "seed26737")
        self.assertEqual(driver.task_spec(7)["replica"], 3)
        self.assertEqual(len({driver.task_spec(i)["velocity_seed"] for i in range(8)}), 8)

    def test_engine_contract_accepts_frozen_authority_dftb_banner(self):
        driver = load_driver()
        observed = {
            "qm_atom_count": [146],
            "qmcharge": [0],
            "spin": [1],
            "link_atom_count": [6],
            "dftb_doubly_occupied_levels": [190],
        }
        self.assertTrue(
            driver._engine_contract_pass(observed),
            "Amber18/DFTB3 authority job 62011285 reports 190 occupied valence levels",
        )

    def test_hit_thresholds_and_temporal_collapse(self):
        driver = load_driver()
        passing = {
            "c12_ow_A": 3.10,
            "o2_c12_ow_deg": 110.0,
            "h_nalpha_A": 2.20,
            "ow_h_nalpha_deg": 155.0,
            "acyl_guard_pass": True,
            "proton_guard_pass": True,
        }
        self.assertTrue(driver.strict_hit(passing))
        for key, value in (
            ("c12_ow_A", 3.60),
            ("o2_c12_ow_deg", 130.0),
            ("h_nalpha_A", 2.60),
            ("ow_h_nalpha_deg", 120.0),
        ):
            row = dict(passing)
            row[key] = value
            self.assertFalse(driver.strict_hit(row), key)
        rows = []
        for frame in range(3):
            row = dict(passing, frame=frame, time_ps=frame * 0.01,
                       water_oxygen_index1=13046, donor_h_index1=13047)
            rows.append(row)
        events = driver.collapse_events(rows)
        self.assertEqual(len(events), 1)
        self.assertEqual(events[0]["consecutive_frames"], 3)
        self.assertEqual(events[0]["water_oxygen_index1"], 13046)
        self.assertEqual(driver.collapse_events(rows[:2]), [])

    def test_near_miss_selection_is_relaxed_ranked_and_water_unique(self):
        driver = load_driver()
        base = {
            "c12_ow_A": 3.10,
            "o2_c12_ow_deg": 110.0,
            "h_nalpha_A": 2.20,
            "ow_h_nalpha_deg": 155.0,
            "acyl_guard_pass": True,
            "proton_guard_pass": True,
            "replica": 0,
            "water_hydrogen_indices1": [13047, 13048],
        }
        rows = [
            dict(base, frame=1, water_oxygen_index1=13046, donor_h_index1=13047),
            dict(base, frame=2, c12_ow_A=3.40, water_oxygen_index1=13046, donor_h_index1=13047),
            dict(base, frame=3, c12_ow_A=4.20, o2_c12_ow_deg=90.0,
                 h_nalpha_A=3.20, ow_h_nalpha_deg=105.0,
                 water_oxygen_index1=14000, donor_h_index1=14001),
            dict(base, frame=4, c12_ow_A=5.10,
                 water_oxygen_index1=15000, donor_h_index1=15001),
        ]
        selected = driver.select_near_misses(rows, limit=4)
        self.assertEqual([row["water_oxygen_index1"] for row in selected], [13046, 14000])
        self.assertEqual(selected[0]["frame"], 1)
        self.assertLessEqual(selected[0]["near_miss_score"], selected[1]["near_miss_score"])

    def test_mdin_runner_and_sbatch_are_unrestrained_and_compact(self):
        driver = load_driver()
        mdin = driver.sampling_input("seed26723", 0, 26723011, "@1,2,3")
        self.assertIn("nstlim=4000", mdin)
        self.assertIn("dt=0.0005", mdin)
        self.assertIn("ifqnt=1", mdin)
        self.assertNotIn("DISANG", mdin)
        self.assertNotIn("nmropt=1", mdin)
        self.assertNotIn("restraintmask", mdin)
        runner = RUNNER.read_text(encoding="utf-8")
        self.assertIn("attempt_${ARRAY_JOB}_${INDEX}", runner)
        self.assertIn('SCRATCH_ROOT="${SLURM_TMPDIR:-/tmp}/nylc_a1_step2_water_reorg_', runner)
        self.assertIn("mpirun --bind-to none -np 8 sander.MPI", runner)
        self.assertIn("sha256sum", runner)
        self.assertIn("flock -x 9", runner)
        self.assertNotIn('cp "$SCRATCH/stage.nc"', runner)
        sbatch = SBATCH.read_text(encoding="utf-8")
        self.assertIn("#SBATCH -n 8", sbatch)
        self.assertIn("#SBATCH --array=0-7", sbatch)
        self.assertNotIn("#SBATCH --array=0-7%", sbatch)
        self.assertIn("PYTHONPYCACHEPREFIX", sbatch)
        self.assertNotIn("tests/test_nylc_a1_step2_water_reorganization_sampling.py", sbatch)
        self.assertNotIn("-m py_compile", sbatch)


if __name__ == "__main__":
    unittest.main()
