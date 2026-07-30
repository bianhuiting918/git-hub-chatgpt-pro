#!/usr/bin/env python3
from __future__ import annotations

import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts/prepare_audit_nylc_a1_step1_tetrahedral_bridge.py"
RUNNER = ROOT / "scripts/run_nylc_a1_step1_tetrahedral_bridge.sh"
SBATCH = ROOT / "slurm/run_nylc_a1_step1_tetrahedral_bridge.sbatch"


def load_driver():
    spec = importlib.util.spec_from_file_location("_a1_tetra_bridge", DRIVER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {DRIVER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class TetrahedralBridgeContract(unittest.TestCase):
    def test_three_sha_fixed_sources_two_schedules_and_carbonyl_coordinate(self):
        self.assertTrue(DRIVER.is_file(), "RED: tetrahedral bridge driver missing")
        self.assertTrue(RUNNER.is_file(), "RED: tetrahedral bridge runner missing")
        self.assertTrue(SBATCH.is_file(), "RED: tetrahedral bridge sbatch missing")
        driver = load_driver()
        contract = driver.describe()
        self.assertEqual(contract["array_task_count"], 6)
        self.assertEqual(contract["source_task_indices"], [1, 3, 5])
        self.assertEqual(contract["schedules"], ["CARBONYL_FIRST", "PT_ASSISTED"])
        self.assertEqual(contract["windows_per_chain"], 8)
        self.assertTrue(contract["explicit_carbonyl_restraint"])
        self.assertEqual(contract["tetrahedral_target"]["c12_o2_A"], 1.38)
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
        self.assertEqual(driver.task_spec(0)["source_task_index"], 1)
        self.assertEqual(driver.task_spec(0)["schedule"], "CARBONYL_FIRST")
        self.assertEqual(driver.task_spec(1)["schedule"], "PT_ASSISTED")
        self.assertEqual(driver.task_spec(5)["source_task_index"], 5)
        self.assertEqual(
            driver.task_spec(0)["source_restart_sha256"],
            "1ed5a5b6a72a43f9f56e41bf427bb95ffc96237c86811a34aec50fa3c89668f1",
        )
        self.assertEqual(
            driver.task_spec(2)["source_restart_sha256"],
            "09e22fb16670fc5fc6427a1aacf80851b97c323fa4eefb900760bde8ebc3be26",
        )
        self.assertEqual(
            driver.task_spec(4)["source_restart_sha256"],
            "740f5fb11de92842af2786ade6af03856dff40351209a613355a9b4dd19b3782",
        )
        text = DRIVER.read_text(encoding="utf-8")
        self.assertIn('REACTIVE["o2"]', text)
        self.assertIn('"c12_o2_A"', text)
        self.assertIn("distance_restraint", text)

    def test_required_response_rejects_opposite_direction_near_target(self):
        driver = load_driver()
        response = driver._response(
            previous=1.2302340709352817,
            current=1.2291623894039534,
            target=1.2676755532014612,
        )
        self.assertTrue(response["required"])
        self.assertFalse(response["same_direction"])
        self.assertFalse(response["pass"])

    def test_tetrahedral_restraints_use_narrow_reaction_coordinate_wells(self):
        driver = load_driver()
        stage = {
            "targets": {
                "attack_A": 1.48,
                "c12_n3_A": 1.62,
                "c12_o2_A": 1.38,
                "nalpha_hg1_A": 1.05,
                "hg1_n3_A": 2.20,
            },
            "force_bond": 24.0,
            "force_carbonyl": 30.0,
            "force_pt": 18.0,
            "proton_coordinate_active": False,
        }
        text = driver.restraints(stage)
        blocks = [block for block in text.split("&rst ") if block.strip()]
        self.assertEqual(len(blocks), 3)
        self.assertIn("r2=1.475, r3=1.485", blocks[0])
        self.assertIn("r2=1.615, r3=1.625", blocks[1])
        self.assertIn("iat=10287,10288", blocks[2])
        self.assertIn("r2=1.375, r3=1.385", blocks[2])

    def test_strict_inheritance_parallel_array_and_no_automatic_downstream(self):
        runner = RUNNER.read_text(encoding="utf-8")
        sbatch = SBATCH.read_text(encoding="utf-8")
        self.assertIn("accepted.rst7", runner)
        self.assertIn("attempt_${ARRAY_JOB}_${INDEX}", runner)
        self.assertIn("mpirun --bind-to none -np 8 sander.MPI", runner)
        self.assertNotIn("62021985", runner)
        self.assertIn("#SBATCH -n 8", sbatch)
        self.assertIn("#SBATCH --array=0-5", sbatch)
        self.assertNotIn("#SBATCH --array=0-5%", sbatch)
        self.assertNotIn("release", sbatch.lower())


if __name__ == "__main__":
    unittest.main()
