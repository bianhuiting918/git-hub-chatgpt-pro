#!/usr/bin/env python3
from __future__ import annotations
import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts/prepare_audit_nylc_a1_step2_water_recruitment.py"
RUNNER = ROOT / "scripts/run_nylc_a1_step2_water_recruitment.sh"
SBATCH = ROOT / "slurm/run_nylc_a1_step2_water_recruitment.sbatch"


def load_driver():
    spec = importlib.util.spec_from_file_location("_step2_water_recruitment", DRIVER)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {DRIVER}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class Step2WaterRecruitmentContract(unittest.TestCase):
    def test_recruitment_is_bounded_and_uses_final_step2_hamiltonian(self):
        self.assertTrue(DRIVER.is_file(), "RED: recruitment driver missing")
        self.assertTrue(RUNNER.is_file(), "RED: recruitment runner missing")
        self.assertTrue(SBATCH.is_file(), "RED: recruitment sbatch missing")
        driver = load_driver()
        self.assertEqual(driver.EXPECTED_NEAREST_WATER, {
            "seed26723": (13046, (13047, 13048)),
            "seed26737": (13046, (13047, 13048)),
        })
        self.assertEqual(driver.GUIDED_TARGETS_A, {
            "c12_ow": 3.00,
            "o2_ow": 3.60,
            "donor_h_nalpha": 2.00,
            "ow_nalpha": 2.87,
        })
        self.assertEqual(driver.GUIDED_FORCE_KCAL_MOL_A2, 5.0)
        self.assertEqual(driver.EXPECTED_CONTRACT, {
            "qm_atom_count": 149,
            "qmcharge": 0,
            "spin": 1,
            "link_atom_count": 6,
            "electron_count": 518,
        })
        self.assertEqual(driver.describe()["array_task_count"], 4)
        self.assertEqual(driver.describe()["release_legs_per_task"], 2)
        self.assertTrue(driver.describe()["guided_restraints_removed_for_release"])
        text = driver.guided_restraints(13046, 13047)
        for pair in ("10287,13046", "10288,13046", "13047,8949", "13046,8949"):
            self.assertIn(f"&rst iat={pair}", text)
        runner = RUNNER.read_text(encoding="utf-8")
        self.assertIn('SCRATCH_ROOT="${SLURM_TMPDIR:-/tmp}/nylc_a1_step2_water_recruit_', runner)
        self.assertIn("srun --exclusive -N 1 -n 8 sander.MPI", runner)
        self.assertIn("for LEG in 0 1", runner)
        sbatch = SBATCH.read_text(encoding="utf-8")
        self.assertIn("#SBATCH -n 16", sbatch)
        self.assertIn("#SBATCH --array=0-3", sbatch)
        self.assertNotIn("#SBATCH --array=0-3%", sbatch)


if __name__ == "__main__":
    unittest.main()
