#!/usr/bin/env python3
from __future__ import annotations
import importlib.util
import pathlib
import tempfile
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
        self.assertEqual(driver.GUIDED_FORCE_KCAL_MOL_A2, 20.0)
        self.assertEqual(driver.GUIDED_MAXCYC, 1500)
        self.assertEqual(driver.GUIDED_NCYC, 450)
        self.assertEqual(driver.format_qmmask([8949, 8950, 13046]), "@8949,8950,13046")
        self.assertNotIn(",@", driver.format_qmmask([8949, 8950, 13046]))
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
        self.assertIn("mpirun --bind-to none -np 8 sander.MPI", runner)
        self.assertNotIn("srun --exclusive", runner)
        self.assertIn('cp "$GUIDED/stage.out" "$OUT/GUIDED_ENGINE.out"', runner)
        self.assertIn('cp "$GUIDED/stage.mdinfo" "$OUT/GUIDED_ENGINE.mdinfo"', runner)
        self.assertIn('cp "$GUIDED/stage.in" "$OUT/GUIDED_ENGINE.in"', runner)
        self.assertIn("for LEG in 0 1", runner)
        sbatch = SBATCH.read_text(encoding="utf-8")
        self.assertIn("#SBATCH -n 16", sbatch)
        self.assertIn("#SBATCH --array=0-3", sbatch)
        self.assertNotIn("#SBATCH --array=0-3%", sbatch)


    def test_release_binary_format_and_same_hamiltonian_banner_authority(self):
        driver = load_driver()
        with tempfile.TemporaryDirectory() as tmp:
            root = pathlib.Path(tmp)
            netcdf = root / "a2.nc"
            ascii_mdcrd = root / "a2.mdcrd"
            netcdf.write_bytes(b"CDF" + bytes([2, 0, 0, 0]))
            ascii_mdcrd.write_text("TITLE\\n", encoding="ascii")
            self.assertEqual(driver.S2.trajectory_format(netcdf), "NETCDF")
            self.assertEqual(driver.S2.trajectory_format(ascii_mdcrd), "AMBER_MDCRD")

        guided = {
            "banner_contract_pass": True,
            "banner_observed": {
                "qm_atom_count": [149],
                "qmcharge": [0],
                "spin": [1],
                "link_atom_count": [6],
                "dftb_doubly_occupied_levels": [194],
                "dftb_valence_electron_count": [388],
            },
            "derived_all_electron_count": 518,
            "restart_sha256": "guided-restart-sha",
        }
        manifest = {
            "qm_contract": {
                "expected": dict(driver.EXPECTED_CONTRACT),
                "qmmask": "@" + ",".join(str(index) for index in range(1, 150)),
            },
            "recruitment": {"guided_result": guided},
        }
        prepared = {
            "input_restart_sha256": "guided-restart-sha",
            "expected_contract": dict(driver.EXPECTED_CONTRACT),
        }
        incomplete_leg_banner = {
            "qm_atom_count": [],
            "qmcharge": [],
            "link_atom_count": [6],
            "electron_count": [],
        }
        effective, source = driver.S2.resolve_leg_banner_contract(
            manifest, prepared, incomplete_leg_banner
        )
        self.assertEqual(source, "GUIDED_ENGINE_SAME_HAMILTONIAN")
        self.assertTrue(driver.S2.banner_pass(effective))
        prepared["input_restart_sha256"] = "different-restart-sha"
        effective, source = driver.S2.resolve_leg_banner_contract(
            manifest, prepared, incomplete_leg_banner
        )
        self.assertEqual(source, "UNVERIFIED")
        self.assertFalse(driver.S2.banner_pass(effective))


    def test_amber18_dftb_banner_contract_uses_engine_and_derived_evidence(self):
        driver = load_driver()
        amber = """
QMMM options:
             ifqnt = True       nquant =      149
              qmgb =        0  qmcharge =        0   adjust_q =        2
              spin =        1     qmcut =  10.0000
QMMM:  nlink =     6                   Link Coords
QMMM: SINGLET STATE CALCULATION
QMMM: RHF CALCULATION, NO. OF DOUBLY OCCUPIED LEVELS =194
   NSTEP       ENERGY          RMS
"""
        observed = driver.parse_recruitment_engine_contract(amber)
        self.assertEqual(observed["qm_atom_count"], [149])
        self.assertEqual(observed["qmcharge"], [0])
        self.assertEqual(observed["spin"], [1])
        self.assertEqual(observed["link_atom_count"], [6])
        self.assertEqual(observed["dftb_doubly_occupied_levels"], [194])
        self.assertEqual(observed["dftb_valence_electron_count"], [388])
        self.assertTrue(driver.recruitment_engine_contract_pass(observed, derived_all_electron_count=518))
        self.assertFalse(driver.recruitment_engine_contract_pass(observed, derived_all_electron_count=510))


if __name__ == "__main__":
    unittest.main()
