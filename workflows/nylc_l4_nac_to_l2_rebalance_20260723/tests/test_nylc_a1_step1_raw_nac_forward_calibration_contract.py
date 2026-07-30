#!/usr/bin/env python3
import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER_PATH = ROOT / "scripts" / "prepare_audit_nylc_a1_step1_raw_nac_forward_calibration.py"


def load_driver():
    spec = importlib.util.spec_from_file_location("raw_nac_forward", DRIVER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class RawNacForwardCalibrationContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_driver() if DRIVER_PATH.exists() else None

    def test_driver_exists(self):
        self.assertTrue(DRIVER_PATH.exists(), f"missing production driver: {DRIVER_PATH}")

    @unittest.skipUnless(DRIVER_PATH.exists(), "production driver not implemented yet")
    def test_matrix_is_two_seeds_times_baseline_and_two_mechanisms_four_scales(self):
        self.assertEqual(self.mod.ARRAY_TASKS, 18)
        self.assertEqual(tuple(self.mod.FORCE_SCALES), (1, 2, 4, 8))
        specs = [self.mod.task_spec(i) for i in range(self.mod.ARRAY_TASKS)]
        for seed in (26723, 26737):
            rows = [row for row in specs if row["seed"] == seed]
            self.assertEqual(len(rows), 9)
            self.assertEqual(sum(row["mode"] == "REACTION_COORDINATE_FREE_BASELINE" for row in rows), 1)
            for mechanism in ("ADDITION_FIRST_RAW_NAC", "FULLY_CONCERTED_RAW_NAC"):
                self.assertEqual(
                    sorted(row["scale"] for row in rows if row["mode"] == mechanism),
                    [1, 2, 4, 8],
                )

    @unittest.skipUnless(DRIVER_PATH.exists(), "production driver not implemented yet")
    def test_sources_are_hash_fixed_raw_unbiased_nac_frames(self):
        expected = {
            26723: {
                "source_gro_sha256": "14477791ce14a35cef0adf9b802b562e091660526ca06de1132f6a74070faf10",
                "start_rst7_sha256": "2440de548c385f092c37f683de7b379ff5b18b6dc16593f7dbc80a9a8a167e14",
            },
            26737: {
                "source_gro_sha256": "a7924184ad3db4c13e0eab4929d46ee99621eacd523e899c0aca39c540350bc8",
                "start_rst7_sha256": "48c3944d3295158b06e96e32e4e07d9bcae9ceba2731f95aae9c1335ff972abe",
            },
        }
        for seed, hashes in expected.items():
            source = self.mod.SOURCES[seed]
            self.assertEqual(source["source_gro_sha256"], hashes["source_gro_sha256"])
            self.assertEqual(source["start_rst7_sha256"], hashes["start_rst7_sha256"])
            self.assertIn("pt2_preorganized_frame_extraction", str(source["source_gro"]))
            self.assertNotIn("attack_inherited", str(source["source_gro"]))

    @unittest.skipUnless(DRIVER_PATH.exists(), "production driver not implemented yet")
    def test_baseline_has_no_reaction_coordinate_restraints(self):
        baseline = self.mod.task_spec(0)
        self.assertEqual(baseline["mode"], "REACTION_COORDINATE_FREE_BASELINE")
        self.assertEqual(self.mod.reactive_restraints(baseline, {}), [])
        self.assertFalse(self.mod.nmropt_for_spec(baseline))

    @unittest.skipUnless(DRIVER_PATH.exists(), "production driver not implemented yet")
    def test_first_window_guard_is_reactant_safe_not_tetrahedral_candidate_gate(self):
        reactant_like = {
            "attack_A": 3.10,
            "cn_A": 1.35,
            "carbonyl_A": 1.23,
            "proton_nearest": "Nalpha",
        }
        overcompressed = dict(reactant_like, attack_A=1.10)
        self.assertTrue(self.mod.first_window_guard(reactant_like)["pass"])
        self.assertFalse(self.mod.first_window_guard(overcompressed)["pass"])
        self.assertNotIn("attack_candidate_upper", self.mod.first_window_guard(reactant_like)["checks"])

    @unittest.skipUnless(DRIVER_PATH.exists(), "production driver not implemented yet")
    def test_frozen_contract_and_no_automatic_downstream(self):
        self.assertEqual(self.mod.QM_CONTRACT["qm_atoms"], 146)
        self.assertEqual(self.mod.QM_CONTRACT["qm_charge"], 0)
        self.assertEqual(self.mod.QM_CONTRACT["electrons"], 510)
        self.assertEqual(self.mod.QM_CONTRACT["link_atoms"], 6)
        self.assertEqual(self.mod.QM_CONTRACT["qm_waters"], 0)
        self.assertEqual(
            self.mod.PRMTOP_SHA256,
            "a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0",
        )
        self.assertFalse(self.mod.AUTO_CONTINUATION)
        self.assertFalse(self.mod.AUTO_RELEASE)
        self.assertFalse(self.mod.AUTO_SHOOTING)
        self.assertFalse(self.mod.AUTO_PMF)


if __name__ == "__main__":
    unittest.main()
