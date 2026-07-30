#!/usr/bin/env python3
import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "prepare_audit_nylc_a1_step1_raw_nac_inherited_chain.py"


def load_driver():
    spec = importlib.util.spec_from_file_location("raw_nac_chain", DRIVER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class RawNacInheritedChainContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_driver() if DRIVER.exists() else None

    def test_driver_exists(self):
        self.assertTrue(DRIVER.exists(), f"missing production driver: {DRIVER}")

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_four_branches_are_selected_from_completed_calibration(self):
        self.assertEqual(self.mod.ARRAY_TASKS, 4)
        self.assertEqual(self.mod.CALIBRATION_JOB, "62471131")
        payload = {
            "status": "PASS_RAW_NAC_FORWARD_FORCE_CALIBRATION_MATRIX",
            "selected_minimum_force_scales": [
                {"seed": 26723, "mode": "ADDITION_FIRST_RAW_NAC", "selected_weakest_scale": 2},
                {"seed": 26723, "mode": "FULLY_CONCERTED_RAW_NAC", "selected_weakest_scale": 4},
                {"seed": 26737, "mode": "ADDITION_FIRST_RAW_NAC", "selected_weakest_scale": 1},
                {"seed": 26737, "mode": "FULLY_CONCERTED_RAW_NAC", "selected_weakest_scale": 8},
            ],
            "per_task": [],
        }
        observed = [self.mod.selection_from_audit(payload, i) for i in range(4)]
        self.assertEqual(
            [(row["seed"], row["mode"], row["scale"]) for row in observed],
            [
                (26723, "ADDITION_FIRST_RAW_NAC", 2),
                (26723, "FULLY_CONCERTED_RAW_NAC", 4),
                (26737, "ADDITION_FIRST_RAW_NAC", 1),
                (26737, "FULLY_CONCERTED_RAW_NAC", 8),
            ],
        )

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_missing_selection_is_not_silently_substituted(self):
        payload = {
            "selected_minimum_force_scales": [
                {"seed": 26723, "mode": "ADDITION_FIRST_RAW_NAC", "selected_weakest_scale": None}
            ],
            "per_task": [],
        }
        with self.assertRaises(ValueError):
            self.mod.selection_from_audit(payload, 0)

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_targets_are_relative_and_mechanism_specific(self):
        geometry = {
            "attack_A": 2.90,
            "c12_n3_A": 1.36,
            "c12_o2_A": 1.24,
            "nalpha_hg1_A": 1.03,
            "hg1_n3_A": 2.10,
        }
        addition = self.mod.propose_targets("ADDITION_FIRST_RAW_NAC", geometry, 1.0)
        concerted = self.mod.propose_targets("FULLY_CONCERTED_RAW_NAC", geometry, 0.5)
        self.assertAlmostEqual(addition["attack_A"], 2.86)
        self.assertAlmostEqual(addition["c12_o2_A"], 1.27)
        self.assertAlmostEqual(addition["c12_n3_A"], geometry["c12_n3_A"])
        self.assertAlmostEqual(addition["nalpha_hg1_A"], geometry["nalpha_hg1_A"])
        self.assertAlmostEqual(concerted["attack_A"], 2.88)
        self.assertAlmostEqual(concerted["c12_o2_A"], 1.255)
        self.assertAlmostEqual(concerted["c12_n3_A"], 1.40)
        self.assertAlmostEqual(concerted["nalpha_hg1_A"], 1.055)
        self.assertAlmostEqual(concerted["hg1_n3_A"], 2.075)

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_strict_inheritance_and_no_automatic_downstream(self):
        description = self.mod.describe()
        self.assertTrue(description["strict_serial_restart_inheritance"])
        self.assertFalse(description["failed_restart_inheritance"])
        self.assertEqual(description["automatic_downstream_action"], "NONE")
        self.assertEqual(self.mod.MAX_WINDOWS, 16)
        self.assertEqual(self.mod.MIN_STEP, 0.125)
        self.assertFalse(self.mod.AUTO_RELEASE)
        self.assertFalse(self.mod.AUTO_SHOOTING)
        self.assertFalse(self.mod.AUTO_PMF)

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_minimization_is_tight_and_restraints_have_physical_lines(self):
        spec = {"task_index": 0, "mode": "ADDITION_FIRST_RAW_NAC", "scale": 1}
        stage = {
            "window_index": 1,
            "force_scale": 1,
            "targets": {
                "attack_A": 2.86,
                "c12_n3_A": 1.36,
                "c12_o2_A": 1.27,
                "nalpha_hg1_A": 1.03,
                "hg1_n3_A": 2.10,
            },
            "active_coordinates": ["attack", "carbonyl"],
            "force_bond": 24.0,
            "force_carbonyl": 30.0,
            "force_pt": 18.0,
            "maxcyc": 2200,
            "ncyc": 550,
        }
        text = self.mod.minimization_input(spec, stage, "@1")
        self.assertIn("drms=0.01", text)
        rst = self.mod.restraints(spec, stage)
        self.assertNotIn(r"\n", rst)
        self.assertEqual(rst.count("&rst"), 2)
        self.assertEqual(rst.count("/\n"), 2)


if __name__ == "__main__":
    unittest.main()
