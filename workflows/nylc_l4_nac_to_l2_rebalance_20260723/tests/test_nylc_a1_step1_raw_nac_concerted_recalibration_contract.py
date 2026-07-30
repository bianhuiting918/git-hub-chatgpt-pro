#!/usr/bin/env python3
import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "prepare_audit_nylc_a1_step1_raw_nac_concerted_recalibration.py"


def load_driver():
    spec = importlib.util.spec_from_file_location("concerted_recal", DRIVER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ConcertedRecalibrationContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_driver() if DRIVER.exists() else None

    def test_driver_exists(self):
        self.assertTrue(DRIVER.exists(), f"missing production driver: {DRIVER}")

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_four_tasks_cover_only_two_missing_concerted_conditions(self):
        self.assertEqual(self.mod.ARRAY_TASKS, 4)
        observed = [self.mod.task_spec(i) for i in range(4)]
        self.assertEqual(
            [(x["seed"], x["variant"], x["force_scale"], x["window_fraction"]) for x in observed],
            [
                (26723, "HALF_WINDOW_SCALE8", 8, 0.5),
                (26723, "FULL_WINDOW_SCALE16", 16, 1.0),
                (26737, "HALF_WINDOW_SCALE8", 8, 0.5),
                (26737, "FULL_WINDOW_SCALE16", 16, 1.0),
            ],
        )
        self.assertEqual([x["baseline_source_task_index"] for x in observed], [0, 0, 9, 9])

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_targets_are_concerted_and_variant_specific(self):
        source = {
            "attack_A": 3.00,
            "c12_n3_A": 1.36,
            "c12_o2_A": 1.22,
            "nalpha_hg1_A": 1.03,
            "hg1_n3_A": 2.00,
        }
        half = self.mod.stage_spec(self.mod.task_spec(0), source)
        full = self.mod.stage_spec(self.mod.task_spec(1), source)
        self.assertEqual(half["active_coordinates"], ["attack", "cn", "carbonyl", "pt"])
        self.assertAlmostEqual(half["targets"]["attack_A"], 2.98)
        self.assertAlmostEqual(half["targets"]["c12_n3_A"], 1.40)
        self.assertAlmostEqual(half["targets"]["c12_o2_A"], 1.235)
        self.assertAlmostEqual(half["targets"]["nalpha_hg1_A"], 1.055)
        self.assertAlmostEqual(half["targets"]["hg1_n3_A"], 1.975)
        self.assertAlmostEqual(full["targets"]["attack_A"], 2.96)
        self.assertAlmostEqual(full["targets"]["c12_n3_A"], 1.44)
        self.assertAlmostEqual(full["targets"]["c12_o2_A"], 1.25)
        self.assertAlmostEqual(full["targets"]["nalpha_hg1_A"], 1.08)
        self.assertAlmostEqual(full["targets"]["hg1_n3_A"], 1.95)
        self.assertEqual(half["force_scale"], 8)
        self.assertEqual(full["force_scale"], 16)

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_weakest_qualified_variant_is_selected_per_seed(self):
        rows = [
            {"task": self.mod.task_spec(0), "eligible_for_inherited_chain_v2": True},
            {"task": self.mod.task_spec(1), "eligible_for_inherited_chain_v2": True},
            {"task": self.mod.task_spec(2), "eligible_for_inherited_chain_v2": False},
            {"task": self.mod.task_spec(3), "eligible_for_inherited_chain_v2": True},
        ]
        selected = self.mod.select_weakest_by_seed(rows)
        self.assertEqual(
            [(x["seed"], x["selected_variant"], x["selected_force_scale"]) for x in selected],
            [(26723, "HALF_WINDOW_SCALE8", 8), (26737, "FULL_WINDOW_SCALE16", 16)],
        )

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_recalibration_has_no_automatic_downstream_action(self):
        description = self.mod.describe()
        self.assertEqual(description["source_calibration_job"], "62471131")
        self.assertEqual(description["automatic_downstream_action"], "NONE")
        self.assertFalse(self.mod.AUTO_CONTINUATION)
        self.assertFalse(self.mod.AUTO_RELEASE)
        self.assertFalse(self.mod.AUTO_SHOOTING)
        self.assertFalse(self.mod.AUTO_PMF)


if __name__ == "__main__":
    unittest.main()
