#!/usr/bin/env python3
import importlib.util
import pathlib
import unittest

HERE = pathlib.Path(__file__).resolve().parent
DRIVER = HERE.parent / "scripts" / "prepare_audit_nylc_a1_step1_reverse_force_calibration.py"


class ReverseForceCalibrationContract(unittest.TestCase):
    def test_frozen_reverse_force_matrix(self):
        self.assertTrue(DRIVER.is_file(), "reverse force calibration driver is missing")
        spec = importlib.util.spec_from_file_location("_reverse_force_cal", DRIVER)
        self.assertIsNotNone(spec)
        self.assertIsNotNone(spec.loader)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)

        description = module.describe()
        self.assertEqual(description["array_task_count"], 12)
        self.assertEqual(description["mechanisms"], ["STEPWISE_REVERSE", "CONCERTED_REVERSE"])
        self.assertEqual(description["force_scales"], [1.0, 2.0, 4.0])
        self.assertTrue(description["first_window_only"])
        self.assertTrue(description["strict_serial_restart_inheritance"])
        self.assertFalse(description["failed_restart_inheritance"])

        expected = [
            ("seed26723", "STEPWISE_REVERSE", 1.0),
            ("seed26723", "STEPWISE_REVERSE", 2.0),
            ("seed26723", "STEPWISE_REVERSE", 4.0),
            ("seed26723", "CONCERTED_REVERSE", 1.0),
            ("seed26723", "CONCERTED_REVERSE", 2.0),
            ("seed26723", "CONCERTED_REVERSE", 4.0),
            ("seed26737", "STEPWISE_REVERSE", 1.0),
            ("seed26737", "STEPWISE_REVERSE", 2.0),
            ("seed26737", "STEPWISE_REVERSE", 4.0),
            ("seed26737", "CONCERTED_REVERSE", 1.0),
            ("seed26737", "CONCERTED_REVERSE", 2.0),
            ("seed26737", "CONCERTED_REVERSE", 4.0),
        ]
        observed = [
            (task["seed"], task["mechanism"], task["force_scale"])
            for task in (module.task_spec(index) for index in range(12))
        ]
        self.assertEqual(observed, expected)

        source = {
            "attack_A": 1.41,
            "c12_n3_A": 3.15,
            "c12_o2_A": 1.24,
            "nalpha_hg1_A": 3.00,
            "hg1_n3_A": 1.02,
        }
        stepwise = module.stage_spec(module.task_spec(0), source)
        concerted = module.stage_spec(module.task_spec(3), source)
        self.assertEqual(stepwise["window_index"], 0)
        self.assertEqual(stepwise["targets"]["attack_A"], source["attack_A"])
        self.assertLess(stepwise["targets"]["c12_n3_A"], source["c12_n3_A"])
        self.assertGreater(stepwise["targets"]["c12_o2_A"], source["c12_o2_A"])
        self.assertLess(
            stepwise["targets"]["nalpha_hg1_A"] - stepwise["targets"]["hg1_n3_A"],
            source["nalpha_hg1_A"] - source["hg1_n3_A"],
        )
        self.assertGreater(concerted["targets"]["attack_A"], source["attack_A"])


if __name__ == "__main__":
    unittest.main()
