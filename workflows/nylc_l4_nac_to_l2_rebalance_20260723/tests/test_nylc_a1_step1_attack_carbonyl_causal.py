#!/usr/bin/env python3
import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "prepare_audit_nylc_a1_step1_attack_carbonyl_causal.py"


def load_driver():
    spec = importlib.util.spec_from_file_location("attack_carbonyl_causal", DRIVER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class AttackCarbonylCausalContract(unittest.TestCase):
    def test_driver_exists(self):
        self.assertTrue(DRIVER.exists(), f"missing production driver: {DRIVER}")

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_matrix_is_two_seeds_times_three_routes(self):
        mod = load_driver()
        self.assertEqual(mod.ARRAY_TASKS, 6)
        specs = [mod.task_spec(i) for i in range(6)]
        self.assertEqual(
            [(x["seed"], x["route"], x["source_kind"]) for x in specs],
            [
                ("seed26723", "ENDPOINT_ATTACK_CARBONYL", "OFF_ANGLE_ENDPOINT"),
                ("seed26723", "ENDPOINT_ATTACK_ONLY", "OFF_ANGLE_ENDPOINT"),
                ("seed26723", "SOURCE_ATTACK_ONLY", "ORIGINAL_NEUTRAL_SOURCE"),
                ("seed26737", "ENDPOINT_ATTACK_CARBONYL", "OFF_ANGLE_ENDPOINT"),
                ("seed26737", "ENDPOINT_ATTACK_ONLY", "OFF_ANGLE_ENDPOINT"),
                ("seed26737", "SOURCE_ATTACK_ONLY", "ORIGINAL_NEUTRAL_SOURCE"),
            ],
        )

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_only_requested_coordinates_are_restrained(self):
        mod = load_driver()
        source = {
            "attack_A": 3.00,
            "c12_n3_A": 1.40,
            "c12_o2_A": 1.22,
            "nalpha_hg1_A": 1.90,
            "hg1_n3_A": 2.10,
        }
        stages = [mod.stage_spec(mod.task_spec(i), source) for i in range(3)]
        self.assertEqual(stages[0]["active_coordinates"], ["attack", "carbonyl"])
        self.assertEqual(stages[1]["active_coordinates"], ["attack"])
        self.assertEqual(stages[2]["active_coordinates"], ["attack"])
        self.assertAlmostEqual(stages[0]["targets"]["attack_A"], 2.96)
        self.assertAlmostEqual(stages[0]["targets"]["c12_o2_A"], 1.25)
        self.assertAlmostEqual(stages[1]["targets"]["c12_o2_A"], 1.22)
        self.assertAlmostEqual(stages[2]["targets"]["c12_o2_A"], 1.22)
        self.assertEqual(stages[0]["force_bond"], 24.0)
        self.assertEqual(stages[0]["force_carbonyl"], 30.0)
        self.assertNotIn("angle", stages[0]["active_coordinates"])
        self.assertNotIn("cn", stages[0]["active_coordinates"])
        self.assertNotIn("pt", stages[0]["active_coordinates"])

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_sources_are_fixed_and_one_window_only(self):
        mod = load_driver()
        self.assertEqual(
            mod.SOURCES["seed26723"]["endpoint_sha256"],
            "58110d8e77c3b6beb15154d74140a847b31aa545187072f9737288e7ae0635f8",
        )
        self.assertEqual(
            mod.SOURCES["seed26737"]["endpoint_sha256"],
            "18472aa8a2770a55431bed09035ea633e10e2cd0950d33979eb812cf3c1cd0ad",
        )
        desc = mod.describe()
        self.assertTrue(desc["first_window_only"])
        self.assertEqual(desc["automatic_downstream_action"], "NONE")
        self.assertEqual(desc["force_attack_kcal_mol_A2"], 24.0)
        self.assertEqual(desc["force_carbonyl_kcal_mol_A2"], 30.0)


if __name__ == "__main__":
    unittest.main()
