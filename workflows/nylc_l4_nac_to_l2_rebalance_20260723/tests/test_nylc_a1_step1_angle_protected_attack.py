#!/usr/bin/env python3
import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "prepare_audit_nylc_a1_step1_angle_protected_attack.py"


def load_driver():
    spec = importlib.util.spec_from_file_location("angle_protected_attack", DRIVER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class AngleProtectedAttackContract(unittest.TestCase):
    def test_driver_exists(self):
        self.assertTrue(DRIVER.exists(), f"missing production driver: {DRIVER}")

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_two_original_neutral_sources_only(self):
        mod = load_driver()
        self.assertEqual(mod.ARRAY_TASKS, 2)
        specs = [mod.task_spec(i) for i in range(2)]
        self.assertEqual([x["seed"] for x in specs], ["seed26723", "seed26737"])
        self.assertTrue(all(x["source_kind"] == "ORIGINAL_NEUTRAL_SOURCE" for x in specs))
        self.assertEqual(
            mod.SOURCES["seed26723"]["restart_sha256"],
            "486e5c3d479e42d78c2d87a670320bbce9d1fcedaa0d5ad3bd99bf4c696506b2",
        )
        self.assertEqual(
            mod.SOURCES["seed26737"]["restart_sha256"],
            "bd1717782acb8b6ec2d2c92a1b4f84d4354baf5797df4577fa80aa57541edccc",
        )

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_attack_only_progress_with_correct_angle_guard(self):
        mod = load_driver()
        source = {
            "attack_A": 3.00,
            "c12_n3_A": 1.40,
            "c12_o2_A": 1.22,
            "nalpha_hg1_A": 0.98,
            "hg1_n3_A": 2.00,
            "attack_angle_deg": 106.0,
        }
        stage = mod.stage_spec(mod.task_spec(0), source)
        self.assertEqual(stage["active_coordinates"], ["attack"])
        self.assertAlmostEqual(stage["targets"]["attack_A"], 2.96)
        self.assertEqual(stage["force_bond"], 24.0)
        self.assertEqual(stage["force_angle_guard"], 2.0)
        self.assertEqual(stage["angle_flat_bottom_deg"], [100.0, 115.0])
        text = mod.restraints(stage)
        rows = [row for row in text.splitlines() if row.strip()]
        self.assertEqual(len(rows), 2)
        self.assertIn(
            f"iat={mod.FWD.REACTIVE['o2']},{mod.FWD.REACTIVE['c12']},{mod.FWD.REACTIVE['og1']}",
            rows[1],
        )
        self.assertIn("r2=100.000", rows[1])
        self.assertIn("r3=115.000", rows[1])
        self.assertNotIn(f"iat={mod.FWD.REACTIVE['c12']},{mod.FWD.REACTIVE['o2']}", text)
        self.assertNotIn(f"iat={mod.FWD.REACTIVE['c12']},{mod.FWD.REACTIVE['n3']}", text)
        self.assertNotIn(f"iat={mod.FWD.REACTIVE['nalpha']},{mod.FWD.REACTIVE['hg1']}", text)

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_one_window_and_strict_stop_contract(self):
        mod = load_driver()
        desc = mod.describe()
        self.assertTrue(desc["first_window_only"])
        self.assertEqual(desc["automatic_downstream_action"], "NONE")
        self.assertEqual(desc["angle_definition"], "O2-C12-OG1")
        self.assertEqual(desc["angle_flat_bottom_deg"], [100.0, 115.0])
        self.assertEqual(desc["stop_conditions"], [
            "technical_or_numerical_failure",
            "thr267_chemical_integrity_failure",
            "attack_response_wrong_direction",
            "attack_angle_below_100_deg",
        ])


if __name__ == "__main__":
    unittest.main()
