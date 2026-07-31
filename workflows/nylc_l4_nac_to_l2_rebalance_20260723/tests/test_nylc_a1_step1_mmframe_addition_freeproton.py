#!/usr/bin/env python3
import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "prepare_audit_nylc_a1_step1_mmframe_addition_freeproton.py"


def load_driver():
    spec = importlib.util.spec_from_file_location("mmframe_addition_freeproton", DRIVER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class MmFrameAdditionFreeProtonContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_driver() if DRIVER.exists() else None

    def test_driver_exists(self):
        self.assertTrue(DRIVER.exists(), f"missing production driver: {DRIVER}")

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_two_approved_mm_sources(self):
        self.assertEqual(self.mod.ARRAY_TASKS, 2)
        self.assertEqual(self.mod.SOURCES[26723]["frame"], 240)
        self.assertEqual(self.mod.SOURCES[26723]["time_ps"], 480.0)
        self.assertEqual(self.mod.SOURCES[26737]["frame"], 103)
        self.assertEqual(self.mod.SOURCES[26737]["time_ps"], 206.0)
        self.assertIn("run.tpr", self.mod.SOURCES[26723]["tpr"])
        self.assertIn("run.xtc", self.mod.SOURCES[26737]["xtc"])

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_six_window_schedule_and_base_forces(self):
        self.assertEqual(self.mod.WINDOWS, 6)
        self.assertEqual(self.mod.ATTACK_DELTA_A, -0.04)
        self.assertEqual(self.mod.CARBONYL_DELTA_A, 0.03)
        self.assertEqual(self.mod.FORCE_ATTACK, 24.0)
        self.assertEqual(self.mod.FORCE_CARBONYL, 30.0)
        source = {"attack_A": 2.972, "c12_o2_A": 1.220}
        last = self.mod.stage_spec(6, source)
        self.assertAlmostEqual(last["attack_target_A"], 2.732)
        self.assertAlmostEqual(last["carbonyl_target_A"], 1.400)

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_restraints_bias_only_attack_and_carbonyl(self):
        stage = self.mod.stage_spec(
            1, {"attack_A": 2.972, "c12_o2_A": 1.220}
        )
        text = self.mod.restraints(stage)
        self.assertNotIn(r"\n", text)
        self.assertEqual(text.count("&rst"), 2)
        self.assertEqual(text.count("/\n"), 2)
        self.assertNotIn(str(self.mod.REACTIVE["hg1"]), text)
        self.assertNotIn(str(self.mod.REACTIVE["n3"]), text)

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_describe_forbids_downstream_and_proton_bias(self):
        payload = self.mod.describe()
        self.assertEqual(payload["array_task_count"], 2)
        self.assertEqual(payload["windows_per_task"], 6)
        self.assertEqual(payload["reaction_coordinate_restraints"], 2)
        self.assertFalse(payload["proton_coordinate_restrained"])
        self.assertFalse(payload["attack_angle_restrained"])
        self.assertEqual(payload["automatic_downstream_action"], "NONE")


if __name__ == "__main__":
    unittest.main()
