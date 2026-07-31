#!/usr/bin/env python3
import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "prepare_audit_nylc_a1_step1_activated_mmframe_addition.py"


def load_driver():
    spec = importlib.util.spec_from_file_location("activated_mmframe_addition", DRIVER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ActivatedMmFrameAdditionContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_driver() if DRIVER.exists() else None

    def test_driver_exists(self):
        self.assertTrue(DRIVER.exists(), f"missing production driver: {DRIVER}")

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_two_sha_frozen_activated_sources(self):
        self.assertEqual(self.mod.ARRAY_TASKS, 2)
        self.assertEqual(self.mod.SOURCES[26723]["frame"], 240)
        self.assertEqual(self.mod.SOURCES[26723]["time_ps"], 480.0)
        self.assertEqual(
            self.mod.SOURCES[26723]["tpr_sha256"],
            "c60078a92c2ace51facde4ef64e453f690177fc88b4d6363427f935944fa2e43",
        )
        self.assertEqual(
            self.mod.SOURCES[26723]["xtc_sha256"],
            "1a54f1b5b9f139b746c22d9e0f7e9a4a94eb8154bf2b881888986eedca933d89",
        )
        self.assertEqual(self.mod.SOURCES[26737]["frame"], 103)
        self.assertEqual(self.mod.SOURCES[26737]["time_ps"], 206.0)
        self.assertEqual(
            self.mod.SOURCES[26737]["tpr_sha256"],
            "dbd19a399547319d10630430ed494d33f6271bab0f30cb0de5a466c6af13ba20",
        )
        self.assertEqual(
            self.mod.SOURCES[26737]["xtc_sha256"],
            "fcba14da98b331368061dcd990f2467628ad77b9b7a9c4ce88090f92e0831b05",
        )

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_activated_a1_topology_contract(self):
        bonds = self.mod.REQUIRED_A1_BONDS
        self.assertIn(
            frozenset((self.mod.THR267["nalpha"], self.mod.THR267["hg1"])),
            bonds,
        )
        self.assertIn(
            frozenset((self.mod.THR267["cb"], self.mod.THR267["og1"])),
            bonds,
        )
        self.assertEqual(
            self.mod.FORBIDDEN_A1_BOND,
            frozenset((self.mod.THR267["og1"], self.mod.THR267["hg1"])),
        )
        self.assertEqual(
            self.mod.describe()["starting_state"], "NALPHA_H3_PLUS_OG1_MINUS"
        )

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_authority_stage_is_reaction_coordinate_free(self):
        payload = self.mod.describe()
        self.assertEqual(self.mod.authority_restraints(), "")
        self.assertFalse(payload["authority_reaction_coordinate_restrained"])
        self.assertFalse(payload["authority_position_restrained"])
        self.assertEqual(payload["authority_stage_count"], 1)

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
    def test_phase1_biases_only_attack_and_carbonyl(self):
        stage = self.mod.stage_spec(
            1, {"attack_A": 2.972, "c12_o2_A": 1.220}
        )
        text = self.mod.reaction_restraints(stage)
        self.assertNotIn(r"\n", text)
        self.assertEqual(text.count("&rst"), 2)
        self.assertEqual(text.count("/\n"), 2)
        self.assertNotIn(str(self.mod.THR267["hg1"]), text)
        self.assertNotIn(str(self.mod.REACTIVE["n3"]), text)

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_describe_forbids_automatic_downstream(self):
        payload = self.mod.describe()
        self.assertEqual(payload["array_task_count"], 2)
        self.assertEqual(payload["windows_per_task"], 6)
        self.assertEqual(payload["reaction_coordinate_restraints"], 2)
        self.assertFalse(payload["proton_coordinate_restrained"])
        self.assertFalse(payload["attack_angle_restrained"])
        self.assertEqual(payload["automatic_downstream_action"], "NONE")


if __name__ == "__main__":
    unittest.main()
