#!/usr/bin/env python3
import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "prepare_audit_nylc_a1_step1_neutral_forward_calibration.py"


def load_driver():
    spec = importlib.util.spec_from_file_location("neutral_forward", DRIVER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class NeutralForwardCalibrationContract(unittest.TestCase):
    def test_driver_exists(self):
        self.assertTrue(DRIVER.exists(), f"missing production driver: {DRIVER}")

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_matrix_is_two_neutral_seeds_times_two_mechanisms(self):
        mod = load_driver()
        self.assertEqual(mod.ARRAY_TASKS, 4)
        specs = [mod.task_spec(i) for i in range(4)]
        self.assertEqual(
            [(x["seed"], x["mechanism"], x["force_scale"]) for x in specs],
            [
                ("seed26723", "ADDITION_FIRST_FORWARD", 1.0),
                ("seed26723", "FULLY_CONCERTED_FORWARD", 1.0),
                ("seed26737", "ADDITION_FIRST_FORWARD", 1.0),
                ("seed26737", "FULLY_CONCERTED_FORWARD", 1.0),
            ],
        )
        self.assertTrue(all(x["starting_state"] == "NALPHA_H2_OG1H" for x in specs))
        self.assertTrue(all(x["qm_contract_key"] == "current" for x in specs))

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_sources_are_fixed_neutral_authority_inputs(self):
        mod = load_driver()
        self.assertEqual(
            mod.SOURCES["seed26723"]["restart_sha256"],
            "486e5c3d479e42d78c2d87a670320bbce9d1fcedaa0d5ad3bd99bf4c696506b2",
        )
        self.assertEqual(
            mod.SOURCES["seed26737"]["restart_sha256"],
            "bd1717782acb8b6ec2d2c92a1b4f84d4354baf5797df4577fa80aa57541edccc",
        )
        self.assertTrue(all("attempt_62533819_" in str(x["restart"]) for x in mod.SOURCES.values()))

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_first_window_only_and_no_automatic_downstream(self):
        mod = load_driver()
        desc = mod.describe()
        self.assertTrue(desc["first_window_only"])
        self.assertEqual(desc["automatic_downstream_action"], "NONE")
        self.assertEqual(desc["qm_contract"], {"qm_atoms": 146, "qm_charge": 0, "electrons": 510, "link_atoms": 6})


if __name__ == "__main__":
    unittest.main()
