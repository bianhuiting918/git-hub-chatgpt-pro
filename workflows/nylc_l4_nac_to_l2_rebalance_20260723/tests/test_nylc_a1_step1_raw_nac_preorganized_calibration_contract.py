#!/usr/bin/env python3
import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "prepare_audit_nylc_a1_step1_raw_nac_preorganized_calibration.py"


def load_driver():
    spec = importlib.util.spec_from_file_location("raw_nac_preorganized", DRIVER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class RawNacPreorganizedCalibrationContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_driver() if DRIVER.exists() else None

    def test_driver_exists(self):
        self.assertTrue(DRIVER.exists(), f"missing production driver: {DRIVER}")

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_matrix_is_two_seeds_times_independent_attack_and_carbonyl_scales(self):
        self.assertEqual(self.mod.ARRAY_TASKS, 18)
        self.assertEqual(tuple(self.mod.ATTACK_FORCE_SCALES), (1, 2, 4))
        self.assertEqual(tuple(self.mod.CARBONYL_FORCE_SCALES), (1, 2, 4))
        specs = [self.mod.task_spec(i) for i in range(self.mod.ARRAY_TASKS)]
        for seed in (26723, 26737):
            rows = [row for row in specs if row["seed"] == seed]
            self.assertEqual(len(rows), 9)
            self.assertEqual(
                {(row["attack_scale"], row["carbonyl_scale"]) for row in rows},
                {(a, c) for a in (1, 2, 4) for c in (1, 2, 4)},
            )

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_sources_are_original_hash_fixed_nac_restarts(self):
        self.assertEqual(
            self.mod.SOURCES[26723]["start_rst7_sha256"],
            "2440de548c385f092c37f683de7b379ff5b18b6dc16593f7dbc80a9a8a167e14",
        )
        self.assertEqual(
            self.mod.SOURCES[26737]["start_rst7_sha256"],
            "48c3944d3295158b06e96e32e4e07d9bcae9ceba2731f95aae9c1335ff972abe",
        )
        self.assertNotIn("62471131", str(self.mod.SOURCES))
        self.assertNotIn("62477118", str(self.mod.SOURCES))

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_attack_and_carbonyl_forces_are_independent_with_fixed_weak_guards(self):
        spec = self.mod.task_spec(5)
        source = {
            "attack_A": 3.00,
            "c12_n3_A": 1.35,
            "c12_o2_A": 1.22,
            "nalpha_hg1_A": 1.02,
            "hg1_n3_A": 2.00,
            "attack_angle_deg": 105.0,
        }
        stage = self.mod.stage_spec(spec, source)
        self.assertNotEqual(stage["force_attack"], stage["force_carbonyl"])
        self.assertEqual(stage["force_angle_guard"], self.mod.WEAK_ANGLE_GUARD_FORCE)
        self.assertEqual(stage["force_pt_guard"], self.mod.WEAK_PT_GUARD_FORCE)
        rst = self.mod.restraints(stage)
        self.assertNotIn(r"\n", rst)
        self.assertEqual(rst.count("&rst"), 5)
        self.assertEqual(rst.count("/\n"), 5)

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_gate_requires_numerical_health_forward_response_and_preorganization(self):
        healthy = {
            "attack_A": 2.96,
            "c12_n3_A": 1.35,
            "c12_o2_A": 1.25,
            "nalpha_hg1_A": 1.02,
            "hg1_n3_A": 2.00,
            "attack_angle_deg": 103.0,
            "hg1_nearest_qm_heavy_atom": "Nalpha",
        }
        bad_pt = dict(healthy, nalpha_hg1_A=1.50, hg1_n3_A=1.20)
        self.assertTrue(self.mod.preorganization_guard(healthy)["pass"])
        self.assertFalse(self.mod.preorganization_guard(bad_pt)["pass"])
        self.assertFalse(self.mod.AUTO_CONTINUATION)
        self.assertFalse(self.mod.AUTO_RELEASE)
        self.assertFalse(self.mod.AUTO_SHOOTING)
        self.assertFalse(self.mod.AUTO_PMF)


if __name__ == "__main__":
    unittest.main()
