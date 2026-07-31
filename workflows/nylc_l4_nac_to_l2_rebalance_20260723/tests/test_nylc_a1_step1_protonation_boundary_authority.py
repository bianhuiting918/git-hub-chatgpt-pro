#!/usr/bin/env python3
import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "prepare_audit_nylc_a1_step1_protonation_boundary_authority.py"


def load_driver():
    spec = importlib.util.spec_from_file_location("a1_protonation_boundary", DRIVER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class ProtonationBoundaryAuthorityContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_driver() if DRIVER.exists() else None

    def test_driver_exists(self):
        self.assertTrue(DRIVER.exists(), f"missing production driver: {DRIVER}")

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_four_task_mapping_and_fixed_sources(self):
        self.assertEqual(self.mod.ARRAY_TASKS, 4)
        expected = [
            ("NEUTRAL_THR267_CURRENT_QM", 26723),
            ("NEUTRAL_THR267_CURRENT_QM", 26737),
            ("ACTIVATED_THR267_EXPANDED_THR268_QM", 26723),
            ("ACTIVATED_THR267_EXPANDED_THR268_QM", 26737),
        ]
        self.assertEqual(
            [(self.mod.task_spec(i)["variant"], self.mod.task_spec(i)["seed"])
             for i in range(4)],
            expected,
        )
        self.assertEqual(self.mod.SOURCES[26723]["frame"], 240)
        self.assertEqual(self.mod.SOURCES[26737]["frame"], 103)

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_qm_masks_and_boundary_move(self):
        self.assertEqual(tuple(self.mod.EXPANDED_THR268_ATOMS), tuple(range(8964, 8978)))
        self.assertNotIn("8964-8977", self.mod.CURRENT_QMMASK)
        self.assertIn("8964-8977", self.mod.EXPANDED_QMMASK)
        self.assertEqual(self.mod.QM_CONTRACTS["current"],
                         {"qm_atoms": 146, "qm_charge": 0, "electrons": 510, "link_atoms": 6})
        self.assertEqual(self.mod.QM_CONTRACTS["expanded_thr268"],
                         {"qm_atoms": 160, "qm_charge": 0, "electrons": 564, "link_atoms": 6})
        self.assertIn((8962, 8964), self.mod.CURRENT_BOUNDARIES)
        self.assertNotIn((8962, 8964), self.mod.EXPANDED_BOUNDARIES)
        self.assertIn((8976, 8978), self.mod.EXPANDED_BOUNDARIES)

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_authority_is_unrestrained_and_strictly_inherited(self):
        payload = self.mod.describe()
        self.assertEqual(self.mod.AUTHORITY_BLOCKS, 4)
        self.assertEqual(self.mod.BLOCK_MAXCYC, 250)
        self.assertTrue(payload["strict_serial_restart_inheritance"])
        self.assertFalse(payload["failed_restart_inheritance"])
        self.assertEqual(payload["automatic_downstream_action"], "NONE")
        text = self.mod.authority_input("@1,2")
        self.assertIn("nmropt=0", text)
        self.assertIn("ntr=0", text)
        self.assertNotIn("DISANG", text)
        self.assertNotIn("restraintmask", text)

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_neutral_hg1_placement_and_proton_classification(self):
        import numpy as np
        coordinates = np.zeros((3, 3), dtype=float)
        coordinates[0] = [0.0, 0.0, 0.0]   # OG1
        coordinates[1] = [2.0, 0.0, 0.0]   # Nalpha
        coordinates[2] = [1.8, 0.0, 0.0]   # HG1 activated source
        moved = self.mod.place_hg1_on_og1(coordinates, 1, 2, 3)
        self.assertAlmostEqual(float(np.linalg.norm(moved[2] - moved[0])),
                               self.mod.NEUTRAL_OH_A, places=6)
        self.assertTrue(np.allclose(moved[:2], coordinates[:2]))
        self.assertEqual(self.mod.classify_hg1_state(0.98, 2.1, 2.8), "OG1H_NEUTRAL")
        self.assertEqual(self.mod.classify_hg1_state(2.1, 1.02, 2.8), "NALPHA_H3_ACTIVATED")
        self.assertEqual(self.mod.classify_hg1_state(2.1, 2.2, 1.04), "N3H_TRANSFERRED")
        self.assertEqual(self.mod.classify_hg1_state(1.5, 1.6, 1.7), "OTHER_OR_SHARED")

    @unittest.skipUnless(DRIVER.exists(), "production driver not implemented yet")
    def test_heavy_skeleton_gate_excludes_proton_identity(self):
        required = set(self.mod.REQUIRED_HEAVY_BONDS)
        self.assertIn(frozenset((8949, 8953)), required)
        self.assertIn(frozenset((8953, 8955)), required)
        self.assertIn(frozenset((8955, 8960)), required)
        self.assertIn(frozenset((8962, 8964)), required)
        self.assertNotIn(frozenset((8949, 8961)), required)
        self.assertNotIn(frozenset((8960, 8961)), required)


if __name__ == "__main__":
    unittest.main()
