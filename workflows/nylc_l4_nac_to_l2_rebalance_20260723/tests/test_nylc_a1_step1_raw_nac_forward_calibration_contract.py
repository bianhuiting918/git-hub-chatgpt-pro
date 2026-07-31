#!/usr/bin/env python3
import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER_PATH = ROOT / "scripts" / "prepare_audit_nylc_a1_step1_raw_nac_forward_calibration.py"


def load_driver():
    spec = importlib.util.spec_from_file_location("raw_nac_forward", DRIVER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class RawNacForwardCalibrationContract(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.mod = load_driver() if DRIVER_PATH.exists() else None

    def test_driver_exists(self):
        self.assertTrue(DRIVER_PATH.exists(), f"missing production driver: {DRIVER_PATH}")

    @unittest.skipUnless(DRIVER_PATH.exists(), "production driver not implemented yet")
    def test_matrix_is_two_seeds_times_baseline_and_two_mechanisms_four_scales(self):
        self.assertEqual(self.mod.ARRAY_TASKS, 18)
        self.assertEqual(tuple(self.mod.FORCE_SCALES), (1, 2, 4, 8))
        specs = [self.mod.task_spec(i) for i in range(self.mod.ARRAY_TASKS)]
        for seed in (26723, 26737):
            rows = [row for row in specs if row["seed"] == seed]
            self.assertEqual(len(rows), 9)
            self.assertEqual(sum(row["mode"] == "REACTION_COORDINATE_FREE_BASELINE" for row in rows), 1)
            for mechanism in ("ADDITION_FIRST_RAW_NAC", "FULLY_CONCERTED_RAW_NAC"):
                self.assertEqual(
                    sorted(row["scale"] for row in rows if row["mode"] == mechanism),
                    [1, 2, 4, 8],
                )

    @unittest.skipUnless(DRIVER_PATH.exists(), "production driver not implemented yet")
    def test_sources_are_hash_fixed_raw_unbiased_nac_frames(self):
        expected = {
            26723: {
                "source_gro_sha256": "14477791ce14a35cef0adf9b802b562e091660526ca06de1132f6a74070faf10",
                "start_rst7_sha256": "2440de548c385f092c37f683de7b379ff5b18b6dc16593f7dbc80a9a8a167e14",
            },
            26737: {
                "source_gro_sha256": "a7924184ad3db4c13e0eab4929d46ee99621eacd523e899c0aca39c540350bc8",
                "start_rst7_sha256": "48c3944d3295158b06e96e32e4e07d9bcae9ceba2731f95aae9c1335ff972abe",
            },
        }
        for seed, hashes in expected.items():
            source = self.mod.SOURCES[seed]
            self.assertEqual(source["source_gro_sha256"], hashes["source_gro_sha256"])
            self.assertEqual(source["start_rst7_sha256"], hashes["start_rst7_sha256"])
            self.assertIn("pt2_preorganized_frame_extraction", str(source["source_gro"]))
            self.assertNotIn("attack_inherited", str(source["source_gro"]))

    @unittest.skipUnless(DRIVER_PATH.exists(), "production driver not implemented yet")
    def test_baseline_has_no_reaction_coordinate_restraints(self):
        baseline = self.mod.task_spec(0)
        self.assertEqual(baseline["mode"], "REACTION_COORDINATE_FREE_BASELINE")
        self.assertEqual(self.mod.reactive_restraints(baseline, {}), [])
        self.assertFalse(self.mod.nmropt_for_spec(baseline))

    @unittest.skipUnless(DRIVER_PATH.exists(), "production driver not implemented yet")
    def test_first_window_guard_is_reactant_safe_not_tetrahedral_candidate_gate(self):
        reactant_like = {
            "attack_A": 3.10,
            "cn_A": 1.35,
            "carbonyl_A": 1.23,
            "proton_nearest": "Nalpha",
        }
        overcompressed = dict(reactant_like, attack_A=1.10)
        self.assertTrue(self.mod.first_window_guard(reactant_like)["pass"])
        self.assertFalse(self.mod.first_window_guard(overcompressed)["pass"])
        self.assertNotIn("attack_candidate_upper", self.mod.first_window_guard(reactant_like)["checks"])

    @unittest.skipUnless(DRIVER_PATH.exists(), "production driver not implemented yet")
    def test_full_geometry_contract_requires_qm_heavy_atom_indices(self):
        self.assertTrue(
            hasattr(self.mod, "validate_full_contract"),
            "missing full QM geometry-contract validator",
        )
        payload = {
            "qm_contract": {
                "qm_atom_count": 146,
                "qmcharge": 0,
                "electron_count_including_link_h": 510,
                "link_atom_count": 6,
                "step1_qm_water_count": 0,
                "qmmask": "@1",
                "qm_heavy_atom_indices": [1],
            }
        }
        contract = self.mod.validate_full_contract(payload, "@1")
        self.assertEqual(contract["qm_heavy_atom_indices"], [1])
        del payload["qm_contract"]["qm_heavy_atom_indices"]
        with self.assertRaises(ValueError):
            self.mod.validate_full_contract(payload, "@1")

    @unittest.skipUnless(DRIVER_PATH.exists(), "production driver not implemented yet")
    def test_frozen_contract_and_no_automatic_downstream(self):
        self.assertEqual(self.mod.QM_CONTRACT["qm_atoms"], 146)
        self.assertEqual(self.mod.QM_CONTRACT["qm_charge"], 0)
        self.assertEqual(self.mod.QM_CONTRACT["electrons"], 510)
        self.assertEqual(self.mod.QM_CONTRACT["link_atoms"], 6)
        self.assertEqual(self.mod.QM_CONTRACT["qm_waters"], 0)
        self.assertEqual(
            self.mod.PRMTOP_SHA256,
            "a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0",
        )
        self.assertFalse(self.mod.AUTO_CONTINUATION)
        self.assertFalse(self.mod.AUTO_RELEASE)
        self.assertFalse(self.mod.AUTO_SHOOTING)
        self.assertFalse(self.mod.AUTO_PMF)

    @unittest.skipUnless(DRIVER_PATH.exists(), "production driver not implemented yet")
    def test_v2_relative_gate_handles_baseline_overshoot_without_sign_error(self):
        item = {
            "source_geometry": {
                "attack_A": 3.00,
                "c12_n3_A": 1.35,
                "c12_o2_A": 1.22,
                "nalpha_hg1_A": 1.02,
                "hg1_n3_A": 2.00,
            },
            "target_geometry": {
                "attack_A": 2.96,
                "c12_n3_A": 1.35,
                "c12_o2_A": 1.25,
                "nalpha_hg1_A": 1.02,
                "hg1_n3_A": 2.00,
            },
            "final_geometry": {
                "attack_A": 2.97,
                "c12_n3_A": 1.35,
                "c12_o2_A": 1.251,
                "nalpha_hg1_A": 1.02,
                "hg1_n3_A": 2.00,
            },
            "response_from_raw_source": {
                "active_coordinates": ["attack", "carbonyl"]
            },
        }
        baseline = {
            "final_geometry": {
                "attack_A": 3.20,
                "c12_n3_A": 1.35,
                "c12_o2_A": 1.27,
                "nalpha_hg1_A": 1.02,
                "hg1_n3_A": 2.00,
            }
        }
        legacy = self.mod.force_effect_vs_baseline(item, baseline)
        revised = self.mod.force_effect_vs_baseline_v2(item, baseline)
        self.assertFalse(legacy["all"])
        self.assertTrue(revised["all"])
        self.assertFalse(
            revised["checks"]["carbonyl"]["forced_minus_baseline_toward_raw_target"]
        )
        self.assertTrue(
            revised["checks"]["carbonyl"]["forced_closer_to_target_than_baseline"]
        )


    @unittest.skipUnless(DRIVER_PATH.exists(), "production driver not implemented yet")
    def test_engine_overflow_is_technical_failure_and_cannot_be_selected(self):
        catastrophic = """
   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
     2200      -1.0000E+08     3.2922E+10     1.2916E+13   C12     10287
 DFTBESCF=**************
 EAMBER  = *************
"""
        healthy = """
   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER
     2200      -5.7000E+05     2.3340E-01     8.6258E+01   C12     10287
 DFTBESCF=   -15545.0964
 EAMBER  =  -569216.7280
"""
        bad = self.mod.engine_numerical_health(catastrophic)
        good = self.mod.engine_numerical_health(healthy)
        self.assertFalse(bad["pass"])
        self.assertIn("ENERGY_FIELD_OVERFLOW", bad["hard_errors"])
        self.assertIn("RMS_NUMERICAL_DIVERGENCE", bad["hard_errors"])
        self.assertIn("GMAX_NUMERICAL_DIVERGENCE", bad["hard_errors"])
        self.assertTrue(good["pass"])
        self.assertEqual(good["hard_errors"], [])


if __name__ == "__main__":
    unittest.main()
