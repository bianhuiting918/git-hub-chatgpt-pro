#!/usr/bin/env python3
import pathlib
import unittest

HERE = pathlib.Path(__file__).resolve().parents[1]
AUDIT = HERE / "scripts" / "audit_nylc_a1_step1_pt2_coordinate.py"


class A1Step1PT2CoordinateContract(unittest.TestCase):
    def test_auditor_exists(self):
        self.assertTrue(AUDIT.is_file(), AUDIT)

    def test_frozen_source_and_atom_identity_contract(self):
        text = AUDIT.read_text(encoding="utf-8")
        for token in (
            "attempt_62011285",
            "attempt_62021985",
            "q06_1p45A/window.rst7",
            "EXPECTED_QM_ATOMS = 146",
            "N_ALPHA = 8949",
            "THR267_OG1 = 8960",
            "TRANSFERRED_HG1 = 8961",
            "L2_C12 = 10287",
            "L2_N3 = 10289",
            "548dca57cb623e69d5fefaa8e270267a5bd599bd5995113219cdf349c346e56c",
            "a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0",
        ):
            self.assertIn(token, text)

    def test_bond_graph_and_qm_membership_are_fail_closed(self):
        text = AUDIT.read_text(encoding="utf-8")
        for token in (
            "nalpha_hg1_bond_present",
            "og1_hg1_bond_absent",
            "c12_n3_bond_present",
            "all_reaction_atoms_in_qm",
            "FAIL_A1_PT2_ATOM_OR_QM_CONTRACT",
        ):
            self.assertIn(token, text)

    def test_preorganization_gate_is_explicit_and_not_mechanism_proof(self):
        text = AUDIT.read_text(encoding="utf-8")
        for token in (
            "nalpha_n3_A",
            "hg1_n3_A",
            "nalpha_hg1_n3_deg",
            "3.5",
            "2.5",
            "135.0",
            "PASS_A1_PT2_DIRECT_RELAY_PREORGANIZED",
            "NOT_EVALUATED_A1_PT2_DIRECT_RELAY_NOT_PREORGANIZED",
            "NOT_EVALUATED_TS_PMF_BARRIER_MECHANISM",
        ):
            self.assertIn(token, text)


if __name__ == "__main__":
    unittest.main()
