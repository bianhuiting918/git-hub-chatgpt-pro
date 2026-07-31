#!/usr/bin/env python3
"""Regression contract for Thr267 covalent integrity in A1 audits."""
from __future__ import annotations

import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = (
    ROOT / "scripts"
    / "prepare_audit_nylc_a1_step1_raw_nac_preorganized_calibration.py"
)
SPEC = importlib.util.spec_from_file_location("_a1_preorg", DRIVER)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot import {DRIVER}")
MODULE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MODULE)


class Thr267ChemicalIntegrityTests(unittest.TestCase):
    def test_accepts_intact_a1_thr267(self) -> None:
        result = MODULE.thr267_chemical_integrity(
            {
                "nalpha_ca_A": 1.46,
                "ca_cb_A": 1.54,
                "cb_og1_A": 1.43,
                "nalpha_h1_A": 1.02,
                "nalpha_h2_A": 1.02,
                "nalpha_hg1_A": 1.03,
            }
        )
        self.assertTrue(result["pass"])
        self.assertTrue(all(result["checks"].values()))

    def test_rejects_broken_heavy_atom_bond(self) -> None:
        result = MODULE.thr267_chemical_integrity(
            {
                "nalpha_ca_A": 1.46,
                "ca_cb_A": 3.20,
                "cb_og1_A": 1.43,
                "nalpha_h1_A": 1.02,
                "nalpha_h2_A": 1.02,
                "nalpha_hg1_A": 1.03,
            }
        )
        self.assertFalse(result["pass"])
        self.assertFalse(result["checks"]["ca_cb_covalent"])

    def test_rejects_migrated_a1_proton(self) -> None:
        result = MODULE.thr267_chemical_integrity(
            {
                "nalpha_ca_A": 1.46,
                "ca_cb_A": 1.54,
                "cb_og1_A": 1.43,
                "nalpha_h1_A": 2.10,
                "nalpha_h2_A": 1.02,
                "nalpha_hg1_A": 1.03,
            }
        )
        self.assertFalse(result["pass"])
        self.assertFalse(result["checks"]["nalpha_h1_covalent"])


if __name__ == "__main__":
    unittest.main()
