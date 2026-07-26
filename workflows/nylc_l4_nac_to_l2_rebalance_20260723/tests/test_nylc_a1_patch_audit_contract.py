#!/usr/bin/env python3
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
AUDITOR = HERE / "scripts" / "audit_nylc_a1_patch.py"
RUNNER = HERE / "scripts" / "run_nylc_a1_parameterization.sh"


class A1PatchAuditBootstrapTests(unittest.TestCase):
    def test_independent_patch_auditor_exists(self):
        self.assertTrue(AUDITOR.is_file())

    def test_rerunnable_parameterization_driver_exists(self):
        self.assertTrue(RUNNER.is_file())


if __name__ == "__main__":
    unittest.main()
