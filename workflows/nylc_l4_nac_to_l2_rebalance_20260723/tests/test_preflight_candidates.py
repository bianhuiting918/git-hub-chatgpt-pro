#!/usr/bin/env python3
import sys
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE / "scripts"))

from preflight_candidates import grompp_command, scan_diagnostics  # noqa: E402


class PreflightTests(unittest.TestCase):
    def test_grompp_is_strict_and_uses_rebuilt_reference(self):
        build = Path("/task/candidates/c/build")
        out = Path("/task/candidates/c/preflight")
        cmd = grompp_command("/gmx", Path("/flow/mdp/em.mdp"), build, out)
        self.assertEqual(cmd[0:2], ["/gmx", "grompp"])
        self.assertIn("-maxwarn", cmd)
        self.assertEqual(cmd[cmd.index("-maxwarn") + 1], "0")
        self.assertEqual(cmd[cmd.index("-c") + 1], str(build / "rebuilt.gro"))
        self.assertEqual(cmd[cmd.index("-r") + 1], str(build / "rebuilt.gro"))
        self.assertEqual(cmd[cmd.index("-p") + 1], str(build / "topol.top"))

    def test_diagnostic_scanner_separates_warning_classes(self):
        text = """LINCS WARNING
Fatal error:
Water molecule cannot be settled
energy is NaN
"""
        observed = scan_diagnostics(text)
        self.assertEqual(observed["lincs_warning"], 1)
        self.assertEqual(observed["fatal"], 1)
        self.assertEqual(observed["settle_problem"], 1)
        self.assertEqual(observed["nan"], 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
