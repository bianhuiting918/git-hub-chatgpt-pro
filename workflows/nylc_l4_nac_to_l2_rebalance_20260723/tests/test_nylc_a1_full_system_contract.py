#!/usr/bin/env python3
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
BUILDER = HERE / "scripts" / "build_nylc_a1_full_system.py"


class A1FullSystemBuilderPresenceTest(unittest.TestCase):
    def test_full_system_builder_exists(self):
        self.assertTrue(BUILDER.is_file(), str(BUILDER))


if __name__ == "__main__":
    unittest.main()
