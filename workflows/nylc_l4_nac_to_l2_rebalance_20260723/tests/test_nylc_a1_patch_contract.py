#!/usr/bin/env python3
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
SCRIPT = HERE / "scripts" / "prepare_nylc_a1_parameter_model.py"


class A1ParameterModelContractTests(unittest.TestCase):
    def test_parameter_model_builder_exists(self):
        self.assertTrue(
            SCRIPT.is_file(),
            "A1 parameter model builder must exist before parameter generation",
        )


if __name__ == "__main__":
    unittest.main()
