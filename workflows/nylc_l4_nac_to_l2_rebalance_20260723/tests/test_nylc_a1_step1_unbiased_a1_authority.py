#!/usr/bin/env python3
"""Contract tests for the two-seed unrestrained A1 authority diagnostic."""
from __future__ import annotations

import importlib.util
import pathlib
import unittest

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "prepare_audit_nylc_a1_step1_unbiased_a1_authority.py"
RUNNER = ROOT / "scripts" / "run_nylc_a1_step1_unbiased_a1_authority.sh"
SBATCH = ROOT / "slurm" / "run_nylc_a1_step1_unbiased_a1_authority.sbatch"


class UnbiasedA1AuthorityContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.assertTrue = unittest.TestCase().assertTrue
        cls.assertTrue(DRIVER.is_file(), f"missing driver {DRIVER}")
        spec = importlib.util.spec_from_file_location("_a1_unbiased", DRIVER)
        if spec is None or spec.loader is None:
            raise RuntimeError(f"cannot import {DRIVER}")
        cls.module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.module)

    def test_two_seed_fixed_source_contract(self) -> None:
        description = self.module.describe()
        self.assertEqual(description["array_task_count"], 2)
        self.assertEqual(description["mpi_ranks_per_task"], 8)
        self.assertEqual(description["unrestrained_minimization_blocks"], 4)
        self.assertEqual(
            description["original_raw_nac_restart_shas"],
            {
                "26723": "2440de548c385f092c37f683de7b379ff5b18b6dc16593f7dbc80a9a8a167e14",
                "26737": "48c3944d3295158b06e96e32e4e07d9bcae9ceba2731f95aae9c1335ff972abe",
            },
        )

    def test_input_is_reaction_coordinate_free(self) -> None:
        text = self.module.unbiased_minimization_input(
            {"task_index": 0}, 0, "@1,@2"
        )
        self.assertIn("ntr=0", text)
        self.assertIn("nmropt=0", text)
        self.assertIn("drms=0.01", text)
        self.assertNotIn("DISANG", text)
        self.assertNotIn("&rst", text)
        self.assertNotIn("restraintmask", text)

    def test_chemical_failure_reclassifies_old_array(self) -> None:
        self.assertEqual(
            self.module.reclassification_status([False] * 18),
            "NOT_EVALUATED_TECHNICAL_CHEMICAL_INTEGRITY",
        )
        self.assertEqual(
            self.module.reclassification_status([True] * 18),
            "PASS_ALL_TASKS_THR267_CHEMICAL_INTEGRITY",
        )

    def test_runner_and_sbatch_resources(self) -> None:
        runner = RUNNER.read_text(encoding="utf-8")
        sbatch = SBATCH.read_text(encoding="utf-8")
        self.assertIn("for STAGE in 0 1 2 3", runner)
        self.assertNotIn("restraints.RST", runner)
        self.assertIn("#SBATCH -n 8", sbatch)
        self.assertIn("#SBATCH --array=0-1", sbatch)
        self.assertNotIn("%", next(
            line for line in sbatch.splitlines() if "--array=" in line
        ))


if __name__ == "__main__":
    unittest.main()
