#!/usr/bin/env python3
import pathlib
import unittest

WORKFLOW_ROOT = pathlib.Path(__file__).resolve().parents[1]
REVERSE = (
    WORKFLOW_ROOT
    / "scripts/run_nylc_a1_step1_product_reverse_boundary.sh"
).read_text(encoding="utf-8")
STEP2 = (
    WORKFLOW_ROOT
    / "scripts/run_nylc_a1_step2_qmwater_endpoint.sh"
).read_text(encoding="utf-8")


class LaunchEnvironmentContract(unittest.TestCase):
    def test_reverse_runner_expands_slurm_and_runtime_variables(self):
        expected = (
            'CODE_ROOT="${A1_PRODUCT_REVERSE_CODE_ROOT:?set immutable code root}"',
            'GITHUB_COMMIT="${A1_PRODUCT_REVERSE_GITHUB_COMMIT:?set immutable GitHub commit}"',
            'INDEX="${SLURM_ARRAY_TASK_ID:?run as array task 0 or 1}"',
            'ARRAY_JOB="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"',
            'SCRATCH_ROOT="${SLURM_TMPDIR:-/tmp}/nylc_a1_product_reverse_$ATTEMPT"',
            '-np "${SLURM_NTASKS:-8}"',
        )
        for contract in expected:
            with self.subTest(contract=contract):
                self.assertIn(contract, REVERSE)
        self.assertNotIn(r"\${", REVERSE)

    def test_step2_runner_uses_scnet_scratch_fallback(self):
        self.assertIn(
            'SCRATCH_ROOT="${SLURM_TMPDIR:-/tmp}/'
            'nylc_a1_step2_qmwater_$ATTEMPT"',
            STEP2,
        )
        self.assertNotIn("${SLURM_TMPDIR:?", STEP2)


if __name__ == "__main__":
    unittest.main()
