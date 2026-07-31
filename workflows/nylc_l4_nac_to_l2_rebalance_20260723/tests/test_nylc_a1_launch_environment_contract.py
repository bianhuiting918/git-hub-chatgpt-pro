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
NEARMISS_RUNNER = (
    WORKFLOW_ROOT
    / "scripts/run_nylc_a1_step2_water_network_nearmiss_continuation.sh"
)
NEARMISS_SBATCH = (
    WORKFLOW_ROOT
    / "slurm/run_nylc_a1_step2_water_network_nearmiss_continuation.sbatch"
)
WATERNET_DRIVER = (
    WORKFLOW_ROOT
    / "scripts/prepare_audit_nylc_a1_step2_water_network_sampling.py"
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

    def test_step2_sparse_hash_manifest_does_not_fail_under_pipefail(self):
        self.assertNotIn(
            '[[ -f "$name" ]] && sha256sum "$name"',
            STEP2,
        )
        self.assertIn(
            'if [[ -f "$name" ]]; then\n'
            '                sha256sum "$name"\n'
            '            fi',
            STEP2,
        )

    def test_step2_nearmiss_continuation_is_sha_fixed_and_unrestrained(self):
        self.assertTrue(NEARMISS_RUNNER.is_file())
        self.assertTrue(NEARMISS_SBATCH.is_file())
        runner = NEARMISS_RUNNER.read_text(encoding="utf-8")
        sbatch = NEARMISS_SBATCH.read_text(encoding="utf-8")
        for token in (
            "SOURCE_ATTEMPT=attempt_62425747_2",
            "309e419b7f14c73c9388fea7325d9e54e2434669019e24b209bf0052981058a4",
            "SOURCE_ATTEMPT=attempt_62425747_6",
            "42c46ebe61ad3016c86a91880ac1d6f94d7f7c7fdcda385bd89935277d378583",
            "$SOURCE_ATTEMPT/direct_near_miss_0.rst7",\n            "a1_step2_water_network_nearmiss_continuation",
            'mpirun --bind-to none -np 8 sander.MPI',
            '-c "$SOURCE_RST7"',
            "prepare_audit_nylc_a1_step2_water_network_sampling.py",
        ):
            with self.subTest(token=token):
                self.assertIn(token, runner)
        self.assertIn("#SBATCH --array=0-7", sbatch)
        self.assertIn("#SBATCH -n 8", sbatch)
        self.assertNotIn("%", sbatch.split("#SBATCH --array=", 1)[1].splitlines()[0])
        self.assertIn("ntr=0, nmropt=0", WATERNET_DRIVER)
        self.assertIn("MIN_CONSECUTIVE_FRAMES = 5", WATERNET_DRIVER)


if __name__ == "__main__":
    unittest.main()
