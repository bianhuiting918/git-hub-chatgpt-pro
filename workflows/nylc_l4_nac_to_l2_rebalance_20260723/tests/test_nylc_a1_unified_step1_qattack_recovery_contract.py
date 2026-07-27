#!/usr/bin/env python3
import pathlib
import unittest

HERE = pathlib.Path(__file__).resolve().parents[1]
PREPARE = HERE / "scripts" / "prepare_nylc_a1_unified_step1_qattack_recovery.py"
AUDIT = HERE / "scripts" / "audit_nylc_a1_unified_step1_qattack_recovery.py"
RUNNER = HERE / "scripts" / "run_nylc_a1_unified_step1_qattack_recovery.sh"
SLURM = HERE / "slurm" / "run_nylc_a1_unified_step1_qattack_recovery.sbatch"


class UnifiedStep1QAttackRecoveryContractTests(unittest.TestCase):
    def test_recovery_files_exist(self):
        for path in (PREPARE, AUDIT, RUNNER, SLURM):
            self.assertTrue(path.is_file(), path)

    def test_common_source_and_frozen_qm_contract(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "attempt_62012919",
            "q03_2p65A/window.rst7",
            "EXPECTED_QM_ATOMS = 146",
            "QMCHARGE = 0",
            "EXPECTED_ELECTRONS = 510",
            '"step1_qm_water_count": 0',
        ):
            self.assertIn(token, text)

    def test_only_force_constant_varies(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "TARGETS_A = (2.65, 2.45, 2.25, 2.05, 1.85, 1.65)",
            "FORCE_CONSTANTS = (50.0, 100.0, 200.0)",
            "r4=4.500",
            "maxcyc=150",
            "ncyc=150",
            "ntmin=2",
            "dx0=0.005",
            "r1=85.0",
            "r2=95.0",
            "r3=115.0",
            "r4=125.0",
            "rk2=20.0",
            "rk3=20.0",
        ):
            self.assertIn(token, text)

    def test_audit_separates_technical_and_scientific_gates(self):
        text = AUDIT.read_text(encoding="utf-8")
        for token in (
            "PASS_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_RECOVERY",
            "FAIL_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_RECOVERY",
            "PASS_CONSTRAINED_ATTACK_BRACKET_SEED",
            "NOT_EVALUATED_ATTACK_BRACKET_NOT_REACHED",
            "NOT_EVALUATED_TS_PMF_BARRIER",
            "Convergence could not be achieved",
            "SANDER BOMB",
            "bond_overflow",
            "start_geometry",
            "endpoint_response",
        ):
            self.assertIn(token, text)

    def test_runner_and_array_are_fail_closed_and_distinct(self):
        runner = RUNNER.read_text(encoding="utf-8")
        slurm = SLURM.read_text(encoding="utf-8")
        for token in (
            "A1_QATTACK_RECOVERY_FORCE",
            "mpirun --bind-to none",
            "SLURM_NTASKS",
            "run_history.tsv",
            "run_history.jsonl",
            "NOT_EVALUATED.json",
            "PASS.json",
        ):
            self.assertIn(token, runner)
        for token in (
            "#SBATCH --array=0-2",
            "#SBATCH -n 8",
            "#SBATCH --mem-per-cpu=2500M",
            "#SBATCH -t 02:00:00",
            "FORCES=(50 100 200)",
            "A1_QATTACK_RECOVERY_CODE_SOURCE",
        ):
            self.assertIn(token, slurm)
        self.assertNotIn("#SBATCH --gres", slurm)
        self.assertNotIn("pmemd.cuda", runner)


if __name__ == "__main__":
    unittest.main()
