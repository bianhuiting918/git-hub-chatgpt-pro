#!/usr/bin/env python3
import pathlib
import unittest

HERE = pathlib.Path(__file__).resolve().parents[1]
PREPARE = HERE / "scripts" / "prepare_nylc_a1_unified_step1_qattack_extension.py"
AUDIT = HERE / "scripts" / "audit_nylc_a1_unified_step1_qattack_extension.py"
RUNNER = HERE / "scripts" / "run_nylc_a1_unified_step1_qattack_extension.sh"
SLURM = HERE / "slurm" / "run_nylc_a1_unified_step1_qattack_extension.sbatch"


class UnifiedStep1QAttackExtensionContractTests(unittest.TestCase):
    def test_extension_files_exist(self):
        for path in (PREPARE, AUDIT, RUNNER, SLURM):
            self.assertTrue(path.is_file(), path)

    def test_extension_is_pinned_to_k200_endpoint(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "attempt_62016429_2_k200",
            "q05_1p65A/window.rst7",
            "80ae4121248f90d8488125ce69d51aa31198d19970f63bbe04fe65ee57503081",
            "6406421d86ae35330cc6a480e5f5b234f9ff8591e97d55743e354095b42a5e74",
            "PASS_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_RECOVERY",
            "NOT_EVALUATED_ATTACK_BRACKET_NOT_REACHED",
        ):
            self.assertIn(token, text)

    def test_extension_changes_only_attack_target(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "TARGET_A = 1.45",
            "FORCE_CONSTANT = 200.0",
            "EXPECTED_QM_ATOMS = 146",
            "EXPECTED_ELECTRONS = 510",
            "QMCHARGE = 0",
            '"step1_qm_water_count": 0',
            "maxcyc=150",
            "ncyc=150",
            "ntmin=2",
            "dx0=0.005",
            "r4=4.500",
            "r2=95.0",
            "r3=115.0",
        ):
            self.assertIn(token, text)

    def test_audit_keeps_technical_and_scientific_status_separate(self):
        text = AUDIT.read_text(encoding="utf-8")
        for token in (
            "PASS_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_EXTENSION",
            "FAIL_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_EXTENSION",
            "PASS_CONSTRAINED_ATTACK_BRACKET_SEED",
            "NOT_EVALUATED_ATTACK_BRACKET_NOT_REACHED",
            "NOT_EVALUATED_TS_PMF_BARRIER",
            "Convergence could not be achieved",
            "SANDER BOMB",
            "bond_overflow",
            "og1_c12_A",
            "c12_o2_A",
            "c12_n3_A",
        ):
            self.assertIn(token, text)

    def test_runner_is_fail_closed_and_uses_separate_output(self):
        runner = RUNNER.read_text(encoding="utf-8")
        slurm = SLURM.read_text(encoding="utf-8")
        for token in (
            "a1_unified_step1_qattack_extension",
            "mpirun --bind-to none",
            "SLURM_NTASKS",
            "run_history.tsv",
            "run_history.jsonl",
            "NOT_EVALUATED.json",
            "PASS.json",
            "sha256sum",
        ):
            self.assertIn(token, runner)
        for token in (
            "#SBATCH -n 8",
            "#SBATCH --mem-per-cpu=2500M",
            "#SBATCH -t 01:00:00",
            "A1_QATTACK_EXTENSION_CODE_SOURCE",
            "A1_QATTACK_EXTENSION_GITHUB_COMMIT",
        ):
            self.assertIn(token, slurm)
        self.assertNotIn("#SBATCH --gres", slurm)
        self.assertNotIn("pmemd.cuda", runner)


if __name__ == "__main__":
    unittest.main()
