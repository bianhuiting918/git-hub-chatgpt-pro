#!/usr/bin/env python3
import pathlib
import unittest

HERE = pathlib.Path(__file__).resolve().parents[1]
PREPARE = HERE / "scripts" / "prepare_nylc_a1_step1_release_test.py"
AUDIT = HERE / "scripts" / "audit_nylc_a1_step1_release_test.py"
RUNNER = HERE / "scripts" / "run_nylc_a1_step1_release_test.sh"
SLURM = HERE / "slurm" / "run_nylc_a1_step1_release_test.sbatch"


class A1Step1ReleaseTestContract(unittest.TestCase):
    def test_release_files_exist(self):
        for path in (PREPARE, AUDIT, RUNNER, SLURM):
            self.assertTrue(path.is_file(), path)

    def test_source_is_exact_constrained_seed(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "attempt_62021985",
            "q06_1p45A/window.rst7",
            "3083663f502f522974f7eac56e1b3a4c7015933c7d949562364d359006c25224",
            "548dca57cb623e69d5fefaa8e270267a5bd599bd5995113219cdf349c346e56c",
            "PASS_CONSTRAINED_ATTACK_BRACKET_SEED",
            "PASS_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_EXTENSION",
        ):
            self.assertIn(token, text)

    def test_release_keeps_hamiltonian_and_removes_reactive_restraints(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "EXPECTED_QM_ATOMS = 146",
            "EXPECTED_ELECTRONS = 510",
            "QMCHARGE = 0",
            '"step1_qm_water_count": 0',
            '"local_release"',
            '"full_release"',
            "maxcyc=500",
            "ncyc=100",
            "ntmin=1",
            "qm_theory='DFTB3'",
            "dftb_telec=200.0",
        ):
            self.assertIn(token, text)
        self.assertNotIn("DISANG=", text)
        self.assertNotIn("nmropt=1", text)

    def test_audit_reports_contact_and_carbonyl_pyramidalization(self):
        text = AUDIT.read_text(encoding="utf-8")
        for token in (
            "C12_C11",
            "carbonyl_neighbor_angle_sum_deg",
            "carbonyl_pyramidalization_deg",
            "PASS_RELEASE_RETAINED_ATTACK_CONTACT",
            "FAIL_RELEASE_RETURNED_TOWARD_REACTANT",
            "NOT_EVALUATED_RELEASE_AMBIGUOUS",
            "NOT_EVALUATED_TS_PMF_BARRIER",
            "Convergence could not be achieved",
            "SANDER BOMB",
            "bond_overflow",
        ):
            self.assertIn(token, text)

    def test_runner_is_fail_closed_and_uses_new_output(self):
        runner = RUNNER.read_text(encoding="utf-8")
        slurm = SLURM.read_text(encoding="utf-8")
        for token in (
            "a1_step1_release_test",
            "local_release",
            "full_release",
            "mpirun --bind-to none",
            "run_history.tsv",
            "run_history.jsonl",
            "NOT_EVALUATED.json",
            "PASS.json",
        ):
            self.assertIn(token, runner)
        for token in (
            "#SBATCH -n 8",
            "#SBATCH --mem-per-cpu=2500M",
            "#SBATCH -t 02:00:00",
            "A1_STEP1_RELEASE_CODE_SOURCE",
            "A1_STEP1_RELEASE_GITHUB_COMMIT",
        ):
            self.assertIn(token, slurm)
        self.assertNotIn("#SBATCH --gres", slurm)
        self.assertNotIn("pmemd.cuda", runner)


if __name__ == "__main__":
    unittest.main()
