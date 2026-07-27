#!/usr/bin/env python3
import pathlib
import unittest

HERE = pathlib.Path(__file__).resolve().parents[1]
PREPARE = HERE / "scripts" / "prepare_nylc_a1_unified_step1_qattack.py"
AUDIT = HERE / "scripts" / "audit_nylc_a1_unified_step1_qattack.py"
RUNNER = HERE / "scripts" / "run_nylc_a1_unified_step1_qattack.sh"
SLURM = HERE / "slurm" / "run_nylc_a1_unified_step1_qattack.sbatch"


class UnifiedStep1QAttackContractTests(unittest.TestCase):
    def test_files_exist(self):
        for path in (PREPARE, AUDIT, RUNNER, SLURM):
            self.assertTrue(path.is_file(), path)

    def test_frozen_source_and_reactive_atoms(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "attempt_62011285",
            "PASS_A1_UNIFIED_CORE_DFTB3_NUMERICAL_PREFLIGHT",
            "THR267_OG1 = 8960",
            "L2_C12 = 10287",
            "L2_O2 = 10288",
            "L2_N3 = 10289",
            "EXPECTED_QM_ATOMS = 146",
            "QMCHARGE = 0",
            "EXPECTED_ELECTRONS = 510",
        ):
            self.assertIn(token, text)

    def test_attack_schedule_is_gradual_and_sd_only(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "3.27",
            "3.05",
            "2.85",
            "2.65",
            "2.45",
            "2.25",
            "2.05",
            "1.85",
            "1.65",
            "ntmin=2",
            "maxcyc=75",
            "ncyc=75",
            "dx0=0.005",
            "nmropt=1",
            "ntr=1",
            "restraint_wt=1.0",
        ):
            self.assertIn(token, text)

    def test_nac_angle_is_flat_bottom_not_exactly_forced(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "r1=85.0",
            "r2=95.0",
            "r3=115.0",
            "r4=125.0",
            "rk2=20.0",
            "rk3=20.0",
        ):
            self.assertIn(token, text)

    def test_auditor_tracks_chemistry_and_hard_errors(self):
        text = AUDIT.read_text(encoding="utf-8")
        for token in (
            "og1_c12_A",
            "o2_c12_og1_deg",
            "c12_o2_A",
            "c12_n3_A",
            "nalpha_h_A",
            "Convergence could not be achieved",
            "SANDER BOMB",
            "bond_overflow",
            "PASS_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_SCAN",
            "FAIL_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_SCAN",
            "NOT_EVALUATED_TS_PMF_BARRIER",
        ):
            self.assertIn(token, text)

    def test_runner_is_sequential_fail_closed_and_mpi_safe(self):
        text = RUNNER.read_text(encoding="utf-8")
        for token in (
            "mpirun --bind-to none",
            "SLURM_NTASKS",
            "-ref",
            "run_history.tsv",
            "run_history.jsonl",
            "NOT_EVALUATED.json",
            "PASS.json",
        ):
            self.assertIn(token, text)
        self.assertNotIn("srun -n", text)
        self.assertNotIn("pmemd.cuda", text)

    def test_slurm_uses_eight_ranks_and_no_dcu(self):
        text = SLURM.read_text(encoding="utf-8")
        for token in (
            "#SBATCH -p xahcnormal",
            "#SBATCH -n 8",
            "#SBATCH -c 1",
            "#SBATCH --mem-per-cpu=2500M",
            "#SBATCH -t 02:00:00",
            "A1_STEP1_QATTACK_CODE_SOURCE",
        ):
            self.assertIn(token, text)
        self.assertNotIn("#SBATCH --gres", text)


if __name__ == "__main__":
    unittest.main()
