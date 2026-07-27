#!/usr/bin/env python3
import pathlib
import unittest

HERE = pathlib.Path(__file__).resolve().parents[1]
PREPARE = HERE / "scripts" / "prepare_nylc_a1_step1_coupled_scout.py"
AUDIT = HERE / "scripts" / "audit_nylc_a1_step1_coupled_scout.py"
RUNNER = HERE / "scripts" / "run_nylc_a1_step1_coupled_scout.sh"
SLURM = HERE / "slurm" / "run_nylc_a1_step1_coupled_scout.sbatch"


class A1Step1CoupledScoutContract(unittest.TestCase):
    def test_files_exist(self):
        for path in (PREPARE, AUDIT, RUNNER, SLURM):
            self.assertTrue(path.is_file(), path)

    def test_source_and_qm_contract_are_frozen(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "attempt_62021985",
            "q06_1p45A/window.rst7",
            "3083663f502f522974f7eac56e1b3a4c7015933c7d949562364d359006c25224",
            "548dca57cb623e69d5fefaa8e270267a5bd599bd5995113219cdf349c346e56c",
            "EXPECTED_QM_ATOMS = 146",
            "EXPECTED_ELECTRONS = 510",
            "QMCHARGE = 0",
            '"step1_qm_water_count": 0',
        ):
            self.assertIn(token, text)

    def test_three_carbonyl_targets_and_release_are_defined(self):
        text = PREPARE.read_text(encoding="utf-8")
        for token in (
            "ATTACK_TARGET_A = 1.45",
            "ATTACK_FORCE = 200.0",
            "CARBONYL_TARGETS_A = (1.30, 1.35, 1.40)",
            "CARBONYL_FORCE = 500.0",
            '"coupled"',
            '"local_release"',
            "r4=4.500",
            "r2=95.0",
            "r3=115.0",
            "maxcyc=300",
            "maxcyc=500",
        ):
            self.assertIn(token, text)

    def test_audit_separates_contact_and_exploratory_geometry_gates(self):
        text = AUDIT.read_text(encoding="utf-8")
        for token in (
            "C12_C11",
            "carbonyl_neighbor_angle_sum_deg",
            "carbonyl_pyramidalization_deg",
            "PASS_RELEASE_RETAINED_ATTACK_CONTACT",
            "FAIL_RELEASE_RETURNED_TOWARD_REACTANT",
            "PASS_PROJECT_EXPLORATORY_TETRAHEDRAL_GEOMETRY",
            "FAIL_PROJECT_EXPLORATORY_TETRAHEDRAL_GEOMETRY",
            "PASS_COUPLED_SCOUT_RELEASED_ATTACK_BASIN_CANDIDATE",
            "NOT_EVALUATED_TS_PMF_BARRIER",
            "Convergence could not be achieved",
            "SANDER BOMB",
            "bond_overflow",
        ):
            self.assertIn(token, text)

    def test_array_and_runner_are_fail_closed(self):
        runner = RUNNER.read_text(encoding="utf-8")
        slurm = SLURM.read_text(encoding="utf-8")
        for token in (
            "A1_COUPLED_SCOUT_CARBONYL_TARGET",
            "coupled",
            "local_release",
            "mpirun --bind-to none",
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
            "TARGETS=(1.30 1.35 1.40)",
            "A1_COUPLED_SCOUT_CODE_SOURCE",
            "A1_COUPLED_SCOUT_GITHUB_COMMIT",
        ):
            self.assertIn(token, slurm)
        self.assertNotIn("#SBATCH --gres", slurm)
        self.assertNotIn("pmemd.cuda", runner)


if __name__ == "__main__":
    unittest.main()
