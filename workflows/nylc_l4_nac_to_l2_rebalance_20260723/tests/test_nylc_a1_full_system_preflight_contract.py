#!/usr/bin/env python3
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
RUNNER = HERE / "scripts" / "run_nylc_a1_full_system_preflight.sh"
AUDITOR = HERE / "scripts" / "audit_nylc_a1_full_system.py"
SBATCH = HERE / "slurm" / "run_nylc_a1_full_system_preflight.sbatch"


class A1FullSystemPreflightContract(unittest.TestCase):
    def test_runner_and_auditor_exist(self):
        self.assertTrue(RUNNER.is_file())
        self.assertTrue(AUDITOR.is_file())
        self.assertTrue(SBATCH.is_file())

    @unittest.skipUnless(RUNNER.is_file(), "runner not implemented")
    def test_runner_enforces_three_independent_scientific_preflights(self):
        text = RUNNER.read_text(encoding="utf-8")
        for token in [
            "nac_evt18_time1206ps", "nac_evt25_time1462ps", "nac_evt08_time1086ps",
            "job_61902961", "PASS_A1_SCREENING_PATCH",
            "build_nylc_a1_full_system.py", "audit_nylc_a1_full_system.py",
            "grompp", "-maxwarn 0", "mdrun", "gmx_single.edr", "potential.xvg",
            "run_history.tsv", "run_history.jsonl",
        ]:
            self.assertIn(token, text)
        self.assertNotIn("/Gromacs-DCU2/", text)

    @unittest.skipUnless(AUDITOR.is_file(), "auditor not implemented")
    def test_auditor_contract_names_technical_gates(self):
        text = AUDITOR.read_text(encoding="utf-8")
        for token in [
            "PASS_A1_FULL_SYSTEM_PREFLIGHT", "minimum_nonbonded_distance_nm",
            "nalpha_hydrogen_count", "ogamma_hydrogen_count",
            "gromacs_potential_energy_kj_mol", "grompp_maxwarn_zero",
        ]:
            self.assertIn(token, text)

    @unittest.skipUnless(SBATCH.is_file(), "sbatch not implemented")
    def test_preflight_is_three_way_scnet_array(self):
        text = SBATCH.read_text(encoding="utf-8")
        for token in [
            "#SBATCH -p xahcnormal", "#SBATCH --array=0-2%3",
            "run_nylc_a1_full_system_preflight.sh", "SLURM_ARRAY_TASK_ID",
        ]:
            self.assertIn(token, text)


if __name__ == "__main__":
    unittest.main()
