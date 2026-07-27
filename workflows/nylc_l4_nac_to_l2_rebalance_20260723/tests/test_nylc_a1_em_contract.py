#!/usr/bin/env python3
import importlib.util
import unittest
from pathlib import Path
HERE=Path(__file__).resolve().parents[1]
RUNNER=HERE/"scripts"/"run_nylc_a1_em.sh"
AUDITOR=HERE/"scripts"/"audit_nylc_a1_em.py"
SBATCH=HERE/"slurm"/"run_nylc_a1_em.sbatch"

class A1EMContract(unittest.TestCase):
    def test_artifacts_exist(self):
        self.assertTrue(RUNNER.is_file())
        self.assertTrue(AUDITOR.is_file())
        self.assertTrue(SBATCH.is_file())

    @unittest.skipUnless(RUNNER.is_file(),"runner missing")
    def test_runner_has_restrained_then_free_em_and_independent_inputs(self):
        text=RUNNER.read_text()
        for token in ["61966928","-DPOSRES","-DPOSRES_L2_1000","em_hrelax",
                      "em_free","-maxwarn 0","run_history.tsv","run_history.jsonl",
                      "audit_nylc_a1_em.py"]:
            self.assertIn(token,text)
        self.assertGreaterEqual(text.count("integrator = steep"),2)

    @unittest.skipUnless(AUDITOR.is_file(),"auditor missing")
    def test_em_auditor_requires_convergence_and_no_overlap(self):
        text=AUDITOR.read_text()
        for token in ["PASS_A1_EM","maximum_force_kj_mol_nm",
                      "minimum_nonbonded_distance_nm","converged to Fmax",
                      "nalpha_hydrogen_count","ogamma_hydrogen_count"]:
            self.assertIn(token,text)


    @unittest.skipUnless(AUDITOR.is_file(),"auditor missing")
    def test_nan_detection_ignores_warnangle_but_rejects_numeric_nan(self):
        spec=importlib.util.spec_from_file_location("a1_em_audit",AUDITOR)
        module=importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        self.assertFalse(module.log_has_nonfinite_token("lincs-warnangle = 30"))
        self.assertTrue(module.log_has_nonfinite_token("Potential Energy = nan"))

    @unittest.skipUnless(SBATCH.is_file(),"sbatch missing")
    def test_em_is_three_way_array(self):
        text=SBATCH.read_text()
        for token in ["#SBATCH -p xahcnormal","#SBATCH --array=0-2%3",
                      "run_nylc_a1_em.sh","SLURM_ARRAY_TASK_ID"]:
            self.assertIn(token,text)

if __name__=="__main__":
    unittest.main()
