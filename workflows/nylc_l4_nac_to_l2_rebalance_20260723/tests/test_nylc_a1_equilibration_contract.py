#!/usr/bin/env python3
import unittest
from pathlib import Path

HERE=Path(__file__).resolve().parents[1]
RUNNER=HERE/"scripts"/"run_nylc_a1_equilibration.sh"
SBATCH=HERE/"slurm"/"run_nylc_a1_equilibration.sbatch"

class A1EquilibrationContract(unittest.TestCase):
    def test_files_exist(self):
        self.assertTrue(RUNNER.is_file())
        self.assertTrue(SBATCH.is_file())

    @unittest.skipUnless(RUNNER.is_file(),"runner missing")
    def test_runner_gates_em_and_runs_release_plus_one_ns_free(self):
        text=RUNNER.read_text()
        for token in ["PASS_A1_EM","nvt50","nvt150","nvt300","npt300r","npt300rel",
                      "npt300free","fully_unrestrained_NPT_1ns","free_tpr_contract.json",
                      "run_history.tsv","run_history.jsonl","refusing to overwrite",
                      "lincs_warning","settle_problem"]:
            self.assertIn(token,text)

    @unittest.skipUnless(SBATCH.is_file(),"sbatch missing")
    def test_sbatch_uses_dcu_array_and_immutable_snapshot(self):
        text=SBATCH.read_text()
        for token in ["#SBATCH -p xahdnormal","#SBATCH --gres=dcu:1",
                      "#SBATCH --array=0-8%9","code_snapshots","sha256sum",
                      "A1_EQ_CODE_ROOT","SLURM_ARRAY_TASK_ID"]:
            self.assertIn(token,text)

if __name__=="__main__":
    unittest.main()
