#!/usr/bin/env python3
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
INSTALLER = HERE / "scripts" / "install_nylc_amberclassic.sh"
SBATCH = HERE / "slurm" / "run_install_nylc_amberclassic.sbatch"


class AmberClassicInstallContractTests(unittest.TestCase):
    def test_installer_is_pinned_and_task_local(self):
        text = INSTALLER.read_text(encoding="utf-8")
        for token in [
            "0b35bfeb96026ffa4e5876391a0828f39b3cfc8d",
            "codeload.github.com/Amber-MD/AmberClassic",
            "CONNECT.TPL",
            "antechamber -h",
            "run_history.tsv",
            "run_history.jsonl",
            "PASS.json",
        ]:
            self.assertIn(token, text)
        self.assertIn("l4_nac_to_l2_rebalance_20260723/tools", text)

    def test_install_runs_on_scnet_compute_node(self):
        text = SBATCH.read_text(encoding="utf-8")
        for token in [
            "#SBATCH -p xahcnormal",
            "#SBATCH -c 8",
            "install_nylc_amberclassic.sh",
            "SLURM_JOB_ID",
        ]:
            self.assertIn(token, text)


if __name__ == "__main__":
    unittest.main()
