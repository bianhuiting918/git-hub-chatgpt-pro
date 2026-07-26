#!/usr/bin/env python3
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
INSTALLER = HERE / "scripts" / "install_nylc_amberclassic.sh"
STAGER = HERE / "scripts" / "stage_nylc_amberclassic_source.sh"
SBATCH = HERE / "slurm" / "run_install_nylc_amberclassic.sbatch"


class AmberClassicInstallContractTests(unittest.TestCase):
    def test_installer_is_pinned_task_local_and_offline(self):
        text = INSTALLER.read_text(encoding="utf-8")
        for token in [
            "0b35bfeb96026ffa4e5876391a0828f39b3cfc8d",
            "STAGED_ARCHIVE",
            "compiler/gcc/11.4.0",
            "gcc --version",
            "gfortran --version",
            'TOOLS_ROOT="$TASK_ROOT/tools"',
            "CONNECT.TPL",
            "antechamber -h",
            "run_history.tsv",
            "run_history.jsonl",
            "PASS.json",
        ]:
            self.assertIn(token, text)
        self.assertIn("l4_nac_to_l2_rebalance_20260723", text)
        self.assertNotIn("codeload.github.com", text)
        self.assertNotIn("curl ", text)

    def test_login_node_stager_is_pinned_audited_and_atomic(self):
        text = STAGER.read_text(encoding="utf-8")
        for token in [
            "codeload.github.com/Amber-MD/AmberClassic",
            "0b35bfeb96026ffa4e5876391a0828f39b3cfc8d",
            "source_archives",
            "mktemp",
            "sha256sum",
            "tar -tzf",
            "mv -n",
            "run_history.tsv",
            "run_history.jsonl",
            "PASS_SOURCE.json",
        ]:
            self.assertIn(token, text)

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
