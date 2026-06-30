import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class RunBlindStage1Stage2GatesTest(unittest.TestCase):
    def test_runner_writes_blocked_status_when_required_inputs_are_missing(self):
        script = Path("work/run_blind_stage1_stage2_gates.py")
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            out_root = tmpdir / "blind_work"

            result = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--project-root",
                    str(tmpdir),
                    "--out-root",
                    str(out_root),
                    "--skip-shell-probes",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 2, result.stderr)
            status_path = out_root / "00_run_status" / "stage1_stage2_gate_status.tsv"
            instructions_path = out_root / "00_run_status" / "stage1_stage2_next_actions.md"
            self.assertTrue(status_path.exists())
            self.assertTrue(instructions_path.exists())

            with status_path.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(rows[0]["gate"], "ligand_smiles")
            self.assertEqual(rows[0]["status"], "blocked")
            self.assertEqual(rows[0]["required_before"], "ligand_build")
            self.assertEqual(rows[-1]["gate"], "stage2_md_queue")
            self.assertEqual(rows[-1]["status"], "blocked")
            self.assertIn("accepted_gs_pose", rows[-1]["note"])

            lowered = instructions_path.read_text(encoding="utf-8").lower()
            self.assertIn("blind", lowered)
            self.assertIn("do not use article", lowered)
            self.assertNotIn("paper ts", lowered)
            self.assertNotIn("paper trajectory", lowered)


if __name__ == "__main__":
    unittest.main()