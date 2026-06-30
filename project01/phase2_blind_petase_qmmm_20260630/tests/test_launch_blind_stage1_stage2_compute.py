import subprocess
import sys
import tempfile
import textwrap
import unittest
from pathlib import Path


class LaunchBlindStage1Stage2ComputeTest(unittest.TestCase):
    def test_launcher_records_runner_evidence_and_preserves_exit_code(self):
        launcher = Path("work/launch_blind_stage1_stage2_compute.py")
        if not launcher.exists():
            launcher = Path("project01/phase2_blind_petase_qmmm_20260630/scripts/launch_blind_stage1_stage2_compute.py")
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            phase_root = tmpdir / "project01" / "phase2_blind_petase_qmmm_20260630"
            phase_root.mkdir(parents=True)
            runner = tmpdir / "fake_runner.py"
            runner.write_text(
                textwrap.dedent(
                    """\
                    import sys
                    print("fake runner stdout")
                    print("fake runner stderr", file=sys.stderr)
                    raise SystemExit(2)
                    """
                ),
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(launcher),
                    "--project-root",
                    str(tmpdir),
                    "--phase-subdir",
                    "project01/phase2_blind_petase_qmmm_20260630",
                    "--runner",
                    str(runner),
                    "--skip-shell-probes",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 2, result.stderr)
            status_dir = phase_root / "blind_work" / "00_run_status"
            summary = status_dir / "compute_launch_summary.md"
            stdout_log = status_dir / "stage1_stage2_runner.stdout.log"
            stderr_log = status_dir / "stage1_stage2_runner.stderr.log"
            command_log = status_dir / "stage1_stage2_runner_command.txt"

            self.assertTrue(summary.exists())
            self.assertTrue(stdout_log.exists())
            self.assertTrue(stderr_log.exists())
            self.assertTrue(command_log.exists())
            self.assertIn("fake runner stdout", stdout_log.read_text(encoding="utf-8"))
            self.assertIn("fake runner stderr", stderr_log.read_text(encoding="utf-8"))

            summary_text = summary.read_text(encoding="utf-8").lower()
            self.assertIn("exit code: 2", summary_text)
            self.assertIn("blind boundary", summary_text)
            self.assertIn("do not use article", summary_text)
            self.assertNotIn("password", summary_text)
            self.assertIn("--skip-shell-probes", command_log.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()