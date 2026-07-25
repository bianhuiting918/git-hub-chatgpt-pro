import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "append_run_history.py"
COLUMNS = [
    "timestamp_local",
    "timestamp_utc",
    "git_commit",
    "run_id",
    "script",
    "parameters",
    "input_manifest",
    "input_sha256",
    "output_path",
    "scheduler_job_id",
    "exit_status",
    "n_input",
    "n_evaluated",
    "n_failed",
    "n_not_evaluated",
    "summary",
]


def base_command(history: Path, run_id: str) -> list[str]:
    return [
        sys.executable,
        str(SCRIPT),
        "--history",
        str(history),
        "--timestamp-local",
        "2026-07-26T04:00:00+08:00",
        "--timestamp-utc",
        "2026-07-25T20:00:00Z",
        "--git-commit",
        "abc123",
        "--run-id",
        run_id,
        "--script",
        "slurm/example.sbatch",
        "--parameters",
        "RUN_ID=" + run_id,
        "--input-manifest",
        "manifest.tsv",
        "--input-sha256",
        "deadbeef",
        "--output-path",
        "results/" + run_id,
        "--scheduler-job-id",
        "123",
        "--exit-status",
        "COMPLETED_0:0",
        "--n-input",
        "3",
        "--n-evaluated",
        "3",
        "--n-failed",
        "0",
        "--n-not-evaluated",
        "0",
        "--summary",
        "technical pass only",
    ]


class AppendRunHistoryTests(unittest.TestCase):
    def test_creates_header_and_appends_distinct_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            history = Path(tmp) / "run_history.tsv"
            first = subprocess.run(base_command(history, "run_a"), text=True, capture_output=True)
            second = subprocess.run(base_command(history, "run_b"), text=True, capture_output=True)
            self.assertEqual(first.returncode, 0, first.stderr)
            self.assertEqual(second.returncode, 0, second.stderr)
            with history.open(newline="", encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(list(rows[0]), COLUMNS)
            self.assertEqual([row["run_id"] for row in rows], ["run_a", "run_b"])

    def test_duplicate_run_id_is_rejected_without_modification(self):
        with tempfile.TemporaryDirectory() as tmp:
            history = Path(tmp) / "run_history.tsv"
            first = subprocess.run(base_command(history, "run_a"), text=True, capture_output=True)
            self.assertEqual(first.returncode, 0, first.stderr)
            before = history.read_bytes()
            duplicate = subprocess.run(base_command(history, "run_a"), text=True, capture_output=True)
            self.assertNotEqual(duplicate.returncode, 0)
            self.assertIn("duplicate run_id", duplicate.stderr)
            self.assertEqual(history.read_bytes(), before)


if __name__ == "__main__":
    unittest.main()
