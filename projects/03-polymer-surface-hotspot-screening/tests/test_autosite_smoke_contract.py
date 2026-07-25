import unittest
from pathlib import Path


SMOKE = (
    Path(__file__).resolve().parents[1]
    / "slurm"
    / "autosite_technical_smoke.sbatch"
)


class AutoSiteSmokeContractTests(unittest.TestCase):
    def test_creates_autosite_output_directory_before_invocation(self):
        text = SMOKE.read_text(encoding="utf-8")
        mkdir_token = 'mkdir -p "$candidate_dir/autosite"'
        invoke_token = '"$ADFR/autosite"'
        self.assertIn(mkdir_token, text)
        self.assertLess(text.index(mkdir_token), text.index(invoke_token))

    def test_run_directory_is_selected_by_explicit_run_id(self):
        text = SMOKE.read_text(encoding="utf-8")
        self.assertIn('RUN_ID=${RUN_ID:?RUN_ID is required}', text)
        self.assertIn('RUN="$ROOT/results/$RUN_ID"', text)


if __name__ == "__main__":
    unittest.main()
