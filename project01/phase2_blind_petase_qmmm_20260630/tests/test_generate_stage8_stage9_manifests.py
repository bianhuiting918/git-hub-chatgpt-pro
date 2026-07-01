import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class GenerateStage8Stage9ManifestsTest(unittest.TestCase):
    def test_generates_free_energy_and_final_validation_manifests(self):
        script = Path("work/generate_stage8_stage9_manifests.py")
        if not script.exists():
            script = Path("project01/phase2_blind_petase_qmmm_20260630/scripts/generate_stage8_stage9_manifests.py")
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            accepted = tmpdir / "accepted_ts_manifest.tsv"
            accepted.write_text(
                "\n".join(
                    [
                        "accepted_ts_id\tensemble_candidate_id\tpath_id\treaction_stage\tstructure_path\tcommittor_mean\tacceptance_status",
                        "ATS_AC_001\tAC_TS_0001\tAC_DIRECT_HIS_GENERAL_BASE\tacylation\twork/ac_ts_001.pdb\t0.51\taccepted",
                        "ATS_DE_001\tDE_TS_0001\tDE_WATER_ATTACK_HIS_GENERAL_BASE\tdeacylation\twork/de_ts_001.pdb\t0.49\taccepted",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            out_root = tmpdir / "blind_work"

            result = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--accepted-ts-manifest",
                    str(accepted),
                    "--out-root",
                    str(out_root),
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)

            ac_pmf = out_root / "08_free_energy" / "acylation_pmf.tsv"
            de_pmf = out_root / "08_free_energy" / "deacylation_pmf.tsv"
            barrier = out_root / "08_free_energy" / "barrier_summary.md"
            compare = out_root / "09_paper_validation" / "blind_vs_paper_comparison.md"
            audit = out_root / "09_paper_validation" / "discrepancy_audit.md"

            for path in (ac_pmf, de_pmf, barrier, compare, audit):
                self.assertTrue(path.exists(), f"missing {path}")
                lowered = path.read_text(encoding="utf-8").lower()
                self.assertNotIn("paper barrier", lowered)
                self.assertNotIn("paper rate", lowered)
                self.assertNotIn("paper ts coordinate", lowered)

            with ac_pmf.open("r", encoding="utf-8", newline="") as handle:
                ac_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(ac_rows[0]["reaction_stage"], "acylation")
            self.assertEqual(ac_rows[0]["sampling_status"], "not_started")
            self.assertEqual(ac_rows[0]["source"], "blind_accepted_ts")

            with de_pmf.open("r", encoding="utf-8", newline="") as handle:
                de_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(de_rows[0]["reaction_stage"], "deacylation")
            self.assertEqual(de_rows[0]["sampling_status"], "not_started")

            self.assertIn("Do not open paper results", compare.read_text(encoding="utf-8"))
            self.assertIn("not_ready", audit.read_text(encoding="utf-8"))


if __name__ == "__main__":
    unittest.main()
