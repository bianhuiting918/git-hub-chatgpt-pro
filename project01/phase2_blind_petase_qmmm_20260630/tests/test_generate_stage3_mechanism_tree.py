import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class GenerateStage3MechanismTreeTest(unittest.TestCase):
    def test_generates_blind_stage3_and_stage4_files(self):
        script = Path("work/generate_stage3_mechanism_tree.py")
        with tempfile.TemporaryDirectory() as tmp:
            out_root = Path(tmp) / "blind_work"

            result = subprocess.run(
                [sys.executable, str(script), "--out-root", str(out_root)],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)

            mech = out_root / "03_mechanism_tree" / "mechanism_hypotheses.yaml"
            cvs = out_root / "03_mechanism_tree" / "candidate_cv_sets.tsv"
            screening = out_root / "04_qmmm_exploration" / "path_screening_table.tsv"
            ts_like = out_root / "04_qmmm_exploration" / "ts_like_guess_manifest.tsv"

            for path in (mech, cvs, screening, ts_like):
                self.assertTrue(path.exists(), f"missing {path}")
                text = path.read_text(encoding="utf-8")
                lowered = text.lower()
                self.assertNotIn("cv77", lowered)
                self.assertNotIn("paper_rc", lowered)
                self.assertNotIn("paper transition state", lowered)

            mech_text = mech.read_text(encoding="utf-8")
            self.assertIn("AC_DIRECT_HIS_GENERAL_BASE", mech_text)
            self.assertIn("DE_WATER_ATTACK_HIS_GENERAL_BASE", mech_text)
            self.assertIn("blind_stage3_boundary", mech_text)

            with cvs.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["reaction_stage"] for row in rows}, {"acylation", "deacylation"})
            self.assertTrue(any(row["cv_role"] == "bond_formation" for row in rows))
            self.assertTrue(any(row["cv_role"] == "proton_transfer" for row in rows))
            self.assertTrue(all(row["source"] == "generic_serine_hydrolase_chemistry" for row in rows))

            with screening.open("r", encoding="utf-8", newline="") as handle:
                screen_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertGreaterEqual(len(screen_rows), 4)
            self.assertTrue(all(row["status"] == "not_started" for row in screen_rows))

            with ts_like.open("r", encoding="utf-8", newline="") as handle:
                ts_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(ts_rows[0]["ts_guess_id"], "pending")
            self.assertEqual(ts_rows[0]["validation_status"], "not_generated")


if __name__ == "__main__":
    unittest.main()
