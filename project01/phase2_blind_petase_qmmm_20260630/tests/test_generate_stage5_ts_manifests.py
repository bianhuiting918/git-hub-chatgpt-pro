import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class GenerateStage5TsManifestsTest(unittest.TestCase):
    def test_generates_ts_refinement_and_ensemble_manifests_from_verified_guesses(self):
        script = Path("work/generate_stage5_ts_manifests.py")
        if not script.exists():
            script = Path("project01/phase2_blind_petase_qmmm_20260630/scripts/generate_stage5_ts_manifests.py")
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            guesses = tmpdir / "ts_like_guess_manifest.tsv"
            guesses.write_text(
                "\n".join(
                    [
                        "ts_guess_id\tpath_id\tsource_pose_id\tsource_method\tstructure_path\timaginary_mode_status\tcommittor_status\tvalidation_status\tnotes",
                        "G_AC_001\tAC_DIRECT_HIS_GENERAL_BASE\tAC_GS_0001\tscan_maximum\twork/guess_ac_001.pdb\tnot_checked\tnot_started\tendpoint_checked\tblind acylation guess",
                        "G_DE_001\tDE_WATER_ATTACK_HIS_GENERAL_BASE\tDE_GS_0001\tstring_image\twork/guess_de_001.pdb\tnot_checked\tnot_started\tendpoint_checked\tblind deacylation guess",
                        "pending\tpending\tpending\tscan_maximum_or_string_image_after_stage4_runs\tpending\tnot_generated\tnot_generated\tnot_generated\tplaceholder",
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
                    "--ts-like-manifest",
                    str(guesses),
                    "--out-root",
                    str(out_root),
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)

            refinement = out_root / "05_ts_refinement" / "ts_refinement_manifest.tsv"
            ac_ensemble = out_root / "06_ts_ensemble" / "acylation_ts_ensemble.tsv"
            de_ensemble = out_root / "06_ts_ensemble" / "deacylation_ts_ensemble.tsv"
            committor = out_root / "07_committor" / "committor_queue.tsv"

            for path in (refinement, ac_ensemble, de_ensemble, committor):
                self.assertTrue(path.exists(), f"missing {path}")
                lowered = path.read_text(encoding="utf-8").lower()
                self.assertNotIn("cv77", lowered)
                self.assertNotIn("paper_rc", lowered)
                self.assertNotIn("paper trajectory", lowered)

            with refinement.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual([row["ts_guess_id"] for row in rows], ["G_AC_001", "G_DE_001"])
            self.assertTrue(all(row["refinement_status"] == "not_started" for row in rows))
            self.assertTrue(all(row["source"] == "blind_stage4_ts_like_guess" for row in rows))

            with ac_ensemble.open("r", encoding="utf-8", newline="") as handle:
                ac_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(ac_rows[0]["reaction_stage"], "acylation")
            self.assertEqual(ac_rows[0]["ensemble_status"], "candidate_pending_refinement")

            with de_ensemble.open("r", encoding="utf-8", newline="") as handle:
                de_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(de_rows[0]["reaction_stage"], "deacylation")
            self.assertEqual(de_rows[0]["ensemble_status"], "candidate_pending_refinement")

            with committor.open("r", encoding="utf-8", newline="") as handle:
                queue_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual({row["committor_status"] for row in queue_rows}, {"not_ready_refinement_required"})


if __name__ == "__main__":
    unittest.main()
