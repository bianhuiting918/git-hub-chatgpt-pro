import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class GenerateStage2ClassicalMdManifestsTest(unittest.TestCase):
    def test_generates_md_queue_only_from_accepted_gs_poses(self):
        script = Path("work/generate_stage2_classical_md_manifests.py")
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            gs_manifest = tmpdir / "gs_pose_manifest.tsv"
            gs_manifest.write_text(
                "\n".join(
                    [
                        "pose_id\ttemplate_pdb\tsubstrate_model\tgeneration_method\trelaxed_structure_path\tser_og_to_c_A\tattack_angle_deg\toxyanion_hbond_count\this_acceptor_distance_A\tleaving_group_relay_distance_A\ttrp_rotamer_label\tpass_fail\trejection_reason\tnext_step",
                        "AC_GS_0001\t6EQE\tPET_dimer_capped\tdocking_md\twork/ac_gs_0001.pdb\t2.9\t102\t2\t2.0\t2.7\trotamer_a\tpass\t\tclassical_md",
                        "AC_GS_0002\t6EQE\tBHET_like\tpending\tpending\t\t\t\t\t\t\tpending\tpending\tpending",
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
                    "--gs-pose-manifest",
                    str(gs_manifest),
                    "--out-root",
                    str(out_root),
                    "--replicates",
                    "3",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)

            queue = out_root / "02_classical_md" / "md_replicate_queue.tsv"
            productive = out_root / "02_classical_md" / "productive_conformer_manifest.tsv"
            rejected = out_root / "02_classical_md" / "rejected_pose_manifest.tsv"
            protocol = out_root / "02_classical_md" / "stage2_classical_md_protocol.md"

            for path in (queue, productive, rejected, protocol):
                self.assertTrue(path.exists(), f"missing {path}")
                lowered = path.read_text(encoding="utf-8").lower()
                self.assertNotIn("paper trajectory", lowered)
                self.assertNotIn("paper ts", lowered)
                self.assertNotIn("paper barrier", lowered)

            with queue.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))

            self.assertEqual(len(rows), 3)
            self.assertEqual({row["source_pose_id"] for row in rows}, {"AC_GS_0001"})
            self.assertEqual(rows[0]["md_job_id"], "MD_AC_GS_0001_R01")
            self.assertEqual(rows[0]["structure_path"], "work/ac_gs_0001.pdb")
            self.assertEqual(rows[0]["stage"], "acylation")
            self.assertEqual(rows[0]["status"], "not_started")

            with productive.open("r", encoding="utf-8", newline="") as handle:
                productive_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(productive_rows[0]["source"], "awaiting_md_replicates")
            self.assertEqual(productive_rows[0]["selection_status"], "not_ready")


if __name__ == "__main__":
    unittest.main()