import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class GenerateStage2AmberMdInputsTest(unittest.TestCase):
    def test_generates_amber_job_inputs_from_md_queue(self):
        script = Path("work/generate_stage2_amber_md_inputs.py")
        if not script.exists():
            script = Path("project01/phase2_blind_petase_qmmm_20260630/scripts/generate_stage2_amber_md_inputs.py")
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            queue = tmpdir / "md_replicate_queue.tsv"
            structure = tmpdir / "seed.pdb"
            structure.write_text("ATOM      1  N   ALA A   1       0.000   0.000   0.000  1.00 20.00           N\nEND\n", encoding="utf-8")
            queue.write_text(
                "\n".join(
                    [
                        "md_job_id\tsource_pose_id\tstage\ttemplate_pdb\tsubstrate_model\tstructure_path\treplicate_index\tensemble\ttemperature_K\tpressure_bar\tplanned_equilibration\tplanned_production\tcluster_input\tstatus\tsource",
                        f"MD_REACTIVE_6EQE_BHET_like_E01_001_R01\tREACTIVE_6EQE_BHET_like_E01_001\tpreacylation_michaelis_complex\t6EQE\tBHET_like\t{structure}\t1\tNPT_after_NVT_equilibration\t300\t1\tminimize_heat_density_equilibrate\tshort_replicate_for_active_site_clustering\tmd_trajectory_after_run\tnot_started\tblind_accepted_gs_pose",
                        "pending\tpending\tpending\tpending\tpending\tpending\tpending\tpending\tpending\tpending\tpending\tpending\tnot_ready_no_accepted_gs_pose\tnot_ready\tblind_workflow_pending",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            out_dir = tmpdir / "amber_jobs"

            result = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--md-queue",
                    str(queue),
                    "--out-dir",
                    str(out_dir),
                    "--max-jobs",
                    "1",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            manifest = out_dir / "amber_md_job_manifest.tsv"
            self.assertTrue(manifest.exists())
            with manifest.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["md_job_id"], "MD_REACTIVE_6EQE_BHET_like_E01_001_R01")
            self.assertEqual(rows[0]["stage"], "preacylation_michaelis_complex")
            self.assertEqual(rows[0]["status"], "inputs_ready_needs_topology_and_md_engine")

            job_dir = Path(rows[0]["job_dir"])
            expected_files = [
                "00_README.md",
                "01_minimize.in",
                "02_heat.in",
                "03_density.in",
                "04_equilibrate.in",
                "05_production.in",
                "run_amber_md.sh",
            ]
            for name in expected_files:
                path = job_dir / name
                self.assertTrue(path.exists(), f"missing {path}")
                lowered = path.read_text(encoding="utf-8").lower()
                self.assertNotIn("paper trajectory", lowered)
                self.assertNotIn("paper reaction coordinate", lowered)
                self.assertNotIn("paper barrier", lowered)
            production = (job_dir / "05_production.in").read_text(encoding="utf-8")
            self.assertIn("nstlim=250000", production)
            runner = (job_dir / "run_amber_md.sh").read_text(encoding="utf-8")
            self.assertIn("sander", runner)
            self.assertIn("pmemd.cuda", runner)


if __name__ == "__main__":
    unittest.main()
