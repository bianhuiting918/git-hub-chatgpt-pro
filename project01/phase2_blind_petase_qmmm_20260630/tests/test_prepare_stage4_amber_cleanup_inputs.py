import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class PrepareStage4AmberCleanupInputsTest(unittest.TestCase):
    def test_generates_staged_mm_and_qmmm_cleanup_inputs(self):
        script = Path("project01/phase2_blind_petase_qmmm_20260630/scripts/prepare_stage4_amber_cleanup_inputs.py")
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            mapped_dir = tmpdir / "mapped_job"
            mapped_dir.mkdir()
            prmtop = tmpdir / "complex.prmtop"
            inpcrd = tmpdir / "complex.inpcrd"
            prmtop.write_text("dummy prmtop\n", encoding="utf-8")
            inpcrd.write_text("dummy inpcrd\n", encoding="utf-8")
            mapped_manifest = tmpdir / "amber_qmmm_topology_mapped_manifest.tsv"
            mapped_manifest.write_text(
                "\n".join(
                    [
                        "mapped_job_id\tqmmm_job_id\ttopology_prep_id\tpose_id\tmapped_atom_count\tmapped_qmmask\tmapped_job_dir\tsmoke_mdin_path\tmin_mdin_path\tequil_mdin_path\trunner_path\tprmtop_path\tinpcrd_path\tstatus",
                        f"MAPPED_AMBER_QMMM_TEST_001\tAMBER_QMMM_TEST_001\tTOPOLOGY_TEST\tPOSE1\t4\t@20,21,22,23\t{mapped_dir}\t{mapped_dir / 'smoke.in'}\t{mapped_dir / 'min.in'}\t{mapped_dir / 'equil.in'}\t{mapped_dir / 'run.sh'}\t{prmtop}\t{inpcrd}\tqmmm_inputs_mapped_to_amber_topology",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            out_dir = tmpdir / "cleanup"

            result = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--mapped-qmmm-manifest",
                    str(mapped_manifest),
                    "--out-dir",
                    str(out_dir),
                    "--mm-solvent-maxcyc",
                    "50",
                    "--mm-all-maxcyc",
                    "25",
                    "--qmmm-short-maxcyc",
                    "2",
                    "--qmmm-full-maxcyc",
                    "100",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            cleanup_manifest = out_dir / "amber_cleanup_manifest.tsv"
            self.assertTrue(cleanup_manifest.exists())
            with cleanup_manifest.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["cleanup_job_id"], "CLEANUP_MAPPED_AMBER_QMMM_TEST_001")
            self.assertEqual(rows[0]["status"], "cleanup_inputs_ready_requires_sander_execution")
            self.assertEqual(rows[0]["mapped_qmmask"], "@20,21,22,23")

            job_dir = Path(rows[0]["cleanup_job_dir"])
            expected = [
                "00_mm_solvent_ion_min.in",
                "01_mm_all_min.in",
                "02_qmmm_short_min.in",
                "03_qmmm_min.in",
                "run_amber_cleanup.sh",
                "00_README.md",
            ]
            for name in expected:
                path = job_dir / name
                self.assertTrue(path.exists(), f"missing {path}")
                lowered = path.read_text(encoding="utf-8").lower()
                self.assertNotIn("paper trajectory", lowered)
                self.assertNotIn("paper reaction coordinate", lowered)
                self.assertNotIn("paper barrier", lowered)

            solvent = (job_dir / "00_mm_solvent_ion_min.in").read_text(encoding="utf-8")
            self.assertIn("imin=1", solvent)
            self.assertIn("maxcyc=50", solvent)
            self.assertIn("ntr=1", solvent)
            self.assertIn("restraintmask='!:WAT,Na+,Cl-'", solvent)
            self.assertNotIn("ifqnt=1", solvent)

            all_mm = (job_dir / "01_mm_all_min.in").read_text(encoding="utf-8")
            self.assertIn("maxcyc=25", all_mm)
            self.assertNotIn("restraintmask", all_mm)

            qmmm_short = (job_dir / "02_qmmm_short_min.in").read_text(encoding="utf-8")
            self.assertIn("ifqnt=1", qmmm_short)
            self.assertIn("maxcyc=2", qmmm_short)
            self.assertIn("qmmask='@20,21,22,23'", qmmm_short)
            self.assertIn("qm_theory='DFTB3'", qmmm_short)

            qmmm_full = (job_dir / "03_qmmm_min.in").read_text(encoding="utf-8")
            self.assertIn("maxcyc=100", qmmm_full)
            self.assertIn("qmmask='@20,21,22,23'", qmmm_full)

            runner = (job_dir / "run_amber_cleanup.sh").read_text(encoding="utf-8")
            self.assertIn("AMBERHOME", runner)
            self.assertIn("00_mm_solvent_ion_min.in", runner)
            self.assertIn("01_mm_all_min.in", runner)
            self.assertIn("02_qmmm_short_min.in", runner)
            self.assertIn("03_qmmm_min.in", runner)
            self.assertIn("03_qmmm_min.rst7", runner)


if __name__ == "__main__":
    unittest.main()
