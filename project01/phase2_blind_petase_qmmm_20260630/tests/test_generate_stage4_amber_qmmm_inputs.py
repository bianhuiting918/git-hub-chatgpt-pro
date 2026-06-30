import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class GenerateStage4AmberQmmmInputsTest(unittest.TestCase):
    def test_generates_sander_dftb3_qmmm_inputs_from_blind_seed_structures(self):
        script = Path("project01/phase2_blind_petase_qmmm_20260630/scripts/generate_stage4_amber_qmmm_inputs.py")
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            repo_root = tmpdir / "repo_root"
            relative_structure = "project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/reactive_pose_seeds/6EQE_BHET_like_E01/seed.pdb"
            manifest = repo_root / "project01" / "phase2_blind_petase_qmmm_20260630" / "blind_work" / "01_system_setup" / "gs_pose_manifest.tsv"
            structure = repo_root / relative_structure
            manifest.parent.mkdir(parents=True, exist_ok=True)
            structure.parent.mkdir(parents=True, exist_ok=True)
            structure.write_text(
                "\n".join(
                    [
                        "ATOM      1  N   TRP A 159      11.000  10.000  10.000  1.00 20.00           N",
                        "ATOM      2  CA  TRP A 159      12.000  10.000  10.000  1.00 20.00           C",
                        "ATOM      3  CB  TRP A 159      13.000  10.000  10.000  1.00 20.00           C",
                        "ATOM      4  CG  TRP A 159      14.000  10.000  10.000  1.00 20.00           C",
                        "ATOM      5  N   SER A 160      11.000  11.000  10.000  1.00 20.00           N",
                        "ATOM      6  CA  SER A 160      12.000  11.000  10.000  1.00 20.00           C",
                        "ATOM      7  CB  SER A 160      13.000  11.000  10.000  1.00 20.00           C",
                        "ATOM      8  OG  SER A 160      14.000  11.000  10.000  1.00 20.00           O",
                        "ATOM      9  N   ASP A 206      11.000  12.000  10.000  1.00 20.00           N",
                        "ATOM     10  CA  ASP A 206      12.000  12.000  10.000  1.00 20.00           C",
                        "ATOM     11  CB  ASP A 206      13.000  12.000  10.000  1.00 20.00           C",
                        "ATOM     12  OD1 ASP A 206      14.000  12.000  10.000  1.00 20.00           O",
                        "ATOM     13  OD2 ASP A 206      14.000  12.500  10.000  1.00 20.00           O",
                        "ATOM     14  N   HIS A 237      11.000  13.000  10.000  1.00 20.00           N",
                        "ATOM     15  CA  HIS A 237      12.000  13.000  10.000  1.00 20.00           C",
                        "ATOM     16  CB  HIS A 237      13.000  13.000  10.000  1.00 20.00           C",
                        "ATOM     17  ND1 HIS A 237      14.000  13.000  10.000  1.00 20.00           N",
                        "ATOM     18  NE2 HIS A 237      14.000  13.500  10.000  1.00 20.00           N",
                        "ATOM     19  N   SER A 238      11.000  14.000  10.000  1.00 20.00           N",
                        "ATOM     20  CA  SER A 238      12.000  14.000  10.000  1.00 20.00           C",
                        "ATOM     21  CB  SER A 238      13.000  14.000  10.000  1.00 20.00           C",
                        "ATOM     22  OG  SER A 238      14.000  14.000  10.000  1.00 20.00           O",
                        "HETATM   23 C001 LIG     1      15.000  11.000  10.000  1.00 20.00           C",
                        "HETATM   24 O002 LIG     1      16.000  11.000  10.000  1.00 20.00           O",
                        "HETATM   25 O003 LIG     1      15.000  12.000  10.000  1.00 20.00           O",
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            manifest.write_text(
                "\n".join(
                    [
                        "pose_id\tpass_fail\tligand_model\trelaxed_structure_path\tgeneration_method",
                        f"REACTIVE_6EQE_BHET_like_E01_001\tpass\tBHET_like\t{relative_structure}\treactive_geometry_seed",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            out_dir = tmpdir / "amber_qmmm_inputs"

            result = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--gs-pose-manifest",
                    str(manifest),
                    "--out-dir",
                    str(out_dir),
                    "--max-jobs",
                    "1",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            job_manifest = out_dir / "amber_qmmm_job_manifest.tsv"
            self.assertTrue(job_manifest.exists())
            with job_manifest.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["qmmm_job_id"], "AMBER_QMMM_REACTIVE_6EQE_BHET_like_E01_001_001")
            self.assertEqual(rows[0]["qm_theory"], "DFTB3")
            self.assertEqual(rows[0]["engine"], "amber_sander")
            self.assertEqual(rows[0]["status"], "inputs_ready_needs_amber_prmtop_inpcrd")

            job_dir = Path(rows[0]["job_dir"])
            for name in (
                "00_README.md",
                "qm_atom_selection.tsv",
                "01_qmmm_min.in",
                "02_qmmm_equil_200ps.in",
                "run_sander_qmmm_smoke.sh",
            ):
                path = job_dir / name
                self.assertTrue(path.exists(), f"missing {path}")
                text = path.read_text(encoding="utf-8")
                lowered = text.lower()
                self.assertNotIn("paper trajectory", lowered)
                self.assertNotIn("paper reaction coordinate", lowered)
                self.assertNotIn("paper barrier", lowered)

            selection = (job_dir / "qm_atom_selection.tsv").read_text(encoding="utf-8")
            self.assertIn("SER\t160", selection)
            self.assertIn("HIS\t237", selection)
            self.assertIn("ASP\t206", selection)
            self.assertIn("LIG\t1", selection)

            min_mdin = (job_dir / "01_qmmm_min.in").read_text(encoding="utf-8")
            self.assertIn("ifqnt=1", min_mdin)
            self.assertIn("&qmmm", min_mdin)
            self.assertIn("qm_theory='DFTB3'", min_mdin)
            self.assertIn("qmmask='@", min_mdin)
            self.assertIn("maxcyc=2000", min_mdin)

            equil_mdin = (job_dir / "02_qmmm_equil_200ps.in").read_text(encoding="utf-8")
            self.assertIn("nstlim=200000", equil_mdin)
            self.assertIn("temp0=310.0", equil_mdin)

            runner = (job_dir / "run_sander_qmmm_smoke.sh").read_text(encoding="utf-8")
            self.assertIn("sander", runner)
            self.assertIn("PRMTOP", runner)
            self.assertIn("INPCRD", runner)


if __name__ == "__main__":
    unittest.main()



