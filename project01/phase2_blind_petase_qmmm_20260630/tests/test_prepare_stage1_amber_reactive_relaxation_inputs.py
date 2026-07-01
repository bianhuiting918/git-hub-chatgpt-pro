import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class PrepareStage1AmberReactiveRelaxationInputsTest(unittest.TestCase):
    def test_generates_restrained_mm_relaxation_inputs_from_no_clash_docking_pose(self):
        script = Path("project01/phase2_blind_petase_qmmm_20260630/scripts/prepare_stage1_amber_reactive_relaxation_inputs.py")
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            complex_pdb = tmpdir / "complex_rank001.pdb"
            complex_pdb.write_text(
                "\n".join(
                    [
                        "ATOM      1  OG  SER A 160       0.000   0.000   0.000  1.00 20.00           O",
                        "ATOM      2  N   MET A 161       2.600   0.800   0.000  1.00 20.00           N",
                        "ATOM      3  NE2 HIS A 237       3.800   0.000   0.000  1.00 20.00           N",
                        "HETATM    4 C005 LIG     1       3.500   0.000   0.000  1.00  0.00           C",
                        "HETATM    5 O006 LIG     1       4.700   0.000   0.000  1.00  0.00           O",
                        "HETATM    6 O004 LIG     1       3.500   1.300   0.000  1.00  0.00           O",
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            score_summary = tmpdir / "score_and_clash_summary.tsv"
            score_summary.write_text(
                "\n".join(
                    [
                        "pose_id\ttemplate_pdb\tsubstrate_model\tgeneration_method\tcomplex_pdb\tser_og_to_c_A\tattack_angle_deg\tmin_ligand_protein_heavy_contact_A",
                        f"DOCK_RANK001\t6EQE\tBHET_like\tvina_box12\t{complex_pdb}\t3.500\t50.0\t2.500",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            labels = tmpdir / "labels.tsv"
            labels.write_text(
                "model_id\tcandidate_id\tscissile_carbonyl_c_atom_name\tscissile_carbonyl_o_atom_name\tleaving_o_atom_name\n"
                "BHET_like\tE01\tC005\tO006\tO004\n",
                encoding="utf-8",
            )
            out_dir = tmpdir / "relax"

            result = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--pose-score-summary",
                    str(score_summary),
                    "--label-tsv",
                    str(labels),
                    "--model-id",
                    "BHET_like",
                    "--candidate-id",
                    "E01",
                    "--ser-og",
                    "A:160:OG",
                    "--his-acceptor",
                    "A:237:NE2",
                    "--oxyanion-donor",
                    "A:161:N",
                    "--out-dir",
                    str(out_dir),
                    "--max-poses",
                    "1",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            manifest = out_dir / "amber_reactive_relaxation_manifest.tsv"
            self.assertTrue(manifest.exists())
            with manifest.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["pose_id"], "DOCK_RANK001")
            self.assertEqual(rows[0]["status"], "reactive_relaxation_inputs_ready_requires_amber_topology")

            job_dir = Path(rows[0]["relaxation_job_dir"])
            complex_for_amber = Path(rows[0]["complex_for_amber_pdb_path"])
            self.assertTrue(complex_for_amber.exists())
            self.assertEqual(
                [line[6:11] for line in complex_for_amber.read_text(encoding="utf-8").splitlines() if line.startswith(("ATOM  ", "HETATM"))],
                ["    1", "    2", "    3", "    4", "    5", "    6"],
            )

            restraints = (job_dir / "reactive_relaxation_restraints.RST").read_text(encoding="utf-8")
            self.assertIn("iat=1,4", restraints)
            self.assertIn("iat=1,4,5", restraints)
            self.assertIn("iat=6,3", restraints)
            self.assertIn("iat=5,2", restraints)

            mdin = (job_dir / "00_restrained_mm_min.in").read_text(encoding="utf-8")
            self.assertIn("imin=1", mdin)
            self.assertIn("nmropt=1", mdin)
            self.assertIn("&wt", mdin)
            self.assertIn("type='END'", mdin)
            self.assertIn("DISANG=reactive_relaxation_restraints.RST", mdin)
            self.assertNotIn("ifqnt=1", mdin)

            runner = (job_dir / "run_reactive_mm_relaxation.sh").read_text(encoding="utf-8")
            self.assertIn("sander", runner)
            self.assertIn("PRMTOP", runner)
            self.assertIn("INPCRD", runner)
            self.assertIn("00_restrained_mm_min.rst7", runner)

            readme = (job_dir / "00_README.md").read_text(encoding="utf-8")
            self.assertIn("complex_for_amber.pdb", readme)


if __name__ == "__main__":
    unittest.main()
