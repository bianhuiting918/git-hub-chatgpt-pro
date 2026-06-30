import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class GenerateStage1ReactivePoseSeedsTest(unittest.TestCase):
    def test_generates_complex_that_passes_generic_geometry_scorer(self):
        script = Path("work/generate_stage1_reactive_pose_seeds.py")
        scorer = Path("work/score_stage1_pose_geometry.py")
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            protein = tmpdir / "protein.pdb"
            ligand = tmpdir / "ligand.pdb"
            labels = tmpdir / "labels.tsv"
            out_dir = tmpdir / "reactive"
            score_tsv = tmpdir / "score.tsv"

            protein.write_text(
                "\n".join(
                    [
                        "ATOM      1  OG  SER A 160       0.000   0.000   0.000  1.00 20.00           O",
                        "ATOM      2  NE2 HIS A 237       2.900   0.000   0.000  1.00 20.00           N",
                        "ATOM      3  ND1 HIS A 237       2.900   1.000   0.000  1.00 20.00           N",
                        "ATOM      4  N   GLY A 161       2.250   2.050   0.000  1.00 20.00           N",
                        "ATOM      5  N   TRP A 185       0.000   4.000   0.000  1.00 20.00           N",
                        "ATOM      6  CA  TRP A 185       0.000   5.000   0.000  1.00 20.00           C",
                        "ATOM      7  CB  TRP A 185       1.000   5.000   0.000  1.00 20.00           C",
                        "ATOM      8  CG  TRP A 185       1.000   6.000   0.000  1.00 20.00           C",
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            ligand.write_text(
                "\n".join(
                    [
                        "HETATM    1 C001 LIG A   1       0.000   0.000   0.000  1.00  0.00           C",
                        "HETATM    2 O002 LIG A   1       1.230   0.000   0.000  1.00  0.00           O",
                        "HETATM    3 O003 LIG A   1       0.000   1.350   0.000  1.00  0.00           O",
                        "HETATM    4 C004 LIG A   1      -1.000   0.000   0.000  1.00  0.00           C",
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            labels.write_text(
                "model_id\tcandidate_id\tscissile_carbonyl_c_atom_name\tscissile_carbonyl_o_atom_name\tleaving_o_atom_name\n"
                "BHET_like\tE01\tC001\tO002\tO003\n",
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--protein-pdb",
                    str(protein),
                    "--ligand-pdb",
                    str(ligand),
                    "--label-tsv",
                    str(labels),
                    "--model-id",
                    "BHET_like",
                    "--candidate-id",
                    "E01",
                    "--template-pdb",
                    "TEST",
                    "--triad-ser",
                    "A:SER160",
                    "--triad-his",
                    "A:HIS237",
                    "--out-dir",
                    str(out_dir),
                    "--max-seeds",
                    "1",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            manifest = out_dir / "reactive_pose_seed_manifest.tsv"
            self.assertTrue(manifest.exists())
            with manifest.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["status"], "ready_for_geometry_scoring")
            self.assertEqual(rows[0]["source"], "blind_reactive_geometry_seed")

            score = subprocess.run(
                [
                    sys.executable,
                    str(scorer),
                    "--complex-pdb",
                    rows[0]["complex_pdb"],
                    "--label-tsv",
                    str(labels),
                    "--model-id",
                    "BHET_like",
                    "--candidate-id",
                    "E01",
                    "--pose-id",
                    rows[0]["pose_id"],
                    "--template-pdb",
                    "TEST",
                    "--generation-method",
                    "reactive_geometry_seed",
                    "--ser-og",
                    "A:160:OG",
                    "--his",
                    "A:237",
                    "--trp",
                    "A:185",
                    "--out-tsv",
                    str(score_tsv),
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(score.returncode, 0, score.stderr)
            with score_tsv.open("r", encoding="utf-8", newline="") as handle:
                score_rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(score_rows[0]["pass_fail"], "pass")


if __name__ == "__main__":
    unittest.main()
