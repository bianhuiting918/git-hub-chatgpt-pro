import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class PrepareStage4AmberTopologyInputsTest(unittest.TestCase):
    def test_generates_ambertools_topology_prep_from_qmmm_manifest(self):
        script = Path("project01/phase2_blind_petase_qmmm_20260630/scripts/prepare_stage4_amber_topology_inputs.py")
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            source_pdb = tmpdir / "seed.pdb"
            source_pdb.write_text(
                "\n".join(
                    [
                        "ATOM      1  N   SER A 160      11.000  10.000  10.000  1.00 20.00           N",
                        "ATOM      2  CA  SER A 160      12.000  10.000  10.000  1.00 20.00           C",
                        "ATOM      3  CB  SER A 160      13.000  10.000  10.000  1.00 20.00           C",
                        "ATOM      4  OG  SER A 160      14.000  10.000  10.000  1.00 20.00           O",
                        "HETATM    5 C001 LIG     1      15.000  11.000  10.000  1.00 20.00           C",
                        "HETATM    6 O002 LIG     1      16.000  11.000  10.000  1.00 20.00           O",
                        "HETATM    7 O003 LIG     1      15.000  12.000  10.000  1.00 20.00           O",
                        "TER",
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            qmmm_manifest = tmpdir / "amber_qmmm_job_manifest.tsv"
            qmmm_manifest.write_text(
                "\n".join(
                    [
                        "qmmm_job_id\tpose_id\tligand_model\tsource_structure_path\tengine\tqm_theory\tqmcharge\tspin\tqm_atom_count\tqmmask\tjob_dir\tqm_selection_path\tmin_mdin_path\tequil_mdin_path\trunner_path\tstatus\tsource",
                        f"AMBER_QMMM_REACTIVE_001_001\tREACTIVE_001\tBHET_like\t{source_pdb}\tamber_sander\tDFTB3\t0\t1\t3\t@5,6,7\t{tmpdir / 'qmmm_job'}\tselection.tsv\tmin.in\tequil.in\trun.sh\tinputs_ready_needs_amber_prmtop_inpcrd\tblind_stage1_accepted_seed",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            out_dir = tmpdir / "topology_prep"

            result = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--amber-qmmm-manifest",
                    str(qmmm_manifest),
                    "--out-dir",
                    str(out_dir),
                    "--max-jobs",
                    "1",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            manifest = out_dir / "amber_topology_prep_manifest.tsv"
            self.assertTrue(manifest.exists())
            with manifest.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["topology_prep_id"], "TOPOLOGY_AMBER_QMMM_REACTIVE_001_001")
            self.assertEqual(rows[0]["status"], "topology_prep_ready_requires_ambertools_execution")
            self.assertEqual(rows[0]["prmtop_path"], str(Path(rows[0]["prep_dir"]) / "complex.prmtop"))
            self.assertEqual(rows[0]["inpcrd_path"], str(Path(rows[0]["prep_dir"]) / "complex.inpcrd"))

            prep_dir = Path(rows[0]["prep_dir"])
            for name in ("complex_for_leap.pdb", "ligand.pdb", "tleap.in", "run_amber_topology_prep.sh", "00_README.md"):
                path = prep_dir / name
                self.assertTrue(path.exists(), f"missing {path}")
                lowered = path.read_text(encoding="utf-8").lower()
                self.assertNotIn("paper trajectory", lowered)
                self.assertNotIn("paper reaction coordinate", lowered)
                self.assertNotIn("paper barrier", lowered)

            ligand = (prep_dir / "ligand.pdb").read_text(encoding="utf-8")
            self.assertIn("HETATM", ligand)
            self.assertIn(" LIG ", ligand)
            self.assertNotIn("ATOM", ligand)

            tleap = (prep_dir / "tleap.in").read_text(encoding="utf-8")
            self.assertIn("source leaprc.protein.ff14SB", tleap)
            self.assertIn("source leaprc.gaff2", tleap)
            self.assertIn("source leaprc.water.tip3p", tleap)
            self.assertIn("LIG = loadmol2 ligand.mol2", tleap)
            self.assertIn("solvatebox complex TIP3PBOX 15.0", tleap)
            self.assertIn("saveamberparm complex complex.prmtop complex.inpcrd", tleap)

            runner = (prep_dir / "run_amber_topology_prep.sh").read_text(encoding="utf-8")
            self.assertIn("antechamber", runner)
            self.assertIn("-rn LIG", runner)
            self.assertIn("-c bcc", runner)
            self.assertIn("parmchk2", runner)
            self.assertIn("tleap", runner)
            self.assertIn("test -s complex.prmtop", runner)
            self.assertIn("test -s complex.inpcrd", runner)


if __name__ == "__main__":
    unittest.main()

