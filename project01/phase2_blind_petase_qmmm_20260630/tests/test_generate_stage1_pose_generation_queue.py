import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class GenerateStage1PoseGenerationQueueTest(unittest.TestCase):
    def test_generates_docking_queue_from_structure_triad_and_ligand_manifest(self):
        script = Path("work/generate_stage1_pose_generation_queue.py")
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            prepared_dir = tmpdir / "prepared"
            ligand_dir = tmpdir / "ligands" / "PET_dimer_capped"
            label_dir = tmpdir / "labels"
            prepared_dir.mkdir(parents=True)
            ligand_dir.mkdir(parents=True)
            label_dir.mkdir(parents=True)

            protein = prepared_dir / "6EQE_chainA_initial_clean_v2.pdb"
            protein.write_text(
                "\n".join(
                    [
                        "ATOM      1  OG  SER A 160      10.000  20.000  30.000  1.00 20.00           O",
                        "ATOM      2  NE2 HIS A 237      12.000  20.000  30.000  1.00 20.00           N",
                        "ATOM      3  OD2 ASP A 206      14.000  20.000  30.000  1.00 20.00           O",
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            ligand = ligand_dir / "PET_dimer_capped_conf001_named_atoms.pdb"
            ligand.write_text("HETATM    1 C001 LIG A   1       0.000   0.000   0.000  1.00  0.00           C\nEND\n", encoding="utf-8")
            labels = label_dir / "PET_dimer_capped_atoms.tsv"
            labels.write_text(
                "model_id\tcandidate_id\tscissile_carbonyl_c_atom_name\tscissile_carbonyl_o_atom_name\tleaving_o_atom_name\tselection_status\n"
                "PET_dimer_capped\tE01\tC001\tO002\tO003\tcandidate_not_selected\n",
                encoding="utf-8",
            )

            prepared_manifest = tmpdir / "prepared_structure_manifest.tsv"
            prepared_manifest.write_text(
                "pdb\tprepared_path\tselected_chain\tprotonation_status\tsha256\tstatus\n"
                f"6EQE\t{protein}\tA\treviewed_ph7\tabc\tprepared\n",
                encoding="utf-8",
            )
            triads = tmpdir / "ser_his_asp_triad_candidates.tsv"
            triads.write_text(
                "pdb\tser\this\tasp\tser_his_min_A\tser_his_atoms\this_asp_min_A\this_asp_atoms\n"
                "6EQE\tA:SER160\tA:HIS237\tA:ASP206\t2.94\tOG-NE2\t2.66\tND1-OD2\n"
                "6EQE\tB:SER160\tB:HIS237\tB:ASP206\t2.94\tOG-NE2\t2.66\tND1-OD2\n",
                encoding="utf-8",
            )
            ligand_manifest = tmpdir / "ligand_build_manifest.tsv"
            ligand_manifest.write_text(
                "model_id\tstatus\tformal_charge\tconformer_count\tsdf_path\tpdb_path\tatom_label_path\tsha256_sdf\tsha256_pdb\tnote\n"
                f"PET_dimer_capped\tbuilt_needs_scissile_candidate_selection\t0\t1\t\t{ligand}\t{labels}\t\tligsha\tforce_field=UFF\n",
                encoding="utf-8",
            )
            out_root = tmpdir / "blind_work"

            result = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--prepared-manifest",
                    str(prepared_manifest),
                    "--triad-manifest",
                    str(triads),
                    "--ligand-build-manifest",
                    str(ligand_manifest),
                    "--out-root",
                    str(out_root),
                    "--box-size",
                    "22",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            queue = out_root / "01_system_setup" / "pose_generation_queue.tsv"
            config = out_root / "01_system_setup" / "docking_inputs" / "POSE_6EQE_PET_dimer_capped_E01.config"
            self.assertTrue(queue.exists())
            self.assertTrue(config.exists())

            with queue.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(rows[0]["pose_job_id"], "POSE_6EQE_PET_dimer_capped_E01")
            self.assertEqual(rows[0]["box_center_x"], "10.000")
            self.assertEqual(rows[0]["box_center_y"], "20.000")
            self.assertEqual(rows[0]["box_center_z"], "30.000")
            self.assertEqual(rows[0]["triad_ser"], "A:SER160")
            self.assertEqual(rows[0]["status"], "ready_for_docking")
            self.assertEqual(rows[0]["source"], "blind_structure_geometry")

            config_text = config.read_text(encoding="utf-8").lower()
            self.assertIn("center_x = 10.000", config_text)
            self.assertIn("size_x = 22.000", config_text)
            self.assertNotIn("paper", config_text)
            self.assertNotIn("article", config_text)


if __name__ == "__main__":
    unittest.main()
