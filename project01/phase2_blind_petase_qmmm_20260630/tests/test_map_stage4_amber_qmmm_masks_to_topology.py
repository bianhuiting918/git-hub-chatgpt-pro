import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


def pdb_line(record, serial, name, resname, chain, resseq, x, y, z, element):
    return f"{record:<6}{serial:5d} {name:<4} {resname:>3} {chain}{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}"


class MapStage4AmberQmmmMasksTest(unittest.TestCase):
    def test_maps_original_seed_qm_atoms_to_tleap_atom_indices_by_occurrence(self):
        script = Path("project01/phase2_blind_petase_qmmm_20260630/scripts/map_stage4_amber_qmmm_masks_to_topology.py")
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            qmmm_dir = tmpdir / "qmmm_job"
            prep_dir = tmpdir / "topology_prep"
            out_dir = tmpdir / "mapped"
            qmmm_dir.mkdir()
            prep_dir.mkdir()

            source_pdb = tmpdir / "source_seed.pdb"
            source_pdb.write_text(
                "\n".join(
                    [
                        pdb_line("ATOM", 5, "CB", "SER", "A", 160, 1.0, 2.0, 3.0, "C"),
                        pdb_line("ATOM", 6, "OG", "SER", "A", 160, 1.5, 2.0, 3.0, "O"),
                        pdb_line("HETATM", 7, "C001", "LIG", " ", 1, 2.0, 2.0, 3.0, "C"),
                        pdb_line("ATOM", 8, "NE2", "HIS", "A", 237, 3.0, 2.0, 3.0, "N"),
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            selection = qmmm_dir / "qm_atom_selection.tsv"
            selection.write_text(
                "\n".join(
                    [
                        "serial\tatom_name\tresname\tresseq\tchain\trecord\trole",
                        "5\tCB\tSER\t160\tA\tATOM\tcatalytic_or_neighbor_sidechain",
                        "6\tOG\tSER\t160\tA\tATOM\tcatalytic_or_neighbor_sidechain",
                        "7\tC001\tLIG\t1\t\tHETATM\tsubstrate_ligand",
                        "8\tNE2\tHIS\t237\tA\tATOM\tcatalytic_or_neighbor_sidechain",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            min_mdin = qmmm_dir / "01_qmmm_min.in"
            equil_mdin = qmmm_dir / "02_qmmm_equil_200ps.in"
            min_mdin.write_text("&cntrl\n  imin=1, maxcyc=2000, ifqnt=1,\n/\n&qmmm\n  qmmask='@5,6,7,8',\n  qm_theory='DFTB3',\n/\n", encoding="utf-8")
            equil_mdin.write_text("&cntrl\n  imin=0, nstlim=200000, ifqnt=1,\n/\n&qmmm\n  qmmask='@5,6,7,8',\n  qm_theory='DFTB3',\n/\n", encoding="utf-8")
            qmmm_manifest = tmpdir / "amber_qmmm_job_manifest.tsv"
            qmmm_manifest.write_text(
                "\n".join(
                    [
                        "qmmm_job_id\tpose_id\tligand_model\tsource_structure_path\tengine\tqm_theory\tqmcharge\tspin\tqm_atom_count\tqmmask\tjob_dir\tqm_selection_path\tmin_mdin_path\tequil_mdin_path\trunner_path\tstatus\tsource",
                        f"AMBER_QMMM_TEST_001\tPOSE1\tBHET_like\t{source_pdb}\tamber_sander\tDFTB3\t0\t1\t4\t@5,6,7,8\t{qmmm_dir}\t{selection}\t{min_mdin}\t{equil_mdin}\t{qmmm_dir / 'run.sh'}\tinputs_ready_needs_amber_prmtop_inpcrd\tblind_stage1_accepted_seed",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            complex_leap = prep_dir / "complex_leap.pdb"
            complex_leap.write_text(
                "\n".join(
                    [
                        pdb_line("ATOM", 20, "CB", "SER", " ", 12, 11.0, 0.0, 8.0, "C"),
                        pdb_line("ATOM", 21, "OG", "SER", " ", 12, 11.5, 0.0, 8.0, "O"),
                        pdb_line("HETATM", 22, "C001", "LIG", " ", 98, 12.0, 0.0, 8.0, "C"),
                        pdb_line("ATOM", 23, "NE2", "HIE", " ", 209, 13.0, 0.0, 8.0, "N"),
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            (prep_dir / "complex.prmtop").write_text("dummy prmtop\n", encoding="utf-8")
            (prep_dir / "complex.inpcrd").write_text("dummy inpcrd\n", encoding="utf-8")
            topology_manifest = tmpdir / "amber_topology_prep_manifest.tsv"
            topology_manifest.write_text(
                "\n".join(
                    [
                        "topology_prep_id\tqmmm_job_id\tpose_id\tligand_model\tsource_structure_path\tprep_dir\tcomplex_pdb_path\tligand_pdb_path\ttleap_input_path\trunner_path\tprmtop_path\tinpcrd_path\tstatus\tsource",
                        f"TOPOLOGY_AMBER_QMMM_TEST_001\tAMBER_QMMM_TEST_001\tPOSE1\tBHET_like\t{source_pdb}\t{prep_dir}\t{prep_dir / 'complex_for_leap.pdb'}\t{prep_dir / 'ligand.pdb'}\t{prep_dir / 'tleap.in'}\t{prep_dir / 'run_amber_topology_prep.sh'}\t{prep_dir / 'complex.prmtop'}\t{prep_dir / 'complex.inpcrd'}\ttopology_prep_ready_requires_ambertools_execution\tblind_stage4_amber_qmmm_manifest",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--amber-qmmm-manifest",
                    str(qmmm_manifest),
                    "--topology-prep-manifest",
                    str(topology_manifest),
                    "--out-dir",
                    str(out_dir),
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            mapped_manifest = out_dir / "amber_qmmm_topology_mapped_manifest.tsv"
            self.assertTrue(mapped_manifest.exists())
            with mapped_manifest.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["mapped_qmmask"], "@20,21,22,23")
            self.assertEqual(rows[0]["status"], "qmmm_inputs_mapped_to_amber_topology")

            job_dir = Path(rows[0]["mapped_job_dir"])
            smoke = (job_dir / "00_qmmm_smoke_min_1step.in").read_text(encoding="utf-8")
            self.assertIn("qmmask='@20,21,22,23'", smoke)
            self.assertIn("maxcyc=1", smoke)
            self.assertIn("ifqnt=1", smoke)
            runner = (job_dir / "run_sander_qmmm_smoke.sh").read_text(encoding="utf-8")
            self.assertIn("complex.prmtop", runner)
            self.assertIn("complex.inpcrd", runner)
            self.assertIn("sander", runner)
            self.assertIn("AMBERHOME", runner)
            self.assertIn("command -v", runner)


if __name__ == "__main__":
    unittest.main()

