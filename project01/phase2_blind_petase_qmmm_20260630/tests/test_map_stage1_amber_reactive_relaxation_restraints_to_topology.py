import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


def pdb_line(record, serial, name, resname, chain, resseq, x, y, z, element):
    return f"{record:<6}{serial:5d} {name:<4} {resname:>3} {chain}{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element:>2}"


class MapStage1AmberReactiveRelaxationRestraintsTest(unittest.TestCase):
    def test_maps_reactive_relaxation_rst_iats_to_tleap_atom_indices(self):
        script = Path(
            "project01/phase2_blind_petase_qmmm_20260630/scripts/"
            "map_stage1_amber_reactive_relaxation_restraints_to_topology.py"
        )
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            relax_dir = tmpdir / "relax_job"
            prep_dir = tmpdir / "topology_prep"
            out_dir = tmpdir / "mapped"
            relax_dir.mkdir()
            prep_dir.mkdir()

            complex_for_amber = relax_dir / "complex_for_amber.pdb"
            complex_for_amber.write_text(
                "\n".join(
                    [
                        pdb_line("ATOM", 1, "OG", "SER", "A", 160, 0.0, 0.0, 0.0, "O"),
                        pdb_line("ATOM", 2, "N", "MET", "A", 161, 2.6, 0.8, 0.0, "N"),
                        pdb_line("ATOM", 3, "NE2", "HIS", "A", 237, 3.8, 0.0, 0.0, "N"),
                        pdb_line("HETATM", 4, "C005", "LIG", " ", 1, 3.5, 0.0, 0.0, "C"),
                        pdb_line("HETATM", 5, "O006", "LIG", " ", 1, 4.7, 0.0, 0.0, "O"),
                        pdb_line("HETATM", 6, "O004", "LIG", " ", 1, 3.5, 1.3, 0.0, "O"),
                        "END",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            source_rst = relax_dir / "reactive_relaxation_restraints.RST"
            source_rst.write_text(
                """# Ser OG to scissile carbonyl carbon
&rst
  iat=1,4,
  r1=2.60, r2=2.70, r3=2.90, r4=3.00,
/
# Ser OG - carbonyl C - carbonyl O attack angle
&rst
  iat=1,4,5,
  r1=90.00, r2=100.00, r3=110.00, r4=120.00,
/
# Leaving oxygen to catalytic His acceptor
&rst
  iat=6,3,
  r1=2.60, r2=2.70, r3=2.90, r4=3.00,
/
# Carbonyl oxygen to oxyanion-hole donor
&rst
  iat=5,2,
  r1=2.60, r2=2.70, r3=2.90, r4=3.00,
/
""",
                encoding="utf-8",
            )
            source_mdin = relax_dir / "00_restrained_mm_min.in"
            source_mdin.write_text(
                """Stage1 reactive-pose restrained MM minimization
&cntrl
  imin=1, maxcyc=500, ncyc=250,
  nmropt=1,
/
&wt
  type='END',
/
DISANG=reactive_relaxation_restraints.RST
LISTOUT=POUT
""",
                encoding="utf-8",
            )
            relaxation_manifest = tmpdir / "amber_reactive_relaxation_manifest.tsv"
            relaxation_manifest.write_text(
                "\n".join(
                    [
                        "relaxation_job_id\tpose_id\ttemplate_pdb\tsubstrate_model\tgeneration_method\tsource_complex_pdb\tcomplex_for_amber_pdb_path\trelaxation_job_dir\trestraint_path\tmm_min_mdin_path\trunner_path\tstatus",
                        (
                            "AMBER_RELAX_TEST\tPOSE1\t6EQE\tBHET_like\tvina\t"
                            f"{complex_for_amber}\t{complex_for_amber}\t{relax_dir}\t"
                            f"{source_rst}\t{source_mdin}\t{relax_dir / 'run.sh'}\t"
                            "reactive_relaxation_inputs_ready_requires_amber_topology"
                        ),
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            complex_leap = prep_dir / "complex_leap.pdb"
            complex_leap.write_text(
                "\n".join(
                    [
                        pdb_line("ATOM", 101, "OG", "SER", " ", 132, 10.0, 0.0, 0.0, "O"),
                        pdb_line("ATOM", 102, "N", "MET", " ", 133, 12.6, 0.8, 0.0, "N"),
                        pdb_line("ATOM", 103, "NE2", "HIE", " ", 209, 13.8, 0.0, 0.0, "N"),
                        pdb_line("HETATM", 104, "C005", "LIG", " ", 266, 13.5, 0.0, 0.0, "C"),
                        pdb_line("HETATM", 105, "O006", "LIG", " ", 266, 14.7, 0.0, 0.0, "O"),
                        pdb_line("HETATM", 106, "O004", "LIG", " ", 266, 13.5, 1.3, 0.0, "O"),
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
                        (
                            "TOPOLOGY_AMBER_RELAX_TEST\tAMBER_TOPO_FOR_AMBER_RELAX_TEST\tPOSE1\tBHET_like\t"
                            f"{complex_for_amber}\t{prep_dir}\t{prep_dir / 'complex_for_leap.pdb'}\t"
                            f"{prep_dir / 'ligand.pdb'}\t{prep_dir / 'tleap.in'}\t"
                            f"{prep_dir / 'run_amber_topology_prep.sh'}\t{prep_dir / 'complex.prmtop'}\t"
                            f"{prep_dir / 'complex.inpcrd'}\ttopology_prep_ready_requires_ambertools_execution\t"
                            "blind_stage1_amber_reactive_relaxation_manifest"
                        ),
                    ]
                )
                + "\n",
                encoding="utf-8",
            )

            result = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--relaxation-manifest",
                    str(relaxation_manifest),
                    "--topology-prep-manifest",
                    str(topology_manifest),
                    "--out-dir",
                    str(out_dir),
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            mapped_manifest = out_dir / "amber_reactive_relaxation_topology_mapped_manifest.tsv"
            self.assertTrue(mapped_manifest.exists())
            with mapped_manifest.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["mapped_iats"], "101,104;101,104,105;106,103;105,102")
            self.assertEqual(rows[0]["status"], "reactive_relaxation_restraints_mapped_to_amber_topology")

            mapped_dir = Path(rows[0]["mapped_relaxation_job_dir"])
            mapped_rst = (mapped_dir / "reactive_relaxation_restraints_mapped.RST").read_text(encoding="utf-8")
            self.assertIn("iat=101,104", mapped_rst)
            self.assertIn("iat=101,104,105", mapped_rst)
            self.assertIn("iat=106,103", mapped_rst)
            self.assertIn("iat=105,102", mapped_rst)
            self.assertNotIn("iat=1,4", mapped_rst)

            mapped_mdin = (mapped_dir / "00_restrained_mm_min_mapped.in").read_text(encoding="utf-8")
            self.assertIn("DISANG=reactive_relaxation_restraints_mapped.RST", mapped_mdin)
            runner = (mapped_dir / "run_reactive_mm_relaxation_mapped.sh").read_text(encoding="utf-8")
            self.assertIn("complex.prmtop", runner)
            self.assertIn("complex.inpcrd", runner)
            self.assertIn("00_restrained_mm_min_mapped.in", runner)


if __name__ == "__main__":
    unittest.main()
