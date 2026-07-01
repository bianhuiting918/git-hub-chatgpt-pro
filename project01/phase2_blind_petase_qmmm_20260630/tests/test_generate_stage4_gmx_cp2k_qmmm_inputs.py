import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class GenerateStage4GmxCp2kQmmmInputsTest(unittest.TestCase):
    def test_generates_gmx_cp2k_single_point_inputs_from_productive_conformers(self):
        script = Path("work/generate_stage4_gmx_cp2k_qmmm_inputs.py")
        if not script.exists():
            script = Path("project01/phase2_blind_petase_qmmm_20260630/scripts/generate_stage4_gmx_cp2k_qmmm_inputs.py")
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            conformers = tmpdir / "productive_conformer_manifest.tsv"
            structure = tmpdir / "conf.gro"
            topology = tmpdir / "topol.top"
            index = tmpdir / "index.ndx"
            structure.write_text("PETase test\n    1\n    1ALA      N    1   0.000   0.000   0.000\n   1.0 1.0 1.0\n", encoding="utf-8")
            topology.write_text("[ system ]\nPETase\n", encoding="utf-8")
            index.write_text("[ QMatoms ]\n1\n", encoding="utf-8")
            conformers.write_text(
                "\n".join(
                    [
                        "conformer_id\tsource_pose_id\tmd_job_id\tstage\tcluster_id\trepresentative_structure_path\tframe_count\tgeometry_pass_fraction\tselection_status\tsource\ttopology_path\tindex_path\tqm_atom_group",
                        f"CONF_001\tREACTIVE_6EQE_BHET_like_E01_001\tMD_REACTIVE_6EQE_BHET_like_E01_001_R01\tpreacylation_michaelis_complex\tCL001\t{structure}\t250\t0.62\tproductive\tblind_stage2_md_cluster\t{topology}\t{index}\tQMatoms",
                    ]
                )
                + "\n",
                encoding="utf-8",
            )
            out_dir = tmpdir / "qmmm_inputs"

            result = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--productive-conformers",
                    str(conformers),
                    "--out-dir",
                    str(out_dir),
                    "--max-jobs",
                    "1",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            manifest = out_dir / "gmx_cp2k_qmmm_job_manifest.tsv"
            self.assertTrue(manifest.exists())
            with manifest.open("r", encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["qmmm_job_id"], "QMMM_CONF_001_SP001")
            self.assertEqual(rows[0]["status"], "inputs_ready_requires_grompp")
            self.assertEqual(rows[0]["energy_source"], "gromacs_log_potential_energy")

            job_dir = Path(rows[0]["job_dir"])
            for name in ("qmmm_sp.mdp", "topol-qmmm.inp", "run_gmx_cp2k_qmmm.sh", "extract_potential_energy.sh", "00_README.md"):
                path = job_dir / name
                self.assertTrue(path.exists(), f"missing {path}")
                lowered = path.read_text(encoding="utf-8").lower()
                self.assertNotIn("paper trajectory", lowered)
                self.assertNotIn("paper reaction coordinate", lowered)
                self.assertNotIn("paper barrier", lowered)

            runner = (job_dir / "run_gmx_cp2k_qmmm.sh").read_text(encoding="utf-8")
            self.assertIn("/Dell/Dell14/bianht/gmx_cp2k_patched.sh", runner)
            self.assertIn("grompp", runner)
            self.assertIn("mdrun", runner)
            self.assertIn("-nsteps 1", runner)

            extractor = (job_dir / "extract_potential_energy.sh").read_text(encoding="utf-8")
            self.assertIn("Potential Energy", extractor)
            self.assertNotIn("Total FORCE_EVAL", extractor)


if __name__ == "__main__":
    unittest.main()
