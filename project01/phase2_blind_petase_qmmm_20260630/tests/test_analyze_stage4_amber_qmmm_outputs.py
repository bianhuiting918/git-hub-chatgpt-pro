import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


class AnalyzeStage4AmberQmmmOutputsTest(unittest.TestCase):
    def test_classifies_completed_nonconverged_and_incomplete_qmmm_outputs(self):
        script = Path("project01/phase2_blind_petase_qmmm_20260630/scripts/analyze_stage4_amber_qmmm_outputs.py")
        with tempfile.TemporaryDirectory() as tmp:
            tmpdir = Path(tmp)
            clean = tmpdir / "clean.out"
            clean.write_text(
                """Header\n   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER\n      5      -1.1000E+05     1.0000E-01     2.0000E+00     C1       10\n VDWAALS =      100.0000  EEL     =     -200.0000  HBOND      =        0.0000\n DFTBESCF=    -5000.2500\n  Maximum number of minimization cycles reached.\n                    FINAL RESULTS\n|           Run   done at 12:00:00.000  on 07/01/2026\n""",
                encoding="utf-8",
            )
            warned = tmpdir / "warned.out"
            warned.write_text(
                """QMMM SCC-DFTB: !!!! ============= WARNING ============= !!!!\n QMMM SCC-DFTB: Convergence could not be achieved in this step.\n   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER\n      1      -1.6217E+05     2.0795E+00     3.7811E+02     NE2      2991\n VDWAALS =    23525.9538  EEL     =  -213867.8381  HBOND      =        0.0000\n DFTBESCF=     2650.4126\n  Maximum number of minimization cycles reached.\n                    FINAL RESULTS\n|           Run   done at 12:02:00.000  on 07/01/2026\n""",
                encoding="utf-8",
            )
            incomplete = tmpdir / "incomplete.out"
            incomplete.write_text(
                """TESTING RELATIVE ERROR over r ranging from 0.0 to cutoff\nQMMM: SINGLET STATE CALCULATION\n   NSTEP       ENERGY          RMS            GMAX         NAME    NUMBER\n      1      -1.6217E+05     2.0795E+00     3.7811E+02     NE2      2991\n VDWAALS =    23525.9538  EEL     =  -213867.8381  HBOND      =        0.0000\n DFTBESCF=     2650.4126\n""",
                encoding="utf-8",
            )
            stderr = tmpdir / "stderr.log"
            stderr.write_text("", encoding="utf-8")
            manifest = tmpdir / "qmmm_health.tsv"

            result = subprocess.run(
                [
                    sys.executable,
                    str(script),
                    "--output-manifest",
                    str(manifest),
                    "--case",
                    f"clean={clean}:{stderr}",
                    "--case",
                    f"warned={warned}:{stderr}",
                    "--case",
                    f"incomplete={incomplete}:{stderr}",
                ],
                text=True,
                capture_output=True,
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            with manifest.open("r", encoding="utf-8", newline="") as handle:
                rows = {row["case_id"]: row for row in csv.DictReader(handle, delimiter="\t")}
            self.assertEqual(rows["clean"]["status"], "completed")
            self.assertEqual(rows["clean"]["scc_convergence_warnings"], "0")
            self.assertEqual(rows["clean"]["last_step"], "5")
            self.assertEqual(rows["clean"]["last_dftbescf"], "-5000.2500")
            self.assertEqual(rows["warned"]["status"], "completed_with_scc_warnings")
            self.assertEqual(rows["warned"]["scc_convergence_warnings"], "1")
            self.assertEqual(rows["warned"]["last_step"], "1")
            self.assertEqual(rows["incomplete"]["status"], "running_or_incomplete")
            self.assertEqual(rows["incomplete"]["run_done"], "0")


if __name__ == "__main__":
    unittest.main()

