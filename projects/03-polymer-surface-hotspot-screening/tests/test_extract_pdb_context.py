import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "extract_pdb_context.py"


def atom_line(serial: int, record: str, atom: str, residue: str, chain: str, resseq: int) -> str:
    return (
        f"{record:<6}{serial:>5} {atom:^4} {residue:>3} {chain:1}{resseq:>4}    "
        f"{serial:>8.3f}{0.0:>8.3f}{0.0:>8.3f}{1.0:>6.2f}{20.0:>6.2f}          C  \n"
    )


class ExtractPdbContextCliTests(unittest.TestCase):
    def test_extracts_only_requested_monomer_chain(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.pdb"
            output = root / "chain_a.pdb"
            source.write_text(
                "HEADER    TEST\n"
                + atom_line(1, "ATOM", "CA", "ALA", "A", 1)
                + atom_line(2, "HETATM", "NA", "NA", "A", 401)
                + atom_line(3, "ATOM", "CA", "ALA", "B", 1)
                + "END\n",
                encoding="ascii",
            )
            run = subprocess.run(
                [sys.executable, str(SCRIPT), "--input", str(source), "--output", str(output), "--chains", "A"],
                text=True,
                capture_output=True,
            )
            self.assertEqual(run.returncode, 0, run.stderr)
            text = output.read_text(encoding="ascii")
            coordinate_lines = [line for line in text.splitlines() if line.startswith(("ATOM  ", "HETATM"))]
            self.assertEqual([line[21] for line in coordinate_lines], ["A", "A"])
            self.assertTrue(text.endswith("END\n"))

    def test_protein_only_excludes_chain_heteroatoms(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.pdb"
            output = root / "protein_only.pdb"
            source.write_text(
                atom_line(1, "ATOM", "CA", "ALA", "A", 1)
                + atom_line(2, "HETATM", "NA", "NA", "A", 401)
                + "END\n",
                encoding="ascii",
            )
            run = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--input",
                    str(source),
                    "--output",
                    str(output),
                    "--chains",
                    "A",
                    "--protein-only",
                ],
                text=True,
                capture_output=True,
            )
            self.assertEqual(run.returncode, 0, run.stderr)
            coordinate_lines = [
                line for line in output.read_text(encoding="ascii").splitlines()
                if line.startswith(("ATOM  ", "HETATM"))
            ]
            self.assertEqual(len(coordinate_lines), 1)
            self.assertTrue(coordinate_lines[0].startswith("ATOM  "))
            self.assertIn("REMARK 999 RECORD POLICY PROTEIN_ONLY", output.read_text(encoding="ascii"))

    def test_extracts_exact_abcd_tetramer_not_fifteen_chain_asymmetric_unit(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.pdb"
            output = root / "abcd.pdb"
            source.write_text(
                "HEADER    TEST\n"
                + "".join(atom_line(index, "ATOM", "CA", "ALA", chain, 1)
                          for index, chain in enumerate("ABCDEO", start=1))
                + "END\n",
                encoding="ascii",
            )
            run = subprocess.run(
                [
                    sys.executable,
                    str(SCRIPT),
                    "--input",
                    str(source),
                    "--output",
                    str(output),
                    "--chains",
                    "A,B,C,D",
                ],
                text=True,
                capture_output=True,
            )
            self.assertEqual(run.returncode, 0, run.stderr)
            coordinate_lines = [
                line for line in output.read_text(encoding="ascii").splitlines()
                if line.startswith(("ATOM  ", "HETATM"))
            ]
            self.assertEqual([line[21] for line in coordinate_lines], list("ABCD"))

    def test_missing_requested_chain_fails_without_output(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            source = root / "source.pdb"
            output = root / "missing.pdb"
            source.write_text(atom_line(1, "ATOM", "CA", "ALA", "A", 1) + "END\n", encoding="ascii")
            run = subprocess.run(
                [sys.executable, str(SCRIPT), "--input", str(source), "--output", str(output), "--chains", "Z"],
                text=True,
                capture_output=True,
            )
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("missing requested chains", run.stderr)
            self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
