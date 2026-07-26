import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


PROJECT = Path(__file__).resolve().parents[1]
CHECKER = PROJECT / "scripts" / "check_receptor.py"
WRAPPER = PROJECT / "scripts" / "prepare_receptor.sh"


def atom_line(
    record,
    serial,
    name,
    altloc,
    residue,
    chain,
    resseq,
    x,
    y,
    z,
    element,
    occupancy=1.0,
):
    return (
        f"{record:<6}{serial:5d} {name:^4}{altloc:1}{residue:>3} {chain:1}"
        f"{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
        f"{occupancy:6.2f}{20.0:6.2f}          {element:>2}\n"
    )


def fixture_pdb():
    lines = [
        "MODEL        1\n",
        atom_line("ATOM", 1, "CA", "", "ALA", "A", 10, 0, 0, 0, "C"),
        atom_line("ATOM", 2, "OG", "A", "SER", "A", 20, 2, 0, 0, "O", 0.6),
        atom_line("ATOM", 3, "OG", "B", "SER", "A", 20, 2.2, 0, 0, "O", 0.4),
        atom_line("HETATM", 4, "ZN", "", "ZN", "A", 401, 4, 0, 0, "ZN"),
        atom_line("HETATM", 5, "C1", "", "EDO", "A", 501, 6, 0, 0, "C"),
        atom_line("ATOM", 6, "CA", "", "GLY", "B", 10, 8, 0, 0, "C"),
        "ENDMDL\n",
        "MODEL        2\n",
        atom_line("ATOM", 7, "CA", "", "ALA", "A", 10, 20, 0, 0, "C"),
        atom_line("ATOM", 8, "OG", "", "SER", "A", 20, 22, 0, 0, "O"),
        "ENDMDL\n",
        "END\n",
    ]
    return "".join(lines)


class ReceptorSelectionTests(unittest.TestCase):
    def run_checker(self, *args):
        return subprocess.run(
            [sys.executable, str(CHECKER), *map(str, args)],
            check=False,
            capture_output=True,
            text=True,
        )

    def test_selects_model_chain_altloc_and_hetero_allowlist_deterministically(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "fixture.pdb"
            source.write_text(fixture_pdb(), encoding="ascii")
            outputs = []
            metadata = []
            for suffix in ("a", "b"):
                output = root / f"selected_{suffix}.pdb"
                meta = root / f"selected_{suffix}.json"
                completed = self.run_checker(
                    "select",
                    "--input",
                    source,
                    "--output",
                    output,
                    "--metadata",
                    meta,
                    "--chains",
                    "A",
                    "--model",
                    "1",
                    "--catalytic-residues",
                    "A:10,A:20",
                    "--retain-hetero",
                    "A:ZN:401",
                )
                self.assertEqual(completed.returncode, 0, completed.stderr)
                outputs.append(output.read_bytes())
                metadata.append(json.loads(meta.read_text(encoding="utf-8")))

            self.assertEqual(outputs[0], outputs[1])
            self.assertEqual(metadata[0], metadata[1])
            text = outputs[0].decode("ascii")
            self.assertIn(" ZN ", text)
            self.assertNotIn(" EDO ", text)
            self.assertNotIn(" B  10", text)
            self.assertNotIn(" 20.000", text)
            self.assertEqual(text.count(" OG "), 1)
            self.assertNotIn("OG  B", text)

            meta = metadata[0]
            self.assertEqual(meta["status"], "SELECTION_PASS")
            self.assertEqual(meta["selected_model"], 1)
            self.assertEqual(meta["selected_chains"], ["A"])
            self.assertEqual(meta["catalytic_mapping_status"], "PASS")
            self.assertEqual(meta["retained_hetero"], ["A:ZN:401"])
            self.assertEqual(meta["removed_hetero"], ["A:EDO:501"])
            self.assertEqual(meta["alternate_location_policy"], "blank_then_A_then_occupancy")

    def test_missing_catalytic_residue_is_not_evaluated_not_biological_failure(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            root = Path(temp_dir)
            source = root / "fixture.pdb"
            output = root / "selected.pdb"
            meta = root / "selected.json"
            source.write_text(fixture_pdb(), encoding="ascii")
            completed = self.run_checker(
                "select",
                "--input",
                source,
                "--output",
                output,
                "--metadata",
                meta,
                "--chains",
                "A",
                "--model",
                "1",
                "--catalytic-residues",
                "A:10,A:999",
            )
            self.assertEqual(completed.returncode, 4)
            self.assertFalse(output.exists())
            payload = json.loads(meta.read_text(encoding="utf-8"))
            self.assertEqual(payload["status"], "NOT_EVALUATED_CATALYTIC_MAPPING")
            self.assertEqual(payload["missing_catalytic_residues"], ["A:999"])

    def test_wrapper_uses_fixed_adfr_preparation_and_refuses_destructive_cleanup(self):
        text = WRAPPER.read_text(encoding="utf-8")
        self.assertIn("check_receptor.py", text)
        self.assertIn("prepare_receptor", text)
        self.assertIn("-A checkhydrogens", text)
        self.assertIn("-U nphs_lps_waters", text)
        self.assertIn("sha256sum", text)
        self.assertIn("metadata", text)
        self.assertNotIn("rm -rf", text)
        self.assertNotIn("nonstdres", text)


if __name__ == "__main__":
    unittest.main()
