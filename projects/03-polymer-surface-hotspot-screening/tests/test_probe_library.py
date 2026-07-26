import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from rdkit import Chem


PROJECT = Path(__file__).resolve().parents[1]
MANIFEST = PROJECT / "config" / "probes.tsv"
BUILDER = PROJECT / "scripts" / "build_probe_library.py"

REQUIRED_COLUMNS = {
    "probe_id",
    "material_family",
    "role",
    "preferred_name",
    "canonical_smiles",
    "formal_charge",
    "heavy_atom_count",
    "total_atom_count",
    "primary_or_control",
    "stage",
    "rationale",
    "source",
    "status",
}


def read_manifest():
    with MANIFEST.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class ProbeLibraryTests(unittest.TestCase):
    def test_manifest_has_required_schema_and_separate_material_families(self):
        with MANIFEST.open(encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            self.assertTrue(REQUIRED_COLUMNS.issubset(reader.fieldnames or []))
            rows = list(reader)

        active = [row for row in rows if row["status"] != "deferred"]
        deferred = [row for row in rows if row["status"] == "deferred"]
        self.assertEqual(len(active), 13)
        self.assertEqual(len(deferred), 3)
        self.assertEqual(
            {row["probe_id"] for row in deferred},
            {"PA6_CAPPED_DIMER", "PA66_CAPPED_REPEAT", "PA66_CAPPED_DIMER"},
        )
        self.assertTrue(all(row["source"].strip() for row in active))
        self.assertEqual(
            len([row for row in active if row["material_family"] == "PET"]),
            6,
        )
        self.assertEqual(
            len([row for row in active if row["material_family"] != "PET"]),
            7,
        )

    def test_declared_active_probe_atom_counts_match_rdkit(self):
        for row in read_manifest():
            if row["status"] == "deferred":
                continue
            with self.subTest(probe_id=row["probe_id"]):
                mol = Chem.MolFromSmiles(row["canonical_smiles"])
                self.assertIsNotNone(mol)
                mol_h = Chem.AddHs(mol)
                self.assertEqual(Chem.GetFormalCharge(mol), int(row["formal_charge"]))
                self.assertEqual(mol.GetNumHeavyAtoms(), int(row["heavy_atom_count"]))
                self.assertEqual(mol_h.GetNumAtoms(), int(row["total_atom_count"]))
                self.assertEqual(
                    Chem.MolToSmiles(mol, canonical=True),
                    row["canonical_smiles"],
                )

    def test_builder_emits_validated_sdf_metadata_gate_and_checksums(self):
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir)
            completed = subprocess.run(
                [
                    sys.executable,
                    str(BUILDER),
                    "--manifest",
                    str(MANIFEST),
                    "--output-dir",
                    str(output_dir),
                ],
                check=False,
                capture_output=True,
                text=True,
            )
            self.assertEqual(completed.returncode, 0, completed.stderr)

            sdf_path = output_dir / "probes.sdf"
            table_path = output_dir / "validated_probes.tsv"
            gate_path = output_dir / "PROBES_PASS.json"
            checksum_path = output_dir / "PROBE_SHA256SUMS"
            for path in (sdf_path, table_path, gate_path, checksum_path):
                self.assertTrue(path.is_file(), path)

            molecules = [
                mol
                for mol in Chem.SDMolSupplier(str(sdf_path), removeHs=False)
                if mol is not None
            ]
            self.assertEqual(len(molecules), 13)
            self.assertTrue(all(mol.GetNumConformers() == 1 for mol in molecules))

            with table_path.open(encoding="utf-8", newline="") as handle:
                validated = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(validated), 13)
            self.assertTrue(all(row["validation_status"] == "PASS" for row in validated))

            gate = json.loads(gate_path.read_text(encoding="utf-8"))
            self.assertEqual(gate["status"], "PROBES_PASS")
            self.assertEqual(gate["n_active"], 13)
            self.assertEqual(gate["n_deferred"], 3)
            self.assertEqual(gate["scientific_scope"], "probe_library_only")

            for line in checksum_path.read_text(encoding="utf-8").splitlines():
                expected, filename = line.split("  ", maxsplit=1)
                observed = hashlib.sha256((output_dir / filename).read_bytes()).hexdigest()
                self.assertEqual(observed, expected)


if __name__ == "__main__":
    unittest.main()
