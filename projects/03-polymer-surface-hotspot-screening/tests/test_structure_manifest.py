import csv
import hashlib
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "build_structure_manifest.py"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_tsv(path: Path, rows):
    fields = [
        "family", "candidate_id", "sequence_md5", "selected_chain",
        "receptor_path", "receptor_sha256", "input_status",
        "structure_provenance",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class StructureManifestCliTests(unittest.TestCase):
    def test_freezes_only_readable_checksum_valid_structures_and_keeps_families_separate(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pet_pdb = root / "pet.pdb"
            nyl_pdb = root / "nyl.pdb"
            extra_pdb = root / "nyl_extra.pdb"
            pet_pdb.write_text("ATOM PET\n", encoding="utf-8")
            nyl_pdb.write_text("ATOM NYL\n", encoding="utf-8")
            extra_pdb.write_text("ATOM NYL EXTRA\n", encoding="utf-8")

            pet = root / "pet.tsv"
            nylon = root / "nylon.tsv"
            recovered = root / "recovered.tsv"
            extra = root / "extra.tsv"
            write_tsv(pet, [
                {"family": "PETase", "candidate_id": "PET_OK", "sequence_md5": "pet_md5",
                 "selected_chain": "A", "receptor_path": str(pet_pdb),
                 "receptor_sha256": sha256(pet_pdb), "input_status": "READY_FOR_STRUCTURE_EVALUATION",
                 "structure_provenance": "PREDICTED_ALPHAFOLD"},
                {"family": "PETase", "candidate_id": "PET_MISSING", "sequence_md5": "pet_missing",
                 "selected_chain": "A", "receptor_path": str(root / "missing.pdb"),
                 "receptor_sha256": "0" * 64, "input_status": "READY_FOR_STRUCTURE_EVALUATION",
                 "structure_provenance": "PREDICTED_ALPHAFOLD"},
            ])
            write_tsv(nylon, [
                {"family": "Nylonase", "candidate_id": "NYL_OK", "sequence_md5": "nyl_md5",
                 "selected_chain": "A", "receptor_path": str(nyl_pdb),
                 "receptor_sha256": sha256(nyl_pdb), "input_status": "READY_FOR_STRUCTURE_EVALUATION",
                 "structure_provenance": "PREDICTED_ESMFOLD"},
            ])
            write_tsv(recovered, [
                {"family": "Nylonase", "candidate_id": "NYL_DUP", "sequence_md5": "nyl_md5",
                 "selected_chain": "A", "receptor_path": str(extra_pdb),
                 "receptor_sha256": sha256(extra_pdb), "input_status": "READY_FOR_STRUCTURE_EVALUATION",
                 "structure_provenance": "EXPERIMENTAL_PDB"},
            ])
            write_tsv(extra, [
                {"family": "Nylonase", "candidate_id": "NYL_EXTRA", "sequence_md5": "nyl_extra",
                 "selected_chain": "A", "receptor_path": str(extra_pdb),
                 "receptor_sha256": sha256(extra_pdb), "input_status": "READY_FOR_STRUCTURE_EVALUATION",
                 "structure_provenance": "EXPERIMENTAL_PDB"},
            ])

            out = root / "out"
            run = subprocess.run([
                sys.executable, str(SCRIPT),
                "--pet-manifest", str(pet),
                "--nylon-manifest", str(nylon),
                "--nylon-recovered-manifest", str(recovered),
                "--nylon-extra-manifest", str(extra),
                "--output-dir", str(out),
            ], text=True, capture_output=True)
            self.assertEqual(run.returncode, 0, run.stderr)

            with (out / "pet_structures.tsv").open(encoding="utf-8") as handle:
                pet_rows = list(csv.DictReader(handle, delimiter="\t"))
            with (out / "nylon_structures.tsv").open(encoding="utf-8") as handle:
                nylon_rows = list(csv.DictReader(handle, delimiter="\t"))
            with (out / "not_evaluated.tsv").open(encoding="utf-8") as handle:
                rejected = list(csv.DictReader(handle, delimiter="\t"))
            summary = json.loads((out / "structure_manifest_summary.json").read_text(encoding="utf-8"))

            self.assertEqual([r["candidate_id"] for r in pet_rows], ["PET_OK"])
            self.assertEqual({r["candidate_id"] for r in nylon_rows}, {"NYL_OK", "NYL_EXTRA"})
            self.assertTrue(all(r["material_family"] == "PET" for r in pet_rows))
            self.assertTrue(all(r["material_family"] == "NYLON" for r in nylon_rows))
            reasons = {r["exclusion_reason"] for r in rejected}
            self.assertIn("NOT_EVALUATED_MISSING_STRUCTURE", reasons)
            self.assertIn("DUPLICATE_LOWER_PRIORITY", reasons)
            self.assertEqual(summary["included"]["PET"], 1)
            self.assertEqual(summary["included"]["NYLON"], 2)

    def test_checksum_mismatch_is_not_included(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pdb = root / "bad.pdb"
            pdb.write_text("ATOM\n", encoding="utf-8")
            pet = root / "pet.tsv"
            empty = root / "empty.tsv"
            write_tsv(pet, [{
                "family": "PETase", "candidate_id": "BAD", "sequence_md5": "bad_md5",
                "selected_chain": "A", "receptor_path": str(pdb),
                "receptor_sha256": "f" * 64, "input_status": "READY_FOR_STRUCTURE_EVALUATION",
                "structure_provenance": "PREDICTED_ALPHAFOLD",
            }])
            write_tsv(empty, [])
            out = root / "out"
            run = subprocess.run([
                sys.executable, str(SCRIPT),
                "--pet-manifest", str(pet),
                "--nylon-manifest", str(empty),
                "--nylon-recovered-manifest", str(empty),
                "--nylon-extra-manifest", str(empty),
                "--output-dir", str(out),
            ], text=True, capture_output=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            with (out / "pet_structures.tsv").open(encoding="utf-8") as handle:
                self.assertEqual(len(list(csv.DictReader(handle, delimiter="\t"))), 0)
            with (out / "not_evaluated.tsv").open(encoding="utf-8") as handle:
                reasons = {r["exclusion_reason"] for r in csv.DictReader(handle, delimiter="\t")}
            self.assertIn("NOT_EVALUATED_CHECKSUM_MISMATCH", reasons)


if __name__ == "__main__":
    unittest.main()
