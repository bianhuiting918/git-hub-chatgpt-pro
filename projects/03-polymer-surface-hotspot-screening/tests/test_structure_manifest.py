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


def write_authority_tsv(path: Path, rows):
    fields = [
        "authority_source", "authority_row_index", "family", "candidate_id",
        "sequence_md5", "status", "reason",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


class StructureManifestCliTests(unittest.TestCase):
    def test_separates_nylon_authority_candidates_from_external_controls(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pet_pdb = root / "pet.pdb"
            nyl_pdb = root / "nyl.pdb"
            recovered_pdb = root / "recovered.pdb"
            control_pdb = root / "control.pdb"
            for path, text in (
                (pet_pdb, "ATOM PET\n"),
                (nyl_pdb, "ATOM NYL\n"),
                (recovered_pdb, "ATOM RECOVERED\n"),
                (control_pdb, "ATOM CONTROL\n"),
            ):
                path.write_text(text, encoding="utf-8")

            pet = root / "pet.tsv"
            nylon = root / "nylon.tsv"
            recovered = root / "recovered.tsv"
            extra = root / "extra.tsv"
            authority = root / "authority.tsv"
            write_tsv(pet, [{
                "family": "PETase", "candidate_id": "PET_OK", "sequence_md5": "pet_md5",
                "selected_chain": "A", "receptor_path": str(pet_pdb),
                "receptor_sha256": sha256(pet_pdb), "input_status": "READY_FOR_STRUCTURE_EVALUATION",
                "structure_provenance": "PREDICTED_ALPHAFOLD",
            }])
            write_tsv(nylon, [{
                "family": "Nylonase", "candidate_id": "AUTH_READY", "sequence_md5": "auth_ready",
                "selected_chain": "A", "receptor_path": str(nyl_pdb),
                "receptor_sha256": sha256(nyl_pdb), "input_status": "READY_FOR_STRUCTURE_EVALUATION",
                "structure_provenance": "PREDICTED_ESMFOLD",
            }])
            write_tsv(recovered, [{
                "family": "Nylonase", "candidate_id": "AUTH_RECOVERED", "sequence_md5": "auth_recovered",
                "selected_chain": "A", "receptor_path": str(recovered_pdb),
                "receptor_sha256": sha256(recovered_pdb), "input_status": "READY_FOR_STRUCTURE_EVALUATION",
                "structure_provenance": "EXPERIMENTAL_PDB",
            }])
            write_tsv(extra, [{
                "family": "Nylonase", "candidate_id": "EXTERNAL_CONTROL", "sequence_md5": "external",
                "selected_chain": "A", "receptor_path": str(control_pdb),
                "receptor_sha256": sha256(control_pdb), "input_status": "READY_FOR_STRUCTURE_EVALUATION",
                "structure_provenance": "EXPERIMENTAL_PDB",
            }])
            write_authority_tsv(authority, [
                {"authority_source": "prepared", "authority_row_index": "1", "family": "Nylonase",
                 "candidate_id": "AUTH_READY", "sequence_md5": "auth_ready",
                 "status": "READY_FOR_STRUCTURE_EVALUATION", "reason": ""},
                {"authority_source": "not_evaluated", "authority_row_index": "2", "family": "Nylonase",
                 "candidate_id": "AUTH_RECOVERED", "sequence_md5": "auth_recovered",
                 "status": "NOT_EVALUATED_NO_EXACT_STRUCTURE", "reason": "old state before recovery"},
                {"authority_source": "dell", "authority_row_index": "3", "family": "Nylonase",
                 "candidate_id": "AUTH_NOT_ON_SUGON", "sequence_md5": "auth_missing",
                 "status": "READY_AFTER_DELL_SYNC", "reason": "not on Sugon"},
            ])

            out = root / "out"
            run = subprocess.run([
                sys.executable, str(SCRIPT),
                "--pet-manifest", str(pet),
                "--nylon-authority-manifest", str(authority),
                "--nylon-manifest", str(nylon),
                "--nylon-recovered-manifest", str(recovered),
                "--nylon-extra-manifest", str(extra),
                "--output-dir", str(out),
            ], text=True, capture_output=True)
            self.assertEqual(run.returncode, 0, run.stderr)

            with (out / "nylon_structures.tsv").open(encoding="utf-8") as handle:
                candidates = list(csv.DictReader(handle, delimiter="\t"))
            with (out / "nylon_external_controls.tsv").open(encoding="utf-8") as handle:
                controls = list(csv.DictReader(handle, delimiter="\t"))
            with (out / "nylon_authority_not_evaluated.tsv").open(encoding="utf-8") as handle:
                missing = list(csv.DictReader(handle, delimiter="\t"))
            summary = json.loads((out / "structure_manifest_summary.json").read_text(encoding="utf-8"))

            self.assertEqual({r["candidate_id"] for r in candidates}, {"AUTH_READY", "AUTH_RECOVERED"})
            self.assertEqual([r["candidate_id"] for r in controls], ["EXTERNAL_CONTROL"])
            self.assertEqual([r["candidate_id"] for r in missing], ["AUTH_NOT_ON_SUGON"])
            self.assertEqual(missing[0]["exclusion_reason"], "NOT_EVALUATED_STRUCTURE_NOT_ON_SUGON")
            self.assertEqual(summary["nylon_authority"]["denominator"], 3)
            self.assertEqual(summary["nylon_authority"]["structure_evaluable"], 2)
            self.assertEqual(summary["nylon_authority"]["not_evaluated"], 1)
            self.assertEqual(summary["nylon_external_controls"]["structure_evaluable"], 1)

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
            authority = root / "authority.tsv"
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
            write_authority_tsv(authority, [
                {"authority_source": "prepared", "authority_row_index": "1", "family": "Nylonase",
                 "candidate_id": "NYL_OK", "sequence_md5": "nyl_md5",
                 "status": "READY_FOR_STRUCTURE_EVALUATION", "reason": ""},
            ])

            out = root / "out"
            run = subprocess.run([
                sys.executable, str(SCRIPT),
                "--pet-manifest", str(pet),
                "--nylon-authority-manifest", str(authority),
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
            with (out / "nylon_external_controls.tsv").open(encoding="utf-8") as handle:
                control_rows = list(csv.DictReader(handle, delimiter="\t"))
            with (out / "not_evaluated.tsv").open(encoding="utf-8") as handle:
                rejected = list(csv.DictReader(handle, delimiter="\t"))
            summary = json.loads((out / "structure_manifest_summary.json").read_text(encoding="utf-8"))

            self.assertEqual([r["candidate_id"] for r in pet_rows], ["PET_OK"])
            self.assertEqual([r["candidate_id"] for r in nylon_rows], ["NYL_OK"])
            self.assertEqual([r["candidate_id"] for r in control_rows], ["NYL_EXTRA"])
            self.assertTrue(all(r["material_family"] == "PET" for r in pet_rows))
            self.assertTrue(all(r["material_family"] == "NYLON" for r in nylon_rows + control_rows))
            reasons = {r["exclusion_reason"] for r in rejected}
            self.assertIn("NOT_EVALUATED_MISSING_STRUCTURE", reasons)
            self.assertIn("DUPLICATE_LOWER_PRIORITY", reasons)
            self.assertEqual(summary["included"]["PET"], 1)
            self.assertEqual(summary["nylon_authority"]["structure_evaluable"], 1)
            self.assertEqual(summary["nylon_external_controls"]["structure_evaluable"], 1)

    def test_checksum_mismatch_is_not_included(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            pdb = root / "bad.pdb"
            pdb.write_text("ATOM\n", encoding="utf-8")
            pet = root / "pet.tsv"
            empty = root / "empty.tsv"
            authority = root / "authority.tsv"
            write_tsv(pet, [{
                "family": "PETase", "candidate_id": "BAD", "sequence_md5": "bad_md5",
                "selected_chain": "A", "receptor_path": str(pdb),
                "receptor_sha256": "f" * 64, "input_status": "READY_FOR_STRUCTURE_EVALUATION",
                "structure_provenance": "PREDICTED_ALPHAFOLD",
            }])
            write_tsv(empty, [])
            write_authority_tsv(authority, [])
            out = root / "out"
            run = subprocess.run([
                sys.executable, str(SCRIPT),
                "--pet-manifest", str(pet),
                "--nylon-authority-manifest", str(authority),
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


    def test_exact_dell_manifest_overrides_older_nylon_structure(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            old_pdb = root / "old.pdb"
            exact_pdb = root / "exact.pdb"
            old_pdb.write_text("ATOM OLD\n", encoding="utf-8")
            exact_pdb.write_text("ATOM EXACT\n", encoding="utf-8")
            empty = root / "empty.tsv"
            write_tsv(empty, [])
            primary = root / "primary.tsv"
            exact = root / "exact.tsv"
            authority = root / "authority.tsv"
            write_tsv(primary, [{
                "family": "Nylonase", "candidate_id": "NYL_OLD", "sequence_md5": "same_md5",
                "selected_chain": "A", "receptor_path": str(old_pdb),
                "receptor_sha256": sha256(old_pdb), "input_status": "READY_FOR_STRUCTURE_EVALUATION",
                "structure_provenance": "PREDICTED_ESMFOLD",
            }])
            write_tsv(exact, [{
                "family": "Nylonase", "candidate_id": "NYL_EXACT", "sequence_md5": "same_md5",
                "selected_chain": "A", "receptor_path": str(exact_pdb),
                "receptor_sha256": sha256(exact_pdb), "input_status": "READY_FOR_STRUCTURE_EVALUATION",
                "structure_provenance": "PREDICTED_ALPHAFOLD",
            }])
            write_authority_tsv(authority, [{
                "authority_source": "dell", "authority_row_index": "1", "family": "Nylonase",
                "candidate_id": "NYL_EXACT", "sequence_md5": "same_md5",
                "status": "READY_AFTER_DELL_SYNC", "reason": "",
            }])
            out = root / "out"
            run = subprocess.run([
                sys.executable, str(SCRIPT),
                "--pet-manifest", str(empty),
                "--nylon-authority-manifest", str(authority),
                "--nylon-exact-manifest", str(exact),
                "--nylon-manifest", str(primary),
                "--nylon-recovered-manifest", str(empty),
                "--nylon-extra-manifest", str(empty),
                "--output-dir", str(out),
            ], text=True, capture_output=True)
            self.assertEqual(run.returncode, 0, run.stderr)
            with (out / "nylon_structures.tsv").open(encoding="utf-8") as handle:
                rows = list(csv.DictReader(handle, delimiter="\t"))
            self.assertEqual(len(rows), 1)
            self.assertEqual(rows[0]["candidate_id"], "NYL_EXACT")
            self.assertEqual(rows[0]["manifest_role"], "nylon_exact_dell")
            self.assertEqual(rows[0]["source_priority"], "-1")
            summary = json.loads((out / "structure_manifest_summary.json").read_text(encoding="utf-8"))
            self.assertEqual(summary["schema_version"], "structure_manifest_v3")
            self.assertEqual(summary["source_priority"]["nylon_exact_dell"], -1)


if __name__ == "__main__":
    unittest.main()
