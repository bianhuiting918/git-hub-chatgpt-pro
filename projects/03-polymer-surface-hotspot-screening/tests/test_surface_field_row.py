import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_surface_field_row.py"


def load_module():
    spec = importlib.util.spec_from_file_location("run_surface_field_row", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def write_manifest(path: Path, receptor: Path, digest: str):
    fields = [
        "material_family", "candidate_universe", "candidate_id", "sequence_md5",
        "selected_chain", "receptor_path", "receptor_sha256", "input_status",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerow({
            "material_family": "PET",
            "candidate_universe": "PET_AUTHORITY",
            "candidate_id": "PET_1",
            "sequence_md5": "a" * 32,
            "selected_chain": "A",
            "receptor_path": str(receptor),
            "receptor_sha256": digest,
            "input_status": "READY_FOR_STRUCTURE_EVALUATION",
        })


def test_loads_zero_based_manifest_row_and_validates_hash(tmp_path):
    module = load_module()
    receptor = tmp_path / "r.pdb"
    receptor.write_text("ATOM\n", encoding="utf-8")
    digest = hashlib.sha256(receptor.read_bytes()).hexdigest()
    manifest = tmp_path / "m.tsv"
    write_manifest(manifest, receptor, digest)
    row = module.load_manifest_row(manifest, 0)
    assert row["candidate_id"] == "PET_1"
    assert module.validate_manifest_row(row, "PET") == ""


def test_bad_hash_is_technical_not_evaluated(tmp_path):
    module = load_module()
    receptor = tmp_path / "r.pdb"
    receptor.write_text("ATOM\n", encoding="utf-8")
    manifest = tmp_path / "m.tsv"
    write_manifest(manifest, receptor, "0" * 64)
    row = module.load_manifest_row(manifest, 0)
    assert module.validate_manifest_row(row, "PET") == "NOT_EVALUATED_RECEPTOR_CHECKSUM_MISMATCH"


def test_existing_field_pass_is_reused(tmp_path):
    module = load_module()
    target = tmp_path / "target"
    target.mkdir()
    (target / "FIELD_PASS.json").write_text(
        json.dumps({"status": "FIELD_PASS"}), encoding="utf-8"
    )
    assert module.existing_status(target) == "PASS"


def test_output_must_stay_under_project_root(tmp_path):
    module = load_module()
    root = tmp_path / "project"
    root.mkdir()
    with pytest.raises(ValueError, match="outside project root"):
        module.require_under_root(tmp_path / "elsewhere", root)
