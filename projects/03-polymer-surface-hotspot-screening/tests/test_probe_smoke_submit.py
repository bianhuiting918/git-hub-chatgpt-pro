from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "submit_probe_match_smoke.py"
SPEC = importlib.util.spec_from_file_location("submit_probe_match_smoke", SCRIPT)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot load {SCRIPT}")
mod = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mod)


def make_row(root: Path, family: str = "PET") -> dict:
    record = "PET_6ILW_A" if family == "PET" else "NYLON_3AXG_A"
    probe = "PET_DMT" if family == "PET" else "NYL_NMA"
    shell = root / "results" / "shell" / record
    patch = root / "results" / "patch" / record
    maps = root / "results" / "maps" / record
    receptor = root / "results" / "receptor" / record / "receptor.pdbqt"
    for gate in [shell / "SHELL_PASS.json", patch / "PATCH_PASS.json", maps / "MAPS_PASS.json"]:
        gate.parent.mkdir(parents=True, exist_ok=True)
        gate.write_text("{}\n", encoding="utf-8")
    receptor.parent.mkdir(parents=True, exist_ok=True)
    receptor.write_text("ATOM\n", encoding="ascii")
    return {
        "material_family": family,
        "record_id": record,
        "probe_ids": probe,
        "shell_dir": str(shell),
        "patch_dir": str(patch),
        "maps_root": str(maps),
        "receptor_pdbqt": str(receptor),
        "output_dir": str(root / "results" / "probe_match_smoke_test" / record),
    }


def test_validate_row_accepts_complete_pet_record(tmp_path):
    row = make_row(tmp_path, "PET")
    validated = mod.validate_row(row, tmp_path)
    assert validated["material_family"] == "PET"
    assert validated["probe_ids"] == "PET_DMT"


def test_validate_row_rejects_output_outside_allowed_root(tmp_path):
    row = make_row(tmp_path, "NYLON")
    row["output_dir"] = str(tmp_path.parent / "escape")
    with pytest.raises(ValueError, match="output"):
        mod.validate_row(row, tmp_path)


def test_validate_row_rejects_family_probe_mismatch(tmp_path):
    row = make_row(tmp_path, "PET")
    row["probe_ids"] = "NYL_NMA"
    with pytest.raises(ValueError, match="probe"):
        mod.validate_row(row, tmp_path)


def test_build_sbatch_arguments_freezes_partition_and_memory(tmp_path):
    row = mod.validate_row(make_row(tmp_path, "PET"), tmp_path)
    arguments = mod.build_sbatch_arguments(
        row,
        root=tmp_path,
        commit="a" * 40,
        wrapper=tmp_path / "wrapper.sh",
        log_dir=tmp_path / "logs",
    )
    assert arguments[0] == "sbatch"
    assert "xahcnormal" in arguments
    assert "3G" in arguments
    assert "--wrap" in arguments
