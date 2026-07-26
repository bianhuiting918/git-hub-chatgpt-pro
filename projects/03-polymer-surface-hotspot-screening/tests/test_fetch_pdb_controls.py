import gzip
import importlib.util
from pathlib import Path

import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "fetch_pdb_controls.py"


def load_module():
    spec = importlib.util.spec_from_file_location("fetch_pdb_controls", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_frozen_sources_are_https_rcsb():
    module = load_module()
    assert set(module.FROZEN_SOURCES) == {"1WYB", "9DYS", "3AXG_BIOASSEMBLY1"}
    assert all(source.url.startswith("https://files.rcsb.org/") for source in module.FROZEN_SOURCES.values())


def test_validate_pdb_requires_coordinates_and_expected_id():
    module = load_module()
    valid = (
        "HEADER    HYDROLASE                         09-FEB-05   1WYB\n"
        "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N\n"
        "END\n"
    ).encode("ascii")
    assert module.validate_pdb_bytes(valid, "1WYB")["n_coordinate_records"] == 1
    with pytest.raises(ValueError, match="expected PDB id"):
        module.validate_pdb_bytes(valid, "9DYS")
    with pytest.raises(ValueError, match="coordinate"):
        module.validate_pdb_bytes(b"HEADER 1WYB\nEND\n", "1WYB")


def test_decode_source_supports_gzip():
    module = load_module()
    raw = (
        "HEADER    HYDROLASE                         09-FEB-05   3AXG\n"
        "ATOM      1  N   ALA A   1      10.000  10.000  10.000  1.00 20.00           N\n"
        "END\n"
    ).encode("ascii")
    assert module.decode_download(gzip.compress(raw), compressed=True) == raw


def test_existing_output_directory_is_never_overwritten(tmp_path):
    module = load_module()
    output = tmp_path / "controls"
    output.mkdir()
    (output / "sentinel").write_text("keep", encoding="utf-8")
    with pytest.raises(FileExistsError):
        module.prepare_output_directory(output)
    assert (output / "sentinel").read_text(encoding="utf-8") == "keep"
