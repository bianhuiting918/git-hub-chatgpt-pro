from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "select_true_thr267_low_energy_nac.sh"


def test_selector_uses_true_nac_geometry_and_potential_energy():
    text = SCRIPT.read_text()
    assert "distance-max-nm 0.35" in text
    assert "angle-min-deg 95" in text
    assert "angle-max-deg 115" in text
    assert "Potential" in text
    assert "lowest_potential_nac_frame" in text


def test_selector_extracts_audited_frame_without_overwriting_parent():
    text = SCRIPT.read_text()
    assert "trjconv" in text
    assert "-dump" in text
    assert "sha256sum" in text
    assert "selected_lowest_potential_nac.json" in text
    assert "NOT_A_SCIENTIFIC_GS_RESTRAINED_SOURCE" in text
    assert "release2.gro" not in text.split("trjconv", 1)[-1]
