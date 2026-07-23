from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "select_true_thr267_unrestrained_gs.sh"


def test_selector_uses_only_unrestrained_nac_frame():
    text = SCRIPT.read_text()
    assert "job_61705307/seed_26703" in text
    assert "lowest_potential_nac_frame" in text
    assert "PASS_UNRESTRAINED_PILOT_NAC_PRESENT" in text
    assert "trjconv" in text
    assert "-dump" in text


def test_selector_preserves_rare_nac_caveat_and_one_ns_result():
    text = SCRIPT.read_text()
    assert "job_61705692/unrestrained_1ns_audit.json" in text
    assert "VERIFIED_RARE_UNRESTRAINED_NAC_GS_CANDIDATE" in text
    assert "NOT_ENSEMBLE_STABLE" in text
    assert "sha256sum" in text
    assert "minimum_ligand_protein_contact" in text


def test_selector_refuses_overwrite():
    text = SCRIPT.read_text()
    assert 'test ! -e "$output_root/source.gro"' in text
    assert 'test ! -e "$output_root/source_manifest.json"' in text
