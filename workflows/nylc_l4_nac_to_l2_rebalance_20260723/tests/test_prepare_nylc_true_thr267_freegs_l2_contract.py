from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "prepare_nylc_true_thr267_freegs_l2.py"


def test_free_gs_builder_uses_correct_true_thr267_and_c18_graph():
    text = SCRIPT.read_text()
    assert '"thr_og1": 8961' in text
    assert '"carbonyl_c": 10297' in text
    assert '"carbonyl_o": 10298' in text
    assert '"amide_n": 10296' in text
    assert '"carbonyl_c": 25' in text
    assert '"carbonyl_o": 26' in text
    assert '"amide_n": 24' in text


def test_builder_uses_audited_l2_and_sha_locked_free_source():
    text = SCRIPT.read_text()
    assert "PA66_L2_GMX.itp" in text
    assert "b0e753c60fd4b71c282d21cc6106a15e73d91d12a20d80e92dd01516162eb301" in text
    assert "0f0c4932e591084feba457819daf1c0bf5cd6ec249778f0673005917b82c4c22" in text
    assert "prepare_one" in text
    assert "source_manifest.json" in text


def test_builder_refuses_overwrite_and_requires_rare_unbiased_status():
    text = SCRIPT.read_text()
    assert "VERIFIED_RARE_UNRESTRAINED_NAC_GS_CANDIDATE" in text
    assert "refusing to overwrite" in text
    assert "BUILD_PASS" in text
