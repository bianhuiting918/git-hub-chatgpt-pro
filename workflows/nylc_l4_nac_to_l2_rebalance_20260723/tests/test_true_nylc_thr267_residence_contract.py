from pathlib import Path

FLOW = Path(__file__).resolve().parents[1]


def _read(relative: str) -> str:
    return (FLOW / relative).read_text()


def test_true_thr267_residence_uses_identity_checked_atom():
    text = _read("scripts/audit_true_nylc_thr267_residence.py")
    assert "TRUE_THR267_OG1 = 8961" in text
    assert "WRONG_THR262_OG1 = 8896" in text
    assert 'expected=(267, "THR", "OG1")' in text
    assert 'expected=(262, "THR", "OG1")' in text


def test_true_thr267_residence_preserves_unrestrained_flags_and_common_gate():
    text = _read("scripts/audit_true_nylc_thr267_residence.py")
    assert "substrate_restrained" in text
    assert "gate_restrained" in text
    assert "distance_nm <= 0.35" in text
    assert "95.0 <= angle_deg <= 115.0" in text
    assert "longest_continuous_nac" in text


def test_corrective_array_covers_only_nylc_c18_and_c23():
    text = _read("slurm/run_true_nylc_thr267_residence_array.sbatch")
    assert "#SBATCH --array=0-1" in text
    assert "C18 10297 10298" in text
    assert "C23 10303 10304" in text
    assert "--workers 4" in text
    assert "SUPERSEDED_WRONG_NUCLEOPHILE" in text
