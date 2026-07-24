from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREP = ROOT / "scripts" / "prepare_nylc_step1_constrained_bracket.py"
AUDIT = ROOT / "scripts" / "audit_nylc_step1_constrained_bracket.py"
SBATCH = ROOT / "slurm" / "run_nylc_step1_constrained_bracket.sbatch"


def test_bracket_uses_authoritative_deprotonated_explicit_water_source():
    assert PREP.exists(), "Step1 constrained-bracket preparer is not implemented"
    text = PREP.read_text()
    assert "step1_dftb3_preflight_post61710861_job61712692" in text
    assert "selected_free_l2_nac.rst7" in text
    assert "WATER_ASP_PAIR = (81218, 81219, 81220)" in text
    assert "EXPECTED_QM_ATOMS = 110" in text
    assert "QMCHARGE = -1" in text


def test_bracket_keeps_od1_and_od2_as_independent_hypotheses():
    text = PREP.read_text() if PREP.exists() else ""
    assert 'choices=("OD1", "OD2")' in text
    assert "ASP306_OD1 = 9572" in text
    assert "ASP306_OD2 = 9573" in text
    assert "THR267_OG1 = 8961" in text
    assert "THR267_HG1 = 8962" in text
    assert "L2_C12 = 10287" in text


def test_bracket_is_gradual_and_does_not_overdetermine_og_h_breaking():
    text = PREP.read_text() if PREP.exists() else ""
    for token in ('"w00"', '"w01"', '"w02"', '"w03"'):
        assert token in text
    assert "iat={THR267_OG1},{L2_C12}" in text
    assert "iat={THR267_HG1},{acceptor_atom}" in text
    assert "iat={THR267_OG1},{THR267_HG1}" not in text
    assert "dftb_telec=200.0" in text
    assert text.count('"proton_A": None') == 2
    assert 'if window["proton_A"] is not None:' in text
    assert "ntr=1" in text and "restraint_wt=1.0" in text
    assert "!@H=" in text


def test_audit_separates_technical_and_scientific_status():
    assert AUDIT.exists(), "Step1 bracket auditor is not implemented"
    text = AUDIT.read_text()
    assert "PASS_TECHNICAL_STEP1_BRACKET" in text
    assert "PASS_SCIENTIFIC_TETRAHEDRAL_SEED_REACHED" in text
    assert "FAIL_SCIENTIFIC_TETRAHEDRAL_SEED_NOT_REACHED" in text
    assert "not a TS, committor, PMF, or barrier" in text
    for token in ("FINAL RESULTS", "Run done", "SANDER BOMB", "NaN", "FATAL"):
        assert token in text


def test_wrapper_is_scnet_cpu_and_failure_isolated_by_acceptor():
    assert SBATCH.exists(), "Step1 bracket Slurm wrapper is not implemented"
    text = SBATCH.read_text()
    assert "#SBATCH -p xahcnormal" in text
    assert "#SBATCH --mem-per-cpu=2500M" in text
    assert 'ACCEPTOR="${ACCEPTOR:?set ACCEPTOR=OD1_or_OD2}"' in text
    assert "run_history.tsv" in text and "run_history.jsonl" in text
    assert "sander" in text
    assert '-ref "$current"' in text
    assert "mdrun" not in text
