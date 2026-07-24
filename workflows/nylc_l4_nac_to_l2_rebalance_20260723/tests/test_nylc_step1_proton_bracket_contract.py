from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREP = ROOT / "scripts" / "prepare_nylc_step1_proton_bracket.py"
AUDIT = ROOT / "scripts" / "audit_nylc_step1_proton_stage.py"
SBATCH = ROOT / "slurm" / "run_nylc_step1_proton_bracket.sbatch"


def test_proton_bracket_locks_passed_attack_seed_and_two_acceptors():
    assert PREP.exists(), "proton bracket preparer is not implemented"
    text = PREP.read_text()
    assert "d5e607dd722be4fd6e5d33016afb4975d18b200a78d1fceedf5bfb9fed1b87bf" in text
    assert 'choices=("OD1", "OD2")' in text
    assert "ASP306_OD1 = 9572" in text and "ASP306_OD2 = 9573" in text
    assert "PASS_ATTACK_PREREQUISITE" in text


def test_proton_bracket_is_gradual_and_does_not_restrain_og_h():
    text = PREP.read_text() if PREP.exists() else ""
    for stage in ('"p00"', '"p01"', '"p02"', '"p03"'):
        assert stage in text
    assert "THR267_HG1, acceptor_atom" in text
    assert "THR267_OG1, THR267_HG1" not in text
    assert "dftb_telec=200.0" in text
    assert "ntr=1" in text and "!@H=" in text


def test_stage_audit_rejects_detached_proton_and_separates_endpoint():
    assert AUDIT.exists(), "proton stage auditor is not implemented"
    text = AUDIT.read_text()
    assert "DETACHED_MAX_A = 1.35" in text
    assert "min(ogh, hacc) <= DETACHED_MAX_A" in text
    assert "PASS_PROTON_BRACKET_STAGE" in text
    assert "PASS_SCIENTIFIC_TETRAHEDRAL_SEED_REACHED" in text
    assert "FAIL_SCIENTIFIC_PROTON_BRACKET_STAGE" in text
    assert "not a TS, committor, PMF, or barrier" in text


def test_wrapper_is_independent_scnet_branch_with_stage_gates():
    assert SBATCH.exists(), "proton bracket Slurm wrapper is not implemented"
    text = SBATCH.read_text()
    assert "#SBATCH -p xahcnormal" in text
    assert "#SBATCH --mem-per-cpu=2500M" in text
    assert 'ACCEPTOR="${ACCEPTOR:?set ACCEPTOR=OD1_or_OD2}"' in text
    assert "for stage in p00 p01 p02 p03" in text
    assert "audit_nylc_step1_proton_stage.py" in text
    assert "SCIENTIFIC_FAIL.json" in text and "tetrahedral_seed.rst7" in text
    assert "run_history.tsv" in text and "run_history.jsonl" in text
    for token in ("FINAL RESULTS", "5.  TIMINGS", "SCC is not converged", "SANDER BOMB", "FATAL", "NaN"):
        assert token in text
    assert "sander" in text and "mdrun" not in text
