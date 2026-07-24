from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SBATCH = ROOT / "slurm" / "run_nylc_step1_attack_continuation.sbatch"


def test_attack_continuation_locks_identical_valid_seed():
    assert SBATCH.exists(), "attack continuation wrapper is not implemented"
    text = SBATCH.read_text()
    assert "step1_constrained_bracket_OD1_job_61723169/run/w01.rst7" in text
    assert "10d9b662cdd0bb2e945b34b798e4688373e43e1cca5737b8c86f7658c1975f38" in text
    assert "sha256sum" in text


def test_attack_continuation_is_attack_only_and_gated_each_round():
    text = SBATCH.read_text() if SBATCH.exists() else ""
    assert "for round in r01 r02 r03" in text
    assert "w01.in" in text and "w01.RST" in text
    assert "w02.in" not in text and "w03.in" not in text
    assert "gate_nylc_step1_attack_prerequisite.py" in text
    assert "PASS_ATTACK_PREREQUISITE" in text
    assert "FAIL_SCIENTIFIC_ATTACK_PREREQUISITE" in text
    assert "attack_seed.rst7" in text


def test_attack_continuation_is_scnet_audited_and_no_trajectory():
    text = SBATCH.read_text() if SBATCH.exists() else ""
    assert "#SBATCH -p xahcnormal" in text
    assert "#SBATCH --mem-per-cpu=2500M" in text
    assert "run_history.tsv" in text and "run_history.jsonl" in text
    assert "FINAL RESULTS" in text and "5.  TIMINGS" in text
    for token in ("SCC is not converged", "SANDER BOMB", "FATAL", "NaN"):
        assert token in text
    assert "sander" in text and "mdrun" not in text
    assert ".xtc" not in text and ".trr" not in text and ".nc" not in text
