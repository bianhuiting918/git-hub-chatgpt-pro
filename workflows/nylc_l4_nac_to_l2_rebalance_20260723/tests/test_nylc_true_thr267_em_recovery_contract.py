from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SBATCH = ROOT / "slurm" / "run_nylc_true_thr267_freegs_em_lbfgs.sbatch"
MDP = ROOT / "mdp" / "em_lbfgs.mdp"


def test_recovery_uses_flexible_water_and_original_rebuilt_start():
    assert "-DFLEXIBLE" in MDP.read_text()
    text = SBATCH.read_text()
    assert 'start="$build/rebuilt.gro"' in text
    assert "runs/em/run.gro" not in text
    assert "em_lbfgs_true" in text


def test_recovery_keeps_strict_force_and_warning_gate():
    text = SBATCH.read_text()
    assert "value<=500.0" in text
    assert "LINCS WARNING" in text
    assert "SETTLE" in text
    assert "NaN" in text
    assert "-maxwarn 0" in text
