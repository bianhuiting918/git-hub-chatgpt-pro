from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SBATCH = ROOT / "slurm" / "run_nylc_true_thr267_freegs_em_cg_flexible.sbatch"
MDP = ROOT / "mdp" / "em_cg_flexible.mdp"


def test_recovery_uses_flexible_water_cg_and_original_start():
    mdp = MDP.read_text()
    assert "-DFLEXIBLE" in mdp
    assert "integrator               = cg" in mdp
    text = SBATCH.read_text()
    assert 'start="$build/rebuilt.gro"' in text
    assert "runs/em/run.gro" not in text
    assert "em_cg_flexible" in text


def test_recovery_keeps_strict_preflight_and_force_gate():
    text = SBATCH.read_text()
    assert "-maxwarn 0" in text
    assert "value<=500.0" in text
    assert "LINCS WARNING" in text
    assert "SETTLE" in text
    assert "NaN" in text
