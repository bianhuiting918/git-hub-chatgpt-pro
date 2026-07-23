from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SBATCH = ROOT / "slurm" / "run_nylc_true_thr267_freegs_em_cg_double.sbatch"


def test_double_recovery_uses_validated_double_binary_and_env():
    text = SBATCH.read_text()
    assert "source /work/home/acshdt1dks/opt/gmx-cp2k/env.sh" in text
    assert "gmx_mpi_d" in text
    assert "2023.1" in text
    assert "em_cg_flexible.mdp" in text
    assert "set +u" in text
    assert "em_cg_flexible_double_retry3" in text
    assert 'start="$build/rebuilt.gro"' in text


def test_double_recovery_is_cpu_and_strict():
    text = SBATCH.read_text()
    assert "-nb cpu" in text
    assert "-pme cpu" in text
    assert "-bonded cpu" in text
    assert "-maxwarn 0" in text
    assert "value<=500.0" in text
    assert "hard_warning_count" in text
