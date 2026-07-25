from pathlib import Path

FLOW = Path(__file__).resolve().parents[1]
SLURM = FLOW / "slurm"


def assert_environment_source_is_nounset_safe(text: str, source_line: str) -> None:
    source_index = text.index(source_line)
    assert text.rfind("set +u", 0, source_index) >= 0
    assert text.find("set -u", source_index + len(source_line)) >= 0


def test_em_array_uses_double_precision_flexible_water_and_independent_failures():
    text = (SLURM / "run_nylc_m1_ensemble_em_array.sbatch").read_text()

    assert "#SBATCH --array=0-11%4" in text
    assert "gromacs-2023.1-cp2k-2024.1/bin/gmx_mpi_d" in text
    assert "em_cg_flexible_m1.mdp" in text
    assert "mpirun -np 1" in text
    assert "fmax" in text and "500.0" in text
    assert "NOT_EVALUATED_BUILD_FAIL" in text
    assert_environment_source_is_nounset_safe(
        text,
        "source /work/home/acshdt1dks/opt/gmx-cp2k/env.sh",
    )


def test_stage_a_array_maps_candidate_and_seed_and_runs_one_dcu():
    text = (SLURM / "run_nylc_m1_ensemble_stageA_array.sbatch").read_text()

    assert "#SBATCH --array=0-35%8" in text
    assert "#SBATCH --gres=dcu:1" in text
    assert "candidate_index=$((SLURM_ARRAY_TASK_ID / 3))" in text
    assert "seed_index=$((SLURM_ARRAY_TASK_ID % 3))" in text
    assert "seeds=(26711 26723 26737)" in text
    for stage in ("nvt50", "nvt150", "nvt300", "npt300r", "npt300rel"):
        assert stage in text
    assert "npt300free_m1_stageA.mdp" in text
    assert '"position_restraints": 0' in text
    assert '"distance_restraints": 0' in text
    assert '"scientific_window": "fully_unrestrained_NPT_100ps"' in text
    assert_environment_source_is_nounset_safe(
        text,
        "source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh",
    )
