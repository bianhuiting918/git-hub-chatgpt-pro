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


def test_stage_a_rank_job_audits_complete_universe_afterany():
    text = (SLURM / "run_nylc_m1_ensemble_rank_stageA.sbatch").read_text()

    assert "#SBATCH -p xahcnormal" in text
    assert "STAGEA_ARRAY_JOB_ID" in text
    assert "for task_id in $(seq 0 35)" in text
    assert "STAGE_A_COMPLETE.json" in text
    assert "NOT_EVALUATED.json" in text
    assert "generate_nylc_m1_ensemble_primitives.py" in text
    assert "audit_nylc_m1_ensemble_replica.py" in text
    assert "rank_nylc_m1_ensemble.py" in text
    assert "expected_stageA_slots=36" in text
    assert "dependency=afterany" in text
    assert "afterok" not in text


def test_free_tpr_contract_only_rejects_nonempty_define_values():
    for name in (
        "run_nylc_m1_ensemble_stageA_array.sbatch",
        "run_nylc_m1_ensemble_stageA_resume_array.sbatch",
    ):
        text = (SLURM / name).read_text()
        assert "define_values = []" in text
        assert 'key.strip().lower() == "define" and value.strip()' in text
        assert "define_present = bool(define_values)" in text


def test_stage_a_resume_uses_immutable_parent_checkpoints_and_new_outputs():
    text = (SLURM / "run_nylc_m1_ensemble_stageA_resume_array.sbatch").read_text()

    assert "#SBATCH --array=0-35%8" in text
    assert "#SBATCH --gres=dcu:1" in text
    assert "SOURCE_STAGEA_ARRAY_JOB_ID" in text
    assert "stageA_resume_job_" in text
    assert '[[ -s "$SOURCE_ROOT/$stage/PASS.json" ]]' in text
    assert 'ln -s "$SOURCE_ROOT/$stage" "$RUN_ROOT/$stage"' in text
    assert "npt300free_m1_stageA.mdp" in text
    assert '"fully_unrestrained": True' in text
    assert "refuse_overwrite=yes" in text
