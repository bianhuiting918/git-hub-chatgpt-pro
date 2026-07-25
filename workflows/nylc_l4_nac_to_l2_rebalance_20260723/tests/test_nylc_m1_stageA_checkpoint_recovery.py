import pathlib


FLOW = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = FLOW / "slurm" / "run_nylc_m1_stageA_checkpoint_recovery.sbatch"


def test_recovery_reuses_only_audited_finished_nvt50_checkpoint():
    text = SCRIPT.read_text()

    assert 'RECOVERY_TASK_ID="${RECOVERY_TASK_ID:?' in text
    assert 'FAILED_STAGEA_ARRAY_JOB_ID="${FAILED_STAGEA_ARRAY_JOB_ID:?' in text
    assert 'grep -q "Finished mdrun"' in text
    assert "nvt50/RECOVERY_AUDIT.json" in text
    assert '"source_exit_code": 132' in text
    assert "run_stage nvt50" not in text
    assert 'run_stage nvt150' in text


def test_recovery_never_overwrites_completed_or_existing_later_stages():
    text = SCRIPT.read_text()

    assert '[[ ! -e "$RUN_ROOT/STAGE_A_COMPLETE.json" ]]' in text
    assert 'for stage in nvt150 nvt300 npt300r npt300rel npt300free_stageA' in text
    assert "refuse_existing_later_stage" in text
    assert "STAGE_A_COMPLETE.json" in text


def test_recovery_keeps_science_not_evaluated_until_free_window_audit():
    text = SCRIPT.read_text()

    assert '"fully_unrestrained": True' in text
    assert '"scientific_status": "NOT_EVALUATED_PENDING_STAGE_A_AUDIT"' in text
    assert "position_restraints != 0" in text
    assert "distance_restraints != 0" in text
