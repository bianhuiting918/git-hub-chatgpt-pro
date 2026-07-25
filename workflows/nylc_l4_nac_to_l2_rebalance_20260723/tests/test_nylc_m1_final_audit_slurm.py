import pathlib


FLOW = pathlib.Path(__file__).resolve().parents[1]
SCRIPT = FLOW / "slurm" / "run_nylc_m1_ensemble_final_audit.sbatch"


def test_final_audit_job_uses_dynamic_stage_b_denominator_and_raw_primitives():
    text = SCRIPT.read_text()

    assert 'STAGEB_ARRAY_JOB_ID="${STAGEB_ARRAY_JOB_ID:?' in text
    assert 'STAGEB_COUNT="${STAGEB_COUNT:?' in text
    assert 'for task_id in $(seq 0 $((STAGEB_COUNT - 1)))' in text
    assert "generate_nylc_m1_ensemble_primitives.py" in text
    assert "assemble_nylc_m1_final_audit_input.py" in text
    assert "independent_audit_nylc_m1_ensemble.py" in text
    assert "stageB_records_manifest.json" in text
    assert "final_audit_input.json" in text


def test_final_audit_job_preserves_failures_and_does_not_make_science_exit_nonzero():
    text = SCRIPT.read_text()

    assert "NOT_EVALUATED_STAGEB_MISSING_SLOT" in text
    assert "NOT_EVALUATED_STAGEB_INCOMPLETE" in text
    assert "FAIL_TECHNICAL_AUDIT" in text
    assert 'scientific_status="$(' in text
    assert "exit 0" in text


def test_final_audit_job_uses_only_unrestrained_stage_b_outputs():
    text = SCRIPT.read_text()

    assert "STAGE_B_COMPLETE.json" in text
    assert "free_tpr_contract.json" in text
    assert "position_restraints" not in text
    assert "stageA_job_61814751" not in text
    assert "61801874" in text
