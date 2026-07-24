from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREPARE = ROOT / "scripts" / "prepare_nylc_step1_ash306_qm_smoke.py"
AUDIT = ROOT / "scripts" / "audit_nylc_step1_ash306_qm_smoke.py"
SBATCH = ROOT / "slurm" / "run_nylc_step1_ash306_qm_smoke.sbatch"


def test_preparer_locks_authoritative_full_system_and_shifted_indices():
    assert PREPARE.exists(), "ASH306 QM smoke preparer is not implemented"
    text = PREPARE.read_text()
    assert "bf6fa73b819da92c1cb0fdc6f97c17bd3db37e1a49c0127ecac39fd9a8d1b8ae" in text
    assert "80a8c9371cf59799c25927c5d1da95e6ecbb4662ba3615f547fa752cf825289d" in text
    assert "THR267_OG1 = 8961" in text
    assert "ASH306_OD1 = 9572" in text
    assert "ASP308_OD1 = 9592" in text
    assert "L2_FIRST = 10274" in text and "L2_LAST = 10352" in text
    assert "WATER_THR_ASN_1 = 50166" in text
    assert "WATER_THR_ASN_2 = 51303" in text
    assert "WATER_ASP_PAIR = 81219" in text


def test_preparer_uses_charge_consistent_even_electron_regions():
    text = PREPARE.read_text() if PREPARE.exists() else ""
    assert '"core": {"qmcharge": 0, "atom_count": 111, "boundary_count": 3, "electrons": 388}' in text
    assert '"network": {"qmcharge": 1, "atom_count": 162, "boundary_count": 7, "electrons": 570}' in text
    assert '"HD2"' in text
    assert "Asp306H_Asp308-" in text
    assert "absolute energy" in text


def test_wrapper_runs_two_sander_smokes_on_cpu_and_records_history():
    assert SBATCH.exists(), "ASH306 QM smoke wrapper is not implemented"
    text = SBATCH.read_text()
    assert "#SBATCH -p xahcnormal" in text and "#SBATCH -c 4" in text
    assert "ash306_full_system_preflight_job_61718715" in text
    assert text.count("sander -O") == 2
    assert "mdrun" not in text
    assert "run_history.tsv" in text and "run_history.jsonl" in text
    assert "trap on_error ERR" in text
    assert "PASS_ASH306_QM_REGION_SMOKE" in text


def test_auditor_hard_fails_scc_and_runtime_errors_without_ts_claim():
    assert AUDIT.exists(), "ASH306 QM smoke auditor is not implemented"
    text = AUDIT.read_text()
    for pattern in ("Convergence could not be achieved", "SANDER BOMB", "NaN", "FATAL"):
        assert pattern in text
    assert "READY_FOR_ASH306_QM_REGION_SMOKE" in text
    assert "PASS_ASH306_QM_REGION_SMOKE" in text
    assert "not a TS, RC, PMF, or barrier" in text
