from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PREP = ROOT / "scripts" / "prepare_nylc_step1_production_qm_smoke.py"
AUDIT = ROOT / "scripts" / "audit_nylc_step1_production_qm_smoke.py"
SBATCH = ROOT / "slurm" / "run_nylc_step1_production_qm_smoke.sbatch"
REAUDIT = ROOT / "slurm" / "run_nylc_step1_production_qm_smoke_reaudit.sbatch"


def test_preparer_defines_audited_core_and_network_regions():
    text = PREP.read_text()
    for atom in (7156, 7768, 8240, 8961, 9572, 9591, 10273, 10351, 50165, 51302, 81218):
        assert str(atom) in text
    assert '"core": {"qmcharge": -1, "atom_count": 110, "boundary_count": 3' in text
    assert '"network": {"qmcharge": 0, "atom_count": 161, "boundary_count": 7' in text
    assert "electron_count_including_link_h" in text
    assert "all_deprotonated_Asp306_Asp308" in text
    assert "DFTB3" in text and "dftb_telec=200.0" in text


def test_smoke_audit_is_numerical_only_and_rejects_scc_warnings():
    text = AUDIT.read_text()
    assert "Convergence could not be achieved" in text
    assert "SANDER BOMB" in text
    assert r're.search(r"Run\s+done"' in text
    assert "FINAL RESULTS" in text
    assert "PASS_PRODUCTION_QM_REGION_SMOKE" in text
    assert "not a TS, RC, PMF, or barrier" in text


def test_wrapper_runs_both_regions_on_scnet_cpu_and_records_history():
    text = SBATCH.read_text()
    assert "#SBATCH -p xahcnormal" in text
    assert "--gres=dcu" not in text
    assert "step1_dftb3_preflight_post61710861_job61712692/system.prmtop" in text
    assert "core_one_step.in" in text and "network_one_step.in" in text
    assert text.count("sander -O") == 2
    assert "step1_production_qm_smoke_job_" in text
    assert "run_history.tsv" in text and "run_history.jsonl" in text
    assert "trap on_error ERR" in text


def test_reaudit_reuses_61716670_outputs_without_rerunning_sander():
    text = REAUDIT.read_text()
    assert "#SBATCH -p xahcnormal" in text
    assert "step1_production_qm_smoke_job_61716670" in text
    assert "audit_nylc_step1_production_qm_smoke.py" in text
    assert "sander -O" not in text
    assert "step1_production_qm_smoke_reaudit_job_" in text
    assert "run_history.tsv" in text and "run_history.jsonl" in text
