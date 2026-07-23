from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
AUDIT = ROOT / "scripts" / "audit_em_double_result.py"
CONT = ROOT / "slurm" / "run_nylc_true_thr267_freegs_rebalance_after_em_double.sbatch"


def test_em_salvage_audits_immutable_completed_result():
    text = AUDIT.read_text()
    assert "Finished mdrun" in text
    assert "converged to Fmax < 500" in text
    assert "fmax_kj_mol_nm" in text
    assert "sha256" in text
    assert "scientific_nac_evidence" in text
    assert "PASS_TECHNICAL_EM" in text


def test_continuation_starts_from_passed_double_em_without_rerunning_em():
    text = CONT.read_text()
    assert "em_cg_flexible_double_retry3/PASS.json" in text
    assert "em_cg_flexible_double_retry3/run.gro" in text
    assert "run_stage em" not in text
    for stage in ("nvt50", "nvt150", "nvt300", "npt300r", "npt300rel", "npt300free"):
        assert stage in text


def test_continuation_free_window_is_unrestrained_and_scientifically_audited():
    text = CONT.read_text()
    assert 'run_stage npt300free' in text
    assert 'no' in text
    assert "npt300free_true_thr267_l2.mdp" in text
    assert "audit_nylc_true_thr267_l2_free.py" in text
    assert "atomnr 8961 plus atomnr 10287" in text
    assert "8961 10287 10288" in text
    assert "residues 261-266" in text
