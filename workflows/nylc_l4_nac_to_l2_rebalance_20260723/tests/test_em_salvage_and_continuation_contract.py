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


def test_em_audit_parses_real_gromacs_force_and_energy_lines():
    import importlib.util
    spec = importlib.util.spec_from_file_location("audit_em_double_result", AUDIT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    sample = """
Potential Energy  = -2.30772586301368e+06
Maximum force     =  4.19978114876091e+02 on atom 3668
"""
    assert module.last_float(r"Maximum force\s*=\s*([+0-9.eE-]+)", sample) == 419.978114876091
    assert module.last_float(r"Potential Energy\s*=\s*([+0-9.eE-]+)", sample) == -2307725.86301368


def test_em_audit_main_accepts_completed_double_cg_log(tmp_path, monkeypatch):
    import importlib.util
    import pytest
    spec = importlib.util.spec_from_file_location("audit_em_double_result_main", AUDIT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    (tmp_path / "run.log").write_text(
        "GROMACS version:    2023.1\n"
        "Precision:          double\n"
        "Polak-Ribiere Conjugate Gradients converged to Fmax < 500 in 153 steps\n"
        "Potential Energy  = -2.30772586301368e+06\n"
        "Maximum force     =  4.19978114876091e+02 on atom 3668\n"
        "Finished mdrun on rank 0\n"
    )
    (tmp_path / "run.processed.mdp").write_text(
        "define = -DPOSRES_L2_1000 -DFLEXIBLE\nintegrator = cg\n"
    )
    (tmp_path / "run.gro").write_text("test\n133589\n")
    (tmp_path / "run.tpr").write_bytes(b"tpr")
    monkeypatch.setattr("sys.argv", ["audit", "--run-root", str(tmp_path)])
    with pytest.raises(SystemExit) as exc:
        module.main()
    assert exc.value.code == 0
    assert (tmp_path / "PASS.json").is_file()
