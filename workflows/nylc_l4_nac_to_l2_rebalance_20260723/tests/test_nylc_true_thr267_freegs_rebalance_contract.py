from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SBATCH = ROOT / "slurm" / "run_nylc_true_thr267_freegs_rebalance.sbatch"
AUDIT = ROOT / "scripts" / "audit_nylc_true_thr267_l2_free.py"


def test_rebalance_requires_build_preflight_and_runs_stages_in_order():
    text = SBATCH.read_text()
    assert "PASS_BUILD_AND_GROMPP" in text
    order = ["run_stage em", "run_stage nvt50", "run_stage nvt150", "run_stage nvt300", "run_stage npt300r", "run_stage npt300rel", "run_stage npt300free"]
    positions = [text.index(item) for item in order]
    assert positions == sorted(positions)


def test_final_window_is_one_ns_fully_unrestrained():
    text = SBATCH.read_text()
    assert 'run_stage npt300free' in text
    assert 'restrained=no' in text
    assert "npt300free.mdp" in text
    assert "scientific audit" in text.lower()


def test_l2_audit_uses_corrected_atoms_and_gate():
    text = SBATCH.read_text()
    assert "atomnr 8961 plus atomnr 10287" in text
    assert "8961 10287 10288" in text
    audit = AUDIT.read_text()
    assert "NylC residues 261-266; Thr267 excluded" in audit
    assert "PASS_UNRESTRAINED_L2_NAC_PRESENT" in audit
    assert "FAIL_UNRESTRAINED_L2_NO_NAC" in audit
    assert "minimum_heavy_atom_contact" in audit
