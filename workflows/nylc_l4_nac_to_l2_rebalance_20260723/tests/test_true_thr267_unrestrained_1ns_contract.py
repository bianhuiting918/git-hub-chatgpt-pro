from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SBATCH = ROOT / "slurm" / "run_true_thr267_unrestrained_1ns.sbatch"
AUDIT = ROOT / "scripts" / "audit_true_thr267_unrestrained_pilot.py"


def test_one_ns_is_checkpoint_continuation_and_unrestrained():
    text = SBATCH.read_text()
    assert "nsteps                   = 500000" in text
    assert "continuation             = yes" in text
    assert "pull                     = no" in text
    assert "define" not in text
    assert "-r " not in text
    assert "free_npt100.cpt" in text
    assert "free_npt100.gro" in text


def test_parent_is_sha_and_scientific_gate_checked():
    text = SBATCH.read_text()
    assert "ab6fb4083f02e9bf1c2af3677500ceccc139ba08882e32750c63136bf3320b6f" in text
    assert "0e8a93774b8c0836bee7a36b041676112c014cd1eec9e96055a2ade971b6d459" in text
    assert "PASS_UNRESTRAINED_PILOT_NAC_PRESENT" in text
    assert "grompp" in text
    assert "-maxwarn 0" in text


def test_one_ns_audit_can_use_generic_stem_and_window_label():
    text = AUDIT.read_text()
    assert '--stem' in text
    assert '--window-type' in text
    assert '--output-name' in text
    sbatch = SBATCH.read_text()
    assert "fully_unrestrained_NPT_1ns" in sbatch
    assert "unrestrained_1ns_audit.json" in sbatch
