from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SBATCH = ROOT / "slurm" / "run_true_thr267_unrestrained_pilot_array.sbatch"
AUDIT = ROOT / "scripts" / "audit_true_thr267_unrestrained_pilot.py"


def test_three_seed_pilot_is_fully_unrestrained():
    text = SBATCH.read_text()
    assert "#SBATCH --array=0-2" in text
    assert "pull                     = no" in text
    assert "define" not in text
    assert "-r " not in text
    assert "gen-vel                  = yes" in text
    assert "free_nvt20" in text
    assert "free_npt100" in text
    assert "Parrinello-Rahman" in text


def test_pilot_uses_true_thr267_identity_and_nac_gate():
    text = SBATCH.read_text()
    assert "atomnr 8961 plus atomnr 10297" in text
    assert "8961 10297 10298" in text
    audit = AUDIT.read_text()
    assert "0.35" in audit
    assert "95.0" in audit
    assert "115.0" in audit
    assert "minimum_heavy_atom_contact" in audit


def test_restrained_source_is_not_called_scientific_gs():
    text = SBATCH.read_text()
    assert "d0f157575193c4662f5463c14f5d9075e817c18142aa74f0f7456b8bc732b4f7" in text
    audit = AUDIT.read_text()
    assert "fully_unrestrained_NPT_100ps" in audit
    assert "PASS_UNRESTRAINED_PILOT_NAC_PRESENT" in audit
    assert "FAIL_UNRESTRAINED_PILOT_NO_NAC" in audit
