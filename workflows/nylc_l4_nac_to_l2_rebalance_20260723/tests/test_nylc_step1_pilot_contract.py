from pathlib import Path

FLOW = Path(__file__).resolve().parents[1]


def _read(relative: str) -> str:
    return (FLOW / relative).read_text()


def test_nylc_pilot_is_limited_to_two_corrected_candidates():
    text = _read("slurm/run_nylc_step1_pilot_array.sbatch")
    assert "#SBATCH --array=0-1%2" in text
    assert "nylc_c18_11854ps" in text
    assert "nylc_c23_29684ps" in text
    assert "nyl50" not in text.lower()
    assert "nyl12" not in text.lower()
    assert "SUPERSEDED_WRONG_NUCLEOPHILE_IDENTITY" in text


def test_contact_relaxation_does_not_require_em_convergence():
    text = _read("slurm/run_nylc_step1_pilot_array.sbatch")
    assert "maximum_force" in text
    assert "EM_CONVERGENCE_NOT_REQUIRED" in text
    assert "value <= 500" not in text
    assert "fmax <= 500" not in text.lower()


def test_low_temperature_pilot_uses_small_timesteps_and_heavy_atom_restraints():
    em = _read("mdp/em_nac_contact.mdp")
    nvt10 = _read("mdp/nvt10_nac.mdp")
    nvt50 = _read("mdp/nvt50_nac.mdp")
    assert "-DPOSRES" in em and "-DPOSRES_L2_1000" in em
    assert "integrator               = steep" in em
    assert "dt                       = 0.0005" in nvt10
    assert "ref-t                    = 10" in nvt10
    assert "-DPOSRES_WEAK" in nvt10 and "-DPOSRES_L2_500" in nvt10
    assert "dt                       = 0.001" in nvt50
    assert "ref-t                    = 50" in nvt50
    assert "-DPOSRES_WEAK" in nvt50 and "-DPOSRES_L2_100" in nvt50


def test_pilot_cannot_claim_scientific_nac_pass():
    text = _read("slurm/run_nylc_step1_pilot_array.sbatch")
    assert "TECHNICAL_PILOT_PASS" in text
    assert "scientific_nac_evidence" in text
    assert '"scientific_nac_evidence": False' in text
    assert "SCIENTIFIC_PASS" not in text
