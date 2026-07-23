from pathlib import Path

FLOW = Path(__file__).resolve().parents[1]


def _read(relative: str) -> str:
    return (FLOW / relative).read_text()


def test_c18_free_1ns_extension_contract():
    text = _read("slurm/run_nylc_c18_step1_free_1ns.sbatch")
    assert "nylc_c18_11854ps" in text
    assert "nylc_c23_29684ps" not in text
    assert "61687591" in text
    assert "bab252dff6f2ce377e023c2cf1d2441da32fb0f2e6795622b782109cc0e26305" in text
    assert "0faa7e20ea5a4f571801c5bb4790777009556f18525df2cfdebf808f37363f57" in text


def test_c18_free_1ns_is_checkpoint_continuation_and_unrestrained():
    text = _read("slurm/run_nylc_c18_step1_free_1ns.sbatch")
    mdp = _read("mdp/npt300free.mdp")
    assert "npt300free.mdp" in text
    assert '-t "$input_cpt"' in text
    assert "-DPOSRES" not in mdp
    assert "nsteps                   = 500000" in mdp
    assert "dt                       = 0.002" in mdp


def test_c18_free_1ns_technical_completion_is_not_scientific_pass():
    text = _read("slurm/run_nylc_c18_step1_free_1ns.sbatch")
    assert "UNRESTRAINED_1NS_COMPLETE" in text
    assert "scientific_status" in text
    assert "NOT_EVALUATED_PENDING_TRAJECTORY_AUDIT" in text
    assert "FINAL_SCIENTIFIC_PASS" not in text
