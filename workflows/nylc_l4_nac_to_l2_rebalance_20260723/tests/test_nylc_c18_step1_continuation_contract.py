from pathlib import Path

FLOW = Path(__file__).resolve().parents[1]


def _read(relative: str) -> str:
    return (FLOW / relative).read_text()


def test_continuation_is_c18_only_and_checks_selected_frame_sha():
    text = _read("slurm/run_nylc_c18_step1_continuation.sbatch")
    assert "nylc_c18_11854ps" in text
    assert "nylc_c23_29684ps" not in text
    assert "c18_50k_late_lowest_potential_nac_9p55ps.gro" in text
    assert "c67cbb1d275863606be62628df6829e8ef15fbad39f3b5a51188c311f95ce235" in text


def test_release_is_staged_before_unrestrained_pilot():
    text = _read("slurm/run_nylc_c18_step1_continuation.sbatch")
    expected = ["nvt100_nac", "nvt150_nac", "nvt300_nac_release", "npt300_nac_free_pilot"]
    positions = [text.index(stage) for stage in expected]
    assert positions == sorted(positions)
    assert "unrestrained=yes" in text


def test_warmer_stage_contracts_reduce_restraints():
    nvt100 = _read("mdp/nvt100_nac.mdp")
    nvt150 = _read("mdp/nvt150_nac.mdp")
    nvt300 = _read("mdp/nvt300_nac_release.mdp")
    free = _read("mdp/npt300_nac_free_pilot.mdp")
    assert "-DPOSRES_WEAK -DPOSRES_L2_100" in nvt100
    assert "-DPOSRES_WEAK -DPOSRES_L2_10" in nvt150
    assert "-DPOSRES_L2_10" in nvt300 and "-DPOSRES_WEAK" not in nvt300
    assert "-DPOSRES" not in free
    assert "pcoupl                   = Parrinello-Rahman" in free


def test_100ps_free_pilot_is_not_final_scientific_pass():
    text = _read("slurm/run_nylc_c18_step1_continuation.sbatch")
    assert "UNRESTRAINED_100PS_PILOT_COMPLETE" in text
    assert "FINAL_SCIENTIFIC_PASS" not in text
    assert "required_final_unrestrained_ps" in text
    assert "1000" in text
