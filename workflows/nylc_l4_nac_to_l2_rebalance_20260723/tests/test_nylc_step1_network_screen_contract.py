from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ANALYZER = ROOT / "scripts" / "analyze_nylc_step1_network_nac.py"
SBATCH = ROOT / "slurm" / "run_nylc_step1_network_screen.sbatch"


def test_exact_c18_reactive_and_active_copy_atoms_are_fixed():
    text = ANALYZER.read_text()
    expected = {
        "TYR146_N": 7156,
        "TYR146_H": 7157,
        "ASN219_ND2": 8240,
        "ASN219_HD21": 8241,
        "ASN219_HD22": 8242,
        "THR267_N": 8949,
        "THR267_OG1": 8961,
        "ASP306_OD1": 9572,
        "ASP306_OD2": 9573,
        "ASP308_OD1": 9591,
        "ASP308_OD2": 9592,
        "L2_REACTIVE_C": 10287,
        "L2_REACTIVE_O": 10288,
        "L2_REACTIVE_N": 10289,
    }
    for name, atom in expected.items():
        assert f"{name} = {atom}" in text


def test_analyzer_requires_strict_nac_and_both_oxyanion_donors():
    text = ANALYZER.read_text()
    assert "distance_nm <= 0.35" in text
    assert "95.0 <= attack_angle_deg <= 115.0" in text
    assert "donor_heavy_max_nm = 0.35" in text
    assert "donor_h_max_nm = 0.25" in text
    assert "donor_angle_min_deg = 135.0" in text
    assert "tyr_oxyanion_ready and asn_oxyanion_ready" in text
    assert "Asp306 protonation microstate" in text


def test_analyzer_never_promotes_energy_only_nac():
    text = ANALYZER.read_text()
    assert "NOT_EVALUATED_NO_FULL_CATALYTIC_PREORGANIZATION" in text
    assert "lowest_potential_full_preorganization_frame" in text
    assert "source.tmp.gro" in text
    assert "selected_step1_network_gs.gro" in text


def test_wrapper_is_cpu_only_analysis_and_never_runs_md():
    text = SBATCH.read_text()
    assert "#SBATCH -p xahcnormal" in text
    assert "--gres=dcu" not in text
    assert "mdrun" not in text
    assert "npt300free/run.xtc" in text
    assert "npt300free/run.tpr" in text
    assert "\\nFREE=" not in text
    assert "\nFREE=" in text
    assert "postprocess_job_61710861" in text
    assert "analyze_nylc_step1_network_nac.py" in text


def test_wrapper_records_history_and_uses_unique_outputs():
    text = SBATCH.read_text()
    assert "run_history.tsv" in text
    assert "run_history.jsonl" in text
    assert "step1_network_job_" in text
    assert "source.tmp.gro" in text
    assert "GMX_MAXBACKUP=-1" in text
