from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ANALYZER = ROOT / "scripts" / "analyze_nylc_step1_water_network.py"
SBATCH = ROOT / "slurm" / "run_nylc_step1_water_network.sbatch"


def test_analyzer_uses_authoritative_434ps_atoms_and_source_sha():
    text = ANALYZER.read_text()
    expected = {
        "THR267_OG1": 8961,
        "LYS189_NZ": 7768,
        "ASN219_OD1": 8239,
        "ASN219_ND2": 8240,
        "ASP306_OD1": 9572,
        "ASP306_OD2": 9573,
        "ASP308_OD1": 9591,
        "ASP308_OD2": 9592,
        "L2_C12": 10287,
        "L2_O2": 10288,
    }
    for name, atom in expected.items():
        assert f"{name} = {atom}" in text
    assert "182f68a41c72c490e0046323c9b8a0f8514ccf7af7dc4d578dcb87342c1c8e2d" in text


def test_analyzer_uses_water_network_not_direct_substrate_oxyanion_gate():
    text = ANALYZER.read_text()
    assert "WATER_EDGE_MAX_NM = 0.40" in text
    assert "thr267_asn219_water_bridges" in text
    assert "asp306_asp308_water_bridges" in text
    assert "lys189_asn219_direct_min_nm" in text
    assert "PASS_STRUCTURAL_GS_HYPOTHESIS" in text
    assert "not proof of a proton-transfer path" in text


def test_wrapper_is_cpu_only_and_records_unique_audit():
    text = SBATCH.read_text()
    assert "#SBATCH -p xahcnormal" in text
    assert "--gres=dcu" not in text
    assert "mdrun" not in text
    assert "step1_dftb3_preflight_post61710861_job61712692/selected_free_l2_nac.gro" in text
    assert "step1_water_network_job_" in text
    assert "run_history.tsv" in text and "run_history.jsonl" in text
    assert "trap on_error ERR" in text
