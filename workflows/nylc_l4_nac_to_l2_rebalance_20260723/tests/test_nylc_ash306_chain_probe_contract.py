from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROBE = ROOT / "scripts" / "prepare_nylc_ash306_chain_probe.py"
SBATCH = ROOT / "slurm" / "run_nylc_ash306_chain_probe.sbatch"


def test_probe_uses_immutable_434ps_source_and_exact_active_chain():
    text = PROBE.read_text()
    assert "182f68a41c72c490e0046323c9b8a0f8514ccf7af7dc4d578dcb87342c1c8e2d" in text
    assert "CHAIN_H_FIRST = 8949" in text
    assert "CHAIN_H_LAST = 10272" in text
    assert 'resid == 306 and resname == "ASP"' in text
    assert 'resname = "ASH"' in text
    assert 'chain_id = "H"' in text


def test_probe_audits_one_added_hd2_and_nterminal_thr267():
    text = PROBE.read_text()
    assert "EXPECTED_CHAIN_ATOMS = 1325" in text
    assert '"ASH", "HD2"' in text
    assert '"THR", {"H1", "H2", "H3", "OG1", "HG1"}' in text
    assert "PROBE_PASS_ASH306_CHAIN" in text
    assert "not yet a full-system topology" in text


def test_wrapper_is_cpu_only_pdb2gmx_probe_with_history():
    text = SBATCH.read_text()
    assert "#SBATCH -p xahcnormal" in text
    assert "pdb2gmx" in text
    assert "-ff amber99sb-ildn" in text
    assert "-water tip3p" in text
    assert "-ignh" in text
    assert "mdrun" not in text
    assert "ash306_chain_probe_job_" in text
    assert "run_history.tsv" in text and "run_history.jsonl" in text
    assert "trap on_error ERR" in text


def test_probe_exports_methionine_sulfur_as_pdb_element_s():
    namespace = {}
    exec(PROBE.read_text(), namespace)
    assert namespace["element"]("SD") == "S"


def test_wrapper_accepts_single_chain_topology_written_directly_to_top():
    text = SBATCH.read_text()
    assert 'itp="$OUT/chainH.top"' in text
