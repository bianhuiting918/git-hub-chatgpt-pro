from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
MDP = ROOT / "mdp" / "nvt100_step1_oxyanion_recapture_response.mdp"
AUDIT = ROOT / "scripts" / "audit_nylc_step1_oxyanion_recapture.py"
SBATCH = ROOT / "slurm" / "run_nylc_step1_oxyanion_recapture_response.sbatch"
MDP2 = ROOT / "mdp" / "nvt100_step1_oxyanion_recapture_response2.mdp"
SBATCH2 = ROOT / "slurm" / "run_nylc_step1_oxyanion_recapture_response2.sbatch"


def test_response_mdp_is_bounded_and_not_scientific_sampling():
    text = MDP.read_text()
    assert "ref-t                    = 100" in text
    assert "nsteps                   = 50000" in text
    assert "dt                       = 0.001" in text
    assert "pcoupl                   = no" in text
    assert "pull-ncoords             = 6" in text
    assert "pull-coord1-init         = 0.330000" in text
    assert "pull-coord2-init         = 105.000000" in text
    assert text.count("pull-coord3-rate         = -0.001000") == 1
    assert text.count("pull-coord5-rate         = -0.001000") == 1
    assert "pull-coord4-groups       = 5 4 5 3" in text
    assert "pull-coord6-groups       = 7 6 7 3" in text
    assert "Gate" not in text


def test_audit_uses_exact_c18_l2_network_and_cannot_scientifically_pass():
    text = AUDIT.read_text()
    expected = {
        "TYR146_N": 7156,
        "TYR146_H": 7157,
        "ASN219_ND2": 8240,
        "ASN219_HD21": 8241,
        "THR267_OG1": 8961,
        "L2_C12": 10287,
        "L2_O2": 10288,
    }
    for name, atom in expected.items():
        assert f"{name} = {atom}" in text
    assert "NOT_EVALUATED_RESTRAINED_OXYANION_RECAPTURE_PILOT" in text
    assert "PASS_TECHNICAL_RESPONSE" in text
    assert "minimum_heavy_ligand_protein_distance_nm" in text


def test_wrapper_extracts_immutable_210ps_source_and_records_terminal_state():
    text = SBATCH.read_text()
    assert "#SBATCH -p xahdnormal" in text
    assert "--gres=dcu:1" in text
    assert "npt300free/run.xtc" in text
    assert "-dump 210" in text
    assert "source.tmp.gro" in text
    assert "source.gro" in text
    assert "grompp" in text and "-maxwarn 0" in text
    assert "run_history.tsv" in text and "run_history.jsonl" in text
    assert "trap on_error ERR" in text
    assert '"$GMX" make_ndx -f "$OUT/source.gro" -o "$OUT/recapture.ndx"' in text
    assert 'cp "$BUILD/source_cycle.ndx"' not in text
    assert r"<<'NDX'\n[" not in text
    assert "<<'NDX'\n[ Thr267_OG1 ]" in text
    assert "[ Tyr146_N ]" in text
    assert "[ Asn219_HD21 ]" in text
    assert "[ Gate ]" not in text


def test_response2_increases_only_bounded_donor_force_and_restarts_same_source():
    mdp = MDP2.read_text()
    assert mdp.count("pull-coord3-k            = 2000") == 1
    assert mdp.count("pull-coord5-k            = 2000") == 1
    assert "pull-coord3-rate         = -0.001000" in mdp
    assert "pull-coord5-rate         = -0.001000" in mdp
    assert "pull-coord1-init         = 0.330000" in mdp
    wrapper = SBATCH2.read_text()
    assert "-dump 210" in wrapper
    assert "step1_oxyanion_recapture_response2_job_" in wrapper
    assert "nvt100_step1_oxyanion_recapture_response2.mdp" in wrapper
    assert "run_history.tsv" in wrapper and "trap on_error ERR" in wrapper
