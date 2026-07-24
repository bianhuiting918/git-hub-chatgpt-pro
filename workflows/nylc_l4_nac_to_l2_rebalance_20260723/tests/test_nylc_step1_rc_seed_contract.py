from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "analyze_nylc_step1_rc_seed.py"
SBATCH = ROOT / "slurm" / "run_nylc_step1_rc_seed_audit.sbatch"


def test_rc_seed_locks_both_microstates_and_reactive_atoms():
    assert SCRIPT.exists(), "Step1 RC seed analyzer is not implemented"
    text = SCRIPT.read_text()
    assert "step1_dftb3_preflight_post61710861_job61712692" in text
    assert "ash306_full_system_preflight_job_61718715" in text
    for token in (
        "THR267_N = 8949", "THR267_OG1 = 8961", "THR267_HG1 = 8962",
        "L2_C12_DEPROT = 10287", "L2_O2_DEPROT = 10288", "L2_N3_DEPROT = 10289",
        "ASP306_OD1 = 9572", "ASP306_OD2 = 9573",
        "ASP308_OD1_DEPROT = 9591", "ASP308_OD2_DEPROT = 9592",
    ):
        assert token in text


def test_rc_seed_keeps_two_asp306_acceptor_hypotheses_and_no_ts_claim():
    text = SCRIPT.read_text() if SCRIPT.exists() else ""
    assert "pt_og1_hg1_to_asp306_od1" in text
    assert "pt_og1_hg1_to_asp306_od2" in text
    assert "attack_og1_to_l2_c12" in text
    assert "amide_c12_n3" in text and "carbonyl_c12_o2" in text
    assert "candidate basin thresholds" in text
    assert "not a validated RC, TS, committor, PMF, or barrier" in text


def test_rc_seed_wrapper_is_light_cpu_only_and_audited():
    assert SBATCH.exists(), "Step1 RC seed wrapper is not implemented"
    text = SBATCH.read_text()
    assert "#SBATCH -p xahcnormal" in text and "#SBATCH -c 1" in text
    assert "analyze_nylc_step1_rc_seed.py" in text
    assert "mdrun" not in text and "sander" not in text
    assert "run_history.tsv" in text and "run_history.jsonl" in text
    assert "PASS_STEP1_RC_SEED_AUDIT" in text


def test_rc_seed_wrapper_respects_partition_memory_per_cpu_limit():
    text = SBATCH.read_text()
    assert "#SBATCH --mem-per-cpu=2500M" in text


def test_rc_seed_distance_uses_full_xyz_vectors_not_scalar_x():
    text = SCRIPT.read_text()
    assert "structure.coordinates[left - 1]" in text
    assert "structure.coordinates[right - 1]" in text
    assert ".xx" not in text
