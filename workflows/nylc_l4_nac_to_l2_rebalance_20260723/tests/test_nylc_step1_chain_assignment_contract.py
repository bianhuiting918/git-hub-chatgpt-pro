from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
ANALYZER = ROOT / "scripts" / "analyze_nylc_step1_chain_assignment.py"
SBATCH = ROOT / "slurm" / "run_nylc_step1_chain_assignment.sbatch"


def test_first_copy_oxyanion_atoms_are_explicit():
    text = ANALYZER.read_text()
    expected = {
        "COPY1_TYR146_N": 2020,
        "COPY1_TYR146_H": 2021,
        "COPY1_ASN219_ND2": 3104,
        "COPY1_ASN219_HD21": 3105,
        "COPY1_ASN219_HD22": 3106,
        "L2_REACTIVE_O": 10288,
    }
    for name, atom in expected.items():
        assert f"{name} = {atom}" in text


def test_all_four_subunit_donor_combinations_are_audited():
    text = ANALYZER.read_text()
    for key in (
        "copy1_tyr_copy1_asn",
        "copy1_tyr_copy2_asn",
        "copy2_tyr_copy1_asn",
        "copy2_tyr_copy2_asn",
    ):
        assert key in text
    assert "NOT_EVALUATED_NO_CROSS_SUBUNIT_OXYANION_PREORGANIZATION" in text
    assert "lowest_potential_preorganized_frame" in text


def test_wrapper_is_cpu_only_reproducible_and_terminally_audited():
    text = SBATCH.read_text()
    assert "#SBATCH -p xahcnormal" in text
    assert "mdrun" not in text
    assert "step1_network_job_61713551/result/step1_network_frames.tsv" in text
    assert "run_history.tsv" in text
    assert "run_history.jsonl" in text
    assert "trap on_error ERR" in text
    angle_command = text.split('| "$GMX" angle', 1)[1].split('> "$OUT/angle.stdout"', 1)[0]
    assert '-s "$FREE/run.tpr"' not in angle_command
