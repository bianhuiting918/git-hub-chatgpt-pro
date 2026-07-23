from pathlib import Path

FLOW = Path(__file__).resolve().parents[1]


def test_all_wrong_nucleophile_routes_are_hard_blocked():
    paths = [
        "slurm/run_nylc_step1_pilot_array.sbatch",
        "slurm/run_nylc_c18_step1_continuation.sbatch",
        "slurm/run_nylc_c18_step1_free_1ns.sbatch",
        "scripts/audit_nylc_c18_free_1ns.sh",
    ]
    for relative in paths:
        text = (FLOW / relative).read_text()
        assert "SUPERSEDED_WRONG_NUCLEOPHILE_IDENTITY" in text
        assert "exit 42" in text
        assert text.index("exit 42") < text.index("TASK_ROOT=")
