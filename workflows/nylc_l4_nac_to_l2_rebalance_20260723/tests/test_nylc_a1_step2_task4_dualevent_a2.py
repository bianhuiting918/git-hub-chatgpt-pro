import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "prepare_audit_nylc_a1_step2_task4_dualevent_a2.py"
SBATCH = ROOT / "slurm" / "run_nylc_a1_step2_task4_dualevent_a2.sbatch"

EVENT_HASHES = (
    "cbe9fbddb0ea8c3a6f12c9ae03d19864b6f4d6e2582901fe35ffa828d0f0c1ca",
    "1b0c3f4_placeholder_until_event1_hash_is_bound",
)


def test_task4_dual_event_driver_contract():
    assert DRIVER.is_file()
    text = DRIVER.read_text(encoding="utf-8")
    assert "attempt_62507122_4" in text
    assert "direct_event_{EVENT_INDEX}.rst7" in text
    assert EVENT_HASHES[0] in text
    assert "SLURM_ARRAY_TASK_ID" in text
    assert "SELECTED_WATER_ATOMS1 = (13046, 13047, 13048)" in text
    assert "EXPECTED_DFTB_DOUBLY_OCCUPIED = 194" in text


def test_task4_dual_event_sbatch_contract():
    assert SBATCH.is_file()
    text = SBATCH.read_text(encoding="utf-8")
    assert "#SBATCH --array=0-1" in text
    assert "#SBATCH -n 16" in text
    assert "mpirun --bind-to none -np 8 sander.MPI" in text
    assert "a1_step2_task4_dualevent_a2_preflight" in text
