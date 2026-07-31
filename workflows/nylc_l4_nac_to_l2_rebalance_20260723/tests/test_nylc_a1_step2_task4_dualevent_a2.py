import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "prepare_audit_nylc_a1_step2_best_direct_pair_a2.py"
SBATCH = ROOT / "slurm" / "run_nylc_a1_step2_best_direct_pair_a2.sbatch"

SOURCE_HASHES = (
    "cbe9fbddb0ea8c3a6f12c9ae03d19864b6f4d6e2582901fe35ffa828d0f0c1ca",
    "5c4227dc8752f1a181e5b17d2890e6bce5a8a790f3ed5f91a181445e72006efe",
)


def test_best_direct_pair_driver_contract():
    assert DRIVER.is_file()
    text = DRIVER.read_text(encoding="utf-8")
    assert "SOURCE_TASKS = (4, 6)" in text
    assert "attempt_62507122_{SOURCE_TASK}" in text
    for digest in SOURCE_HASHES:
        assert digest in text
    assert "SELECTED_WATER_ATOMS1 = (13046, 13047, 13048)" in text
    assert "EXPECTED_DFTB_DOUBLY_OCCUPIED = 194" in text


def test_best_direct_pair_sbatch_contract():
    assert SBATCH.is_file()
    text = SBATCH.read_text(encoding="utf-8")
    assert "#SBATCH --array=0-1" in text
    assert "#SBATCH -n 16" in text
    assert "mpirun --bind-to none -np 8 sander.MPI" in text
    assert "a1_step2_best_direct_pair_a2_preflight" in text
