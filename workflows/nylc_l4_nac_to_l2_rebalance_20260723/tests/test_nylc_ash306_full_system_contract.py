import runpy
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BUILDER = ROOT / "scripts" / "build_nylc_ash306_full_system.py"


def load_builder():
    assert BUILDER.exists(), "full-system builder is not implemented"
    return runpy.run_path(str(BUILDER))


def test_builder_locks_immutable_source_probe_and_atom_ranges():
    ns = load_builder()
    assert ns["SOURCE_GRO_SHA256"] == "182f68a41c72c490e0046323c9b8a0f8514ccf7af7dc4d578dcb87342c1c8e2d"
    assert ns["PROBE_GRO_SHA256"] == "a3998b6b5f5df5a40c78861dfd4b65983dff5d96bba1cd7ed238c7d57528ea44"
    assert ns["PROBE_TOP_SHA256"] == "630b45eb20713b1f604ae11a193ded778edcfb52f44b70c9377d6858bd7690f9"
    assert (ns["CHAIN_H_FIRST"], ns["CHAIN_H_LAST"]) == (8949, 10272)
    assert (ns["L2_FIRST"], ns["L2_LAST"]) == (10273, 10351)
    assert ns["SOURCE_ATOMS"] == ns["FINAL_ATOMS"] == 133589


def test_builder_selects_farthest_sodium_under_periodic_boundaries():
    ns = load_builder()
    cell = ((10.0, 0.0, 0.0), (0.0, 10.0, 0.0), (0.0, 0.0, 10.0))
    solute = [(0.0, 0.0, 0.0), (1.0, 0.0, 0.0)]
    sodiums = [
        {"source_index": 100, "xyz_nm": (9.9, 0.0, 0.0)},
        {"source_index": 101, "xyz_nm": (5.0, 5.0, 5.0)},
    ]
    selected = ns["select_farthest_sodium"](sodiums, solute, cell)
    assert selected["source_index"] == 101
    assert selected["minimum_solute_distance_nm"] > 8.0


def test_builder_updates_only_na_molecule_count():
    ns = load_builder()
    top = """[ molecules ]\nProtein_chain_H 1\nPA66_L2 1\nSOL 40990\nNA 144\nCL 124\n"""
    updated = ns["update_na_molecule_count"](top)
    assert "NA 143" in updated
    assert "SOL 40990" in updated and "CL 124" in updated
    assert updated.count("NA 143") == 1


def test_builder_extracts_only_chain_molecule_topology():
    ns = load_builder()
    top = """#include \"forcefield.itp\"\n[ moleculetype ]\nProtein_chain_H 3\n[ atoms ]\n1 N 267 THR N 1 0.1 14\n#ifdef POSRES\n#include \"posre.itp\"\n#endif\n#include \"tip3p.itp\"\n[ system ]\nProtein\n[ molecules ]\nProtein_chain_H 1\n"""
    extracted = ns["extract_chain_itp"](top)
    assert extracted.startswith("[ moleculetype ]")
    assert "[ atoms ]" in extracted
    assert "forcefield.itp" not in extracted
    assert "tip3p.itp" not in extracted
    assert "[ system ]" not in extracted
    assert "POSRES" not in extracted


def test_builder_contract_keeps_gate_and_scientific_boundaries_explicit():
    text = BUILDER.read_text() if BUILDER.exists() else ""
    assert "GATE_RESIDUES = tuple(range(261, 267))" in text
    assert "Thr267 is excluded from the gate group" in text
    assert "not a TS, RC, PMF, or barrier" in text
    assert "minimum_ligand_protein_heavy_distance_nm" in text
    assert "nac_distance_nm" in text and "nac_angle_deg" in text


SBATCH = ROOT / "slurm" / "run_nylc_ash306_full_system_preflight.sbatch"
AUDITOR = ROOT / "scripts" / "audit_nylc_ash306_full_system_preflight.py"


def test_wrapper_is_cpu_only_grompp_without_dynamics_and_records_history():
    assert SBATCH.exists(), "full-system preflight wrapper is not implemented"
    text = SBATCH.read_text()
    assert "#SBATCH -p xahcnormal" in text
    assert "#SBATCH -c 4" in text
    assert "gmx_mpi" in text and "grompp" in text
    assert "mdrun" not in text
    assert "ash306_full_system_preflight_job_" in text
    assert "61717760" in text
    assert "trap on_error ERR" in text
    assert "run_history.tsv" in text and "run_history.jsonl" in text
    assert "PASS_ASH306_FULL_SYSTEM_PREFLIGHT" in text


def test_auditor_requires_neutral_133589_atom_parmed_system():
    assert AUDITOR.exists(), "full-system preflight auditor is not implemented"
    text = AUDITOR.read_text()
    assert "import parmed as pmd" in text
    assert "EXPECTED_ATOMS = 133589" in text
    assert "abs(total_charge)" in text
    assert "PASS_ASH306_FULL_SYSTEM_PREFLIGHT" in text
    assert "PROBE_PASS_ASH306_CHAIN" in text
    assert "not a TS, RC, PMF, or barrier" in text
