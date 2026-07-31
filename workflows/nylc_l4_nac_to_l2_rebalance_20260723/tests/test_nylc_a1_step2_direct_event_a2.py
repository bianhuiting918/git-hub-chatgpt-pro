import importlib.util
import os
from pathlib import Path


WORKFLOW_ROOT = Path(os.environ["A1_WORKFLOW_ROOT"])
ADAPTER = (
    WORKFLOW_ROOT
    / "scripts"
    / "prepare_audit_nylc_a1_step2_direct_event_a2.py"
)
SBATCH = WORKFLOW_ROOT / "slurm" / "run_nylc_a1_step2_direct_event_a2.sbatch"


def load_adapter():
    assert ADAPTER.is_file(), "missing task6 direct-event A2 adapter"
    spec = importlib.util.spec_from_file_location("direct_event_a2", ADAPTER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_task6_source_and_water_are_immutable():
    module = load_adapter()
    assert module.SOURCE_RESTART_SHA256 == (
        "5c4227dc8752f1a181e5b17d2890e6bce5a8a790f3ed5f91a181445e72006efe"
    )
    assert module.SELECTED_WATER_ATOMS1 == (13046, 13047, 13048)
    assert module.SELECTED_DONOR_H1 == 13048
    assert module.ALLOWED_REACTED_TOPOLOGY_BONDS1 == {
        (8949, 8961),
        (10287, 10289),
    }


def test_adapter_exposes_single_seed_base_driver_contract():
    module = load_adapter()
    source = module.source_for_index(0)
    assert source["seed"] == "seed26737"
    assert source["restart_sha256"] == module.SOURCE_RESTART_SHA256
    assert source["velocity_seeds"] == (26737621, 26737622)
    assert callable(module.fixed_select_water)
    assert callable(module.validate_task6_authority)


def test_adapter_normalizes_repeated_atom_mask_prefixes_for_amber18():
    module = load_adapter()
    assert module.normalize_amber_atom_mask("@1,@2,@3") == "@1,2,3"
    assert module.normalize_amber_atom_mask("@1,2,3") == "@1,2,3"


def test_sbatch_uses_cluster_working_openmpi_launcher_for_both_legs():
    text = SBATCH.read_text(encoding="utf-8")
    assert "srun --exclusive" not in text
    assert text.count("mpirun --bind-to none -np 8 sander.MPI") == 1


def test_direct_engine_banner_uses_qmmm_options_and_dftb_valence_authority():
    module = load_adapter()
    text = """
QMMM options:
             ifqnt = True       nquant =      149
              qmgb =        0  qmcharge =        0   adjust_q =        2
|QMMM: Running QMMM calculation in parallel mode on    8 threads.
QMMM:  nlink =     6                   Link Coords
QMMM: RHF CALCULATION, NO. OF DOUBLY OCCUPIED LEVELS =194
 NSTEP =        1
"""
    observed = module.parse_direct_engine_banner(text)
    assert observed["qm_atom_count"] == [149]
    assert observed["qmcharge"] == [0]
    assert observed["link_atom_count"] == [6]
    assert observed["doubly_occupied_levels"] == [194]
    manifest = {
        "qm_contract": {
            "expected": dict(module.BASE.EXPECTED_CONTRACT),
            "qmmask": "@1-149",
        }
    }
    prepared = {"expected_contract": dict(module.BASE.EXPECTED_CONTRACT)}
    effective, source = module.resolve_direct_leg_banner_contract(
        manifest, prepared, observed
    )
    assert effective == {
        "qm_atom_count": [149],
        "qmcharge": [0],
        "link_atom_count": [6],
        "electron_count": [518],
    }
    assert source == "LEG_ENGINE_DFTB_VALENCE_PLUS_FIXED_COMPOSITION"


def test_sbatch_persists_reauditable_engine_and_geometry_artifacts():
    text = SBATCH.read_text(encoding="utf-8")
    for name in (
        "stage.out",
        "stage.mdinfo",
        "stage.rst7",
        "a2.mdcrd",
        "engine.rc",
    ):
        assert name in text
    assert "A2_LEG_${LEG}.stage.out" in text
    assert "A2_LEG_${LEG}.a2.mdcrd" in text


def test_adapter_adds_thr267_heavy_skeleton_integrity_audit():
    module = load_adapter()
    assert callable(module.audit_thr267_heavy_skeleton)
    assert callable(module.audit_leg_with_persisted_thr_integrity)
