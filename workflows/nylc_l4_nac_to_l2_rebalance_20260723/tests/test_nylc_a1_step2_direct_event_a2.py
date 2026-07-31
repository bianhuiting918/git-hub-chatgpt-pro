import importlib.util
import os
from pathlib import Path


WORKFLOW_ROOT = Path(os.environ["A1_WORKFLOW_ROOT"])
ADAPTER = (
    WORKFLOW_ROOT
    / "scripts"
    / "prepare_audit_nylc_a1_step2_direct_event_a2.py"
)


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
