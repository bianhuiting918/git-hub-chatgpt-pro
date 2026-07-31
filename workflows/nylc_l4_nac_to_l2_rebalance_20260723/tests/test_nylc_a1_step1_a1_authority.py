#!/usr/bin/env python3
"""Regression tests for A1 covalent-integrity and no-bias authority."""
from __future__ import annotations

import importlib.util
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
DRIVER = ROOT / "scripts" / "prepare_audit_nylc_a1_step1_a1_authority.py"


def load_driver():
    spec = importlib.util.spec_from_file_location("_a1_authority", DRIVER)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_covalent_integrity_accepts_intact_a1_and_rejects_observed_breaks():
    module = load_driver()
    intact = {
        "nalpha_ca_A": 1.47,
        "ca_cb_A": 1.54,
        "cb_og1_A": 1.43,
        "nalpha_h1_A": 1.02,
        "nalpha_h2_A": 1.02,
        "nalpha_hg1_A": 1.03,
    }
    expected = {"H1": "Nalpha", "H2": "Nalpha", "HG1": "Nalpha"}
    good = module.chemical_integrity_from_metrics(intact, expected)
    assert good["pass"] is True
    assert good["criterion"] == "A1_THR267_COVALENT_INTEGRITY_V1"

    broken = dict(intact)
    broken["nalpha_ca_A"] = 3.70
    broken["ca_cb_A"] = 2.80
    broken["cb_og1_A"] = 2.05
    bad = module.chemical_integrity_from_metrics(broken, expected)
    assert bad["pass"] is False
    assert bad["checks"]["nalpha_ca_covalent"] is False
    assert bad["checks"]["ca_cb_covalent"] is False
    assert bad["checks"]["cb_og1_covalent"] is False


def test_wrong_proton_acceptor_is_technical_integrity_failure():
    module = load_driver()
    intact = {
        "nalpha_ca_A": 1.47,
        "ca_cb_A": 1.54,
        "cb_og1_A": 1.43,
        "nalpha_h1_A": 1.02,
        "nalpha_h2_A": 1.02,
        "nalpha_hg1_A": 1.03,
    }
    nearest = {"H1": "CA", "H2": "Nalpha", "HG1": "Nalpha"}
    result = module.chemical_integrity_from_metrics(intact, nearest)
    assert result["pass"] is False
    assert result["checks"]["all_nalpha_hydrogens_attached"] is False


def test_authority_contract_is_two_seed_and_reaction_coordinate_free():
    module = load_driver()
    description = module.describe()
    assert description["array_task_count"] == 2
    assert description["mpi_ranks_per_task"] == 8
    assert description["reaction_coordinate_restraints"] == 0
    assert description["source_restart_shas"] == {
        "26723": "2440de548c385f092c37f683de7b379ff5b18b6dc16593f7dbc80a9a8a167e14",
        "26737": "48c3944d3295158b06e96e32e4e07d9bcae9ceba2731f95aae9c1335ff972abe",
    }


def test_reclassification_never_mutates_legacy_rows():
    module = load_driver()
    legacy = {
        "status": "PASS_TECHNICAL_A1_RAW_NAC_PREORGANIZED_CALIBRATION",
        "technical_complete": True,
        "eligible_for_inherited_chain": True,
    }
    integrity = {"pass": False, "criterion": "A1_THR267_COVALENT_INTEGRITY_V1"}
    revised = module.reclassify_legacy_result(legacy, integrity)
    assert legacy["technical_complete"] is True
    assert legacy["eligible_for_inherited_chain"] is True
    assert revised["status"] == "NOT_EVALUATED_TECHNICAL_CHEMICAL_INTEGRITY"
    assert revised["technical_complete"] is False
    assert revised["eligible_for_inherited_chain"] is False
    assert revised["legacy_result_preserved"] is True
