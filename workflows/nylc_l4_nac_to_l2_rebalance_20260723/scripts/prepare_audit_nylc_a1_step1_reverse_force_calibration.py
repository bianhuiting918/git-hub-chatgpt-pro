#!/usr/bin/env python3
"""Run the frozen NylC A1 product-to-reactant first-window force calibration."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import os
import pathlib
import re
import shutil
from typing import Any, Mapping

HERE = pathlib.Path(__file__).resolve().parent
BRIDGE_PATH = HERE / "prepare_audit_nylc_a1_step1_tetrahedral_bridge.py"
_SPEC = importlib.util.spec_from_file_location("_a1_tetra_bridge", BRIDGE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"cannot import {BRIDGE_PATH}")
BRIDGE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(BRIDGE)
BASE = BRIDGE.BASE
REACTIVE = BASE.REACTIVE
TASK_ROOT = BASE.TASK_ROOT
PRMTOP = BASE.PRMTOP

ARRAY_TASKS = 12
MPI_RANKS = 8
MECHANISMS = ("STEPWISE_REVERSE", "CONCERTED_REVERSE")
FORCE_SCALES = (1.0, 2.0, 4.0)
OUTPUT_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/a1_step1_reverse_force_calibration"
)
SOURCE_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/a1_step1_acyl_release_md_continuation"
)
SOURCES = (
    {
        "seed": "seed26723",
        "attempt": "attempt_62216380_0",
        "restart_sha256": "5d8f76d2c90e3e8c707b640c55a93938f53e18dc30d6f93d54adda467da25f41",
    },
    {
        "seed": "seed26737",
        "attempt": "attempt_62216380_1",
        "restart_sha256": "4cc60ad4d7be099b3f76040f1ee8c49b91deecebb20ed98570bc192b3d511132",
    },
)
BASE_FORCES = {
    "force_bond": 24.0,
    "force_carbonyl": 30.0,
    "force_pt": 18.0,
}
FIRST_WINDOW_DELTAS_A = {
    "attack": 0.04,
    "cn": -0.08,
    "carbonyl": 0.03,
    "nalpha_hg1": -0.05,
    "hg1_n3": 0.05,
}
TECHNICAL_PASS = "PASS_TECHNICAL_A1_REVERSE_FORCE_CALIBRATION"
TECHNICAL_FAIL = "NOT_EVALUATED_TECHNICAL_A1_REVERSE_FORCE_CALIBRATION"
SCIENTIFIC_STATUS = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..11")
    seed_index = index // 6
    local = index % 6
    mechanism = MECHANISMS[local // 3]
    force_scale = FORCE_SCALES[local % 3]
    source = SOURCES[seed_index]
    return {
        "task_index": index,
        "seed_index": seed_index,
        "seed": source["seed"],
        "source_attempt": source["attempt"],
        "source_restart_sha256": source["restart_sha256"],
        "source_basin": "P",
        "destination": "FIRST_REVERSE_WINDOW",
        "mechanism": mechanism,
        "force_scale": force_scale,
    }


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "mechanisms": list(MECHANISMS),
        "force_scales": list(FORCE_SCALES),
        "mapping": "2 seeds x 2 mechanisms x 3 force scales",
        "first_window_only": True,
        "strict_serial_restart_inheritance": True,
        "failed_restart_inheritance": False,
        "stepwise_active_groups": ["C12-N3", "qPT", "C12-O2"],
        "concerted_active_groups": ["OG1-C12", "C12-N3", "qPT", "C12-O2"],
        "first_window_deltas_A": dict(FIRST_WINDOW_DELTAS_A),
        "frozen_step1_contract": BASE.describe()["frozen_step1_contract"],
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }


def stage_spec(task: Mapping[str, Any], source: Mapping[str, Any]) -> dict[str, Any]:
    mechanism = str(task["mechanism"])
    if mechanism not in MECHANISMS:
        raise ValueError(mechanism)
    scale = float(task["force_scale"])
    target = {
        "attack_A": float(source["attack_A"]),
        "c12_n3_A": float(source["c12_n3_A"]) + FIRST_WINDOW_DELTAS_A["cn"],
        "c12_o2_A": float(source["c12_o2_A"]) + FIRST_WINDOW_DELTAS_A["carbonyl"],
        "nalpha_hg1_A": float(source["nalpha_hg1_A"])
        + FIRST_WINDOW_DELTAS_A["nalpha_hg1"],
        "hg1_n3_A": float(source["hg1_n3_A"]) + FIRST_WINDOW_DELTAS_A["hg1_n3"],
    }
    if mechanism == "CONCERTED_REVERSE":
        target["attack_A"] = (
            float(source["attack_A"]) + FIRST_WINDOW_DELTAS_A["attack"]
        )
    return {
        "window_index": 0,
        "mechanism": mechanism,
        "force_scale": scale,
        "targets": target,
        "attack_coordinate_active": mechanism == "CONCERTED_REVERSE",
        "force_bond": BASE_FORCES["force_bond"] * scale,
        "force_carbonyl": BASE_FORCES["force_carbonyl"] * scale,
        "force_pt": BASE_FORCES["force_pt"] * scale,
        "maxcyc": 2200,
        "ncyc": 550,
    }


def _source_path(task: Mapping[str, Any]) -> pathlib.Path:
    return SOURCE_ROOT / str(task["source_attempt"]) / "release_md_endpoint.rst7"


def _source_manifest(task: Mapping[str, Any]) -> pathlib.Path:
    return SOURCE_ROOT / str(task["source_attempt"]) / "ENDPOINT_MANIFEST.json"


def _validate_contract(payload: Mapping[str, Any]) -> dict[str, Any]:
    contract = dict(payload["qm_contract"])
    checks = (
        int(contract["qm_atom_count"]) == 146,
        int(contract["qmcharge"]) == 0,
        int(contract["electron_count_including_link_h"]) == 510,
        int(contract["link_atom_count"]) == 6,
        int(contract["step1_qm_water_count"]) == 0,
        bool(contract.get("qmmask")),
    )
    if not all(checks):
        raise ValueError("frozen Step1 QM contract changed")
    return contract


def initialize(task_index: int, root: pathlib.Path, commit: str) -> dict[str, Any]:
    if root.exists():
        raise FileExistsError(root)
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("runtime is not bound to a full Git commit")
    task = task_spec(task_index)
    restart = _source_path(task)
    source_manifest = _source_manifest(task)
    if not restart.is_file() or sha256(restart) != task["source_restart_sha256"]:
        raise ValueError("fixed product-side restart SHA failed")
    if not source_manifest.is_file():
        raise FileNotFoundError(source_manifest)
    endpoint = BASE.read_json(source_manifest)
    contract = _validate_contract(endpoint)
    if sha256(PRMTOP) != BASE.AUTH.BASE.EXPECTED_PRMTOP_SHA256:
        raise ValueError("frozen prmtop SHA changed")
    geometry = BASE.AUTH._geometry(restart, {"qm_contract": contract})
    root.mkdir(parents=True, exist_ok=False)
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_REVERSE_FORCE_CALIBRATION",
        "github_commit": commit,
        "scientific_status": SCIENTIFIC_STATUS,
        "task": task,
        "source": {
            "restart": str(restart),
            "restart_sha256": task["source_restart_sha256"],
            "authority": "PROVISIONAL_ACYL_CONNECTIVITY_SENSITIVITY_V1",
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": sha256(PRMTOP),
        "qm_contract": contract,
        "source_geometry": geometry,
        "first_window_only": True,
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }
    BASE.write_json(root / "MANIFEST.json", manifest)
    return manifest


def _qpt(geometry: Mapping[str, Any]) -> float:
    return float(geometry["nalpha_hg1_A"]) - float(geometry["hg1_n3_A"])


def restraints(stage: Mapping[str, Any]) -> str:
    target = stage["targets"]
    rows = []
    if stage["attack_coordinate_active"]:
        rows.append(
            BRIDGE._tight_distance_restraint(
                REACTIVE["og1"], REACTIVE["c12"], target["attack_A"],
                stage["force_bond"],
            )
        )
    rows.extend(
        (
            BRIDGE._tight_distance_restraint(
                REACTIVE["c12"], REACTIVE["n3"], target["c12_n3_A"],
                stage["force_bond"],
            ),
            BRIDGE._tight_distance_restraint(
                REACTIVE["c12"], REACTIVE["o2"], target["c12_o2_A"],
                stage["force_carbonyl"],
            ),
            BRIDGE._tight_distance_restraint(
                REACTIVE["nalpha"], REACTIVE["hg1"], target["nalpha_hg1_A"],
                stage["force_pt"],
            ),
            BRIDGE._tight_distance_restraint(
                REACTIVE["hg1"], REACTIVE["n3"], target["hg1_n3_A"],
                stage["force_pt"],
            ),
        )
    )
    return "".join(rows)


def prepare(root: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest = BASE.read_json(root / "MANIFEST.json")
    if scratch.exists():
        raise FileExistsError(scratch)
    scratch.mkdir(parents=True, exist_ok=False)
    stage = stage_spec(manifest["task"], manifest["source_geometry"])
    (scratch / "stage.in").write_text(
        BRIDGE.tetra_minimization_input(
            manifest["task"], stage, manifest["qm_contract"]["qmmask"]
        ),
        encoding="utf-8",
    )
    (scratch / "restraints.RST").write_text(restraints(stage), encoding="utf-8")
    prepared = {
        "schema_version": 1,
        "stage": stage,
        "input_restart": manifest["source"]["restart"],
        "input_restart_sha256": manifest["source"]["restart_sha256"],
        "first_window_only": True,
        "restrained": True,
        "shooting_forbidden": True,
    }
    BASE.write_json(root / "WINDOW_MANIFEST.json", prepared)
    return prepared


def chemical_guard(geometry: Mapping[str, Any]) -> dict[str, Any]:
    values = (
        geometry["attack_A"], geometry["c12_n3_A"], geometry["c12_o2_A"],
        geometry["nalpha_hg1_A"], geometry["hg1_n3_A"],
    )
    checks = {
        "finite_geometry": all(math.isfinite(float(value)) for value in values),
        "attack_product_bond_not_overcompressed_or_lost": (
            1.25 <= float(geometry["attack_A"]) <= 1.85
        ),
        "cn_not_overcompressed": float(geometry["c12_n3_A"]) >= 1.25,
        "carbonyl_safe": 1.15 <= float(geometry["c12_o2_A"]) <= 1.55,
        "proton_not_detached_or_misrouted": (
            geometry.get("hg1_nearest_qm_heavy_atom")
            in (REACTIVE["nalpha"], REACTIVE["n3"])
        ),
    }
    checks["all"] = all(checks.values())
    return checks


def _write_hashes(root: pathlib.Path) -> None:
    rows = []
    for path in sorted(root.iterdir()):
        if path.is_file() and path.name != "SHA256.tsv":
            rows.append(f"{sha256(path)}  {path.name}")
    (root / "SHA256.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")


def audit(root: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest = BASE.read_json(root / "MANIFEST.json")
    prepared = BASE.read_json(root / "WINDOW_MANIFEST.json")
    stage = prepared["stage"]
    technical, geometry, diagnostics = BASE.AUTH.TETRA._technical(
        scratch / "stage.out", scratch / "stage.rst7", manifest
    )
    previous = manifest["source_geometry"]
    target = stage["targets"]
    response = {"all": False}
    guard = {"all": False}
    if technical:
        response = {
            "attack": BRIDGE._response(
                previous["attack_A"], geometry["attack_A"], target["attack_A"]
            ),
            "cn": BRIDGE._response(
                previous["c12_n3_A"], geometry["c12_n3_A"], target["c12_n3_A"]
            ),
            "carbonyl": BRIDGE._response(
                previous["c12_o2_A"], geometry["c12_o2_A"], target["c12_o2_A"]
            ),
            "pt": BRIDGE._response(
                _qpt(previous), _qpt(geometry),
                target["nalpha_hg1_A"] - target["hg1_n3_A"],
            ),
        }
        active = ["cn", "carbonyl", "pt"]
        if stage["attack_coordinate_active"]:
            active.append("attack")
        response["active_coordinates"] = active
        response["all"] = all(response[name]["pass"] for name in active)
        guard = chemical_guard(geometry)
        shutil.copy2(scratch / "stage.rst7", root / "stage.rst7")
    for name in ("stage.in", "restraints.RST"):
        shutil.copy2(scratch / name, root / name)
    stage_out = scratch / "stage.out"
    if stage_out.is_file():
        tail = stage_out.read_text(encoding="utf-8", errors="replace").splitlines()[-200:]
        (root / "ENGINE_TAIL.txt").write_text("\n".join(tail) + "\n", encoding="utf-8")
    qualified = bool(technical and response["all"] and guard["all"])
    gate = (
        "NOT_EVALUATED_TECHNICAL_A1_REVERSE_FORCE_CALIBRATION"
        if not technical
        else "PASS_FIRST_WINDOW_FORCE_RESPONSE"
        if qualified
        else "FAIL_FIRST_WINDOW_ACTUAL_RESPONSE_OR_GUARD"
    )
    result = {
        "schema_version": 1,
        "status": TECHNICAL_PASS if technical else TECHNICAL_FAIL,
        "technical_complete": technical,
        "scientific_status": SCIENTIFIC_STATUS,
        "scientific_gate": gate,
        "task": manifest["task"],
        "source_geometry": previous,
        "target_geometry": target,
        "final_geometry": geometry if technical else None,
        "source_qPT_A": _qpt(previous),
        "target_qPT_A": target["nalpha_hg1_A"] - target["hg1_n3_A"],
        "final_qPT_A": _qpt(geometry) if technical else None,
        "response_gate": response,
        "chemical_guard": guard,
        "eligible_for_inherited_chain": qualified,
        "first_window_only": True,
        "automatic_downstream_action": "NONE",
        "diagnostics": diagnostics,
        "restrained_structures_are_not_ts": True,
    }
    BASE.write_json(root / "RESULT.json", result)
    BASE.write_json(
        root / ("PASS.json" if technical else "NOT_EVALUATED.json"), result
    )
    _write_hashes(root)
    return result


def merge_if_ready(output_root: pathlib.Path, array_job: str) -> bool:
    paths = [
        output_root / f"attempt_{array_job}_{index}" / "RESULT.json"
        for index in range(ARRAY_TASKS)
    ]
    if not all(path.is_file() for path in paths):
        return False
    results = [BASE.read_json(path) for path in paths]
    selected = []
    for seed in ("seed26723", "seed26737"):
        for mechanism in MECHANISMS:
            candidates = [
                item for item in results
                if item["task"]["seed"] == seed
                and item["task"]["mechanism"] == mechanism
                and item.get("eligible_for_inherited_chain") is True
            ]
            selected.append({
                "seed": seed,
                "mechanism": mechanism,
                "selected_force_scale": (
                    min(item["task"]["force_scale"] for item in candidates)
                    if candidates else None
                ),
            })
    technical = sum(item.get("technical_complete") is True for item in results)
    qualified = sum(
        item.get("eligible_for_inherited_chain") is True for item in results
    )
    complete_selection = all(
        item["selected_force_scale"] is not None for item in selected
    )
    payload = {
        "schema_version": 1,
        "status": (
            "NOT_EVALUATED_TECHNICAL_A1_REVERSE_FORCE_CALIBRATION"
            if technical != ARRAY_TASKS
            else "PASS_REVERSE_FORCE_CALIBRATION_MATRIX"
            if complete_selection
            else "PARTIAL_REVERSE_FORCE_CALIBRATION_MATRIX"
        ),
        "denominator_tasks": ARRAY_TASKS,
        "technical_pass_tasks": technical,
        "qualified_first_window_tasks": qualified,
        "selected_minimum_force_scales": selected,
        "per_task": sorted(results, key=lambda item: item["task"]["task_index"]),
        "scientific_status": SCIENTIFIC_STATUS,
        "automatic_downstream_action": "NONE",
    }
    audit = output_root / "audit" / f"nylc_a1_reverse_force_calibration_{array_job}.json"
    audit.parent.mkdir(parents=True, exist_ok=True)
    if audit.exists() and BASE.read_json(audit) != payload:
        raise FileExistsError(audit)
    if not audit.exists():
        BASE.write_json(audit, payload)
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode", required=True,
        choices=("describe", "initialize", "prepare", "audit", "merge-if-ready"),
    )
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--root", type=pathlib.Path)
    parser.add_argument("--scratch", type=pathlib.Path)
    parser.add_argument("--github-commit", default="unknown")
    parser.add_argument("--output-root", type=pathlib.Path)
    parser.add_argument("--array-job")
    args = parser.parse_args()
    if args.mode == "describe":
        print(json.dumps(describe(), sort_keys=True))
    elif args.mode == "initialize":
        initialize(args.task_index, args.root.resolve(), args.github_commit)
    elif args.mode == "prepare":
        prepare(args.root.resolve(), args.scratch.resolve())
    elif args.mode == "audit":
        audit(args.root.resolve(), args.scratch.resolve())
    else:
        merge_if_ready(args.output_root.resolve(), args.array_job)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
