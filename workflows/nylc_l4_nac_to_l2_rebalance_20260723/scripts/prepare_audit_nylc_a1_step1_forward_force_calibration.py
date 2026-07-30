#!/usr/bin/env python3
"""Run matched NylC A1 reactant-to-product first-window force calibrations."""
from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import re
import shutil
from typing import Any, Mapping

HERE = pathlib.Path(__file__).resolve().parent
CAL_PATH = HERE / "prepare_audit_nylc_a1_step1_reverse_force_calibration.py"
_SPEC = importlib.util.spec_from_file_location("_a1_reverse_force", CAL_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"cannot import {CAL_PATH}")
CAL = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(CAL)
BASE = CAL.BASE
REACTIVE = CAL.REACTIVE
TASK_ROOT = CAL.TASK_ROOT
PRMTOP = CAL.PRMTOP

ARRAY_TASKS = 12
MPI_RANKS = 8
MECHANISMS = ("ADDITION_FIRST_FORWARD", "FULLY_CONCERTED_FORWARD")
FORCE_SCALES = (1.0, 2.0, 4.0)
OUTPUT_ROOT = TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_step1_forward_force_calibration"
SOURCE_ROOT = TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_attack_inherited_chain"
CONTRACT_ROOT = TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_step1_acyl_release_md_continuation"
SOURCES = (
    {
        "seed": "seed26723",
        "attempt": "attempt_62136360_0",
        "contract_attempt": "attempt_62216380_0",
        "restart_sha256": "9d93b60aee9e8d6fa97493396d757f3bf4d0a47595591968e3debaebb2640e86",
    },
    {
        "seed": "seed26737",
        "attempt": "attempt_62156801_1",
        "contract_attempt": "attempt_62216380_1",
        "restart_sha256": "a12c018e195c32b34239004aea36ebc80de8275c9f9332b44084ad8d88d4db0e",
    },
)
FIRST_WINDOW_DELTAS_A = {
    "attack": -0.04,
    "cn": 0.08,
    "carbonyl": 0.03,
    "nalpha_hg1": 0.05,
    "hg1_n3": -0.05,
}
SCIENTIFIC_STATUS = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..11")
    seed_index = index // 6
    local = index % 6
    source = SOURCES[seed_index]
    return {
        "task_index": index,
        "seed_index": seed_index,
        "seed": source["seed"],
        "source_attempt": source["attempt"],
        "contract_attempt": source["contract_attempt"],
        "source_restart_sha256": source["restart_sha256"],
        "mechanism": MECHANISMS[local // 3],
        "force_scale": FORCE_SCALES[local % 3],
        "source_basin": "R",
        "destination": "FIRST_FORWARD_WINDOW",
    }


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "mechanisms": list(MECHANISMS),
        "force_scales": list(FORCE_SCALES),
        "mapping": "2 seeds x 2 mechanisms x 3 force scales",
        "active_groups": {
            "ADDITION_FIRST_FORWARD": ["OG1-C12", "C12-O2"],
            "FULLY_CONCERTED_FORWARD": ["OG1-C12", "C12-N3", "C12-O2", "qPT"],
        },
        "monitored_unrestrained_groups": {
            "ADDITION_FIRST_FORWARD": ["C12-N3", "qPT"],
            "FULLY_CONCERTED_FORWARD": [],
        },
        "first_window_deltas_A": dict(FIRST_WINDOW_DELTAS_A),
        "frozen_step1_contract": BASE.describe()["frozen_step1_contract"],
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }


def _qpt(geometry: Mapping[str, Any]) -> float:
    return float(geometry["nalpha_hg1_A"]) - float(geometry["hg1_n3_A"])


def stage_spec(task: Mapping[str, Any], source: Mapping[str, Any]) -> dict[str, Any]:
    mechanism = str(task["mechanism"])
    if mechanism not in MECHANISMS:
        raise ValueError(mechanism)
    scale = float(task["force_scale"])
    target = {
        "attack_A": float(source["attack_A"]) + FIRST_WINDOW_DELTAS_A["attack"],
        "c12_n3_A": float(source["c12_n3_A"]),
        "c12_o2_A": float(source["c12_o2_A"]) + FIRST_WINDOW_DELTAS_A["carbonyl"],
        "nalpha_hg1_A": float(source["nalpha_hg1_A"]),
        "hg1_n3_A": float(source["hg1_n3_A"]),
    }
    active = ["attack", "carbonyl"]
    if mechanism == "FULLY_CONCERTED_FORWARD":
        target["c12_n3_A"] += FIRST_WINDOW_DELTAS_A["cn"]
        target["nalpha_hg1_A"] += FIRST_WINDOW_DELTAS_A["nalpha_hg1"]
        target["hg1_n3_A"] += FIRST_WINDOW_DELTAS_A["hg1_n3"]
        active = ["attack", "cn", "carbonyl", "pt"]
    return {
        "window_index": 0,
        "mechanism": mechanism,
        "force_scale": scale,
        "targets": target,
        "active_coordinates": active,
        "force_bond": CAL.BASE_FORCES["force_bond"] * scale,
        "force_carbonyl": CAL.BASE_FORCES["force_carbonyl"] * scale,
        "force_pt": CAL.BASE_FORCES["force_pt"] * scale,
        "maxcyc": 2200,
        "ncyc": 550,
    }


def _source_path(task: Mapping[str, Any]) -> pathlib.Path:
    return SOURCE_ROOT / str(task["source_attempt"]) / "terminal_endpoint.rst7"


def initialize(task_index: int, root: pathlib.Path, commit: str) -> dict[str, Any]:
    if root.exists():
        raise FileExistsError(root)
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("runtime is not bound to a full Git commit")
    task = task_spec(task_index)
    restart = _source_path(task)
    chain = BASE.read_json(SOURCE_ROOT / str(task["source_attempt"]) / "CHAIN_MANIFEST.json")
    contract_manifest = BASE.read_json(
        CONTRACT_ROOT / str(task["contract_attempt"]) / "ENDPOINT_MANIFEST.json"
    )
    if not all((
        chain.get("status") == "PASS_TECHNICAL_A1_ATTACK_INHERITED_CHAIN",
        chain.get("inheritance_sha256_verified") is True,
        restart.is_file(),
        CAL.sha256(restart) == task["source_restart_sha256"],
        chain["persistent_restarts"]["terminal_endpoint"]["sha256"]
        == task["source_restart_sha256"],
    )):
        raise ValueError("fixed reactant-side source authority failed")
    contract = CAL._validate_contract(contract_manifest)
    if CAL.sha256(PRMTOP) != BASE.AUTH.BASE.EXPECTED_PRMTOP_SHA256:
        raise ValueError("frozen prmtop SHA changed")
    geometry = BASE.AUTH._geometry(restart, {"qm_contract": contract})
    root.mkdir(parents=True, exist_ok=False)
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_FORWARD_FORCE_CALIBRATION",
        "github_commit": commit,
        "scientific_status": SCIENTIFIC_STATUS,
        "task": task,
        "source": {
            "restart": str(restart),
            "restart_sha256": task["source_restart_sha256"],
            "chain_manifest": str(
                SOURCE_ROOT / str(task["source_attempt"]) / "CHAIN_MANIFEST.json"
            ),
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": CAL.sha256(PRMTOP),
        "qm_contract": contract,
        "source_geometry": geometry,
        "first_window_only": True,
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }
    BASE.write_json(root / "MANIFEST.json", manifest)
    return manifest


def restraints(stage: Mapping[str, Any]) -> str:
    target = stage["targets"]
    rows = [
        CAL.BRIDGE._tight_distance_restraint(
            REACTIVE["og1"], REACTIVE["c12"], target["attack_A"], stage["force_bond"]
        ),
        CAL.BRIDGE._tight_distance_restraint(
            REACTIVE["c12"], REACTIVE["o2"], target["c12_o2_A"],
            stage["force_carbonyl"],
        ),
    ]
    if stage["mechanism"] == "FULLY_CONCERTED_FORWARD":
        rows.extend((
            CAL.BRIDGE._tight_distance_restraint(
                REACTIVE["c12"], REACTIVE["n3"], target["c12_n3_A"],
                stage["force_bond"],
            ),
            CAL.BRIDGE._tight_distance_restraint(
                REACTIVE["nalpha"], REACTIVE["hg1"], target["nalpha_hg1_A"],
                stage["force_pt"],
            ),
            CAL.BRIDGE._tight_distance_restraint(
                REACTIVE["hg1"], REACTIVE["n3"], target["hg1_n3_A"],
                stage["force_pt"],
            ),
        ))
    return "".join(rows)


def prepare(root: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest = BASE.read_json(root / "MANIFEST.json")
    if scratch.exists():
        raise FileExistsError(scratch)
    scratch.mkdir(parents=True, exist_ok=False)
    stage = stage_spec(manifest["task"], manifest["source_geometry"])
    (scratch / "stage.in").write_text(
        CAL.BRIDGE.tetra_minimization_input(
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


def _write_hashes(root: pathlib.Path) -> None:
    rows = []
    for path in sorted(root.iterdir()):
        if path.is_file() and path.name != "SHA256.tsv":
            rows.append(f"{CAL.sha256(path)}  {path.name}")
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
            "attack": CAL.BRIDGE._response(
                previous["attack_A"], geometry["attack_A"], target["attack_A"]
            ),
            "cn": CAL.BRIDGE._response(
                previous["c12_n3_A"], geometry["c12_n3_A"], target["c12_n3_A"]
            ),
            "carbonyl": CAL.BRIDGE._response(
                previous["c12_o2_A"], geometry["c12_o2_A"], target["c12_o2_A"]
            ),
            "pt": CAL.BRIDGE._response(
                _qpt(previous), _qpt(geometry),
                target["nalpha_hg1_A"] - target["hg1_n3_A"],
            ),
            "active_coordinates": list(stage["active_coordinates"]),
        }
        response["all"] = all(response[name]["pass"] for name in response["active_coordinates"])
        guard = CAL.chemical_guard(geometry)
        shutil.copy2(scratch / "stage.rst7", root / "stage.rst7")
    for name in ("stage.in", "restraints.RST", "stage.mdinfo"):
        source = scratch / name
        if source.is_file():
            shutil.copy2(source, root / name)
    stage_out = scratch / "stage.out"
    if stage_out.is_file():
        tail = stage_out.read_text(encoding="utf-8", errors="replace").splitlines()[-200:]
        (root / "ENGINE_TAIL.txt").write_text("\n".join(tail) + "\n", encoding="utf-8")
    qualified = bool(technical and response["all"] and guard["all"])
    gate = (
        "NOT_EVALUATED_TECHNICAL_A1_FORWARD_FORCE_CALIBRATION"
        if not technical else
        "PASS_FIRST_FORWARD_WINDOW_FORCE_RESPONSE"
        if qualified else
        "FAIL_FIRST_FORWARD_WINDOW_ACTUAL_RESPONSE_OR_GUARD"
    )
    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_A1_FORWARD_FORCE_CALIBRATION"
            if technical else "NOT_EVALUATED_TECHNICAL_A1_FORWARD_FORCE_CALIBRATION"
        ),
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
        "eligible_for_forward_continuation": qualified,
        "first_window_only": True,
        "automatic_downstream_action": "NONE",
        "diagnostics": diagnostics,
        "restrained_structures_are_not_ts": True,
    }
    BASE.write_json(root / "RESULT.json", result)
    BASE.write_json(root / ("PASS.json" if technical else "NOT_EVALUATED.json"), result)
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
                and item.get("eligible_for_forward_continuation") is True
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
    payload = {
        "schema_version": 1,
        "status": (
            "NOT_EVALUATED_TECHNICAL_A1_FORWARD_FORCE_CALIBRATION"
            if technical != ARRAY_TASKS else
            "PASS_FORWARD_FORCE_CALIBRATION_MATRIX"
            if all(item["selected_force_scale"] is not None for item in selected) else
            "PARTIAL_FORWARD_FORCE_CALIBRATION_MATRIX"
        ),
        "denominator_tasks": ARRAY_TASKS,
        "technical_pass_tasks": technical,
        "qualified_first_window_tasks": sum(
            item.get("eligible_for_forward_continuation") is True for item in results
        ),
        "selected_minimum_force_scales": selected,
        "per_task": sorted(results, key=lambda item: item["task"]["task_index"]),
        "scientific_status": SCIENTIFIC_STATUS,
        "automatic_downstream_action": "NONE",
    }
    audit = output_root / "audit" / f"nylc_a1_forward_force_calibration_{array_job}.json"
    audit.parent.mkdir(parents=True, exist_ok=True)
    if audit.exists() and BASE.read_json(audit) != payload:
        raise FileExistsError(audit)
    if not audit.exists():
        BASE.write_json(audit, payload)
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=(
        "describe", "initialize", "prepare", "audit", "merge-if-ready",
    ))
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
