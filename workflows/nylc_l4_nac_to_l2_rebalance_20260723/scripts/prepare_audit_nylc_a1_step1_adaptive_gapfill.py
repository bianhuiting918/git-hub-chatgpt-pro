#!/usr/bin/env python3
"""Adaptive inherited continuation of the eight incomplete NylC A1 Step1 routes.

This reuses the frozen four-anchor Hamiltonian, targets, chemistry guards, and
engine audit.  Restrained routes are not TS, committor, PMF, barrier, or
mechanism evidence.
"""
from __future__ import annotations

import argparse
import fcntl
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
BASE_PATH = HERE / "prepare_audit_nylc_a1_step1_four_anchor_bidirectional_gapfill.py"
_SPEC = importlib.util.spec_from_file_location("_a1_four_anchor_authority", BASE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"cannot import {BASE_PATH}")
BASE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(BASE)

TASK_ROOT = BASE.TASK_ROOT
PRMTOP = BASE.PRMTOP
SOURCE_JOB = "62267272"
SOURCE_TASKS = (2, 3, 4, 5, 10, 11, 12, 13)
SOURCE_ROOT = (
    TASK_ROOT / "a1_activated_nac_20260726/qmmm/"
    "a1_step1_four_anchor_bidirectional_gapfill"
)
ARRAY_TASKS = 8
MPI_RANKS = 8
MAXIMUM_NEW_WINDOWS = 8
MINIMUM_BRACKET_FRACTION = 1.0 / 16.0
TECHNICAL_PASS = "PASS_TECHNICAL_ADAPTIVE_GAPFILL"
TECHNICAL_FAIL = "NOT_EVALUATED_TECHNICAL_ADAPTIVE_GAPFILL"
SCIENTIFIC_STATUS = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: pathlib.Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def midpoint(
    low: Mapping[str, float], high: Mapping[str, float]
) -> dict[str, float]:
    if set(low) != set(high):
        raise ValueError("midpoint mappings must have identical keys")
    return {key: (float(low[key]) + float(high[key])) / 2.0 for key in low}


def next_bracket(low: float, high: float, passed: bool) -> tuple[float, float]:
    middle = (float(low) + float(high)) / 2.0
    return (middle, float(high)) if passed else (float(low), middle)


def interpolate(
    low: Mapping[str, float], high: Mapping[str, float], fraction: float
) -> dict[str, float]:
    if set(low) != set(high):
        raise ValueError("target mappings must have identical keys")
    f = float(fraction)
    if not 0.0 <= f <= 1.0:
        raise ValueError("fraction must be between 0 and 1")
    return {
        key: float(low[key]) + f * (float(high[key]) - float(low[key]))
        for key in low
    }


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..7")
    source_task_index = SOURCE_TASKS[index]
    source_task = BASE.task_spec(source_task_index)
    return {
        "task_index": index,
        "source_task_index": source_task_index,
        "seed_index": source_task["seed_index"],
        "seed": source_task["seed"],
        "anchor": source_task["anchor"],
        "direction": source_task["direction"],
    }


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "source_job": SOURCE_JOB,
        "source_tasks": list(SOURCE_TASKS),
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "maximum_new_windows": MAXIMUM_NEW_WINDOWS,
        "minimum_bracket_fraction": MINIMUM_BRACKET_FRACTION,
        "strict_serial_restart_inheritance": True,
        "failed_restart_inheritance": False,
        "frozen_step1_contract": {
            "qm_atom_count": 146,
            "qmcharge": 0,
            "electron_count": 510,
            "link_atom_count": 6,
            "qm_water_count": 0,
        },
        "restrained_structures_are_not_ts": True,
        "automatic_downstream_action": "NONE",
    }


def _load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def _source_authority(task: Mapping[str, Any]) -> dict[str, Any]:
    root = SOURCE_ROOT / f"attempt_{SOURCE_JOB}_{task['source_task_index']}"
    manifest_path = root / "MANIFEST.json"
    result_path = root / "RESULT.json"
    if not manifest_path.is_file() or not result_path.is_file():
        raise ValueError("source attempt lacks MANIFEST.json or RESULT.json")
    source_manifest = _load_json(manifest_path)
    source_result = _load_json(result_path)
    if (
        source_result.get("technical_complete") is not True
        or source_result.get("status") != BASE.TECHNICAL_PASS
        or source_result.get("stop_reason") != "SCIENTIFIC_GUARD_STOP"
        or len(source_manifest.get("stages", [])) != 4
    ):
        raise ValueError("source attempt is not the frozen four-stage guard-stop authority")
    source_task = source_result.get("task", {})
    for key in ("seed", "anchor", "direction"):
        if source_task.get(key) != task[key]:
            raise ValueError(f"source task mismatch for {key}")
    safe_indices = [
        index for index, stage in enumerate(source_manifest["stages"])
        if stage.get("technical_pass") is True and stage.get("guard_stop") is False
    ]
    if not safe_indices:
        raise ValueError("source route has no accepted guard-safe stage")
    safe_index = max(safe_indices)
    rejected_index = safe_index + 1
    if rejected_index >= len(source_manifest["stages"]):
        raise ValueError("source route does not contain the first rejected stage")
    rejected = source_manifest["stages"][rejected_index]
    if rejected.get("technical_pass") is not True or rejected.get("guard_stop") is not True:
        raise ValueError("source next stage is not the frozen guard rejection")
    safe = source_manifest["stages"][safe_index]
    safe_restart = root / f"stage_{safe_index:02d}" / "stage.rst7"
    if not safe_restart.is_file() or sha256(safe_restart) != safe["output_restart_sha256"]:
        raise ValueError("persisted last safe restart SHA authority failed")
    if sha256(PRMTOP) != BASE.BASE.EXPECTED_PRMTOP_SHA256:
        raise ValueError("frozen prmtop SHA changed")
    return {
        "root": root,
        "manifest": source_manifest,
        "result": source_result,
        "safe_index": safe_index,
        "rejected_index": rejected_index,
        "safe_stage": safe,
        "rejected_stage": rejected,
        "safe_restart": safe_restart,
        "safe_restart_sha256": safe["output_restart_sha256"],
    }


def _terminal_state(manifest: Mapping[str, Any]) -> tuple[bool, str | None]:
    state = manifest["adaptive_state"]
    windows = manifest["windows"]
    if manifest.get("technical_failure"):
        return True, "NOT_EVALUATED_TECHNICAL"
    if state.get("destination_reached"):
        return True, "DESTINATION_REACHED"
    if len(windows) >= MAXIMUM_NEW_WINDOWS:
        return True, "MAXIMUM_WINDOWS"
    width = float(state["high_fraction"]) - float(state["low_fraction"])
    if width <= MINIMUM_BRACKET_FRACTION and windows and not windows[-1].get("accepted_for_inheritance"):
        return True, "MINIMUM_BRACKET_REFINED"
    return False, None


def initialize(task_index: int, root: pathlib.Path, commit: str) -> dict[str, Any]:
    if root.exists():
        raise FileExistsError(root)
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("runtime is not bound to a full Git commit")
    task = task_spec(task_index)
    authority = _source_authority(task)
    source_manifest = authority["manifest"]
    contract = source_manifest["qm_contract"]
    observed_contract = {
        "qm_atom_count": contract.get("qm_atom_count"),
        "qmcharge": contract.get("qmcharge"),
        "electron_count": contract.get("electron_count_including_link_h"),
        "link_atom_count": contract.get("link_atom_count"),
        "qm_water_count": contract.get("step1_qm_water_count"),
    }
    if observed_contract != describe()["frozen_step1_contract"]:
        raise ValueError(f"frozen Step1 QM contract mismatch: {observed_contract}")
    original_specs = source_manifest["stage_specs"]
    next_original_index = authority["rejected_index"]
    low_targets = dict(authority["safe_stage"]["stage_spec"]["targets"])
    high_targets = dict(original_specs[next_original_index]["targets"])
    root.mkdir(parents=True, exist_ok=False)
    accepted = root / "accepted.rst7"
    shutil.copy2(authority["safe_restart"], accepted)
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_ADAPTIVE_GAPFILL",
        "github_commit": commit,
        "scientific_status": SCIENTIFIC_STATUS,
        "task": task,
        "source": {
            "job": SOURCE_JOB,
            "attempt": str(authority["root"]),
            "safe_stage_index": authority["safe_index"],
            "rejected_stage_index": authority["rejected_index"],
            "safe_restart": str(authority["safe_restart"]),
            "safe_restart_sha256": authority["safe_restart_sha256"],
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": sha256(PRMTOP),
        "qm_contract": contract,
        "original_stage_specs": original_specs,
        "accepted_restart": str(accepted),
        "accepted_restart_sha256": sha256(accepted),
        "adaptive_state": {
            "original_stage_index": next_original_index,
            "low_targets": low_targets,
            "high_targets": high_targets,
            "low_fraction": 0.0,
            "high_fraction": 1.0,
            "destination_reached": False,
        },
        "windows": [],
        "technical_failure": False,
        "restrained_structures_are_not_ts": True,
        "automatic_downstream_action": "NONE",
    }
    write_json(root / "MANIFEST.json", manifest)
    return manifest


def _candidate_fraction(state: Mapping[str, Any]) -> float:
    low = float(state["low_fraction"])
    high = float(state["high_fraction"])
    width = high - low
    if width <= MINIMUM_BRACKET_FRACTION:
        return high
    return (low + high) / 2.0


def prepare_window(
    root: pathlib.Path,
    output: pathlib.Path,
    scratch: pathlib.Path,
    input_rst7: pathlib.Path,
    window_index: int,
) -> dict[str, Any]:
    manifest = _load_json(root / "MANIFEST.json")
    terminal, reason = _terminal_state(manifest)
    if terminal:
        raise ValueError(f"route is terminal: {reason}")
    if int(window_index) != len(manifest["windows"]):
        raise ValueError("window order is not strict")
    if sha256(input_rst7) != manifest["accepted_restart_sha256"]:
        raise ValueError("window input breaks accepted restart SHA inheritance")
    state = manifest["adaptive_state"]
    fraction = _candidate_fraction(state)
    target = interpolate(state["low_targets"], state["high_targets"], fraction)
    original_index = int(state["original_stage_index"])
    stage = dict(manifest["original_stage_specs"][original_index])
    stage["targets"] = target
    stage["adaptive_fraction"] = fraction
    stage["adaptive_window_index"] = int(window_index)
    stage["stage_index"] = int(window_index)
    stage["target_anchor"] = (
        manifest["original_stage_specs"][original_index].get("target_anchor")
        if math.isclose(fraction, 1.0, abs_tol=1.0e-12) else None
    )
    output.mkdir(parents=True, exist_ok=False)
    scratch.mkdir(parents=True, exist_ok=False)
    (scratch / "stage.in").write_text(
        BASE.minimization_input(manifest["task"], stage, manifest["qm_contract"]["qmmask"]),
        encoding="utf-8",
    )
    (scratch / "restraints.RST").write_text(BASE.restraints(stage), encoding="utf-8")
    prepared = {
        "schema_version": 1,
        "window_index": int(window_index),
        "original_stage_index": original_index,
        "candidate_fraction": fraction,
        "stage_spec": stage,
        "input_restart": str(input_rst7),
        "input_restart_sha256": manifest["accepted_restart_sha256"],
        "restrained": True,
        "shooting_forbidden": True,
    }
    write_json(output / "WINDOW_MANIFEST.json", prepared)
    return prepared


def audit_window(
    root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path
) -> dict[str, Any]:
    manifest_path = root / "MANIFEST.json"
    manifest = _load_json(manifest_path)
    prepared = _load_json(output / "WINDOW_MANIFEST.json")
    restart = scratch / "stage.rst7"
    technical, geometry, diagnostics = BASE.TETRA._technical(
        scratch / "stage.out", restart, manifest
    )
    stage = prepared["stage_spec"]
    target_anchor = stage.get("target_anchor")
    chemical = BASE.REVERSE.chemical_guard(geometry) if technical else {"all": False}
    checkpoint = (
        BASE.anchor_gate(target_anchor, geometry)
        if technical and target_anchor else None
    )
    guard_pass = bool(
        technical
        and chemical.get("all") is True
        and (checkpoint is None or checkpoint.get("all") is True)
    )
    persistent_restart = output / "stage.rst7"
    if technical:
        shutil.copy2(restart, persistent_restart)
    accepted = bool(technical and guard_pass)
    if accepted:
        shutil.copy2(persistent_restart, root / "accepted.rst7")
        manifest["accepted_restart_sha256"] = sha256(root / "accepted.rst7")
    state = manifest["adaptive_state"]
    fraction = float(prepared["candidate_fraction"])
    if technical:
        if accepted:
            state["low_fraction"] = fraction
            if math.isclose(fraction, 1.0, abs_tol=1.0e-12):
                original_index = int(state["original_stage_index"])
                if original_index == len(manifest["original_stage_specs"]) - 1:
                    state["destination_reached"] = True
                else:
                    state["original_stage_index"] = original_index + 1
                    state["low_targets"] = dict(state["high_targets"])
                    state["high_targets"] = dict(
                        manifest["original_stage_specs"][original_index + 1]["targets"]
                    )
                    state["low_fraction"] = 0.0
                    state["high_fraction"] = 1.0
        else:
            state["high_fraction"] = fraction
    else:
        manifest["technical_failure"] = True
    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_A1_ADAPTIVE_GAPFILL_WINDOW"
            if technical else "NOT_EVALUATED_TECHNICAL_A1_ADAPTIVE_GAPFILL_WINDOW"
        ),
        "technical_pass": technical,
        "window_index": prepared["window_index"],
        "original_stage_index": prepared["original_stage_index"],
        "fraction": fraction,
        "target": stage["targets"],
        "input_restart_sha256": prepared["input_restart_sha256"],
        "output_restart": str(persistent_restart) if technical else None,
        "output_restart_sha256": sha256(persistent_restart) if persistent_restart.is_file() else None,
        "geometry": geometry,
        "chemical_guard": chemical,
        "anchor_gate": checkpoint,
        "guard_pass": guard_pass,
        "accepted_for_inheritance": accepted,
        "restrained": True,
        "shooting_forbidden": True,
        "diagnostics": diagnostics,
    }
    manifest["windows"].append(result)
    terminal, reason = _terminal_state(manifest)
    manifest["terminal"] = terminal
    manifest["terminal_reason"] = reason
    write_json(output / "RESULT.json", result)
    write_json(output / ("PASS.json" if technical else "NOT_EVALUATED.json"), result)
    write_json(manifest_path, manifest)
    return result


def finalize(root: pathlib.Path) -> dict[str, Any]:
    manifest = _load_json(root / "MANIFEST.json")
    terminal, reason = _terminal_state(manifest)
    if not terminal:
        raise ValueError("cannot finalize a non-terminal adaptive route")
    windows = manifest["windows"]
    technical = bool(
        windows
        and not manifest["technical_failure"]
        and all(window.get("technical_pass") is True for window in windows)
    )
    inheritance = True
    previous = manifest["source"]["safe_restart_sha256"]
    for window in windows:
        inheritance = inheritance and window.get("input_restart_sha256") == previous
        if window.get("accepted_for_inheritance"):
            previous = window.get("output_restart_sha256")
    destination = "R" if manifest["task"]["direction"] == "REACTANT" else "P"
    final_geometry = next(
        (window["geometry"] for window in reversed(windows) if window.get("accepted_for_inheritance")),
        {},
    )
    destination_gate = BASE.anchor_gate(destination, final_geometry) if final_geometry else {"all": False}
    reached = bool(
        technical
        and manifest["adaptive_state"].get("destination_reached")
        and destination_gate.get("all") is True
    )
    scientific_gate = (
        "NOT_EVALUATED_TECHNICAL"
        if not technical
        else "PASS_GUIDED_ROUTE_REACHED_DESTINATION"
        if reached
        else "PARTIAL_BRACKET_REFINED"
    )
    result = {
        "schema_version": 1,
        "status": TECHNICAL_PASS if technical else TECHNICAL_FAIL,
        "technical_complete": technical,
        "scientific_status": SCIENTIFIC_STATUS,
        "scientific_gate": scientific_gate,
        "task": manifest["task"],
        "source_task_index": manifest["task"]["source_task_index"],
        "windows_attempted": len(windows),
        "accepted_windows": sum(window.get("accepted_for_inheritance") is True for window in windows),
        "terminal_reason": reason,
        "inheritance_sha256_verified": inheritance,
        "destination_anchor": destination,
        "destination_gate": destination_gate,
        "remaining_bracket_fraction": (
            manifest["adaptive_state"]["high_fraction"]
            - manifest["adaptive_state"]["low_fraction"]
        ),
        "restrained_structures_are_not_ts": True,
        "automatic_downstream_action": "NONE",
    }
    write_json(root / "RESULT.json", result)
    write_json(root / ("PASS.json" if technical else "NOT_EVALUATED.json"), result)
    rows = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.name != "SHA256.tsv" and path.suffix in (".json", ".rst7"):
            rows.append(f"{sha256(path)}  {path.relative_to(root)}")
    (root / "SHA256.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    return result


def merge_if_ready(output_root: pathlib.Path, array_job: str) -> bool:
    audit = output_root / "audit" / f"nylc_a1_adaptive_gapfill_{array_job}.json"
    audit.parent.mkdir(parents=True, exist_ok=True)
    with audit.with_suffix(".lock").open("a+", encoding="utf-8") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        paths = [
            output_root / f"attempt_{array_job}_{index}" / "RESULT.json"
            for index in range(ARRAY_TASKS)
        ]
        if not all(path.is_file() for path in paths):
            return False
        results = [_load_json(path) for path in paths]
        technical = all(result.get("technical_complete") is True for result in results)
        reached = sum(
            result.get("scientific_gate") == "PASS_GUIDED_ROUTE_REACHED_DESTINATION"
            for result in results
        )
        payload = {
            "schema_version": 1,
            "status": (
                "NOT_EVALUATED_TECHNICAL_A1_ADAPTIVE_GAPFILL"
                if not technical
                else "PASS_ALL_EIGHT_GUIDED_ROUTES_REACHED_DESTINATION"
                if reached == ARRAY_TASKS
                else "PARTIAL_ADAPTIVE_GAPFILL"
            ),
            "denominator_tasks": ARRAY_TASKS,
            "technical_pass_tasks": sum(
                result.get("technical_complete") is True for result in results
            ),
            "guided_destination_pass_tasks": reached,
            "per_task": sorted(results, key=lambda item: item["task"]["task_index"]),
            "scientific_status": SCIENTIFIC_STATUS,
            "restrained_structures_are_not_ts": True,
            "automatic_downstream_action": "NONE",
        }
        if audit.exists():
            if _load_json(audit) != payload:
                raise FileExistsError(audit)
        else:
            write_json(audit, payload)
        return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=(
        "describe", "initialize", "prepare-window", "audit-window",
        "finalize", "merge-if-ready",
    ))
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--window-index", type=int)
    parser.add_argument("--root", type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--output-root", type=pathlib.Path)
    parser.add_argument("--scratch", type=pathlib.Path)
    parser.add_argument("--input-rst7", type=pathlib.Path)
    parser.add_argument("--github-commit", default="unknown")
    parser.add_argument("--array-job")
    args = parser.parse_args()
    if args.mode == "describe":
        print(json.dumps(describe(), sort_keys=True))
    elif args.mode == "initialize":
        initialize(args.task_index, args.root.resolve(), args.github_commit)
    elif args.mode == "prepare-window":
        prepare_window(
            args.root.resolve(), args.output.resolve(), args.scratch.resolve(),
            args.input_rst7.resolve(), args.window_index,
        )
    elif args.mode == "audit-window":
        audit_window(args.root.resolve(), args.output.resolve(), args.scratch.resolve())
    elif args.mode == "finalize":
        finalize(args.root.resolve())
    else:
        merge_if_ready(args.output_root.resolve(), args.array_job)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
