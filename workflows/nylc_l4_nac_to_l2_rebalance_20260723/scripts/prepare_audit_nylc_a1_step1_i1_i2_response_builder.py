#!/usr/bin/env python3
"""Build true restrained I1/I2 anchors from observed geometry response.

The prior equal-target interpolation did not move the actual reactive geometry.
This continuation therefore advances only from the last accepted restart,
limits each target increment, requires measured coordinate response, and never
inherits a failed restart. Restrained anchors are not TS or barrier evidence.
"""
from __future__ import annotations

import argparse
import fcntl
import hashlib
import importlib.util
import json
import os
import pathlib
import re
import shutil
from typing import Any, Mapping

HERE = pathlib.Path(__file__).resolve().parent
BASE_PATH = HERE / "prepare_audit_nylc_a1_step1_four_anchor_bidirectional_gapfill.py"
_SPEC = importlib.util.spec_from_file_location("_i1_i2_base", BASE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"cannot import {BASE_PATH}")
BASE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(BASE)

TASK_ROOT = BASE.TASK_ROOT
PRMTOP = BASE.PRMTOP
SOURCE_ROOT = TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_step1_adaptive_gapfill"
OUTPUT_NAME = "a1_step1_i1_i2_response_builder"
SOURCE_SPECS = (
    {"seed_index": 0, "seed": "seed26723", "anchor": "I1", "job": "62289293", "task": 1},
    {"seed_index": 0, "seed": "seed26723", "anchor": "I2", "job": "62289293", "task": 2},
    {"seed_index": 1, "seed": "seed26737", "anchor": "I1", "job": "62296486", "task": 5},
    {"seed_index": 1, "seed": "seed26737", "anchor": "I2", "job": "62289293", "task": 6},
)
ARRAY_TASKS = 4
MPI_RANKS = 8
MAXIMUM_WINDOWS = 12
MINIMUM_STEP_SCALE = 0.25
PRIMARY_STEP_A = 0.10
CARBONYL_STEP_A = 0.02
TECHNICAL_PASS = "PASS_TECHNICAL_A1_I1_I2_RESPONSE_BUILDER"
TECHNICAL_FAIL = "NOT_EVALUATED_TECHNICAL_A1_I1_I2_RESPONSE_BUILDER"
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


def load_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def step_toward(current: float, target: float, maximum_step: float) -> float:
    current = float(current)
    target = float(target)
    step = abs(float(maximum_step))
    if abs(target - current) <= step:
        return target
    return current + step if target > current else current - step


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..3")
    return {"task_index": index, **SOURCE_SPECS[index]}


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "maximum_windows": MAXIMUM_WINDOWS,
        "minimum_step_scale": MINIMUM_STEP_SCALE,
        "primary_step_A": PRIMARY_STEP_A,
        "carbonyl_step_A": CARBONYL_STEP_A,
        "actual_target_residual_is_required": True,
        "failed_restart_inheritance": False,
        "strict_serial_restart_inheritance": True,
        "source_specs": list(SOURCE_SPECS),
        "frozen_step1_contract": {
            "qm_atom_count": 146,
            "qmcharge": 0,
            "electron_count": 510,
            "link_atom_count": 6,
            "qm_water_count": 0,
        },
        "automatic_downstream_action": "NONE",
    }


def _source_authority(task: Mapping[str, Any]) -> dict[str, Any]:
    root = SOURCE_ROOT / f"attempt_{task['job']}_{task['task']}"
    manifest_path = root / "MANIFEST.json"
    result_path = root / "RESULT.json"
    restart = root / "accepted.rst7"
    if not manifest_path.is_file() or not result_path.is_file() or not restart.is_file():
        raise ValueError("source adaptive attempt is not terminal")
    manifest = load_json(manifest_path)
    result = load_json(result_path)
    if (
        result.get("technical_complete") is not True
        or result.get("status") != "PASS_TECHNICAL_ADAPTIVE_GAPFILL"
        or result.get("inheritance_sha256_verified") is not True
    ):
        raise ValueError("source adaptive attempt is not technical/SHA authority")
    if manifest.get("accepted_restart_sha256") != sha256(restart):
        raise ValueError("source accepted restart SHA mismatch")
    accepted = [window for window in manifest.get("windows", []) if window.get("accepted_for_inheritance") is True]
    if not accepted:
        raise ValueError("source has no accepted window")
    source_task = manifest.get("task", {})
    if source_task.get("seed") != task["seed"] or source_task.get("anchor") != task["anchor"]:
        raise ValueError("source seed/anchor mismatch")
    contract = manifest.get("qm_contract", {})
    observed = {
        "qm_atom_count": contract.get("qm_atom_count"),
        "qmcharge": contract.get("qmcharge"),
        "electron_count": contract.get("electron_count_including_link_h"),
        "link_atom_count": contract.get("link_atom_count"),
        "qm_water_count": contract.get("step1_qm_water_count"),
    }
    if observed != describe()["frozen_step1_contract"]:
        raise ValueError(f"frozen Step1 contract mismatch: {observed}")
    if sha256(PRMTOP) != BASE.BASE.EXPECTED_PRMTOP_SHA256:
        raise ValueError("frozen prmtop SHA changed")
    return {
        "root": root,
        "restart": restart,
        "restart_sha256": sha256(restart),
        "geometry": dict(accepted[-1]["geometry"]),
        "qm_contract": contract,
    }


def _target_for(anchor: str, geometry: Mapping[str, Any], scale: float) -> dict[str, float]:
    desired = BASE.ANCHOR_TARGETS[anchor]
    return {
        "attack_A": step_toward(geometry["attack_A"], desired["attack_A"], PRIMARY_STEP_A * scale),
        "cn_A": step_toward(geometry["c12_n3_A"], desired["cn_A"], PRIMARY_STEP_A * scale),
        "c12_o2_A": step_toward(geometry["c12_o2_A"], desired["c12_o2_A"], CARBONYL_STEP_A * scale),
        "nalpha_hg1_A": float(geometry["nalpha_hg1_A"]),
        "hg1_n3_A": float(geometry["hg1_n3_A"]),
    }


def _terminal(manifest: Mapping[str, Any]) -> tuple[bool, str | None]:
    state = manifest["adaptive_state"]
    if manifest.get("technical_failure"):
        return True, "NOT_EVALUATED_TECHNICAL"
    if state.get("anchor_reached"):
        return True, "ACTUAL_ANCHOR_REACHED"
    if state.get("stalled_at_target"):
        return True, "TARGET_DISTANCES_WITHOUT_TETRAHEDRAL_RESPONSE"
    if float(state["step_scale"]) < MINIMUM_STEP_SCALE:
        return True, "MINIMUM_ACTUAL_RESPONSE_STEP_FAILED"
    if len(manifest.get("windows", [])) >= MAXIMUM_WINDOWS:
        return True, "MAXIMUM_WINDOWS"
    return False, None


def initialize(task_index: int, root: pathlib.Path, commit: str) -> dict[str, Any]:
    if root.exists():
        raise FileExistsError(root)
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("runtime is not bound to a full Git commit")
    task = task_spec(task_index)
    source = _source_authority(task)
    root.mkdir(parents=True, exist_ok=False)
    accepted = root / "accepted.rst7"
    shutil.copy2(source["restart"], accepted)
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_I1_I2_ACTUAL_RESPONSE_BUILDER",
        "github_commit": commit,
        "scientific_status": SCIENTIFIC_STATUS,
        "task": task,
        "source": {
            "attempt": str(source["root"]),
            "accepted_restart": str(source["restart"]),
            "accepted_restart_sha256": source["restart_sha256"],
            "observed_geometry": source["geometry"],
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": sha256(PRMTOP),
        "qm_contract": source["qm_contract"],
        "accepted_restart": str(accepted),
        "accepted_restart_sha256": sha256(accepted),
        "accepted_geometry": source["geometry"],
        "adaptive_state": {
            "step_scale": 1.0,
            "anchor_reached": False,
            "stalled_at_target": False,
        },
        "windows": [],
        "technical_failure": False,
        "restrained_structures_are_not_ts": True,
        "automatic_downstream_action": "NONE",
    }
    write_json(root / "MANIFEST.json", manifest)
    return manifest


def prepare_window(
    root: pathlib.Path,
    output: pathlib.Path,
    scratch: pathlib.Path,
    input_rst7: pathlib.Path,
    window_index: int,
) -> dict[str, Any]:
    manifest = load_json(root / "MANIFEST.json")
    terminal, reason = _terminal(manifest)
    if terminal:
        raise ValueError(f"route is terminal: {reason}")
    if int(window_index) != len(manifest["windows"]):
        raise ValueError("window order is not strict")
    if sha256(input_rst7) != manifest["accepted_restart_sha256"]:
        raise ValueError("input breaks accepted restart SHA inheritance")
    geometry = manifest["accepted_geometry"]
    scale = float(manifest["adaptive_state"]["step_scale"])
    targets = _target_for(manifest["task"]["anchor"], geometry, scale)
    stage = {
        "stage_index": int(window_index),
        "targets": targets,
        "maxcyc": 800,
        "ncyc": 200,
        "force_bonds": 50.0,
        "force_proton": 20.0,
        "force_carbonyl": 30.0,
    }
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
        "anchor": manifest["task"]["anchor"],
        "step_scale": scale,
        "before_geometry": geometry,
        "target": targets,
        "input_restart_sha256": manifest["accepted_restart_sha256"],
        "restrained": True,
        "shooting_forbidden": True,
    }
    write_json(output / "WINDOW_MANIFEST.json", prepared)
    return prepared


def audit_window(root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest_path = root / "MANIFEST.json"
    manifest = load_json(manifest_path)
    prepared = load_json(output / "WINDOW_MANIFEST.json")
    restart = scratch / "stage.rst7"
    technical, geometry, diagnostics = BASE.TETRA._technical(
        scratch / "stage.out", restart, manifest
    )
    target = prepared["target"]
    before = prepared["before_geometry"]
    residuals = {}
    response = {}
    anchor = manifest["task"]["anchor"]
    desired = BASE.ANCHOR_TARGETS[anchor]
    if technical:
        residuals = {
            "attack_A": abs(float(geometry["attack_A"]) - target["attack_A"]),
            "cn_A": abs(float(geometry["c12_n3_A"]) - target["cn_A"]),
            "c12_o2_A": abs(float(geometry["c12_o2_A"]) - target["c12_o2_A"]),
            "nalpha_hg1_A": abs(float(geometry["nalpha_hg1_A"]) - target["nalpha_hg1_A"]),
            "hg1_n3_A": abs(float(geometry["hg1_n3_A"]) - target["hg1_n3_A"]),
        }
        primary_key = "attack_A" if anchor == "I1" else "c12_n3_A"
        target_key = "attack_A" if anchor == "I1" else "cn_A"
        before_gap = abs(float(before[primary_key]) - float(desired[target_key]))
        after_gap = abs(float(geometry[primary_key]) - float(desired[target_key]))
        requested = abs(float(before[primary_key]) - float(target[target_key]))
        improvement = before_gap - after_gap
        response = {
            "primary_coordinate": primary_key,
            "requested_A": requested,
            "improvement_A": improvement,
            "minimum_required_A": max(0.005, 0.20 * requested),
            "pass": requested <= 1.0e-6 or improvement >= max(0.005, 0.20 * requested),
        }
    residual_pass = bool(
        technical
        and residuals["attack_A"] <= 0.15
        and residuals["cn_A"] <= 0.15
        and residuals["c12_o2_A"] <= 0.06
        and residuals["nalpha_hg1_A"] <= 0.15
        and residuals["hg1_n3_A"] <= 0.15
    )
    chemical = BASE.REVERSE.chemical_guard(geometry) if technical else {"all": False}
    anchor_gate = BASE.anchor_gate(anchor, geometry) if technical else {"all": False}
    accepted = bool(technical and chemical.get("all") is True and residual_pass and response.get("pass") is True)
    persistent = output / "stage.rst7"
    if technical:
        shutil.copy2(restart, persistent)
    if accepted:
        shutil.copy2(persistent, root / "accepted.rst7")
        manifest["accepted_restart_sha256"] = sha256(root / "accepted.rst7")
        manifest["accepted_geometry"] = geometry
        manifest["adaptive_state"]["step_scale"] = min(
            1.0, float(manifest["adaptive_state"]["step_scale"]) * 1.25
        )
        if anchor_gate.get("all") is True:
            manifest["adaptive_state"]["anchor_reached"] = True
        targets_at_desired = all(
            abs(float(target[key]) - float(desired[key])) <= 1.0e-8
            for key in ("attack_A", "cn_A", "c12_o2_A")
        )
        if targets_at_desired and anchor_gate.get("all") is not True:
            manifest["adaptive_state"]["stalled_at_target"] = True
    elif technical:
        manifest["adaptive_state"]["step_scale"] = (
            float(manifest["adaptive_state"]["step_scale"]) * 0.5
        )
    else:
        manifest["technical_failure"] = True
    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_A1_I1_I2_RESPONSE_WINDOW"
            if technical else "NOT_EVALUATED_TECHNICAL_A1_I1_I2_RESPONSE_WINDOW"
        ),
        "technical_pass": technical,
        "window_index": prepared["window_index"],
        "anchor": anchor,
        "step_scale": prepared["step_scale"],
        "target": target,
        "before_geometry": before,
        "geometry": geometry,
        "target_residuals": residuals,
        "actual_response": response,
        "residual_pass": residual_pass,
        "chemical_guard": chemical,
        "anchor_gate": anchor_gate,
        "guard_pass": bool(chemical.get("all") is True),
        "accepted_for_inheritance": accepted,
        "input_restart_sha256": prepared["input_restart_sha256"],
        "output_restart": str(persistent) if technical else None,
        "output_restart_sha256": sha256(persistent) if persistent.is_file() else None,
        "restrained": True,
        "shooting_forbidden": True,
        "diagnostics": diagnostics,
    }
    manifest["windows"].append(result)
    terminal, reason = _terminal(manifest)
    manifest["terminal"] = terminal
    manifest["terminal_reason"] = reason
    write_json(output / "RESULT.json", result)
    write_json(output / ("PASS.json" if technical else "NOT_EVALUATED.json"), result)
    write_json(manifest_path, manifest)
    return result


def finalize(root: pathlib.Path) -> dict[str, Any]:
    manifest = load_json(root / "MANIFEST.json")
    terminal, reason = _terminal(manifest)
    if not terminal:
        raise ValueError("cannot finalize non-terminal route")
    windows = manifest["windows"]
    technical = bool(
        windows
        and not manifest["technical_failure"]
        and all(window.get("technical_pass") is True for window in windows)
    )
    inheritance = True
    previous = manifest["source"]["accepted_restart_sha256"]
    for window in windows:
        inheritance = inheritance and window.get("input_restart_sha256") == previous
        if window.get("accepted_for_inheritance"):
            previous = window.get("output_restart_sha256")
    reached = bool(technical and manifest["adaptive_state"].get("anchor_reached"))
    result = {
        "schema_version": 1,
        "status": TECHNICAL_PASS if technical else TECHNICAL_FAIL,
        "technical_complete": technical,
        "scientific_status": SCIENTIFIC_STATUS,
        "scientific_gate": (
            f"PASS_TRUE_{manifest['task']['anchor']}_RESTRAINED_ANCHOR"
            if reached else "PARTIAL_ACTUAL_RESPONSE_I1_I2_ANCHOR_BUILD"
        ),
        "task": manifest["task"],
        "windows_attempted": len(windows),
        "accepted_windows": sum(window.get("accepted_for_inheritance") is True for window in windows),
        "inheritance_sha256_verified": inheritance,
        "terminal_reason": reason,
        "final_geometry": manifest["accepted_geometry"],
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
    audit = output_root / "audit" / f"nylc_a1_i1_i2_response_builder_{array_job}.json"
    audit.parent.mkdir(parents=True, exist_ok=True)
    with audit.with_suffix(".lock").open("a+", encoding="utf-8") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        paths = [output_root / f"attempt_{array_job}_{index}" / "RESULT.json" for index in range(ARRAY_TASKS)]
        if not all(path.is_file() for path in paths):
            return False
        results = [load_json(path) for path in paths]
        technical = all(item.get("technical_complete") is True for item in results)
        reached = sum(str(item.get("scientific_gate", "")).startswith("PASS_TRUE_") for item in results)
        payload = {
            "schema_version": 1,
            "status": (
                "NOT_EVALUATED_TECHNICAL_A1_I1_I2_RESPONSE_BUILDER"
                if not technical else
                "PASS_TRUE_I1_I2_ANCHORS_REPRODUCED"
                if reached == ARRAY_TASKS else
                "PARTIAL_ACTUAL_RESPONSE_I1_I2_ANCHORS"
            ),
            "denominator_tasks": ARRAY_TASKS,
            "technical_pass_tasks": sum(item.get("technical_complete") is True for item in results),
            "true_anchor_pass_tasks": reached,
            "per_task": sorted(results, key=lambda item: item["task"]["task_index"]),
            "scientific_status": SCIENTIFIC_STATUS,
            "automatic_downstream_action": "NONE",
        }
        if audit.exists():
            if load_json(audit) != payload:
                raise FileExistsError(audit)
        else:
            write_json(audit, payload)
        return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=(
        "describe", "initialize", "prepare-window", "audit-window", "finalize", "merge-if-ready"
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
        prepare_window(args.root.resolve(), args.output.resolve(), args.scratch.resolve(), args.input_rst7.resolve(), args.window_index)
    elif args.mode == "audit-window":
        audit_window(args.root.resolve(), args.output.resolve(), args.scratch.resolve())
    elif args.mode == "finalize":
        finalize(args.root.resolve())
    else:
        merge_if_ready(args.output_root.resolve(), args.array_job)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
