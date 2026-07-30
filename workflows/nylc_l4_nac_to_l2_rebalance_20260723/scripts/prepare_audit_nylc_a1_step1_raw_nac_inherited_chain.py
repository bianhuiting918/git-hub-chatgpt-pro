#!/usr/bin/env python3
"""Strict inherited Step1 chains from matched raw-NAC first-window calibration."""
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
RAW_PATH = HERE / "prepare_audit_nylc_a1_step1_raw_nac_forward_calibration.py"
_SPEC = importlib.util.spec_from_file_location("_a1_raw_nac_calibration", RAW_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"cannot import {RAW_PATH}")
RAW = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(RAW)

BASE = RAW.BASE
REACTIVE = RAW.REACTIVE
TASK_ROOT = RAW.TASK_ROOT
PRMTOP = RAW.PRMTOP
CALIBRATION_JOB = "62471131"
CALIBRATION_ROOT = RAW.OUTPUT_ROOT
CALIBRATION_AUDIT = (
    CALIBRATION_ROOT / "audit"
    / f"nylc_a1_raw_nac_forward_calibration_{CALIBRATION_JOB}.json"
)
OUTPUT_ROOT = (
    TASK_ROOT / "a1_activated_nac_20260726/qmmm"
    / "a1_step1_raw_nac_inherited_chain"
)
ARRAY_TASKS = 4
MAX_WINDOWS = 16
INITIAL_STEP = 1.0
MIN_STEP = 0.125
MAX_STEP = 4.0
MIN_ACTUAL_A = 0.005
MIN_RESPONSE_FRACTION = 0.10
AUTO_RELEASE = False
AUTO_SHOOTING = False
AUTO_PMF = False
SCIENTIFIC_STATUS = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"
BRANCHES = (
    (26723, "ADDITION_FIRST_RAW_NAC"),
    (26723, "FULLY_CONCERTED_RAW_NAC"),
    (26737, "ADDITION_FIRST_RAW_NAC"),
    (26737, "FULLY_CONCERTED_RAW_NAC"),
)


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: pathlib.Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def _source_task_index(seed: int, mode: str, scale: int) -> int:
    if scale not in RAW.FORCE_SCALES:
        raise ValueError("selected force scale is outside the calibrated matrix")
    base = 0 if int(seed) == 26723 else 9 if int(seed) == 26737 else None
    if base is None:
        raise ValueError("unknown seed")
    offset = RAW.FORCE_SCALES.index(scale)
    if mode == "ADDITION_FIRST_RAW_NAC":
        return base + 1 + offset
    if mode == "FULLY_CONCERTED_RAW_NAC":
        return base + 5 + offset
    raise ValueError("unknown calibrated mechanism")


def selection_from_audit(payload: Mapping[str, Any], task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..3")
    seed, mode = BRANCHES[index]
    matches = [
        row for row in payload.get("selected_minimum_force_scales", [])
        if int(row.get("seed", -1)) == seed and row.get("mode") == mode
    ]
    if len(matches) != 1 or matches[0].get("selected_weakest_scale") is None:
        raise ValueError(f"calibration has no unique selected scale for {seed}/{mode}")
    scale = int(matches[0]["selected_weakest_scale"])
    return {
        "task_index": index,
        "seed": seed,
        "seed_index": 0 if seed == 26723 else 1,
        "mode": mode,
        "mechanism": mode,
        "scale": scale,
        "force_scale": scale,
        "source_calibration_task_index": _source_task_index(seed, mode, scale),
        "source_calibration_job": CALIBRATION_JOB,
    }


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": 8,
        "calibration_job": CALIBRATION_JOB,
        "branches": [{"seed": seed, "mode": mode} for seed, mode in BRANCHES],
        "maximum_windows_per_chain": MAX_WINDOWS,
        "initial_step_multiplier": INITIAL_STEP,
        "minimum_step_multiplier": MIN_STEP,
        "maximum_step_multiplier": MAX_STEP,
        "strict_serial_restart_inheritance": True,
        "failed_restart_inheritance": False,
        "selected_weakest_calibrated_force_is_fixed_within_chain": True,
        "frozen_step1_contract": dict(RAW.QM_CONTRACT),
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }


def _qpt(geometry: Mapping[str, Any]) -> float:
    if "qPT_A" in geometry:
        return float(geometry["qPT_A"])
    return float(geometry["nalpha_hg1_A"]) - float(geometry["hg1_n3_A"])


def propose_targets(
    mode: str, geometry: Mapping[str, Any], step: float
) -> dict[str, float]:
    multiplier = float(step)
    if multiplier <= 0.0:
        raise ValueError("step multiplier must be positive")
    target = {
        "attack_A": float(geometry["attack_A"]) + RAW.FIRST_WINDOW_DELTAS_A["attack"] * multiplier,
        "c12_n3_A": float(geometry["c12_n3_A"]),
        "c12_o2_A": float(geometry["c12_o2_A"]) + RAW.FIRST_WINDOW_DELTAS_A["carbonyl"] * multiplier,
        "nalpha_hg1_A": float(geometry["nalpha_hg1_A"]),
        "hg1_n3_A": float(geometry["hg1_n3_A"]),
    }
    if mode == "FULLY_CONCERTED_RAW_NAC":
        target["c12_n3_A"] += RAW.FIRST_WINDOW_DELTAS_A["cn"] * multiplier
        target["nalpha_hg1_A"] += RAW.FIRST_WINDOW_DELTAS_A["nalpha_hg1"] * multiplier
        target["hg1_n3_A"] += RAW.FIRST_WINDOW_DELTAS_A["hg1_n3"] * multiplier
    elif mode != "ADDITION_FIRST_RAW_NAC":
        raise ValueError("unknown mechanism")
    return target


def active_coordinates(mode: str) -> list[str]:
    if mode == "ADDITION_FIRST_RAW_NAC":
        return ["attack", "carbonyl"]
    if mode == "FULLY_CONCERTED_RAW_NAC":
        return ["attack", "cn", "carbonyl", "pt"]
    raise ValueError("unknown mechanism")


def stage_spec(manifest: Mapping[str, Any], window_index: int) -> dict[str, Any]:
    task = manifest["task"]
    scale = int(task["scale"])
    return {
        "window_index": int(window_index),
        "mode": task["mode"],
        "mechanism": task["mode"],
        "force_scale": scale,
        "step_multiplier": float(manifest["state"]["step"]),
        "targets": propose_targets(
            task["mode"], manifest["accepted_geometry"], manifest["state"]["step"]
        ),
        "active_coordinates": active_coordinates(task["mode"]),
        "force_bond": RAW.FWD.CAL.BASE_FORCES["force_bond"] * scale,
        "force_carbonyl": RAW.FWD.CAL.BASE_FORCES["force_carbonyl"] * scale,
        "force_pt": RAW.FWD.CAL.BASE_FORCES["force_pt"] * scale,
        "maxcyc": 2200,
        "ncyc": 550,
    }


def restraints(spec: Mapping[str, Any], stage: Mapping[str, Any]) -> str:
    target = stage["targets"]
    rows = [
        RAW.FWD.CAL.BRIDGE._tight_distance_restraint(
            REACTIVE["og1"], REACTIVE["c12"], target["attack_A"], stage["force_bond"]
        ),
        RAW.FWD.CAL.BRIDGE._tight_distance_restraint(
            REACTIVE["c12"], REACTIVE["o2"], target["c12_o2_A"],
            stage["force_carbonyl"],
        ),
    ]
    if spec["mode"] == "FULLY_CONCERTED_RAW_NAC":
        rows.extend((
            RAW.FWD.CAL.BRIDGE._tight_distance_restraint(
                REACTIVE["c12"], REACTIVE["n3"], target["c12_n3_A"],
                stage["force_bond"],
            ),
            RAW.FWD.CAL.BRIDGE._tight_distance_restraint(
                REACTIVE["nalpha"], REACTIVE["hg1"], target["nalpha_hg1_A"],
                stage["force_pt"],
            ),
            RAW.FWD.CAL.BRIDGE._tight_distance_restraint(
                REACTIVE["hg1"], REACTIVE["n3"], target["hg1_n3_A"],
                stage["force_pt"],
            ),
        ))
    return "".join(rows)


def minimization_input(
    spec: Mapping[str, Any], stage: Mapping[str, Any], qmmask: str
) -> str:
    text = RAW.FWD.CAL.BRIDGE.tetra_minimization_input(spec, stage, qmmask)
    if "drms=0.01" not in text:
        raise ValueError("tight minimization drms authority changed")
    return text


def _coordinate(geometry: Mapping[str, Any], name: str) -> float:
    if name == "attack":
        return float(geometry["attack_A"])
    if name == "cn":
        return float(geometry["c12_n3_A"])
    if name == "carbonyl":
        return float(geometry["c12_o2_A"])
    if name == "pt":
        return _qpt(geometry)
    raise ValueError(name)


def actual_response_gate(
    previous: Mapping[str, Any], current: Mapping[str, Any],
    target: Mapping[str, Any], active: list[str]
) -> dict[str, Any]:
    checks: dict[str, Any] = {}
    for name in active:
        old = _coordinate(previous, name)
        new = _coordinate(current, name)
        goal = _coordinate(target, name)
        expected = goal - old
        actual = new - old
        same_direction = expected * actual > 0.0
        at_target = abs(new - goal) <= 0.01
        minimum = max(MIN_ACTUAL_A, MIN_RESPONSE_FRACTION * abs(expected))
        passed = same_direction and (at_target or abs(actual) >= minimum)
        checks[name] = {
            "expected_delta_A": expected,
            "actual_delta_A": actual,
            "minimum_required_actual_A": minimum,
            "same_direction": same_direction,
            "at_target": at_target,
            "pass": passed,
        }
    return {"active_coordinates": list(active), "checks": checks,
            "all": bool(active) and all(row["pass"] for row in checks.values())}


def restrained_hint(mode: str, geometry: Mapping[str, Any]) -> dict[str, Any]:
    attack = float(geometry["attack_A"])
    cn = float(geometry["c12_n3_A"])
    carbonyl = float(geometry["c12_o2_A"])
    qpt = _qpt(geometry)
    angle = float(geometry.get("attack_angle_deg", 0.0))
    oop = max(
        float(geometry.get("product_out_of_plane_A", 0.0)),
        float(geometry.get("c12_reactant_plane_out_of_plane_A", 0.0)),
    )
    angle_sum = float(geometry.get("product_angle_sum_deg", 360.0))
    tetra_response = carbonyl >= 1.28 or oop >= 0.15 or angle_sum <= 350.0
    if mode == "ADDITION_FIRST_RAW_NAC":
        checks = {
            "attack_bond_formed": 1.30 <= attack <= 1.75,
            "cn_bond_retained": 1.25 <= cn <= 1.90,
            "carbonyl_lengthened": 1.28 <= carbonyl <= 1.50,
            "proton_not_product_side": qpt <= 0.20,
            "attack_angle_nac": 90.0 <= angle <= 130.0,
            "tetrahedral_response": tetra_response,
        }
    elif mode == "FULLY_CONCERTED_RAW_NAC":
        checks = {
            "attack_boundary": 1.45 <= attack <= 2.20,
            "cn_boundary": 1.45 <= cn <= 2.20,
            "bond_balance": abs(attack - cn) <= 0.60,
            "pt_boundary": abs(qpt) <= 0.50,
            "attack_angle_nac": 90.0 <= angle <= 130.0,
            "tetrahedral_response": tetra_response,
        }
    else:
        raise ValueError("unknown mechanism")
    checks["all"] = all(checks.values())
    return checks


def initialize(task_index: int, root: pathlib.Path, commit: str) -> dict[str, Any]:
    if root.exists():
        raise FileExistsError(root)
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("runtime is not bound to a full Git commit")
    if not CALIBRATION_AUDIT.is_file():
        raise FileNotFoundError(CALIBRATION_AUDIT)
    audit = read_json(CALIBRATION_AUDIT)
    task = selection_from_audit(audit, task_index)
    source_root = CALIBRATION_ROOT / (
        f"attempt_{CALIBRATION_JOB}_{task['source_calibration_task_index']}"
    )
    source_result = read_json(source_root / "RESULT.json")
    source_manifest = read_json(source_root / "MANIFEST.json")
    restart = source_root / "stage.rst7"
    checks = (
        source_result.get("technical_complete") is True,
        source_result.get("eligible_by_absolute_response") is True,
        int(source_result["task"]["seed"]) == task["seed"],
        source_result["task"]["mode"] == task["mode"],
        int(source_result["task"]["scale"]) == task["scale"],
        restart.is_file(),
        sha256(PRMTOP) == RAW.PRMTOP_SHA256,
    )
    audit_rows = [
        row for row in audit.get("per_task", [])
        if int(row["task"]["task_index"]) == task["source_calibration_task_index"]
    ]
    if (
        not all(checks) or len(audit_rows) != 1
        or audit_rows[0].get("eligible_for_inherited_chain") is not True
    ):
        raise ValueError("selected calibration source authority failed")
    contract = RAW.validate_full_contract(
        source_manifest, source_manifest["qm_contract"]["qmmask"]
    )
    root.mkdir(parents=True, exist_ok=False)
    accepted = root / "accepted.rst7"
    shutil.copy2(restart, accepted)
    geometry = BASE.AUTH._geometry(accepted, {"qm_contract": contract})
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_RAW_NAC_INHERITED_CHAIN",
        "github_commit": commit,
        "scientific_status": SCIENTIFIC_STATUS,
        "task": task,
        "source": {
            "attempt": str(source_root),
            "result": str(source_root / "RESULT.json"),
            "restart": str(restart),
            "restart_sha256": sha256(restart),
            "calibration_audit": str(CALIBRATION_AUDIT),
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": sha256(PRMTOP),
        "qm_contract": contract,
        "accepted_restart": str(accepted),
        "accepted_restart_sha256": sha256(accepted),
        "accepted_geometry": geometry,
        "state": {"step": INITIAL_STEP},
        "windows": [],
        "restrained_hints": [],
        "technical_failure": False,
        "terminal": False,
        "terminal_reason": None,
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }
    write_json(root / "MANIFEST.json", manifest)
    return manifest


def prepare_window(
    root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path,
    input_rst7: pathlib.Path, window_index: int
) -> dict[str, Any]:
    manifest = read_json(root / "MANIFEST.json")
    if manifest.get("terminal"):
        raise ValueError(f"chain is terminal: {manifest.get('terminal_reason')}")
    if int(window_index) != len(manifest["windows"]):
        raise ValueError("window order is not strict")
    if sha256(input_rst7) != manifest["accepted_restart_sha256"]:
        raise ValueError("window input breaks accepted restart SHA inheritance")
    if output.exists() or scratch.exists():
        raise FileExistsError("window output or scratch already exists")
    stage = stage_spec(manifest, window_index)
    output.mkdir(parents=True, exist_ok=False)
    scratch.mkdir(parents=True, exist_ok=False)
    (scratch / "stage.in").write_text(
        minimization_input(manifest["task"], stage, manifest["qm_contract"]["qmmask"]),
        encoding="utf-8",
    )
    (scratch / "restraints.RST").write_text(
        restraints(manifest["task"], stage), encoding="utf-8"
    )
    prepared = {
        "schema_version": 1,
        "stage": stage,
        "input_restart": str(input_rst7),
        "input_restart_sha256": manifest["accepted_restart_sha256"],
        "strict_serial_restart_inheritance": True,
        "failed_restart_inheritance": False,
        "restrained": True,
        "shooting_forbidden": True,
    }
    write_json(output / "WINDOW_MANIFEST.json", prepared)
    return prepared


def audit_window(
    root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path
) -> dict[str, Any]:
    manifest_path = root / "MANIFEST.json"
    manifest = read_json(manifest_path)
    prepared = read_json(output / "WINDOW_MANIFEST.json")
    stage = prepared["stage"]
    restart = scratch / "stage.rst7"
    technical, geometry, diagnostics = BASE.AUTH.TETRA._technical(
        scratch / "stage.out", restart, manifest
    )
    previous = manifest["accepted_geometry"]
    response = {"all": False, "active_coordinates": [], "checks": {}}
    guard = {"pass": False, "checks": {}}
    hint = {"all": False}
    persistent = output / "stage.rst7"
    if technical:
        shutil.copy2(restart, persistent)
        response = actual_response_gate(
            previous, geometry, stage["targets"], stage["active_coordinates"]
        )
        guard = RAW.first_window_guard(geometry)
        hint = restrained_hint(manifest["task"]["mode"], geometry)
    accepted = bool(technical and response["all"] and guard["pass"])
    hint_pass = bool(accepted and hint["all"])
    if accepted:
        shutil.copy2(persistent, root / "accepted.rst7")
        manifest["accepted_restart_sha256"] = sha256(root / "accepted.rst7")
        manifest["accepted_geometry"] = geometry
        manifest["state"]["step"] = min(
            MAX_STEP, float(manifest["state"]["step"]) * 1.5
        )
    else:
        manifest["state"]["step"] = float(manifest["state"]["step"]) / 2.0
    if hint_pass:
        manifest["restrained_hints"].append({
            "window_index": stage["window_index"],
            "restart": str(persistent),
            "restart_sha256": sha256(persistent),
            "geometry": geometry,
            "gate": hint,
        })
        manifest["terminal"] = True
        manifest["terminal_reason"] = "RESTRAINED_HINT_REACHED"
    if not technical:
        manifest["technical_failure"] = True
        manifest["terminal"] = True
        manifest["terminal_reason"] = "NOT_EVALUATED_TECHNICAL"
    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_A1_RAW_NAC_INHERITED_WINDOW"
            if technical else
            "NOT_EVALUATED_TECHNICAL_A1_RAW_NAC_INHERITED_WINDOW"
        ),
        "scientific_status": SCIENTIFIC_STATUS,
        "window_index": stage["window_index"],
        "technical_pass": technical,
        "input_restart_sha256": prepared["input_restart_sha256"],
        "output_restart": str(persistent) if persistent.is_file() else None,
        "output_restart_sha256": sha256(persistent) if persistent.is_file() else None,
        "stage": stage,
        "previous_geometry": previous,
        "geometry": geometry if technical else None,
        "actual_response_gate": response,
        "chemical_guard": guard,
        "restrained_hint_gate": hint,
        "accepted_for_inheritance": accepted,
        "restrained_hint": hint_pass,
        "next_step_multiplier": manifest["state"]["step"],
        "diagnostics": diagnostics,
        "shooting_forbidden": True,
    }
    manifest["windows"].append(result)
    if (
        not manifest["terminal"]
        and len(manifest["windows"]) >= MAX_WINDOWS
    ):
        manifest["terminal"] = True
        manifest["terminal_reason"] = "MAXIMUM_WINDOWS"
    if (
        not manifest["terminal"]
        and float(manifest["state"]["step"]) < MIN_STEP
    ):
        manifest["terminal"] = True
        manifest["terminal_reason"] = "MINIMUM_RESPONSE_STEP_FAILED"
    write_json(output / "RESULT.json", result)
    write_json(output / ("PASS.json" if technical else "NOT_EVALUATED.json"), result)
    write_json(manifest_path, manifest)
    return result


def _write_hashes(root: pathlib.Path) -> None:
    rows = []
    for path in sorted(root.rglob("*")):
        if (
            path.is_file() and path.name != "SHA256.tsv"
            and path.suffix in (".json", ".rst7", ".in", ".RST", ".txt")
        ):
            rows.append(f"{sha256(path)}  {path.relative_to(root)}")
    (root / "SHA256.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")


def finalize(root: pathlib.Path) -> dict[str, Any]:
    manifest = read_json(root / "MANIFEST.json")
    if not manifest.get("terminal"):
        raise ValueError("cannot finalize a non-terminal chain")
    windows = manifest["windows"]
    technical = bool(
        windows and not manifest["technical_failure"]
        and all(row.get("technical_pass") is True for row in windows)
    )
    previous = manifest["source"]["restart_sha256"]
    inheritance = True
    for row in windows:
        inheritance = inheritance and row.get("input_restart_sha256") == previous
        if row.get("accepted_for_inheritance") is True:
            previous = row.get("output_restart_sha256")
    accepted_count = sum(row.get("accepted_for_inheritance") is True for row in windows)
    hints = manifest["restrained_hints"]
    gate = (
        "NOT_EVALUATED_TECHNICAL_A1_RAW_NAC_INHERITED_CHAIN"
        if not technical else
        "PASS_RESTRAINED_FORWARD_BOUNDARY_HINT"
        if hints else
        "PARTIAL_FORWARD_ACTUAL_RESPONSE"
        if accepted_count else
        "FAIL_NO_ACTUAL_FORWARD_RESPONSE"
    )
    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_A1_RAW_NAC_INHERITED_CHAIN"
            if technical else
            "NOT_EVALUATED_TECHNICAL_A1_RAW_NAC_INHERITED_CHAIN"
        ),
        "technical_complete": technical,
        "scientific_status": SCIENTIFIC_STATUS,
        "scientific_gate": gate,
        "task": manifest["task"],
        "terminal_reason": manifest["terminal_reason"],
        "windows_attempted": len(windows),
        "accepted_windows": accepted_count,
        "inheritance_sha256_verified": inheritance,
        "final_accepted_restart": manifest["accepted_restart"],
        "final_accepted_restart_sha256": manifest["accepted_restart_sha256"],
        "final_geometry": manifest["accepted_geometry"],
        "restrained_hint_count": len(hints),
        "restrained_hints": hints,
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }
    write_json(root / "RESULT.json", result)
    write_json(root / ("PASS.json" if technical else "NOT_EVALUATED.json"), result)
    _write_hashes(root)
    return result


def merge_if_ready(output_root: pathlib.Path, array_job: str) -> bool:
    paths = [
        output_root / f"attempt_{array_job}_{index}" / "RESULT.json"
        for index in range(ARRAY_TASKS)
    ]
    if not all(path.is_file() for path in paths):
        return False
    results = [read_json(path) for path in paths]
    technical = sum(row.get("technical_complete") is True for row in results)
    per_mode = {}
    for mode in RAW.MODES[1:]:
        rows = [row for row in results if row["task"]["mode"] == mode]
        per_mode[mode] = {
            "denominator_seeds": 2,
            "technical_pass_seeds": sum(row.get("technical_complete") is True for row in rows),
            "restrained_hint_seeds": sum(row.get("restrained_hint_count", 0) > 0 for row in rows),
        }
    payload = {
        "schema_version": 1,
        "status": (
            "NOT_EVALUATED_TECHNICAL_A1_RAW_NAC_INHERITED_CHAIN"
            if technical != ARRAY_TASKS else
            "PASS_BOTH_MECHANISMS_RESTRAINED_HINTS_REPRODUCED"
            if all(row["restrained_hint_seeds"] == 2 for row in per_mode.values()) else
            "PARTIAL_A1_RAW_NAC_INHERITED_CHAIN"
        ),
        "denominator_tasks": ARRAY_TASKS,
        "technical_pass_tasks": technical,
        "per_mechanism": per_mode,
        "per_task": results,
        "next_action": "REACTION_COORDINATE_FREE_RELEASE_OF_RESTRAINED_HINTS_ONLY",
        "scientific_status": SCIENTIFIC_STATUS,
        "automatic_downstream_action": "NONE",
    }
    audit = (
        output_root / "audit"
        / f"nylc_a1_raw_nac_inherited_chain_{array_job}.json"
    )
    if audit.exists() and read_json(audit) != payload:
        raise FileExistsError(audit)
    if not audit.exists():
        write_json(audit, payload)
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=(
        "describe", "initialize", "prepare-window", "audit-window",
        "finalize", "merge-if-ready",
    ))
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--root", type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--scratch", type=pathlib.Path)
    parser.add_argument("--input-rst7", type=pathlib.Path)
    parser.add_argument("--window-index", type=int)
    parser.add_argument("--github-commit", default="unknown")
    parser.add_argument("--output-root", type=pathlib.Path)
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
