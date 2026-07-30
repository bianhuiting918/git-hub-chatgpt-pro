#!/usr/bin/env python3
"""Build product-side Step1 tetrahedral bridges with an explicit carbonyl coordinate."""
from __future__ import annotations

import argparse
import fcntl
import importlib.util
import json
import math
import os
import pathlib
import re
import shutil
from typing import Any, Mapping

HERE = pathlib.Path(__file__).resolve().parent
CONT_PATH = HERE / "prepare_audit_nylc_a1_step1_2cv_continuation.py"
_SPEC = importlib.util.spec_from_file_location("_a1_2cv_cont", CONT_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"cannot import {CONT_PATH}")
CONT = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(CONT)
BASE = CONT.BASE
REACTIVE = BASE.REACTIVE
TASK_ROOT = BASE.TASK_ROOT
PRMTOP = BASE.PRMTOP

ARRAY_TASKS = 6
MPI_RANKS = 8
WINDOWS = 8
SCHEDULES = ("CARBONYL_FIRST", "PT_ASSISTED")
SOURCE_TASKS = (1, 3, 5)
SOURCE_ROOT = TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_step1_2cv_continuation"
OUTPUT_ROOT = TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_step1_tetrahedral_bridge"
TECHNICAL_PASS = "PASS_TECHNICAL_A1_TETRAHEDRAL_BRIDGE"
TECHNICAL_FAIL = "NOT_EVALUATED_TECHNICAL_A1_TETRAHEDRAL_BRIDGE"
SCIENTIFIC_STATUS = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"

SOURCES = {
    1: {
        "seed": "seed26723",
        "restart_sha256": "1ed5a5b6a72a43f9f56e41bf427bb95ffc96237c86811a34aec50fa3c89668f1",
    },
    3: {
        "seed": "seed26723",
        "restart_sha256": "09e22fb16670fc5fc6427a1aacf80851b97c323fa4eefb900760bde8ebc3be26",
    },
    5: {
        "seed": "seed26737",
        "restart_sha256": "740f5fb11de92842af2786ade6af03856dff40351209a613355a9b4dd19b3782",
    },
}
TETRA_TARGET = {
    "attack_A": 1.48,
    "c12_n3_A": 1.62,
    "c12_o2_A": 1.38,
    "nalpha_hg1_A": 1.05,
    "hg1_n3_A": 2.20,
}
FORCES = {
    "force_bond": 24.0,
    "force_carbonyl": 30.0,
    "force_pt": 18.0,
    "maxcyc": 2200,
    "ncyc": 550,
}


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..5")
    source_task = SOURCE_TASKS[index // 2]
    source = SOURCES[source_task]
    return {
        "task_index": index,
        "source_array_job": "62388613",
        "source_task_index": source_task,
        "seed": source["seed"],
        "source_restart_sha256": source["restart_sha256"],
        "source_basin": "P",
        "destination": "TETRAHEDRAL_BRIDGE",
        "schedule": SCHEDULES[index % 2],
    }


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "source_task_indices": list(SOURCE_TASKS),
        "schedules": list(SCHEDULES),
        "windows_per_chain": WINDOWS,
        "strict_serial_restart_inheritance": True,
        "failed_restart_inheritance": False,
        "explicit_carbonyl_restraint": True,
        "tetrahedral_target": dict(TETRA_TARGET),
        "frozen_step1_contract": BASE.describe()["frozen_step1_contract"],
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }


def _linear(start: float, end: float, progress: float) -> float:
    return float(start) + float(progress) * (float(end) - float(start))


def stage_targets(schedule: str, source: Mapping[str, Any], window_index: int) -> dict[str, Any]:
    if schedule not in SCHEDULES:
        raise ValueError(schedule)
    p = (int(window_index) + 1.0) / WINDOWS
    if schedule == "CARBONYL_FIRST":
        carbonyl_p = min(p / 0.50, 1.0)
        bond_p = max(0.0, (p - 0.125) / 0.875)
        pt_p = max(0.0, (p - 0.50) / 0.50)
    else:
        carbonyl_p = p
        bond_p = max(0.0, (p - 0.125) / 0.875)
        pt_p = min(p / 0.50, 1.0)
    target = {
        "attack_A": _linear(source["attack_A"], TETRA_TARGET["attack_A"], p),
        "c12_n3_A": _linear(source["c12_n3_A"], TETRA_TARGET["c12_n3_A"], bond_p),
        "c12_o2_A": _linear(source["c12_o2_A"], TETRA_TARGET["c12_o2_A"], carbonyl_p),
        "nalpha_hg1_A": _linear(
            source["nalpha_hg1_A"], TETRA_TARGET["nalpha_hg1_A"], pt_p
        ),
        "hg1_n3_A": _linear(source["hg1_n3_A"], TETRA_TARGET["hg1_n3_A"], pt_p),
    }
    return {
        "window_index": int(window_index),
        "progress": p,
        "schedule": schedule,
        "targets": target,
        "proton_coordinate_active": pt_p > 0.0,
        **FORCES,
    }


def _contract(source_manifest: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    return BASE._contract(source_manifest)


def initialize(task_index: int, root: pathlib.Path, commit: str) -> dict[str, Any]:
    if root.exists():
        raise FileExistsError(root)
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("runtime is not bound to a full Git commit")
    task = task_spec(task_index)
    source_root = SOURCE_ROOT / f"attempt_62388613_{task['source_task_index']}"
    result = BASE.read_json(source_root / "RESULT.json")
    source_manifest = BASE.read_json(source_root / "MANIFEST.json")
    restart = source_root / "accepted.rst7"
    checks = (
        result.get("technical_complete") is True,
        result.get("inheritance_sha256_verified") is True,
        int(result.get("accepted_windows", 0)) > 0,
        restart.is_file(),
        BASE.sha256(restart) == task["source_restart_sha256"],
    )
    if not all(checks):
        raise ValueError("fixed accepted tetrahedral-bridge source authority failed")
    contract, observed = _contract(source_manifest)
    if BASE.sha256(PRMTOP) != BASE.AUTH.BASE.EXPECTED_PRMTOP_SHA256:
        raise ValueError("frozen prmtop SHA changed")
    root.mkdir(parents=True, exist_ok=False)
    accepted = root / "accepted.rst7"
    shutil.copy2(restart, accepted)
    geometry = BASE.AUTH._geometry(accepted, {"qm_contract": contract})
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_TETRAHEDRAL_BRIDGE",
        "github_commit": commit,
        "scientific_status": SCIENTIFIC_STATUS,
        "task": task,
        "source": {
            "attempt": str(source_root),
            "restart": str(restart),
            "restart_sha256": task["source_restart_sha256"],
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": BASE.sha256(PRMTOP),
        "qm_contract": contract,
        "observed_frozen_contract": observed,
        "source_geometry": geometry,
        "accepted_geometry": geometry,
        "accepted_restart": str(accepted),
        "accepted_restart_sha256": BASE.sha256(accepted),
        "windows": [],
        "technical_failure": False,
        "terminal": False,
        "terminal_reason": None,
        "explicit_carbonyl_restraint": True,
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }
    BASE.write_json(root / "MANIFEST.json", manifest)
    return manifest


def restraints(stage: Mapping[str, Any]) -> str:
    target = stage["targets"]
    text = "".join(
        (
            BASE.AUTH.AC.distance_restraint(
                REACTIVE["og1"], REACTIVE["c12"], target["attack_A"], stage["force_bond"]
            ),
            BASE.AUTH.AC.distance_restraint(
                REACTIVE["c12"], REACTIVE["n3"], target["c12_n3_A"], stage["force_bond"]
            ),
            BASE.AUTH.AC.distance_restraint(
                REACTIVE["c12"], REACTIVE["o2"], target["c12_o2_A"],
                stage["force_carbonyl"],
            ),
        )
    )
    if stage["proton_coordinate_active"]:
        text += "".join(
            (
                BASE.AUTH.AC.distance_restraint(
                    REACTIVE["nalpha"], REACTIVE["hg1"],
                    target["nalpha_hg1_A"], stage["force_pt"],
                ),
                BASE.AUTH.AC.distance_restraint(
                    REACTIVE["hg1"], REACTIVE["n3"],
                    target["hg1_n3_A"], stage["force_pt"],
                ),
            )
        )
    return text


def _terminal(manifest: Mapping[str, Any]) -> tuple[bool, str | None]:
    if manifest.get("technical_failure"):
        return True, "NOT_EVALUATED_TECHNICAL"
    if manifest.get("terminal"):
        return True, manifest.get("terminal_reason")
    if len(manifest["windows"]) >= WINDOWS:
        return True, "MAXIMUM_WINDOWS"
    return False, None


def prepare_window(
    root: pathlib.Path,
    output: pathlib.Path,
    scratch: pathlib.Path,
    input_rst7: pathlib.Path,
    window_index: int,
) -> dict[str, Any]:
    manifest = BASE.read_json(root / "MANIFEST.json")
    terminal, reason = _terminal(manifest)
    if terminal:
        raise ValueError(f"chain is terminal: {reason}")
    if int(window_index) != len(manifest["windows"]):
        raise ValueError("window order is not strict")
    if BASE.sha256(input_rst7) != manifest["accepted_restart_sha256"]:
        raise ValueError("window input breaks accepted restart SHA inheritance")
    stage = stage_targets(
        manifest["task"]["schedule"], manifest["source_geometry"], window_index
    )
    output.mkdir(parents=True, exist_ok=False)
    scratch.mkdir(parents=True, exist_ok=False)
    (scratch / "stage.in").write_text(
        BASE.minimization_input(
            manifest["task"], stage, manifest["qm_contract"]["qmmask"]
        ),
        encoding="utf-8",
    )
    (scratch / "restraints.RST").write_text(restraints(stage), encoding="utf-8")
    prepared = {
        "schema_version": 1,
        "stage": stage,
        "input_restart": str(input_rst7),
        "input_restart_sha256": manifest["accepted_restart_sha256"],
        "strict_inheritance": True,
        "restrained": True,
        "shooting_forbidden": True,
    }
    BASE.write_json(output / "WINDOW_MANIFEST.json", prepared)
    return prepared


def _response(previous: float, current: float, target: float, tolerance: float = 0.06) -> dict[str, Any]:
    expected = float(target) - float(previous)
    actual = float(current) - float(previous)
    required = abs(expected) >= 0.005
    same_direction = expected * actual > 0.0
    at_target = abs(float(current) - float(target)) <= tolerance
    passed = (not required) or at_target or (same_direction and abs(actual) >= 0.005)
    return {
        "required": required,
        "expected_delta_A": expected,
        "actual_delta_A": actual,
        "same_direction": same_direction,
        "at_target": at_target,
        "pass": passed,
    }


def _qpt(geometry: Mapping[str, Any]) -> float:
    return float(geometry["nalpha_hg1_A"]) - float(geometry["hg1_n3_A"])


def chemical_guard(geometry: Mapping[str, Any]) -> dict[str, Any]:
    nearest = geometry.get("hg1_nearest_qm_heavy_atom")
    values = [
        geometry["attack_A"], geometry["c12_n3_A"], geometry["c12_o2_A"],
        geometry["nalpha_hg1_A"], geometry["hg1_n3_A"],
    ]
    checks = {
        "finite_geometry": all(math.isfinite(float(value)) for value in values),
        "attack_not_overcompressed": 1.25 <= float(geometry["attack_A"]) <= 1.80,
        "cn_not_overcompressed": float(geometry["c12_n3_A"]) >= 1.25,
        "carbonyl_safe": 1.15 <= float(geometry["c12_o2_A"]) <= 1.55,
        "proton_not_detached_or_misrouted": nearest in (
            REACTIVE["nalpha"], REACTIVE["n3"]
        ),
    }
    checks["all"] = all(checks.values())
    return checks


def tetrahedral_gate(geometry: Mapping[str, Any]) -> dict[str, Any]:
    angle_sum = float(geometry.get("product_angle_sum_deg", 360.0))
    product_oop = float(geometry.get("product_out_of_plane_A", 0.0))
    checks = {
        "attack_bond_retained": 1.30 <= float(geometry["attack_A"]) <= 1.70,
        "cn_bond_substantially_reformed": float(geometry["c12_n3_A"]) <= 1.90,
        "carbonyl_lengthened": 1.32 <= float(geometry["c12_o2_A"]) <= 1.48,
        "proton_moved_off_n3_side": _qpt(geometry) <= 0.20,
        "pyramidalized": angle_sum <= 350.0 or product_oop >= 0.15,
        "proton_not_misrouted": geometry.get("hg1_nearest_qm_heavy_atom")
        in (REACTIVE["nalpha"], REACTIVE["n3"]),
    }
    checks["all"] = all(checks.values())
    return checks


def audit_window(root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest_path = root / "MANIFEST.json"
    manifest = BASE.read_json(manifest_path)
    prepared = BASE.read_json(output / "WINDOW_MANIFEST.json")
    stage = prepared["stage"]
    restart = scratch / "stage.rst7"
    technical, geometry, diagnostics = BASE.AUTH.TETRA._technical(
        scratch / "stage.out", restart, manifest
    )
    persistent = output / "stage.rst7"
    if technical:
        shutil.copy2(restart, persistent)
    previous = manifest["accepted_geometry"]
    target = stage["targets"]
    response = {"all": False}
    guard = {"all": False}
    if technical:
        response = {
            "attack": _response(previous["attack_A"], geometry["attack_A"], target["attack_A"]),
            "cn": _response(previous["c12_n3_A"], geometry["c12_n3_A"], target["c12_n3_A"]),
            "carbonyl": _response(
                previous["c12_o2_A"], geometry["c12_o2_A"], target["c12_o2_A"]
            ),
            "pt": _response(_qpt(previous), _qpt(geometry),
                            target["nalpha_hg1_A"] - target["hg1_n3_A"]),
        }
        response["all"] = all(item["pass"] for key, item in response.items() if key != "all")
        guard = chemical_guard(geometry)
    accepted = bool(technical and response["all"] and guard["all"])
    if accepted:
        shutil.copy2(persistent, root / "accepted.rst7")
        manifest["accepted_restart_sha256"] = BASE.sha256(root / "accepted.rst7")
        manifest["accepted_geometry"] = geometry
    else:
        manifest["terminal"] = True
        manifest["terminal_reason"] = (
            "NOT_EVALUATED_TECHNICAL" if not technical
            else "ACTUAL_TETRAHEDRAL_RESPONSE_FAILED"
        )
    if not technical:
        manifest["technical_failure"] = True
    candidate = tetrahedral_gate(geometry) if technical else {"all": False}
    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_A1_TETRAHEDRAL_BRIDGE_WINDOW"
            if technical else "NOT_EVALUATED_TECHNICAL_A1_TETRAHEDRAL_BRIDGE_WINDOW"
        ),
        "technical_pass": technical,
        "window_index": stage["window_index"],
        "schedule": stage["schedule"],
        "targets": target,
        "geometry": geometry,
        "qPT_A": _qpt(geometry) if technical else None,
        "response_gate": response,
        "chemical_guard": guard,
        "tetrahedral_gate": candidate,
        "accepted_for_inheritance": accepted,
        "input_restart_sha256": prepared["input_restart_sha256"],
        "output_restart": str(persistent) if persistent.is_file() else None,
        "output_restart_sha256": BASE.sha256(persistent) if persistent.is_file() else None,
        "diagnostics": diagnostics,
        "restrained": True,
        "shooting_forbidden": True,
    }
    manifest["windows"].append(result)
    if len(manifest["windows"]) >= WINDOWS:
        manifest["terminal"] = True
        manifest["terminal_reason"] = "MAXIMUM_WINDOWS"
    BASE.write_json(output / "RESULT.json", result)
    BASE.write_json(
        output / ("PASS.json" if technical else "NOT_EVALUATED.json"), result
    )
    BASE.write_json(manifest_path, manifest)
    return result


def _write_hashes(root: pathlib.Path) -> None:
    rows = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.name != "SHA256.tsv" and path.suffix in (".json", ".rst7"):
            rows.append(f"{BASE.sha256(path)}  {path.relative_to(root)}")
    (root / "SHA256.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")


def finalize(root: pathlib.Path) -> dict[str, Any]:
    manifest = BASE.read_json(root / "MANIFEST.json")
    terminal, reason = _terminal(manifest)
    if not terminal:
        raise ValueError("cannot finalize a non-terminal chain")
    windows = manifest["windows"]
    technical = bool(
        windows
        and not manifest["technical_failure"]
        and all(item.get("technical_pass") is True for item in windows)
    )
    previous = manifest["source"]["restart_sha256"]
    inheritance = True
    for window in windows:
        inheritance = inheritance and window["input_restart_sha256"] == previous
        if window.get("accepted_for_inheritance"):
            previous = window["output_restart_sha256"]
    accepted = sum(item.get("accepted_for_inheritance") is True for item in windows)
    final_geometry = manifest["accepted_geometry"]
    candidate = tetrahedral_gate(final_geometry)
    gate = (
        "NOT_EVALUATED_TECHNICAL_A1_TETRAHEDRAL_BRIDGE"
        if not technical else
        "PASS_RESTRAINED_TETRAHEDRAL_BRIDGE_CANDIDATE"
        if candidate["all"] else
        "PARTIAL_TETRAHEDRALIZATION_RESPONSE"
        if accepted else
        "FAIL_NO_ACTUAL_TETRAHEDRAL_RESPONSE"
    )
    result = {
        "schema_version": 1,
        "status": TECHNICAL_PASS if technical else TECHNICAL_FAIL,
        "technical_complete": technical,
        "scientific_status": SCIENTIFIC_STATUS,
        "scientific_gate": gate,
        "task": manifest["task"],
        "terminal_reason": reason,
        "windows_attempted": len(windows),
        "accepted_windows": accepted,
        "inheritance_sha256_verified": inheritance,
        "final_geometry": final_geometry,
        "final_qPT_A": _qpt(final_geometry),
        "tetrahedral_gate": candidate,
        "accepted_restart_sha256": manifest["accepted_restart_sha256"],
        "restrained_structures_are_not_ts": True,
        "automatic_downstream_action": "NONE",
    }
    BASE.write_json(root / "RESULT.json", result)
    BASE.write_json(root / ("PASS.json" if technical else "NOT_EVALUATED.json"), result)
    _write_hashes(root)
    return result


def merge_if_ready(output_root: pathlib.Path, array_job: str) -> bool:
    audit = output_root / "audit" / f"nylc_a1_tetrahedral_bridge_{array_job}.json"
    audit.parent.mkdir(parents=True, exist_ok=True)
    with audit.with_suffix(".lock").open("a+", encoding="utf-8") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        paths = [
            output_root / f"attempt_{array_job}_{index}" / "RESULT.json"
            for index in range(ARRAY_TASKS)
        ]
        if not all(path.is_file() for path in paths):
            return False
        results = [BASE.read_json(path) for path in paths]
        technical = sum(item.get("technical_complete") is True for item in results)
        candidates = sum(
            item.get("tetrahedral_gate", {}).get("all") is True for item in results
        )
        payload = {
            "schema_version": 1,
            "status": (
                "NOT_EVALUATED_TECHNICAL_A1_TETRAHEDRAL_BRIDGE"
                if technical != ARRAY_TASKS else
                "PASS_RESTRAINED_TETRAHEDRAL_BRIDGE_CANDIDATES"
                if candidates else
                "FAIL_NO_RESTRAINED_TETRAHEDRAL_BRIDGE_CANDIDATE"
            ),
            "denominator_tasks": ARRAY_TASKS,
            "technical_pass_tasks": technical,
            "restrained_tetrahedral_candidate_tasks": candidates,
            "per_task": sorted(results, key=lambda item: item["task"]["task_index"]),
            "scientific_status": SCIENTIFIC_STATUS,
            "automatic_downstream_action": "NONE",
        }
        if audit.exists():
            if BASE.read_json(audit) != payload:
                raise FileExistsError(audit)
        else:
            BASE.write_json(audit, payload)
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
