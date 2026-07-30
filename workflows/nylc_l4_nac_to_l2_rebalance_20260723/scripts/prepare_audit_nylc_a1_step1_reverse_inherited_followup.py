#!/usr/bin/env python3
"""Continue the two qualified NylC A1 reverse-force sources with strict inheritance."""
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
TASK_ROOT = CAL.TASK_ROOT
PRMTOP = CAL.PRMTOP

ARRAY_TASKS = 2
MPI_RANKS = 8
MAX_WINDOWS = 8
OUTPUT_ROOT = TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_step1_reverse_inherited_followup"
SOURCE_ROOT = TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_step1_reverse_force_calibration"
SOURCE_ARRAY_JOB = "62441145"
PER_WINDOW_DELTAS_A = {
    "cn": -0.08,
    "carbonyl": 0.03,
    "nalpha_hg1": -0.05,
    "hg1_n3": 0.05,
}
SOURCES = (
    {
        "seed": "seed26737",
        "source_task_index": 6,
        "force_scale": 1.0,
        "route_label": "LOW_BIAS",
        "restart_sha256": "54ecabec3cbb993a81b36ce86630355e2143797481d71726ca83edd4d65862df",
    },
    {
        "seed": "seed26737",
        "source_task_index": 8,
        "force_scale": 4.0,
        "route_label": "TARGET_CLOSE",
        "restart_sha256": "2f30ffd0746cc7c53f36ac102b446256ea575e24a2d5bfa9da15a1455db96a7f",
    },
)
SCIENTIFIC_STATUS = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..1")
    source = SOURCES[index]
    return {
        "task_index": index,
        "seed": source["seed"],
        "source_array_job": SOURCE_ARRAY_JOB,
        "source_task_index": source["source_task_index"],
        "source_restart_sha256": source["restart_sha256"],
        "force_scale": source["force_scale"],
        "route_label": source["route_label"],
        "mechanism": "STEPWISE_REVERSE",
        "source_basin": "P",
        "destination": "REACTANT_SIDE",
    }


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "max_windows_per_chain": MAX_WINDOWS,
        "strict_serial_restart_inheritance": True,
        "failed_restart_inheritance": False,
        "active_groups": ["C12-N3", "C12-O2", "qPT"],
        "monitored_unrestrained_groups": ["OG1-C12"],
        "per_window_deltas_A": dict(PER_WINDOW_DELTAS_A),
        "sources": [task_spec(i) for i in range(ARRAY_TASKS)],
        "frozen_step1_contract": BASE.describe()["frozen_step1_contract"],
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }


def _qpt(geometry: Mapping[str, Any]) -> float:
    return float(geometry["nalpha_hg1_A"]) - float(geometry["hg1_n3_A"])


def _source_root(task: Mapping[str, Any]) -> pathlib.Path:
    return SOURCE_ROOT / f"attempt_{SOURCE_ARRAY_JOB}_{task['source_task_index']}"


def initialize(task_index: int, root: pathlib.Path, commit: str) -> dict[str, Any]:
    if root.exists():
        raise FileExistsError(root)
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("runtime is not bound to a full Git commit")
    task = task_spec(task_index)
    source_root = _source_root(task)
    source_result = BASE.read_json(source_root / "RESULT.json")
    source_manifest = BASE.read_json(source_root / "MANIFEST.json")
    restart = source_root / "stage.rst7"
    if not all((
        source_result.get("technical_complete") is True,
        source_result.get("eligible_for_inherited_chain") is True,
        restart.is_file(),
        CAL.sha256(restart) == task["source_restart_sha256"],
    )):
        raise ValueError("qualified reverse-force source authority failed")
    contract = CAL._validate_contract(source_manifest)
    if CAL.sha256(PRMTOP) != BASE.AUTH.BASE.EXPECTED_PRMTOP_SHA256:
        raise ValueError("frozen prmtop SHA changed")
    root.mkdir(parents=True, exist_ok=False)
    accepted = root / "accepted.rst7"
    shutil.copy2(restart, accepted)
    geometry = BASE.AUTH._geometry(accepted, {"qm_contract": contract})
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_REVERSE_INHERITED_FOLLOWUP",
        "github_commit": commit,
        "scientific_status": SCIENTIFIC_STATUS,
        "task": task,
        "source": {
            "attempt": str(source_root),
            "restart": str(restart),
            "restart_sha256": task["source_restart_sha256"],
            "response_gate": source_result["response_gate"],
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": CAL.sha256(PRMTOP),
        "qm_contract": contract,
        "source_geometry": geometry,
        "accepted_geometry": geometry,
        "accepted_restart_sha256": CAL.sha256(accepted),
        "windows": [],
        "technical_failure": False,
        "terminal": False,
        "terminal_reason": None,
        "strict_serial_restart_inheritance": True,
        "failed_restart_inheritance": False,
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }
    BASE.write_json(root / "CHAIN_MANIFEST.json", manifest)
    return manifest


def stage_spec(task: Mapping[str, Any], source: Mapping[str, Any], window_index: int) -> dict[str, Any]:
    scale = float(task["force_scale"])
    target = {
        "attack_A": float(source["attack_A"]),
        "c12_n3_A": float(source["c12_n3_A"]) + PER_WINDOW_DELTAS_A["cn"],
        "c12_o2_A": float(source["c12_o2_A"]) + PER_WINDOW_DELTAS_A["carbonyl"],
        "nalpha_hg1_A": float(source["nalpha_hg1_A"]) + PER_WINDOW_DELTAS_A["nalpha_hg1"],
        "hg1_n3_A": float(source["hg1_n3_A"]) + PER_WINDOW_DELTAS_A["hg1_n3"],
    }
    return {
        "window_index": int(window_index),
        "mechanism": "STEPWISE_REVERSE",
        "force_scale": scale,
        "targets": target,
        "attack_coordinate_active": False,
        "force_bond": CAL.BASE_FORCES["force_bond"] * scale,
        "force_carbonyl": CAL.BASE_FORCES["force_carbonyl"] * scale,
        "force_pt": CAL.BASE_FORCES["force_pt"] * scale,
        "maxcyc": 2200,
        "ncyc": 550,
    }


def _terminal(manifest: Mapping[str, Any]) -> bool:
    return bool(
        manifest.get("terminal")
        or manifest.get("technical_failure")
        or len(manifest["windows"]) >= MAX_WINDOWS
    )


def prepare_window(
    root: pathlib.Path,
    output: pathlib.Path,
    scratch: pathlib.Path,
    input_rst7: pathlib.Path,
    window_index: int,
) -> dict[str, Any]:
    manifest = BASE.read_json(root / "CHAIN_MANIFEST.json")
    if _terminal(manifest):
        raise ValueError("chain is terminal")
    if int(window_index) != len(manifest["windows"]):
        raise ValueError("window order is not strict")
    if CAL.sha256(input_rst7) != manifest["accepted_restart_sha256"]:
        raise ValueError("input breaks accepted restart SHA inheritance")
    output.mkdir(parents=True, exist_ok=False)
    scratch.mkdir(parents=True, exist_ok=False)
    stage = stage_spec(manifest["task"], manifest["accepted_geometry"], window_index)
    (scratch / "stage.in").write_text(
        CAL.BRIDGE.tetra_minimization_input(
            manifest["task"], stage, manifest["qm_contract"]["qmmask"]
        ),
        encoding="utf-8",
    )
    (scratch / "restraints.RST").write_text(CAL.restraints(stage), encoding="utf-8")
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


def _persist_small_files(output: pathlib.Path, scratch: pathlib.Path) -> None:
    for name in ("stage.in", "restraints.RST", "stage.mdinfo"):
        source = scratch / name
        if source.is_file():
            shutil.copy2(source, output / name)
    source = scratch / "stage.out"
    if source.is_file():
        tail = source.read_text(encoding="utf-8", errors="replace").splitlines()[-200:]
        (output / "ENGINE_TAIL.txt").write_text("\n".join(tail) + "\n", encoding="utf-8")


def audit_window(root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest_path = root / "CHAIN_MANIFEST.json"
    manifest = BASE.read_json(manifest_path)
    prepared = BASE.read_json(output / "WINDOW_MANIFEST.json")
    stage = prepared["stage"]
    technical, geometry, diagnostics = BASE.AUTH.TETRA._technical(
        scratch / "stage.out", scratch / "stage.rst7", manifest
    )
    previous = manifest["accepted_geometry"]
    target = stage["targets"]
    response = {"all": False}
    guard = {"all": False}
    if technical:
        response = {
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
            "active_coordinates": ["cn", "carbonyl", "pt"],
        }
        response["all"] = all(response[name]["pass"] for name in response["active_coordinates"])
        guard = CAL.chemical_guard(geometry)
        shutil.copy2(scratch / "stage.rst7", output / "stage.rst7")
    accepted = bool(technical and response["all"] and guard["all"])
    if accepted:
        shutil.copy2(output / "stage.rst7", root / "accepted.rst7")
        manifest["accepted_restart_sha256"] = CAL.sha256(root / "accepted.rst7")
        manifest["accepted_geometry"] = geometry
    else:
        manifest["terminal"] = True
        manifest["terminal_reason"] = (
            "NOT_EVALUATED_TECHNICAL" if not technical
            else "ACTUAL_RESPONSE_OR_GUARD_FAILED"
        )
    if not technical:
        manifest["technical_failure"] = True
    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_A1_REVERSE_INHERITED_WINDOW"
            if technical else "NOT_EVALUATED_TECHNICAL_A1_REVERSE_INHERITED_WINDOW"
        ),
        "technical_pass": technical,
        "window_index": stage["window_index"],
        "source_geometry": previous,
        "target_geometry": target,
        "geometry": geometry if technical else None,
        "source_qPT_A": _qpt(previous),
        "target_qPT_A": target["nalpha_hg1_A"] - target["hg1_n3_A"],
        "qPT_A": _qpt(geometry) if technical else None,
        "response_gate": response,
        "chemical_guard": guard,
        "accepted_for_inheritance": accepted,
        "input_restart_sha256": prepared["input_restart_sha256"],
        "output_restart_sha256": (
            CAL.sha256(output / "stage.rst7") if (output / "stage.rst7").is_file() else None
        ),
        "diagnostics": diagnostics,
        "restrained": True,
        "shooting_forbidden": True,
    }
    manifest["windows"].append(result)
    if len(manifest["windows"]) >= MAX_WINDOWS:
        manifest["terminal"] = True
        manifest["terminal_reason"] = "MAXIMUM_WINDOWS"
    BASE.write_json(output / "RESULT.json", result)
    BASE.write_json(output / ("PASS.json" if technical else "NOT_EVALUATED.json"), result)
    BASE.write_json(manifest_path, manifest)
    _persist_small_files(output, scratch)
    return result


def _write_hashes(root: pathlib.Path) -> None:
    rows = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.name != "SHA256.tsv" and path.suffix in (".json", ".rst7"):
            rows.append(f"{CAL.sha256(path)}  {path.relative_to(root)}")
    (root / "SHA256.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")


def finalize(root: pathlib.Path) -> dict[str, Any]:
    manifest = BASE.read_json(root / "CHAIN_MANIFEST.json")
    if not _terminal(manifest):
        raise ValueError("cannot finalize non-terminal chain")
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
    gate = (
        "NOT_EVALUATED_TECHNICAL_A1_REVERSE_INHERITED_FOLLOWUP"
        if not technical else
        "PASS_EIGHT_INHERITED_REVERSE_WINDOWS"
        if accepted == MAX_WINDOWS else
        "PARTIAL_INHERITED_REVERSE_RESPONSE"
        if accepted else
        "FAIL_NO_ADDITIONAL_REVERSE_RESPONSE"
    )
    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_A1_REVERSE_INHERITED_FOLLOWUP"
            if technical else "NOT_EVALUATED_TECHNICAL_A1_REVERSE_INHERITED_FOLLOWUP"
        ),
        "technical_complete": technical,
        "scientific_status": SCIENTIFIC_STATUS,
        "scientific_gate": gate,
        "task": manifest["task"],
        "windows_attempted": len(windows),
        "accepted_windows": accepted,
        "terminal_reason": manifest.get("terminal_reason"),
        "inheritance_sha256_verified": inheritance,
        "accepted_restart_sha256": manifest["accepted_restart_sha256"],
        "final_geometry": manifest["accepted_geometry"],
        "final_qPT_A": _qpt(manifest["accepted_geometry"]),
        "automatic_downstream_action": "NONE",
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
    technical = sum(item.get("technical_complete") is True for item in results)
    payload = {
        "schema_version": 1,
        "status": (
            "NOT_EVALUATED_TECHNICAL_A1_REVERSE_INHERITED_FOLLOWUP"
            if technical != ARRAY_TASKS else
            "PASS_REVERSE_INHERITED_FOLLOWUP_MATRIX"
        ),
        "denominator_tasks": ARRAY_TASKS,
        "technical_pass_tasks": technical,
        "per_task": sorted(results, key=lambda item: item["task"]["task_index"]),
        "scientific_status": SCIENTIFIC_STATUS,
        "automatic_downstream_action": "NONE",
    }
    audit = output_root / "audit" / f"nylc_a1_reverse_inherited_followup_{array_job}.json"
    audit.parent.mkdir(parents=True, exist_ok=True)
    if audit.exists() and BASE.read_json(audit) != payload:
        raise FileExistsError(audit)
    if not audit.exists():
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
