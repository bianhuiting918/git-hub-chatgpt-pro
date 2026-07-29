#!/usr/bin/env python3
"""Bidirectional adaptive two-coordinate Step1 scout for NylC A1.

Eight parallel chains combine two seeds, two endpoint sources, and two
mechanistic schedules.  Within a chain, only accepted windows are inherited.
The bond-exchange pair and proton-transfer pair are coordinated surrogates for
two collective coordinates.  Carbonyl length and out-of-plane displacement are
observables, never restrained targets.  Restrained chains are not TS,
committor, PMF, barrier, or mechanism evidence.
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
AUTHORITY_PATH = HERE / "prepare_audit_nylc_a1_step1_four_anchor_bidirectional_gapfill.py"
_SPEC = importlib.util.spec_from_file_location("_a1_two_cv_authority", AUTHORITY_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"cannot import {AUTHORITY_PATH}")
AUTH = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(AUTH)

TASK_ROOT = AUTH.TASK_ROOT
PRMTOP = AUTH.PRMTOP
REACTIVE = AUTH.REACTIVE
MECHANISMS = ("ADDITION_FIRST", "CONCERTED")
SOURCE_BASINS = ("R", "P")
ARRAY_TASKS = 8
MPI_RANKS = 8
MAX_WINDOWS = 20
INITIAL_STEP = 0.125
MIN_STEP = 1.0 / 64.0
MIN_RESPONSE_FRACTION = 0.10
MIN_RESPONSE_A = 0.01
TECHNICAL_PASS = "PASS_TECHNICAL_A1_BIDIRECTIONAL_2CV"
TECHNICAL_FAIL = "NOT_EVALUATED_TECHNICAL_A1_BIDIRECTIONAL_2CV"
SCIENTIFIC_STATUS = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"

R_TARGET = {
    "attack_A": 2.18,
    "cn_A": 1.38,
    "nalpha_hg1_A": 1.04,
    "hg1_n3_A": 2.45,
}
P_TARGET = {
    "attack_A": 1.42,
    "cn_A": 3.10,
    "nalpha_hg1_A": 2.70,
    "hg1_n3_A": 1.02,
}


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: pathlib.Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def read_json(path: pathlib.Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..7")
    seed_index = index // 4
    local = index % 4
    source_basin = SOURCE_BASINS[local // 2]
    mechanism = MECHANISMS[local % 2]
    return {
        "task_index": index,
        "seed_index": seed_index,
        "seed": AUTH.REACTANT_SOURCES[seed_index]["seed"],
        "source_basin": source_basin,
        "destination_basin": "P" if source_basin == "R" else "R",
        "mechanism": mechanism,
    }


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "mechanisms": list(MECHANISMS),
        "source_basins": list(SOURCE_BASINS),
        "mapping": "task=seed*4+source_basin*2+mechanism",
        "maximum_windows_per_chain": MAX_WINDOWS,
        "initial_lambda_step": INITIAL_STEP,
        "minimum_lambda_step": MIN_STEP,
        "strict_serial_restart_inheritance": True,
        "failed_restart_inheritance": False,
        "coordinate_model": {
            "xi_bond": "d(OG1,C12)-d(C12,N3)",
            "xi_pt": "d(Nalpha,HG1)-d(HG1,N3)",
            "amber_implementation": "paired_distance_surrogate_with_shared_lambda",
            "carbonyl_and_oop": "observed_not_restrained",
        },
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


def _linear(a: float, b: float, progress: float) -> float:
    return float(a) + float(progress) * (float(b) - float(a))


def path_targets(mechanism: str, fraction: float) -> dict[str, float]:
    if mechanism not in MECHANISMS:
        raise ValueError(f"unknown mechanism {mechanism}")
    f = float(fraction)
    if not 0.0 <= f <= 1.0:
        raise ValueError("path fraction must be between 0 and 1")
    if math.isclose(f, 0.0, abs_tol=1.0e-12):
        return dict(R_TARGET)
    if math.isclose(f, 1.0, abs_tol=1.0e-12):
        return dict(P_TARGET)
    if mechanism == "CONCERTED":
        attack_progress = cn_progress = pt_progress = f
    else:
        attack_progress = min(f / 0.45, 1.0)
        cn_progress = 0.0 if f <= 0.30 else min((f - 0.30) / 0.70, 1.0)
        pt_progress = 0.0 if f <= 0.50 else min((f - 0.50) / 0.50, 1.0)
    return {
        "attack_A": _linear(R_TARGET["attack_A"], P_TARGET["attack_A"], attack_progress),
        "cn_A": _linear(R_TARGET["cn_A"], P_TARGET["cn_A"], cn_progress),
        "nalpha_hg1_A": _linear(
            R_TARGET["nalpha_hg1_A"], P_TARGET["nalpha_hg1_A"], pt_progress
        ),
        "hg1_n3_A": _linear(
            R_TARGET["hg1_n3_A"], P_TARGET["hg1_n3_A"], pt_progress
        ),
    }


def progress_coordinates(geometry: Mapping[str, Any]) -> dict[str, float]:
    return {
        "xi_bond_A": float(geometry["attack_A"]) - float(geometry["c12_n3_A"]),
        "xi_pt_A": float(geometry["nalpha_hg1_A"]) - float(geometry["hg1_n3_A"]),
    }


def _target_coordinates(mechanism: str, fraction: float) -> dict[str, float]:
    target = path_targets(mechanism, fraction)
    return {
        "xi_bond_A": target["attack_A"] - target["cn_A"],
        "xi_pt_A": target["nalpha_hg1_A"] - target["hg1_n3_A"],
    }


def _response_check(previous: float, current: float, old_target: float,
                    new_target: float, required: bool) -> dict[str, Any]:
    expected = float(new_target) - float(old_target)
    actual = float(current) - float(previous)
    same_direction = expected * actual > 0.0
    ratio = actual / expected if same_direction and not math.isclose(expected, 0.0) else 0.0
    at_target = abs(float(current) - float(new_target)) <= 0.10
    passed = (not required) or at_target or (
        same_direction
        and abs(actual) >= MIN_RESPONSE_A
        and ratio >= MIN_RESPONSE_FRACTION
    )
    return {
        "required": required,
        "expected_delta_A": expected,
        "actual_delta_A": actual,
        "response_fraction": ratio,
        "same_direction": same_direction,
        "at_proposed_target": at_target,
        "pass": passed,
    }


def actual_response_gate(
    task: Mapping[str, Any],
    previous_geometry: Mapping[str, Any],
    current_geometry: Mapping[str, Any],
    accepted_fraction: float,
    proposed_fraction: float,
) -> dict[str, Any]:
    previous = progress_coordinates(previous_geometry)
    current = progress_coordinates(current_geometry)
    old_target = _target_coordinates(task["mechanism"], accepted_fraction)
    new_target = _target_coordinates(task["mechanism"], proposed_fraction)
    bond = _response_check(
        previous["xi_bond_A"], current["xi_bond_A"],
        old_target["xi_bond_A"], new_target["xi_bond_A"], True,
    )
    pt_required = abs(new_target["xi_pt_A"] - old_target["xi_pt_A"]) > 1.0e-8
    proton = _response_check(
        previous["xi_pt_A"], current["xi_pt_A"],
        old_target["xi_pt_A"], new_target["xi_pt_A"], pt_required,
    )
    nearest = current_geometry.get("hg1_nearest_qm_heavy_atom")
    finite = all(
        math.isfinite(float(current_geometry[key]))
        for key in ("attack_A", "c12_n3_A", "nalpha_hg1_A", "hg1_n3_A", "c12_o2_A")
    )
    chemical_sanity = bool(
        finite
        and float(current_geometry["attack_A"]) >= 1.25
        and float(current_geometry["c12_n3_A"]) >= 1.25
        and 1.10 <= float(current_geometry["c12_o2_A"]) <= 1.60
        and nearest in (REACTIVE["nalpha"], REACTIVE["n3"])
    )
    wrong_direction = bool(
        (bond["required"] and not bond["same_direction"] and not bond["at_proposed_target"])
        or (
            proton["required"]
            and not proton["same_direction"]
            and not proton["at_proposed_target"]
        )
    )
    return {
        "previous_coordinates": previous,
        "current_coordinates": current,
        "old_target_coordinates": old_target,
        "new_target_coordinates": new_target,
        "bond_response": bond,
        "proton_response": proton,
        "chemical_sanity": chemical_sanity,
        "wrong_direction": wrong_direction,
        "all": bool(chemical_sanity and bond["pass"] and proton["pass"]),
    }


def _source_for(task: Mapping[str, Any]) -> dict[str, Any]:
    collection = (
        AUTH.REACTANT_SOURCES
        if task["source_basin"] == "R"
        else AUTH.PRODUCT_SOURCES
    )
    return dict(collection[int(task["seed_index"])])


def _contract(source_manifest: Mapping[str, Any]) -> tuple[dict[str, Any], dict[str, Any]]:
    contract = dict(source_manifest["qm_contract"])
    observed = {
        "qm_atom_count": contract.get("qm_atom_count"),
        "qmcharge": contract.get("qmcharge"),
        "electron_count": contract.get("electron_count_including_link_h"),
        "link_atom_count": contract.get("link_atom_count"),
        "qm_water_count": contract.get("step1_qm_water_count"),
    }
    if observed != describe()["frozen_step1_contract"]:
        raise ValueError(f"frozen Step1 contract mismatch: {observed}")
    return contract, observed


def _terminal(manifest: Mapping[str, Any]) -> tuple[bool, str | None]:
    state = manifest["state"]
    if manifest.get("technical_failure"):
        return True, "NOT_EVALUATED_TECHNICAL"
    if state.get("destination_reached"):
        return True, "DESTINATION_REACHED"
    if len(manifest["windows"]) >= MAX_WINDOWS:
        return True, "MAXIMUM_WINDOWS"
    if float(state["step"]) < MIN_STEP:
        return True, "MINIMUM_ACTUAL_RESPONSE_STEP_FAILED"
    return False, None


def initialize(task_index: int, root: pathlib.Path, commit: str) -> dict[str, Any]:
    if root.exists():
        raise FileExistsError(root)
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("runtime is not bound to a full Git commit")
    task = task_spec(task_index)
    source = _source_for(task)
    restart = pathlib.Path(source["restart"])
    if not restart.is_file() or sha256(restart) != source["restart_sha256"]:
        raise ValueError("fixed source restart SHA authority failed")
    source_manifest, qmmask, _ = AUTH.REVERSE.validate_authority(task["seed_index"])
    contract, observed = _contract(source_manifest)
    if sha256(PRMTOP) != AUTH.BASE.EXPECTED_PRMTOP_SHA256:
        raise ValueError("frozen prmtop SHA changed")
    root.mkdir(parents=True, exist_ok=False)
    accepted = root / "accepted.rst7"
    shutil.copy2(restart, accepted)
    geometry = AUTH._geometry(accepted, {"qm_contract": contract})
    initial_fraction = 0.0 if task["source_basin"] == "R" else 1.0
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_BIDIRECTIONAL_2CV",
        "github_commit": commit,
        "scientific_status": SCIENTIFIC_STATUS,
        "task": task,
        "source": {
            "restart": str(restart),
            "restart_sha256": source["restart_sha256"],
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": sha256(PRMTOP),
        "qm_contract": contract,
        "observed_frozen_contract": observed,
        "accepted_restart": str(accepted),
        "accepted_restart_sha256": sha256(accepted),
        "accepted_geometry": geometry,
        "state": {
            "accepted_fraction": initial_fraction,
            "destination_fraction": 1.0 - initial_fraction,
            "direction": 1.0 if initial_fraction == 0.0 else -1.0,
            "step": INITIAL_STEP,
            "destination_reached": False,
        },
        "windows": [],
        "technical_failure": False,
        "carbonyl_and_oop_are_observables_only": True,
        "restrained_structures_are_not_ts": True,
        "automatic_downstream_action": "NONE",
    }
    write_json(root / "MANIFEST.json", manifest)
    return manifest


def restraints(stage: Mapping[str, Any]) -> str:
    target = stage["targets"]
    text = "".join((
        AUTH.AC.distance_restraint(
            REACTIVE["og1"], REACTIVE["c12"], target["attack_A"], stage["force_bond"]
        ),
        AUTH.AC.distance_restraint(
            REACTIVE["c12"], REACTIVE["n3"], target["cn_A"], stage["force_bond"]
        ),
    ))
    if stage["proton_coordinate_active"]:
        text += "".join((
            AUTH.AC.distance_restraint(
                REACTIVE["nalpha"], REACTIVE["hg1"],
                target["nalpha_hg1_A"], stage["force_pt"],
            ),
            AUTH.AC.distance_restraint(
                REACTIVE["hg1"], REACTIVE["n3"],
                target["hg1_n3_A"], stage["force_pt"],
            ),
        ))
    return text


def minimization_input(task: Mapping[str, Any], stage: Mapping[str, Any],
                       qmmask: str) -> str:
    return f"""NylC A1 bidirectional 2CV task {task['task_index']} window {stage['window_index']}
&cntrl
  imin=1, ntmin=2, maxcyc={stage['maxcyc']}, ncyc={stage['ncyc']}, dx0=0.005,
  ntb=1, cut=10.0, ntpr=5, ntxo=1, ifqnt=1,
  ntr=1, nmropt=1, restraint_wt=1.0,
  restraintmask='{AUTH.BASE.NON_QM_SOLUTE_HEAVY_MASK}', drms=0.10,
/
{AUTH.BASE.qmmm_block(qmmask)}&wt type='END' /
DISANG=restraints.RST
DUMPAVE=restraint.dat
"""


def prepare_window(root: pathlib.Path, output: pathlib.Path,
                   scratch: pathlib.Path, input_rst7: pathlib.Path,
                   window_index: int) -> dict[str, Any]:
    manifest = read_json(root / "MANIFEST.json")
    terminal, reason = _terminal(manifest)
    if terminal:
        raise ValueError(f"chain is terminal: {reason}")
    if int(window_index) != len(manifest["windows"]):
        raise ValueError("window order is not strict")
    if sha256(input_rst7) != manifest["accepted_restart_sha256"]:
        raise ValueError("window input breaks accepted restart SHA inheritance")
    state = manifest["state"]
    old_fraction = float(state["accepted_fraction"])
    direction = float(state["direction"])
    proposed = old_fraction + direction * float(state["step"])
    destination = float(state["destination_fraction"])
    proposed = min(proposed, destination) if direction > 0 else max(proposed, destination)
    old_target = path_targets(manifest["task"]["mechanism"], old_fraction)
    target = path_targets(manifest["task"]["mechanism"], proposed)
    old_qpt = old_target["nalpha_hg1_A"] - old_target["hg1_n3_A"]
    new_qpt = target["nalpha_hg1_A"] - target["hg1_n3_A"]
    stage = {
        "window_index": int(window_index),
        "accepted_fraction": old_fraction,
        "proposed_fraction": proposed,
        "targets": target,
        "target_coordinates": _target_coordinates(
            manifest["task"]["mechanism"], proposed
        ),
        "proton_coordinate_active": abs(new_qpt - old_qpt) > 1.0e-8,
        "force_bond": 8.0,
        "force_pt": 6.0,
        "maxcyc": 800,
        "ncyc": 200,
    }
    output.mkdir(parents=True, exist_ok=False)
    scratch.mkdir(parents=True, exist_ok=False)
    (scratch / "stage.in").write_text(
        minimization_input(manifest["task"], stage, manifest["qm_contract"]["qmmask"]),
        encoding="utf-8",
    )
    (scratch / "restraints.RST").write_text(restraints(stage), encoding="utf-8")
    prepared = {
        "schema_version": 1,
        "stage": stage,
        "input_restart": str(input_rst7),
        "input_restart_sha256": manifest["accepted_restart_sha256"],
        "carbonyl_and_oop_are_observables_only": True,
        "restrained": True,
        "shooting_forbidden": True,
    }
    write_json(output / "WINDOW_MANIFEST.json", prepared)
    return prepared


def audit_window(root: pathlib.Path, output: pathlib.Path,
                 scratch: pathlib.Path) -> dict[str, Any]:
    manifest_path = root / "MANIFEST.json"
    manifest = read_json(manifest_path)
    prepared = read_json(output / "WINDOW_MANIFEST.json")
    stage = prepared["stage"]
    restart = scratch / "stage.rst7"
    technical, geometry, diagnostics = AUTH.TETRA._technical(
        scratch / "stage.out", restart, manifest
    )
    persistent = output / "stage.rst7"
    if technical:
        shutil.copy2(restart, persistent)
    response = (
        actual_response_gate(
            manifest["task"], manifest["accepted_geometry"], geometry,
            stage["accepted_fraction"], stage["proposed_fraction"],
        )
        if technical else {"all": False}
    )
    chemical = (
        AUTH.REVERSE.chemical_guard(geometry) if technical else {"all": False}
    )
    destination_basin = manifest["task"]["destination_basin"]
    at_destination = math.isclose(
        float(stage["proposed_fraction"]),
        float(manifest["state"]["destination_fraction"]),
        abs_tol=1.0e-12,
    )
    endpoint_gate = (
        AUTH.anchor_gate(destination_basin, geometry)
        if technical and at_destination else None
    )
    accepted = bool(
        technical
        and response.get("all") is True
        and chemical.get("all") is True
        and (endpoint_gate is None or endpoint_gate.get("all") is True)
    )
    if accepted:
        shutil.copy2(persistent, root / "accepted.rst7")
        manifest["accepted_restart_sha256"] = sha256(root / "accepted.rst7")
        manifest["accepted_geometry"] = geometry
        manifest["state"]["accepted_fraction"] = stage["proposed_fraction"]
        if at_destination:
            manifest["state"]["destination_reached"] = True
    else:
        manifest["state"]["step"] = float(manifest["state"]["step"]) / 2.0
    if not technical:
        manifest["technical_failure"] = True
    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_A1_BIDIRECTIONAL_2CV_WINDOW"
            if technical else
            "NOT_EVALUATED_TECHNICAL_A1_BIDIRECTIONAL_2CV_WINDOW"
        ),
        "technical_pass": technical,
        "window_index": stage["window_index"],
        "accepted_fraction": stage["accepted_fraction"],
        "proposed_fraction": stage["proposed_fraction"],
        "input_restart_sha256": prepared["input_restart_sha256"],
        "output_restart": str(persistent) if persistent.is_file() else None,
        "output_restart_sha256": sha256(persistent) if persistent.is_file() else None,
        "targets": stage["targets"],
        "target_coordinates": stage["target_coordinates"],
        "geometry": geometry,
        "actual_coordinates": progress_coordinates(geometry) if technical else None,
        "actual_response_gate": response,
        "chemical_guard": chemical,
        "destination_gate": endpoint_gate,
        "accepted_for_inheritance": accepted,
        "next_step": manifest["state"]["step"],
        "carbonyl_and_oop_are_observables_only": True,
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
    manifest = read_json(root / "MANIFEST.json")
    terminal, reason = _terminal(manifest)
    if not terminal:
        raise ValueError("cannot finalize a non-terminal chain")
    windows = manifest["windows"]
    technical = bool(
        windows
        and not manifest["technical_failure"]
        and all(item.get("technical_pass") is True for item in windows)
    )
    inheritance = True
    previous = manifest["source"]["restart_sha256"]
    for window in windows:
        inheritance = inheritance and window.get("input_restart_sha256") == previous
        if window.get("accepted_for_inheritance"):
            previous = window.get("output_restart_sha256")
    accepted_count = sum(item.get("accepted_for_inheritance") is True for item in windows)
    reached = bool(
        technical
        and inheritance
        and manifest["state"].get("destination_reached") is True
    )
    scientific_gate = (
        "NOT_EVALUATED_TECHNICAL_A1_BIDIRECTIONAL_2CV"
        if not technical else
        "PASS_RESTRAINED_2CV_CHAIN_REACHED_OPPOSITE_BASIN"
        if reached else
        "PARTIAL_2CV_ACTUAL_RESPONSE"
        if accepted_count else
        "FAIL_NO_ACTUAL_2CV_RESPONSE"
    )
    result = {
        "schema_version": 1,
        "status": TECHNICAL_PASS if technical else TECHNICAL_FAIL,
        "technical_complete": technical,
        "scientific_status": SCIENTIFIC_STATUS,
        "scientific_gate": scientific_gate,
        "task": manifest["task"],
        "terminal_reason": reason,
        "windows_attempted": len(windows),
        "accepted_windows": accepted_count,
        "final_accepted_fraction": manifest["state"]["accepted_fraction"],
        "final_step": manifest["state"]["step"],
        "inheritance_sha256_verified": inheritance,
        "destination_reached": reached,
        "final_geometry": manifest["accepted_geometry"],
        "final_coordinates": progress_coordinates(manifest["accepted_geometry"]),
        "carbonyl_and_oop_are_observables_only": True,
        "restrained_structures_are_not_ts": True,
        "automatic_downstream_action": "NONE",
    }
    write_json(root / "RESULT.json", result)
    write_json(root / ("PASS.json" if technical else "NOT_EVALUATED.json"), result)
    rows = []
    for path in sorted(root.rglob("*")):
        if (
            path.is_file()
            and path.name != "SHA256.tsv"
            and path.suffix in (".json", ".rst7")
        ):
            rows.append(f"{sha256(path)}  {path.relative_to(root)}")
    (root / "SHA256.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")
    return result


def merge_if_ready(output_root: pathlib.Path, array_job: str) -> bool:
    audit = output_root / "audit" / f"nylc_a1_bidirectional_2cv_{array_job}.json"
    audit.parent.mkdir(parents=True, exist_ok=True)
    with audit.with_suffix(".lock").open("a+", encoding="utf-8") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        paths = [
            output_root / f"attempt_{array_job}_{index}" / "RESULT.json"
            for index in range(ARRAY_TASKS)
        ]
        if not all(path.is_file() for path in paths):
            return False
        results = [read_json(path) for path in paths]
        technical_count = sum(
            result.get("technical_complete") is True for result in results
        )
        reached_count = sum(
            result.get("destination_reached") is True for result in results
        )
        payload = {
            "schema_version": 1,
            "status": (
                "NOT_EVALUATED_TECHNICAL_A1_BIDIRECTIONAL_2CV"
                if technical_count != ARRAY_TASKS else
                "PASS_ALL_EIGHT_RESTRAINED_2CV_CHAINS_REACHED_OPPOSITE_BASIN"
                if reached_count == ARRAY_TASKS else
                "PARTIAL_A1_BIDIRECTIONAL_2CV"
            ),
            "denominator_tasks": ARRAY_TASKS,
            "technical_pass_tasks": technical_count,
            "opposite_basin_reached_tasks": reached_count,
            "per_task": sorted(results, key=lambda item: item["task"]["task_index"]),
            "scientific_status": SCIENTIFIC_STATUS,
            "restrained_structures_are_not_ts": True,
            "automatic_downstream_action": "NONE",
        }
        if audit.exists():
            if read_json(audit) != payload:
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
