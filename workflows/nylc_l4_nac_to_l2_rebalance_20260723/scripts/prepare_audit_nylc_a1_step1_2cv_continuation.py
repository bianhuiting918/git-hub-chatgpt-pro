#!/usr/bin/env python3
"""Continue accepted product-to-reactant Step1 two-CV chains with bounded force tiers."""
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
BASE_PATH = HERE / "prepare_audit_nylc_a1_step1_bidirectional_2cv.py"
_SPEC = importlib.util.spec_from_file_location("_a1_2cv_base", BASE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"cannot import {BASE_PATH}")
BASE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(BASE)

ARRAY_TASKS = 8
ANCHOR_TASKS = (2, 3, 6, 7)
MAX_WINDOWS = 12
INITIAL_STEP = 1.0 / 32.0
MIN_STEP = 1.0 / 256.0
TECHNICAL_PASS = "PASS_TECHNICAL_A1_2CV_CONTINUATION"
TECHNICAL_FAIL = "NOT_EVALUATED_TECHNICAL_A1_2CV_CONTINUATION"

BASE.MAX_WINDOWS = MAX_WINDOWS
BASE.MIN_STEP = MIN_STEP

SOURCE_ROOT = (
    BASE.TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/a1_step1_bidirectional_2cv"
)
ANCHORS = {
    2: {
        "seed": "seed26723",
        "seed_index": 0,
        "mechanism": "ADDITION_FIRST",
        "source_fraction": 0.75,
        "restart_sha256": "cda206065902568f1406fa6f227f2d924b581a0140b12a2c0ba1bc639a037823",
    },
    3: {
        "seed": "seed26723",
        "seed_index": 0,
        "mechanism": "CONCERTED",
        "source_fraction": 0.75,
        "restart_sha256": "edb2a3a244cd2cefd90bfc3eb7ca9eb93a2b4bbda7f8fed3b15dc8e16a826b98",
    },
    6: {
        "seed": "seed26737",
        "seed_index": 1,
        "mechanism": "ADDITION_FIRST",
        "source_fraction": 0.625,
        "restart_sha256": "e4f2489b3e64117ff87363dcb107c44ca4e2b12e2d50b2b92e347e525037c564",
    },
    7: {
        "seed": "seed26737",
        "seed_index": 1,
        "mechanism": "CONCERTED",
        "source_fraction": 0.625,
        "restart_sha256": "a8037eb79eed6e3007f6ca14b51f5bcfe4d25459a306a77510562a22de838626",
    },
}
FORCE_TIERS = {
    "MODERATE": {"force_bond": 16.0, "force_pt": 12.0, "maxcyc": 1600, "ncyc": 400},
    "FIRM": {"force_bond": 24.0, "force_pt": 18.0, "maxcyc": 2000, "ncyc": 500},
}


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..7")
    source_task = ANCHOR_TASKS[index // 2]
    tier = ("MODERATE", "FIRM")[index % 2]
    anchor = ANCHORS[source_task]
    return {
        "task_index": index,
        "source_task_index": source_task,
        "source_array_job": "62306108",
        "seed": anchor["seed"],
        "seed_index": anchor["seed_index"],
        "source_basin": "P",
        "destination_basin": "R",
        "mechanism": anchor["mechanism"],
        "force_tier": tier,
        "source_fraction": anchor["source_fraction"],
        "source_restart_sha256": anchor["restart_sha256"],
    }


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": 8,
        "anchor_task_indices": list(ANCHOR_TASKS),
        "force_tiers": list(FORCE_TIERS),
        "maximum_windows_per_chain": MAX_WINDOWS,
        "initial_lambda_step": INITIAL_STEP,
        "minimum_lambda_step": MIN_STEP,
        "strict_serial_restart_inheritance": True,
        "failed_restart_inheritance": False,
        "carbonyl_and_oop": "observed_not_restrained",
        "frozen_step1_contract": BASE.describe()["frozen_step1_contract"],
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }


def initialize(task_index: int, root: pathlib.Path, commit: str) -> dict[str, Any]:
    if root.exists():
        raise FileExistsError(root)
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("runtime is not bound to a full Git commit")
    task = task_spec(task_index)
    source_root = SOURCE_ROOT / f"attempt_62306108_{task['source_task_index']}"
    source_result = BASE.read_json(source_root / "RESULT.json")
    source_manifest = BASE.read_json(source_root / "MANIFEST.json")
    restart = source_root / "accepted.rst7"
    expected_fraction = float(task["source_fraction"])
    checks = (
        source_result.get("technical_complete") is True,
        source_result.get("inheritance_sha256_verified") is True,
        int(source_result.get("accepted_windows", 0)) > 0,
        int(source_result["task"]["task_index"]) == int(task["source_task_index"]),
        math.isclose(
            float(source_result["final_accepted_fraction"]),
            expected_fraction,
            abs_tol=1.0e-12,
        ),
        restart.is_file(),
        BASE.sha256(restart) == task["source_restart_sha256"],
    )
    if not all(checks):
        raise ValueError("fixed accepted continuation source authority failed")
    contract, observed = BASE._contract(source_manifest)
    if BASE.sha256(BASE.PRMTOP) != BASE.AUTH.BASE.EXPECTED_PRMTOP_SHA256:
        raise ValueError("frozen prmtop SHA changed")
    root.mkdir(parents=True, exist_ok=False)
    accepted = root / "accepted.rst7"
    shutil.copy2(restart, accepted)
    geometry = BASE.AUTH._geometry(accepted, {"qm_contract": contract})
    tier = dict(FORCE_TIERS[task["force_tier"]])
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_2CV_CONTINUATION",
        "github_commit": commit,
        "scientific_status": BASE.SCIENTIFIC_STATUS,
        "task": task,
        "source": {
            "attempt": str(source_root),
            "result": str(source_root / "RESULT.json"),
            "restart": str(restart),
            "restart_sha256": task["source_restart_sha256"],
            "source_fraction": expected_fraction,
        },
        "prmtop": str(BASE.PRMTOP),
        "prmtop_sha256": BASE.sha256(BASE.PRMTOP),
        "qm_contract": contract,
        "observed_frozen_contract": observed,
        "accepted_restart": str(accepted),
        "accepted_restart_sha256": BASE.sha256(accepted),
        "accepted_geometry": geometry,
        "state": {
            "accepted_fraction": expected_fraction,
            "destination_fraction": 0.0,
            "direction": -1.0,
            "step": INITIAL_STEP,
            "destination_reached": False,
        },
        "continuation": {
            "force_tier": task["force_tier"],
            **tier,
            "maximum_windows": MAX_WINDOWS,
            "minimum_lambda_step": MIN_STEP,
        },
        "windows": [],
        "technical_failure": False,
        "terminal": False,
        "terminal_reason": None,
        "carbonyl_and_oop_are_observables_only": True,
        "restrained_structures_are_not_ts": True,
        "automatic_downstream_action": "NONE",
    }
    BASE.write_json(root / "MANIFEST.json", manifest)
    return manifest


def prepare_window(
    root: pathlib.Path,
    output: pathlib.Path,
    scratch: pathlib.Path,
    input_rst7: pathlib.Path,
    window_index: int,
) -> dict[str, Any]:
    prepared = BASE.prepare_window(root, output, scratch, input_rst7, window_index)
    manifest = BASE.read_json(root / "MANIFEST.json")
    tier = manifest["continuation"]
    stage = prepared["stage"]
    for key in ("force_bond", "force_pt", "maxcyc", "ncyc"):
        stage[key] = tier[key]
    BASE.write_json(output / "WINDOW_MANIFEST.json", prepared)
    (scratch / "stage.in").write_text(
        BASE.minimization_input(
            manifest["task"], stage, manifest["qm_contract"]["qmmask"]
        ),
        encoding="utf-8",
    )
    (scratch / "restraints.RST").write_text(
        BASE.restraints(stage), encoding="utf-8"
    )
    return prepared


def audit_window(
    root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path
) -> dict[str, Any]:
    return BASE.audit_window(root, output, scratch)


def _write_hashes(root: pathlib.Path) -> None:
    rows = []
    for path in sorted(root.rglob("*")):
        if (
            path.is_file()
            and path.name != "SHA256.tsv"
            and path.suffix in (".json", ".rst7")
        ):
            rows.append(f"{BASE.sha256(path)}  {path.relative_to(root)}")
    (root / "SHA256.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")


def finalize(root: pathlib.Path) -> dict[str, Any]:
    result = BASE.finalize(root)
    result["status"] = TECHNICAL_PASS if result["technical_complete"] else TECHNICAL_FAIL
    result["continuation_from_job"] = "62306108"
    result["continuation_force_tier"] = result["task"]["force_tier"]
    BASE.write_json(root / "RESULT.json", result)
    marker = "PASS.json" if result["technical_complete"] else "NOT_EVALUATED.json"
    BASE.write_json(root / marker, result)
    _write_hashes(root)
    return result


def merge_if_ready(output_root: pathlib.Path, array_job: str) -> bool:
    audit = output_root / "audit" / f"nylc_a1_2cv_continuation_{array_job}.json"
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
        reached = sum(item.get("destination_reached") is True for item in results)
        responsive = sum(int(item.get("accepted_windows", 0)) > 0 for item in results)
        payload = {
            "schema_version": 1,
            "status": (
                "NOT_EVALUATED_TECHNICAL_A1_2CV_CONTINUATION"
                if technical != ARRAY_TASKS
                else "PASS_ALL_CONTINUATIONS_REACHED_REACTANT_BASIN"
                if reached == ARRAY_TASKS
                else "PARTIAL_A1_2CV_CONTINUATION"
            ),
            "denominator_tasks": ARRAY_TASKS,
            "technical_pass_tasks": technical,
            "actual_response_tasks": responsive,
            "reactant_basin_reached_tasks": reached,
            "per_task": sorted(results, key=lambda item: item["task"]["task_index"]),
            "scientific_status": BASE.SCIENTIFIC_STATUS,
            "restrained_structures_are_not_ts": True,
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
    parser.add_argument(
        "--mode",
        required=True,
        choices=(
            "describe",
            "initialize",
            "prepare-window",
            "audit-window",
            "finalize",
            "merge-if-ready",
        ),
    )
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
            args.root.resolve(),
            args.output.resolve(),
            args.scratch.resolve(),
            args.input_rst7.resolve(),
            args.window_index,
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
