#!/usr/bin/env python3
"""Two-seed neutral-Thr attack-only first window with the correct attack-angle guard."""
from __future__ import annotations

import argparse
import importlib.util
import json
import math
import pathlib
from typing import Any, Mapping

HERE = pathlib.Path(__file__).resolve().parent


def _load(name: str, path: pathlib.Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CAUSAL = _load(
    "_angle_guard_causal",
    HERE / "prepare_audit_nylc_a1_step1_attack_carbonyl_causal.py",
)
NEUTRAL = CAUSAL.NEUTRAL
FWD = CAUSAL.FWD
PB = CAUSAL.PB
BASE = CAUSAL.BASE
PRMTOP = CAUSAL.PRMTOP
TASK_ROOT = CAUSAL.TASK_ROOT

ARRAY_TASKS = 2
MPI_RANKS = 8
ATTACK_DELTA_A = -0.04
FORCE_ATTACK = 24.0
FORCE_ANGLE_GUARD = 2.0
ANGLE_FLAT_BOTTOM_DEG = (100.0, 115.0)
OUTPUT_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/a1_step1_angle_protected_attack"
)
SOURCES = {
    seed: {
        "task_index": int(source["task_index"]),
        "attempt": pathlib.Path(source["attempt"]),
        "restart": pathlib.Path(source["restart"]),
        "restart_sha256": str(source["restart_sha256"]),
        "baseline_result": pathlib.Path(source["baseline_result"]),
        "baseline_result_sha256": str(source["baseline_result_sha256"]),
    }
    for seed, source in NEUTRAL.SOURCES.items()
}
SCIENTIFIC_STATUS = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..1")
    seed = ("seed26723", "seed26737")[index]
    return {
        "task_index": index,
        "seed_index": SOURCES[seed]["task_index"],
        "seed": seed,
        "mechanism": "ADDITION_FIRST_FORWARD",
        "route": "ANGLE_PROTECTED_ATTACK_ONLY",
        "source_kind": "ORIGINAL_NEUTRAL_SOURCE",
        "starting_state": "NALPHA_H2_OG1H",
        "qm_contract_key": "current",
        "source_basin": "NEUTRAL_RAW_NAC",
        "destination": "ANGLE_PROTECTED_FIRST_FORWARD_WINDOW",
        "force_scale": 1.0,
    }


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "mapping": "2 original healthy neutral-Thr seeds x one attack-only first window",
        "first_window_only": True,
        "attack_delta_A": ATTACK_DELTA_A,
        "force_attack_kcal_mol_A2": FORCE_ATTACK,
        "force_angle_guard": FORCE_ANGLE_GUARD,
        "angle_definition": "O2-C12-OG1",
        "angle_flat_bottom_deg": list(ANGLE_FLAT_BOTTOM_DEG),
        "unrestrained_coordinates": ["C12-O2", "C12-N3", "qPT"],
        "stop_conditions": [
            "technical_or_numerical_failure",
            "thr267_chemical_integrity_failure",
            "attack_response_wrong_direction",
            "attack_angle_below_100_deg",
        ],
        "source_authority_job": 62533819,
        "qm_contract": dict(NEUTRAL.QM_CONTRACT),
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }


def stage_spec(task: Mapping[str, Any], source: Mapping[str, Any]) -> dict[str, Any]:
    target = {
        "attack_A": float(source["attack_A"]) + ATTACK_DELTA_A,
        "c12_n3_A": float(source["c12_n3_A"]),
        "c12_o2_A": float(source["c12_o2_A"]),
        "nalpha_hg1_A": float(source["nalpha_hg1_A"]),
        "hg1_n3_A": float(source["hg1_n3_A"]),
    }
    return {
        "window_index": 0,
        "mechanism": str(task["mechanism"]),
        "route": str(task["route"]),
        "force_scale": 1.0,
        "targets": target,
        "active_coordinates": ["attack"],
        "force_bond": FORCE_ATTACK,
        "force_angle_guard": FORCE_ANGLE_GUARD,
        "angle_definition": "O2-C12-OG1",
        "angle_flat_bottom_deg": list(ANGLE_FLAT_BOTTOM_DEG),
        "force_carbonyl": 0.0,
        "force_pt": 0.0,
        "maxcyc": 2200,
        "ncyc": 550,
    }


def _angle_restraint(
    atom1: int,
    atom2: int,
    atom3: int,
    lower: float,
    upper: float,
    force: float,
) -> str:
    return (
        f"&rst iat={atom1},{atom2},{atom3}, "
        f"r1={max(0.0, lower - 25.0):.3f}, "
        f"r2={lower:.3f}, r3={upper:.3f}, "
        f"r4={min(180.0, upper + 25.0):.3f}, "
        f"rk2={force:.1f}, rk3={force:.1f}, /\n"
    )


def restraints(stage: Mapping[str, Any]) -> str:
    target = stage["targets"]
    lower, upper = stage["angle_flat_bottom_deg"]
    return "".join(
        (
            FWD.CAL.BRIDGE._tight_distance_restraint(
                FWD.REACTIVE["og1"],
                FWD.REACTIVE["c12"],
                target["attack_A"],
                stage["force_bond"],
            ),
            _angle_restraint(
                FWD.REACTIVE["o2"],
                FWD.REACTIVE["c12"],
                FWD.REACTIVE["og1"],
                lower,
                upper,
                stage["force_angle_guard"],
            ),
        )
    )


def initialize(task_index: int, root: pathlib.Path, commit: str) -> dict[str, Any]:
    index = int(task_index)
    manifest = NEUTRAL.initialize(index * 2, root, commit)
    task = task_spec(index)
    if manifest["task"]["seed"] != task["seed"]:
        raise ValueError("neutral authority task mapping changed")
    manifest["schema_version"] = 1
    manifest["status"] = "READY_A1_ANGLE_PROTECTED_ATTACK"
    manifest["scientific_status"] = SCIENTIFIC_STATUS
    manifest["task"] = task
    manifest["angle_definition"] = "O2-C12-OG1"
    manifest["angle_flat_bottom_deg"] = list(ANGLE_FLAT_BOTTOM_DEG)
    manifest["first_window_only"] = True
    manifest["automatic_downstream_action"] = "NONE"
    BASE.write_json(root / "MANIFEST.json", manifest)
    return manifest


def prepare(root: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest = BASE.read_json(root / "MANIFEST.json")
    if scratch.exists():
        raise FileExistsError(scratch)
    scratch.mkdir(parents=True, exist_ok=False)
    stage = stage_spec(manifest["task"], manifest["source_geometry"])
    (scratch / "stage.in").write_text(
        FWD.CAL.BRIDGE.tetra_minimization_input(
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


def audit(root: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    result = FWD.audit(root, scratch)
    integrity = {"pass": False}
    if result.get("technical_complete") is True and (root / "stage.rst7").is_file():
        integrity = PB.measure_integrity(PRMTOP, root / "stage.rst7")

    final_geometry = result.get("final_geometry") or {}
    final_angle = final_geometry.get("attack_angle_deg")
    attack_gate = result.get("response_gate", {}).get("attack", {})
    attack_pass = attack_gate.get("pass") is True and attack_gate.get("same_direction") is True
    angle_pass = bool(
        final_angle is not None
        and math.isfinite(float(final_angle))
        and float(final_angle) >= ANGLE_FLAT_BOTTOM_DEG[0]
    )
    technical = result.get("technical_complete") is True
    healthy = bool(technical and integrity.get("pass") is True)
    response_pass = bool(healthy and attack_pass and angle_pass)

    result["thr267_chemical_integrity"] = integrity
    result["attack_angle_definition"] = "O2-C12-OG1"
    result["angle_flat_bottom_deg"] = list(ANGLE_FLAT_BOTTOM_DEG)
    result["angle_guard_pass"] = angle_pass
    result["angle_protected_response_pass"] = response_pass
    result["eligible_for_forward_continuation"] = False
    result["automatic_downstream_action"] = "NONE"
    result["scientific_gate"] = (
        "NOT_EVALUATED_TECHNICAL_A1_ANGLE_PROTECTED_ATTACK"
        if not healthy
        else "PASS_ANGLE_PROTECTED_ATTACK_RESPONSE"
        if response_pass
        else "FAIL_ANGLE_PROTECTED_ATTACK_RESPONSE"
    )
    if technical and integrity.get("pass") is not True:
        result["status"] = "NOT_EVALUATED_TECHNICAL_CHEMICAL_INTEGRITY"
        result["technical_complete"] = False
        pass_path = root / "PASS.json"
        if pass_path.exists():
            pass_path.unlink()
        BASE.write_json(root / "NOT_EVALUATED.json", result)
    else:
        BASE.write_json(root / "PASS.json", result)
    BASE.write_json(root / "RESULT.json", result)
    FWD._write_hashes(root)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode", required=True, choices=("describe", "initialize", "prepare", "audit")
    )
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--root", type=pathlib.Path)
    parser.add_argument("--scratch", type=pathlib.Path)
    parser.add_argument("--github-commit", default="unknown")
    args = parser.parse_args()
    if args.mode == "describe":
        print(json.dumps(describe(), sort_keys=True))
    elif args.mode == "initialize":
        initialize(args.task_index, args.root.resolve(), args.github_commit)
    elif args.mode == "prepare":
        prepare(args.root.resolve(), args.scratch.resolve())
    else:
        audit(args.root.resolve(), args.scratch.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
