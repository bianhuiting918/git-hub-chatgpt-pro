#!/usr/bin/env python3
"""Two-seed one-window test of whether carbonyl pulling drives attack-angle collapse."""
from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import re
import shutil
from typing import Any, Mapping

HERE = pathlib.Path(__file__).resolve().parent


def _load(name: str, path: pathlib.Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


NEUTRAL = _load(
    "_neutral_forward_causal_base",
    HERE / "prepare_audit_nylc_a1_step1_neutral_forward_calibration.py",
)
FWD = NEUTRAL.FWD
PB = NEUTRAL.PB
BASE = NEUTRAL.BASE
TASK_ROOT = NEUTRAL.TASK_ROOT
PRMTOP = NEUTRAL.PRMTOP

ARRAY_TASKS = 6
MPI_RANKS = 8
ROUTES = (
    "ENDPOINT_ATTACK_CARBONYL",
    "ENDPOINT_ATTACK_ONLY",
    "SOURCE_ATTACK_ONLY",
)
OUTPUT_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/a1_step1_attack_carbonyl_causal"
)
PREVIOUS_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/a1_step1_neutral_forward_calibration"
)
SOURCES = {
    "seed26723": {
        "authority": NEUTRAL.SOURCES["seed26723"],
        "endpoint": PREVIOUS_ROOT / "attempt_62539741_0/stage.rst7",
        "endpoint_sha256": "58110d8e77c3b6beb15154d74140a847b31aa545187072f9737288e7ae0635f8",
        "endpoint_result": PREVIOUS_ROOT / "attempt_62539741_0/RESULT.json",
        "endpoint_result_sha256": "d5efa25c391cb362702c797735e9f77ac2cfc31ca8cab2a288975bb39f8a93da",
    },
    "seed26737": {
        "authority": NEUTRAL.SOURCES["seed26737"],
        "endpoint": PREVIOUS_ROOT / "attempt_62539741_2/stage.rst7",
        "endpoint_sha256": "18472aa8a2770a55431bed09035ea633e10e2cd0950d33979eb812cf3c1cd0ad",
        "endpoint_result": PREVIOUS_ROOT / "attempt_62539741_2/RESULT.json",
        "endpoint_result_sha256": "b2a368c251a6ea625644c089271093c05694b6787ebd2f59c9875ab6d873848c",
    },
}
ATTACK_DELTA_A = -0.04
CARBONYL_DELTA_A = 0.03
FORCE_ATTACK = 24.0
FORCE_CARBONYL = 30.0
SCIENTIFIC_STATUS = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..5")
    seed = ("seed26723", "seed26737")[index // 3]
    route = ROUTES[index % 3]
    return {
        "task_index": index,
        "seed_index": index // 3,
        "seed": seed,
        "route": route,
        "mechanism": route,
        "source_kind": (
            "ORIGINAL_NEUTRAL_SOURCE"
            if route == "SOURCE_ATTACK_ONLY"
            else "OFF_ANGLE_ENDPOINT"
        ),
        "starting_state": "NALPHA_H2_OG1H",
        "qm_contract_key": "current",
        "source_basin": "NEUTRAL_RAW_NAC",
        "destination": "ONE_WINDOW_CAUSAL_DIAGNOSTIC",
        "force_scale": 1.0,
    }


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "mapping": "2 neutral seeds x 3 one-window causal routes",
        "routes": list(ROUTES),
        "force_attack_kcal_mol_A2": FORCE_ATTACK,
        "force_carbonyl_kcal_mol_A2": FORCE_CARBONYL,
        "attack_delta_A": ATTACK_DELTA_A,
        "carbonyl_delta_A": CARBONYL_DELTA_A,
        "angle_cn_qpt_restraints": False,
        "first_window_only": True,
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }


def stage_spec(task: Mapping[str, Any], source: Mapping[str, Any]) -> dict[str, Any]:
    route = str(task["route"])
    if route not in ROUTES:
        raise ValueError(route)
    target = {
        "attack_A": float(source["attack_A"]) + ATTACK_DELTA_A,
        "c12_n3_A": float(source["c12_n3_A"]),
        "c12_o2_A": float(source["c12_o2_A"]),
        "nalpha_hg1_A": float(source["nalpha_hg1_A"]),
        "hg1_n3_A": float(source["hg1_n3_A"]),
    }
    active = ["attack"]
    if route == "ENDPOINT_ATTACK_CARBONYL":
        target["c12_o2_A"] += CARBONYL_DELTA_A
        active.append("carbonyl")
    return {
        "window_index": 0,
        "mechanism": route,
        "force_scale": 1.0,
        "targets": target,
        "active_coordinates": active,
        "force_bond": FORCE_ATTACK,
        "force_carbonyl": FORCE_CARBONYL,
        "force_pt": FWD.CAL.BASE_FORCES["force_pt"],
        "maxcyc": 2200,
        "ncyc": 550,
    }


def restraints(stage: Mapping[str, Any]) -> str:
    target = stage["targets"]
    rows = [
        FWD.CAL.BRIDGE._tight_distance_restraint(
            FWD.REACTIVE["og1"],
            FWD.REACTIVE["c12"],
            target["attack_A"],
            stage["force_bond"],
        )
    ]
    if "carbonyl" in stage["active_coordinates"]:
        rows.append(
            FWD.CAL.BRIDGE._tight_distance_restraint(
                FWD.REACTIVE["c12"],
                FWD.REACTIVE["o2"],
                target["c12_o2_A"],
                stage["force_carbonyl"],
            )
        )
    return "".join(rows)


def _validated_contract(seed: str) -> dict[str, Any]:
    authority = SOURCES[seed]["authority"]
    attempt = pathlib.Path(authority["attempt"])
    manifest = BASE.read_json(attempt / "SOURCE_MANIFEST.json")
    result = BASE.read_json(attempt / "RESULT.json")
    if not all(
        (
            result.get("status") == "PASS_A1_PROTONATION_BOUNDARY_AUTHORITY",
            result.get("accepted_blocks") == 4,
            manifest["task"].get("starting_state") == "NALPHA_H2_OG1H",
            manifest["task"].get("qm_contract_key") == "current",
            manifest.get("source_integrity", {}).get("pass") is True,
        )
    ):
        raise ValueError("neutral authority source contract failed")
    contract = manifest["qm_contract"]
    if int(contract["qm_atom_count"]) != 146 or int(contract["qmcharge"]) != 0:
        raise ValueError("neutral current-QM contract changed")
    return contract


def initialize(task_index: int, root: pathlib.Path, commit: str) -> dict[str, Any]:
    if root.exists():
        raise FileExistsError(root)
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("runtime is not bound to a full Git commit")
    task = task_spec(task_index)
    source = SOURCES[task["seed"]]
    contract = _validated_contract(task["seed"])
    if task["source_kind"] == "OFF_ANGLE_ENDPOINT":
        restart = pathlib.Path(source["endpoint"])
        expected_sha = source["endpoint_sha256"]
        result_path = pathlib.Path(source["endpoint_result"])
        if (
            not result_path.is_file()
            or NEUTRAL.sha256(result_path) != source["endpoint_result_sha256"]
        ):
            raise ValueError("fixed endpoint result SHA failed")
        prior = BASE.read_json(result_path)
        if not all(
            (
                prior.get("technical_complete") is True,
                prior.get("thr267_chemical_integrity", {}).get("pass") is True,
                prior.get("task", {}).get("seed") == task["seed"],
            )
        ):
            raise ValueError("endpoint source authority failed")
    else:
        restart = pathlib.Path(source["authority"]["restart"])
        expected_sha = source["authority"]["restart_sha256"]
        result_path = pathlib.Path(source["authority"]["attempt"]) / "RESULT.json"
    if not restart.is_file() or NEUTRAL.sha256(restart) != expected_sha:
        raise ValueError("fixed source restart SHA failed")
    if NEUTRAL.sha256(PRMTOP) != PB.PRMTOP_SHA256:
        raise ValueError("frozen prmtop SHA changed")
    integrity = PB.measure_integrity(PRMTOP, restart)
    if integrity.get("pass") is not True:
        raise ValueError("source Thr267 chemical integrity failed")
    geometry = BASE.AUTH._geometry(restart, {"qm_contract": contract})

    root.mkdir(parents=True, exist_ok=False)
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_ATTACK_CARBONYL_CAUSAL",
        "github_commit": commit,
        "scientific_status": SCIENTIFIC_STATUS,
        "task": task,
        "source": {
            "restart": str(restart),
            "restart_sha256": expected_sha,
            "authority_result": str(result_path),
            "source_kind": task["source_kind"],
        },
        "source_geometry": geometry,
        "source_thr267_chemical_integrity": integrity,
        "prmtop": str(PRMTOP),
        "prmtop_sha256": NEUTRAL.sha256(PRMTOP),
        "qm_contract": contract,
        "first_window_only": True,
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }
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
    source_angle = result.get("source_geometry", {}).get("attack_angle_deg")
    final_angle = (result.get("final_geometry") or {}).get("attack_angle_deg")
    response_pass = result.get("response_gate", {}).get("all") is True
    technical = result.get("technical_complete") is True
    healthy = bool(technical and integrity.get("pass") is True)
    result["thr267_chemical_integrity"] = integrity
    result["attack_angle_change_deg"] = (
        float(final_angle) - float(source_angle)
        if source_angle is not None and final_angle is not None
        else None
    )
    result["causal_response_pass"] = bool(healthy and response_pass)
    result["eligible_for_forward_continuation"] = False
    result["scientific_gate"] = (
        "NOT_EVALUATED_TECHNICAL_A1_ATTACK_CARBONYL_CAUSAL"
        if not healthy
        else "PASS_ONE_WINDOW_CAUSAL_RESPONSE"
        if response_pass
        else "FAIL_ONE_WINDOW_CAUSAL_RESPONSE"
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
