#!/usr/bin/env python3
"""Two-seed neutral-Thr267 Step1 first-window mechanism calibration."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import pathlib
import re
from typing import Any, Mapping

HERE = pathlib.Path(__file__).resolve().parent


def _load(name: str, path: pathlib.Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FWD = _load(
    "_neutral_forward_base",
    HERE / "prepare_audit_nylc_a1_step1_forward_force_calibration.py",
)
PB = _load(
    "_neutral_forward_integrity",
    HERE / "prepare_audit_nylc_a1_step1_protonation_boundary_authority.py",
)
BASE = FWD.BASE
TASK_ROOT = FWD.TASK_ROOT
PRMTOP = FWD.PRMTOP

ARRAY_TASKS = 4
MPI_RANKS = 8
MECHANISMS = ("ADDITION_FIRST_FORWARD", "FULLY_CONCERTED_FORWARD")
FORCE_SCALE = 1.0
OUTPUT_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/a1_step1_neutral_forward_calibration"
)
AUTHORITY_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/a1_step1_protonation_boundary_authority"
)
SOURCES = {
    "seed26723": {
        "task_index": 0,
        "attempt": AUTHORITY_ROOT / "attempt_62533819_0",
        "restart": AUTHORITY_ROOT / "attempt_62533819_0/source.rst7",
        "restart_sha256": "486e5c3d479e42d78c2d87a670320bbce9d1fcedaa0d5ad3bd99bf4c696506b2",
        "baseline_result": AUTHORITY_ROOT
        / "attempt_62533819_0/block_01/RESULT.json",
        "baseline_result_sha256": "8769352cc5af3513ef3ac89d257e9f5d5653d35f079a25615bcea4e9d5e75161",
    },
    "seed26737": {
        "task_index": 1,
        "attempt": AUTHORITY_ROOT / "attempt_62533819_1",
        "restart": AUTHORITY_ROOT / "attempt_62533819_1/source.rst7",
        "restart_sha256": "bd1717782acb8b6ec2d2c92a1b4f84d4354baf5797df4577fa80aa57541edccc",
        "baseline_result": AUTHORITY_ROOT
        / "attempt_62533819_1/block_01/RESULT.json",
        "baseline_result_sha256": "044e54aa8a2e83a653e544fea13d42589047d43d6464e16bc4fb2c752761cbb8",
    },
}
QM_CONTRACT = {"qm_atoms": 146, "qm_charge": 0, "electrons": 510, "link_atoms": 6}
SCIENTIFIC_STATUS = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..3")
    seed = ("seed26723", "seed26737")[index // 2]
    return {
        "task_index": index,
        "seed_index": SOURCES[seed]["task_index"],
        "seed": seed,
        "mechanism": MECHANISMS[index % 2],
        "force_scale": FORCE_SCALE,
        "starting_state": "NALPHA_H2_OG1H",
        "qm_contract_key": "current",
        "source_basin": "NEUTRAL_RAW_NAC",
        "destination": "FIRST_FORWARD_WINDOW",
    }


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "mapping": "2 neutral seeds x 2 mechanisms x one 1x first window",
        "mechanisms": list(MECHANISMS),
        "force_scale": FORCE_SCALE,
        "first_window_only": True,
        "qm_contract": dict(QM_CONTRACT),
        "source_authority_job": 62533819,
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }


def _qpt(geometry: Mapping[str, Any]) -> float:
    return float(geometry["nalpha_hg1_A"]) - float(geometry["hg1_n3_A"])


def _coord(geometry: Mapping[str, Any], name: str) -> float:
    if name == "attack":
        return float(geometry["attack_A"])
    if name == "cn":
        return float(geometry["c12_n3_A"])
    if name == "carbonyl":
        return float(geometry["c12_o2_A"])
    if name == "pt":
        return _qpt(geometry)
    raise KeyError(name)


def initialize(task_index: int, root: pathlib.Path, commit: str) -> dict[str, Any]:
    if root.exists():
        raise FileExistsError(root)
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("runtime is not bound to a full Git commit")
    task = task_spec(task_index)
    source = SOURCES[task["seed"]]
    restart = pathlib.Path(source["restart"])
    baseline_path = pathlib.Path(source["baseline_result"])
    if not restart.is_file() or sha256(restart) != source["restart_sha256"]:
        raise ValueError("fixed neutral source restart SHA failed")
    if not baseline_path.is_file() or sha256(baseline_path) != source["baseline_result_sha256"]:
        raise ValueError("matched unrestrained baseline SHA failed")

    authority = BASE.read_json(pathlib.Path(source["attempt"]) / "SOURCE_MANIFEST.json")
    authority_result = BASE.read_json(pathlib.Path(source["attempt"]) / "RESULT.json")
    baseline = BASE.read_json(baseline_path)
    if not all((
        authority_result.get("status") == "PASS_A1_PROTONATION_BOUNDARY_AUTHORITY",
        authority_result.get("accepted_blocks") == 4,
        authority["task"].get("starting_state") == "NALPHA_H2_OG1H",
        authority["task"].get("qm_contract_key") == "current",
        authority.get("source_integrity", {}).get("pass") is True,
        baseline.get("accepted") is True,
        baseline.get("input_restart_sha256") == source["restart_sha256"],
    )):
        raise ValueError("neutral authority source contract failed")
    contract = authority["qm_contract"]
    if int(contract["qm_atom_count"]) != 146 or int(contract["qmcharge"]) != 0:
        raise ValueError("neutral current-QM contract changed")
    if sha256(PRMTOP) != PB.PRMTOP_SHA256:
        raise ValueError("frozen prmtop SHA changed")

    root.mkdir(parents=True, exist_ok=False)
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_NEUTRAL_FORWARD_CALIBRATION",
        "github_commit": commit,
        "scientific_status": SCIENTIFIC_STATUS,
        "task": task,
        "source": {
            "restart": str(restart),
            "restart_sha256": source["restart_sha256"],
            "authority_attempt": str(source["attempt"]),
        },
        "matched_unrestrained_baseline": {
            "result": str(baseline_path),
            "result_sha256": source["baseline_result_sha256"],
            "geometry": baseline["geometry"],
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": sha256(PRMTOP),
        "qm_contract": contract,
        "source_geometry": authority["source_geometry"],
        "first_window_only": True,
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }
    BASE.write_json(root / "MANIFEST.json", manifest)
    return manifest


def prepare(root: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    return FWD.prepare(root, scratch)


def audit(root: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    result = FWD.audit(root, scratch)
    manifest = BASE.read_json(root / "MANIFEST.json")
    integrity = {"pass": False}
    if result.get("technical_complete") is True and (root / "stage.rst7").is_file():
        integrity = PB.measure_integrity(PRMTOP, root / "stage.rst7")

    active = list(result.get("response_gate", {}).get("active_coordinates", []))
    baseline_geometry = manifest["matched_unrestrained_baseline"]["geometry"]
    final_geometry = result.get("final_geometry")
    target_geometry = result.get("target_geometry")
    checks: dict[str, Any] = {}
    if final_geometry is not None and target_geometry is not None:
        for name in active:
            baseline_error = abs(_coord(baseline_geometry, name) - _coord(target_geometry, name))
            forced_error = abs(_coord(final_geometry, name) - _coord(target_geometry, name))
            checks[name] = {
                "baseline_abs_error_A": baseline_error,
                "forced_abs_error_A": forced_error,
                "forced_closer_to_target": forced_error + 1.0e-4 < baseline_error,
            }
    baseline_gate = {"checks": checks, "all": bool(active) and all(
        item["forced_closer_to_target"] for item in checks.values()
    )}
    eligible = bool(
        result.get("technical_complete") is True
        and integrity.get("pass") is True
        and result.get("response_gate", {}).get("all") is True
        and result.get("chemical_guard", {}).get("all") is True
        and baseline_gate["all"] is True
    )
    result["thr267_chemical_integrity"] = integrity
    result["response_vs_matched_unrestrained_baseline"] = baseline_gate
    result["eligible_for_forward_continuation"] = eligible
    result["scientific_gate"] = (
        "PASS_NEUTRAL_FIRST_WINDOW_RESPONSE"
        if eligible else "FAIL_NEUTRAL_FIRST_WINDOW_RESPONSE_OR_INTEGRITY"
    )
    if result.get("technical_complete") is True and integrity.get("pass") is not True:
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
    parser.add_argument("--mode", required=True, choices=("describe", "initialize", "prepare", "audit"))
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
