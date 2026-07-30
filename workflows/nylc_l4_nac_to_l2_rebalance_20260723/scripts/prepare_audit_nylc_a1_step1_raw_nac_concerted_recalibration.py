#!/usr/bin/env python3
"""Local concerted first-window recalibration after the raw-NAC v2 matrix."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import pathlib
import re
import shutil
from typing import Any, Mapping

HERE = pathlib.Path(__file__).resolve().parent
RAW_PATH = HERE / "prepare_audit_nylc_a1_step1_raw_nac_forward_calibration.py"
_SPEC = importlib.util.spec_from_file_location("_raw_nac_forward", RAW_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"cannot import {RAW_PATH}")
RAW = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(RAW)

TASK_ROOT = RAW.TASK_ROOT
PRMTOP = RAW.PRMTOP
PRMTOP_SHA256 = RAW.PRMTOP_SHA256
CALIBRATION_JOB = "62471131"
CALIBRATION_ROOT = RAW.OUTPUT_ROOT
CALIBRATION_AUDIT_V2 = (
    CALIBRATION_ROOT / "audit"
    / f"nylc_a1_raw_nac_forward_calibration_{CALIBRATION_JOB}_v2.json"
)
OUTPUT_ROOT = (
    TASK_ROOT / "a1_activated_nac_20260726/qmmm"
    / "a1_step1_raw_nac_concerted_recalibration"
)
ARRAY_TASKS = 4
MPI_RANKS = 8
VARIANTS = (
    ("HALF_WINDOW_SCALE8", 8, 0.5),
    ("FULL_WINDOW_SCALE16", 16, 1.0),
)
AUTO_CONTINUATION = False
AUTO_RELEASE = False
AUTO_SHOOTING = False
AUTO_PMF = False
SCIENTIFIC_STATUS = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"


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


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..3")
    seed = (26723, 26737)[index // 2]
    variant, scale, fraction = VARIANTS[index % 2]
    baseline_index = 0 if seed == 26723 else 9
    return {
        "task_index": index,
        "seed": seed,
        "seed_index": 0 if seed == 26723 else 1,
        "mode": "FULLY_CONCERTED_RAW_NAC",
        "mechanism": "FULLY_CONCERTED_RAW_NAC",
        "variant": variant,
        "force_scale": scale,
        "scale": scale,
        "window_fraction": fraction,
        "baseline_source_task_index": baseline_index,
        "source_calibration_job": CALIBRATION_JOB,
        "source_basin": "RAW_UNBIASED_NAC",
        "first_window_only": True,
    }


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "mapping": "2 seeds x (half-window scale8 plus full-window scale16)",
        "source_calibration_job": CALIBRATION_JOB,
        "source_calibration_audit_v2": str(CALIBRATION_AUDIT_V2),
        "variants": [
            {"variant": name, "force_scale": scale, "window_fraction": fraction}
            for name, scale, fraction in VARIANTS
        ],
        "matched_baseline_reused": True,
        "addition_first_recomputed": False,
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
        "qm_contract": dict(RAW.QM_CONTRACT),
    }


def stage_spec(spec: Mapping[str, Any], source: Mapping[str, Any]) -> dict[str, Any]:
    fraction = float(spec["window_fraction"])
    scale = int(spec["force_scale"])
    target = {
        "attack_A": float(source["attack_A"]) + RAW.FIRST_WINDOW_DELTAS_A["attack"] * fraction,
        "c12_n3_A": float(source["c12_n3_A"]) + RAW.FIRST_WINDOW_DELTAS_A["cn"] * fraction,
        "c12_o2_A": float(source["c12_o2_A"]) + RAW.FIRST_WINDOW_DELTAS_A["carbonyl"] * fraction,
        "nalpha_hg1_A": float(source["nalpha_hg1_A"]) + RAW.FIRST_WINDOW_DELTAS_A["nalpha_hg1"] * fraction,
        "hg1_n3_A": float(source["hg1_n3_A"]) + RAW.FIRST_WINDOW_DELTAS_A["hg1_n3"] * fraction,
    }
    return {
        "window_index": 0,
        "mode": spec["mode"],
        "mechanism": spec["mechanism"],
        "variant": spec["variant"],
        "force_scale": scale,
        "window_fraction": fraction,
        "targets": target,
        "active_coordinates": ["attack", "cn", "carbonyl", "pt"],
        "force_bond": RAW.FWD.CAL.BASE_FORCES["force_bond"] * scale,
        "force_carbonyl": RAW.FWD.CAL.BASE_FORCES["force_carbonyl"] * scale,
        "force_pt": RAW.FWD.CAL.BASE_FORCES["force_pt"] * scale,
        "maxcyc": 2200,
        "ncyc": 550,
    }


def restraints(stage: Mapping[str, Any]) -> str:
    target = stage["targets"]
    rows = [
        RAW.FWD.CAL.BRIDGE._tight_distance_restraint(
            RAW.REACTIVE["og1"], RAW.REACTIVE["c12"],
            target["attack_A"], stage["force_bond"],
        ),
        RAW.FWD.CAL.BRIDGE._tight_distance_restraint(
            RAW.REACTIVE["c12"], RAW.REACTIVE["o2"],
            target["c12_o2_A"], stage["force_carbonyl"],
        ),
        RAW.FWD.CAL.BRIDGE._tight_distance_restraint(
            RAW.REACTIVE["c12"], RAW.REACTIVE["n3"],
            target["c12_n3_A"], stage["force_bond"],
        ),
        RAW.FWD.CAL.BRIDGE._tight_distance_restraint(
            RAW.REACTIVE["nalpha"], RAW.REACTIVE["hg1"],
            target["nalpha_hg1_A"], stage["force_pt"],
        ),
        RAW.FWD.CAL.BRIDGE._tight_distance_restraint(
            RAW.REACTIVE["hg1"], RAW.REACTIVE["n3"],
            target["hg1_n3_A"], stage["force_pt"],
        ),
    ]
    return "".join(rows)


def initialize(task_index: int, root: pathlib.Path, commit: str) -> dict[str, Any]:
    if root.exists():
        raise FileExistsError(root)
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("runtime is not bound to a full Git commit")
    if not CALIBRATION_AUDIT_V2.is_file():
        raise FileNotFoundError(CALIBRATION_AUDIT_V2)
    spec = task_spec(task_index)
    audit = read_json(CALIBRATION_AUDIT_V2)
    selections = [
        row for row in audit.get("selected_minimum_force_scales", [])
        if int(row.get("seed", -1)) == spec["seed"]
        and row.get("mode") == "FULLY_CONCERTED_RAW_NAC"
    ]
    if len(selections) != 1 or selections[0].get("selected_weakest_scale") is not None:
        raise ValueError("concerted branch is not uniquely missing in v2 calibration")

    baseline_root = CALIBRATION_ROOT / (
        f"attempt_{CALIBRATION_JOB}_{spec['baseline_source_task_index']}"
    )
    baseline_manifest = read_json(baseline_root / "MANIFEST.json")
    baseline_result = read_json(baseline_root / "RESULT.json")
    source_restart = baseline_root / "source.rst7"
    checks = (
        baseline_result.get("technical_complete") is True,
        baseline_result.get("matched_baseline") is True,
        int(baseline_result["task"]["seed"]) == spec["seed"],
        source_restart.is_file(),
        sha256(source_restart) == RAW.SOURCES[spec["seed"]]["start_rst7_sha256"],
        sha256(PRMTOP) == PRMTOP_SHA256,
    )
    if not all(checks):
        raise ValueError("matched raw-NAC source/baseline authority failed")

    root.mkdir(parents=True, exist_ok=False)
    persistent_source = root / "source.rst7"
    shutil.copy2(source_restart, persistent_source)
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_RAW_NAC_CONCERTED_RECALIBRATION",
        "github_commit": commit,
        "scientific_status": SCIENTIFIC_STATUS,
        "task": spec,
        "source": {
            "restart": str(persistent_source),
            "restart_sha256": sha256(persistent_source),
            "original_restart": str(source_restart),
            "matched_baseline_result": str(baseline_root / "RESULT.json"),
            "matched_baseline_result_sha256": sha256(baseline_root / "RESULT.json"),
            "calibration_audit_v2": str(CALIBRATION_AUDIT_V2),
            "calibration_audit_v2_sha256": sha256(CALIBRATION_AUDIT_V2),
        },
        "source_geometry": baseline_manifest["source_geometry"],
        "qm_contract": baseline_manifest["qm_contract"],
        "prmtop": str(PRMTOP),
        "prmtop_sha256": sha256(PRMTOP),
        "first_window_only": True,
        "matched_baseline_reused": True,
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }
    write_json(root / "MANIFEST.json", manifest)
    return manifest


def prepare(root: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest = read_json(root / "MANIFEST.json")
    if scratch.exists():
        raise FileExistsError(scratch)
    scratch.mkdir(parents=True, exist_ok=False)
    stage = stage_spec(manifest["task"], manifest["source_geometry"])
    stage_input = RAW.FWD.CAL.BRIDGE.tetra_minimization_input(
        manifest["task"], stage, manifest["qm_contract"]["qmmask"]
    )
    if "drms=0.01" not in stage_input:
        raise ValueError("tight minimization authority changed")
    restraint_text = restraints(stage)
    if restraint_text.count("&rst") != 5 or r"\n" in restraint_text:
        raise ValueError("concerted restraint authority is not five physical records")
    (scratch / "stage.in").write_text(stage_input, encoding="utf-8")
    (scratch / "restraints.RST").write_text(restraint_text, encoding="utf-8")
    prepared = {
        "schema_version": 1,
        "stage": stage,
        "input_restart": manifest["source"]["restart"],
        "input_restart_sha256": manifest["source"]["restart_sha256"],
        "reaction_coordinate_restraints": 5,
        "nmropt": 1,
        "first_window_only": True,
        "shooting_forbidden": True,
    }
    write_json(root / "WINDOW_MANIFEST.json", prepared)
    return prepared


def _write_hashes(root: pathlib.Path) -> None:
    rows = []
    for path in sorted(root.iterdir()):
        if path.is_file() and path.name != "SHA256.tsv":
            rows.append(f"{sha256(path)}  {path.name}")
    (root / "SHA256.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")


def audit(root: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest = read_json(root / "MANIFEST.json")
    stage = read_json(root / "WINDOW_MANIFEST.json")["stage"]
    technical, geometry, diagnostics = RAW.BASE.AUTH.TETRA._technical(
        scratch / "stage.out", scratch / "stage.rst7", manifest
    )
    source = manifest["source_geometry"]
    response: dict[str, Any] = {"all": False, "active_coordinates": []}
    guard: dict[str, Any] = {"pass": False, "checks": {}}
    effect: dict[str, Any] | None = None
    absolute = False
    eligible = False
    if technical:
        response = RAW._source_response(
            source, geometry, stage["targets"], stage["active_coordinates"]
        )
        guard = RAW.first_window_guard(geometry)
        absolute = bool(response["all"] and guard["pass"])
        item = {
            "source_geometry": source,
            "target_geometry": stage["targets"],
            "final_geometry": geometry,
            "response_from_raw_source": response,
        }
        baseline = read_json(pathlib.Path(manifest["source"]["matched_baseline_result"]))
        effect = RAW.force_effect_vs_baseline_v2(item, baseline)
        eligible = bool(absolute and effect["all"])
        shutil.copy2(scratch / "stage.rst7", root / "stage.rst7")
    for name in ("stage.in", "restraints.RST", "stage.mdinfo"):
        path = scratch / name
        if path.is_file():
            shutil.copy2(path, root / name)
    stage_out = scratch / "stage.out"
    if stage_out.is_file():
        tail = stage_out.read_text(encoding="utf-8", errors="replace").splitlines()[-200:]
        (root / "ENGINE_TAIL.txt").write_text("\n".join(tail) + "\n", encoding="utf-8")

    gate = (
        "NOT_EVALUATED_TECHNICAL_A1_RAW_NAC_CONCERTED_RECALIBRATION"
        if not technical else
        "PASS_CONCERTED_FIRST_WINDOW_RESPONSE_V2"
        if eligible else
        "FAIL_CONCERTED_FIRST_WINDOW_RESPONSE_OR_GUARD_V2"
    )
    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_A1_RAW_NAC_CONCERTED_RECALIBRATION"
            if technical else
            "NOT_EVALUATED_TECHNICAL_A1_RAW_NAC_CONCERTED_RECALIBRATION"
        ),
        "technical_complete": technical,
        "scientific_status": SCIENTIFIC_STATUS,
        "scientific_gate": gate,
        "task": manifest["task"],
        "source_geometry": source,
        "target_geometry": stage["targets"],
        "final_geometry": geometry if technical else None,
        "source_qPT_A": RAW._qpt(source),
        "target_qPT_A": RAW._qpt(stage["targets"]),
        "final_qPT_A": RAW._qpt(geometry) if technical else None,
        "response_from_raw_source": response,
        "first_window_guard": guard,
        "eligible_by_absolute_response": absolute,
        "force_effect_vs_matched_baseline_v2": effect,
        "eligible_for_inherited_chain_v2": eligible,
        "first_window_only": True,
        "automatic_downstream_action": "NONE",
        "diagnostics": diagnostics,
        "restrained_structures_are_not_ts": True,
    }
    write_json(root / "RESULT.json", result)
    write_json(root / ("PASS.json" if technical else "NOT_EVALUATED.json"), result)
    _write_hashes(root)
    return result


def select_weakest_by_seed(rows: list[Mapping[str, Any]]) -> list[dict[str, Any]]:
    selected = []
    for seed in (26723, 26737):
        candidates = [
            row for row in rows
            if int(row["task"]["seed"]) == seed
            and row.get("eligible_for_inherited_chain_v2") is True
        ]
        candidates.sort(
            key=lambda row: (
                int(row["task"]["force_scale"]),
                float(row["task"]["window_fraction"]),
            )
        )
        winner = candidates[0] if candidates else None
        selected.append({
            "seed": seed,
            "mode": "FULLY_CONCERTED_RAW_NAC",
            "selected_variant": winner["task"]["variant"] if winner else None,
            "selected_force_scale": winner["task"]["force_scale"] if winner else None,
            "selected_window_fraction": winner["task"]["window_fraction"] if winner else None,
            "selected_task_index": winner["task"]["task_index"] if winner else None,
        })
    return selected


def merge_if_ready(output_root: pathlib.Path, array_job: str) -> bool:
    paths = [
        output_root / f"attempt_{array_job}_{index}" / "RESULT.json"
        for index in range(ARRAY_TASKS)
    ]
    if not all(path.is_file() for path in paths):
        return False
    rows = [read_json(path) for path in paths]
    selected = select_weakest_by_seed(rows)
    technical = sum(row.get("technical_complete") is True for row in rows)
    selected_seeds = sum(row["selected_variant"] is not None for row in selected)
    payload = {
        "schema_version": 1,
        "status": (
            "NOT_EVALUATED_TECHNICAL_A1_RAW_NAC_CONCERTED_RECALIBRATION"
            if technical != ARRAY_TASKS else
            "PASS_CONCERTED_RECALIBRATION_REPRODUCED"
            if selected_seeds == 2 else
            "PARTIAL_CONCERTED_RECALIBRATION"
        ),
        "denominator_tasks": ARRAY_TASKS,
        "technical_pass_tasks": technical,
        "qualified_tasks_v2": sum(
            row.get("eligible_for_inherited_chain_v2") is True for row in rows
        ),
        "selected_seeds": selected_seeds,
        "selected_by_seed": selected,
        "per_task": rows,
        "next_action": (
            "STRICT_INHERITED_CONCERTED_WINDOWS_FROM_SELECTED_VARIANTS"
            if selected_seeds == 2 else
            "NO_CONCERTED_INHERITANCE_UNDER_TESTED_LOCAL_RECALIBRATION"
        ),
        "scientific_status": SCIENTIFIC_STATUS,
        "automatic_downstream_action": "NONE",
    }
    audit_path = (
        output_root / "audit"
        / f"nylc_a1_raw_nac_concerted_recalibration_{array_job}.json"
    )
    if audit_path.exists() and read_json(audit_path) != payload:
        raise FileExistsError(audit_path)
    if not audit_path.exists():
        write_json(audit_path, payload)
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode", required=True,
        choices=("describe", "initialize", "prepare", "audit", "merge-if-ready"),
    )
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--root", type=pathlib.Path)
    parser.add_argument("--scratch", type=pathlib.Path)
    parser.add_argument("--github-commit", default="unknown")
    parser.add_argument("--output-root", type=pathlib.Path)
    parser.add_argument("--array-job")
    args = parser.parse_args()
    if args.mode == "describe":
        print(json.dumps(describe(), sort_keys=True))
    elif args.mode == "initialize":
        initialize(args.task_index, args.root.resolve(), args.github_commit)
    elif args.mode == "prepare":
        prepare(args.root.resolve(), args.scratch.resolve())
    elif args.mode == "audit":
        audit(args.root.resolve(), args.scratch.resolve())
    else:
        merge_if_ready(args.output_root.resolve(), args.array_job)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
