#!/usr/bin/env python3
"""Four-anchor, bidirectional restrained gap fill for NylC A1 Step1.

R, I1, I2, and P are chemical anchors.  Every anchor is propagated toward
both endpoint basins for two independent seeds.  The resulting 16 inherited
chains are restrained path probes only, never TS, committor, PMF, barrier, or
mechanism evidence.
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
REVERSE_PATH = HERE / "prepare_audit_nylc_a1_step1_product_reverse_boundary.py"


def _load(name: str, path: pathlib.Path) -> Any:
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


REVERSE = _load("_a1_reverse_for_four_anchor", REVERSE_PATH)
AC = REVERSE.AC
BASE = REVERSE.BASE
TETRA = REVERSE.TETRA
PRMTOP = REVERSE.PRMTOP
REACTIVE = AC.REACTIVE_ATOMS

TASK_ROOT = pathlib.Path(
    "/work/home/acshdt1dks/nylon_pa66_scnet_20260708/"
    "l4_nac_to_l2_rebalance_20260723"
)
REACTANT_ROOT = (
    TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_attack_inherited_chain"
)
REACTANT_SOURCES = (
    {
        "seed_index": 0,
        "seed": "seed26723",
        "restart": REACTANT_ROOT / "attempt_62136360_0/terminal_endpoint.rst7",
        "restart_sha256": "9d93b60aee9e8d6fa97493396d757f3bf4d0a47595591968e3debaebb2640e86",
    },
    {
        "seed_index": 1,
        "seed": "seed26737",
        "restart": REACTANT_ROOT / "attempt_62156801_1/terminal_endpoint.rst7",
        "restart_sha256": "a12c018e195c32b34239004aea36ebc80de8275c9f9332b44084ad8d88d4db0e",
    },
)
PRODUCT_SOURCES = tuple(
    {
        "seed_index": item["seed_index"],
        "seed": item["seed"],
        "restart": item["restart"],
        "restart_sha256": item["restart_sha256"],
    }
    for item in REVERSE.SOURCES
)

ANCHORS = ("R", "I1", "I2", "P")
DIRECTIONS = ("REACTANT", "PRODUCT")
ANCHOR_TARGETS = {
    "R": {
        "attack_A": 2.18, "cn_A": 1.38, "nalpha_hg1_A": 1.04,
        "hg1_n3_A": 2.45, "c12_o2_A": 1.26,
        "proton_site": "NALPHA", "tetrahedral_required": False,
    },
    "I1": {
        "attack_A": 1.60, "cn_A": 1.50, "nalpha_hg1_A": 1.05,
        "hg1_n3_A": 2.10, "c12_o2_A": 1.34,
        "proton_site": "NALPHA", "tetrahedral_required": True,
    },
    "I2": {
        "attack_A": 1.60, "cn_A": 1.55, "nalpha_hg1_A": 2.10,
        "hg1_n3_A": 1.05, "c12_o2_A": 1.34,
        "proton_site": "N3", "tetrahedral_required": True,
    },
    "P": {
        "attack_A": 1.42, "cn_A": 3.10, "nalpha_hg1_A": 2.70,
        "hg1_n3_A": 1.02, "c12_o2_A": 1.23,
        "proton_site": "N3", "tetrahedral_required": False,
    },
}
SOURCE_ANCHOR = {"R": "R", "I1": "R", "I2": "P", "P": "P"}
POINTS_PER_SEGMENT = 4
SCIENTIFIC_STATUS = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"
TECHNICAL_PASS = "PASS_TECHNICAL_A1_FOUR_ANCHOR_BIDIRECTIONAL_GAPFILL"
TECHNICAL_FAIL = "NOT_EVALUATED_TECHNICAL_A1_FOUR_ANCHOR_BIDIRECTIONAL_GAPFILL"


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: pathlib.Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < 16:
        raise ValueError("task index must be 0..15")
    direction_index = index % 2
    anchor_index = (index // 2) % 4
    seed_index = index // 8
    return {
        "task_index": index,
        "seed_index": seed_index,
        "seed": REACTANT_SOURCES[seed_index]["seed"],
        "anchor": ANCHORS[anchor_index],
        "direction": DIRECTIONS[direction_index],
    }


def _interpolate(left: str, right: str, fraction: float) -> dict[str, float]:
    result = {}
    for key in ("attack_A", "cn_A", "nalpha_hg1_A", "hg1_n3_A", "c12_o2_A"):
        a = float(ANCHOR_TARGETS[left][key])
        b = float(ANCHOR_TARGETS[right][key])
        result[key] = a + (b - a) * float(fraction)
    return result


def _segment(left: str, right: str, role: str) -> list[dict[str, Any]]:
    records = []
    for step in range(1, POINTS_PER_SEGMENT + 1):
        fraction = step / POINTS_PER_SEGMENT
        records.append({
            "from_anchor": left,
            "to_anchor": right,
            "fraction": fraction,
            "role": role,
            "anchor_checkpoint": fraction == 1.0,
            "targets": _interpolate(left, right, fraction),
            "maxcyc": 600,
            "ncyc": 150,
            "force_bonds": 20.0,
            "force_proton": 15.0,
            "force_carbonyl": 12.0,
        })
    return records


def stage_specs(task_index: int) -> list[dict[str, Any]]:
    task = task_spec(task_index)
    anchor = task["anchor"]
    source_anchor = SOURCE_ANCHOR[anchor]
    stages: list[dict[str, Any]] = []
    if source_anchor != anchor:
        direction = 1 if ANCHORS.index(anchor) > ANCHORS.index(source_anchor) else -1
        current = ANCHORS.index(source_anchor)
        while current != ANCHORS.index(anchor):
            nxt = current + direction
            stages.extend(_segment(ANCHORS[current], ANCHORS[nxt], "BUILD_START_ANCHOR"))
            current = nxt
    destination = "R" if task["direction"] == "REACTANT" else "P"
    if anchor == destination:
        stages.append({
            "from_anchor": anchor,
            "to_anchor": anchor,
            "fraction": 1.0,
            "role": "ENDPOINT_BASIN_CONTROL",
            "anchor_checkpoint": True,
            "targets": {
                key: float(ANCHOR_TARGETS[anchor][key])
                for key in ("attack_A", "cn_A", "nalpha_hg1_A", "hg1_n3_A", "c12_o2_A")
            },
            "maxcyc": 400,
            "ncyc": 100,
            "force_bonds": 10.0,
            "force_proton": 8.0,
            "force_carbonyl": 6.0,
        })
    else:
        direction = 1 if ANCHORS.index(destination) > ANCHORS.index(anchor) else -1
        current = ANCHORS.index(anchor)
        while current != ANCHORS.index(destination):
            nxt = current + direction
            stages.extend(_segment(ANCHORS[current], ANCHORS[nxt], "PROPAGATE_TO_ENDPOINT"))
            current = nxt
    for stage_index, stage in enumerate(stages):
        stage["stage_index"] = stage_index
        stage["target_anchor"] = stage["to_anchor"] if stage["anchor_checkpoint"] else None
    return stages


def _source_for(task: Mapping[str, Any]) -> dict[str, Any]:
    source_anchor = SOURCE_ANCHOR[task["anchor"]]
    collection = REACTANT_SOURCES if source_anchor == "R" else PRODUCT_SOURCES
    return dict(collection[int(task["seed_index"])])


def _geometry(path: pathlib.Path, manifest: Mapping[str, Any]) -> dict[str, Any]:
    import parmed as pmd
    structure = pmd.load_file(str(PRMTOP), xyz=str(path))
    AC._validate_structure(structure)
    return AC.geometry_from_coordinates(
        AC._coordinate_map(structure),
        manifest["qm_contract"]["qm_heavy_atom_indices"],
    )


def _finite_geometry(geometry: Mapping[str, Any]) -> bool:
    return bool(geometry) and all(
        isinstance(value, (int, float)) and math.isfinite(float(value))
        for key, value in geometry.items()
        if key != "hg1_nearest_qm_heavy_atom"
    )


def anchor_gate(anchor: str, geometry: Mapping[str, Any]) -> dict[str, bool]:
    target = ANCHOR_TARGETS[anchor]
    nearest = geometry.get("hg1_nearest_qm_heavy_atom")
    checks = {
        "finite": _finite_geometry(geometry),
        "attack_residual_le_0p30": abs(geometry.get("attack_A", 99.0) - target["attack_A"]) <= 0.30,
        "cn_residual_le_0p35": abs(geometry.get("c12_n3_A", 99.0) - target["cn_A"]) <= 0.35,
        "carbonyl_residual_le_0p15": abs(geometry.get("c12_o2_A", 99.0) - target["c12_o2_A"]) <= 0.15,
        "proton_site": (
            nearest == REACTIVE["nalpha"]
            if target["proton_site"] == "NALPHA"
            else nearest == REACTIVE["n3"]
        ),
        "proton_bonded": (
            geometry.get("nalpha_hg1_A", 99.0) <= 1.25
            if target["proton_site"] == "NALPHA"
            else geometry.get("hg1_n3_A", 99.0) <= 1.25
        ),
        "not_overcompressed": (
            geometry.get("attack_A", 0.0) >= 1.25
            and geometry.get("c12_n3_A", 0.0) >= 1.25
        ),
    }
    tetra = (
        geometry.get("c12_o2_A", 0.0) >= 1.28
        or geometry.get("c12_reactant_plane_out_of_plane_A", 0.0) >= 0.15
    )
    checks["tetrahedral_response"] = (
        tetra if target["tetrahedral_required"] else True
    )
    checks["all"] = all(checks.values())
    return checks


def restraints(stage: Mapping[str, Any]) -> str:
    t = stage["targets"]
    return "".join((
        AC.distance_restraint(REACTIVE["og1"], REACTIVE["c12"], t["attack_A"], stage["force_bonds"]),
        AC.distance_restraint(REACTIVE["c12"], REACTIVE["n3"], t["cn_A"], stage["force_bonds"]),
        AC.distance_restraint(REACTIVE["nalpha"], REACTIVE["hg1"], t["nalpha_hg1_A"], stage["force_proton"]),
        AC.distance_restraint(REACTIVE["hg1"], REACTIVE["n3"], t["hg1_n3_A"], stage["force_proton"]),
        AC.distance_restraint(REACTIVE["c12"], REACTIVE["o2"], t["c12_o2_A"], stage["force_carbonyl"]),
    ))


def minimization_input(task: Mapping[str, Any], stage: Mapping[str, Any], qmmask: str) -> str:
    return f"""NylC A1 four-anchor gapfill task {task['task_index']} stage {stage['stage_index']}
&cntrl
  imin=1, ntmin=2, maxcyc={stage['maxcyc']}, ncyc={stage['ncyc']}, dx0=0.005,
  ntb=1, cut=10.0, ntpr=5, ntxo=1, ifqnt=1,
  ntr=1, nmropt=1, restraint_wt=1.0,
  restraintmask='{BASE.NON_QM_SOLUTE_HEAVY_MASK}', drms=0.10,
/
{BASE.qmmm_block(qmmask)}&wt type='END' /
DISANG=restraints.RST
DUMPAVE=restraint.dat
"""


def initialize(task_index: int, root: pathlib.Path, start: pathlib.Path, commit: str) -> None:
    if root.exists():
        raise FileExistsError(root)
    task = task_spec(task_index)
    source = _source_for(task)
    if sha256(source["restart"]) != source["restart_sha256"]:
        raise ValueError("fixed source restart SHA changed")
    source_manifest, qmmask, _ = REVERSE.validate_authority(task["seed_index"])
    if sha256(PRMTOP) != BASE.EXPECTED_PRMTOP_SHA256:
        raise ValueError("frozen prmtop SHA changed")
    shutil.copy2(source["restart"], start)
    root.mkdir(parents=True, exist_ok=False)
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_FOUR_ANCHOR_BIDIRECTIONAL_GAPFILL",
        "scientific_status": SCIENTIFIC_STATUS,
        "github_commit": commit,
        "task": task,
        "source": {
            "restart": str(source["restart"]),
            "restart_sha256": source["restart_sha256"],
            "source_anchor": SOURCE_ANCHOR[task["anchor"]],
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": sha256(PRMTOP),
        "qm_contract": dict(source_manifest["qm_contract"]),
        "start_restart": str(start),
        "start_restart_sha256": sha256(start),
        "source_geometry": _geometry(start, {
            "qm_contract": source_manifest["qm_contract"]
        }),
        "stage_specs": stage_specs(task_index),
        "stages": [],
        "guard_stop": False,
        "restrained_structures_are_not_ts": True,
        "automatic_downstream_action": "NONE",
    }
    write_json(root / "MANIFEST.json", manifest)


def prepare_stage(root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path,
                  input_rst7: pathlib.Path, stage_index: int) -> None:
    manifest = json.loads((root / "MANIFEST.json").read_text(encoding="utf-8"))
    specs = manifest["stage_specs"]
    if stage_index != len(manifest["stages"]) or not 0 <= stage_index < len(specs):
        raise ValueError("stage order is not strict")
    previous_sha = (
        manifest["start_restart_sha256"]
        if stage_index == 0 else manifest["stages"][-1]["output_restart_sha256"]
    )
    if sha256(input_rst7) != previous_sha:
        raise ValueError("stage input breaks SHA inheritance")
    output.mkdir(parents=True, exist_ok=False)
    scratch.mkdir(parents=True, exist_ok=False)
    stage = specs[stage_index]
    (scratch / "stage.in").write_text(
        minimization_input(manifest["task"], stage, manifest["qm_contract"]["qmmask"]),
        encoding="utf-8",
    )
    (scratch / "restraints.RST").write_text(restraints(stage), encoding="utf-8")
    write_json(output / "STAGE_MANIFEST.json", {
        "schema_version": 1,
        "stage_spec": stage,
        "input_restart": str(input_rst7),
        "input_restart_sha256": previous_sha,
        "restrained": True,
        "shooting_forbidden": True,
    })


def audit_stage(root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest_path = root / "MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    prepared = json.loads((output / "STAGE_MANIFEST.json").read_text(encoding="utf-8"))
    scratch_restart = scratch / "stage.rst7"
    technical, geometry, diagnostics = TETRA._technical(
        scratch / "stage.out", scratch_restart, manifest
    )
    persistent_restart = output / "stage.rst7"
    if technical:
        shutil.copy2(scratch_restart, persistent_restart)
    stage = prepared["stage_spec"]
    target_anchor = stage.get("target_anchor")
    checkpoint = anchor_gate(target_anchor, geometry) if technical and target_anchor else None
    chemical = REVERSE.chemical_guard(geometry) if technical else {"all": False}
    guard_stop = bool(
        technical and (
            not chemical["all"]
            or (target_anchor in ("I1", "I2") and not checkpoint["all"])
        )
    )
    result = {
        "schema_version": 1,
        "status": "PASS_TECHNICAL_A1_FOUR_ANCHOR_STAGE" if technical else "NOT_EVALUATED_TECHNICAL_A1_FOUR_ANCHOR_STAGE",
        "technical_pass": technical,
        "scientific_status": SCIENTIFIC_STATUS,
        "stage_spec": stage,
        "input_restart_sha256": prepared["input_restart_sha256"],
        "output_restart": str(persistent_restart),
        "output_restart_sha256": sha256(persistent_restart) if persistent_restart.is_file() else None,
        "geometry": geometry,
        "chemical_guard": chemical,
        "anchor_gate": checkpoint,
        "guard_stop": guard_stop,
        "restrained": True,
        "shooting_forbidden": True,
        "diagnostics": diagnostics,
    }
    write_json(output / "RESULT.json", result)
    write_json(output / ("PASS.json" if technical else "NOT_EVALUATED.json"), result)
    manifest["stages"].append(result)
    manifest["guard_stop"] = guard_stop
    write_json(manifest_path, manifest)
    if not technical:
        raise RuntimeError("stage did not technically pass")
    return result


def finalize(root: pathlib.Path, stop_reason: str, technical_failure: bool = False) -> dict[str, Any]:
    manifest = json.loads((root / "MANIFEST.json").read_text(encoding="utf-8"))
    stages = manifest["stages"]
    inheritance = True
    previous = manifest["start_restart_sha256"]
    for stage in stages:
        inheritance = inheritance and stage.get("input_restart_sha256") == previous
        previous = stage.get("output_restart_sha256")
    expected = len(manifest["stage_specs"])
    technical = bool(
        not technical_failure
        and stages
        and all(stage.get("technical_pass") is True for stage in stages)
        and inheritance
        and (
            (stop_reason == "COMPLETED_ALL_STAGES" and len(stages) == expected)
            or (stop_reason == "SCIENTIFIC_GUARD_STOP" and manifest["guard_stop"])
        )
    )
    destination = "R" if manifest["task"]["direction"] == "REACTANT" else "P"
    destination_gate = (
        anchor_gate(destination, stages[-1]["geometry"]) if stages else {"all": False}
    )
    scientific_gate = (
        "NOT_EVALUATED_TECHNICAL_FOUR_ANCHOR_GAPFILL"
        if not technical else
        "PASS_GUIDED_ROUTE_REACHED_DESTINATION"
        if stop_reason == "COMPLETED_ALL_STAGES" and destination_gate["all"] else
        "FAIL_GUIDED_ROUTE_DID_NOT_REACH_DESTINATION"
    )
    result = {
        "schema_version": 1,
        "status": TECHNICAL_PASS if technical else TECHNICAL_FAIL,
        "technical_complete": technical,
        "scientific_status": SCIENTIFIC_STATUS,
        "scientific_gate": scientific_gate,
        "task": manifest["task"],
        "stop_reason": stop_reason,
        "stage_denominator_expected": expected,
        "stage_denominator_completed": len(stages),
        "inheritance_sha256_verified": inheritance,
        "destination_anchor": destination,
        "destination_gate": destination_gate,
        "restrained_structures_are_not_ts": True,
        "automatic_downstream_action": "NONE",
    }
    write_json(root / "RESULT.json", result)
    write_json(root / ("PASS.json" if technical else "NOT_EVALUATED.json"), result)
    records = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.name != "SHA256.tsv" and path.suffix in (".json", ".rst7"):
            records.append(f"{sha256(path)}  {path.relative_to(root)}")
    (root / "SHA256.tsv").write_text("\n".join(records) + "\n", encoding="utf-8")
    return result


def merge_if_ready(output_root: pathlib.Path, array_job: str) -> bool:
    audit = output_root / "audit" / f"nylc_a1_four_anchor_gapfill_{array_job}.json"
    audit.parent.mkdir(parents=True, exist_ok=True)
    lock_path = audit.with_suffix(".lock")
    with lock_path.open("a+", encoding="utf-8") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        try:
            paths = [output_root / f"attempt_{array_job}_{i}" / "RESULT.json" for i in range(16)]
            if not all(path.is_file() for path in paths):
                return False
            results = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
            technical = all(item.get("technical_complete") is True for item in results)
            reached = sum(item.get("scientific_gate") == "PASS_GUIDED_ROUTE_REACHED_DESTINATION" for item in results)
            payload = {
                "schema_version": 1,
                "status": (
                    "NOT_EVALUATED_TECHNICAL_A1_FOUR_ANCHOR_GAPFILL"
                    if not technical else
                    "PASS_GUIDED_FOUR_ANCHOR_GAP_COVERAGE"
                    if reached == 16 else
                    "PARTIAL_GUIDED_FOUR_ANCHOR_GAP_COVERAGE"
                ),
                "denominator_tasks": 16,
                "technical_pass_tasks": sum(item.get("technical_complete") is True for item in results),
                "guided_destination_pass_tasks": reached,
                "per_task": sorted(results, key=lambda item: item["task"]["task_index"]),
                "scientific_status": SCIENTIFIC_STATUS,
                "restrained_structures_are_not_ts": True,
                "automatic_downstream_action": "NONE",
            }
            if audit.exists() and json.loads(audit.read_text(encoding="utf-8")) != payload:
                raise FileExistsError(audit)
            if not audit.exists():
                write_json(audit, payload)
            return True
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "anchors": list(ANCHORS),
        "directions": list(DIRECTIONS),
        "seed_count": 2,
        "array_task_count": 16,
        "array_mapping": "task=((seed*4+anchor)*2+direction)",
        "points_per_adjacent_segment": POINTS_PER_SEGMENT,
        "stage_restart_persisted": True,
        "restrained_structures_are_not_ts": True,
        "automatic_downstream_action": "NONE",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=(
        "describe", "describe-task", "initialize", "prepare-stage",
        "audit-stage", "finalize", "merge-if-ready",
    ))
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--stage-index", type=int)
    parser.add_argument("--root", type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--output-root", type=pathlib.Path)
    parser.add_argument("--scratch", type=pathlib.Path)
    parser.add_argument("--start-rst7", type=pathlib.Path)
    parser.add_argument("--input-rst7", type=pathlib.Path)
    parser.add_argument("--github-commit", default="unknown")
    parser.add_argument("--stop-reason", default="UNKNOWN")
    parser.add_argument("--array-job")
    parser.add_argument("--technical-failure", action="store_true")
    args = parser.parse_args()
    if args.mode == "describe":
        print(json.dumps(describe(), sort_keys=True))
    elif args.mode == "describe-task":
        payload = task_spec(args.task_index)
        payload["stage_specs"] = stage_specs(args.task_index)
        print(json.dumps(payload, sort_keys=True))
    elif args.mode == "initialize":
        initialize(args.task_index, args.root.resolve(), args.start_rst7.resolve(), args.github_commit)
    elif args.mode == "prepare-stage":
        prepare_stage(args.root.resolve(), args.output.resolve(), args.scratch.resolve(), args.input_rst7.resolve(), args.stage_index)
    elif args.mode == "audit-stage":
        audit_stage(args.root.resolve(), args.output.resolve(), args.scratch.resolve())
    elif args.mode == "finalize":
        finalize(args.root.resolve(), args.stop_reason, args.technical_failure)
    else:
        merge_if_ready(args.output_root.resolve(), args.array_job)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
