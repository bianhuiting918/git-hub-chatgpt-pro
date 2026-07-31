#!/usr/bin/env python3
"""Seed-specific raw-NAC Step1 attack/carbonyl force calibration."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import pathlib
import re
import shutil
from typing import Any, Mapping

HERE = pathlib.Path(__file__).resolve().parent
RAW_FORWARD_PATH = HERE / "prepare_audit_nylc_a1_step1_raw_nac_forward_calibration.py"


def _load(name: str, path: pathlib.Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RAW_FORWARD = _load("_a1_raw_forward", RAW_FORWARD_PATH)
FWD = RAW_FORWARD.FWD
RAW = RAW_FORWARD.RAW
BASE = RAW_FORWARD.BASE
REACTIVE = RAW_FORWARD.REACTIVE
TASK_ROOT = RAW_FORWARD.TASK_ROOT
PRMTOP = RAW_FORWARD.PRMTOP
PRMTOP_SHA256 = RAW_FORWARD.PRMTOP_SHA256
QM_CONTRACT = dict(RAW_FORWARD.QM_CONTRACT)
SOURCES = {
    seed: {
        key: value for key, value in source.items()
        if key != "contract_attempt"
    }
    for seed, source in RAW_FORWARD.SOURCES.items()
}

ARRAY_TASKS = 18
MPI_RANKS = 8
ATTACK_FORCE_SCALES = (1, 2, 4)
CARBONYL_FORCE_SCALES = (1, 2, 4)
WEAK_ANGLE_GUARD_FORCE = 2.0
WEAK_PT_GUARD_FORCE = 3.0
FIRST_WINDOW_ATTACK_DELTA_A = -0.04
FIRST_WINDOW_CARBONYL_DELTA_A = 0.03
OUTPUT_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/qmmm"
    / "a1_step1_raw_nac_preorganized_calibration"
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


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..17")
    seed = (26723, 26737)[index // 9]
    local = index % 9
    attack_scale = ATTACK_FORCE_SCALES[local // 3]
    carbonyl_scale = CARBONYL_FORCE_SCALES[local % 3]
    return {
        "task_index": index,
        "seed": seed,
        "seed_index": RAW_FORWARD.SOURCES[seed]["seed_index"],
        "candidate": RAW_FORWARD.SOURCES[seed]["candidate"],
        "mode": "ADDITION_FIRST_PREORGANIZED_RAW_NAC",
        "mechanism": "ADDITION_FIRST",
        "attack_scale": attack_scale,
        "carbonyl_scale": carbonyl_scale,
        "source_basin": "HASH_FIXED_RAW_UNBIASED_NAC",
        "first_window_only": True,
    }


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "mapping": "2 seeds x 3 attack-force scales x 3 carbonyl-force scales",
        "attack_force_scales": list(ATTACK_FORCE_SCALES),
        "carbonyl_force_scales": list(CARBONYL_FORCE_SCALES),
        "weak_angle_guard_force": WEAK_ANGLE_GUARD_FORCE,
        "weak_pt_guard_force": WEAK_PT_GUARD_FORCE,
        "original_raw_nac_restart_shas": {
            str(seed): source["start_rst7_sha256"]
            for seed, source in SOURCES.items()
        },
        "qm_contract": dict(QM_CONTRACT),
        "prmtop_sha256": PRMTOP_SHA256,
        "strict_serial_restart_inheritance": False,
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }


def _qpt(geometry: Mapping[str, Any]) -> float:
    return float(geometry["nalpha_hg1_A"]) - float(geometry["hg1_n3_A"])


def stage_spec(
    spec: Mapping[str, Any], source: Mapping[str, Any]
) -> dict[str, Any]:
    targets = {
        "attack_A": float(source["attack_A"]) + FIRST_WINDOW_ATTACK_DELTA_A,
        "c12_n3_A": float(source["c12_n3_A"]),
        "c12_o2_A": float(source["c12_o2_A"]) + FIRST_WINDOW_CARBONYL_DELTA_A,
        "nalpha_hg1_A": float(source["nalpha_hg1_A"]),
        "hg1_n3_A": float(source["hg1_n3_A"]),
        "attack_angle_deg": float(source["attack_angle_deg"]),
    }
    return {
        "window_index": 0,
        "targets": targets,
        "active_coordinates": ["attack", "carbonyl"],
        "force_attack": (
            FWD.CAL.BASE_FORCES["force_bond"] * float(spec["attack_scale"])
        ),
        "force_carbonyl": (
            FWD.CAL.BASE_FORCES["force_carbonyl"]
            * float(spec["carbonyl_scale"])
        ),
        "force_angle_guard": WEAK_ANGLE_GUARD_FORCE,
        "force_pt_guard": WEAK_PT_GUARD_FORCE,
        "force_bond": (
            FWD.CAL.BASE_FORCES["force_bond"] * float(spec["attack_scale"])
        ),
        "force_pt": WEAK_PT_GUARD_FORCE,
        "maxcyc": 2200,
        "ncyc": 550,
    }


def _angle_restraint(
    atom1: int, atom2: int, atom3: int, target: float, force: float
) -> str:
    return (
        f"&rst iat={atom1},{atom2},{atom3}, "
        f"r1={max(0.0, target - 25.0):.3f}, "
        f"r2={target - 2.0:.3f}, r3={target + 2.0:.3f}, "
        f"r4={min(180.0, target + 25.0):.3f}, "
        f"rk2={force:.1f}, rk3={force:.1f}, /\n"
    )


def restraints(stage: Mapping[str, Any]) -> str:
    target = stage["targets"]
    rows = [
        FWD.CAL.BRIDGE._tight_distance_restraint(
            REACTIVE["og1"], REACTIVE["c12"], target["attack_A"],
            stage["force_attack"],
        ),
        FWD.CAL.BRIDGE._tight_distance_restraint(
            REACTIVE["c12"], REACTIVE["o2"], target["c12_o2_A"],
            stage["force_carbonyl"],
        ),
        _angle_restraint(
            REACTIVE["og1"], REACTIVE["c12"], REACTIVE["n3"],
            target["attack_angle_deg"], stage["force_angle_guard"],
        ),
        FWD.CAL.BRIDGE._tight_distance_restraint(
            REACTIVE["nalpha"], REACTIVE["hg1"], target["nalpha_hg1_A"],
            stage["force_pt_guard"],
        ),
        FWD.CAL.BRIDGE._tight_distance_restraint(
            REACTIVE["hg1"], REACTIVE["n3"], target["hg1_n3_A"],
            stage["force_pt_guard"],
        ),
    ]
    return "".join(rows)


def preorganization_guard(geometry: Mapping[str, Any]) -> dict[str, Any]:
    nearest = geometry.get(
        "hg1_nearest_qm_heavy_atom", geometry.get("proton_nearest")
    )
    qpt = _qpt(geometry)
    values = (
        geometry["attack_A"], geometry["c12_n3_A"], geometry["c12_o2_A"],
        geometry["attack_angle_deg"], qpt,
    )
    checks = {
        "finite_geometry": all(math.isfinite(float(value)) for value in values),
        "attack_not_overcompressed": float(geometry["attack_A"]) >= 1.25,
        "cn_not_overcompressed": float(geometry["c12_n3_A"]) >= 1.25,
        "carbonyl_safe": 1.15 <= float(geometry["c12_o2_A"]) <= 1.55,
        "attack_angle_preorganized": (
            80.0 <= float(geometry["attack_angle_deg"]) <= 140.0
        ),
        "reactant_side_qpt": qpt <= -0.30,
        "proton_nearest_nalpha": nearest in (
            "Nalpha", REACTIVE["nalpha"]
        ),
    }
    return {"pass": all(checks.values()), "checks": checks, "qPT_A": qpt}


def _write_hashes(root: pathlib.Path) -> None:
    rows = []
    for path in sorted(root.iterdir()):
        if path.is_file() and path.name != "SHA256.tsv":
            rows.append(f"{sha256(path)}  {path.name}")
    (root / "SHA256.tsv").write_text(
        "\n".join(rows) + "\n", encoding="utf-8"
    )


def initialize(
    task_index: int, root: pathlib.Path, commit: str
) -> dict[str, Any]:
    if root.exists():
        raise FileExistsError(root)
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("runtime is not bound to a full Git commit")
    spec = task_spec(task_index)
    source = RAW_FORWARD.SOURCES[spec["seed"]]
    authority = dict(RAW.BASE.SOURCES[spec["seed_index"]])
    qmmask, authority_gro = RAW.BASE.validate_authority(authority)
    source_gro = pathlib.Path(source["source_gro"])
    if authority_gro.resolve() != source_gro.resolve():
        raise ValueError("raw NAC source path differs from frozen authority")
    if sha256(source_gro) != source["source_gro_sha256"]:
        raise ValueError("raw NAC source GRO SHA changed")
    if sha256(PRMTOP) != PRMTOP_SHA256:
        raise ValueError("frozen prmtop SHA changed")

    import parmed as pmd
    structure = pmd.load_file(str(PRMTOP), xyz=str(source_gro))
    if len(structure.atoms) != RAW.BASE.EXPECTED_SYSTEM_ATOMS:
        raise ValueError("coordinate transplant changed atom count")
    root.mkdir(parents=True, exist_ok=False)
    start = root / "source.rst7"
    structure.save(str(start), overwrite=False)
    if sha256(start) != source["start_rst7_sha256"]:
        raise ValueError("deterministic raw NAC restart SHA changed")

    contract_manifest = (
        RAW_FORWARD.CONTRACT_ROOT
        / source["contract_attempt"]
        / "ENDPOINT_MANIFEST.json"
    )
    contract = RAW_FORWARD.validate_full_contract(
        BASE.read_json(contract_manifest), qmmask
    )
    geometry = BASE.AUTH._geometry(start, {"qm_contract": contract})
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_RAW_NAC_PREORGANIZED_CALIBRATION",
        "github_commit": commit,
        "task": spec,
        "source": {
            "source_gro": str(source_gro),
            "source_gro_sha256": source["source_gro_sha256"],
            "restart": str(start),
            "restart_sha256": source["start_rst7_sha256"],
            "calibration_restart_reuse": False,
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": sha256(PRMTOP),
        "qm_contract": contract,
        "source_geometry": geometry,
        "scientific_status": SCIENTIFIC_STATUS,
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
    spec = manifest["task"]
    stage = stage_spec(spec, manifest["source_geometry"])
    stage_input = FWD.CAL.BRIDGE.tetra_minimization_input(
        spec, stage, manifest["qm_contract"]["qmmask"]
    )
    (scratch / "stage.in").write_text(stage_input, encoding="utf-8")
    rst = restraints(stage)
    if rst.count("&rst") != 5 or r"\n" in rst:
        raise ValueError("preorganization restraints lost physical-line authority")
    (scratch / "restraints.RST").write_text(rst, encoding="utf-8")
    prepared = {
        "schema_version": 1,
        "stage": stage,
        "input_restart": manifest["source"]["restart"],
        "input_restart_sha256": manifest["source"]["restart_sha256"],
        "reaction_coordinate_restraints": 2,
        "weak_preorganization_restraints": 3,
        "first_window_only": True,
        "automatic_downstream_action": "NONE",
    }
    BASE.write_json(root / "WINDOW_MANIFEST.json", prepared)
    return prepared


def audit(root: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest = BASE.read_json(root / "MANIFEST.json")
    prepared = BASE.read_json(root / "WINDOW_MANIFEST.json")
    stage_out = scratch / "stage.out"
    technical, geometry, diagnostics = BASE.AUTH.TETRA._technical(
        stage_out, scratch / "stage.rst7", manifest
    )
    engine_text = (
        stage_out.read_text(encoding="utf-8", errors="replace")
        if stage_out.is_file() else ""
    )
    numerical = RAW_FORWARD.engine_numerical_health(engine_text)
    technical = bool(technical and numerical["pass"])
    diagnostics = dict(diagnostics)
    diagnostics["numerical_health"] = numerical
    source = manifest["source_geometry"]
    response = {"all": False, "active_coordinates": ["attack", "carbonyl"]}
    guard = {"pass": False, "checks": {}}
    if technical:
        response = {
            "attack": FWD.CAL.BRIDGE._response(
                source["attack_A"], geometry["attack_A"],
                prepared["stage"]["targets"]["attack_A"],
            ),
            "carbonyl": FWD.CAL.BRIDGE._response(
                source["c12_o2_A"], geometry["c12_o2_A"],
                prepared["stage"]["targets"]["c12_o2_A"],
            ),
            "active_coordinates": ["attack", "carbonyl"],
        }
        response["all"] = (
            response["attack"]["pass"] and response["carbonyl"]["pass"]
        )
        guard = preorganization_guard(geometry)
        shutil.copy2(scratch / "stage.rst7", root / "stage.rst7")
    for name in ("stage.in", "restraints.RST", "stage.mdinfo"):
        path = scratch / name
        if path.is_file():
            shutil.copy2(path, root / name)
    if stage_out.is_file():
        tail = engine_text.splitlines()[-240:]
        (root / "ENGINE_TAIL.txt").write_text(
            "\n".join(tail) + "\n", encoding="utf-8"
        )

    qualified = bool(technical and response["all"] and guard["pass"])
    if not technical:
        gate = "NOT_EVALUATED_TECHNICAL_NUMERICAL_OR_ENGINE_FAILURE"
    elif qualified:
        gate = "PASS_PREORGANIZED_FIRST_WINDOW_RESPONSE"
    else:
        gate = "FAIL_PREORGANIZED_FIRST_WINDOW_RESPONSE_OR_GUARD"
    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_A1_RAW_NAC_PREORGANIZED_CALIBRATION"
            if technical else
            "NOT_EVALUATED_TECHNICAL_A1_RAW_NAC_PREORGANIZED_CALIBRATION"
        ),
        "technical_complete": technical,
        "scientific_gate": gate,
        "scientific_status": SCIENTIFIC_STATUS,
        "task": manifest["task"],
        "source_geometry": source,
        "target_geometry": prepared["stage"]["targets"],
        "final_geometry": geometry if technical else None,
        "response_from_raw_source": response,
        "preorganization_guard": guard,
        "eligible_for_inherited_chain": qualified,
        "diagnostics": diagnostics,
        "first_window_only": True,
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }
    BASE.write_json(root / "RESULT.json", result)
    BASE.write_json(
        root / ("PASS.json" if technical else "NOT_EVALUATED.json"), result
    )
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
    selected = []
    for seed in (26723, 26737):
        candidates = [
            row for row in results
            if row["task"]["seed"] == seed
            and row.get("eligible_for_inherited_chain") is True
        ]
        candidates.sort(key=lambda row: (
            row["task"]["attack_scale"] + row["task"]["carbonyl_scale"],
            max(row["task"]["attack_scale"], row["task"]["carbonyl_scale"]),
            row["task"]["attack_scale"],
            row["task"]["carbonyl_scale"],
        ))
        winner = candidates[0] if candidates else None
        selected.append({
            "seed": seed,
            "selected_task_index": (
                winner["task"]["task_index"] if winner else None
            ),
            "selected_attack_scale": (
                winner["task"]["attack_scale"] if winner else None
            ),
            "selected_carbonyl_scale": (
                winner["task"]["carbonyl_scale"] if winner else None
            ),
        })
    technical = sum(
        row.get("technical_complete") is True for row in results
    )
    selected_seeds = sum(
        row["selected_task_index"] is not None for row in selected
    )
    payload = {
        "schema_version": 1,
        "status": (
            "NOT_EVALUATED_TECHNICAL_A1_RAW_NAC_PREORGANIZED_CALIBRATION"
            if technical != ARRAY_TASKS else
            "PASS_BOTH_SEEDS_PREORGANIZED_FORCE_PAIRS"
            if selected_seeds == 2 else
            "PARTIAL_A1_RAW_NAC_PREORGANIZED_CALIBRATION"
        ),
        "denominator_tasks": ARRAY_TASKS,
        "technical_pass_tasks": technical,
        "qualified_tasks": sum(
            row.get("eligible_for_inherited_chain") is True for row in results
        ),
        "selected_seed_specific_force_pairs": selected,
        "per_task": results,
        "next_action": (
            "STRICT_INHERITED_WINDOWS_FROM_SELECTED_SEED_SPECIFIC_FORCE_PAIRS"
            if technical == ARRAY_TASKS and selected_seeds == 2 else
            "NO_INHERITANCE_REVIEW_NUMERICAL_AND_COORDINATE_RESPONSE"
        ),
        "automatic_downstream_action": "NONE",
        "scientific_status": SCIENTIFIC_STATUS,
        "restrained_structures_are_not_ts": True,
    }
    audit_path = (
        output_root / "audit"
        / f"nylc_a1_raw_nac_preorganized_calibration_{array_job}.json"
    )
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    if audit_path.exists() and BASE.read_json(audit_path) != payload:
        raise FileExistsError(audit_path)
    if not audit_path.exists():
        BASE.write_json(audit_path, payload)
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
