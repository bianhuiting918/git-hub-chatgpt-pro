#!/usr/bin/env python3
"""Matched raw-NAC baseline and first-window Step1 forward calibrations."""
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
FWD_PATH = HERE / "prepare_audit_nylc_a1_step1_forward_force_calibration.py"
RAW_PATH = HERE / "prepare_audit_nylc_a1_step1_attack_inherited.py"


def _load(name: str, path: pathlib.Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


FWD = _load("_a1_existing_forward", FWD_PATH)
RAW = _load("_a1_raw_source", RAW_PATH)
BASE = FWD.BASE
REACTIVE = FWD.REACTIVE
TASK_ROOT = FWD.TASK_ROOT
PRMTOP = FWD.PRMTOP

ARRAY_TASKS = 18
MPI_RANKS = 8
FORCE_SCALES = (1, 2, 4, 8)
MODES = (
    "REACTION_COORDINATE_FREE_BASELINE",
    "ADDITION_FIRST_RAW_NAC",
    "FULLY_CONCERTED_RAW_NAC",
)
OUTPUT_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/a1_step1_raw_nac_forward_calibration"
)
SOURCE_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/pt2_preorganized_frame_extraction"
    / "attempt_62112503"
)
CONTRACT_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/a1_step1_acyl_release_md_continuation"
)
SOURCES = {
    26723: {
        "seed_index": 0,
        "candidate": "seed26723_t378_f189",
        "contract_attempt": "attempt_62216380_0",
        "source_gro": SOURCE_ROOT / "seed26723_t378_f189/source.gro",
        "source_gro_sha256": "14477791ce14a35cef0adf9b802b562e091660526ca06de1132f6a74070faf10",
        "start_rst7_sha256": "2440de548c385f092c37f683de7b379ff5b18b6dc16593f7dbc80a9a8a167e14",
    },
    26737: {
        "seed_index": 1,
        "candidate": "seed26737_t676_f338",
        "contract_attempt": "attempt_62216380_1",
        "source_gro": SOURCE_ROOT / "seed26737_t676_f338/source.gro",
        "source_gro_sha256": "a7924184ad3db4c13e0eab4929d46ee99621eacd523e899c0aca39c540350bc8",
        "start_rst7_sha256": "48c3944d3295158b06e96e32e4e07d9bcae9ceba2731f95aae9c1335ff972abe",
    },
}
PRMTOP_SHA256 = "a61d15bf0bf78675be93275d45f274e808ed6ae450fc1ca21a8e14aee8c12ca0"
QM_CONTRACT = {
    "qm_atoms": 146,
    "qm_charge": 0,
    "electrons": 510,
    "link_atoms": 6,
    "qm_waters": 0,
}
FIRST_WINDOW_DELTAS_A = {
    "attack": -0.04,
    "cn": 0.08,
    "carbonyl": 0.03,
    "nalpha_hg1": 0.05,
    "hg1_n3": -0.05,
}
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
    if local == 0:
        mode = MODES[0]
        scale = None
    elif local <= 4:
        mode = MODES[1]
        scale = FORCE_SCALES[local - 1]
    else:
        mode = MODES[2]
        scale = FORCE_SCALES[local - 5]
    return {
        "task_index": index,
        "seed": seed,
        "seed_label": f"seed{seed}",
        "seed_index": SOURCES[seed]["seed_index"],
        "candidate": SOURCES[seed]["candidate"],
        "mode": mode,
        "mechanism": mode,
        "scale": scale,
        "force_scale": scale,
        "source_basin": "RAW_UNBIASED_NAC",
        "destination": "MATCHED_BASELINE_OR_FIRST_FORWARD_WINDOW",
    }


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "mapping": "2 seeds x (1 matched baseline + 2 mechanisms x 4 force scales)",
        "modes": list(MODES),
        "force_scales": list(FORCE_SCALES),
        "source_authority": "hash-fixed strict-NAC source.gro transplanted onto frozen prmtop",
        "matched_baseline": (
            "same Hamiltonian, positional restraints and minimization; "
            "no reaction-coordinate restraints"
        ),
        "first_window_deltas_A": dict(FIRST_WINDOW_DELTAS_A),
        "qm_contract": dict(QM_CONTRACT),
        "prmtop_sha256": PRMTOP_SHA256,
        "automatic_downstream_action": "NONE",
        "restrained_structures_are_not_ts": True,
    }


def _qpt(geometry: Mapping[str, Any]) -> float:
    return float(geometry["nalpha_hg1_A"]) - float(geometry["hg1_n3_A"])


def stage_spec(spec: Mapping[str, Any], source: Mapping[str, Any]) -> dict[str, Any]:
    mode = str(spec["mode"])
    target = {
        "attack_A": float(source["attack_A"]),
        "c12_n3_A": float(source["c12_n3_A"]),
        "c12_o2_A": float(source["c12_o2_A"]),
        "nalpha_hg1_A": float(source["nalpha_hg1_A"]),
        "hg1_n3_A": float(source["hg1_n3_A"]),
    }
    active: list[str] = []
    scale = spec.get("scale")
    if mode != MODES[0]:
        target["attack_A"] += FIRST_WINDOW_DELTAS_A["attack"]
        target["c12_o2_A"] += FIRST_WINDOW_DELTAS_A["carbonyl"]
        active = ["attack", "carbonyl"]
        if mode == MODES[2]:
            target["c12_n3_A"] += FIRST_WINDOW_DELTAS_A["cn"]
            target["nalpha_hg1_A"] += FIRST_WINDOW_DELTAS_A["nalpha_hg1"]
            target["hg1_n3_A"] += FIRST_WINDOW_DELTAS_A["hg1_n3"]
            active = ["attack", "cn", "carbonyl", "pt"]
    return {
        "window_index": 0,
        "mode": mode,
        "mechanism": mode,
        "force_scale": scale,
        "targets": target,
        "active_coordinates": active,
        "force_bond": FWD.CAL.BASE_FORCES["force_bond"] * float(scale or 0.0),
        "force_carbonyl": FWD.CAL.BASE_FORCES["force_carbonyl"] * float(scale or 0.0),
        "force_pt": FWD.CAL.BASE_FORCES["force_pt"] * float(scale or 0.0),
        "maxcyc": 2200,
        "ncyc": 550,
    }


def reactive_restraints(
    spec: Mapping[str, Any], source_geometry: Mapping[str, Any]
) -> list[str]:
    if spec["mode"] == MODES[0]:
        return []
    stage = stage_spec(spec, source_geometry)
    target = stage["targets"]
    rows = [
        FWD.CAL.BRIDGE._tight_distance_restraint(
            REACTIVE["og1"], REACTIVE["c12"], target["attack_A"], stage["force_bond"]
        ),
        FWD.CAL.BRIDGE._tight_distance_restraint(
            REACTIVE["c12"], REACTIVE["o2"], target["c12_o2_A"],
            stage["force_carbonyl"],
        ),
    ]
    if spec["mode"] == MODES[2]:
        rows.extend((
            FWD.CAL.BRIDGE._tight_distance_restraint(
                REACTIVE["c12"], REACTIVE["n3"], target["c12_n3_A"],
                stage["force_bond"],
            ),
            FWD.CAL.BRIDGE._tight_distance_restraint(
                REACTIVE["nalpha"], REACTIVE["hg1"], target["nalpha_hg1_A"],
                stage["force_pt"],
            ),
            FWD.CAL.BRIDGE._tight_distance_restraint(
                REACTIVE["hg1"], REACTIVE["n3"], target["hg1_n3_A"],
                stage["force_pt"],
            ),
        ))
    return rows


def nmropt_for_spec(spec: Mapping[str, Any]) -> bool:
    return spec["mode"] != MODES[0]


def first_window_guard(geometry: Mapping[str, Any]) -> dict[str, Any]:
    nearest = geometry.get("hg1_nearest_qm_heavy_atom", geometry.get("proton_nearest"))
    values = (
        geometry["attack_A"], geometry.get("c12_n3_A", geometry.get("cn_A")),
        geometry.get("c12_o2_A", geometry.get("carbonyl_A")),
    )
    checks = {
        "finite_geometry": all(math.isfinite(float(value)) for value in values),
        "attack_not_overcompressed": float(geometry["attack_A"]) >= 1.25,
        "cn_not_overcompressed": float(
            geometry.get("c12_n3_A", geometry.get("cn_A"))
        ) >= 1.25,
        "carbonyl_safe": 1.15 <= float(
            geometry.get("c12_o2_A", geometry.get("carbonyl_A"))
        ) <= 1.55,
        "proton_not_detached_or_misrouted": nearest in (
            "Nalpha", "N3", REACTIVE["nalpha"], REACTIVE["n3"]
        ),
    }
    return {"pass": all(checks.values()), "checks": checks}


def _baseline_input(spec: Mapping[str, Any], stage: Mapping[str, Any], qmmask: str) -> str:
    text = FWD.CAL.BRIDGE.tetra_minimization_input(spec, stage, qmmask)
    if text.count("nmropt=1") != 1:
        raise ValueError("base minimization nmropt authority changed")
    text = text.replace("nmropt=1", "nmropt=0")
    kept = [
        line for line in text.splitlines()
        if not line.strip().startswith(("&wt", "DISANG=", "DUMPAVE="))
    ]
    return "\n".join(kept) + "\n"


def validate_full_contract(
    payload: Mapping[str, Any], qmmask: str
) -> dict[str, Any]:
    contract = FWD.CAL._validate_contract(payload)
    heavy = contract.get("qm_heavy_atom_indices")
    if contract.get("qmmask") != qmmask:
        raise ValueError("raw NAC qmmask differs from frozen full contract")
    if not isinstance(heavy, list) or not heavy or not all(
        isinstance(index, int) and index > 0 for index in heavy
    ):
        raise ValueError("full QM geometry contract lacks qm_heavy_atom_indices")
    return contract


def initialize(task_index: int, root: pathlib.Path, commit: str) -> dict[str, Any]:
    if root.exists():
        raise FileExistsError(root)
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("runtime is not bound to a full Git commit")
    spec = task_spec(task_index)
    source = SOURCES[spec["seed"]]
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
    if len(structure.atoms) != RAW.BASE.EXPECTED_SYSTEM_ATOMS or structure.box is None:
        raise ValueError("coordinate transplant changed atom count or box")
    for atom_index, resname, atom_name in (
        (REACTIVE["nalpha"], "THR", "N"),
        (REACTIVE["og1"], "THR", "OG1"),
        (REACTIVE["hg1"], "THR", "HG1"),
        (REACTIVE["c12"], "L2", "C12"),
        (REACTIVE["o2"], "L2", "O2"),
        (REACTIVE["n3"], "L2", "N3"),
    ):
        atom = structure.atoms[atom_index - 1]
        if atom.residue.name != resname or atom.name != atom_name:
            raise ValueError(f"frozen atom identity changed at {atom_index}")

    root.mkdir(parents=True, exist_ok=False)
    start = root / "source.rst7"
    structure.save(str(start), overwrite=False)
    if sha256(start) != source["start_rst7_sha256"]:
        raise ValueError("deterministic raw NAC restart SHA changed")
    contract_manifest = (
        CONTRACT_ROOT / source["contract_attempt"] / "ENDPOINT_MANIFEST.json"
    )
    if not contract_manifest.is_file():
        raise FileNotFoundError(contract_manifest)
    contract = validate_full_contract(BASE.read_json(contract_manifest), qmmask)
    geometry = BASE.AUTH._geometry(start, {"qm_contract": contract})
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_RAW_NAC_FORWARD_CALIBRATION",
        "github_commit": commit,
        "scientific_status": SCIENTIFIC_STATUS,
        "task": spec,
        "source": {
            "source_gro": str(source_gro),
            "source_gro_sha256": source["source_gro_sha256"],
            "restart": str(start),
            "restart_sha256": source["start_rst7_sha256"],
            "coordinate_transplant": "ParMed frozen prmtop plus same-order raw NAC GRO",
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": sha256(PRMTOP),
        "qm_contract": contract,
        "source_geometry": geometry,
        "matched_baseline": spec["mode"] == MODES[0],
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
    spec = manifest["task"]
    stage = stage_spec(spec, manifest["source_geometry"])
    if nmropt_for_spec(spec):
        stage_input = FWD.CAL.BRIDGE.tetra_minimization_input(
            spec, stage, manifest["qm_contract"]["qmmask"]
        )
        rows = reactive_restraints(spec, manifest["source_geometry"])
        (scratch / "restraints.RST").write_text("".join(rows), encoding="utf-8")
    else:
        stage_input = _baseline_input(
            spec, stage, manifest["qm_contract"]["qmmask"]
        )
        rows = []
    (scratch / "stage.in").write_text(stage_input, encoding="utf-8")
    prepared = {
        "schema_version": 1,
        "stage": stage,
        "input_restart": manifest["source"]["restart"],
        "input_restart_sha256": manifest["source"]["restart_sha256"],
        "reaction_coordinate_restraints": len(rows),
        "nmropt": int(nmropt_for_spec(spec)),
        "matched_baseline": spec["mode"] == MODES[0],
        "first_window_only": True,
        "shooting_forbidden": True,
    }
    BASE.write_json(root / "WINDOW_MANIFEST.json", prepared)
    return prepared


def _source_response(
    previous: Mapping[str, Any], current: Mapping[str, Any],
    target: Mapping[str, Any], active: list[str]
) -> dict[str, Any]:
    response = {
        "attack": FWD.CAL.BRIDGE._response(
            previous["attack_A"], current["attack_A"], target["attack_A"]
        ),
        "cn": FWD.CAL.BRIDGE._response(
            previous["c12_n3_A"], current["c12_n3_A"], target["c12_n3_A"]
        ),
        "carbonyl": FWD.CAL.BRIDGE._response(
            previous["c12_o2_A"], current["c12_o2_A"], target["c12_o2_A"]
        ),
        "pt": FWD.CAL.BRIDGE._response(
            _qpt(previous), _qpt(current),
            target["nalpha_hg1_A"] - target["hg1_n3_A"],
        ),
        "active_coordinates": list(active),
    }
    response["all"] = all(response[name]["pass"] for name in active)
    return response


def _write_hashes(root: pathlib.Path) -> None:
    rows = []
    for path in sorted(root.iterdir()):
        if path.is_file() and path.name != "SHA256.tsv":
            rows.append(f"{sha256(path)}  {path.name}")
    (root / "SHA256.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")


def audit(root: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest = BASE.read_json(root / "MANIFEST.json")
    prepared = BASE.read_json(root / "WINDOW_MANIFEST.json")
    stage = prepared["stage"]
    technical, geometry, diagnostics = BASE.AUTH.TETRA._technical(
        scratch / "stage.out", scratch / "stage.rst7", manifest
    )
    previous = manifest["source_geometry"]
    baseline = manifest["task"]["mode"] == MODES[0]
    response: dict[str, Any] = {"all": False, "active_coordinates": []}
    guard: dict[str, Any] = {"pass": False, "checks": {}}
    if technical:
        guard = first_window_guard(geometry)
        if not baseline:
            response = _source_response(
                previous, geometry, stage["targets"], stage["active_coordinates"]
            )
        shutil.copy2(scratch / "stage.rst7", root / "stage.rst7")
    for name in ("stage.in", "restraints.RST", "stage.mdinfo"):
        path = scratch / name
        if path.is_file():
            shutil.copy2(path, root / name)
    stage_out = scratch / "stage.out"
    if stage_out.is_file():
        tail = stage_out.read_text(encoding="utf-8", errors="replace").splitlines()[-200:]
        (root / "ENGINE_TAIL.txt").write_text("\n".join(tail) + "\n", encoding="utf-8")

    qualified = bool(technical and not baseline and response["all"] and guard["pass"])
    if not technical:
        gate = "NOT_EVALUATED_TECHNICAL_A1_RAW_NAC_FORWARD_CALIBRATION"
    elif baseline:
        gate = "PASS_REACTION_COORDINATE_FREE_BASELINE_TECHNICAL"
    elif qualified:
        gate = "PASS_FIRST_FORWARD_WINDOW_ABSOLUTE_RESPONSE"
    else:
        gate = "FAIL_FIRST_FORWARD_WINDOW_ABSOLUTE_RESPONSE_OR_GUARD"
    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_A1_RAW_NAC_FORWARD_CALIBRATION"
            if technical else
            "NOT_EVALUATED_TECHNICAL_A1_RAW_NAC_FORWARD_CALIBRATION"
        ),
        "technical_complete": technical,
        "scientific_status": SCIENTIFIC_STATUS,
        "scientific_gate": gate,
        "task": manifest["task"],
        "source_geometry": previous,
        "target_geometry": stage["targets"],
        "final_geometry": geometry if technical else None,
        "source_qPT_A": _qpt(previous),
        "target_qPT_A": (
            stage["targets"]["nalpha_hg1_A"] - stage["targets"]["hg1_n3_A"]
        ),
        "final_qPT_A": _qpt(geometry) if technical else None,
        "response_from_raw_source": response,
        "first_window_guard": guard,
        "eligible_by_absolute_response": qualified,
        "matched_baseline": baseline,
        "first_window_only": True,
        "automatic_downstream_action": "NONE",
        "diagnostics": diagnostics,
        "restrained_structures_are_not_ts": True,
    }
    BASE.write_json(root / "RESULT.json", result)
    BASE.write_json(root / ("PASS.json" if technical else "NOT_EVALUATED.json"), result)
    _write_hashes(root)
    return result


def _coordinate_value(geometry: Mapping[str, Any], name: str) -> float:
    if name == "attack":
        return float(geometry["attack_A"])
    if name == "cn":
        return float(geometry["c12_n3_A"])
    if name == "carbonyl":
        return float(geometry["c12_o2_A"])
    if name == "pt":
        return _qpt(geometry)
    raise ValueError(name)


def force_effect_vs_baseline(
    item: Mapping[str, Any], baseline: Mapping[str, Any]
) -> dict[str, Any]:
    active = list(item["response_from_raw_source"]["active_coordinates"])
    checks = {}
    for name in active:
        source = _coordinate_value(item["source_geometry"], name)
        target = _coordinate_value(item["target_geometry"], name)
        final = _coordinate_value(item["final_geometry"], name)
        control = _coordinate_value(baseline["final_geometry"], name)
        expected = target - source
        effect = final - control
        checks[name] = {
            "expected_delta_from_raw_A": expected,
            "forced_minus_baseline_A": effect,
            "toward_target_relative_to_baseline": expected * effect > 0.0,
            "closer_to_target_than_baseline": abs(final - target) < abs(control - target),
        }
        checks[name]["pass"] = (
            checks[name]["toward_target_relative_to_baseline"]
            and checks[name]["closer_to_target_than_baseline"]
        )
    return {"all": bool(active) and all(row["pass"] for row in checks.values()),
            "active_coordinates": active, "checks": checks}


def merge_if_ready(output_root: pathlib.Path, array_job: str) -> bool:
    paths = [
        output_root / f"attempt_{array_job}_{index}" / "RESULT.json"
        for index in range(ARRAY_TASKS)
    ]
    if not all(path.is_file() for path in paths):
        return False
    results = [BASE.read_json(path) for path in paths]
    baselines = {
        item["task"]["seed"]: item for item in results if item["matched_baseline"]
    }
    per_task = []
    selected = []
    for item in sorted(results, key=lambda row: row["task"]["task_index"]):
        row = dict(item)
        if not row["matched_baseline"] and row.get("technical_complete") is True:
            effect = force_effect_vs_baseline(row, baselines[row["task"]["seed"]])
            row["force_effect_vs_matched_baseline"] = effect
            row["eligible_for_inherited_chain"] = bool(
                row["eligible_by_absolute_response"] and effect["all"]
            )
        else:
            row["force_effect_vs_matched_baseline"] = None
            row["eligible_for_inherited_chain"] = False
        per_task.append(row)
    for seed in (26723, 26737):
        for mode in MODES[1:]:
            candidates = [
                row for row in per_task
                if row["task"]["seed"] == seed
                and row["task"]["mode"] == mode
                and row["eligible_for_inherited_chain"]
            ]
            selected.append({
                "seed": seed,
                "mode": mode,
                "selected_weakest_scale": (
                    min(row["task"]["scale"] for row in candidates)
                    if candidates else None
                ),
            })
    technical = sum(row.get("technical_complete") is True for row in results)
    baseline_pass = sum(
        row.get("technical_complete") is True and row["matched_baseline"]
        for row in results
    )
    complete_selection = all(row["selected_weakest_scale"] is not None for row in selected)
    payload = {
        "schema_version": 1,
        "status": (
            "NOT_EVALUATED_TECHNICAL_A1_RAW_NAC_FORWARD_CALIBRATION"
            if technical != ARRAY_TASKS or baseline_pass != 2 else
            "PASS_RAW_NAC_FORWARD_FORCE_CALIBRATION_MATRIX"
            if complete_selection else
            "PARTIAL_RAW_NAC_FORWARD_FORCE_CALIBRATION_MATRIX"
        ),
        "denominator_tasks": ARRAY_TASKS,
        "technical_pass_tasks": technical,
        "matched_baseline_pass_seeds": baseline_pass,
        "qualified_force_tasks": sum(row["eligible_for_inherited_chain"] for row in per_task),
        "selected_minimum_force_scales": selected,
        "per_task": per_task,
        "next_action": (
            "STRICT_INHERITED_WINDOWS_FROM_SELECTED_SCALES"
            if complete_selection else
            "HIGHER_FORCE_OR_SMALLER_WINDOW_CALIBRATION"
        ),
        "scientific_status": SCIENTIFIC_STATUS,
        "automatic_downstream_action": "NONE",
    }
    audit = output_root / "audit" / f"nylc_a1_raw_nac_forward_calibration_{array_job}.json"
    audit.parent.mkdir(parents=True, exist_ok=True)
    if audit.exists() and BASE.read_json(audit) != payload:
        raise FileExistsError(audit)
    if not audit.exists():
        BASE.write_json(audit, payload)
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=(
        "describe", "initialize", "prepare", "audit", "merge-if-ready",
    ))
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
