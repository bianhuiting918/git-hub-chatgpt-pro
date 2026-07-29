#!/usr/bin/env python3
"""Bounded NylC A1 Step2 water-recruitment A2 scout.

One deterministic nearest complete H2O is guided into an attack pose with four
independent distance restraints.  All guide restraints are removed before two
independent final-Hamiltonian A2 release legs.  This is not a product, TS, path,
PMF, barrier, or mechanism calculation.
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
from typing import Any, Mapping, Sequence

HERE = pathlib.Path(__file__).resolve().parent
S2_PATH = HERE / "prepare_audit_nylc_a1_step2_qmwater_endpoint.py"
_SPEC = importlib.util.spec_from_file_location("_nylc_step2_base", S2_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"cannot load {S2_PATH}")
S2 = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(S2)
AC = S2.ACYL

TASK_ROOT = S2.TASK_ROOT
PRMTOP = S2.PRMTOP
SOURCES = S2.SOURCES
REACTIVE = S2.REACTIVE
EXPECTED_CONTRACT = dict(S2.EXPECTED_CONTRACT)
EXPECTED_NEAREST_WATER = {
    "seed26723": (13046, (13047, 13048)),
    "seed26737": (13046, (13047, 13048)),
}
GUIDED_TARGETS_A = {
    "c12_ow": 3.00,
    "o2_ow": 3.60,
    "donor_h_nalpha": 2.00,
    "ow_nalpha": 2.87,
}
GUIDED_FORCE_KCAL_MOL_A2 = 5.0
GUIDED_MAXCYC = 600
GUIDED_NCYC = 200
ARRAY_TASKS = 4
SCOPE = "STEP2_BOUNDED_WATER_RECRUITMENT_A2_ONLY_NOT_PRODUCT_TS_PATH_PMF_BARRIER_OR_MECHANISM"
NEXT = "STEP2_PRODUCT_ENDPOINT_BLOCKED_PENDING_RECRUITED_A2_PASS"


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: pathlib.Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def task_spec(task_index: int) -> dict[str, Any]:
    if not 0 <= int(task_index) < ARRAY_TASKS:
        raise ValueError("task index must be 0..3")
    seed_index = int(task_index) // 2
    donor_branch = int(task_index) % 2
    source = S2.source_for_index(seed_index)
    return {
        "task_index": int(task_index),
        "seed_index": seed_index,
        "donor_branch": donor_branch,
        "source": source,
    }


def guided_restraints(oxygen_index1: int, donor_h_index1: int) -> str:
    return "".join(
        (
            AC.distance_restraint(REACTIVE["c12"], oxygen_index1, GUIDED_TARGETS_A["c12_ow"], GUIDED_FORCE_KCAL_MOL_A2),
            AC.distance_restraint(REACTIVE["o2"], oxygen_index1, GUIDED_TARGETS_A["o2_ow"], GUIDED_FORCE_KCAL_MOL_A2),
            AC.distance_restraint(donor_h_index1, REACTIVE["nalpha"], GUIDED_TARGETS_A["donor_h_nalpha"], GUIDED_FORCE_KCAL_MOL_A2),
            AC.distance_restraint(oxygen_index1, REACTIVE["nalpha"], GUIDED_TARGETS_A["ow_nalpha"], GUIDED_FORCE_KCAL_MOL_A2),
        )
    )


def guided_input(seed: str, branch: int, qmmask: str) -> str:
    return f"""NylC A1 Step2 water recruitment {seed} donor branch {branch}
&cntrl
  imin=1, ntmin=2, maxcyc={GUIDED_MAXCYC}, ncyc={GUIDED_NCYC}, dx0=0.005,
  ntb=1, cut=10.0, ntpr=5, ntxo=1, ifqnt=1,
  ntr=1, nmropt=1, restraint_wt=1.0,
  restraintmask='{AC.BASE.NON_QM_SOLUTE_HEAVY_MASK}', drms=0.10,
/
{S2.qmmm_block(qmmask)}&wt type='END' /
DISANG=restraints.RST
DUMPAVE=restraint.dat
"""


def _nearest_complete_water(structure: Any) -> tuple[Any, tuple[Any, Any], float, int]:
    cell = S2._periodic_cell(structure.box)
    c12 = structure.atoms[REACTIVE["c12"] - 1]
    rows = []
    complete = 0
    for residue in structure.residues:
        parsed = S2._complete_h2o(residue)
        if parsed is None:
            continue
        complete += 1
        oxygen, hydrogens = parsed
        rows.append((S2._distance(c12, oxygen, cell), residue.idx, oxygen, hydrogens))
    if not rows:
        raise ValueError("no complete water exists in the full-system restart")
    rows.sort(key=lambda item: (item[0], item[1], item[2].idx))
    distance, _, oxygen, hydrogens = rows[0]
    return oxygen, hydrogens, float(distance), complete


def _derive_qm_contract(structure: Any, source_manifest: Mapping[str, Any], water_indices: Sequence[int]) -> tuple[list[int], dict[str, Any], str]:
    base_indices = S2._parse_qmmask(source_manifest["qm_contract"]["qmmask"], 146)
    if set(base_indices).intersection(water_indices):
        raise ValueError("nearest water overlaps the frozen Step1 QM mask")
    qm_indices = base_indices + [int(value) for value in water_indices]
    if len(qm_indices) != 149 or len(set(qm_indices)) != 149:
        raise ValueError("Step2 QM mask is not 149 unique atoms")
    boundaries = S2._boundary_bonds(structure, qm_indices)
    derived = {
        "qm_atom_count": len(qm_indices),
        "qmcharge": 0,
        "spin": 1,
        "link_atom_count": len(boundaries),
        "electron_count": 518,
    }
    if any(derived[key] != value for key, value in EXPECTED_CONTRACT.items()):
        raise ValueError(f"derived Step2 contract {derived} != {EXPECTED_CONTRACT}")
    qmmask = ",".join(f"@{index}" for index in qm_indices)
    contract = dict(source_manifest["qm_contract"])
    contract.update(
        {
            "expected": EXPECTED_CONTRACT,
            "derived": derived,
            "boundary_bonds": boundaries,
            "qmmask": qmmask,
            "base_qm_atom_count": 146,
            "selected_complete_water_atom_indices1": list(water_indices),
            "step2_qm_water_count": 1,
        }
    )
    return qm_indices, contract, qmmask


def initialize(task_index: int, output: pathlib.Path, scratch: pathlib.Path, code_root: pathlib.Path, github_commit: str) -> None:
    if output.exists() or scratch.exists():
        raise FileExistsError("output or scratch already exists")
    if not re.fullmatch(r"[0-9a-f]{40}", github_commit):
        raise ValueError("runtime code is not bound to a full Git commit")
    task = task_spec(task_index)
    source = task["source"]
    source_manifest, _ = S2.validate_authority(source, code_root)
    import parmed as pmd

    structure = pmd.load_file(str(PRMTOP), xyz=str(source["restart"]))
    AC._validate_structure(structure)
    oxygen, hydrogens, nearest_distance, complete_count = _nearest_complete_water(structure)
    observed = (oxygen.idx + 1, tuple(atom.idx + 1 for atom in hydrogens))
    if observed != EXPECTED_NEAREST_WATER[source["seed"]]:
        raise ValueError(f"nearest-water identity changed: {observed}")
    donor = hydrogens[task["donor_branch"]]
    water_indices = [oxygen.idx + 1, *(atom.idx + 1 for atom in hydrogens)]
    _, contract, qmmask = _derive_qm_contract(structure, source_manifest, water_indices)

    output.mkdir(parents=True, exist_ok=False)
    guided = scratch / "guided"
    guided.mkdir(parents=True, exist_ok=False)
    (guided / "stage.in").write_text(
        guided_input(source["seed"], task["donor_branch"], qmmask), encoding="utf-8"
    )
    (guided / "restraints.RST").write_text(
        guided_restraints(oxygen.idx + 1, donor.idx + 1), encoding="utf-8"
    )
    selected = {
        "residue_index1": oxygen.residue.idx + 1,
        "residue_name": oxygen.residue.name,
        "oxygen_index1": oxygen.idx + 1,
        "hydrogen_indices1": [atom.idx + 1 for atom in hydrogens],
        "donor_h_index1": donor.idx + 1,
        "initial_c12_ow_A": nearest_distance,
        "complete_water_denominator": complete_count,
        "selection_policy": "nearest_complete_H2O_by_triclinic_minimum_image_then_two_donor_H_branches",
    }
    manifest = {
        "schema_version": 1,
        "status": "READY_STEP2_WATER_RECRUITMENT_GUIDED_MIN",
        "scientific_scope": SCOPE,
        "github_commit": github_commit,
        "task_index": task["task_index"],
        "seed_index": task["seed_index"],
        "donor_branch": task["donor_branch"],
        "seed": source["seed"],
        "source": {
            "restart": str(source["restart"]),
            "restart_sha256": source["restart_sha256"],
            "attempt": source["attempt"],
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": S2.EXPECTED_PRMTOP_SHA256,
        "selected_water": selected,
        "qm_contract": contract,
        "box_lengths_A": [float(value) for value in structure.box[:6]],
        "recruitment": {
            "guided_targets_A": GUIDED_TARGETS_A,
            "guided_force_kcal_mol_A2": GUIDED_FORCE_KCAL_MOL_A2,
            "restraint_representation": "four_independent_distance_restraints",
            "all_guided_restraints_removed_for_release": True,
            "guided_result": None,
        },
        "protocol": {
            "guided_minimization_maxcyc": GUIDED_MAXCYC,
            "release_legs": 2,
            "release_execution": "PARALLEL_TWO_INDEPENDENT_8_RANK_LEGS",
            "release_reactive_restraints": False,
            "release_position_restraints": False,
            "release_trajectory_policy": "scratch_only",
        },
        "NEXT": NEXT,
        "product_endpoint_implemented": False,
    }
    write_json(output / "A2_MANIFEST.json", manifest)
    write_json(
        output / "READY_GUIDED.json",
        {
            "status": "READY_STEP2_WATER_RECRUITMENT_GUIDED_MIN",
            "seed": source["seed"],
            "task_index": task["task_index"],
            "selected_water": selected,
            "expected_contract": EXPECTED_CONTRACT,
            "NEXT": NEXT,
        },
    )


def _guided_geometry(restart: pathlib.Path, manifest: Mapping[str, Any]) -> dict[str, Any]:
    import parmed as pmd

    structure = pmd.load_file(str(PRMTOP), xyz=str(restart))
    cell = S2._periodic_cell(structure.box)
    water = manifest["selected_water"]
    c12 = structure.atoms[REACTIVE["c12"] - 1]
    o2 = structure.atoms[REACTIVE["o2"] - 1]
    nalpha = structure.atoms[REACTIVE["nalpha"] - 1]
    ow = structure.atoms[int(water["oxygen_index1"]) - 1]
    donor = structure.atoms[int(water["donor_h_index1"]) - 1]
    acyl = AC._read_structure_geometry(restart, manifest)
    result = dict(acyl)
    result.update(
        {
            "water_c12_ow_A": S2._distance(c12, ow, cell),
            "water_attack_angle_deg": S2._angle(o2, c12, ow, cell),
            "water_h_nalpha_A": S2._distance(donor, nalpha, cell),
            "water_hbond_angle_deg": S2._angle(ow, donor, nalpha, cell),
        }
    )
    result["acyl_sensitivity_like"] = S2._acyl_sensitivity_like(result)
    result["guided_attack_pose"] = bool(
        2.40 <= result["water_c12_ow_A"] <= 3.50
        and 95.0 <= result["water_attack_angle_deg"] <= 125.0
        and result["water_h_nalpha_A"] <= 2.40
        and result["water_hbond_angle_deg"] >= 135.0
        and result["acyl_sensitivity_like"]
    )
    return result


def audit_guided(output: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest_path = output / "A2_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    stage = scratch / "guided"
    text = (stage / "stage.out").read_text(encoding="utf-8", errors="replace") if (stage / "stage.out").is_file() else ""
    engine_rc = int((stage / "engine.rc").read_text().strip()) if (stage / "engine.rc").is_file() else 999
    restart = stage / "stage.rst7"
    banner = S2.parse_banner(text)
    hard = {name: len(re.findall(pattern, text, re.I)) for name, pattern in S2.HARD_PATTERNS.items()}
    geometry: dict[str, Any] = {}
    geometry_error = None
    if restart.is_file() and restart.stat().st_size:
        try:
            geometry = _guided_geometry(restart, manifest)
        except Exception as error:
            geometry_error = f"{type(error).__name__}: {error}"
    technical = bool(
        engine_rc == 0
        and AC.stage_output_complete("guided_min", text)
        and restart.is_file()
        and restart.stat().st_size > 0
        and S2.banner_pass(banner)
        and sum(hard.values()) == 0
        and geometry
        and geometry_error is None
    )
    guided_pass = bool(technical and geometry.get("guided_attack_pose"))
    result = {
        "schema_version": 1,
        "status": "PASS_TECHNICAL_STEP2_WATER_RECRUITMENT_GUIDED_MIN" if technical else "NOT_EVALUATED_TECHNICAL_STEP2_WATER_RECRUITMENT_GUIDED_MIN",
        "technical_pass": technical,
        "guided_attack_pose_pass": guided_pass,
        "classification": "PASS_GUIDED_ATTACK_POSE" if guided_pass else "FAIL_GUIDED_ATTACK_POSE" if technical else "NOT_EVALUATED_GUIDED_TECHNICAL",
        "engine_exit_code": engine_rc,
        "banner_contract_pass": S2.banner_pass(banner),
        "banner_observed": banner,
        "geometry": geometry,
        "restart_path": str(restart),
        "restart_sha256": sha256(restart) if restart.is_file() and restart.stat().st_size else "",
        "diagnostics": {"hard_error_hits": hard, "geometry_error": geometry_error},
        "NEXT": NEXT,
    }
    write_json(output / "GUIDED_RESULT.json", result)
    manifest["recruitment"]["guided_result"] = result
    write_json(manifest_path, manifest)
    if guided_pass:
        source = S2.source_for_index(int(manifest["seed_index"]))
        for leg, velocity_seed in enumerate(source["velocity_seeds"]):
            leg_root = scratch / f"a2_leg{leg}"
            leg_root.mkdir()
            (leg_root / "stage.in").write_text(
                S2.md_input(manifest["seed"], leg, velocity_seed, manifest["qm_contract"]["qmmask"]),
                encoding="utf-8",
            )
            write_json(
                leg_root / "PREPARED.json",
                {
                    "seed": manifest["seed"],
                    "leg": leg,
                    "velocity_seed": velocity_seed,
                    "input_restart": str(restart),
                    "input_restart_sha256": result["restart_sha256"],
                    "expected_contract": EXPECTED_CONTRACT,
                    "reactive_restraints": False,
                },
            )
    return result


def finalize(output: pathlib.Path) -> dict[str, Any]:
    manifest = json.loads((output / "A2_MANIFEST.json").read_text(encoding="utf-8"))
    guided = json.loads((output / "GUIDED_RESULT.json").read_text(encoding="utf-8"))
    if not guided.get("guided_attack_pose_pass"):
        status = "FAIL_STEP2_WATER_RECRUITMENT_GUIDED_POSE" if guided.get("technical_pass") else "NOT_EVALUATED_TECHNICAL_STEP2_WATER_RECRUITMENT"
        result = {
            "schema_version": 1,
            "status": status,
            "classification": status,
            "scientific_scope": SCOPE,
            "task_index": manifest["task_index"],
            "seed_index": manifest["seed_index"],
            "donor_branch": manifest["donor_branch"],
            "seed": manifest["seed"],
            "selected_water": manifest["selected_water"],
            "guided_result": guided,
            "technical_complete": bool(guided.get("technical_pass")),
            "candidate": None,
            "NEXT": NEXT,
            "automatic_downstream_action": "NONE",
        }
    else:
        S2.finalize(output)
        result = json.loads((output / "RESULT.json").read_text(encoding="utf-8"))
        result.update(
            {
                "scientific_scope": SCOPE,
                "task_index": manifest["task_index"],
                "seed_index": manifest["seed_index"],
                "donor_branch": manifest["donor_branch"],
                "guided_result": guided,
                "technical_complete": all(item.get("technical_pass") is True for item in result.get("per_leg", [])),
                "candidate": result if result.get("status") == "PASS_EXPLORATORY_STEP2_A2_ACYL_BASIN_REVALIDATED" else None,
                "NEXT": NEXT,
                "automatic_downstream_action": "NONE",
            }
        )
    for name in ("PASS.json", "FAIL.json", "NOT_EVALUATED.json"):
        path = output / name
        if path.exists():
            path.unlink()
    write_json(output / "RESULT.json", result)
    marker = "PASS.json" if str(result["status"]).startswith("PASS") else "NOT_EVALUATED.json" if str(result["status"]).startswith("NOT_EVALUATED") else "FAIL.json"
    write_json(output / marker, result)
    _write_hashes(output)
    return result


def _write_hashes(output: pathlib.Path) -> None:
    names = (
        "A2_MANIFEST.json", "READY_GUIDED.json", "GUIDED_RESULT.json",
        "A2_LEG_0.json", "A2_LEG_1.json", "RESULT.json", "PASS.json",
        "FAIL.json", "NOT_EVALUATED.json", "a2_leg0_endpoint.rst7",
        "a2_leg1_endpoint.rst7",
    )
    rows = [f"{sha256(output / name)}  {name}" for name in names if (output / name).is_file()]
    (output / "SHA256.tsv").write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")


def merge_if_ready(output_root: pathlib.Path, array_job: str, seed_index: int) -> bool:
    audit = output_root / "audit" / f"nylc_a1_step2_water_recruitment_{array_job}_seed{seed_index}.json"
    audit.parent.mkdir(parents=True, exist_ok=True)
    with (audit.with_suffix(".json.lock")).open("a+") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        try:
            paths = [output_root / f"attempt_{array_job}_{seed_index * 2 + branch}" / "RESULT.json" for branch in (0, 1)]
            if not all(path.is_file() for path in paths):
                return False
            results = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
            technical = all(item.get("technical_complete") is True for item in results)
            passing = [item for item in sorted(results, key=lambda item: item["donor_branch"]) if str(item.get("status", "")).startswith("PASS")]
            selected = passing[0] if technical and passing else None
            payload = {
                "schema_version": 1,
                "status": "NOT_EVALUATED_TECHNICAL_STEP2_WATER_RECRUITMENT" if not technical else "PASS_RECRUITED_A2_WATER_SEED" if selected else "FAIL_NO_STABLE_RECRUITED_A2_WATER",
                "seed_index": seed_index,
                "denominator_donor_branches": 2,
                "per_branch": sorted(results, key=lambda item: item["donor_branch"]),
                "selected_candidate": selected,
                "NEXT": NEXT,
                "automatic_downstream_action": "NONE",
            }
            if audit.exists():
                if json.loads(audit.read_text(encoding="utf-8")) != payload:
                    raise FileExistsError(f"existing audit differs: {audit}")
            else:
                write_json(audit, payload)
            return True
        finally:
            fcntl.flock(handle.fileno(), fcntl.LOCK_UN)


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "seed_count": 2,
        "donor_branches_per_seed": 2,
        "array_task_count": ARRAY_TASKS,
        "nearest_water": EXPECTED_NEAREST_WATER,
        "guided_targets_A": GUIDED_TARGETS_A,
        "guided_force_kcal_mol_A2": GUIDED_FORCE_KCAL_MOL_A2,
        "expected_contract": EXPECTED_CONTRACT,
        "release_legs_per_task": 2,
        "guided_restraints_removed_for_release": True,
        "automatic_downstream_action": "NONE",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=("describe", "initialize", "audit-guided", "audit-leg", "finalize", "merge-if-ready"))
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--seed-index", type=int)
    parser.add_argument("--leg", type=int, choices=(0, 1))
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--output-root", type=pathlib.Path)
    parser.add_argument("--scratch", type=pathlib.Path)
    parser.add_argument("--code-root", type=pathlib.Path)
    parser.add_argument("--github-commit", default="unknown")
    parser.add_argument("--array-job")
    args = parser.parse_args()
    if args.mode == "describe":
        print(json.dumps(describe(), sort_keys=True))
    elif args.mode == "initialize":
        initialize(args.task_index, args.output.resolve(), args.scratch.resolve(), args.code_root.resolve(), args.github_commit)
    elif args.mode == "audit-guided":
        audit_guided(args.output.resolve(), args.scratch.resolve())
    elif args.mode == "audit-leg":
        S2.audit_leg(args.output.resolve(), args.scratch.resolve(), args.leg)
    elif args.mode == "finalize":
        finalize(args.output.resolve())
    else:
        merge_if_ready(args.output_root.resolve(), args.array_job, args.seed_index)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
