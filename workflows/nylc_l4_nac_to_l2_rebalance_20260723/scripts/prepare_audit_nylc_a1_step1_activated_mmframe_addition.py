#!/usr/bin/env python3
"""Two-seed activated-A1 MM-frame addition-first pilot."""
from __future__ import annotations

import argparse
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
PRE_PATH = HERE / "prepare_audit_nylc_a1_step1_raw_nac_preorganized_calibration.py"


def _load(name: str, path: pathlib.Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PRE = _load("_a1_preorganized", PRE_PATH)
BASE = PRE.BASE
FWD = PRE.FWD
RAW_FORWARD = PRE.RAW_FORWARD
REACTIVE = PRE.REACTIVE
THR267 = dict(PRE.THR267)
TASK_ROOT = PRE.TASK_ROOT
PRMTOP = PRE.PRMTOP
PRMTOP_SHA256 = PRE.PRMTOP_SHA256

ARRAY_TASKS = 2
MPI_RANKS = 8
WINDOWS = 6
ATTACK_DELTA_A = -0.04
CARBONYL_DELTA_A = 0.03
FORCE_ATTACK = 24.0
FORCE_CARBONYL = 30.0
MAXCYC = 2200
NCYC = 550
OUTPUT_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/a1_step1_activated_mmframe_addition"
)
CONTRACT_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/a1_step1_acyl_release_md_continuation"
)
SOURCE_BASE = (
    TASK_ROOT
    / "a1_activated_nac_20260726/equilibration/nac_evt25_time1462ps"
)
SOURCES = {
    26723: {
        "task_index": 0,
        "frame": 240,
        "time_ps": 480.0,
        "tpr": SOURCE_BASE / "seed26723/attempt_61970146_4_61970151/npt300free/run.tpr",
        "xtc": SOURCE_BASE / "seed26723/attempt_61970146_4_61970151/npt300free/run.xtc",
        "tpr_sha256": "c60078a92c2ace51facde4ef64e453f690177fc88b4d6363427f935944fa2e43",
        "xtc_sha256": "1a54f1b5b9f139b746c22d9e0f7e9a4a94eb8154bf2b881888986eedca933d89",
        "contract_attempt": "attempt_62216380_0",
    },
    26737: {
        "task_index": 1,
        "frame": 103,
        "time_ps": 206.0,
        "tpr": SOURCE_BASE / "seed26737/attempt_61970146_5_61970146/npt300free/run.tpr",
        "xtc": SOURCE_BASE / "seed26737/attempt_61970146_5_61970146/npt300free/run.xtc",
        "tpr_sha256": "dbd19a399547319d10630430ed494d33f6271bab0f30cb0de5a466c6af13ba20",
        "xtc_sha256": "fcba14da98b331368061dcd990f2467628ad77b9b7a9c4ce88090f92e0831b05",
        "contract_attempt": "attempt_62216380_1",
    },
}
REQUIRED_A1_BONDS = {
    frozenset((THR267["nalpha"], THR267["h1"])),
    frozenset((THR267["nalpha"], THR267["h2"])),
    frozenset((THR267["nalpha"], THR267["hg1"])),
    frozenset((THR267["nalpha"], THR267["ca"])),
    frozenset((THR267["ca"], THR267["cb"])),
    frozenset((THR267["cb"], THR267["og1"])),
}
FORBIDDEN_A1_BOND = frozenset((THR267["og1"], THR267["hg1"]))
SCIENTIFIC_STATUS = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: pathlib.Path, payload: Mapping[str, Any]) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(temporary, path)


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if index not in (0, 1):
        raise ValueError("task index must be 0 or 1")
    seed = (26723, 26737)[index]
    return {
        "task_index": index,
        "seed": seed,
        "mechanism": "ACTIVATED_A1_ADDITION_FIRST",
        "starting_state": "NALPHA_H3_PLUS_OG1_MINUS",
    }


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "starting_state": "NALPHA_H3_PLUS_OG1_MINUS",
        "authority_stage_count": 1,
        "authority_reaction_coordinate_restrained": False,
        "authority_position_restrained": False,
        "windows_per_task": WINDOWS,
        "reaction_coordinate_restraints": 2,
        "proton_coordinate_restrained": False,
        "attack_angle_restrained": False,
        "strict_serial_restart_inheritance": True,
        "failed_restart_inheritance": False,
        "automatic_downstream_action": "NONE",
        "sources": {
            str(seed): {
                key: str(value) if isinstance(value, pathlib.Path) else value
                for key, value in source.items()
            }
            for seed, source in SOURCES.items()
        },
    }


def authority_restraints() -> str:
    return ""


def stage_spec(window_index: int, source: Mapping[str, Any]) -> dict[str, Any]:
    index = int(window_index)
    if not 1 <= index <= WINDOWS:
        raise ValueError("window index must be 1..6")
    return {
        "window_index": index,
        "attack_target_A": float(source["attack_A"]) + index * ATTACK_DELTA_A,
        "carbonyl_target_A": float(source["c12_o2_A"]) + index * CARBONYL_DELTA_A,
        "force_attack": FORCE_ATTACK,
        "force_carbonyl": FORCE_CARBONYL,
    }


def reaction_restraints(stage: Mapping[str, Any]) -> str:
    return "".join(
        (
            FWD.CAL.BRIDGE._tight_distance_restraint(
                REACTIVE["og1"], REACTIVE["c12"],
                float(stage["attack_target_A"]), FORCE_ATTACK,
            ),
            FWD.CAL.BRIDGE._tight_distance_restraint(
                REACTIVE["c12"], REACTIVE["o2"],
                float(stage["carbonyl_target_A"]), FORCE_CARBONYL,
            ),
        )
    )


def _validate_prmtop_a1_graph(structure: Any) -> dict[str, Any]:
    bond_set = {
        frozenset((bond.atom1.idx + 1, bond.atom2.idx + 1))
        for bond in structure.bonds
    }
    missing = sorted(
        [sorted(pair) for pair in REQUIRED_A1_BONDS if pair not in bond_set]
    )
    return {
        "pass": not missing and FORBIDDEN_A1_BOND not in bond_set,
        "missing_required_bonds": missing,
        "forbidden_og1_hg1_present": FORBIDDEN_A1_BOND in bond_set,
        "required_bonds": sorted([sorted(pair) for pair in REQUIRED_A1_BONDS]),
    }


def _source_geometry_gate(geometry: Mapping[str, Any]) -> dict[str, Any]:
    checks = {
        "strict_nac_attack": float(geometry["attack_A"]) <= 3.5,
        "strict_nac_angle": 95.0 <= float(geometry["attack_angle_deg"]) <= 115.0,
    }
    return {"pass": all(checks.values()), "checks": checks}


def _contract(source: Mapping[str, Any]) -> dict[str, Any]:
    path = CONTRACT_ROOT / str(source["contract_attempt"]) / "ENDPOINT_MANIFEST.json"
    payload = BASE.read_json(path)
    return RAW_FORWARD.validate_full_contract(payload, payload["qm_contract"]["qmmask"])


def initialize(task_index: int, root: pathlib.Path, commit: str) -> dict[str, Any]:
    if root.exists():
        raise FileExistsError(root)
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("runtime is not bound to a full Git commit")
    task = task_spec(task_index)
    source = SOURCES[task["seed"]]
    for key in ("tpr", "xtc"):
        path = pathlib.Path(source[key])
        if not path.is_file() or sha256(path) != source[f"{key}_sha256"]:
            raise ValueError(f"fixed {key} authority failed")
    if sha256(PRMTOP) != PRMTOP_SHA256:
        raise ValueError("frozen prmtop SHA changed")

    import MDAnalysis as mda
    import parmed as pmd

    universe = mda.Universe(str(source["tpr"]), str(source["xtc"]))
    universe.trajectory[int(source["frame"])]
    structure = pmd.load_file(str(PRMTOP))
    if len(universe.atoms) != len(structure.atoms):
        raise ValueError("TPR/XTC atom count differs from frozen Amber topology")
    structure.coordinates = universe.atoms.positions.copy()
    dims = universe.dimensions
    if dims is None or len(dims) != 6:
        raise ValueError("source frame lacks triclinic periodic box")
    structure.box = [float(value) for value in dims]

    root.mkdir(parents=True, exist_ok=False)
    restart = root / "source.rst7"
    structure.save(str(restart), overwrite=False)
    contract = _contract(source)
    geometry = BASE.AUTH._geometry(restart, {"qm_contract": contract})
    topology_gate = _validate_prmtop_a1_graph(structure)
    integrity = PRE.measure_thr267_chemical_integrity(restart)
    geometry_gate = _source_geometry_gate(geometry)
    ready = bool(topology_gate["pass"] and integrity["pass"] and geometry_gate["pass"])
    manifest = {
        "schema_version": 1,
        "status": "READY_ACTIVATED_A1_MMFRAME" if ready else
                  "NOT_EVALUATED_TECHNICAL_SOURCE_INTEGRITY",
        "github_commit": commit,
        "task": task,
        "source": {
            "tpr": str(source["tpr"]),
            "xtc": str(source["xtc"]),
            "tpr_sha256": source["tpr_sha256"],
            "xtc_sha256": source["xtc_sha256"],
            "frame": source["frame"],
            "time_ps": source["time_ps"],
            "restart": str(restart),
            "restart_sha256": sha256(restart),
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": sha256(PRMTOP),
        "qm_contract": contract,
        "source_geometry": geometry,
        "source_topology_gate": topology_gate,
        "source_thr267_integrity": integrity,
        "source_geometry_gate": geometry_gate,
        "ready_for_authority": ready,
        "accepted_restart": None,
        "accepted_restart_sha256": None,
        "accepted_geometry": None,
        "authority": None,
        "windows": [],
        "terminal": not ready,
        "terminal_reason": None if ready else
                           "NOT_EVALUATED_TECHNICAL_SOURCE_INTEGRITY",
        "scientific_status": SCIENTIFIC_STATUS,
        "automatic_downstream_action": "NONE",
    }
    write_json(root / "SOURCE_MANIFEST.json", manifest)
    return manifest


def _qmmm_block(qmmask: str) -> str:
    return f"""&qmmm
  qmmask='{qmmask}',
  qmcharge=0,
  spin=1,
  qm_theory='DFTB3',
  dftb_telec=200.0,
  qmshake=0,
/
"""


def _authority_input(qmmask: str) -> str:
    return f"""NylC activated A1 reaction-coordinate-free authority
&cntrl
  imin=1, ntmin=2, maxcyc={MAXCYC}, ncyc={NCYC}, dx0=0.005, drms=0.01,
  ntb=1, cut=10.0, ntpr=5, ntxo=1,
  ifqnt=1, nmropt=0, ntr=0,
/
{_qmmm_block(qmmask)}"""


def _reaction_input(qmmask: str, index: int) -> str:
    mask = FWD.RAW.NON_QM_SOLUTE_HEAVY_MASK
    return f"""NylC activated A1 addition window {index}
&cntrl
  imin=1, ntmin=2, maxcyc={MAXCYC}, ncyc={NCYC}, dx0=0.005, drms=0.01,
  ntb=1, cut=10.0, ntpr=5, ntxo=1,
  ifqnt=1, nmropt=1, ntr=1,
  restraint_wt=1.0,
  restraintmask='{mask}',
/
{_qmmm_block(qmmask)}&wt type='END' /
DISANG=restraints.RST
DUMPAVE=restraint.dat
"""


def prepare_authority(root: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest = BASE.read_json(root / "SOURCE_MANIFEST.json")
    if not manifest["ready_for_authority"] or manifest["terminal"]:
        raise ValueError("source is not eligible for authority minimization")
    scratch.mkdir(parents=True, exist_ok=False)
    text = _authority_input(manifest["qm_contract"]["qmmask"])
    if "nmropt=0" not in text or "ntr=0" not in text or "DISANG" in text:
        raise ValueError("authority input is not reaction-coordinate free")
    (scratch / "authority.in").write_text(text, encoding="utf-8")
    prepared = {
        "input_restart": manifest["source"]["restart"],
        "input_restart_sha256": manifest["source"]["restart_sha256"],
        "reaction_coordinate_restraints": 0,
        "position_restrained": False,
    }
    write_json(root / "AUTHORITY_MANIFEST.json", prepared)
    return prepared


def _engine_evidence(
    stage_out: pathlib.Path, restart: pathlib.Path, manifest: Mapping[str, Any]
) -> tuple[bool, dict[str, Any], dict[str, Any], dict[str, Any]]:
    technical, geometry, diagnostics = BASE.AUTH.TETRA._technical(
        stage_out, restart, manifest
    )
    text = (
        stage_out.read_text(encoding="utf-8", errors="replace")
        if stage_out.is_file() else ""
    )
    numerical = RAW_FORWARD.engine_numerical_health(text)
    return bool(technical and numerical["pass"]), geometry, diagnostics, numerical


def _copy_engine_evidence(
    target: pathlib.Path, scratch: pathlib.Path, prefix: str
) -> None:
    for name in (f"{prefix}.in", f"{prefix}.out", f"{prefix}.mdinfo",
                 f"{prefix}.rst7", "restraints.RST", "restraint.dat"):
        source = scratch / name
        if source.is_file():
            shutil.copy2(source, target / name)
    stage_out = scratch / f"{prefix}.out"
    if stage_out.is_file():
        tail = stage_out.read_text(
            encoding="utf-8", errors="replace"
        ).splitlines()[-240:]
        (target / "ENGINE_TAIL.txt").write_text(
            "\n".join(tail) + "\n", encoding="utf-8"
        )


def audit_authority(root: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest = BASE.read_json(root / "SOURCE_MANIFEST.json")
    restart = scratch / "authority.rst7"
    engine_ok, geometry, diagnostics, numerical = _engine_evidence(
        scratch / "authority.out", restart, manifest
    )
    integrity = {"pass": False}
    geometry_gate = {"pass": False}
    if engine_ok:
        integrity = PRE.measure_thr267_chemical_integrity(restart)
        geometry_gate = _source_geometry_gate(geometry)
    output_sha = sha256(restart) if restart.is_file() else None
    new_restart = output_sha not in (None, manifest["source"]["restart_sha256"])
    passed = bool(
        engine_ok and integrity["pass"] and geometry_gate["pass"] and new_restart
    )
    result = {
        "schema_version": 1,
        "status": "PASS_ACTIVATED_A1_AUTHORITY" if passed else
                  "NOT_EVALUATED_TECHNICAL_ACTIVATED_A1_AUTHORITY",
        "technical_complete": engine_ok,
        "numerical_health": numerical,
        "thr267_chemical_integrity": integrity,
        "geometry_gate": geometry_gate,
        "source_geometry": manifest["source_geometry"],
        "final_geometry": geometry if engine_ok else None,
        "output_restart_sha256": output_sha,
        "new_restart_sha": new_restart,
        "eligible_for_windows": passed,
        "diagnostics": diagnostics,
    }
    _copy_engine_evidence(root, scratch, "authority")
    if passed:
        persistent = root / "authority.rst7"
        if not persistent.is_file():
            shutil.copy2(restart, persistent)
        manifest["accepted_restart"] = str(persistent)
        manifest["accepted_restart_sha256"] = sha256(persistent)
        manifest["accepted_geometry"] = geometry
    else:
        manifest["terminal"] = True
        manifest["terminal_reason"] = result["status"]
    manifest["authority"] = result
    write_json(root / "AUTHORITY_RESULT.json", result)
    write_json(root / "SOURCE_MANIFEST.json", manifest)
    return result


def prepare_window(
    root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path, index: int
) -> dict[str, Any]:
    manifest = BASE.read_json(root / "SOURCE_MANIFEST.json")
    if manifest["terminal"] or not manifest["authority"]["eligible_for_windows"]:
        raise ValueError("chain is terminal before requested window")
    if int(index) != len(manifest["windows"]) + 1:
        raise ValueError("window order is not strict")
    input_restart = pathlib.Path(manifest["accepted_restart"])
    if sha256(input_restart) != manifest["accepted_restart_sha256"]:
        raise ValueError("accepted restart SHA inheritance failed")
    output.mkdir(parents=True, exist_ok=False)
    scratch.mkdir(parents=True, exist_ok=False)
    stage = stage_spec(index, manifest["authority"]["final_geometry"])
    (scratch / "stage.in").write_text(
        _reaction_input(manifest["qm_contract"]["qmmask"], index),
        encoding="utf-8",
    )
    restraints = reaction_restraints(stage)
    if restraints.count("&rst") != 2 or r"\n" in restraints:
        raise ValueError("reaction restraint authority failed")
    (scratch / "restraints.RST").write_text(restraints, encoding="utf-8")
    prepared = {
        "window_index": index,
        "stage": stage,
        "input_restart": str(input_restart),
        "input_restart_sha256": manifest["accepted_restart_sha256"],
        "strict_inheritance": True,
    }
    write_json(output / "WINDOW_MANIFEST.json", prepared)
    return prepared


def audit_window(
    root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path
) -> dict[str, Any]:
    manifest = BASE.read_json(root / "SOURCE_MANIFEST.json")
    prepared = BASE.read_json(output / "WINDOW_MANIFEST.json")
    restart = scratch / "stage.rst7"
    engine_ok, geometry, diagnostics, numerical = _engine_evidence(
        scratch / "stage.out", restart, manifest
    )
    integrity = {"pass": False}
    if engine_ok:
        integrity = PRE.measure_thr267_chemical_integrity(restart)
    angle_ok = bool(
        engine_ok and 95.0 <= float(geometry["attack_angle_deg"]) <= 115.0
    )
    previous = manifest["accepted_geometry"]
    stage = prepared["stage"]
    response = {"all": False}
    if engine_ok:
        response = {
            "attack": FWD.CAL.BRIDGE._response(
                previous["attack_A"], geometry["attack_A"],
                stage["attack_target_A"],
            ),
            "carbonyl": FWD.CAL.BRIDGE._response(
                previous["c12_o2_A"], geometry["c12_o2_A"],
                stage["carbonyl_target_A"],
            ),
        }
        response["all"] = response["attack"]["pass"] and response["carbonyl"]["pass"]
    accepted = bool(engine_ok and integrity["pass"] and angle_ok and response["all"])
    _copy_engine_evidence(output, scratch, "stage")
    output_sha = sha256(restart) if restart.is_file() else None
    result = {
        "schema_version": 1,
        "window_index": prepared["window_index"],
        "status": "PASS_ACTIVATED_A1_ADDITION_WINDOW" if accepted else
                  "FAIL_ACTIVATED_A1_ADDITION_WINDOW",
        "technical_complete": engine_ok,
        "numerical_health": numerical,
        "thr267_chemical_integrity": integrity,
        "attack_angle_guard": angle_ok,
        "input_restart_sha256": prepared["input_restart_sha256"],
        "output_restart_sha256": output_sha,
        "source_geometry": previous,
        "target_geometry": stage,
        "final_geometry": geometry if engine_ok else None,
        "response": response,
        "accepted": accepted,
        "diagnostics": diagnostics,
    }
    if accepted:
        persistent = output / "stage.rst7"
        if not persistent.is_file():
            shutil.copy2(restart, persistent)
        manifest["accepted_restart"] = str(persistent)
        manifest["accepted_restart_sha256"] = sha256(persistent)
        manifest["accepted_geometry"] = geometry
    else:
        manifest["terminal"] = True
        manifest["terminal_reason"] = result["status"]
    manifest["windows"].append({
        "window_index": prepared["window_index"],
        "result": str(output / "RESULT.json"),
        "accepted": accepted,
        "input_restart_sha256": prepared["input_restart_sha256"],
        "output_restart_sha256": output_sha,
    })
    write_json(output / "RESULT.json", result)
    write_json(root / "SOURCE_MANIFEST.json", manifest)
    return result


def _write_hashes(root: pathlib.Path) -> None:
    rows = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.name != "SHA256.tsv":
            rows.append(f"{sha256(path)}  {path.relative_to(root)}")
    (root / "SHA256.tsv").write_text(
        "\n".join(rows) + "\n", encoding="utf-8"
    )


def finalize(root: pathlib.Path) -> dict[str, Any]:
    manifest = BASE.read_json(root / "SOURCE_MANIFEST.json")
    accepted = sum(row["accepted"] for row in manifest["windows"])
    authority_pass = bool(
        manifest.get("authority")
        and manifest["authority"].get("eligible_for_windows")
    )
    if not authority_pass:
        status = manifest.get("terminal_reason") or "NOT_EVALUATED_TECHNICAL"
    elif accepted:
        status = "ACTIVATED_A1_EARLY_TETRAHEDRAL_RESPONSE"
    else:
        status = "NO_TETRAHEDRAL_RESPONSE_UNDER_TESTED_BIAS"
    result = {
        "schema_version": 1,
        "status": status,
        "technical_authority_pass": authority_pass,
        "attempted_windows": len(manifest["windows"]),
        "accepted_windows": accepted,
        "terminal_reason": manifest.get("terminal_reason"),
        "terminal_geometry": manifest.get("accepted_geometry"),
        "strict_inheritance": True,
        "scientific_status": SCIENTIFIC_STATUS,
        "automatic_downstream_action": "NONE",
    }
    write_json(root / "RESULT.json", result)
    marker = "PASS.json" if authority_pass else "NOT_EVALUATED.json"
    write_json(root / marker, result)
    _write_hashes(root)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=(
        "describe", "initialize", "prepare-authority", "audit-authority",
        "prepare-window", "audit-window", "finalize",
    ))
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--root", type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--scratch", type=pathlib.Path)
    parser.add_argument("--window-index", type=int)
    parser.add_argument("--github-commit", default="unknown")
    args = parser.parse_args()
    if args.mode == "describe":
        print(json.dumps(describe(), sort_keys=True))
    elif args.mode == "initialize":
        initialize(args.task_index, args.root.resolve(), args.github_commit)
    elif args.mode == "prepare-authority":
        prepare_authority(args.root.resolve(), args.scratch.resolve())
    elif args.mode == "audit-authority":
        audit_authority(args.root.resolve(), args.scratch.resolve())
    elif args.mode == "prepare-window":
        prepare_window(
            args.root.resolve(), args.output.resolve(), args.scratch.resolve(),
            args.window_index,
        )
    elif args.mode == "audit-window":
        audit_window(args.root.resolve(), args.output.resolve(), args.scratch.resolve())
    else:
        finalize(args.root.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
