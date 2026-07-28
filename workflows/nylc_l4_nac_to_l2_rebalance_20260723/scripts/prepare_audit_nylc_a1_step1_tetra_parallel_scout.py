#!/usr/bin/env python3
"""Parallel NylC A1 tetrahedralization scouts; restrained geometry probes, not paths."""

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

_AC_PATH = pathlib.Path(__file__).with_name("prepare_audit_nylc_a1_step1_acyl_endpoint_stability.py")
_spec = importlib.util.spec_from_file_location("_a1_acyl_base", _AC_PATH)
if _spec is None or _spec.loader is None:
    raise RuntimeError(f"cannot import immutable base {_AC_PATH}")
AC = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(AC)
BASE = AC.BASE

A_WINDOWS = (
    {"attack": 2.20, "co": 1.27, "k_attack": 25.0, "k_co": 10.0, "maxcyc": 400, "ncyc": 100},
    {"attack": 1.95, "co": 1.31, "k_attack": 40.0, "k_co": 15.0, "maxcyc": 600, "ncyc": 150},
    {"attack": 1.75, "co": 1.36, "k_attack": 60.0, "k_co": 20.0, "maxcyc": 600, "ncyc": 150},
    {"attack": 1.58, "co": 1.42, "k_attack": 80.0, "k_co": 25.0, "maxcyc": 600, "ncyc": 150},
)
B_WINDOWS = (
    {"attack": 1.85, "k_attack": 40.0, "k_torsion": 5.0, "maxcyc": 600, "ncyc": 150},
    {"attack": 1.65, "k_attack": 60.0, "k_torsion": 10.0, "maxcyc": 800, "ncyc": 200},
)
SCIENTIFIC_BOUNDARY = "NOT_EVALUATED_TS_PMF_BARRIER_MECHANISM"
NEXT_ACTION = "DO_NOT_START_PMF_NEB_STRING"
SCOPE = "bounded_tetrahedralization_geometry_scout_not_path_or_barrier"


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if index in (0, 1):
        return {"task_index": index, "branch": "A_ATTACK_CO", "seed_index": index,
                "sign": None, "windows": A_WINDOWS}
    if 2 <= index <= 5:
        offset = index - 2
        return {"task_index": index, "branch": "B_ATTACK_OOP", "seed_index": offset % 2,
                "sign": "plus" if offset < 2 else "minus", "windows": B_WINDOWS}
    raise ValueError(f"task index {index} is outside 0..5")


def _vec_sub(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(x - y for x, y in zip(a, b))


def _vec_add(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return tuple(x + y for x, y in zip(a, b))


def _scale(a: tuple[float, float, float], value: float) -> tuple[float, float, float]:
    return tuple(value * x for x in a)


def _dot(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return sum(x * y for x, y in zip(a, b))


def _cross(a: tuple[float, float, float], b: tuple[float, float, float]) -> tuple[float, float, float]:
    return (a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0])


def _unit(a: tuple[float, float, float]) -> tuple[float, float, float]:
    norm = math.sqrt(_dot(a, a))
    if norm == 0:
        raise ValueError("degenerate geometry")
    return _scale(a, 1.0 / norm)


def _dihedral(a: tuple[float, float, float], b: tuple[float, float, float],
              c: tuple[float, float, float], d: tuple[float, float, float]) -> float:
    b0 = _vec_sub(a, b)
    b1 = _vec_sub(c, b)
    b2 = _vec_sub(d, c)
    axis = _unit(b1)
    v = _vec_sub(b0, _scale(axis, _dot(b0, axis)))
    w = _vec_sub(b2, _scale(axis, _dot(b2, axis)))
    return math.degrees(math.atan2(_dot(_cross(axis, v), w), _dot(v, w)))


def _xyz(structure: Any, atom_index: int) -> tuple[float, float, float]:
    atom = structure.atoms[atom_index - 1]
    return float(atom.xx), float(atom.xy), float(atom.xz)


def torsion_target(structure: Any, sign: str) -> float:
    c12 = _xyz(structure, AC.REACTIVE_ATOMS["c12"])
    o2 = _xyz(structure, AC.REACTIVE_ATOMS["o2"])
    n3 = _xyz(structure, AC.REACTIVE_ATOMS["n3"])
    c11 = _xyz(structure, AC.REACTIVE_ATOMS["c11"])
    normal = _unit(_cross(_vec_sub(n3, o2), _vec_sub(c11, o2)))
    signed = _dot(_vec_sub(c12, o2), normal)
    projection = _vec_sub(c12, _scale(normal, signed))
    shifted = _vec_add(projection, _scale(normal, 0.20 if sign == "plus" else -0.20))
    return _dihedral(shifted, o2, n3, c11)


def torsion_restraint(target: float, force: float) -> str:
    return (
        f"&rst iat=10287,10288,10289,10286, r1={target - 15.0:.3f}, "
        f"r2={target - 2.0:.3f}, r3={target + 2.0:.3f}, r4={target + 15.0:.3f}, "
        f"rk2={force:.1f}, rk3={force:.1f}, /\n"
    )


def amber_min_input(title: str, qmmask: str, maxcyc: int, ncyc: int,
                    restrained: bool) -> str:
    disang = "\nDISANG=restraints.RST\nDUMPAVE=restraint.dat" if restrained else ""
    return f"""{title}
&cntrl
  imin=1, ntmin=2, maxcyc={maxcyc}, ncyc={ncyc}, dx0=0.005,
  ntb=1, cut=10.0, ntpr=5, ntxo=1, ifqnt=1,
  ntr={1 if restrained else 0}, nmropt={1 if restrained else 0},
  restraint_wt=1.0,
  restraintmask='{BASE.NON_QM_SOLUTE_HEAVY_MASK if restrained else ""}',
  drms=0.10,
/
{BASE.qmmm_block(qmmask)}&wt type='END' /{disang}
"""


def initialize(task_index: int, output: pathlib.Path, start_rst7: pathlib.Path,
               github_commit: str) -> None:
    task = task_spec(task_index)
    AC.initialize(task["seed_index"], output, start_rst7, github_commit)
    manifest_path = output / "ENDPOINT_MANIFEST.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    import parmed as pmd
    structure = pmd.load_file(str(AC.PRMTOP), xyz=str(start_rst7))
    manifest.update({
        "status": "READY_A1_TETRA_PARALLEL_SCOUT",
        "scope": SCOPE, "scientific_status": SCIENTIFIC_BOUNDARY,
        "next_action": NEXT_ACTION, "task": {k: v for k, v in task.items() if k != "windows"},
        "window_specs": list(task["windows"]),
        "baseline_geometry": AC._read_structure_geometry(start_rst7, manifest),
        "oop_torsion_definition": {
            "iat": [10287, 10288, 10289, 10286],
            "meaning": "torsion(C12,O2,N3,C11); actual C12-plane OOP is the scientific metric",
            "target_deg": torsion_target(structure, task["sign"]) if task["branch"] == "B_ATTACK_OOP" else None,
        },
        "results": [],
    })
    AC.write_json(manifest_path, manifest)


def prepare(task_index: int, window_index: int, root: pathlib.Path,
            output: pathlib.Path, input_rst7: pathlib.Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    if not input_rst7.is_file() or input_rst7.stat().st_size == 0:
        raise ValueError("missing inherited restart")
    manifest = json.loads((root / "ENDPOINT_MANIFEST.json").read_text(encoding="utf-8"))
    task = task_spec(task_index)
    if manifest["task"]["task_index"] != task["task_index"]:
        raise ValueError("task identity mismatch")
    windows = task["windows"]
    if not 0 <= window_index < len(windows):
        raise ValueError("window index outside branch")
    spec = dict(windows[window_index])
    output.mkdir(parents=True, exist_ok=False)
    restraints = AC.distance_restraint(AC.REACTIVE_ATOMS["og1"], AC.REACTIVE_ATOMS["c12"],
                                       spec["attack"], spec["k_attack"])
    if task["branch"] == "A_ATTACK_CO":
        restraints += AC.distance_restraint(AC.REACTIVE_ATOMS["c12"], AC.REACTIVE_ATOMS["o2"],
                                            spec["co"], spec["k_co"])
    else:
        target = float(manifest["oop_torsion_definition"]["target_deg"])
        restraints += torsion_restraint(target, spec["k_torsion"])
    (output / "stage.in").write_text(
        amber_min_input(f"NylC A1 tetra scout {task['branch']} task={task_index} window={window_index}",
                        manifest["qm_contract"]["qmmask"], spec["maxcyc"], spec["ncyc"], True),
        encoding="utf-8")
    (output / "restraints.RST").write_text(restraints, encoding="utf-8")
    AC.write_json(output / "WINDOW_MANIFEST.json", {
        "schema_version": 1, "status": "READY_A1_TETRA_SCOUT_WINDOW",
        "github_commit": manifest["github_commit"], "scope": SCOPE,
        "task": manifest["task"], "window_index": window_index, "window_spec": spec,
        "input_restart": str(input_rst7), "input_restart_sha256": AC.sha256(input_rst7),
        "prmtop": manifest["prmtop"], "prmtop_sha256": manifest["prmtop_sha256"],
        "qm_contract": manifest["qm_contract"],
        "oop_torsion_definition": manifest["oop_torsion_definition"],
        "reaction_restraints": (
            ["attack_OG1_C12", "carbonyl_C12_O2"] if task["branch"] == "A_ATTACK_CO"
            else ["attack_OG1_C12", "torsion_C12_O2_N3_C11"]
        ),
        "audit_only": ["PT2", "C12_N3", "attack_angle", "actual_C12_plane_OOP"],
    })


def _technical(stage_out: pathlib.Path, restart: pathlib.Path,
               root_manifest: Mapping[str, Any]) -> tuple[bool, dict[str, Any], dict[str, Any]]:
    text = stage_out.read_text(encoding="utf-8", errors="replace") if stage_out.is_file() else ""
    hard = AC._hard_hits(text)
    complete = "FINAL RESULTS" in text and bool(re.search(r"Run\s+done", text))
    scc = len(re.findall(r"Convergence could not be achieved", text, re.I))
    vlimit = len(re.findall(r"vlimit\s+exceeded", text, re.I))
    overflow = len(re.findall(r"BOND\s*=\s*\*+", text, re.I))
    geometry: dict[str, Any] = {}
    if restart.is_file() and restart.stat().st_size:
        geometry = AC._read_structure_geometry(restart, root_manifest)
    finite = bool(geometry) and all(
        isinstance(value, (int, float)) and math.isfinite(float(value))
        for key, value in geometry.items() if key != "hg1_nearest_qm_heavy_atom"
    )
    passed = bool(complete and restart.is_file() and restart.stat().st_size
                  and finite and scc == 0 and vlimit == 0 and overflow == 0
                  and sum(hard.values()) == 0)
    return passed, geometry, {
        "complete": complete, "scc_warnings": scc, "vlimit_warnings": vlimit,
        "bond_overflow": overflow, "hard_error_hits": hard,
    }


def geometry_gates(geometry: Mapping[str, Any], attack_target: float | None = None) -> dict[str, Any]:
    technical_geometry = bool(geometry)
    residual = (abs(float(geometry["attack_A"]) - attack_target)
                if technical_geometry and attack_target is not None else None)
    reactant_guard = bool(
        technical_geometry
        and 1.30 <= float(geometry["c12_n3_A"]) <= 1.60
        and float(geometry["nalpha_hg1_A"]) <= 1.20
        and float(geometry["qPT_A"]) <= -0.40
        and int(geometry["hg1_nearest_qm_heavy_atom"]) == 8949
    )
    tetra = bool(
        reactant_guard and float(geometry["attack_A"]) <= 1.70
        and 1.28 <= float(geometry["c12_o2_A"]) <= 1.45
        and float(geometry["c12_reactant_plane_out_of_plane_A"]) >= 0.20
        and 90.0 <= float(geometry["attack_angle_deg"]) <= 130.0
        and (residual is None or residual <= 0.20)
    )
    return {
        "attack_target_residual_A": residual,
        "attack_target_residual_le_0p20_A": bool(residual is None or residual <= 0.20),
        "attack_le_1p70_A": bool(technical_geometry and float(geometry["attack_A"]) <= 1.70),
        "c12_o2_1p28_to_1p45_A": bool(technical_geometry and 1.28 <= float(geometry["c12_o2_A"]) <= 1.45),
        "c12_plane_oop_ge_0p20_A": bool(technical_geometry and float(geometry["c12_reactant_plane_out_of_plane_A"]) >= 0.20),
        "c12_n3_1p30_to_1p60_A": bool(technical_geometry and 1.30 <= float(geometry["c12_n3_A"]) <= 1.60),
        "attack_angle_90_to_130_deg": bool(technical_geometry and 90.0 <= float(geometry["attack_angle_deg"]) <= 130.0),
        "nalpha_hg1_le_1p20_A": bool(technical_geometry and float(geometry["nalpha_hg1_A"]) <= 1.20),
        "qPT_le_minus_0p40_A": bool(technical_geometry and float(geometry["qPT_A"]) <= -0.40),
        "hg1_nearest_qm_heavy_is_nalpha": bool(technical_geometry and int(geometry["hg1_nearest_qm_heavy_atom"]) == 8949),
        "reactant_guard": reactant_guard,
        "TETRAHEDRAL_LIKE_RESTRAINED": tetra,
    }


def audit(root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path) -> None:
    root_manifest_path = root / "ENDPOINT_MANIFEST.json"
    root_manifest = json.loads(root_manifest_path.read_text(encoding="utf-8"))
    window_manifest = json.loads((output / "WINDOW_MANIFEST.json").read_text(encoding="utf-8"))
    restart = scratch / "stage.rst7"
    technical, geometry, diagnostics = _technical(scratch / "stage.out", restart, root_manifest)
    gates = geometry_gates(geometry, float(window_manifest["window_spec"]["attack"]) if technical else None)
    guard_stop = bool(technical and not gates["reactant_guard"])
    candidate = bool(technical and gates["TETRAHEDRAL_LIKE_RESTRAINED"])
    close = bool(technical and gates["attack_le_1p70_A"] and gates["attack_target_residual_le_0p20_A"])
    classification = (
        "OUT_OF_SCOPE_COUPLED_PT_OR_CN_RESPONSE" if guard_stop
        else "TETRAHEDRAL_LIKE_RESTRAINED" if candidate
        else "FORCED_CLOSE_CONTACT" if close
        else "NO_TETRAHEDRAL_LIKE_RESPONSE"
    )
    result = {
        "schema_version": 1,
        "status": "PASS_TECHNICAL_A1_TETRA_SCOUT_WINDOW" if technical else "NOT_EVALUATED_TECHNICAL_A1_TETRA_SCOUT_WINDOW",
        "scientific_status": SCIENTIFIC_BOUNDARY, "next_action": NEXT_ACTION,
        "classification": classification, "technical_pass": technical,
        "stopped_by_reactant_guard": guard_stop, "gates": gates,
        "observed": geometry, "diagnostics": diagnostics,
        "input_restart_sha256": window_manifest["input_restart_sha256"],
        "output_restart": str(restart),
        "output_restart_sha256": AC.sha256(restart) if restart.is_file() and restart.stat().st_size else None,
        "task": window_manifest["task"], "window_index": window_manifest["window_index"],
        "window_spec": window_manifest["window_spec"],
    }
    AC.write_json(output / "RESULT.json", result)
    AC.write_json(output / ("PASS.json" if technical else "NOT_EVALUATED.json"), result)
    root_manifest["results"].append({
        "window_index": result["window_index"], "classification": classification,
        "technical_pass": technical, "input_restart_sha256": result["input_restart_sha256"],
        "output_restart_sha256": result["output_restart_sha256"],
    })
    AC.write_json(root_manifest_path, root_manifest)
    if not technical:
        raise RuntimeError("tetra scout window did not technically PASS")


def prepare_release(root: pathlib.Path, output: pathlib.Path, input_rst7: pathlib.Path) -> None:
    if output.exists():
        raise FileExistsError(output)
    manifest = json.loads((root / "ENDPOINT_MANIFEST.json").read_text(encoding="utf-8"))
    output.mkdir(parents=True, exist_ok=False)
    (output / "stage.in").write_text(
        amber_min_input("NylC A1 tetrahedral candidate full release",
                        manifest["qm_contract"]["qmmask"], 800, 200, False),
        encoding="utf-8")
    AC.write_json(output / "RELEASE_MANIFEST.json", {
        "schema_version": 1, "status": "READY_FULL_RESTRAINT_FREE_RELEASE",
        "input_restart": str(input_rst7), "input_restart_sha256": AC.sha256(input_rst7),
        "ntr": 0, "nmropt": 0, "reactive_restraints": [],
        "environment_restraint": False, "scientific_status": SCIENTIFIC_BOUNDARY,
    })


def audit_release(root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path) -> None:
    root_manifest = json.loads((root / "ENDPOINT_MANIFEST.json").read_text(encoding="utf-8"))
    release_manifest = json.loads((output / "RELEASE_MANIFEST.json").read_text(encoding="utf-8"))
    restart = scratch / "stage.rst7"
    technical, geometry, diagnostics = _technical(scratch / "stage.out", restart, root_manifest)
    gates = geometry_gates(geometry, None)
    persists = bool(technical and gates["TETRAHEDRAL_LIKE_RESTRAINED"])
    result = {
        "schema_version": 1,
        "status": "PASS_TECHNICAL_FULL_RELEASE" if technical else "NOT_EVALUATED_TECHNICAL_FULL_RELEASE",
        "scientific_gate": (
            "PERSISTS_TETRAHEDRAL_AFTER_FULL_RELEASE" if persists
            else "COLLAPSES_AFTER_FULL_RELEASE" if technical
            else "NOT_EVALUATED_TECHNICAL_FULL_RELEASE"
        ),
        "technical_pass": technical, "gates": gates, "observed": geometry,
        "diagnostics": diagnostics,
        "input_restart_sha256": release_manifest["input_restart_sha256"],
        "output_restart": str(restart),
        "output_restart_sha256": AC.sha256(restart) if restart.is_file() and restart.stat().st_size else None,
        "scientific_status": SCIENTIFIC_BOUNDARY, "next_action": NEXT_ACTION,
    }
    AC.write_json(output / "RESULT.json", result)
    AC.write_json(output / ("PASS.json" if technical else "NOT_EVALUATED.json"), result)
    if not technical:
        raise RuntimeError("full release did not technically PASS")


def finalize(root: pathlib.Path, stop_reason: str, terminal_rst7: pathlib.Path | None) -> None:
    manifest = json.loads((root / "ENDPOINT_MANIFEST.json").read_text(encoding="utf-8"))
    results = [json.loads(path.read_text(encoding="utf-8"))
               for path in sorted(root.glob("window_*/RESULT.json"))]
    release_path = root / "full_release" / "RESULT.json"
    release = json.loads(release_path.read_text(encoding="utf-8")) if release_path.is_file() else None
    technical = bool(results) and all(item.get("technical_pass") is True for item in results)
    candidate = next((item for item in results if item["gates"]["TETRAHEDRAL_LIKE_RESTRAINED"]), None)
    persisted: dict[str, Any] = {}
    for label, source in (
        ("restrained_candidate", pathlib.Path(candidate["output_restart"]) if candidate else None),
        ("released_endpoint", pathlib.Path(release["output_restart"]) if release else None),
        ("terminal_endpoint", terminal_rst7),
    ):
        if source is not None and source.is_file() and source.stat().st_size:
            target = root / f"{label}.rst7"
            if not target.exists():
                shutil.copy2(source, target)
            persisted[label] = {"path": str(target), "sha256": AC.sha256(target)}
    previous = manifest["start_restart_sha256"]
    inheritance_records = []
    inheritance_ok = True
    for item in results:
        if item["input_restart_sha256"] != previous or not item["output_restart_sha256"]:
            inheritance_ok = False
        inheritance_records.append({
            "window_index": item["window_index"],
            "input_sha256": item["input_restart_sha256"],
            "output_sha256": item["output_restart_sha256"],
        })
        previous = item["output_restart_sha256"]
    scientific_gate = (
        release["scientific_gate"] if release is not None
        else "RESTRAINED_TETRAHEDRAL_CANDIDATE_NOT_RELEASED" if candidate
        else "NO_TETRAHEDRAL_LIKE_RESPONSE" if technical
        else "NOT_EVALUATED_TECHNICAL_A1_TETRA_SCOUT"
    )
    chain = {
        "schema_version": 1,
        "status": "PASS_TECHNICAL_A1_TETRA_SCOUT" if technical else "NOT_EVALUATED_TECHNICAL_A1_TETRA_SCOUT",
        "scientific_gate": scientific_gate, "scientific_status": SCIENTIFIC_BOUNDARY,
        "next_action": NEXT_ACTION, "scope": SCOPE, "task": manifest["task"],
        "stop_reason": stop_reason, "window_count": len(results),
        "inheritance_sha256_verified": inheritance_ok,
        "inheritance": inheritance_records, "windows": results,
        "full_release": release, "persisted_restarts": persisted,
    }
    AC.write_json(root / "CHAIN_MANIFEST.json", chain)
    hashes = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and (
            path.suffix == ".json" or path.suffix == ".rst7"
        ):
            hashes.append(f"{AC.sha256(path)}  {path.relative_to(root)}")
    (root / "SHA256.tsv").write_text("\n".join(hashes) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True,
                        choices=("initialize", "prepare", "audit", "prepare-release",
                                 "audit-release", "finalize", "describe"))
    parser.add_argument("--task-index", type=int, required=True)
    parser.add_argument("--window-index", type=int)
    parser.add_argument("--root", type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--scratch", type=pathlib.Path)
    parser.add_argument("--input-rst7", type=pathlib.Path)
    parser.add_argument("--start-rst7", type=pathlib.Path)
    parser.add_argument("--github-commit", default="unknown")
    parser.add_argument("--stop-reason", default="UNKNOWN")
    parser.add_argument("--terminal-rst7", type=pathlib.Path)
    args = parser.parse_args()
    if args.mode == "describe":
        task = task_spec(args.task_index)
        print(json.dumps({"task": {k: v for k, v in task.items() if k != "windows"},
                          "window_count": len(task["windows"])}))
    elif args.mode == "initialize":
        initialize(args.task_index, args.root.resolve(), args.start_rst7.resolve(), args.github_commit)
    elif args.mode == "prepare":
        prepare(args.task_index, args.window_index, args.root.resolve(),
                args.output.resolve(), args.input_rst7.resolve())
    elif args.mode == "audit":
        audit(args.root.resolve(), args.output.resolve(), args.scratch.resolve())
    elif args.mode == "prepare-release":
        prepare_release(args.root.resolve(), args.output.resolve(), args.input_rst7.resolve())
    elif args.mode == "audit-release":
        audit_release(args.root.resolve(), args.output.resolve(), args.scratch.resolve())
    else:
        finalize(args.root.resolve(), args.stop_reason,
                 args.terminal_rst7.resolve() if args.terminal_rst7 else None)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
