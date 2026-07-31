#!/usr/bin/env python3
"""Two-axis unrestrained A1 authority control: Thr267 protonation vs QM boundary."""
from __future__ import annotations

import argparse
import copy
import importlib.util
import json
import os
import pathlib
import re
import shutil
from typing import Any, Mapping

HERE = pathlib.Path(__file__).resolve().parent
ACT_PATH = HERE / "prepare_audit_nylc_a1_step1_activated_mmframe_addition.py"


def _load(name: str, path: pathlib.Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ACT = _load("_a1_activated_mmframe", ACT_PATH)
BASE = ACT.BASE
RAW_FORWARD = ACT.RAW_FORWARD
THR267 = dict(ACT.THR267)
REACTIVE = dict(ACT.REACTIVE)
TASK_ROOT = ACT.TASK_ROOT
PRMTOP = ACT.PRMTOP
PRMTOP_SHA256 = ACT.PRMTOP_SHA256
SOURCES = ACT.SOURCES

ARRAY_TASKS = 4
MPI_RANKS = 8
AUTHORITY_BLOCKS = 4
BLOCK_MAXCYC = 250
BLOCK_NCYC = 75
NEUTRAL_OH_A = 0.96
OUTPUT_ROOT = (
    TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/a1_step1_protonation_boundary_authority"
)

CURRENT_QMMASK = (
    "@7160-7174,7756-7771,8235-8242,8949-8963,"
    "9567-9573,9587-9592,10273-10351"
)
EXPANDED_THR268_ATOMS = tuple(range(8964, 8978))
EXPANDED_QMMASK = (
    "@7160-7174,7756-7771,8235-8242,8949-8963,8964-8977,"
    "9567-9573,9587-9592,10273-10351"
)
QM_CONTRACTS = {
    "current": {"qm_atoms": 146, "qm_charge": 0, "electrons": 510, "link_atoms": 6},
    "expanded_thr268": {
        "qm_atoms": 160, "qm_charge": 0, "electrons": 564, "link_atoms": 6
    },
}
CURRENT_BOUNDARIES = (
    (7158, 7160), (7754, 7756), (8233, 8235),
    (8962, 8964), (9565, 9567), (9585, 9587),
)
EXPANDED_BOUNDARIES = (
    (7158, 7160), (7754, 7756), (8233, 8235),
    (8976, 8978), (9565, 9567), (9585, 9587),
)
REQUIRED_HEAVY_BONDS = (
    frozenset((8949, 8953)),  # Thr267 N-CA
    frozenset((8953, 8955)),  # CA-CB
    frozenset((8955, 8960)),  # CB-OG1
    frozenset((8962, 8964)),  # Thr267 C - Thr268 N
)
REQUIRED_NH_BONDS = (
    frozenset((8949, 8950)),
    frozenset((8949, 8951)),
)
SCIENTIFIC_STATUS = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..3")
    seed = (26723, 26737)[index % 2]
    if index < 2:
        return {
            "task_index": index,
            "seed": seed,
            "variant": "NEUTRAL_THR267_CURRENT_QM",
            "starting_state": "NALPHA_H2_OG1H",
            "qm_contract_key": "current",
            "qmmask": CURRENT_QMMASK,
        }
    return {
        "task_index": index,
        "seed": seed,
        "variant": "ACTIVATED_THR267_EXPANDED_THR268_QM",
        "starting_state": "NALPHA_H3_PLUS_OG1_MINUS",
        "qm_contract_key": "expanded_thr268",
        "qmmask": EXPANDED_QMMASK,
    }


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "variants": [
            "NEUTRAL_THR267_CURRENT_QM",
            "ACTIVATED_THR267_EXPANDED_THR268_QM",
        ],
        "authority_blocks": AUTHORITY_BLOCKS,
        "maxcyc_per_block": BLOCK_MAXCYC,
        "reaction_coordinate_restrained": False,
        "position_restrained": False,
        "strict_serial_restart_inheritance": True,
        "failed_restart_inheritance": False,
        "automatic_downstream_action": "NONE",
    }


def place_hg1_on_og1(
    coordinates: Any, og1_atom: int, nalpha_atom: int, hg1_atom: int
):
    """Place HG1 at a neutral O-H distance; atom numbers are one based."""
    import numpy as np

    moved = np.asarray(coordinates, dtype=float).copy()
    og = moved[int(og1_atom) - 1]
    hg = moved[int(hg1_atom) - 1]
    direction = hg - og
    norm = float(np.linalg.norm(direction))
    if norm < 1.0e-8:
        direction = moved[int(nalpha_atom) - 1] - og
        norm = float(np.linalg.norm(direction))
    if norm < 1.0e-8:
        raise ValueError("cannot determine neutral HG1 placement direction")
    moved[int(hg1_atom) - 1] = og + NEUTRAL_OH_A * direction / norm
    return moved


def classify_hg1_state(
    og1_hg1_A: float, nalpha_hg1_A: float, n3_hg1_A: float
) -> str:
    distances = {
        "OG1H_NEUTRAL": float(og1_hg1_A),
        "NALPHA_H3_ACTIVATED": float(nalpha_hg1_A),
        "N3H_TRANSFERRED": float(n3_hg1_A),
    }
    state, distance = min(distances.items(), key=lambda item: item[1])
    return state if distance <= 1.25 else "OTHER_OR_SHARED"


def _mask_atoms(mask: str) -> set[int]:
    text = mask.strip()
    if not text.startswith("@"):
        raise ValueError("qmmask must be atom based")
    atoms: set[int] = set()
    for token in text[1:].split(","):
        if "-" in token:
            left, right = (int(value) for value in token.split("-", 1))
            atoms.update(range(left, right + 1))
        else:
            atoms.add(int(token))
    return atoms


def _boundary_pairs(structure: Any, qmmask: str) -> tuple[tuple[int, int], ...]:
    qm = _mask_atoms(qmmask)
    pairs = []
    for bond in structure.bonds:
        left, right = bond.atom1.idx + 1, bond.atom2.idx + 1
        if (left in qm) != (right in qm):
            pairs.append((left, right) if left not in qm else (right, left))
    return tuple(sorted(pairs))


def _contract(task: Mapping[str, Any], source: Mapping[str, Any]) -> dict[str, Any]:
    current = ACT._contract(source)
    expected = QM_CONTRACTS[str(task["qm_contract_key"])]
    contract = copy.deepcopy(current)
    contract["qmmask"] = str(task["qmmask"])
    contract["qm_atom_count"] = expected["qm_atoms"]
    contract["qmcharge"] = expected["qm_charge"]
    contract["electron_count_including_link_h"] = expected["electrons"]
    contract["link_atom_count"] = expected["link_atoms"]
    contract["step1_qm_water_count"] = 0
    if task["qm_contract_key"] == "expanded_thr268":
        heavy = set(int(value) for value in contract.get("qm_heavy_atom_indices", []))
        heavy.update((8964, 8966, 8968, 8970, 8974, 8976, 8977))
        contract["qm_heavy_atom_indices"] = sorted(heavy)
    return contract


def _minimum_distance(coordinates: Any, left: int, right: int, box: Any) -> float:
    import numpy as np
    from MDAnalysis.lib.distances import minimize_vectors

    delta = np.asarray(coordinates[right - 1] - coordinates[left - 1], dtype=float)
    if box is not None:
        delta = minimize_vectors(delta.reshape(1, 3), box)[0]
    return float(np.linalg.norm(delta))


def measure_integrity(prmtop: pathlib.Path, restart: pathlib.Path) -> dict[str, Any]:
    import parmed as pmd

    structure = pmd.load_file(str(prmtop), str(restart))
    coordinates = structure.coordinates
    box = structure.box
    pairs = {
        "Nalpha-CA": (8949, 8953, 1.85),
        "CA-CB": (8953, 8955, 1.95),
        "CB-OG1": (8955, 8960, 1.85),
        "C267-N268": (8962, 8964, 1.85),
        "Nalpha-H1": (8949, 8950, 1.35),
        "Nalpha-H2": (8949, 8951, 1.35),
    }
    distances = {
        name: _minimum_distance(coordinates, left, right, box)
        for name, (left, right, _limit) in pairs.items()
    }
    checks = {
        name: distances[name] <= limit
        for name, (_left, _right, limit) in pairs.items()
    }
    proton = {
        "OG1-HG1": _minimum_distance(coordinates, 8960, 8961, box),
        "Nalpha-HG1": _minimum_distance(coordinates, 8949, 8961, box),
        "N3-HG1": _minimum_distance(
            coordinates, int(REACTIVE["n3"]), 8961, box
        ),
    }
    return {
        "pass": all(checks.values()),
        "distances_A": distances,
        "checks": checks,
        "proton_distances_A": proton,
        "proton_state": classify_hg1_state(
            proton["OG1-HG1"], proton["Nalpha-HG1"], proton["N3-HG1"]
        ),
    }


def _source_geometry_gate(geometry: Mapping[str, Any]) -> dict[str, Any]:
    checks = {
        "attack_distance": float(geometry["attack_A"]) <= 3.5,
        "attack_angle": 95.0 <= float(geometry["attack_angle_deg"]) <= 115.0,
    }
    return {"pass": all(checks.values()), "checks": checks}


def initialize(
    task_index: int, root: pathlib.Path, github_commit: str
) -> dict[str, Any]:
    if root.exists():
        raise FileExistsError(root)
    if not re.fullmatch(r"[0-9a-f]{40}", github_commit):
        raise ValueError("runtime is not bound to a full Git commit")
    if ACT.sha256(PRMTOP) != PRMTOP_SHA256:
        raise ValueError("frozen prmtop SHA changed")

    task = task_spec(task_index)
    source = SOURCES[int(task["seed"])]
    for key in ("tpr", "xtc"):
        path = pathlib.Path(source[key])
        if not path.is_file() or ACT.sha256(path) != source[f"{key}_sha256"]:
            raise ValueError(f"fixed {key} authority failed")

    import MDAnalysis as mda
    import parmed as pmd

    universe = mda.Universe(str(source["tpr"]), str(source["xtc"]))
    universe.trajectory[int(source["frame"])]
    structure = pmd.load_file(str(PRMTOP))
    if len(universe.atoms) != len(structure.atoms):
        raise ValueError("TPR/XTC atom count differs from frozen Amber topology")
    dimensions = universe.dimensions
    if dimensions is None or len(dimensions) != 6:
        raise ValueError("source frame lacks triclinic periodic box")
    bonds = [(bond.atom1.idx, bond.atom2.idx) for bond in structure.bonds]
    coordinates = ACT.make_bonded_fragments_whole(
        universe.atoms.positions.copy(), bonds, dimensions
    )
    if task["variant"] == "NEUTRAL_THR267_CURRENT_QM":
        coordinates = place_hg1_on_og1(
            coordinates, THR267["og1"], THR267["nalpha"], THR267["hg1"]
        )
    structure.coordinates = coordinates
    structure.box = [float(value) for value in dimensions]

    expected_boundaries = (
        CURRENT_BOUNDARIES
        if task["qm_contract_key"] == "current"
        else EXPANDED_BOUNDARIES
    )
    actual_boundaries = _boundary_pairs(structure, str(task["qmmask"]))
    boundary_ok = actual_boundaries == tuple(sorted(expected_boundaries))
    contract = _contract(task, source)

    root.mkdir(parents=True, exist_ok=False)
    restart = root / "source.rst7"
    structure.save(str(restart), overwrite=False)
    geometry = BASE.AUTH._geometry(restart, {"qm_contract": contract})
    integrity = measure_integrity(PRMTOP, restart)
    geometry_gate = _source_geometry_gate(geometry)
    ready = bool(boundary_ok and integrity["pass"] and geometry_gate["pass"])
    manifest = {
        "schema_version": 1,
        "status": "READY_A1_PROTONATION_BOUNDARY_AUTHORITY"
                  if ready else "NOT_EVALUATED_TECHNICAL_SOURCE_INTEGRITY",
        "github_commit": github_commit,
        "task": task,
        "source": {
            "tpr": str(source["tpr"]),
            "xtc": str(source["xtc"]),
            "tpr_sha256": source["tpr_sha256"],
            "xtc_sha256": source["xtc_sha256"],
            "frame": source["frame"],
            "time_ps": source["time_ps"],
            "restart": str(restart),
            "restart_sha256": ACT.sha256(restart),
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": ACT.sha256(PRMTOP),
        "qm_contract": contract,
        "expected_boundary_pairs": expected_boundaries,
        "actual_boundary_pairs": actual_boundaries,
        "boundary_gate": boundary_ok,
        "source_geometry": geometry,
        "source_geometry_gate": geometry_gate,
        "source_integrity": integrity,
        "ready": ready,
        "accepted_restart": str(restart) if ready else None,
        "accepted_restart_sha256": ACT.sha256(restart) if ready else None,
        "blocks": [],
        "terminal": not ready,
        "terminal_reason": None if ready else
                           "NOT_EVALUATED_TECHNICAL_SOURCE_INTEGRITY",
        "scientific_status": SCIENTIFIC_STATUS,
        "automatic_downstream_action": "NONE",
    }
    ACT.write_json(root / "SOURCE_MANIFEST.json", manifest)
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


def authority_input(qmmask: str) -> str:
    return f"""NylC A1 protonation-boundary unrestrained authority
&cntrl
  imin=1, ntmin=2, maxcyc={BLOCK_MAXCYC}, ncyc={BLOCK_NCYC},
  dx0=0.005, drms=0.01,
  ntb=1, cut=10.0, ntpr=5, ntxo=1,
  ifqnt=1, nmropt=0, ntr=0,
/
{_qmmm_block(qmmask)}"""


def prepare_block(
    root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path, index: int
) -> dict[str, Any]:
    manifest = BASE.read_json(root / "SOURCE_MANIFEST.json")
    block = int(index)
    if not 1 <= block <= AUTHORITY_BLOCKS:
        raise ValueError("authority block index must be 1..4")
    if manifest["terminal"] or not manifest["ready"]:
        raise ValueError("authority chain is terminal")
    if block != len(manifest["blocks"]) + 1:
        raise ValueError("authority block order is not strict")
    input_restart = pathlib.Path(manifest["accepted_restart"])
    if ACT.sha256(input_restart) != manifest["accepted_restart_sha256"]:
        raise ValueError("accepted restart SHA inheritance failed")
    output.mkdir(parents=True, exist_ok=False)
    scratch.mkdir(parents=True, exist_ok=False)
    text = authority_input(manifest["qm_contract"]["qmmask"])
    if "nmropt=0" not in text or "ntr=0" not in text or "DISANG" in text:
        raise ValueError("authority input is not unrestrained")
    (scratch / "block.in").write_text(text, encoding="utf-8")
    prepared = {
        "block_index": block,
        "input_restart": str(input_restart),
        "input_restart_sha256": manifest["accepted_restart_sha256"],
        "strict_inheritance": True,
        "reaction_coordinate_restraints": 0,
        "position_restrained": False,
    }
    ACT.write_json(output / "BLOCK_MANIFEST.json", prepared)
    return prepared


def audit_block(
    root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path
) -> dict[str, Any]:
    manifest = BASE.read_json(root / "SOURCE_MANIFEST.json")
    prepared = BASE.read_json(output / "BLOCK_MANIFEST.json")
    restart = scratch / "block.rst7"
    engine_ok, geometry, diagnostics, numerical = ACT._engine_evidence(
        scratch / "block.out", restart, manifest
    )
    integrity = {"pass": False}
    if engine_ok:
        integrity = measure_integrity(pathlib.Path(manifest["prmtop"]), restart)
    output_sha = ACT.sha256(restart) if restart.is_file() else None
    new_restart = output_sha not in (None, prepared["input_restart_sha256"])
    accepted = bool(engine_ok and integrity["pass"] and new_restart)
    result = {
        "schema_version": 1,
        "block_index": prepared["block_index"],
        "status": "PASS_A1_UNRESTRAINED_AUTHORITY_BLOCK"
                  if accepted else "NOT_EVALUATED_TECHNICAL_A1_AUTHORITY_BLOCK",
        "technical_complete": engine_ok,
        "numerical_health": numerical,
        "chemical_integrity": integrity,
        "geometry": geometry if engine_ok else None,
        "input_restart_sha256": prepared["input_restart_sha256"],
        "output_restart_sha256": output_sha,
        "new_restart_sha": new_restart,
        "accepted": accepted,
        "diagnostics": diagnostics,
    }
    ACT._copy_engine_evidence(output, scratch, "block")
    if accepted:
        persistent = output / "block.rst7"
        shutil.copy2(restart, persistent)
        manifest["accepted_restart"] = str(persistent)
        manifest["accepted_restart_sha256"] = ACT.sha256(persistent)
    else:
        manifest["terminal"] = True
        manifest["terminal_reason"] = result["status"]
    manifest["blocks"].append({
        "block_index": prepared["block_index"],
        "result": str(output / "RESULT.json"),
        "accepted": accepted,
        "input_restart_sha256": prepared["input_restart_sha256"],
        "output_restart_sha256": output_sha,
    })
    ACT.write_json(output / "RESULT.json", result)
    ACT.write_json(root / "SOURCE_MANIFEST.json", manifest)
    return result


def finalize(root: pathlib.Path) -> dict[str, Any]:
    manifest = BASE.read_json(root / "SOURCE_MANIFEST.json")
    accepted = sum(bool(row["accepted"]) for row in manifest["blocks"])
    complete = bool(manifest["ready"] and accepted == AUTHORITY_BLOCKS)
    last = None
    if manifest["blocks"]:
        last = BASE.read_json(pathlib.Path(manifest["blocks"][-1]["result"]))
    status = (
        "PASS_A1_PROTONATION_BOUNDARY_AUTHORITY"
        if complete else
        manifest.get("terminal_reason") or "NOT_EVALUATED_TECHNICAL_INCOMPLETE"
    )
    result = {
        "schema_version": 1,
        "status": status,
        "task": manifest["task"],
        "technical_complete": complete,
        "accepted_blocks": accepted,
        "expected_blocks": AUTHORITY_BLOCKS,
        "strict_inheritance": True,
        "source_integrity": manifest["source_integrity"],
        "terminal_geometry": last.get("geometry") if last else None,
        "terminal_chemical_integrity": (
            last.get("chemical_integrity") if last else None
        ),
        "terminal_reason": manifest.get("terminal_reason"),
        "scientific_status": SCIENTIFIC_STATUS,
        "automatic_downstream_action": "NONE",
    }
    ACT.write_json(root / "RESULT.json", result)
    ACT.write_json(
        root / ("PASS.json" if complete else "NOT_EVALUATED.json"), result
    )
    ACT._write_hashes(root)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode", required=True,
        choices=("describe", "initialize", "prepare-block", "audit-block", "finalize")
    )
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--root", type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--scratch", type=pathlib.Path)
    parser.add_argument("--block-index", type=int)
    parser.add_argument("--github-commit", default="unknown")
    args = parser.parse_args()
    if args.mode == "describe":
        print(json.dumps(describe(), sort_keys=True))
    elif args.mode == "initialize":
        initialize(args.task_index, args.root.resolve(), args.github_commit)
    elif args.mode == "prepare-block":
        prepare_block(
            args.root.resolve(), args.output.resolve(), args.scratch.resolve(),
            args.block_index,
        )
    elif args.mode == "audit-block":
        audit_block(
            args.root.resolve(), args.output.resolve(), args.scratch.resolve()
        )
    else:
        finalize(args.root.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
