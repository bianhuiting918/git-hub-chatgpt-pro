#!/usr/bin/env python3
"""A1 covalent-integrity reclassification and unrestrained authority diagnostic."""
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
PREORG_PATH = HERE / "prepare_audit_nylc_a1_step1_raw_nac_preorganized_calibration.py"
PBC_PATH = HERE / "prepare_audit_nylc_a1_step2_qmwater_endpoint.py"


def _load(name: str, path: pathlib.Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot import {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


PREORG = _load("_a1_preorg", PREORG_PATH)
PBC = _load("_a1_pbc", PBC_PATH)
FWD = PREORG.RAW_FORWARD
BASE = PREORG.BASE
TASK_ROOT = PREORG.TASK_ROOT
PRMTOP = PREORG.PRMTOP
PRMTOP_SHA256 = PREORG.PRMTOP_SHA256
SOURCES = PREORG.SOURCES
QM_CONTRACT = PREORG.QM_CONTRACT

ARRAY_TASKS = 2
MPI_RANKS = 8
OUTPUT_ROOT = (
    TASK_ROOT / "a1_activated_nac_20260726/qmmm"
    / "a1_step1_a1_unbiased_authority"
)
LEGACY_OUTPUT_ROOT = PREORG.OUTPUT_ROOT
LEGACY_ARRAY_JOB = "62500360"
N_ALPHA = 8949
H1 = 8950
H2 = 8951
CA = 8952
CB = 8954
OG1 = 8960
HG1 = 8961
HYDROGENS = {"H1": H1, "H2": H2, "HG1": HG1}
ATOM_IDENTITIES = {
    N_ALPHA: ("THR", "N"),
    H1: ("THR", "H1"),
    H2: ("THR", "H2"),
    CA: ("THR", "CA"),
    CB: ("THR", "CB"),
    OG1: ("THR", "OG1"),
    HG1: ("THR", "HG1"),
}
SCIENTIFIC_STATUS = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "mapping": "2 original hash-fixed raw-NAC seeds",
        "source_restart_shas": {
            str(seed): source["start_rst7_sha256"]
            for seed, source in SOURCES.items()
        },
        "qm_contract": dict(QM_CONTRACT),
        "reaction_coordinate_restraints": 0,
        "attack_restraints": 0,
        "carbonyl_restraints": 0,
        "angle_restraints": 0,
        "qpt_restraints": 0,
        "criterion": "A1_THR267_COVALENT_INTEGRITY_V1",
        "legacy_result_mutation": False,
        "automatic_downstream_action": "NONE",
    }


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..1")
    seed = (26723, 26737)[index]
    return {
        "task_index": index,
        "seed": seed,
        "seed_index": PREORG.RAW_FORWARD.SOURCES[seed]["seed_index"],
        "candidate": PREORG.RAW_FORWARD.SOURCES[seed]["candidate"],
        "mode": "A1_REACTION_COORDINATE_FREE_AUTHORITY",
        "source_basin": "HASH_FIXED_RAW_UNBIASED_NAC",
        "progression_restraints": 0,
    }


def chemical_integrity_from_metrics(
    metrics: Mapping[str, Any], nearest_acceptors: Mapping[str, Any]
) -> dict[str, Any]:
    values = [float(value) for value in metrics.values()]
    checks = {
        "finite_metrics": all(math.isfinite(value) for value in values),
        "nalpha_ca_covalent": 1.20 <= float(metrics["nalpha_ca_A"]) <= 1.75,
        "ca_cb_covalent": 1.20 <= float(metrics["ca_cb_A"]) <= 1.85,
        "cb_og1_covalent": 1.15 <= float(metrics["cb_og1_A"]) <= 1.75,
        "nalpha_h1_covalent": 0.80 <= float(metrics["nalpha_h1_A"]) <= 1.30,
        "nalpha_h2_covalent": 0.80 <= float(metrics["nalpha_h2_A"]) <= 1.30,
        "nalpha_hg1_covalent": 0.80 <= float(metrics["nalpha_hg1_A"]) <= 1.30,
        "all_nalpha_hydrogens_attached": all(
            nearest_acceptors.get(name) == "Nalpha" for name in HYDROGENS
        ),
    }
    return {
        "pass": all(checks.values()),
        "criterion": "A1_THR267_COVALENT_INTEGRITY_V1",
        "checks": checks,
        "metrics_A": {key: float(value) for key, value in metrics.items()},
        "nearest_qm_heavy_acceptors": dict(nearest_acceptors),
        "thresholds_A": {
            "nalpha_ca": [1.20, 1.75],
            "ca_cb": [1.20, 1.85],
            "cb_og1": [1.15, 1.75],
            "nalpha_h": [0.80, 1.30],
        },
    }


def measure_covalent_integrity(
    prmtop: pathlib.Path, restart: pathlib.Path, contract: Mapping[str, Any]
) -> dict[str, Any]:
    import parmed as pmd

    structure = pmd.load_file(str(prmtop), xyz=str(restart))
    if structure.box is None:
        raise ValueError("restart lacks periodic box")
    for index1, (resname, atom_name) in ATOM_IDENTITIES.items():
        atom = structure.atoms[index1 - 1]
        if atom.residue.name != resname or atom.name != atom_name:
            raise ValueError(
                f"atom identity mismatch at {index1}: "
                f"{atom.residue.name}:{atom.name} != {resname}:{atom_name}"
            )
    cell = PBC._periodic_cell(structure.box)
    atoms = structure.atoms

    def distance(left: int, right: int) -> float:
        return PBC._distance(atoms[left - 1], atoms[right - 1], cell)

    metrics = {
        "nalpha_ca_A": distance(N_ALPHA, CA),
        "ca_cb_A": distance(CA, CB),
        "cb_og1_A": distance(CB, OG1),
        "nalpha_h1_A": distance(N_ALPHA, H1),
        "nalpha_h2_A": distance(N_ALPHA, H2),
        "nalpha_hg1_A": distance(N_ALPHA, HG1),
    }
    heavy_indices = [
        int(index) for index in contract.get("qm_heavy_atom_indices", [])
        if int(index) not in HYDROGENS.values()
    ]
    if N_ALPHA not in heavy_indices:
        raise ValueError("QM heavy-atom authority lacks Nalpha")
    nearest = {}
    nearest_details = {}
    for name, hydrogen in HYDROGENS.items():
        ranked = sorted(
            (distance(hydrogen, heavy), heavy) for heavy in heavy_indices
        )
        nearest_distance, nearest_index = ranked[0]
        nearest[name] = "Nalpha" if nearest_index == N_ALPHA else (
            f"{atoms[nearest_index - 1].residue.name}:"
            f"{atoms[nearest_index - 1].name}@{nearest_index}"
        )
        nearest_details[name] = {
            "index1": nearest_index,
            "label": nearest[name],
            "distance_A": nearest_distance,
        }
    result = chemical_integrity_from_metrics(metrics, nearest)
    result["nearest_qm_heavy_details"] = nearest_details
    result["restart"] = str(restart)
    result["restart_sha256"] = sha256(restart)
    return result


def reclassify_legacy_result(
    legacy: Mapping[str, Any], integrity: Mapping[str, Any]
) -> dict[str, Any]:
    row = json.loads(json.dumps(legacy))
    row["legacy_status_before_chemical_integrity_reclassification"] = legacy.get(
        "status"
    )
    row["legacy_scientific_gate_before_chemical_integrity_reclassification"] = (
        legacy.get("scientific_gate")
    )
    row["legacy_result_preserved"] = True
    row["chemical_integrity_v1"] = dict(integrity)
    if integrity.get("pass") is not True:
        row["status"] = "NOT_EVALUATED_TECHNICAL_CHEMICAL_INTEGRITY"
        row["scientific_gate"] = "NOT_EVALUATED_TECHNICAL_CHEMICAL_INTEGRITY"
        row["technical_complete"] = False
        row["eligible_by_absolute_response"] = False
        row["eligible_for_inherited_chain"] = False
        row["eligible_for_inherited_chain_v2"] = False
    return row


def initialize(task_index: int, root: pathlib.Path, commit: str) -> dict[str, Any]:
    if not re.fullmatch(r"[0-9a-f]{40}", commit):
        raise ValueError("runtime is not bound to a full Git commit")
    spec = task_spec(task_index)
    mapped_preorg_index = 0 if spec["seed"] == 26723 else 9
    manifest = PREORG.initialize(mapped_preorg_index, root, commit)
    manifest["status"] = "READY_A1_REACTION_COORDINATE_FREE_AUTHORITY"
    manifest["task"] = spec
    manifest["reaction_coordinate_restraints"] = 0
    manifest["automatic_downstream_action"] = "NONE"
    BASE.write_json(root / "MANIFEST.json", manifest)
    return manifest


def prepare(root: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest = BASE.read_json(root / "MANIFEST.json")
    if scratch.exists():
        raise FileExistsError(scratch)
    scratch.mkdir(parents=True, exist_ok=False)
    seed = manifest["task"]["seed"]
    mapped = PREORG.task_spec(0 if seed == 26723 else 9)
    stage = PREORG.stage_spec(mapped, manifest["source_geometry"])
    stage_input = FWD._baseline_input(
        mapped, stage, manifest["qm_contract"]["qmmask"]
    )
    forbidden = ("nmropt=1", "DISANG=", "DUMPAVE=", "&wt")
    if any(token in stage_input for token in forbidden):
        raise ValueError("reaction-coordinate restraint authority leaked into A1 control")
    if "nmropt=0" not in stage_input:
        raise ValueError("A1 control is not explicitly reaction-coordinate free")
    (scratch / "stage.in").write_text(stage_input, encoding="utf-8")
    payload = {
        "schema_version": 1,
        "mode": "A1_REACTION_COORDINATE_FREE_AUTHORITY",
        "input_restart": manifest["source"]["restart"],
        "input_restart_sha256": manifest["source"]["restart_sha256"],
        "maxcyc": stage["maxcyc"],
        "ncyc": stage["ncyc"],
        "reaction_coordinate_restraints": 0,
        "nmropt": 0,
        "automatic_downstream_action": "NONE",
    }
    BASE.write_json(root / "WINDOW_MANIFEST.json", payload)
    return payload


def _write_hashes(root: pathlib.Path) -> None:
    rows = []
    for path in sorted(root.iterdir()):
        if path.is_file() and path.name != "SHA256.tsv":
            rows.append(f"{sha256(path)}  {path.name}")
    (root / "SHA256.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")


def audit(root: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest = BASE.read_json(root / "MANIFEST.json")
    stage_out = scratch / "stage.out"
    technical, geometry, diagnostics = BASE.AUTH.TETRA._technical(
        stage_out, scratch / "stage.rst7", manifest
    )
    engine_text = (
        stage_out.read_text(encoding="utf-8", errors="replace")
        if stage_out.is_file() else ""
    )
    numerical = FWD.engine_numerical_health(engine_text)
    technical = bool(technical and numerical["pass"])
    diagnostics = dict(diagnostics)
    diagnostics["numerical_health"] = numerical
    integrity = None
    if technical:
        integrity = measure_covalent_integrity(
            PRMTOP, scratch / "stage.rst7", manifest["qm_contract"]
        )
        shutil.copy2(scratch / "stage.rst7", root / "stage.rst7")
    for name in ("stage.in", "stage.mdinfo"):
        path = scratch / name
        if path.is_file():
            shutil.copy2(path, root / name)
    if stage_out.is_file():
        (root / "ENGINE_TAIL.txt").write_text(
            "\n".join(engine_text.splitlines()[-240:]) + "\n",
            encoding="utf-8",
        )

    integrity_pass = bool(technical and integrity and integrity["pass"])
    if not technical:
        status = "NOT_EVALUATED_TECHNICAL_ENGINE_OR_NUMERICAL_FAILURE"
        gate = status
    elif not integrity_pass:
        status = "NOT_EVALUATED_TECHNICAL_CHEMICAL_INTEGRITY"
        gate = status
    else:
        status = "PASS_A1_UNBIASED_CHEMICAL_INTEGRITY_AUTHORITY"
        gate = status
    result = {
        "schema_version": 1,
        "status": status,
        "technical_execution_complete": technical,
        "chemical_integrity_pass": integrity_pass,
        "scientific_gate": gate,
        "scientific_status": SCIENTIFIC_STATUS,
        "task": manifest["task"],
        "source_geometry": manifest["source_geometry"],
        "final_reactive_geometry": geometry if technical else None,
        "chemical_integrity_v1": integrity,
        "diagnostics": diagnostics,
        "reaction_coordinate_restraints": 0,
        "eligible_for_progression_calibration": integrity_pass,
        "automatic_downstream_action": "NONE",
    }
    BASE.write_json(root / "RESULT.json", result)
    BASE.write_json(
        root / ("PASS.json" if integrity_pass else "NOT_EVALUATED.json"), result
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
    integrity_pass = sum(
        row.get("chemical_integrity_pass") is True for row in results
    )
    payload = {
        "schema_version": 1,
        "criterion": "A1_THR267_COVALENT_INTEGRITY_V1",
        "status": (
            "PASS_BOTH_SEEDS_A1_UNBIASED_CHEMICAL_INTEGRITY_AUTHORITY"
            if integrity_pass == ARRAY_TASKS else
            "NOT_EVALUATED_TECHNICAL_A1_UNBIASED_CHEMICAL_INTEGRITY"
        ),
        "denominator_seeds": ARRAY_TASKS,
        "chemical_integrity_pass_seeds": integrity_pass,
        "per_seed": results,
        "next_action": (
            "REVIEW_PROGRESSION_RESTRAINTS_AND_WINDOW_SIZE"
            if integrity_pass == ARRAY_TASKS else
            "STOP_PREFORMED_A1_AND_VERSION_NEUTRAL_THR_CONCERTED_QPT_ATTACK"
        ),
        "automatic_downstream_action": "NONE",
        "scientific_status": SCIENTIFIC_STATUS,
    }
    audit_path = (
        output_root / "audit"
        / f"nylc_a1_unbiased_chemical_integrity_authority_{array_job}.json"
    )
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    if audit_path.exists() and BASE.read_json(audit_path) != payload:
        raise FileExistsError(audit_path)
    if not audit_path.exists():
        BASE.write_json(audit_path, payload)
    return True


def reclassify_legacy(output_root: pathlib.Path, array_job: str) -> pathlib.Path:
    legacy_audit = (
        output_root / "audit"
        / f"nylc_a1_raw_nac_preorganized_calibration_{array_job}.json"
    )
    if not legacy_audit.is_file():
        raise FileNotFoundError(legacy_audit)
    rows = []
    for index in range(18):
        attempt = output_root / f"attempt_{array_job}_{index}"
        result_path = attempt / "RESULT.json"
        restart = attempt / "stage.rst7"
        manifest_path = attempt / "MANIFEST.json"
        if not all(path.is_file() for path in (result_path, restart, manifest_path)):
            raise FileNotFoundError(f"incomplete legacy attempt {index}")
        legacy = BASE.read_json(result_path)
        manifest = BASE.read_json(manifest_path)
        integrity = measure_covalent_integrity(
            PRMTOP, restart, manifest["qm_contract"]
        )
        row = reclassify_legacy_result(legacy, integrity)
        row["legacy_task_index"] = index
        rows.append(row)
    failed = [
        row["legacy_task_index"] for row in rows
        if row["chemical_integrity_v1"]["pass"] is not True
    ]
    payload = {
        "schema_version": 2,
        "criterion": "A1_THR267_COVALENT_INTEGRITY_V1",
        "status": (
            "NOT_EVALUATED_TECHNICAL_CHEMICAL_INTEGRITY"
            if failed else
            "PASS_LEGACY_A1_CHEMICAL_INTEGRITY_SENSITIVITY"
        ),
        "legacy_audit": str(legacy_audit),
        "legacy_audit_sha256": sha256(legacy_audit),
        "legacy_manifest_mutation": False,
        "denominator_tasks": 18,
        "chemical_integrity_pass_tasks": 18 - len(failed),
        "chemical_integrity_failed_task_indices": failed,
        "per_task": rows,
        "eligible_for_inherited_chain": False,
        "automatic_downstream_action": "NONE",
    }
    audit_path = (
        output_root / "audit"
        / f"nylc_a1_raw_nac_preorganized_calibration_{array_job}"
        "_v2_chemical_integrity.json"
    )
    if audit_path.exists() and BASE.read_json(audit_path) != payload:
        raise FileExistsError(audit_path)
    if not audit_path.exists():
        BASE.write_json(audit_path, payload)
    return audit_path


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", required=True, choices=(
        "describe", "initialize", "prepare", "audit", "merge-if-ready",
        "reclassify-legacy",
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
    elif args.mode == "merge-if-ready":
        merge_if_ready(args.output_root.resolve(), args.array_job)
    else:
        print(reclassify_legacy(args.output_root.resolve(), args.array_job))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
