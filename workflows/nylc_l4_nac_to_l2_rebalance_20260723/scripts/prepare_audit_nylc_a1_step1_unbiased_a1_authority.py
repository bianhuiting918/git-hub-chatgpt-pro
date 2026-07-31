#!/usr/bin/env python3
"""Two-seed unrestrained A1 authority diagnostic and immutable reclassification."""
from __future__ import annotations

import argparse
import importlib.util
import json
import pathlib
import re
import shutil
from typing import Any, Mapping

HERE = pathlib.Path(__file__).resolve().parent
PRE_PATH = HERE / "prepare_audit_nylc_a1_step1_raw_nac_preorganized_calibration.py"
SPEC = importlib.util.spec_from_file_location("_a1_preorg_integrity", PRE_PATH)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError(f"cannot import {PRE_PATH}")
PRE = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(PRE)

ARRAY_TASKS = 2
MPI_RANKS = 8
STAGES = 4
MAXCYC_PER_STAGE = 250
NCYC_PER_STAGE = 100
OUTPUT_ROOT = (
    PRE.TASK_ROOT
    / "a1_activated_nac_20260726/qmmm/a1_step1_unbiased_a1_authority"
)
OLD_ARRAY_JOB = "62500360"
OLD_OUTPUT_ROOT = PRE.OUTPUT_ROOT
SCIENTIFIC_STATUS = "NOT_EVALUATED_TS_COMMITTOR_PMF_BARRIER_MECHANISM"


def task_spec(index: int) -> dict[str, Any]:
    task_index = int(index)
    if task_index not in (0, 1):
        raise ValueError("task index must be 0 or 1")
    seed = (26723, 26737)[task_index]
    return {
        "task_index": task_index,
        "seed": seed,
        "seed_label": f"seed{seed}",
        "source_preorg_task_index": (0, 9)[task_index],
        "mode": "UNRESTRAINED_A1_AUTHORITY_DIAGNOSTIC",
    }


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "unrestrained_minimization_blocks": STAGES,
        "maxcyc_per_block": MAXCYC_PER_STAGE,
        "ncyc_per_block": NCYC_PER_STAGE,
        "reaction_coordinate_restraints": False,
        "position_restraints": False,
        "strict_serial_restart_inheritance": True,
        "original_raw_nac_restart_shas": {
            str(seed): source["start_rst7_sha256"]
            for seed, source in PRE.SOURCES.items()
        },
        "qm_contract": dict(PRE.QM_CONTRACT),
        "automatic_downstream_action": "NONE",
        "scientific_status": SCIENTIFIC_STATUS,
    }


def qmmm_block(qmmask: str) -> str:
    return f"""&qmmm
  qmmask='{qmmask}',
  qmcharge=0,
  spin=1,
  qm_theory='DFTB3',
  dftb_telec=200.0,
  qmshake=0,
/
"""


def unbiased_minimization_input(
    task: Mapping[str, Any], stage_index: int, qmmask: str
) -> str:
    return f"""NylC A1 unrestrained authority task {task['task_index']} block {stage_index}
&cntrl
  imin=1, ntmin=2, maxcyc={MAXCYC_PER_STAGE}, ncyc={NCYC_PER_STAGE}, dx0=0.005,
  ntb=1, cut=10.0, ntpr=5, ntxo=1, ifqnt=1,
  ntr=0, nmropt=0, drms=0.01,
/
{qmmm_block(qmmask)}"""


def initialize(index: int, root: pathlib.Path, commit: str) -> dict[str, Any]:
    task = task_spec(index)
    PRE.initialize(task["source_preorg_task_index"], root, commit)
    manifest = PRE.BASE.read_json(root / "MANIFEST.json")
    source_integrity = PRE.measure_thr267_chemical_integrity(root / "source.rst7")
    manifest.update(
        {
            "schema_version": 1,
            "status": "READY_A1_UNRESTRAINED_AUTHORITY_DIAGNOSTIC",
            "task": task,
            "source_thr267_chemical_integrity": source_integrity,
            "accepted_restart": str(root / "source.rst7"),
            "accepted_restart_sha256": PRE.sha256(root / "source.rst7"),
            "stages": [],
            "terminal": not source_integrity["pass"],
            "terminal_reason": (
                None
                if source_integrity["pass"]
                else "NOT_EVALUATED_SOURCE_THR267_CHEMICAL_INTEGRITY"
            ),
            "reaction_coordinate_restraints": False,
            "position_restraints": False,
            "automatic_downstream_action": "NONE",
            "scientific_status": SCIENTIFIC_STATUS,
        }
    )
    PRE.BASE.write_json(root / "MANIFEST.json", manifest)
    return manifest


def prepare_stage(
    root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path, stage_index: int
) -> dict[str, Any]:
    manifest = PRE.BASE.read_json(root / "MANIFEST.json")
    if manifest.get("terminal"):
        raise ValueError(f"diagnostic already terminal: {manifest['terminal_reason']}")
    if int(stage_index) != len(manifest["stages"]):
        raise ValueError("unrestrained blocks must execute in strict order")
    if output.exists() or scratch.exists():
        raise FileExistsError(output if output.exists() else scratch)
    input_restart = pathlib.Path(manifest["accepted_restart"])
    expected_sha = manifest["accepted_restart_sha256"]
    if not input_restart.is_file() or PRE.sha256(input_restart) != expected_sha:
        raise ValueError("unrestrained block input breaks SHA inheritance")
    output.mkdir(parents=True)
    scratch.mkdir(parents=True)
    text = unbiased_minimization_input(
        manifest["task"], stage_index, manifest["qm_contract"]["qmmask"]
    )
    if any(token in text for token in ("DISANG", "&rst", "restraintmask")):
        raise ValueError("unrestrained authority input contains a restraint")
    (scratch / "stage.in").write_text(text, encoding="utf-8")
    payload = {
        "schema_version": 1,
        "stage_index": int(stage_index),
        "input_restart": str(input_restart),
        "input_restart_sha256": expected_sha,
        "reaction_coordinate_restraints": False,
        "position_restraints": False,
    }
    PRE.BASE.write_json(output / "STAGE_MANIFEST.json", payload)
    return payload


def audit_stage(
    root: pathlib.Path, output: pathlib.Path, scratch: pathlib.Path
) -> dict[str, Any]:
    manifest_path = root / "MANIFEST.json"
    manifest = PRE.BASE.read_json(manifest_path)
    prepared = PRE.BASE.read_json(output / "STAGE_MANIFEST.json")
    stage_out = scratch / "stage.out"
    restart = scratch / "stage.rst7"
    engine_ok, geometry, diagnostics = PRE.BASE.AUTH.TETRA._technical(
        stage_out, restart, manifest
    )
    text = (
        stage_out.read_text(encoding="utf-8", errors="replace")
        if stage_out.is_file()
        else ""
    )
    numerical = PRE.RAW_FORWARD.engine_numerical_health(text)
    engine_technical = bool(engine_ok and numerical["pass"])
    integrity = {"pass": False, "checks": {}, "distances_A": {}}
    if engine_technical:
        integrity = PRE.measure_thr267_chemical_integrity(restart)
    technical = bool(engine_technical and integrity["pass"])
    if not engine_technical:
        classification = "NOT_EVALUATED_TECHNICAL_ENGINE_OR_NUMERICAL"
    elif not integrity["pass"]:
        classification = "NOT_EVALUATED_TECHNICAL_CHEMICAL_INTEGRITY"
    else:
        classification = "PASS_UNRESTRAINED_A1_AUTHORITY_STAGE"

    for name in ("stage.in", "stage.out", "stage.mdinfo", "ENGINE_RC.txt"):
        source = scratch / name
        if source.is_file():
            shutil.copy2(source, output / name)
    if restart.is_file():
        shutil.copy2(restart, output / "stage.rst7")

    diagnostics = dict(diagnostics)
    diagnostics["numerical_health"] = numerical
    result = {
        "schema_version": 1,
        "status": classification,
        "technical_complete": technical,
        "engine_technical_complete": engine_technical,
        "stage_index": prepared["stage_index"],
        "input_restart_sha256": prepared["input_restart_sha256"],
        "output_restart_sha256": (
            PRE.sha256(output / "stage.rst7")
            if (output / "stage.rst7").is_file()
            else None
        ),
        "thr267_chemical_integrity": integrity,
        "geometry": geometry if engine_technical else None,
        "diagnostics": diagnostics,
        "reaction_coordinate_restraints": False,
        "position_restraints": False,
        "scientific_status": SCIENTIFIC_STATUS,
    }
    PRE.BASE.write_json(output / "RESULT.json", result)
    PRE.BASE.write_json(
        output / ("PASS.json" if technical else "NOT_EVALUATED.json"), result
    )
    manifest["stages"].append(result)
    if technical:
        manifest["accepted_restart"] = str(output / "stage.rst7")
        manifest["accepted_restart_sha256"] = result["output_restart_sha256"]
    else:
        manifest["terminal"] = True
        manifest["terminal_reason"] = classification
    if len(manifest["stages"]) >= STAGES:
        manifest["terminal"] = True
        manifest["terminal_reason"] = (
            "COMPLETED_UNRESTRAINED_A1_AUTHORITY_BLOCKS"
            if all(stage["technical_complete"] for stage in manifest["stages"])
            else manifest["terminal_reason"]
        )
    PRE.BASE.write_json(manifest_path, manifest)
    return result


def _write_hashes(root: pathlib.Path) -> None:
    rows = []
    for path in sorted(root.rglob("*")):
        if path.is_file() and path.name != "SHA256.tsv":
            rows.append(f"{PRE.sha256(path)}  {path.relative_to(root)}")
    (root / "SHA256.tsv").write_text("\n".join(rows) + "\n", encoding="utf-8")


def finalize(root: pathlib.Path) -> dict[str, Any]:
    manifest = PRE.BASE.read_json(root / "MANIFEST.json")
    stages = manifest["stages"]
    passed = (
        len(stages) == STAGES
        and all(stage.get("technical_complete") is True for stage in stages)
    )
    result = {
        "schema_version": 1,
        "status": (
            "PASS_UNBIASED_A1_CHEMICAL_INTEGRITY_AUTHORITY"
            if passed
            else "NOT_EVALUATED_TECHNICAL_UNBIASED_A1_AUTHORITY"
        ),
        "technical_complete": passed,
        "seed": manifest["task"]["seed_label"],
        "source_restart_sha256": manifest["source"]["restart_sha256"],
        "source_thr267_chemical_integrity": manifest[
            "source_thr267_chemical_integrity"
        ],
        "stages_attempted": len(stages),
        "stages": stages,
        "strict_sha_inheritance": all(
            stage["input_restart_sha256"]
            == (
                manifest["source"]["restart_sha256"]
                if index == 0
                else stages[index - 1]["output_restart_sha256"]
            )
            for index, stage in enumerate(stages)
        ),
        "reaction_coordinate_restraints": False,
        "position_restraints": False,
        "automatic_downstream_action": "NONE",
        "scientific_status": SCIENTIFIC_STATUS,
    }
    PRE.BASE.write_json(root / "RESULT.json", result)
    PRE.BASE.write_json(
        root / ("PASS.json" if passed else "NOT_EVALUATED.json"), result
    )
    _write_hashes(root)
    return result


def reclassification_status(integrity_passes: list[bool]) -> str:
    return (
        "PASS_ALL_TASKS_THR267_CHEMICAL_INTEGRITY"
        if integrity_passes and all(integrity_passes)
        else "NOT_EVALUATED_TECHNICAL_CHEMICAL_INTEGRITY"
    )


def reclassify_existing(array_job: str = OLD_ARRAY_JOB) -> dict[str, Any]:
    rows = []
    for index in range(PRE.ARRAY_TASKS):
        root = OLD_OUTPUT_ROOT / f"attempt_{array_job}_{index}"
        old_result = PRE.BASE.read_json(root / "RESULT.json")
        restart = root / "stage.rst7"
        integrity = (
            PRE.measure_thr267_chemical_integrity(restart)
            if restart.is_file()
            else {"pass": False, "checks": {}, "distances_A": {}}
        )
        rows.append(
            {
                "task_index": index,
                "seed": old_result["task"]["seed"],
                "old_status": old_result["status"],
                "old_eligible_for_inherited_chain": old_result.get(
                    "eligible_for_inherited_chain"
                ),
                "stage_restart_sha256": (
                    PRE.sha256(restart) if restart.is_file() else None
                ),
                "thr267_chemical_integrity": integrity,
                "reclassified_status": (
                    "UNCHANGED_CHEMICAL_INTEGRITY_PASS"
                    if integrity["pass"]
                    else "NOT_EVALUATED_TECHNICAL_CHEMICAL_INTEGRITY"
                ),
            }
        )
    passes = [row["thr267_chemical_integrity"]["pass"] for row in rows]
    payload = {
        "schema_version": 1,
        "status": reclassification_status(passes),
        "source_array_job": array_job,
        "denominator_tasks": len(rows),
        "chemical_integrity_pass_tasks": sum(passes),
        "chemical_integrity_fail_tasks": len(rows) - sum(passes),
        "old_manifests_untouched": True,
        "old_restart_inheritance_revoked_on_failure": True,
        "per_task": rows,
        "scientific_status": SCIENTIFIC_STATUS,
    }
    audit = (
        OLD_OUTPUT_ROOT
        / "audit"
        / f"nylc_a1_raw_nac_preorganized_calibration_{array_job}_chemical_integrity_v1.json"
    )
    audit.parent.mkdir(parents=True, exist_ok=True)
    if audit.exists():
        if PRE.BASE.read_json(audit) != payload:
            raise FileExistsError(audit)
    else:
        PRE.BASE.write_json(audit, payload)
    return payload


def merge_if_ready(output_root: pathlib.Path, array_job: str) -> bool:
    paths = [
        output_root / f"attempt_{array_job}_{index}" / "RESULT.json"
        for index in range(ARRAY_TASKS)
    ]
    if not all(path.is_file() for path in paths):
        return False
    results = [PRE.BASE.read_json(path) for path in paths]
    technical = sum(result.get("technical_complete") is True for result in results)
    payload = {
        "schema_version": 1,
        "status": (
            "PASS_BOTH_SEEDS_UNBIASED_A1_CHEMICAL_INTEGRITY_AUTHORITY"
            if technical == ARRAY_TASKS
            else "NOT_EVALUATED_TECHNICAL_UNBIASED_A1_AUTHORITY"
        ),
        "denominator_seeds": ARRAY_TASKS,
        "technical_pass_seeds": technical,
        "per_seed": results,
        "automatic_downstream_action": "NONE",
        "scientific_status": SCIENTIFIC_STATUS,
    }
    audit = (
        output_root / "audit"
        / f"nylc_a1_unbiased_a1_authority_{array_job}.json"
    )
    audit.parent.mkdir(parents=True, exist_ok=True)
    if audit.exists():
        if PRE.BASE.read_json(audit) != payload:
            raise FileExistsError(audit)
    else:
        PRE.BASE.write_json(audit, payload)
    return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        required=True,
        choices=(
            "describe", "initialize", "prepare-stage", "audit-stage",
            "finalize", "reclassify-existing", "merge-if-ready",
        ),
    )
    parser.add_argument("--task-index", type=int)
    parser.add_argument("--stage-index", type=int)
    parser.add_argument("--root", type=pathlib.Path)
    parser.add_argument("--output", type=pathlib.Path)
    parser.add_argument("--scratch", type=pathlib.Path)
    parser.add_argument("--github-commit", default="unknown")
    parser.add_argument("--array-job", default=OLD_ARRAY_JOB)
    parser.add_argument("--output-root", type=pathlib.Path, default=OUTPUT_ROOT)
    args = parser.parse_args()
    if args.mode == "describe":
        print(json.dumps(describe(), sort_keys=True))
    elif args.mode == "initialize":
        initialize(args.task_index, args.root.resolve(), args.github_commit)
    elif args.mode == "prepare-stage":
        prepare_stage(
            args.root.resolve(), args.output.resolve(), args.scratch.resolve(),
            args.stage_index,
        )
    elif args.mode == "audit-stage":
        audit_stage(args.root.resolve(), args.output.resolve(), args.scratch.resolve())
    elif args.mode == "finalize":
        finalize(args.root.resolve())
    elif args.mode == "reclassify-existing":
        print(json.dumps(reclassify_existing(args.array_job), sort_keys=True))
    else:
        merge_if_ready(args.output_root.resolve(), args.array_job)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
