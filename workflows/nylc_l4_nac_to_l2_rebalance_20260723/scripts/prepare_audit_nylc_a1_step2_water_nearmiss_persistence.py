#!/usr/bin/env python3
"""Unrestrained persistence sampling from two Step2 water near-miss frames."""
from __future__ import annotations

import argparse
import fcntl
import importlib.util
import json
import pathlib
from typing import Any, Mapping

HERE = pathlib.Path(__file__).resolve().parent
BASE_PATH = HERE / "prepare_audit_nylc_a1_step2_water_reorganization_sampling.py"
_SPEC = importlib.util.spec_from_file_location("_nylc_step2_water_reorg_base", BASE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"cannot load {BASE_PATH}")
BASE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(BASE)

ARRAY_TASKS = 8
REPLICAS_PER_SEED = 4
COMBINED_ARRAY_TASKS = 9
MPI_RANKS = 8
VELOCITY_SEEDS = {
    "seed26723": (26723101, 26723102, 26723103, 26723104),
    "seed26737": (26737101, 26737102, 26737103, 26737104),
}
SOURCES = {
    "seed26723": {
        "seed_index": 0,
        "parent_task": 3,
        "restart": BASE.TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_step2_water_reorganization_sampling/attempt_62300430_3/selected_near_miss_0.rst7",
        "restart_sha256": "a65834ae6c0f82291fd0b8c38ece41983b501b2444deb98de1794b92d79622a7",
    },
    "seed26737": {
        "seed_index": 1,
        "parent_task": 7,
        "restart": BASE.TASK_ROOT / "a1_activated_nac_20260726/qmmm/a1_step2_water_reorganization_sampling/attempt_62300430_7/selected_near_miss_0.rst7",
        "restart_sha256": "34ad5af01f5a2cd7e7ce0342590ab66f201902a31fd49236546441afa82f5ea6",
    },
}
SCOPE = "STEP2_WATER_NEARMISS_PERSISTENCE_ONLY_NOT_PRODUCT_TS_PATH_PMF_BARRIER_OR_MECHANISM"
NEXT = "PROMOTE_TO_149QM_A2_ONLY_AFTER_CROSS_SEED_PERSISTENT_WATER_PASS"


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("persistence task index must be 0..7")
    seed = "seed26723" if index < REPLICAS_PER_SEED else "seed26737"
    replica = index % REPLICAS_PER_SEED
    return {
        "task_index": index,
        "seed_index": SOURCES[seed]["seed_index"],
        "replica": replica,
        "seed": seed,
        "velocity_seed": VELOCITY_SEEDS[seed][replica],
        "source": SOURCES[seed],
    }


def sampling_input(seed: str, replica: int, velocity_seed: int, qmmask: str) -> str:
    return BASE.sampling_input(seed, replica, velocity_seed, qmmask)


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "seed_count": 2,
        "replicas_per_seed": REPLICAS_PER_SEED,
        "persistence_task_count": ARRAY_TASKS,
        "combined_array_task_count": COMBINED_ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "steps": BASE.MD_STEPS,
        "dt_ps": BASE.MD_DT_PS,
        "duration_ps": BASE.MD_STEPS * BASE.MD_DT_PS,
        "expected_frames": BASE.EXPECTED_FRAMES,
        "minimum_consecutive_frames": BASE.MIN_CONSECUTIVE_FRAMES,
        "reactive_restraints": False,
        "water_position_restraints": False,
        "scan_all_complete_waters": True,
        "step1_qm_contract": BASE.STEP1_QM_CONTRACT,
        "source_restart_sha256": {
            seed: source["restart_sha256"] for seed, source in SOURCES.items()
        },
        "velocity_seeds": VELOCITY_SEEDS,
        "automatic_downstream_action": "NONE",
    }


def initialize(
    task_index: int,
    output: pathlib.Path,
    scratch: pathlib.Path,
    code_root: pathlib.Path,
    github_commit: str,
) -> dict[str, Any]:
    task = task_spec(task_index)
    source = task["source"]
    restart = pathlib.Path(source["restart"])
    if not restart.is_file():
        raise FileNotFoundError(restart)
    if BASE.sha256(restart) != source["restart_sha256"]:
        raise ValueError("near-miss restart SHA256 mismatch")
    base_task_index = 0 if task["seed_index"] == 0 else 4
    manifest = BASE.initialize(
        base_task_index, output, scratch, code_root, github_commit
    )

    import parmed as pmd

    structure = pmd.load_file(str(BASE.PRMTOP), xyz=str(restart))
    box = [float(value) for value in structure.box[:6]]
    BASE.S2._periodic_cell(box)
    (scratch / "stage.in").write_text(
        sampling_input(
            task["seed"], task["replica"], task["velocity_seed"], manifest["qmmask"]
        ),
        encoding="utf-8",
    )
    manifest.update(
        {
            "status": "READY_STEP2_WATER_NEARMISS_PERSISTENCE",
            "scientific_scope": SCOPE,
            "task_index": task["task_index"],
            "seed_index": task["seed_index"],
            "replica": task["replica"],
            "seed": task["seed"],
            "velocity_seed": task["velocity_seed"],
            "box": box,
            "source": {
                "parent_job": "62300430",
                "parent_task": source["parent_task"],
                "restart": str(restart),
                "restart_sha256": source["restart_sha256"],
                "selection": "best_persisted_near_miss_per_seed",
            },
            "protocol": describe(),
            "automatic_downstream_action": "NONE",
            "NEXT": NEXT,
        }
    )
    BASE.write_json(output / "SOURCE_MANIFEST.json", manifest)
    return manifest


def audit(output: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    result = BASE.audit(output, scratch)
    manifest = json.loads(
        (output / "SOURCE_MANIFEST.json").read_text(encoding="utf-8")
    )
    technical = result.get("technical_complete") is True
    result.update(
        {
            "status": (
                "PASS_TECHNICAL_STEP2_WATER_NEARMISS_PERSISTENCE"
                if technical
                else "NOT_EVALUATED_TECHNICAL_STEP2_WATER_NEARMISS_PERSISTENCE"
            ),
            "scientific_scope": SCOPE,
            "parent_near_miss": manifest["source"],
            "automatic_downstream_action": "NONE",
            "NEXT": NEXT,
        }
    )
    BASE.write_json(output / "RESULT.json", result)
    BASE.write_json(
        output / ("PASS.json" if technical else "NOT_EVALUATED.json"), result
    )
    BASE._write_hashes(output)
    return result


def merge_if_ready(output_root: pathlib.Path, array_job: str) -> bool:
    audit_path = (
        output_root
        / "audit"
        / f"nylc_a1_step2_water_nearmiss_persistence_{array_job}.json"
    )
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    lock = audit_path.with_suffix(".json.lock")
    with lock.open("a+") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        paths = [
            output_root / f"attempt_{array_job}_{index}" / "RESULT.json"
            for index in range(ARRAY_TASKS)
        ]
        if not all(path.is_file() for path in paths):
            return False
        results = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
        technical_count = sum(
            item.get("technical_complete") is True for item in results
        )
        selected_events: dict[str, list[dict[str, Any]]] = {}
        selected_near_misses: dict[str, list[dict[str, Any]]] = {}
        for seed in ("seed26723", "seed26737"):
            events = [
                dict(event)
                for result in results
                if result["seed"] == seed
                for event in result.get("events", [])
            ]
            selected_events[seed] = sorted(
                events,
                key=lambda item: (
                    item["score"],
                    item["task_index"],
                    item["start_frame"],
                ),
            )[: BASE.MAX_EVENTS_PER_SEED]
            near = [
                dict(row)
                for result in results
                if result["seed"] == seed
                for row in result.get("near_misses", [])
            ]
            selected_near_misses[seed] = sorted(
                near,
                key=lambda row: (
                    row["near_miss_score"],
                    row["task_index"],
                    row["frame"],
                ),
            )[: BASE.MAX_NEAR_MISSES_PER_SEED]
        reproduced = (
            technical_count == ARRAY_TASKS
            and all(selected_events[seed] for seed in selected_events)
        )
        status = (
            "NOT_EVALUATED_TECHNICAL_STEP2_WATER_NEARMISS_PERSISTENCE"
            if technical_count != ARRAY_TASKS
            else "PASS_STEP2_PERSISTENT_WATER_REPRODUCED"
            if reproduced
            else "FAIL_NO_REPRODUCED_STEP2_PERSISTENT_WATER"
        )
        payload = {
            "schema_version": 1,
            "status": status,
            "technical_denominator": f"{technical_count}/{ARRAY_TASKS}",
            "per_seed_selected_events": selected_events,
            "per_seed_selected_near_misses": selected_near_misses,
            "per_task": sorted(results, key=lambda item: item["task_index"]),
            "automatic_downstream_action": "NONE",
            "NEXT": NEXT,
        }
        if audit_path.exists():
            if json.loads(audit_path.read_text(encoding="utf-8")) != payload:
                raise FileExistsError(f"existing audit differs: {audit_path}")
        else:
            BASE.write_json(audit_path, payload)
        return True


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode",
        required=True,
        choices=("describe", "initialize", "audit", "merge-if-ready"),
    )
    parser.add_argument("--task-index", type=int)
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
        initialize(
            args.task_index,
            args.output.resolve(),
            args.scratch.resolve(),
            args.code_root.resolve(),
            args.github_commit,
        )
    elif args.mode == "audit":
        audit(args.output.resolve(), args.scratch.resolve())
    else:
        merge_if_ready(args.output_root.resolve(), args.array_job)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
