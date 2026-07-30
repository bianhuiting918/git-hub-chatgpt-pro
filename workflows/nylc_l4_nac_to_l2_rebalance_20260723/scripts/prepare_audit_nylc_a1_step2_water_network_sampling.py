#!/usr/bin/env python3
"""Unrestrained 5 ps Step2 water-network sampling from stable acyl endpoints.

All waters remain MM. Direct-water and two-water-relay geometries are
recorded separately. This protocol generates A2 candidates only; it does not
construct a Step2 product, TS, committor, PMF, barrier, or mechanism claim.
"""
from __future__ import annotations

import argparse
import fcntl
import importlib.util
import json
import math
import pathlib
import re
import shutil
from typing import Any, Mapping, Sequence

HERE = pathlib.Path(__file__).resolve().parent
BASE_PATH = HERE / "prepare_audit_nylc_a1_step2_water_reorganization_sampling.py"
_SPEC = importlib.util.spec_from_file_location("_nylc_step2_water_reorg_base", BASE_PATH)
if _SPEC is None or _SPEC.loader is None:
    raise RuntimeError(f"cannot load {BASE_PATH}")
BASE = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(BASE)

ARRAY_TASKS = 8
REPLICAS_PER_SEED = 4
MPI_RANKS = 8
MD_STEPS = 10000
MD_DT_PS = 0.0005
MD_NTWX = 20
EXPECTED_FRAMES = MD_STEPS // MD_NTWX
MIN_CONSECUTIVE_FRAMES = 5
LOCAL_WATER_RADIUS_A = 6.0
MAX_SELECTED_PER_TASK = 1
MAX_SELECTED_PER_SEED = 2
VELOCITY_SEEDS = {
    "seed26723": (26723201, 26723202, 26723203, 26723204),
    "seed26737": (26737201, 26737202, 26737203, 26737204),
}
DIRECT_FILTER = dict(BASE.HIT_FILTER)
DIRECT_NEAR_FILTER = dict(BASE.NEAR_MISS_FILTER)
RELAY_FILTER = {
    "c12_ow_A": (2.70, 3.50),
    "o2_c12_ow_deg": (95.0, 125.0),
    "attack_h_relay_o_A_max": 2.50,
    "attack_ow_attack_h_relay_o_deg_min": 130.0,
    "relay_h_nalpha_A_max": 2.50,
    "relay_ow_relay_h_nalpha_deg_min": 130.0,
}
RELAY_NEAR_FILTER = {
    "c12_ow_A": (2.50, 4.50),
    "o2_c12_ow_deg": (80.0, 145.0),
    "attack_h_relay_o_A_max": 3.20,
    "attack_ow_attack_h_relay_o_deg_min": 100.0,
    "relay_h_nalpha_A_max": 3.50,
    "relay_ow_relay_h_nalpha_deg_min": 100.0,
}
SCOPE = "STEP2_WATER_NETWORK_SAMPLING_ONLY_NOT_PRODUCT_TS_PATH_PMF_BARRIER_OR_MECHANISM"
NEXT = "PROMOTE_ROUTE_SPECIFIC_CROSS_SEED_WATER_NETWORK_CANDIDATES_TO_FINAL_QM_A2_PREFLIGHT"


def task_spec(task_index: int) -> dict[str, Any]:
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..7")
    seed_index = index // REPLICAS_PER_SEED
    replica = index % REPLICAS_PER_SEED
    source = BASE.S2.source_for_index(seed_index)
    return {
        "task_index": index,
        "seed_index": seed_index,
        "replica": replica,
        "seed": source["seed"],
        "velocity_seed": VELOCITY_SEEDS[source["seed"]][replica],
        "source": source,
    }


def sampling_input(seed: str, replica: int, velocity_seed: int, qmmask: str) -> str:
    return f"""NylC A1 Step2 water network {seed} replica{replica}
&cntrl
  imin=0, irest=0, ntx=1, nstlim={MD_STEPS}, dt={MD_DT_PS},
  ntb=1, cut=10.0, ntt=3, gamma_ln=2.0,
  temp0=300.0, tempi=300.0, ig={velocity_seed},
  ntc=1, ntf=1, ntwx={MD_NTWX}, ntpr=100, ntwr={MD_STEPS},
  ioutfm=1, ntxo=1, ifqnt=1, ntr=0, nmropt=0,
/
{BASE.S2.qmmm_block(qmmask)}
"""


def describe() -> dict[str, Any]:
    return {
        "schema_version": 1,
        "seed_count": 2,
        "replicas_per_seed": REPLICAS_PER_SEED,
        "array_task_count": ARRAY_TASKS,
        "mpi_ranks_per_task": MPI_RANKS,
        "steps": MD_STEPS,
        "dt_ps": MD_DT_PS,
        "duration_ps": MD_STEPS * MD_DT_PS,
        "trajectory_stride_steps": MD_NTWX,
        "expected_frames": EXPECTED_FRAMES,
        "minimum_consecutive_frames": MIN_CONSECUTIVE_FRAMES,
        "local_water_radius_A": LOCAL_WATER_RADIUS_A,
        "step1_qm_contract": dict(BASE.STEP1_QM_CONTRACT),
        "scan_all_complete_waters": True,
        "scan_direct_and_two_water_relay": True,
        "direct_filter": dict(DIRECT_FILTER),
        "relay_filter": dict(RELAY_FILTER),
        "water_identity_restraints": False,
        "water_position_restraints": False,
        "reactive_restraints": False,
        "automatic_downstream_action": "NONE",
        "candidate_is_not_ts": True,
    }


def direct_water_candidate(row: Mapping[str, Any]) -> bool:
    return bool(
        DIRECT_FILTER["c12_ow_A"][0] <= float(row["c12_ow_A"]) <= DIRECT_FILTER["c12_ow_A"][1]
        and DIRECT_FILTER["o2_c12_ow_deg"][0] <= float(row["o2_c12_ow_deg"]) <= DIRECT_FILTER["o2_c12_ow_deg"][1]
        and float(row["h_nalpha_A"]) <= DIRECT_FILTER["h_nalpha_A_max"]
        and float(row["ow_h_nalpha_deg"]) >= DIRECT_FILTER["ow_h_nalpha_deg_min"]
        and row.get("acyl_guard_pass") is True
        and row.get("proton_guard_pass") is True
    )


def direct_near_candidate(row: Mapping[str, Any]) -> bool:
    return bool(
        DIRECT_NEAR_FILTER["c12_ow_A"][0] <= float(row["c12_ow_A"]) <= DIRECT_NEAR_FILTER["c12_ow_A"][1]
        and DIRECT_NEAR_FILTER["o2_c12_ow_deg"][0] <= float(row["o2_c12_ow_deg"]) <= DIRECT_NEAR_FILTER["o2_c12_ow_deg"][1]
        and float(row["h_nalpha_A"]) <= DIRECT_NEAR_FILTER["h_nalpha_A_max"]
        and float(row["ow_h_nalpha_deg"]) >= DIRECT_NEAR_FILTER["ow_h_nalpha_deg_min"]
        and row.get("acyl_guard_pass") is True
        and row.get("proton_guard_pass") is True
    )


def relay_water_candidate(row: Mapping[str, Any]) -> bool:
    return bool(
        RELAY_FILTER["c12_ow_A"][0] <= float(row["c12_ow_A"]) <= RELAY_FILTER["c12_ow_A"][1]
        and RELAY_FILTER["o2_c12_ow_deg"][0] <= float(row["o2_c12_ow_deg"]) <= RELAY_FILTER["o2_c12_ow_deg"][1]
        and float(row["attack_h_relay_o_A"]) <= RELAY_FILTER["attack_h_relay_o_A_max"]
        and float(row["attack_ow_attack_h_relay_o_deg"]) >= RELAY_FILTER["attack_ow_attack_h_relay_o_deg_min"]
        and float(row["relay_h_nalpha_A"]) <= RELAY_FILTER["relay_h_nalpha_A_max"]
        and float(row["relay_ow_relay_h_nalpha_deg"]) >= RELAY_FILTER["relay_ow_relay_h_nalpha_deg_min"]
        and row.get("acyl_guard_pass") is True
        and row.get("proton_guard_pass") is True
    )


def relay_near_candidate(row: Mapping[str, Any]) -> bool:
    return bool(
        RELAY_NEAR_FILTER["c12_ow_A"][0] <= float(row["c12_ow_A"]) <= RELAY_NEAR_FILTER["c12_ow_A"][1]
        and RELAY_NEAR_FILTER["o2_c12_ow_deg"][0] <= float(row["o2_c12_ow_deg"]) <= RELAY_NEAR_FILTER["o2_c12_ow_deg"][1]
        and float(row["attack_h_relay_o_A"]) <= RELAY_NEAR_FILTER["attack_h_relay_o_A_max"]
        and float(row["attack_ow_attack_h_relay_o_deg"]) >= RELAY_NEAR_FILTER["attack_ow_attack_h_relay_o_deg_min"]
        and float(row["relay_h_nalpha_A"]) <= RELAY_NEAR_FILTER["relay_h_nalpha_A_max"]
        and float(row["relay_ow_relay_h_nalpha_deg"]) >= RELAY_NEAR_FILTER["relay_ow_relay_h_nalpha_deg_min"]
        and row.get("acyl_guard_pass") is True
        and row.get("proton_guard_pass") is True
    )


def _direct_score(row: Mapping[str, Any]) -> float:
    return (
        ((float(row["c12_ow_A"]) - 3.05) / 0.35) ** 2
        + ((float(row["o2_c12_ow_deg"]) - 110.0) / 15.0) ** 2
        + ((float(row["h_nalpha_A"]) - 1.90) / 0.60) ** 2
        + ((float(row["ow_h_nalpha_deg"]) - 180.0) / 50.0) ** 2
    )


def _relay_score(row: Mapping[str, Any]) -> float:
    return (
        ((float(row["c12_ow_A"]) - 3.05) / 0.35) ** 2
        + ((float(row["o2_c12_ow_deg"]) - 110.0) / 15.0) ** 2
        + ((float(row["attack_h_relay_o_A"]) - 1.85) / 0.65) ** 2
        + ((float(row["attack_ow_attack_h_relay_o_deg"]) - 180.0) / 50.0) ** 2
        + ((float(row["relay_h_nalpha_A"]) - 1.90) / 0.60) ** 2
        + ((float(row["relay_ow_relay_h_nalpha_deg"]) - 180.0) / 50.0) ** 2
    )


def _event_key(row: Mapping[str, Any], route: str) -> tuple[int, ...]:
    if route == "direct":
        return (
            int(row["replica"]),
            int(row["attack_water_oxygen_index1"]),
            int(row["attack_donor_h_index1"]),
        )
    return (
        int(row["replica"]),
        int(row["attack_water_oxygen_index1"]),
        int(row["attack_donor_h_index1"]),
        int(row["relay_water_oxygen_index1"]),
        int(row["relay_donor_h_index1"]),
    )


def _collapse_events(rows: Sequence[Mapping[str, Any]], route: str) -> list[dict[str, Any]]:
    grouped: dict[tuple[int, ...], list[dict[str, Any]]] = {}
    for source in rows:
        row = dict(source)
        passed = direct_water_candidate(row) if route == "direct" else relay_water_candidate(row)
        if passed:
            grouped.setdefault(_event_key(row, route), []).append(row)
    events: list[dict[str, Any]] = []
    score = _direct_score if route == "direct" else _relay_score
    for key, values in grouped.items():
        ordered = sorted(values, key=lambda item: int(item["frame"]))
        run: list[dict[str, Any]] = []
        for row in ordered:
            if run and int(row["frame"]) != int(run[-1]["frame"]) + 1:
                if len(run) >= MIN_CONSECUTIVE_FRAMES:
                    events.append(_event_from_run(route, key, run, score))
                run = []
            run.append(row)
        if len(run) >= MIN_CONSECUTIVE_FRAMES:
            events.append(_event_from_run(route, key, run, score))
    return sorted(events, key=lambda item: (item["score"], item["replica"], item["start_frame"]))


def _event_from_run(route: str, key: tuple[int, ...], run: Sequence[Mapping[str, Any]], score_fn: Any) -> dict[str, Any]:
    representative = min(run, key=lambda row: (score_fn(row), int(row["frame"])))
    payload = {
        "route": route,
        "replica": key[0],
        "start_frame": int(run[0]["frame"]),
        "end_frame": int(run[-1]["frame"]),
        "consecutive_frames": len(run),
        "representative_frame": int(representative["frame"]),
        "score": score_fn(representative),
        "representative_geometry": dict(representative),
    }
    return payload


def _select_near(rows: Sequence[Mapping[str, Any]], route: str, limit: int = 4) -> list[dict[str, Any]]:
    score_fn = _direct_score if route == "direct" else _relay_score
    predicate = direct_near_candidate if route == "direct" else relay_near_candidate
    best: dict[tuple[int, ...], dict[str, Any]] = {}
    for source in rows:
        if not predicate(source):
            continue
        row = dict(source)
        row["route"] = route
        row["near_miss_score"] = score_fn(row)
        identity = (
            (int(row["attack_water_oxygen_index1"]),)
            if route == "direct"
            else (int(row["attack_water_oxygen_index1"]), int(row["relay_water_oxygen_index1"]))
        )
        incumbent = best.get(identity)
        if incumbent is None or (row["near_miss_score"], int(row["frame"])) < (
            incumbent["near_miss_score"], int(incumbent["frame"])
        ):
            best[identity] = row
    return sorted(best.values(), key=lambda row: (row["near_miss_score"], int(row["frame"])))[: int(limit)]


def _scan_frames_network(trajectory: pathlib.Path, manifest: Mapping[str, Any]):
    import parmed as pmd

    topology = pmd.load_file(str(BASE.PRMTOP))
    frames, warnings = BASE.S2._read_md_frames(trajectory, len(topology.atoms))
    waters = BASE._complete_water_indices(topology)
    cell = BASE.S2._periodic_cell(manifest["box"])
    heavy = [int(value) for value in manifest["qm_heavy_atom_indices"]]
    direct_rows: list[dict[str, Any]] = []
    relay_rows: list[dict[str, Any]] = []
    guards: list[dict[str, Any]] = []
    for frame_index, frame in enumerate(frames):
        coords = {index + 1: frame[index] for index in range(len(frame))}
        acyl_geometry = BASE.AC.geometry_from_coordinates(coords, heavy)
        acyl_guard, proton_guard = BASE._guards(acyl_geometry)
        guards.append({"frame": frame_index, "acyl_guard_pass": acyl_guard, "proton_guard_pass": proton_guard})
        if not (acyl_guard and proton_guard):
            continue
        c12 = coords[BASE.REACTIVE["c12"]]
        o2 = coords[BASE.REACTIVE["o2"]]
        nalpha = coords[BASE.REACTIVE["nalpha"]]
        local = []
        for oxygen_index1, hydrogen_indices1 in waters:
            c12_ow = BASE.S2._distance(c12, coords[oxygen_index1], cell)
            if c12_ow <= LOCAL_WATER_RADIUS_A:
                local.append((oxygen_index1, hydrogen_indices1, c12_ow))
        for attack_o, attack_hs, c12_ow in local:
            if not DIRECT_NEAR_FILTER["c12_ow_A"][0] <= c12_ow <= DIRECT_NEAR_FILTER["c12_ow_A"][1]:
                continue
            attack_angle = BASE.S2._angle(o2, c12, coords[attack_o], cell)
            if not DIRECT_NEAR_FILTER["o2_c12_ow_deg"][0] <= attack_angle <= DIRECT_NEAR_FILTER["o2_c12_ow_deg"][1]:
                continue
            direct_options = []
            for attack_h in attack_hs:
                h_nalpha = BASE.S2._distance(coords[attack_h], nalpha, cell)
                h_angle = BASE.S2._angle(coords[attack_o], coords[attack_h], nalpha, cell)
                direct_options.append((h_nalpha, -h_angle, attack_h, h_angle))
            h_nalpha, _, direct_h, h_angle = min(direct_options)
            direct_rows.append({
                "frame": frame_index,
                "time_ps": (frame_index + 1) * MD_NTWX * MD_DT_PS,
                "replica": manifest["replica"],
                "attack_water_oxygen_index1": attack_o,
                "attack_water_hydrogen_indices1": list(attack_hs),
                "attack_donor_h_index1": direct_h,
                "c12_ow_A": c12_ow,
                "o2_c12_ow_deg": attack_angle,
                "h_nalpha_A": h_nalpha,
                "ow_h_nalpha_deg": h_angle,
                "acyl_guard_pass": acyl_guard,
                "proton_guard_pass": proton_guard,
            })
            relay_options = []
            for attack_h in attack_hs:
                for relay_o, relay_hs, _ in local:
                    if relay_o == attack_o:
                        continue
                    attack_h_relay_o = BASE.S2._distance(coords[attack_h], coords[relay_o], cell)
                    if attack_h_relay_o > RELAY_NEAR_FILTER["attack_h_relay_o_A_max"]:
                        continue
                    bridge_angle = BASE.S2._angle(coords[attack_o], coords[attack_h], coords[relay_o], cell)
                    if bridge_angle < RELAY_NEAR_FILTER["attack_ow_attack_h_relay_o_deg_min"]:
                        continue
                    for relay_h in relay_hs:
                        relay_h_nalpha = BASE.S2._distance(coords[relay_h], nalpha, cell)
                        relay_angle = BASE.S2._angle(coords[relay_o], coords[relay_h], nalpha, cell)
                        row = {
                            "frame": frame_index,
                            "time_ps": (frame_index + 1) * MD_NTWX * MD_DT_PS,
                            "replica": manifest["replica"],
                            "attack_water_oxygen_index1": attack_o,
                            "attack_water_hydrogen_indices1": list(attack_hs),
                            "attack_donor_h_index1": attack_h,
                            "relay_water_oxygen_index1": relay_o,
                            "relay_water_hydrogen_indices1": list(relay_hs),
                            "relay_donor_h_index1": relay_h,
                            "c12_ow_A": c12_ow,
                            "o2_c12_ow_deg": attack_angle,
                            "attack_h_relay_o_A": attack_h_relay_o,
                            "attack_ow_attack_h_relay_o_deg": bridge_angle,
                            "relay_h_nalpha_A": relay_h_nalpha,
                            "relay_ow_relay_h_nalpha_deg": relay_angle,
                            "acyl_guard_pass": acyl_guard,
                            "proton_guard_pass": proton_guard,
                        }
                        if relay_near_candidate(row):
                            relay_options.append(row)
            if relay_options:
                relay_rows.append(min(relay_options, key=_relay_score))
    return frames, direct_rows, relay_rows, guards, warnings


def initialize(task_index: int, output: pathlib.Path, scratch: pathlib.Path, code_root: pathlib.Path, github_commit: str) -> dict[str, Any]:
    task = task_spec(task_index)
    base_index = 0 if task["seed_index"] == 0 else 4
    manifest = BASE.initialize(base_index, output, scratch, code_root, github_commit)
    (scratch / "stage.in").write_text(
        sampling_input(task["seed"], task["replica"], task["velocity_seed"], manifest["qmmask"]),
        encoding="utf-8",
    )
    manifest.update({
        "status": "READY_STEP2_WATER_NETWORK_SAMPLING",
        "scientific_scope": SCOPE,
        "task_index": task["task_index"],
        "seed_index": task["seed_index"],
        "replica": task["replica"],
        "seed": task["seed"],
        "velocity_seed": task["velocity_seed"],
        "protocol": describe(),
        "automatic_downstream_action": "NONE",
        "NEXT": NEXT,
    })
    BASE.write_json(output / "SOURCE_MANIFEST.json", manifest)
    return manifest


def _persist_candidate(output: pathlib.Path, frames: Sequence[Any], box: Sequence[float], item: dict[str, Any], prefix: str) -> None:
    frame_index = int(item.get("representative_frame", item.get("frame", 0)))
    destination = output / f"{prefix}.rst7"
    BASE._write_selected_restart(destination, frames[frame_index], box)
    item["selected_restart"] = destination.name
    item["selected_restart_sha256"] = BASE.sha256(destination)


def audit(output: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest = json.loads((output / "SOURCE_MANIFEST.json").read_text(encoding="utf-8"))
    stage_out = scratch / "stage.out"
    trajectory = scratch / "stage.nc"
    restart = scratch / "stage.rst7"
    engine_rc = int((scratch / "engine.rc").read_text().strip()) if (scratch / "engine.rc").is_file() else 999
    text = stage_out.read_text(encoding="utf-8", errors="replace") if stage_out.is_file() else ""
    banner = BASE._parse_engine_contract(text)
    hard = {name: len(re.findall(pattern, text, re.I)) for name, pattern in BASE.HARD_PATTERNS.items()}
    frames: list[Any] = []
    direct_rows: list[dict[str, Any]] = []
    relay_rows: list[dict[str, Any]] = []
    guards: list[dict[str, Any]] = []
    warnings: list[str] = []
    parse_error = None
    if trajectory.is_file() and trajectory.stat().st_size:
        try:
            frames, direct_rows, relay_rows, guards, warnings = _scan_frames_network(trajectory, manifest)
        except Exception as error:
            parse_error = f"{type(error).__name__}: {error}"
    technical = bool(
        engine_rc == 0
        and BASE._last_nstep(text) == MD_STEPS
        and restart.is_file() and restart.stat().st_size > 0
        and trajectory.is_file() and trajectory.stat().st_size > 0
        and len(frames) == EXPECTED_FRAMES
        and BASE._engine_contract_pass(banner)
        and sum(hard.values()) == 0
        and not warnings
        and parse_error is None
    )
    direct_hits = [row for row in direct_rows if direct_water_candidate(row)]
    relay_hits = [row for row in relay_rows if relay_water_candidate(row)]
    direct_events = _collapse_events(direct_hits, "direct") if technical else []
    relay_events = _collapse_events(relay_hits, "relay") if technical else []
    direct_near = _select_near(direct_rows, "direct") if technical else []
    relay_near = _select_near(relay_rows, "relay") if technical else []
    for rank, item in enumerate(direct_events[:MAX_SELECTED_PER_TASK]):
        _persist_candidate(output, frames, manifest["box"], item, f"direct_event_{rank}")
    for rank, item in enumerate(relay_events[:MAX_SELECTED_PER_TASK]):
        _persist_candidate(output, frames, manifest["box"], item, f"relay_event_{rank}")
    if not direct_events:
        for rank, item in enumerate(direct_near[:MAX_SELECTED_PER_TASK]):
            _persist_candidate(output, frames, manifest["box"], item, f"direct_near_miss_{rank}")
    if not relay_events:
        for rank, item in enumerate(relay_near[:MAX_SELECTED_PER_TASK]):
            _persist_candidate(output, frames, manifest["box"], item, f"relay_near_miss_{rank}")
    for name, rows in (
        ("DIRECT_WATER_EVENTS.jsonl", direct_events),
        ("RELAY_WATER_EVENTS.jsonl", relay_events),
        ("DIRECT_NEAR_MISSES.jsonl", direct_near),
        ("RELAY_NEAR_MISSES.jsonl", relay_near),
    ):
        with (output / name).open("w", encoding="utf-8") as handle:
            for row in rows:
                handle.write(json.dumps(row, sort_keys=True) + "\n")
    status = "PASS_TECHNICAL_STEP2_WATER_NETWORK_SAMPLING" if technical else "NOT_EVALUATED_TECHNICAL_STEP2_WATER_NETWORK_SAMPLING"
    result = {
        "schema_version": 1,
        "status": status,
        "technical_complete": technical,
        "scientific_scope": SCOPE,
        "task_index": manifest["task_index"],
        "seed_index": manifest["seed_index"],
        "replica": manifest["replica"],
        "seed": manifest["seed"],
        "frame_count": len(frames),
        "direct_strict_hit_frame_count": len({int(row["frame"]) for row in direct_hits}),
        "relay_strict_hit_frame_count": len({int(row["frame"]) for row in relay_hits}),
        "direct_event_count": len(direct_events),
        "relay_event_count": len(relay_events),
        "direct_events": direct_events,
        "relay_events": relay_events,
        "direct_near_misses": direct_near,
        "relay_near_misses": relay_near,
        "guard_occupancy": {
            "acyl": sum(item["acyl_guard_pass"] for item in guards) / len(guards) if guards else 0.0,
            "proton": sum(item["proton_guard_pass"] for item in guards) / len(guards) if guards else 0.0,
        },
        "engine": {
            "exit_code": engine_rc,
            "last_nstep": BASE._last_nstep(text),
            "banner": banner,
            "banner_contract_pass": BASE._engine_contract_pass(banner),
            "hard_error_hits": hard,
            "parse_warnings": warnings,
            "parse_error": parse_error,
        },
        "automatic_downstream_action": "NONE",
        "NEXT": NEXT,
    }
    BASE.write_json(output / "RESULT.json", result)
    BASE.write_json(output / ("PASS.json" if technical else "NOT_EVALUATED.json"), result)
    if technical:
        shutil.copy2(restart, output / "final_endpoint.rst7")
    BASE._write_hashes(output)
    return result


def merge_if_ready(output_root: pathlib.Path, array_job: str) -> bool:
    audit_path = output_root / "audit" / f"nylc_a1_step2_water_network_{array_job}.json"
    audit_path.parent.mkdir(parents=True, exist_ok=True)
    lock = audit_path.with_suffix(".json.lock")
    with lock.open("a+") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        paths = [output_root / f"attempt_{array_job}_{index}" / "RESULT.json" for index in range(ARRAY_TASKS)]
        if not all(path.is_file() for path in paths):
            return False
        results = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
        technical_count = sum(item.get("technical_complete") is True for item in results)
        by_seed: dict[str, dict[str, list[dict[str, Any]]]] = {}
        for seed in ("seed26723", "seed26737"):
            seed_results = [item for item in results if item["seed"] == seed]
            by_seed[seed] = {}
            for route in ("direct", "relay"):
                events = [
                    dict(event)
                    for item in seed_results
                    for event in item.get(f"{route}_events", [])
                ]
                by_seed[seed][f"{route}_events"] = sorted(
                    events, key=lambda item: (item["score"], item["task_index"] if "task_index" in item else 0, item["start_frame"])
                )[:MAX_SELECTED_PER_SEED]
                near = [
                    dict(row)
                    for item in seed_results
                    for row in item.get(f"{route}_near_misses", [])
                ]
                by_seed[seed][f"{route}_near_misses"] = sorted(
                    near, key=lambda item: (item["near_miss_score"], item["frame"])
                )[:MAX_SELECTED_PER_SEED]
        direct_reproduced = technical_count == ARRAY_TASKS and all(by_seed[seed]["direct_events"] for seed in by_seed)
        relay_reproduced = technical_count == ARRAY_TASKS and all(by_seed[seed]["relay_events"] for seed in by_seed)
        status = (
            "NOT_EVALUATED_TECHNICAL_STEP2_WATER_NETWORK"
            if technical_count != ARRAY_TASKS
            else "PASS_STEP2_DIRECT_WATER_NETWORK_REPRODUCED"
            if direct_reproduced
            else "PASS_STEP2_TWO_WATER_RELAY_REPRODUCED"
            if relay_reproduced
            else "FAIL_NO_REPRODUCED_STEP2_WATER_NETWORK"
        )
        payload = {
            "schema_version": 1,
            "status": status,
            "technical_denominator": f"{technical_count}/{ARRAY_TASKS}",
            "direct_reproduced": direct_reproduced,
            "relay_reproduced": relay_reproduced,
            "per_seed": by_seed,
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
    parser.add_argument("--mode", required=True, choices=("describe", "initialize", "audit", "merge-if-ready"))
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
        initialize(args.task_index, args.output.resolve(), args.scratch.resolve(), args.code_root.resolve(), args.github_commit)
    elif args.mode == "audit":
        audit(args.output.resolve(), args.scratch.resolve())
    else:
        merge_if_ready(args.output_root.resolve(), args.array_job)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
