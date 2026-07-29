#!/usr/bin/env python3
"""Unrestrained Step1-QM sampling for transient Step2 attack-water organization.

All waters remain MM.  The driver scans every complete water under triclinic
PBC and records only compact, reproducible geometry evidence.  It does not
construct a Step2 product, TS, path, PMF, barrier, or mechanism claim.
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
import shutil
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
EXPECTED_PRMTOP_SHA256 = S2.EXPECTED_PRMTOP_SHA256
SOURCES = S2.SOURCES
REACTIVE = S2.REACTIVE
ARRAY_TASKS = 8
REPLICAS_PER_SEED = 4
MPI_RANKS = 8
MD_STEPS = 4000
MD_DT_PS = 0.0005
MD_NTWX = 20
EXPECTED_FRAMES = MD_STEPS // MD_NTWX
VELOCITY_SEEDS = {
    "seed26723": (26723011, 26723012, 26723013, 26723014),
    "seed26737": (26737011, 26737012, 26737013, 26737014),
}
STEP1_QM_CONTRACT = {
    "qm_atom_count": 146,
    "qmcharge": 0,
    "electron_count": 510,
    "link_atom_count": 6,
    "qm_water_count": 0,
}
HIT_FILTER = {
    "c12_ow_A": (2.70, 3.50),
    "o2_c12_ow_deg": (95.0, 125.0),
    "h_nalpha_A_max": 2.50,
    "ow_h_nalpha_deg_min": 130.0,
}
MIN_CONSECUTIVE_FRAMES = 3
MAX_EVENTS_PER_SEED = 2
EXPECTED_DFTB_DOUBLY_OCCUPIED = 191
SCOPE = "STEP2_WATER_REORGANIZATION_SAMPLING_ONLY_NOT_PRODUCT_TS_PATH_PMF_BARRIER_OR_MECHANISM"
NEXT = "PROMOTE_STRICT_WATER_HITS_TO_149QM_A2_ONLY_AFTER_CROSS_SEED_PASS"
HARD_PATTERNS = {
    "sander_bomb": r"SANDER BOMB",
    "segmentation": r"segmentation",
    "forrtl": r"forrtl",
    "nan": r"\bnan\b",
    "fatal": r"FATAL",
    "scc": r"Convergence could not be achieved",
    "vlimit": r"vlimit\s+exceeded",
    "overflow": r"BOND\s*=\s*\*+",
}


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
    index = int(task_index)
    if not 0 <= index < ARRAY_TASKS:
        raise ValueError("task index must be 0..7")
    seed_index = index // REPLICAS_PER_SEED
    replica = index % REPLICAS_PER_SEED
    source = S2.source_for_index(seed_index)
    return {
        "task_index": index,
        "seed_index": seed_index,
        "replica": replica,
        "seed": source["seed"],
        "velocity_seed": VELOCITY_SEEDS[source["seed"]][replica],
        "source": source,
    }


def sampling_input(seed: str, replica: int, velocity_seed: int, qmmask: str) -> str:
    return f"""NylC A1 Step2 water reorganization {seed} replica{replica}
&cntrl
  imin=0, irest=0, ntx=1, nstlim={MD_STEPS}, dt={MD_DT_PS},
  ntb=1, cut=10.0, ntt=3, gamma_ln=2.0,
  temp0=300.0, tempi=300.0, ig={velocity_seed},
  ntc=1, ntf=1, ntwx={MD_NTWX}, ntpr=100, ntwr={MD_STEPS},
  ioutfm=1, ntxo=1, ifqnt=1, ntr=0, nmropt=0,
/
{S2.qmmm_block(qmmask)}
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
        "step1_qm_contract": STEP1_QM_CONTRACT,
        "reactive_restraints": False,
        "water_position_restraints": False,
        "selected_water_identity_restraint": False,
        "scan_all_complete_waters": True,
        "triclinic_minimum_image": True,
        "hit_filter": HIT_FILTER,
        "minimum_consecutive_frames": MIN_CONSECUTIVE_FRAMES,
        "maximum_events_per_seed": MAX_EVENTS_PER_SEED,
        "velocity_seeds": VELOCITY_SEEDS,
        "source_restart_sha256": {
            item["seed"]: item["restart_sha256"] for item in SOURCES
        },
        "prmtop_sha256": EXPECTED_PRMTOP_SHA256,
        "automatic_downstream_action": "NONE",
    }


def _parse_engine_contract(text: str) -> dict[str, list[int]]:
    start = text.find("QMMM options:")
    region = text[start:] if start >= 0 else ""
    stop = region.find("\n   NSTEP")
    if stop >= 0:
        region = region[:stop]
    patterns = {
        "qm_atom_count": r"\bnquant\s*[=:]\s*(\d+)",
        "qmcharge": r"\bqmcharge\s*[=:]\s*(-?\d+)",
        "spin": r"\bspin\s*[=:]\s*(\d+)",
        "link_atom_count": r"\bnlink\s*[=:]\s*(\d+)",
        "dftb_doubly_occupied_levels": (
            r"RHF\s+CALCULATION,\s+NO\.\s+OF\s+DOUBLY\s+OCCUPIED\s+LEVELS\s*=\s*(\d+)"
        ),
    }
    return {
        key: sorted(set(int(value) for value in re.findall(pattern, region, re.I)))
        for key, pattern in patterns.items()
    }


def _engine_contract_pass(observed: Mapping[str, Sequence[int]]) -> bool:
    expected = {
        "qm_atom_count": [146],
        "qmcharge": [0],
        "spin": [1],
        "link_atom_count": [6],
        "dftb_doubly_occupied_levels": [EXPECTED_DFTB_DOUBLY_OCCUPIED],
    }
    return all(list(observed.get(key, [])) == value for key, value in expected.items())


def _last_nstep(text: str) -> int | None:
    values = re.findall(r"\bNSTEP\s*=\s*(\d+)", text, re.I)
    return int(values[-1]) if values else None


def _guards(geometry: Mapping[str, Any]) -> tuple[bool, bool]:
    gate = AC.PRODUCT_GATE
    acyl = bool(
        geometry["attack_A"] <= gate["attack_A"][1]
        and geometry["c12_n3_A"] >= gate["c12_n3_A_min"]
        and gate["c12_o2_A"][0] <= geometry["c12_o2_A"] <= gate["c12_o2_A"][1]
        and geometry["product_out_of_plane_A"] <= gate["product_out_of_plane_A_max"]
        and geometry["product_angle_sum_deg"] >= gate["product_angle_sum_deg_min"]
    )
    proton = bool(
        gate["hg1_n3_A"][0] <= geometry["hg1_n3_A"] <= gate["hg1_n3_A"][1]
        and geometry["nalpha_hg1_A"] >= gate["nalpha_hg1_A_min"]
        and geometry["qPT_A"] >= gate["qPT_A_min"]
        and geometry["hg1_nearest_qm_heavy_atom"] == REACTIVE["n3"]
    )
    return acyl, proton


def strict_hit(row: Mapping[str, Any]) -> bool:
    return bool(
        HIT_FILTER["c12_ow_A"][0] <= float(row["c12_ow_A"]) <= HIT_FILTER["c12_ow_A"][1]
        and HIT_FILTER["o2_c12_ow_deg"][0]
        <= float(row["o2_c12_ow_deg"])
        <= HIT_FILTER["o2_c12_ow_deg"][1]
        and float(row["h_nalpha_A"]) <= HIT_FILTER["h_nalpha_A_max"]
        and float(row["ow_h_nalpha_deg"]) >= HIT_FILTER["ow_h_nalpha_deg_min"]
        and row.get("acyl_guard_pass") is True
        and row.get("proton_guard_pass") is True
    )


def _hit_score(row: Mapping[str, Any]) -> float:
    return (
        ((float(row["c12_ow_A"]) - 3.05) / 0.35) ** 2
        + ((float(row["o2_c12_ow_deg"]) - 110.0) / 15.0) ** 2
        + ((float(row["h_nalpha_A"]) - 1.90) / 0.60) ** 2
        + ((float(row["ow_h_nalpha_deg"]) - 180.0) / 50.0) ** 2
    )


def collapse_events(rows: Sequence[Mapping[str, Any]]) -> list[dict[str, Any]]:
    grouped: dict[tuple[int, int, int], list[dict[str, Any]]] = {}
    for item in rows:
        row = dict(item)
        if not strict_hit(row):
            continue
        key = (
            int(row.get("replica", 0)),
            int(row["water_oxygen_index1"]),
            int(row["donor_h_index1"]),
        )
        grouped.setdefault(key, []).append(row)
    events: list[dict[str, Any]] = []
    for key, values in grouped.items():
        ordered = sorted(values, key=lambda item: int(item["frame"]))
        run: list[dict[str, Any]] = []
        for row in ordered:
            if run and int(row["frame"]) != int(run[-1]["frame"]) + 1:
                if len(run) >= MIN_CONSECUTIVE_FRAMES:
                    events.append(_event_from_run(key, run))
                run = []
            run.append(row)
        if len(run) >= MIN_CONSECUTIVE_FRAMES:
            events.append(_event_from_run(key, run))
    return sorted(events, key=lambda item: (item["score"], item["replica"], item["start_frame"]))


def _event_from_run(key: tuple[int, int, int], run: Sequence[Mapping[str, Any]]) -> dict[str, Any]:
    representative = min(run, key=lambda row: (_hit_score(row), int(row["frame"])))
    return {
        "replica": key[0],
        "water_oxygen_index1": key[1],
        "donor_h_index1": key[2],
        "start_frame": int(run[0]["frame"]),
        "end_frame": int(run[-1]["frame"]),
        "start_time_ps": float(run[0]["time_ps"]),
        "end_time_ps": float(run[-1]["time_ps"]),
        "consecutive_frames": len(run),
        "representative_frame": int(representative["frame"]),
        "score": _hit_score(representative),
        "representative_geometry": {
            name: representative[name]
            for name in (
                "c12_ow_A",
                "o2_c12_ow_deg",
                "h_nalpha_A",
                "ow_h_nalpha_deg",
                "acyl_guard_pass",
                "proton_guard_pass",
            )
        },
    }


def _complete_water_indices(topology: Any) -> list[tuple[int, tuple[int, int]]]:
    waters: list[tuple[int, tuple[int, int]]] = []
    for residue in topology.residues:
        parsed = S2._complete_h2o(residue)
        if parsed is None:
            continue
        oxygen, hydrogens = parsed
        waters.append((oxygen.idx + 1, tuple(atom.idx + 1 for atom in hydrogens)))
    return waters


def _scan_frames(
    trajectory: pathlib.Path, manifest: Mapping[str, Any]
) -> tuple[list[Any], list[dict[str, Any]], list[dict[str, Any]], list[str]]:
    import parmed as pmd

    topology = pmd.load_file(str(PRMTOP))
    frames, warnings = S2._read_md_frames(trajectory, len(topology.atoms))
    waters = _complete_water_indices(topology)
    cell = S2._periodic_cell(manifest["box"])
    heavy = [int(value) for value in manifest["qm_heavy_atom_indices"]]
    hits: list[dict[str, Any]] = []
    frame_guards: list[dict[str, Any]] = []
    for frame_index, frame in enumerate(frames):
        coords = {index + 1: frame[index] for index in range(len(frame))}
        acyl_geometry = AC.geometry_from_coordinates(coords, heavy)
        acyl_guard, proton_guard = _guards(acyl_geometry)
        frame_guards.append({
            "frame": frame_index,
            "acyl_guard_pass": acyl_guard,
            "proton_guard_pass": proton_guard,
        })
        if not (acyl_guard and proton_guard):
            continue
        c12 = coords[REACTIVE["c12"]]
        o2 = coords[REACTIVE["o2"]]
        nalpha = coords[REACTIVE["nalpha"]]
        for oxygen_index1, hydrogen_indices1 in waters:
            oxygen = coords[oxygen_index1]
            c12_ow = S2._distance(c12, oxygen, cell)
            if not HIT_FILTER["c12_ow_A"][0] <= c12_ow <= HIT_FILTER["c12_ow_A"][1]:
                continue
            attack_angle = S2._angle(o2, c12, oxygen, cell)
            if not HIT_FILTER["o2_c12_ow_deg"][0] <= attack_angle <= HIT_FILTER["o2_c12_ow_deg"][1]:
                continue
            donors = []
            for hydrogen_index1 in hydrogen_indices1:
                hydrogen = coords[hydrogen_index1]
                distance = S2._distance(hydrogen, nalpha, cell)
                angle = S2._angle(oxygen, hydrogen, nalpha, cell)
                donors.append((distance, -angle, hydrogen_index1, angle))
            h_distance, _, donor_h, h_angle = min(donors)
            row = {
                "frame": frame_index,
                "time_ps": (frame_index + 1) * MD_NTWX * MD_DT_PS,
                "replica": manifest["replica"],
                "water_oxygen_index1": oxygen_index1,
                "water_hydrogen_indices1": list(hydrogen_indices1),
                "donor_h_index1": donor_h,
                "c12_ow_A": c12_ow,
                "o2_c12_ow_deg": attack_angle,
                "h_nalpha_A": h_distance,
                "ow_h_nalpha_deg": h_angle,
                "acyl_guard_pass": acyl_guard,
                "proton_guard_pass": proton_guard,
            }
            if strict_hit(row):
                hits.append(row)
    return frames, hits, frame_guards, warnings


def _write_selected_restart(
    destination: pathlib.Path, frame: Any, box: Sequence[float]
) -> None:
    import parmed as pmd

    structure = pmd.load_file(str(PRMTOP))
    structure.coordinates = frame
    structure.box = box
    structure.save(str(destination), overwrite=False)


def initialize(
    task_index: int,
    output: pathlib.Path,
    scratch: pathlib.Path,
    code_root: pathlib.Path,
    github_commit: str,
) -> dict[str, Any]:
    if output.exists() or scratch.exists():
        raise FileExistsError("output or scratch already exists")
    if not re.fullmatch(r"[0-9a-f]{40}", github_commit):
        raise ValueError("runtime is not bound to a full Git commit")
    task = task_spec(task_index)
    source = task["source"]
    source_manifest, _ = S2.validate_authority(source, code_root)
    contract = source_manifest["qm_contract"]
    checks = {
        "qm_atom_count": contract.get("qm_atom_count"),
        "qmcharge": contract.get("qmcharge"),
        "electron_count": contract.get("electron_count_including_link_h"),
        "link_atom_count": contract.get("link_atom_count"),
        "qm_water_count": contract.get("step1_qm_water_count"),
    }
    if checks != STEP1_QM_CONTRACT:
        raise ValueError(f"Step1 QM contract changed: {checks}")
    if sha256(PRMTOP) != EXPECTED_PRMTOP_SHA256:
        raise ValueError("prmtop SHA256 mismatch")
    import parmed as pmd

    structure = pmd.load_file(str(PRMTOP), xyz=str(source["restart"]))
    AC._validate_structure(structure)
    box = [float(value) for value in structure.box[:6]]
    S2._periodic_cell(box)
    qmmask = contract["qmmask"]
    base_indices = S2._parse_qmmask(qmmask, 146)
    heavy = [
        index for index in base_indices
        if int(getattr(structure.atoms[index - 1], "atomic_number", 0) or 0) > 1
    ]
    output.mkdir(parents=True, exist_ok=False)
    scratch.mkdir(parents=True, exist_ok=False)
    (scratch / "stage.in").write_text(
        sampling_input(task["seed"], task["replica"], task["velocity_seed"], qmmask),
        encoding="utf-8",
    )
    manifest = {
        "schema_version": 1,
        "status": "READY_STEP2_WATER_REORGANIZATION_SAMPLING",
        "scientific_scope": SCOPE,
        "github_commit": github_commit,
        "task_index": task["task_index"],
        "seed_index": task["seed_index"],
        "replica": task["replica"],
        "seed": task["seed"],
        "velocity_seed": task["velocity_seed"],
        "source": {
            "attempt": source["attempt"],
            "restart": str(source["restart"]),
            "restart_sha256": source["restart_sha256"],
        },
        "prmtop": str(PRMTOP),
        "prmtop_sha256": EXPECTED_PRMTOP_SHA256,
        "box": box,
        "qm_contract": checks,
        "qmmask": qmmask,
        "qm_heavy_atom_indices": heavy,
        "water_universe": {
            "complete_water_count": len(_complete_water_indices(structure)),
            "identity_fixed": False,
            "all_waters_remain_mm": True,
        },
        "protocol": describe(),
        "automatic_downstream_action": "NONE",
        "NEXT": NEXT,
    }
    write_json(output / "SOURCE_MANIFEST.json", manifest)
    return manifest


def audit(output: pathlib.Path, scratch: pathlib.Path) -> dict[str, Any]:
    manifest = json.loads((output / "SOURCE_MANIFEST.json").read_text(encoding="utf-8"))
    stage_out = scratch / "stage.out"
    trajectory = scratch / "stage.nc"
    restart = scratch / "stage.rst7"
    engine_rc = int((scratch / "engine.rc").read_text().strip()) if (scratch / "engine.rc").is_file() else 999
    text = stage_out.read_text(encoding="utf-8", errors="replace") if stage_out.is_file() else ""
    banner = _parse_engine_contract(text)
    hard = {name: len(re.findall(pattern, text, re.I)) for name, pattern in HARD_PATTERNS.items()}
    frames: list[Any] = []
    hits: list[dict[str, Any]] = []
    guards: list[dict[str, Any]] = []
    parse_warnings: list[str] = []
    parse_error = None
    if trajectory.is_file() and trajectory.stat().st_size:
        try:
            frames, hits, guards, parse_warnings = _scan_frames(trajectory, manifest)
        except Exception as error:
            parse_error = f"{type(error).__name__}: {error}"
    technical = bool(
        engine_rc == 0
        and _last_nstep(text) == MD_STEPS
        and restart.is_file()
        and restart.stat().st_size > 0
        and trajectory.is_file()
        and trajectory.stat().st_size > 0
        and len(frames) == EXPECTED_FRAMES
        and _engine_contract_pass(banner)
        and sum(hard.values()) == 0
        and not parse_warnings
        and parse_error is None
    )
    events = collapse_events(hits)[:MAX_EVENTS_PER_SEED] if technical else []
    trajectory_sha = sha256(trajectory) if trajectory.is_file() else ""
    for rank, event in enumerate(events):
        event["trajectory_sha256"] = trajectory_sha
        event["source_restart_sha256"] = manifest["source"]["restart_sha256"]
        event["seed"] = manifest["seed"]
        event["task_index"] = manifest["task_index"]
        selected = output / f"selected_water_event_{rank}.rst7"
        _write_selected_restart(selected, frames[event["representative_frame"]], manifest["box"])
        event["selected_restart"] = selected.name
        event["selected_restart_sha256"] = sha256(selected)
    with (output / "WATER_HITS.jsonl").open("w", encoding="utf-8") as handle:
        for event in events:
            handle.write(json.dumps(event, sort_keys=True) + "\n")
    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_STEP2_WATER_REORGANIZATION_SAMPLING"
            if technical else "NOT_EVALUATED_TECHNICAL_STEP2_WATER_REORGANIZATION_SAMPLING"
        ),
        "scientific_scope": SCOPE,
        "technical_complete": technical,
        "task_index": manifest["task_index"],
        "seed_index": manifest["seed_index"],
        "replica": manifest["replica"],
        "seed": manifest["seed"],
        "frame_count": len(frames),
        "strict_hit_frame_count": len(hits),
        "event_count": len(events),
        "events": events,
        "guard_occupancy": {
            "acyl": sum(item["acyl_guard_pass"] for item in guards) / len(guards) if guards else 0.0,
            "proton": sum(item["proton_guard_pass"] for item in guards) / len(guards) if guards else 0.0,
        },
        "engine": {
            "exit_code": engine_rc,
            "last_nstep": _last_nstep(text),
            "banner": banner,
            "banner_contract_pass": _engine_contract_pass(banner),
            "hard_error_hits": hard,
            "trajectory_sha256": trajectory_sha,
            "parse_warnings": parse_warnings,
            "parse_error": parse_error,
        },
        "automatic_downstream_action": "NONE",
        "NEXT": NEXT,
    }
    write_json(output / "RESULT.json", result)
    marker = output / ("PASS.json" if technical else "NOT_EVALUATED.json")
    write_json(marker, result)
    if technical:
        shutil.copy2(restart, output / "final_endpoint.rst7")
    _write_hashes(output)
    return result


def _write_hashes(output: pathlib.Path) -> None:
    names = [
        path.name for path in output.iterdir()
        if path.is_file() and path.name != "SHA256.tsv"
    ]
    rows = [f"{sha256(output / name)}  {name}" for name in sorted(names)]
    (output / "SHA256.tsv").write_text("\n".join(rows) + ("\n" if rows else ""), encoding="utf-8")


def merge_if_ready(output_root: pathlib.Path, array_job: str) -> bool:
    audit = output_root / "audit" / f"nylc_a1_step2_water_reorganization_{array_job}.json"
    audit.parent.mkdir(parents=True, exist_ok=True)
    lock = audit.with_suffix(".json.lock")
    with lock.open("a+") as handle:
        fcntl.flock(handle.fileno(), fcntl.LOCK_EX)
        paths = [output_root / f"attempt_{array_job}_{index}" / "RESULT.json" for index in range(ARRAY_TASKS)]
        if not all(path.is_file() for path in paths):
            return False
        results = [json.loads(path.read_text(encoding="utf-8")) for path in paths]
        technical = all(item.get("technical_complete") is True for item in results)
        selected_by_seed: dict[str, list[dict[str, Any]]] = {}
        for seed in ("seed26723", "seed26737"):
            candidates = [
                dict(event)
                for result in results if result["seed"] == seed
                for event in result.get("events", [])
            ]
            selected_by_seed[seed] = sorted(
                candidates,
                key=lambda item: (item["score"], item["task_index"], item["start_frame"]),
            )[:MAX_EVENTS_PER_SEED]
        reproduced = technical and all(selected_by_seed[seed] for seed in selected_by_seed)
        status = (
            "NOT_EVALUATED_TECHNICAL_STEP2_WATER_REORGANIZATION"
            if not technical
            else "PASS_STEP2_PREORGANIZED_WATER_REPRODUCED"
            if reproduced
            else "FAIL_NO_REPRODUCED_STEP2_PREORGANIZED_WATER"
        )
        payload = {
            "schema_version": 1,
            "status": status,
            "technical_denominator": f"{sum(item.get('technical_complete') is True for item in results)}/8",
            "per_seed_selected_events": selected_by_seed,
            "per_task": sorted(results, key=lambda item: item["task_index"]),
            "automatic_downstream_action": "NONE",
            "NEXT": NEXT,
        }
        if audit.exists():
            if json.loads(audit.read_text(encoding="utf-8")) != payload:
                raise FileExistsError(f"existing audit differs: {audit}")
        else:
            write_json(audit, payload)
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
