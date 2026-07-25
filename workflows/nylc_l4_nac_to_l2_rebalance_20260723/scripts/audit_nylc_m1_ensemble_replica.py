#!/usr/bin/env python3
"""Audit one fully unrestrained NylC M1 ensemble replica from primitive metrics."""

from __future__ import annotations

import argparse
import json
import math
import pathlib
import re
import statistics
from typing import Iterable

from analyze_nac_series import audit_series, read_xvg

HARD_PATTERNS = {
    "fatal": r"(?i)fatal(?:\s+error)?",
    "lincs_warning": r"(?i)lincs\s+warning",
    "settle_problem": r"(?i)(?:cannot\s+be\s+settled|settle[^\n]*(?:warning|error|failed))",
    "nan": r"(?i)\bnan\b",
}
CHANNELS = (
    "OgH_to_Nalpha_preorg",
    "OgH_to_Asp306",
    "OgH_to_Asp308",
    "OgH_to_water",
    "OgH_via_water_to_Asp306",
    "OgH_via_water_to_Asp308",
)


def _read_table(path: pathlib.Path, minimum_columns: int) -> list[list[float]]:
    rows = []
    for line_number, raw in enumerate(path.read_text(errors="strict").splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith(("#", "@")):
            continue
        fields = [float(value) for value in line.split()]
        if len(fields) < minimum_columns:
            raise ValueError(f"{path}:{line_number}: expected {minimum_columns} columns")
        if not all(math.isfinite(value) for value in fields):
            raise ValueError(f"{path}:{line_number}: non-finite value")
        rows.append(fields)
    if not rows:
        raise ValueError(f"{path}: no data rows")
    return rows


def _window(rows: Iterable[list[float]], start: float, end: float) -> list[list[float]]:
    selected = [row for row in rows if start - 1.0e-6 <= row[0] <= end + 1.0e-6]
    if not selected:
        raise ValueError("analysis window contains no frames")
    return selected


def _validate_times(reference: list[list[float]], other: list[list[float]], label: str) -> None:
    if len(reference) != len(other):
        raise ValueError(f"{label}: frame count differs from NAC series")
    for left, right in zip(reference, other):
        if abs(left[0] - right[0]) > 1.0e-6:
            raise ValueError(f"{label}: frame times differ from NAC series")


def _stats(values: Iterable[float]) -> dict:
    values = [float(value) for value in values]
    if not values or not all(math.isfinite(value) for value in values):
        raise ValueError("statistics require finite values")
    return {
        "frame_count": len(values),
        "mean": statistics.fmean(values),
        "stdev": statistics.pstdev(values),
        "min": min(values),
        "max": max(values),
    }


def _median_spacing(rows: list[list[float]]) -> float:
    if len(rows) < 2:
        raise ValueError("at least two frames are required")
    gaps = [right[0] - left[0] for left, right in zip(rows, rows[1:])]
    if not all(gap > 0 and math.isfinite(gap) for gap in gaps):
        raise ValueError("frame times must be strictly increasing")
    return statistics.median(gaps)


def _scan_log(path: pathlib.Path) -> dict:
    text = path.read_text(errors="replace")
    counts = {name: len(re.findall(pattern, text)) for name, pattern in HARD_PATTERNS.items()}
    counts["finished_mdrun"] = bool(re.search(r"Finished mdrun", text))
    return counts


def _thermodynamics(rows: list[list[float]]) -> dict:
    temperature = _stats(row[1] for row in rows)
    pressure = _stats(row[2] for row in rows)
    volume = _stats(row[3] for row in rows)
    mean_volume = abs(volume["mean"])
    volume_cv = volume["stdev"] / mean_volume if mean_volume else math.inf
    midpoint = len(rows) // 2
    first = rows[:midpoint] or rows[:1]
    second = rows[midpoint:] or rows[-1:]
    volume_drift = abs(
        statistics.fmean(row[3] for row in first)
        - statistics.fmean(row[3] for row in second)
    ) / mean_volume if mean_volume else math.inf
    stable = (
        295.0 <= temperature["mean"] <= 305.0
        and temperature["min"] >= 270.0
        and temperature["max"] <= 330.0
        and abs(pressure["mean"]) <= 100.0
        and volume_cv <= 0.02
        and volume_drift <= 0.02
    )
    return {
        "temperature_K": temperature,
        "pressure_bar": pressure,
        "volume_nm3": volume,
        "volume_coefficient_of_variation": volume_cv,
        "volume_half_drift_fraction": volume_drift,
        "stable": stable,
        "pressure_variance_is_not_a_failure": True,
    }


def _proton_summary(path: pathlib.Path, start: float, end: float, nac_by_time: dict) -> dict:
    records = []
    for line_number, raw in enumerate(path.read_text(errors="strict").splitlines(), 1):
        if not raw.strip():
            continue
        record = json.loads(raw)
        time = float(record["time_ps"])
        if start - 1.0e-6 <= time <= end + 1.0e-6:
            records.append(record)
    if not records:
        raise ValueError(f"{path}: no proton-network rows in analysis window")
    result = {}
    for channel in CHANNELS:
        favorable = []
        nac_favorable = []
        for record in records:
            value = bool(record.get("channels", {}).get(channel, {}).get("favorable", False))
            favorable.append(value)
            nac_favorable.append(value and nac_by_time.get(round(float(record["time_ps"]), 6), False))
        nac_denominator = sum(
            bool(nac_by_time.get(round(float(record["time_ps"]), 6), False))
            for record in records
        )
        result[channel] = {
            "denominator_frames": len(records),
            "favorable_frames": sum(favorable),
            "favorable_occupancy": sum(favorable) / len(records),
            "nac_denominator_frames": nac_denominator,
            "nac_favorable_frames": sum(nac_favorable),
            "nac_favorable_occupancy": (
                sum(nac_favorable) / nac_denominator if nac_denominator else None
            ),
            "interpretation": "fixed_topology_geometry_preorganization_only",
        }
    return result


def audit_replica(
    run_root: pathlib.Path,
    manifest: dict,
    window: tuple[float, float],
) -> dict:
    root = pathlib.Path(run_root)
    start, end = map(float, window)
    if start > end or not all(math.isfinite(value) for value in (start, end)):
        raise ValueError("invalid analysis window")
    if not manifest.get("fully_unrestrained"):
        raise ValueError("replica manifest is not fully unrestrained")

    distance_all = read_xvg(root / "nac_distance.xvg")
    angle_all = read_xvg(root / "nac_angle.xvg")
    sample_interval = statistics.median(
        right[0] - left[0] for left, right in zip(distance_all, distance_all[1:])
    )
    nac = audit_series(
        distance_all,
        angle_all,
        0.35,
        95.0,
        115.0,
        analysis_start_ps=start,
        analysis_end_ps=end,
        sample_interval_ps=sample_interval,
    )
    after_start = start + 20.0
    after = audit_series(
        distance_all,
        angle_all,
        0.35,
        95.0,
        115.0,
        analysis_start_ps=after_start,
        analysis_end_ps=end,
        sample_interval_ps=sample_interval,
    )
    nac["analysis_frame_count_after_20ps"] = after["frame_count"]
    nac["nac_frame_count_after_20ps"] = after["nac_frame_count"]
    nac["nac_occupancy_after_20ps"] = after["nac_occupancy"]
    nac["longest_event"] = nac.pop("longest_continuous_nac") or {
        "duration_ps": 0.0,
        "frame_count": 0,
    }
    nac["longest_event_after_20ps"] = after["longest_continuous_nac"] or {
        "duration_ps": 0.0,
        "frame_count": 0,
    }

    distance = _window(_read_table(root / "nac_distance.xvg", 2), start, end)
    angle = _window(_read_table(root / "nac_angle.xvg", 2), start, end)
    gate = _window(_read_table(root / "gate_opening.xvg", 2), start, end)
    thermo = _window(_read_table(root / "thermo.xvg", 4), start, end)
    pocket = _window(_read_table(root / "pocket_state.xvg", 3), start, end)
    for label, rows in (("angle", angle), ("gate", gate), ("thermo", thermo), ("pocket", pocket)):
        _validate_times(distance, rows, label)

    retained = [row[1] >= 3.0 and row[2] <= 1.2 for row in pocket]
    contact = json.loads((root / "minimum_contact.json").read_text())
    minimum_protein = float(contact["minimum_ligand_protein_heavy_nm"])
    minimum_water = float(contact["minimum_ligand_water_heavy_nm"])
    if not all(math.isfinite(value) for value in (minimum_protein, minimum_water)):
        raise ValueError("minimum contact contains non-finite value")
    severe_clash = minimum_protein < 0.18

    nac_by_time = {
        round(row[0], 6): bool(row[1] <= 0.35 and 95.0 <= other[1] <= 115.0)
        for row, other in zip(distance, angle)
    }
    numerical = _scan_log(root / "run.log")
    numerical_pass = numerical.pop("finished_mdrun") and not any(numerical.values())
    thermo_audit = _thermodynamics(thermo)

    if not numerical_pass:
        scientific = "NOT_EVALUATED_NUMERICAL_FAIL"
    elif severe_clash:
        scientific = "FAIL_SEVERE_CLASH"
    elif not retained[-1]:
        scientific = "FAIL_UNBOUND"
    elif not thermo_audit["stable"]:
        scientific = "NOT_EVALUATED_THERMODYNAMIC_FAIL"
    elif nac["nac_frame_count"] == 0:
        scientific = "FAIL_REPLICA_NO_NAC"
    else:
        scientific = "PASS_REPLICA_NAC_PRESENT"

    potential_path = root / "potential_energy.xvg"
    potential = None
    if potential_path.is_file():
        potential_rows = _window(_read_table(potential_path, 2), start, end)
        _validate_times(distance, potential_rows, "potential")
        potential = _stats(row[1] for row in potential_rows)

    return {
        "schema_version": 1,
        "candidate_id": manifest["candidate_id"],
        "velocity_seed": int(manifest["velocity_seed"]),
        "microstate": manifest["microstate"],
        "gate_definition": manifest["gate_definition"],
        "fully_unrestrained": True,
        "analysis": {
            "window_ps": [start, end],
            "frame_count": len(distance),
            "sampling_interval_ps": _median_spacing(distance),
            "actual_denominator": len(distance),
        },
        "technical_status": "PASS" if numerical_pass else "FAIL",
        "scientific_status": scientific,
        "nac": nac,
        "gate_opening_nm": _stats(row[1] for row in gate),
        "thermodynamics": thermo_audit,
        "potential_energy_kj_mol": potential,
        "bound_state": {
            "definition": "contacts_le_0.45nm>=3 and ligand_pocket_com<=1.2nm",
            "retained_frame_count": sum(retained),
            "retained_occupancy": sum(retained) / len(retained),
            "final_frame_retained": retained[-1],
            "minimum_ligand_protein_heavy_nm": minimum_protein,
            "minimum_ligand_water_heavy_nm": minimum_water,
            "severe_clash": severe_clash,
            "pocket_contact_count": _stats(row[1] for row in pocket),
            "ligand_pocket_com_nm": _stats(row[2] for row in pocket),
        },
        "proton_preorganization": _proton_summary(
            root / "proton_network.jsonl", start, end, nac_by_time
        ),
        "numerical_issue_counts": numerical,
        "scientific_scope": "fixed_topology_preorganization_not_proton_transfer",
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", type=pathlib.Path, required=True)
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--window-start-ps", type=float, required=True)
    parser.add_argument("--window-end-ps", type=float, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    audit = audit_replica(
        args.run_root,
        json.loads(args.manifest.read_text()),
        (args.window_start_ps, args.window_end_ps),
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    print(json.dumps(audit, sort_keys=True))
    return 0 if audit["technical_status"] == "PASS" else 2


if __name__ == "__main__":
    raise SystemExit(main())
