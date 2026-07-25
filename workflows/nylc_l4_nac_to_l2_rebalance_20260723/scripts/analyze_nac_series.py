#!/usr/bin/env python3
"""Audit paired GROMACS distance/angle series for strict NAC residence."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path


def read_xvg(path: Path):
    rows = []
    for line_number, raw in enumerate(path.read_text().splitlines(), 1):
        line = raw.strip()
        if not line or line.startswith(("#", "@")):
            continue
        fields = line.split()
        if len(fields) < 2:
            raise ValueError(f"{path}:{line_number}: expected at least two columns")
        time, value = float(fields[0]), float(fields[1])
        if not (math.isfinite(time) and math.isfinite(value)):
            raise ValueError(f"{path}:{line_number}: non-finite value")
        rows.append((time, value))
    if not rows:
        raise ValueError(f"{path}: no data rows")
    return rows


def audit_series(
    distance_rows,
    angle_rows,
    distance_max,
    angle_min,
    angle_max,
    potential_rows=None,
    analysis_start_ps=None,
    analysis_end_ps=None,
    sample_interval_ps=2.0,
):
    if len(distance_rows) != len(angle_rows):
        raise ValueError("time series do not align: different frame counts")
    if potential_rows is not None and len(distance_rows) != len(potential_rows):
        raise ValueError("time series do not align: different frame counts")

    all_rows = []
    for index, ((distance_time, distance), (angle_time, angle)) in enumerate(
        zip(distance_rows, angle_rows)
    ):
        if abs(distance_time - angle_time) > 1e-6:
            raise ValueError("time series do not align: frame times differ")
        row = {
            "time_ps": distance_time,
            "distance_nm": distance,
            "angle_deg": angle,
            "nac": distance <= distance_max and angle_min <= angle <= angle_max,
        }
        if potential_rows is not None:
            potential_time, potential = potential_rows[index]
            if abs(distance_time - potential_time) > 1e-6:
                raise ValueError("time series do not align: frame times differ")
            row["potential_energy_kj_mol"] = potential
        all_rows.append(row)

    if not math.isfinite(sample_interval_ps) or sample_interval_ps <= 0:
        raise ValueError("sample interval must be a positive finite value")
    window_start = all_rows[0]["time_ps"] if analysis_start_ps is None else analysis_start_ps
    window_end = all_rows[-1]["time_ps"] if analysis_end_ps is None else analysis_end_ps
    if not (math.isfinite(window_start) and math.isfinite(window_end)):
        raise ValueError("analysis window must be finite")
    if window_start > window_end:
        raise ValueError("analysis window start exceeds end")

    rows = [
        row
        for row in all_rows
        if window_start - 1e-6 <= row["time_ps"] <= window_end + 1e-6
    ]
    if not rows:
        raise ValueError("analysis window contains no frames")

    events = []
    current = []

    def finish_event():
        if not current:
            return
        event = {
            "event_id": f"event_{len(events) + 1:04d}",
            "start_ps": current[0]["time_ps"],
            "end_ps": current[-1]["time_ps"],
            "duration_ps": current[-1]["time_ps"] - current[0]["time_ps"],
            "frame_count": len(current),
            "member_times_ps": [item["time_ps"] for item in current],
            "mean_distance_nm": sum(item["distance_nm"] for item in current) / len(current),
            "mean_angle_deg": sum(item["angle_deg"] for item in current) / len(current),
        }
        if potential_rows is not None:
            event["minimum_potential_energy_kj_mol"] = min(
                item["potential_energy_kj_mol"] for item in current
            )
        events.append(event)
        current.clear()

    for row in rows:
        if not row["nac"]:
            finish_event()
            continue
        if (
            current
            and row["time_ps"] - current[-1]["time_ps"] > sample_interval_ps + 1e-6
        ):
            finish_event()
        current.append(row)
    finish_event()

    longest = None
    if events:
        longest = max(
            events,
            key=lambda event: (
                event["duration_ps"],
                event["frame_count"],
                -event["mean_distance_nm"],
                -abs(event["mean_angle_deg"] - 105.0),
            ),
        )

    nac_rows = [row for row in rows if row["nac"]]
    lowest_potential = None
    if potential_rows is not None and nac_rows:
        selected = min(nac_rows, key=lambda row: row["potential_energy_kj_mol"])
        lowest_potential = {
            "time_ps": selected["time_ps"],
            "potential_energy_kj_mol": selected["potential_energy_kj_mol"],
            "distance_nm": selected["distance_nm"],
            "angle_deg": selected["angle_deg"],
        }

    return {
        "schema_version": 2,
        "analysis_window_ps": [window_start, window_end],
        "sample_interval_ps": sample_interval_ps,
        "frame_count": len(rows),
        "nac_frame_count": len(nac_rows),
        "nac_occupancy": len(nac_rows) / len(rows),
        "gates": {
            "distance_max_nm": distance_max,
            "angle_min_deg": angle_min,
            "angle_max_deg": angle_max,
        },
        "nac_events": events,
        "longest_continuous_nac": longest,
        "lowest_potential_nac_frame": lowest_potential,
        "nac_run_count": len(events),
        "time_start_ps": rows[0]["time_ps"],
        "time_end_ps": rows[-1]["time_ps"],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--distance-xvg", type=Path, required=True)
    parser.add_argument("--angle-xvg", type=Path, required=True)
    parser.add_argument("--potential-xvg", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--distance-max-nm", type=float, default=0.35)
    parser.add_argument("--angle-min-deg", type=float, default=95.0)
    parser.add_argument("--angle-max-deg", type=float, default=115.0)
    parser.add_argument("--analysis-start-ps", type=float)
    parser.add_argument("--analysis-end-ps", type=float)
    parser.add_argument("--sample-interval-ps", type=float, default=2.0)
    args = parser.parse_args()
    try:
        result = audit_series(
            read_xvg(args.distance_xvg),
            read_xvg(args.angle_xvg),
            args.distance_max_nm,
            args.angle_min_deg,
            args.angle_max_deg,
            read_xvg(args.potential_xvg) if args.potential_xvg else None,
            analysis_start_ps=args.analysis_start_ps,
            analysis_end_ps=args.analysis_end_ps,
            sample_interval_ps=args.sample_interval_ps,
        )
    except ValueError as error:
        print(str(error), file=sys.stderr)
        return 2
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
