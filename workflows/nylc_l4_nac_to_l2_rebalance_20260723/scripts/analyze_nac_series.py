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
):
    if len(distance_rows) != len(angle_rows):
        raise ValueError("time series do not align: different frame counts")
    if potential_rows is not None and len(distance_rows) != len(potential_rows):
        raise ValueError("time series do not align: different frame counts")
    rows = []
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
        rows.append(row)

    runs = []
    start = None
    for index, row in enumerate(rows + [{"nac": False}]):
        if row["nac"] and start is None:
            start = index
        elif not row["nac"] and start is not None:
            stop = index - 1
            subset = rows[start : stop + 1]
            runs.append(
                {
                    "start_ps": subset[0]["time_ps"],
                    "end_ps": subset[-1]["time_ps"],
                    "duration_ps": subset[-1]["time_ps"] - subset[0]["time_ps"],
                    "frame_count": len(subset),
                    "mean_distance_nm": sum(item["distance_nm"] for item in subset) / len(subset),
                    "mean_angle_deg": sum(item["angle_deg"] for item in subset) / len(subset),
                }
            )
            start = None

    longest = None
    if runs:
        longest = max(
            runs,
            key=lambda run: (
                run["duration_ps"],
                run["frame_count"],
                -run["mean_distance_nm"],
                -abs(run["mean_angle_deg"] - 105.0),
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
        "schema_version": 1,
        "frame_count": len(rows),
        "nac_frame_count": len(nac_rows),
        "nac_occupancy": len(nac_rows) / len(rows),
        "gates": {
            "distance_max_nm": distance_max,
            "angle_min_deg": angle_min,
            "angle_max_deg": angle_max,
        },
        "longest_continuous_nac": longest,
        "lowest_potential_nac_frame": lowest_potential,
        "nac_run_count": len(runs),
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
    args = parser.parse_args()
    try:
        result = audit_series(
            read_xvg(args.distance_xvg),
            read_xvg(args.angle_xvg),
            args.distance_max_nm,
            args.angle_min_deg,
            args.angle_max_deg,
            read_xvg(args.potential_xvg) if args.potential_xvg else None,
        )
    except ValueError as error:
        print(str(error), file=sys.stderr)
        return 2
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
