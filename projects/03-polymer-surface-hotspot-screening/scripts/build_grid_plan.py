#!/usr/bin/env python3
"""Build deterministic whole-receptor AutoGrid tile plans with shared A/C/OA/HD geometry."""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
import sys
from pathlib import Path


CHANNELS = ["A", "C", "OA", "HD"]


def canonical_geometry(geometry: dict) -> str:
    return json.dumps(geometry, sort_keys=True, separators=(",", ":"))


def geometry_sha256(geometry: dict) -> str:
    return hashlib.sha256(canonical_geometry(geometry).encode("utf-8")).hexdigest()


def validate_channel_geometries(geometries: dict[str, dict]) -> None:
    missing = [channel for channel in CHANNELS if channel not in geometries]
    if missing:
        raise ValueError(f"missing channel geometry: {','.join(missing)}")
    fingerprints = {
        canonical_geometry({
            "center": geometries[channel]["center"],
            "npts": geometries[channel]["npts"],
            "spacing": geometries[channel]["spacing"],
        })
        for channel in CHANNELS
    }
    if len(fingerprints) != 1:
        raise ValueError("nonidentical channel geometry")


def even_ceiling(value: float) -> int:
    result = int(math.ceil(value - 1e-12))
    return result if result % 2 == 0 else result + 1


def axis_windows(
    lower: float,
    upper: float,
    spacing: float,
    overlap: float,
    max_npts: int,
) -> list[dict]:
    extent = upper - lower
    if extent <= 0:
        raise ValueError("bounding box extent must be positive")
    if max_npts < 2 or max_npts % 2:
        raise ValueError("max_npts must be a positive even integer")
    if spacing <= 0 or overlap < 0:
        raise ValueError("spacing must be positive and overlap nonnegative")

    required_npts = even_ceiling(extent / spacing)
    if required_npts <= max_npts:
        npts = max(2, required_npts)
        span = npts * spacing
        center = (lower + upper) / 2.0
        return [{
            "center": center,
            "npts": npts,
            "lower": center - span / 2.0,
            "upper": center + span / 2.0,
        }]

    span = max_npts * spacing
    stride = span - overlap
    if stride <= 0:
        raise ValueError("overlap must be smaller than maximum tile span")

    starts = [lower]
    final_start = upper - span
    while starts[-1] < final_start - 1e-9:
        candidate = starts[-1] + stride
        if candidate >= final_start:
            candidate = final_start
        if candidate <= starts[-1] + 1e-9:
            break
        starts.append(candidate)

    return [{
        "center": start + span / 2.0,
        "npts": max_npts,
        "lower": start,
        "upper": start + span,
    } for start in starts]


def rounded(values):
    return [round(float(value), 6) for value in values]


def build_plan(
    bbox_min: list[float],
    bbox_max: list[float],
    spacing: float,
    margin: float,
    overlap: float,
    max_npts: int,
) -> dict:
    if len(bbox_min) != 3 or len(bbox_max) != 3:
        raise ValueError("bounding boxes must contain three coordinates")
    if margin < 0:
        raise ValueError("margin must be nonnegative")

    target_lower = [float(value) - margin for value in bbox_min]
    target_upper = [float(value) + margin for value in bbox_max]
    windows = [
        axis_windows(target_lower[i], target_upper[i], spacing, overlap, max_npts)
        for i in range(3)
    ]

    tiles = []
    for tile_index, axes in enumerate(itertools.product(*windows), start=1):
        center = rounded(axis["center"] for axis in axes)
        npts = [int(axis["npts"]) for axis in axes]
        lower = rounded(axis["lower"] for axis in axes)
        upper = rounded(axis["upper"] for axis in axes)
        geometry = {"center": center, "npts": npts, "spacing": float(spacing)}
        fingerprint = geometry_sha256(geometry)
        channel_geometry = {channel: dict(geometry) for channel in CHANNELS}
        validate_channel_geometries(channel_geometry)
        tiles.append({
            "tile_id": f"tile_{tile_index:04d}",
            "center": center,
            "npts": npts,
            "spacing": float(spacing),
            "lower": lower,
            "upper": upper,
            "geometry_sha256": fingerprint,
            "channel_geometry_sha256": {
                channel: fingerprint for channel in CHANNELS
            },
        })

    return {
        "status": "GRID_PLAN_PASS",
        "schema_version": "whole_receptor_grid_plan_v1",
        "channels": CHANNELS,
        "spacing_angstrom": float(spacing),
        "margin_angstrom": float(margin),
        "overlap_angstrom": float(overlap),
        "max_npts": int(max_npts),
        "target_lower": rounded(target_lower),
        "target_upper": rounded(target_upper),
        "tile_count": len(tiles),
        "tiles": tiles,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptor-metadata", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--spacing", type=float, default=0.75)
    parser.add_argument("--margin", type=float, default=5.0)
    parser.add_argument("--overlap", type=float, default=8.0)
    parser.add_argument("--max-npts", type=int, default=126)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    if args.output.exists():
        print(f"refusing to overwrite grid plan: {args.output}", file=sys.stderr)
        return 2
    metadata = json.loads(args.receptor_metadata.read_text(encoding="utf-8"))
    if metadata.get("status") != "RECEPTOR_PASS":
        print("receptor metadata is not RECEPTOR_PASS", file=sys.stderr)
        return 3
    prepared_bbox = metadata.get("prepared_bbox")
    if not isinstance(prepared_bbox, dict):
        print("receptor metadata lacks prepared_bbox", file=sys.stderr)
        return 3

    try:
        plan = build_plan(
            prepared_bbox["min"],
            prepared_bbox["max"],
            args.spacing,
            args.margin,
            args.overlap,
            args.max_npts,
        )
    except (KeyError, TypeError, ValueError) as error:
        print(f"invalid grid plan input: {error}", file=sys.stderr)
        return 4

    plan["receptor_metadata_sha256"] = hashlib.sha256(
        args.receptor_metadata.read_bytes()
    ).hexdigest()
    plan["prepared_pdbqt_sha256"] = metadata.get("prepared_pdbqt_sha256")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as handle:
        json.dump(plan, handle, indent=2, sort_keys=True)
        handle.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
