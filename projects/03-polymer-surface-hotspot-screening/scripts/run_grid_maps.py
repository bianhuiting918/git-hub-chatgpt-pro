#!/usr/bin/env python3
"""Generate and audit whole-receptor A/C/OA/HD AutoGrid maps from a frozen plan."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import subprocess
import sys
from pathlib import Path


CHANNELS = ("A", "C", "OA", "HD")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def render_gpf(
    receptor_name: str,
    parameter_name: str,
    prefix: str,
    receptor_types: list[str],
    center: list[float],
    npts: list[int],
    spacing: float,
) -> str:
    lines = [
        "outlev 1",
        f"parameter_file {parameter_name}",
        f"npts {npts[0]} {npts[1]} {npts[2]}",
        f"gridfld {prefix}.maps.fld",
        f"spacing {spacing:.6f}",
        f"receptor_types {' '.join(sorted(set(receptor_types)))}",
        "ligand_types A C OA HD",
        f"receptor {receptor_name}",
        f"gridcenter {center[0]:.6f} {center[1]:.6f} {center[2]:.6f}",
        "smooth 0.500000",
    ]
    lines.extend(f"map {prefix}.{channel}.map" for channel in CHANNELS)
    lines.extend([
        f"elecmap {prefix}.e.map",
        f"dsolvmap {prefix}.d.map",
        "dielectric -0.145600",
    ])
    return "\n".join(lines) + "\n"


def parse_map_header(path: Path) -> dict:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"missing map: {path}")
    spacing = npts = center = None
    with path.open(encoding="ascii", errors="strict") as handle:
        for _ in range(12):
            line = handle.readline()
            if not line:
                break
            fields = line.split()
            if not fields:
                continue
            if fields[0] == "SPACING":
                spacing = float(fields[1])
            elif fields[0] == "NELEMENTS":
                npts = [int(value) for value in fields[1:4]]
            elif fields[0] == "CENTER":
                center = [float(value) for value in fields[1:4]]
    if spacing is None or npts is None or center is None:
        raise ValueError(f"incomplete map header: {path}")
    return {"spacing": spacing, "npts": npts, "center": center}


def validate_map_set(paths: dict[str, Path]) -> dict:
    headers = {}
    for channel in CHANNELS:
        path = paths.get(channel)
        if path is None or not path.is_file():
            raise ValueError(f"missing map: {channel}")
        headers[channel] = parse_map_header(path)
    canonical = {
        json.dumps(headers[channel], sort_keys=True, separators=(",", ":"))
        for channel in CHANNELS
    }
    if len(canonical) != 1:
        raise ValueError("nonidentical map geometry")
    return headers["A"]


def receptor_atom_types(path: Path) -> list[str]:
    atom_types = set()
    with path.open(encoding="ascii", errors="strict") as handle:
        for line in handle:
            if line[:6] in {"ATOM  ", "HETATM"}:
                fields = line.split()
                if fields:
                    atom_types.add(fields[-1])
    if not atom_types:
        raise ValueError("receptor PDBQT has no atom types")
    return sorted(atom_types)


def run_command(command: list[str], cwd: Path, stdout_path: Path, stderr_path: Path) -> int:
    with stdout_path.open("x", encoding="utf-8") as stdout, stderr_path.open(
        "x", encoding="utf-8"
    ) as stderr:
        completed = subprocess.run(command, cwd=cwd, stdout=stdout, stderr=stderr)
    return completed.returncode


def write_json_exclusive(path: Path, payload: dict) -> None:
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--receptor-pdbqt", type=Path, required=True)
    parser.add_argument("--grid-plan", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--record-id", required=True)
    parser.add_argument("--autogrid", type=Path, required=True)
    parser.add_argument("--autosite", type=Path, required=True)
    parser.add_argument("--parameter-file", type=Path, required=True)
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path("/work/home/acshdt1dks/polymer_surface_hotspot_screen_20260725"),
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    required = [args.receptor_pdbqt, args.grid_plan, args.autogrid, args.autosite, args.parameter_file]
    if any(not path.is_file() for path in required):
        print("one or more required inputs/executables are missing", file=sys.stderr)
        return 2
    root = args.project_root.resolve()
    output = args.output_dir.resolve()
    if root != output and root not in output.parents:
        print("refusing output outside project root", file=sys.stderr)
        return 2
    if args.output_dir.exists():
        pass_gate = args.output_dir / "MAPS_PASS.json"
        if pass_gate.is_file():
            print(f"existing MAPS_PASS retained: {pass_gate}")
            return 0
        print("refusing partial or unverified map output", file=sys.stderr)
        return 3

    plan = json.loads(args.grid_plan.read_text(encoding="utf-8"))
    if plan.get("status") != "GRID_PLAN_PASS" or plan.get("channels") != list(CHANNELS):
        print("grid plan is not a compatible GRID_PLAN_PASS", file=sys.stderr)
        return 4
    if plan.get("prepared_pdbqt_sha256") != sha256(args.receptor_pdbqt):
        print("receptor checksum differs from grid plan", file=sys.stderr)
        return 4

    args.output_dir.mkdir(parents=True)
    try:
        types = receptor_atom_types(args.receptor_pdbqt)
        tile_results = []
        for tile in plan["tiles"]:
            tile_dir = args.output_dir / tile["tile_id"]
            tile_dir.mkdir()
            receptor = tile_dir / "receptor.pdbqt"
            parameter = tile_dir / "AD4.1_bound.dat"
            shutil.copyfile(args.receptor_pdbqt, receptor)
            shutil.copyfile(args.parameter_file, parameter)
            prefix = "whole_receptor"
            gpf = tile_dir / f"{prefix}.gpf"
            gpf.write_text(render_gpf(
                receptor.name, parameter.name, prefix, types,
                tile["center"], tile["npts"], tile["spacing"],
            ), encoding="ascii", newline="\n")
            autogrid_rc = run_command(
                [str(args.autogrid), "-p", gpf.name, "-l", f"{prefix}.glg"],
                tile_dir,
                tile_dir / "autogrid.stdout.log",
                tile_dir / "autogrid.stderr.log",
            )
            if autogrid_rc != 0:
                raise RuntimeError(f"autogrid failed for {tile['tile_id']} rc={autogrid_rc}")
            map_paths = {
                channel: tile_dir / f"{prefix}.{channel}.map" for channel in CHANNELS
            }
            geometry = validate_map_set(map_paths)
            expected = {
                "spacing": float(tile["spacing"]),
                "npts": list(tile["npts"]),
                "center": [round(float(value), 3) for value in tile["center"]],
            }
            observed = {
                "spacing": geometry["spacing"],
                "npts": geometry["npts"],
                "center": [round(float(value), 3) for value in geometry["center"]],
            }
            if observed != expected:
                raise ValueError(f"map geometry differs from plan for {tile['tile_id']}")
            boxdim = [float(n) * float(tile["spacing"]) for n in tile["npts"]]
            autosite_dir = tile_dir / "autosite"
            autosite_rc = run_command(
                [
                    str(args.autosite), "-r", str(receptor), "-o", str(autosite_dir),
                    "-n", "10", "--spacing", str(tile["spacing"]),
                    "--boxcenter", *[str(value) for value in tile["center"]],
                    "--boxdim", *[str(value) for value in boxdim],
                ],
                tile_dir,
                tile_dir / "autosite.stdout.log",
                tile_dir / "autosite.stderr.log",
            )
            if autosite_rc != 0:
                raise RuntimeError(f"autosite failed for {tile['tile_id']} rc={autosite_rc}")
            tile_results.append({
                "tile_id": tile["tile_id"],
                "geometry": geometry,
                "map_sha256": {channel: sha256(path) for channel, path in map_paths.items()},
                "autogrid_exit_code": autogrid_rc,
                "autosite_exit_code": autosite_rc,
            })
        payload = {
            "status": "MAPS_PASS",
            "scientific_scope": "map_generation_only",
            "record_id": args.record_id,
            "receptor_pdbqt_sha256": sha256(args.receptor_pdbqt),
            "grid_plan_sha256": sha256(args.grid_plan),
            "channels": list(CHANNELS),
            "tile_count": len(tile_results),
            "tiles": tile_results,
        }
        write_json_exclusive(args.output_dir / "MAPS_PASS.json", payload)
        checksums = []
        for path in sorted(args.output_dir.rglob("*")):
            if path.is_file() and path.name != "SHA256SUMS":
                checksums.append(f"{sha256(path)}  {path.relative_to(args.output_dir)}")
        (args.output_dir / "SHA256SUMS").write_text("\n".join(checksums) + "\n", encoding="ascii")
    except Exception as error:
        write_json_exclusive(args.output_dir / "MAPS_FAIL.json", {
            "status": "NOT_EVALUATED_GRID",
            "record_id": args.record_id,
            "error": str(error),
        })
        print(str(error), file=sys.stderr)
        return 5
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
