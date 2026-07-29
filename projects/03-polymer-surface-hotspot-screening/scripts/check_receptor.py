#!/usr/bin/env python3
"""Deterministically select PDB receptor context and audit prepared PDBQT output."""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
from collections import defaultdict
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_csv(value: str) -> list[str]:
    return list(dict.fromkeys(part.strip() for part in value.split(",") if part.strip()))


def coordinate(line: str) -> tuple[float, float, float]:
    return float(line[30:38]), float(line[38:46]), float(line[46:54])


def residue_token(line: str) -> str:
    chain = line[21].strip() or "_"
    number = line[22:26].strip()
    insertion = line[26].strip()
    return f"{chain}:{number}{insertion}"


def hetero_token(line: str) -> str:
    chain = line[21].strip() or "_"
    residue = line[17:20].strip()
    number = line[22:26].strip()
    insertion = line[26].strip()
    return f"{chain}:{residue}:{number}{insertion}"


def atom_key(line: str) -> tuple[str, str, str, str, str]:
    return (
        line[:6],
        line[21],
        line[22:26],
        line[26],
        line[12:16],
    )


def occupancy(line: str) -> float:
    try:
        return float(line[54:60])
    except ValueError:
        return 0.0


def altloc_rank(line: str) -> tuple[int, float, str]:
    altloc = line[16].strip()
    preference = 2 if altloc == "" else 1 if altloc == "A" else 0
    return preference, occupancy(line), "" if altloc == "" else altloc


def write_json_exclusive(path: Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)
        handle.write("\n")


def parse_pdb(path: Path, selected_model: int, chains: set[str]):
    current_model = 1
    saw_model = False
    model_ids: set[int] = set()
    candidates: list[tuple[int, str]] = []
    found_chains: set[str] = set()
    all_hetero: set[str] = set()

    with path.open(encoding="ascii", errors="strict") as handle:
        for order, raw in enumerate(handle):
            line = raw.rstrip("\r\n")
            record = line[:6]
            if record == "MODEL ":
                saw_model = True
                try:
                    current_model = int(line[10:14].strip())
                except ValueError:
                    current_model = int(line.split()[1])
                model_ids.add(current_model)
                continue
            if record == "ENDMDL":
                continue
            if record not in {"ATOM  ", "HETATM"}:
                continue
            if not saw_model:
                model_ids.add(1)
                current_model = 1
            if current_model != selected_model:
                continue
            chain = line[21].strip() or "_"
            if chain not in chains:
                continue
            found_chains.add(chain)
            if record == "HETATM":
                all_hetero.add(hetero_token(line))
            candidates.append((order, line))

    return candidates, found_chains, sorted(model_ids), all_hetero


def select_records(candidates, retained_hetero: set[str]):
    grouped: dict[tuple[str, str, str, str, str], list[tuple[int, str]]] = defaultdict(list)
    for order, line in candidates:
        if line[:6] == "HETATM" and hetero_token(line) not in retained_hetero:
            continue
        grouped[atom_key(line)].append((order, line))

    chosen: list[tuple[int, str]] = []
    altloc_warnings: list[str] = []
    for key, options in grouped.items():
        best_order, best_line = max(
            options,
            key=lambda item: (altloc_rank(item[1]), -item[0]),
        )
        if len(options) > 1:
            altloc_warnings.append(
                f"{residue_token(best_line)}:{best_line[12:16].strip()}"
            )
        normalized = best_line[:16] + " " + best_line[17:]
        chosen.append((best_order, normalized))
    chosen.sort(key=lambda item: item[0])
    return [line for _, line in chosen], sorted(set(altloc_warnings))


def bbox(lines: list[str]) -> dict[str, list[float]]:
    points = [coordinate(line) for line in lines]
    return {
        "min": [min(point[i] for point in points) for i in range(3)],
        "max": [max(point[i] for point in points) for i in range(3)],
    }


def select_command(args: argparse.Namespace) -> int:
    if args.output.exists() or args.metadata.exists():
        print("refusing to overwrite receptor selection output", file=sys.stderr)
        return 2

    chains = set(parse_csv(args.chains))
    catalytic = set(parse_csv(args.catalytic_residues))
    retained_requested = set(parse_csv(args.retain_hetero or ""))
    if not chains:
        print("at least one chain is required", file=sys.stderr)
        return 2

    candidates, found_chains, model_ids, all_hetero = parse_pdb(
        args.input,
        args.model,
        chains,
    )
    missing_chains = sorted(chains.difference(found_chains))
    base = {
        "input_sha256": sha256(args.input),
        "selected_model": args.model,
        "available_models": model_ids,
        "selected_chains": sorted(chains),
        "alternate_location_policy": "blank_then_A_then_occupancy",
        "retained_hetero": sorted(all_hetero.intersection(retained_requested)),
        "removed_hetero": sorted(all_hetero.difference(retained_requested)),
    }
    if missing_chains:
        write_json_exclusive(
            args.metadata,
            {
                **base,
                "status": "NOT_EVALUATED_CHAIN_MAPPING",
                "missing_chains": missing_chains,
            },
        )
        return 3
    if not candidates:
        write_json_exclusive(
            args.metadata,
            {**base, "status": "NOT_EVALUATED_STRUCTURE_PARSE"},
        )
        return 3

    selected, altloc_warnings = select_records(candidates, retained_requested)
    protein_residues = {
        residue_token(line) for line in selected if line[:6] == "ATOM  "
    }
    missing_catalytic = sorted(catalytic.difference(protein_residues))
    if missing_catalytic:
        write_json_exclusive(
            args.metadata,
            {
                **base,
                "status": "NOT_EVALUATED_CATALYTIC_MAPPING",
                "catalytic_mapping_status": "NOT_EVALUATED",
                "missing_catalytic_residues": missing_catalytic,
            },
        )
        return 4

    args.output.parent.mkdir(parents=True, exist_ok=True)
    header = [
        f"REMARK 999 INPUT_SHA256 {base['input_sha256']}\n",
        f"REMARK 999 SELECTED_MODEL {args.model}\n",
        f"REMARK 999 SELECTED_CHAINS {','.join(sorted(chains))}\n",
    ]
    with args.output.open("x", encoding="ascii", newline="\n") as handle:
        handle.writelines(header)
        for line in selected:
            handle.write(line + "\n")
        handle.write("END\n")

    atom_count = sum(line[:6] == "ATOM  " for line in selected)
    hetero_count = sum(line[:6] == "HETATM" for line in selected)
    write_json_exclusive(
        args.metadata,
        {
            **base,
            "status": "SELECTION_PASS",
            "catalytic_mapping_status": "PASS",
            "missing_catalytic_residues": [],
            "selected_atom_count": atom_count,
            "selected_hetero_atom_count": hetero_count,
            "selected_coordinate_count": len(selected),
            "selected_bbox": bbox(selected),
            "selected_pdb_sha256": sha256(args.output),
            "warnings": [
                f"ALTLOC_RESOLVED:{token}" for token in altloc_warnings
            ],
        },
    )
    return 0


def finalize_command(args: argparse.Namespace) -> int:
    if args.output_metadata.exists():
        print("refusing to overwrite final receptor metadata", file=sys.stderr)
        return 2
    selection = json.loads(args.selection_metadata.read_text(encoding="utf-8"))
    if selection.get("status") != "SELECTION_PASS":
        print("selection metadata is not PASS", file=sys.stderr)
        return 3

    coordinates = []
    with args.prepared_pdbqt.open(encoding="ascii", errors="strict") as handle:
        for raw in handle:
            if raw[:6] in {"ATOM  ", "HETATM"}:
                coordinates.append(raw.rstrip("\r\n"))
    if not coordinates:
        print("prepared PDBQT has no coordinate records", file=sys.stderr)
        return 4

    payload = {
        **selection,
        "status": "RECEPTOR_PASS",
        "prepared_pdbqt_sha256": sha256(args.prepared_pdbqt),
        "prepared_atom_count": len(coordinates),
        "prepared_bbox": bbox(coordinates),
    }
    write_json_exclusive(args.output_metadata, payload)
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest="command", required=True)

    select = subparsers.add_parser("select")
    select.add_argument("--input", type=Path, required=True)
    select.add_argument("--output", type=Path, required=True)
    select.add_argument("--metadata", type=Path, required=True)
    select.add_argument("--chains", required=True)
    select.add_argument("--model", type=int, default=1)
    select.add_argument("--catalytic-residues", required=True)
    select.add_argument("--retain-hetero", default="")
    select.set_defaults(func=select_command)

    finalize = subparsers.add_parser("finalize")
    finalize.add_argument("--selection-metadata", type=Path, required=True)
    finalize.add_argument("--prepared-pdbqt", type=Path, required=True)
    finalize.add_argument("--output-metadata", type=Path, required=True)
    finalize.set_defaults(func=finalize_command)
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    return args.func(args)


if __name__ == "__main__":
    raise SystemExit(main())
