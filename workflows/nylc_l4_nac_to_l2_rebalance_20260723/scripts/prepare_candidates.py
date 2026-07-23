#!/usr/bin/env python3
"""Prepare five independent audited PA66-L2 systems from immutable L4 NAC frames."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import re
import shutil
from pathlib import Path
from typing import Any, Dict, Iterable, Tuple

from build_l4_to_l2 import (
    _minimum_image_vector,
    build_l2_coordinates,
    discover_heavy_mapping,
    extract_ligand_coordinates,
    minimum_ligand_environment_distance,
    parse_gro,
    parse_itp,
    replace_ligand,
    write_gro,
)

Vector = Tuple[float, float, float]
POSRE_LEVELS = (1000, 500, 100, 10)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def make_posre(path: Path, force_constant: int) -> None:
    lines = [
        "; PA66_L2 heavy-atom position restraints; target atoms 1-33 are heavy",
        "[ position_restraints ]",
        "; atom  type      fx      fy      fz",
    ]
    lines.extend(
        f"{atom_id:6d}     1 {force_constant:7d} {force_constant:7d} {force_constant:7d}"
        for atom_id in range(1, 34)
    )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _l2_include_contract() -> str:
    lines = ['#include "PA66_L2_GMX.itp"']
    for force_constant in POSRE_LEVELS:
        lines.extend(
            [
                f"#ifdef POSRES_L2_{force_constant}",
                f'#include "posre_l2_{force_constant}.itp"',
                "#endif",
            ]
        )
    return "\n".join(lines)


def rewrite_topology(source: str, source_itp_name: str, source_molecule_type: str) -> str:
    old_blocks = re.compile(
        r"(?ms)^\s*#ifdef\s+POSRES_LIG(?:_WEAK)?\s*$.*?^\s*#endif\s*$\n?"
    )
    rewritten = old_blocks.sub("", source)
    include_pattern = re.compile(
        r'(?m)^\s*#include\s+"' + re.escape(source_itp_name) + r'"\s*$'
    )
    rewritten, include_count = include_pattern.subn(_l2_include_contract(), rewritten)
    if include_count != 1:
        raise ValueError(
            f"Expected one include for {source_itp_name}, found {include_count}"
        )
    molecule_pattern = re.compile(
        r"(?m)^(\s*)" + re.escape(source_molecule_type) + r"(\s+)1(\s*(?:;.*)?)$"
    )
    rewritten, molecule_count = molecule_pattern.subn(
        r"\1PA66_L2\g<2>1\g<3>", rewritten
    )
    if molecule_count != 1:
        raise ValueError(
            f"Expected one molecule row for {source_molecule_type}, found {molecule_count}"
        )
    return rewritten.rstrip() + "\n"


def _angle_degrees(left: Vector, right: Vector) -> float:
    left_norm = math.sqrt(sum(value * value for value in left))
    right_norm = math.sqrt(sum(value * value for value in right))
    if left_norm == 0.0 or right_norm == 0.0:
        raise ValueError("Cannot calculate angle using a zero-length vector")
    cosine = sum(a * b for a, b in zip(left, right)) / (left_norm * right_norm)
    cosine = max(-1.0, min(1.0, cosine))
    return math.degrees(math.acos(cosine))


def _reaction_geometry(
    system: Any, carbonyl_global: int, amide_global: int, thr_global: int
) -> Dict[str, float]:
    carbon = system.atoms[carbonyl_global - 1].coordinate
    amide = system.atoms[amide_global - 1].coordinate
    thr = system.atoms[thr_global - 1].coordinate
    carbon_to_thr = _minimum_image_vector(
        tuple(a - b for a, b in zip(thr, carbon)), system.box
    )
    carbon_to_amide = _minimum_image_vector(
        tuple(a - b for a, b in zip(amide, carbon)), system.box
    )
    return {
        "distance_thr_og1_to_carbonyl_c_nm": math.sqrt(
            sum(value * value for value in carbon_to_thr)
        ),
        "angle_thr_og1_carbonyl_c_amide_n_deg": _angle_degrees(
            carbon_to_thr, carbon_to_amide
        ),
    }


def _copy_topology_inputs(candidate: Dict[str, Any], build_dir: Path) -> None:
    branch = Path(candidate["branch_dir"])
    source_ligand = candidate["source_itp"]
    for source in sorted(branch.glob("*.itp")):
        if source.name == source_ligand or source.name.startswith("posre_P66"):
            continue
        shutil.copy2(source, build_dir / source.name)
    shutil.copy2(branch / "cycle.ndx", build_dir / "source_cycle.ndx")


def _mapping_payload(
    source: Any, target: Any, mapping: Dict[int, int | None], heavy: Any
) -> Dict[str, Any]:
    target_to_source = []
    source_to_l2 = []
    for target_id, source_id in sorted(mapping.items()):
        row = {
            "l2_atom_id": target_id,
            "l2_atom_name": target.atoms[target_id].atom_name,
            "element": target.atoms[target_id].element,
            "source_atom_id": source_id,
            "source_atom_name": (
                source.atoms[source_id].atom_name if source_id is not None else None
            ),
            "coordinate_origin": "source_exact" if source_id is not None else "rebuilt_endpoint",
        }
        target_to_source.append(row)
        if source_id is not None:
            source_to_l2.append(
                {
                    "source_atom_id": source_id,
                    "source_atom_name": source.atoms[source_id].atom_name,
                    "l2_atom_id": target_id,
                    "l2_atom_name": target.atoms[target_id].atom_name,
                }
            )
    return {
        "target_to_source": target_to_source,
        "source_to_l2": sorted(source_to_l2, key=lambda row: row["source_atom_id"]),
        "cut_edges": [
            {
                "l2_endpoint": cut.target_endpoint,
                "source_endpoint": cut.source_endpoint,
                "removed_source_neighbor": cut.removed_source_neighbor,
            }
            for cut in heavy.cut_edges
        ],
        "target_reactive_triplet": {
            "carbonyl_c": heavy.target_reactive_triplet[0],
            "carbonyl_o": heavy.target_reactive_triplet[1],
            "amide_n": heavy.target_reactive_triplet[2],
        },
    }


def prepare_one(
    candidate: Dict[str, Any], task_root: Path, l2_itp_path: Path
) -> Dict[str, Any]:
    candidate_root = task_root / "candidates" / candidate["id"]
    build_dir = candidate_root / "build"
    audit_dir = candidate_root / "audit"
    build_dir.mkdir(parents=True, exist_ok=True)
    audit_dir.mkdir(parents=True, exist_ok=True)

    source_gro_path = task_root / "inputs" / candidate["id"] / "source.gro"
    source_itp_path = Path(candidate["branch_dir"]) / candidate["source_itp"]
    source_top_path = Path(candidate["branch_dir"]) / "topol.top"
    source_topology = parse_itp(source_itp_path)
    target_topology = parse_itp(l2_itp_path)
    source_system = parse_gro(source_gro_path)

    reactive = candidate["reactive_local_atoms"]
    heavy = discover_heavy_mapping(
        source_topology,
        target_topology,
        (reactive["carbonyl_c"], reactive["carbonyl_o"], reactive["amide_n"]),
    )
    source_coordinates = extract_ligand_coordinates(
        source_system,
        candidate["ligand_first_global_atom"],
        candidate["source_atom_count"],
    )
    target_coordinates, atom_mapping = build_l2_coordinates(
        source_topology, target_topology, heavy, source_coordinates
    )
    rebuilt = replace_ligand(
        source_system,
        candidate["ligand_first_global_atom"],
        candidate["source_atom_count"],
        target_topology,
        target_coordinates,
    )
    rebuilt_path = build_dir / "rebuilt.gro"
    write_gro(rebuilt, rebuilt_path)
    rebuilt_roundtrip = parse_gro(rebuilt_path)

    _copy_topology_inputs(candidate, build_dir)
    shutil.copy2(l2_itp_path, build_dir / "PA66_L2_GMX.itp")
    for force_constant in POSRE_LEVELS:
        make_posre(build_dir / f"posre_l2_{force_constant}.itp", force_constant)
    rewritten_topology = rewrite_topology(
        source_top_path.read_text(encoding="utf-8"),
        candidate["source_itp"],
        candidate["source_molecule_type"],
    )
    (build_dir / "topol.top").write_text(rewritten_topology, encoding="utf-8")

    mapping_payload = _mapping_payload(
        source_topology, target_topology, atom_mapping, heavy
    )
    mapping_path = audit_dir / "atom_mapping.json"
    mapping_path.write_text(json.dumps(mapping_payload, indent=2) + "\n", encoding="utf-8")

    target_c, target_o, target_n = heavy.target_reactive_triplet
    first = candidate["ligand_first_global_atom"]
    target_globals = {
        "carbonyl_c": first + target_c - 1,
        "carbonyl_o": first + target_o - 1,
        "amide_n": first + target_n - 1,
        "thr_og1": candidate["reactive_global_atoms"]["thr_og1"],
    }
    preservation = {}
    max_delta = 0.0
    for role, target_id in (
        ("carbonyl_c", target_c),
        ("carbonyl_o", target_o),
        ("amide_n", target_n),
    ):
        source_id = heavy.target_to_source[target_id]
        source_coord = source_coordinates[source_id]
        written_coord = rebuilt_roundtrip.atoms[first + target_id - 2].coordinate
        delta = max(abs(a - b) for a, b in zip(source_coord, written_coord))
        max_delta = max(max_delta, delta)
        preservation[role] = {
            "source_local_atom": source_id,
            "l2_local_atom": target_id,
            "source_coordinate_nm": source_coord,
            "rebuilt_coordinate_nm": written_coord,
            "max_abs_delta_nm": delta,
        }

    atom_count_expected = (
        len(source_system.atoms) - candidate["source_atom_count"] + len(target_topology.atoms)
    )
    charge = sum(atom.charge for atom in target_topology.atoms.values())
    min_environment = minimum_ligand_environment_distance(
        rebuilt_roundtrip, first, len(target_topology.atoms)
    )
    geometry = _reaction_geometry(
        rebuilt_roundtrip,
        target_globals["carbonyl_c"],
        target_globals["amide_n"],
        target_globals["thr_og1"],
    )
    checks = {
        "l2_atom_count_79": len(target_topology.atoms) == 79,
        "l2_heavy_atom_count_33": sum(
            atom.element != "H" for atom in target_topology.atoms.values()
        ) == 33,
        "l2_charge_neutral": abs(charge) < 1e-5,
        "two_cut_edges": len(heavy.cut_edges) == 2,
        "whole_system_atom_count": len(rebuilt_roundtrip.atoms) == atom_count_expected,
        "reactive_coordinates_preserved_to_gro_precision": max_delta <= 5e-4,
        "minimum_ligand_environment_distance_gt_0p05_nm": min_environment > 0.05,
        "reaction_geometry_finite": all(math.isfinite(value) for value in geometry.values()),
    }
    audit = {
        "schema_version": 1,
        "candidate_id": candidate["id"],
        "state": "BUILD_PASS" if all(checks.values()) else "BUILD_FAIL",
        "source": {
            "gro": str(source_gro_path),
            "gro_sha256": sha256(source_gro_path),
            "itp": str(source_itp_path),
            "itp_sha256": sha256(source_itp_path),
            "topology": str(source_top_path),
            "topology_sha256": sha256(source_top_path),
            "atom_count": len(source_system.atoms),
            "molecule_atom_count": candidate["source_atom_count"],
        },
        "target": {
            "gro": str(rebuilt_path),
            "gro_sha256": sha256(rebuilt_path),
            "itp": str(build_dir / "PA66_L2_GMX.itp"),
            "itp_sha256": sha256(build_dir / "PA66_L2_GMX.itp"),
            "topology": str(build_dir / "topol.top"),
            "topology_sha256": sha256(build_dir / "topol.top"),
            "molecule_type": target_topology.molecule_type,
            "atom_count": len(rebuilt_roundtrip.atoms),
            "molecule_atom_count": len(target_topology.atoms),
            "charge_e": charge,
        },
        "atom_count_expected": atom_count_expected,
        "reactive_global_atoms": target_globals,
        "reactive_coordinate_preservation": preservation,
        "max_reactive_coordinate_delta_nm": max_delta,
        "reaction_geometry": geometry,
        "source_reported_geometry": candidate["source_geometry"],
        "minimum_ligand_environment_distance_nm": min_environment,
        "mapping_file": str(mapping_path),
        "checks": checks,
    }
    audit_path = audit_dir / "build_audit.json"
    audit_path.write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    return audit


def prepare_all(manifest_path: Path, task_root: Path) -> Dict[str, Any]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("formal_candidate_count") != 5 or len(manifest.get("candidates", [])) != 5:
        raise ValueError("Formal candidate universe must contain exactly five candidates")
    l2 = manifest["audited_l2_itp"]
    l2_itp_path = Path(l2["path"])
    if sha256(l2_itp_path) != l2["sha256"]:
        raise ValueError("Audited L2 ITP SHA256 mismatch")

    rows = []
    for candidate in manifest["candidates"]:
        try:
            audit = prepare_one(candidate, task_root, l2_itp_path)
            rows.append(
                {
                    "candidate_id": candidate["id"],
                    "state": audit["state"],
                    "reason": "" if audit["state"] == "BUILD_PASS" else "one_or_more_build_checks_failed",
                    "audit": str(task_root / "candidates" / candidate["id"] / "audit" / "build_audit.json"),
                }
            )
        except Exception as exc:
            audit_dir = task_root / "candidates" / candidate["id"] / "audit"
            audit_dir.mkdir(parents=True, exist_ok=True)
            failure = {
                "schema_version": 1,
                "candidate_id": candidate["id"],
                "state": "BUILD_FAIL",
                "exception_type": type(exc).__name__,
                "reason": str(exc),
            }
            (audit_dir / "build_audit.json").write_text(
                json.dumps(failure, indent=2) + "\n", encoding="utf-8"
            )
            rows.append(
                {
                    "candidate_id": candidate["id"],
                    "state": "BUILD_FAIL",
                    "reason": f"{type(exc).__name__}: {exc}",
                    "audit": str(audit_dir / "build_audit.json"),
                }
            )

    summary_path = task_root / "audit" / "prepare_summary.tsv"
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["candidate_id\tstate\treason\taudit"]
    lines.extend(
        "\t".join(
            str(row[key]).replace("\t", " ").replace("\n", " ")
            for key in ("candidate_id", "state", "reason", "audit")
        )
        for row in rows
    )
    summary_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return {
        "formal_candidate_count": 5,
        "build_pass_count": sum(row["state"] == "BUILD_PASS" for row in rows),
        "rows": rows,
        "summary": str(summary_path),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--task-root", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(prepare_all(args.manifest, args.task_root), indent=2))


if __name__ == "__main__":
    main()
