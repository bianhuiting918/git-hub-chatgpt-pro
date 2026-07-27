#!/usr/bin/env python3
"""Independent structural and topology audit for a full-system NylC A1 preflight."""

import argparse
import importlib.util
import json
import math
import pathlib
from collections import deque

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
BUILDER_PATH = HERE / "build_nylc_a1_full_system.py"


class AuditError(RuntimeError):
    pass


def load_builder():
    spec = importlib.util.spec_from_file_location("a1_builder_for_audit", BUILDER_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def box_matrix(values):
    if len(values) == 3:
        return np.diag(values)
    if len(values) != 9:
        raise AuditError(f"unsupported GRO box length {len(values)}")
    v1x, v2y, v3z, v1y, v1z, v2x, v2z, v3x, v3y = values
    return np.array([[v1x, v2x, v3x], [v1y, v2y, v3y], [v1z, v2z, v3z]])


def local_exclusions(bonds, source, depth=3):
    graph = {}
    for first, second in bonds:
        graph.setdefault(first, set()).add(second)
        graph.setdefault(second, set()).add(first)
    distance = {source: 0}
    queue = deque([source])
    while queue:
        atom = queue.popleft()
        if distance[atom] == depth:
            continue
        for neighbor in graph.get(atom, ()):
            if neighbor not in distance:
                distance[neighbor] = distance[atom] + 1
                queue.append(neighbor)
    return set(distance)


def minimum_nonbonded(gro, chain_itp_text, transferred_global_atom_id):
    builder = load_builder()
    _, blocks = builder._sections(chain_itp_text)
    atoms = builder._atoms(blocks)
    bonds = builder._bonds(blocks)
    thr_by_name = {row["name"]: atom_id for atom_id, row in atoms.items() if row["resnr"] == 267}
    hg_local = thr_by_name["HG1"]
    chain_first = transferred_global_atom_id - hg_local + 1
    excluded_local = local_exclusions(bonds, hg_local, 3)
    excluded_global = {chain_first + atom_id - 1 for atom_id in excluded_local}

    parsed = builder.parse_gro(gro)
    coordinates = np.asarray([parsed["atoms"][i]["xyz_nm"] for i in range(1, parsed["count"] + 1)])
    target = coordinates[transferred_global_atom_id - 1]
    box = box_matrix([float(x) for x in parsed["box"].split()])
    inv_box = np.linalg.inv(box)
    delta = coordinates - target
    fractional = delta @ inv_box.T
    fractional -= np.rint(fractional)
    minimum_image = fractional @ box.T
    distances = np.linalg.norm(minimum_image, axis=1)
    for atom_id in excluded_global:
        if 1 <= atom_id <= len(distances):
            distances[atom_id - 1] = np.inf
    index = int(np.argmin(distances))
    value = float(distances[index])
    if not math.isfinite(value):
        raise AuditError("minimum nonbonded distance is not finite")
    return value, index + 1, sorted(excluded_global)


def audit(gro_text, chain_itp_text, build_audit, potential_energy, grompp_maxwarn_zero):
    if build_audit.get("status") != "PASS_A1_FULL_SYSTEM_BUILD":
        raise AuditError("full-system build audit is not PASS")
    builder = load_builder()
    _, blocks = builder._sections(chain_itp_text)
    atoms = builder._atoms(blocks)
    bonds = builder._bonds(blocks)
    thr = {i: row for i, row in atoms.items() if row["resnr"] == 267}
    by_name = {row["name"]: i for i, row in thr.items()}
    neighbor = {i: set() for i in atoms}
    for first, second in bonds:
        neighbor[first].add(second)
        neighbor[second].add(first)
    n_hydrogens = sorted(
        atoms[i]["name"] for i in neighbor[by_name["N"]] if atoms[i]["mass"] < 2.0
    )
    og_hydrogens = sorted(
        atoms[i]["name"] for i in neighbor[by_name["OG1"]] if atoms[i]["mass"] < 2.0
    )
    if n_hydrogens != ["H1", "H2", "HG1"]:
        raise AuditError(f"Nalpha hydrogen contract failed: {n_hydrogens}")
    if og_hydrogens:
        raise AuditError(f"Ogamma remains protonated: {og_hydrogens}")
    if not grompp_maxwarn_zero:
        raise AuditError("grompp was not run with -maxwarn 0")
    if not math.isfinite(potential_energy):
        raise AuditError("GROMACS potential energy is not finite")

    transferred = build_audit["coordinate_patch"]["transferred_global_atom_id"]
    min_distance, closest, excluded = minimum_nonbonded(
        gro_text, chain_itp_text, transferred
    )
    if min_distance < 0.08:
        raise AuditError(f"nonbonded overlap {min_distance:.6f} nm")
    return {
        "schema_version": 1,
        "status": "PASS_A1_FULL_SYSTEM_PREFLIGHT",
        "nalpha_hydrogen_count": len(n_hydrogens),
        "nalpha_hydrogens": n_hydrogens,
        "ogamma_hydrogen_count": len(og_hydrogens),
        "ogamma_hydrogens": og_hydrogens,
        "minimum_nonbonded_distance_nm": min_distance,
        "minimum_nonbonded_closest_global_atom_id": closest,
        "topological_exclusion_count": len(excluded),
        "gromacs_potential_energy_kj_mol": potential_energy,
        "grompp_maxwarn_zero": True,
        "scientific_scope": "technical preflight only; not an unconstrained NAC stability pass",
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gro", required=True)
    ap.add_argument("--chain-itp", required=True)
    ap.add_argument("--build-audit", required=True)
    ap.add_argument("--potential-energy-kj-mol", required=True, type=float)
    ap.add_argument("--grompp-maxwarn-zero", action="store_true")
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    result = audit(
        pathlib.Path(a.gro).read_text(),
        pathlib.Path(a.chain_itp).read_text(),
        json.loads(pathlib.Path(a.build_audit).read_text()),
        a.potential_energy_kj_mol,
        a.grompp_maxwarn_zero,
    )
    pathlib.Path(a.output).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")


if __name__ == "__main__":
    main()
