#!/usr/bin/env python3
"""Gate entry from attack-only relaxation into proton-transfer windows."""
import argparse
import json
import math
import pathlib

import numpy as np
import parmed as pmd

THR267_OG1 = 8961
THR267_HG1 = 8962
L2_C12 = 10287
L2_O2 = 10288
L2_N3 = 10289
ATTACK_MAX_A = 2.80
OGH_MAX_A = 1.20


def lattice(box):
    a, b, c, alpha, beta, gamma = [float(v) for v in box]
    alpha, beta, gamma = map(math.radians, (alpha, beta, gamma))
    va = np.array([a, 0.0, 0.0])
    vb = np.array([b * math.cos(gamma), b * math.sin(gamma), 0.0])
    cx = c * math.cos(beta)
    cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
    cz = math.sqrt(max(0.0, c * c - cx * cx - cy * cy))
    return np.array([va, vb, [cx, cy, cz]])


def distance(structure, left, right):
    cell = lattice(structure.box)
    delta = structure.coordinates[left - 1] - structure.coordinates[right - 1]
    frac = delta @ np.linalg.inv(cell)
    frac -= np.rint(frac)
    return float(np.linalg.norm(frac @ cell))


def require(structure, number, residue, atom):
    value = structure.atoms[number - 1]
    if value.residue.name != residue or value.name != atom:
        raise ValueError(f"identity mismatch at {number}: {value.residue.name}:{value.name}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prmtop", type=pathlib.Path, required=True)
    parser.add_argument("--rst7", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    structure = pmd.load_file(str(args.prmtop), xyz=str(args.rst7))
    require(structure, THR267_OG1, "THR", "OG1")
    require(structure, THR267_HG1, "THR", "HG1")
    require(structure, L2_C12, "L2", "C12")
    require(structure, L2_O2, "L2", "O2")
    require(structure, L2_N3, "L2", "N3")
    metrics = {
        "attack_og1_c12_A": distance(structure, THR267_OG1, L2_C12),
        "thr_og1_hg1_A": distance(structure, THR267_OG1, THR267_HG1),
        "carbonyl_c12_o2_A": distance(structure, L2_C12, L2_O2),
        "amide_c12_n3_A": distance(structure, L2_C12, L2_N3),
    }
    passed = (
        metrics["attack_og1_c12_A"] <= ATTACK_MAX_A
        and metrics["thr_og1_hg1_A"] <= OGH_MAX_A
        and 1.10 <= metrics["carbonyl_c12_o2_A"] <= 1.45
        and metrics["amide_c12_n3_A"] <= 1.70
    )
    result = {
        "schema_version": 1,
        "status": "PASS_ATTACK_PREREQUISITE" if passed else "FAIL_SCIENTIFIC_ATTACK_PREREQUISITE",
        "distances_A": metrics,
        "thresholds_A": {"attack_og1_c12_max": ATTACK_MAX_A, "thr_og1_hg1_max": OGH_MAX_A, "carbonyl_c12_o2_range": [1.10, 1.45], "amide_c12_n3_max": 1.70},
        "interpretation": "Geometry gate for entering proton-transfer seed windows only; not a TS, committor, PMF, or barrier.",
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    raise SystemExit(main())
