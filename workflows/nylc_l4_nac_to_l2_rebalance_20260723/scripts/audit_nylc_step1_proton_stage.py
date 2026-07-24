#!/usr/bin/env python3
"""Audit one staged Step1 proton-bracket restart."""
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
DETACHED_MAX_A = 1.35


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=pathlib.Path, required=True)
    parser.add_argument("--prmtop", type=pathlib.Path, required=True)
    parser.add_argument("--rst7", type=pathlib.Path, required=True)
    parser.add_argument("--stage", choices=("p00", "p01", "p02", "p03"), required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    spec = next(value for value in manifest["schedule"] if value["name"] == args.stage)
    structure = pmd.load_file(str(args.prmtop), xyz=str(args.rst7))
    acceptor = int(manifest["acceptor_atom"])
    attack = distance(structure, THR267_OG1, L2_C12)
    ogh = distance(structure, THR267_OG1, THR267_HG1)
    hacc = distance(structure, THR267_HG1, acceptor)
    co = distance(structure, L2_C12, L2_O2)
    cn = distance(structure, L2_C12, L2_N3)
    proton_attached = min(ogh, hacc) <= DETACHED_MAX_A
    chemistry_valid = 1.10 <= co <= 1.55 and cn <= 1.70
    progress = attack <= spec["attack_gate_A"] and hacc <= spec["hacc_gate_A"]
    endpoint = args.stage == "p03" and attack <= 1.95 and hacc <= 1.30 and ogh >= 1.20
    passed = proton_attached and chemistry_valid and (endpoint if args.stage == "p03" else progress)
    if passed and args.stage == "p03":
        status = "PASS_SCIENTIFIC_TETRAHEDRAL_SEED_REACHED"
    elif passed:
        status = "PASS_PROTON_BRACKET_STAGE"
    else:
        status = "FAIL_SCIENTIFIC_PROTON_BRACKET_STAGE"
    result = {
        "schema_version": 1,
        "status": status,
        "candidate_id": manifest["candidate_id"],
        "acceptor_hypothesis": manifest["acceptor_hypothesis"],
        "stage": args.stage,
        "distances_A": {"attack_og1_c12": attack, "thr_og1_hg1": ogh, "hg1_acceptor": hacc, "carbonyl_c12_o2": co, "amide_c12_n3": cn},
        "gates": {"proton_attached": proton_attached, "chemistry_valid": chemistry_valid, "stage_progress": progress, "endpoint": endpoint},
        "thresholds_A": {"detached_proton_max_nearest_bond": DETACHED_MAX_A, "attack_stage_max": spec["attack_gate_A"], "hacc_stage_max": spec["hacc_gate_A"]},
        "interpretation": "A passing stage is a constrained bracket seed only; not a TS, committor, PMF, or barrier.",
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    raise SystemExit(main())
