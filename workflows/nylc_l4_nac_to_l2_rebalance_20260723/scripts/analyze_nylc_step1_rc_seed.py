#!/usr/bin/env python3
"""Audit NylC Step1 reaction-coordinate seed geometries in two microstates."""
import argparse
import json
import math
import pathlib

import numpy as np
import parmed as pmd

DEPROT_SOURCE = "step1_dftb3_preflight_post61710861_job61712692"
ASH_SOURCE = "ash306_full_system_preflight_job_61718715"
THR267_N = 8949
THR267_OG1 = 8961
THR267_HG1 = 8962
L2_C12_DEPROT = 10287
L2_O2_DEPROT = 10288
L2_N3_DEPROT = 10289
ASP306_OD1 = 9572
ASP306_OD2 = 9573
ASP308_OD1_DEPROT = 9591
ASP308_OD2_DEPROT = 9592


def cell_from_box(box):
    a, b, c, alpha, beta, gamma = [float(value) for value in box]
    alpha, beta, gamma = map(math.radians, (alpha, beta, gamma))
    va = np.array([a, 0.0, 0.0])
    vb = np.array([b * math.cos(gamma), b * math.sin(gamma), 0.0])
    cx = c * math.cos(beta)
    cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
    cz = math.sqrt(max(0.0, c * c - cx * cx - cy * cy))
    return np.array([va, vb, [cx, cy, cz]])


def minimum_image_vector(left, right, cell):
    delta = np.asarray(left) - np.asarray(right)
    fractional = delta @ np.linalg.inv(cell)
    fractional -= np.rint(fractional)
    shifts = np.asarray([(i, j, k) for i in (-1.0, 0.0, 1.0) for j in (-1.0, 0.0, 1.0) for k in (-1.0, 0.0, 1.0)])
    candidates = (fractional[None, :] + shifts) @ cell
    return candidates[np.argmin(np.linalg.norm(candidates, axis=1))]


def distance(structure, left, right):
    cell = cell_from_box(structure.box)
    vector = minimum_image_vector(structure.atoms[left - 1].xx, structure.atoms[right - 1].xx, cell)
    return float(np.linalg.norm(vector))


def identity(structure, number, residue, atom):
    value = structure.atoms[number - 1]
    if value.residue.name != residue or value.name != atom:
        raise ValueError(f"identity mismatch at {number}: {value.residue.name}:{value.name}")
    return number


def load_system(prmtop, rst7):
    structure = pmd.load_file(str(prmtop), xyz=str(rst7))
    if structure.box is None:
        raise ValueError("periodic box is missing")
    return structure


def metrics_deprotonated(structure):
    atoms = {
        "thr_n": identity(structure, THR267_N, "THR", "N"),
        "thr_og1": identity(structure, THR267_OG1, "THR", "OG1"),
        "thr_hg1": identity(structure, THR267_HG1, "THR", "HG1"),
        "l2_c12": identity(structure, L2_C12_DEPROT, "L2", "C12"),
        "l2_o2": identity(structure, L2_O2_DEPROT, "L2", "O2"),
        "l2_n3": identity(structure, L2_N3_DEPROT, "L2", "N3"),
        "asp306_od1": identity(structure, ASP306_OD1, "ASP", "OD1"),
        "asp306_od2": identity(structure, ASP306_OD2, "ASP", "OD2"),
        "asp308_od1": identity(structure, ASP308_OD1_DEPROT, "ASP", "OD1"),
        "asp308_od2": identity(structure, ASP308_OD2_DEPROT, "ASP", "OD2"),
    }
    d = lambda left, right: distance(structure, atoms[left], atoms[right])
    oh = d("thr_og1", "thr_hg1")
    return {
        "atom_mapping": atoms,
        "distances_A": {
            "attack_og1_to_l2_c12": d("thr_og1", "l2_c12"),
            "thr_og1_hg1": oh,
            "thr_hg1_to_thr_n": d("thr_hg1", "thr_n"),
            "thr_hg1_to_asp306_od1": d("thr_hg1", "asp306_od1"),
            "thr_hg1_to_asp306_od2": d("thr_hg1", "asp306_od2"),
            "carbonyl_c12_o2": d("l2_c12", "l2_o2"),
            "amide_c12_n3": d("l2_c12", "l2_n3"),
            "asp306_od1_to_asp308_od1": d("asp306_od1", "asp308_od1"),
            "asp306_od2_to_asp308_od2": d("asp306_od2", "asp308_od2"),
        },
        "candidate_cvs_A": {
            "attack_og1_to_l2_c12": d("thr_og1", "l2_c12"),
            "pt_og1_hg1_to_asp306_od1": oh - d("thr_hg1", "asp306_od1"),
            "pt_og1_hg1_to_asp306_od2": oh - d("thr_hg1", "asp306_od2"),
            "carbonyl_c12_o2": d("l2_c12", "l2_o2"),
            "amide_c12_n3": d("l2_c12", "l2_n3"),
        },
    }


def metrics_ash306(structure):
    mapping = {
        "thr_og1": identity(structure, 8961, "THR", "OG1"),
        "thr_hg1": identity(structure, 8962, "THR", "HG1"),
        "ash306_od1": identity(structure, 9572, "ASH", "OD1"),
        "ash306_od2": identity(structure, 9573, "ASH", "OD2"),
        "ash306_hd2": identity(structure, 9574, "ASH", "HD2"),
        "asp308_od1": identity(structure, 9592, "ASP", "OD1"),
        "asp308_od2": identity(structure, 9593, "ASP", "OD2"),
        "l2_c12": identity(structure, 10288, "L2", "C12"),
        "l2_o2": identity(structure, 10289, "L2", "O2"),
        "l2_n3": identity(structure, 10290, "L2", "N3"),
        "bridge_water_ow": identity(structure, 81219, "SOL", "OW"),
    }
    d = lambda left, right: distance(structure, mapping[left], mapping[right])
    return {
        "atom_mapping": mapping,
        "distances_A": {
            "attack_og1_to_l2_c12": d("thr_og1", "l2_c12"),
            "thr_og1_hg1": d("thr_og1", "thr_hg1"),
            "ash306_od2_hd2": d("ash306_od2", "ash306_hd2"),
            "ash306_hd2_to_bridge_water_ow": d("ash306_hd2", "bridge_water_ow"),
            "ash306_hd2_to_asp308_od1": d("ash306_hd2", "asp308_od1"),
            "ash306_hd2_to_asp308_od2": d("ash306_hd2", "asp308_od2"),
            "carbonyl_c12_o2": d("l2_c12", "l2_o2"),
            "amide_c12_n3": d("l2_c12", "l2_n3"),
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--task-root", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    task = args.task_root.resolve()
    deprot = task / "candidates/nylc_C18_trueT267_freeGS/qmmm" / DEPROT_SOURCE
    ash = task / "candidates/nylc_C18_trueT267_freeGS/qmmm" / ASH_SOURCE
    deprot_pass = json.loads((deprot / "PASS.json").read_text())
    ash_pass = json.loads((ash / "PASS.json").read_text())
    if deprot_pass.get("status") != "PASS_DFTB3_PREFLIGHT":
        raise ValueError("all-deprotonated source is not PASS_DFTB3_PREFLIGHT")
    if ash_pass.get("status") != "PASS_ASH306_FULL_SYSTEM_PREFLIGHT":
        raise ValueError("ASH306 source is not PASS_ASH306_FULL_SYSTEM_PREFLIGHT")
    deprot_structure = load_system(deprot / "system.prmtop", deprot / "selected_free_l2_nac.rst7")
    ash_structure = load_system(ash / "amber/system.prmtop", ash / "amber/system.rst7")
    deprot_metrics = metrics_deprotonated(deprot_structure)
    ash_metrics = metrics_ash306(ash_structure)
    d = deprot_metrics["distances_A"]
    closer = "OD1" if d["thr_hg1_to_asp306_od1"] <= d["thr_hg1_to_asp306_od2"] else "OD2"
    result = {
        "schema_version": 1,
        "candidate_id": "nylc_C18_trueT267_freeGS",
        "status": "PASS_STEP1_RC_SEED_AUDIT",
        "all_deprotonated_reactant_candidate": deprot_metrics,
        "ash306_protonated_reference": ash_metrics,
        "geometrically_closer_asp306_oxygen_to_thr_hg1": closer,
        "candidate basin thresholds": {
            "reactant_seed_A": {
                "attack_og1_c12_gt": 2.6,
                "thr_og1_hg1_lt": 1.2,
                "asp306_o_hg1_gt": 1.3,
                "amide_c12_n3_lt": 1.7,
            },
            "tetrahedral_seed_A": {
                "attack_og1_c12_lt": 1.8,
                "thr_og1_hg1_gt": 1.3,
                "asp306_o_hg1_lt": 1.2,
                "amide_c12_n3_lt": 1.7,
            },
            "scope": "find_ts bracketing seeds only; both Asp306 OD1 and OD2 hypotheses remain active until constrained-path evidence resolves them",
        },
        "recommended_next": "Run short, independently audited constrained DFTB3 minimization brackets for attack plus proton transfer to OD1 and OD2 in the all-deprotonated core region.",
        "interpretation": "Geometry and candidate basin thresholds only; not a validated RC, TS, committor, PMF, or barrier.",
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
