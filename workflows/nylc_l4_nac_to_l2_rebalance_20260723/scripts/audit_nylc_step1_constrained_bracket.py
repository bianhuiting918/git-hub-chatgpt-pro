#!/usr/bin/env python3
"""Audit technical completion and scientific endpoint of a NylC Step1 bracket."""
import argparse
import json
import math
import pathlib
import re

import numpy as np
import parmed as pmd

THR267_OG1 = 8961
THR267_HG1 = 8962
L2_C12 = 10287
L2_O2 = 10288
L2_N3 = 10289
BAD_PATTERNS = ("SANDER BOMB", "NaN", "FATAL", "forrtl", "segmentation")
WINDOWS = ("w00", "w01", "w02", "w03")


def cell(box):
    a, b, c, alpha, beta, gamma = [float(v) for v in box]
    alpha, beta, gamma = map(math.radians, (alpha, beta, gamma))
    va = np.array([a, 0.0, 0.0])
    vb = np.array([b * math.cos(gamma), b * math.sin(gamma), 0.0])
    cx = c * math.cos(beta)
    cy = c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma)
    cz = math.sqrt(max(0.0, c * c - cx * cx - cy * cy))
    return np.array([va, vb, [cx, cy, cz]])


def distance(structure, left, right):
    lattice = cell(structure.box)
    delta = structure.coordinates[left - 1] - structure.coordinates[right - 1]
    frac = delta @ np.linalg.inv(lattice)
    frac -= np.rint(frac)
    return float(np.linalg.norm(frac @ lattice))


def last_energy(text):
    values = []
    for line in text.splitlines():
        match = re.match(r"^\s*\d+\s+(-?\d+(?:\.\d*)?(?:[Ee][+-]?\d+)?)\s+", line)
        if match:
            values.append(float(match.group(1)))
    return values[-1] if values else None


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--prmtop", type=pathlib.Path, required=True)
    parser.add_argument("--input-dir", type=pathlib.Path, required=True)
    parser.add_argument("--run-dir", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    manifest = json.loads((args.input_dir / "manifest.json").read_text())
    acceptor = int(manifest["acceptor_atom"])
    records = []
    technical_failures = []
    for name in WINDOWS:
        out = args.run_dir / f"{name}.out"
        rst = args.run_dir / f"{name}.rst7"
        if not out.exists() or not rst.exists():
            technical_failures.append(f"{name}:missing_output")
            continue
        text = out.read_text(errors="replace")
        hits = {token: text.count(token) for token in BAD_PATTERNS}
        final_results = "FINAL RESULTS" in text
        run_done = "Run done" in text or "5.  TIMINGS" in text
        scc_warnings = text.lower().count("scc is not converged")
        if not final_results or not run_done or any(hits.values()) or scc_warnings:
            technical_failures.append(f"{name}:amber_completion_or_error_gate")
        structure = pmd.load_file(str(args.prmtop), xyz=str(rst))
        values = {
            "attack_og1_c12_A": distance(structure, THR267_OG1, L2_C12),
            "thr_og1_hg1_A": distance(structure, THR267_OG1, THR267_HG1),
            "hg1_acceptor_A": distance(structure, THR267_HG1, acceptor),
            "carbonyl_c12_o2_A": distance(structure, L2_C12, L2_O2),
            "amide_c12_n3_A": distance(structure, L2_C12, L2_N3),
        }
        records.append({"window": name, "final_energy_kcal_mol": last_energy(text), "final_results": final_results, "run_done": run_done, "bad_pattern_hits": hits, "scc_warnings": scc_warnings, "distances_A": values})

    technical_status = "PASS_TECHNICAL_STEP1_BRACKET" if not technical_failures and len(records) == len(WINDOWS) else "FAIL_TECHNICAL_STEP1_BRACKET"
    endpoint = records[-1]["distances_A"] if records else {}
    reached = bool(endpoint and endpoint["attack_og1_c12_A"] <= 1.90 and endpoint["hg1_acceptor_A"] <= 1.25 and endpoint["thr_og1_hg1_A"] >= 1.25 and endpoint["amide_c12_n3_A"] <= 1.70)
    scientific_status = "PASS_SCIENTIFIC_TETRAHEDRAL_SEED_REACHED" if technical_status.startswith("PASS") and reached else "FAIL_SCIENTIFIC_TETRAHEDRAL_SEED_NOT_REACHED"
    monotonic_attack = len(records) == len(WINDOWS) and all(records[i + 1]["distances_A"]["attack_og1_c12_A"] <= records[i]["distances_A"]["attack_og1_c12_A"] + 0.10 for i in range(len(records) - 1))
    monotonic_proton = len(records) == len(WINDOWS) and all(records[i + 1]["distances_A"]["hg1_acceptor_A"] <= records[i]["distances_A"]["hg1_acceptor_A"] + 0.10 for i in range(len(records) - 1))
    result = {
        "schema_version": 1,
        "candidate_id": manifest["candidate_id"],
        "acceptor_hypothesis": manifest["acceptor_hypothesis"],
        "technical_status": technical_status,
        "scientific_status": scientific_status,
        "technical_failures": technical_failures,
        "monotonic_checks": {"attack_nonincreasing_with_tolerance": monotonic_attack, "proton_approach_nonincreasing_with_tolerance": monotonic_proton},
        "windows": records,
        "endpoint_thresholds_A": {"attack_og1_c12_max": 1.90, "hg1_acceptor_max": 1.25, "thr_og1_hg1_min": 1.25, "amide_c12_n3_max": 1.70},
        "interpretation": "A reached endpoint is only a constrained tetrahedral-side seed; not a TS, committor, PMF, or barrier.",
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))
    return 0 if technical_status.startswith("PASS") else 31


if __name__ == "__main__":
    raise SystemExit(main())
