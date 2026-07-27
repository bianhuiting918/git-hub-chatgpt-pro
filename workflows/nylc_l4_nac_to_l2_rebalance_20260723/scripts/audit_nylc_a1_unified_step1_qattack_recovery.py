#!/usr/bin/env python3
"""Audit one corrected unified-core Step1 q_attack recovery replica."""
import argparse
import json
import math
import pathlib
import re

import parmed as pmd

THR267_N = 8949
THR267_OG1 = 8960
N_ALPHA_H = (8950, 8951, 8961)
L2_C12 = 10287
L2_O2 = 10288
L2_N3 = 10289
HARD_PATTERNS = {
    "sander_bomb": r"SANDER BOMB",
    "segmentation": r"segmentation",
    "forrtl": r"forrtl",
    "nan": r"\bnan\b",
    "fatal": r"FATAL",
}


def distance(structure, i, j):
    a, b = structure.atoms[i - 1], structure.atoms[j - 1]
    return math.sqrt((a.xx-b.xx)**2 + (a.xy-b.xy)**2 + (a.xz-b.xz)**2)


def angle(structure, i, j, k):
    a, b, c = [structure.atoms[index - 1] for index in (i, j, k)]
    u = (a.xx-b.xx, a.xy-b.xy, a.xz-b.xz)
    v = (c.xx-b.xx, c.xy-b.xy, c.xz-b.xz)
    denominator = math.sqrt(sum(x*x for x in u)) * math.sqrt(sum(x*x for x in v))
    cosine = sum(x*y for x, y in zip(u, v)) / denominator
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def geometry(prmtop, restart):
    structure = pmd.load_file(str(prmtop), xyz=str(restart))
    return {
        "og1_c12_A": distance(structure, THR267_OG1, L2_C12),
        "o2_c12_og1_deg": angle(structure, L2_O2, L2_C12, THR267_OG1),
        "c12_o2_A": distance(structure, L2_C12, L2_O2),
        "c12_n3_A": distance(structure, L2_C12, L2_N3),
        "nalpha_h_A": [distance(structure, THR267_N, h) for h in N_ALPHA_H],
    }


def inspect_window(prmtop, directory, target):
    output = directory / "window.out"
    restart = directory / "window.rst7"
    text = output.read_text(encoding="utf-8", errors="replace") if output.is_file() else ""
    hard = {name: len(re.findall(pattern, text, re.I)) for name, pattern in HARD_PATTERNS.items()}
    complete = "FINAL RESULTS" in text and bool(re.search(r"Run\s+done", text))
    scc = len(re.findall(r"Convergence could not be achieved", text, re.I))
    vlimit = len(re.findall(r"vlimit\s+exceeded", text, re.I))
    bond_overflow = len(re.findall(r"BOND\s*=\s*\*+", text, re.I))
    records = re.findall(
        r"^\s*(\d+)\s+(-?\d+\.\d+E[+-]\d+)\s+(\d+\.\d+E[+-]\d+)\s+(\d+\.\d+E[+-]\d+)",
        text,
        re.M,
    )
    observed = geometry(prmtop, restart) if restart.is_file() and restart.stat().st_size else {}
    values = []
    if observed:
        values = [
            observed["og1_c12_A"],
            observed["o2_c12_og1_deg"],
            observed["c12_o2_A"],
            observed["c12_n3_A"],
            *observed["nalpha_h_A"],
        ]
    plausible = bool(values) and all(math.isfinite(value) for value in values) and (
        1.1 <= observed["og1_c12_A"] <= 5.0
        and 0.9 <= observed["c12_o2_A"] <= 2.5
        and 1.0 <= observed["c12_n3_A"] <= 3.0
        and 0.0 <= observed["o2_c12_og1_deg"] <= 180.0
    )
    technical = (
        complete
        and scc == 0
        and vlimit == 0
        and bond_overflow == 0
        and sum(hard.values()) == 0
        and plausible
    )
    return {
        "target_A": target,
        "path": str(directory),
        "complete": complete,
        "scc_warnings": scc,
        "vlimit_warnings": vlimit,
        "bond_overflow": bond_overflow,
        "hard_error_hits": hard,
        "geometry_plausible": plausible,
        "geometry": observed,
        "last_minimization_record": (
            {
                "nstep": int(records[-1][0]),
                "energy_kcal_mol": float(records[-1][1]),
                "rms": float(records[-1][2]),
                "gmax": float(records[-1][3]),
            }
            if records else None
        ),
        "technical_pass": technical,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan-root", type=pathlib.Path, required=True)
    parser.add_argument("--source-prmtop", type=pathlib.Path, required=True)
    parser.add_argument("--start-restart", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    scan_root = args.scan_root.resolve()
    manifest = json.loads((scan_root / "SCAN_MANIFEST.json").read_text(encoding="utf-8"))
    if manifest.get("status") != "READY_A1_UNIFIED_STEP1_QATTACK_RECOVERY":
        raise SystemExit("recovery manifest is not READY")
    start_geometry = geometry(args.source_prmtop, args.start_restart)
    windows = [
        inspect_window(args.source_prmtop, scan_root / item["name"], item["target_A"])
        for item in manifest["windows"]
    ]
    technical = len(windows) == 6 and all(window["technical_pass"] for window in windows)
    endpoint = windows[-1]["geometry"] if windows else {}
    bracket = technical and bool(endpoint) and (
        endpoint["og1_c12_A"] <= 1.90
        and 95.0 <= endpoint["o2_c12_og1_deg"] <= 115.0
        and endpoint["c12_o2_A"] <= 1.70
        and endpoint["c12_n3_A"] <= 1.80
    )
    endpoint_response = (
        {
            "delta_og1_c12_A": endpoint["og1_c12_A"] - start_geometry["og1_c12_A"],
            "delta_c12_o2_A": endpoint["c12_o2_A"] - start_geometry["c12_o2_A"],
            "delta_c12_n3_A": endpoint["c12_n3_A"] - start_geometry["c12_n3_A"],
            "delta_nalpha_h_A": [
                value - start
                for value, start in zip(endpoint["nalpha_h_A"], start_geometry["nalpha_h_A"])
            ],
        }
        if endpoint else {}
    )
    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_RECOVERY"
            if technical else "FAIL_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_RECOVERY"
        ),
        "scientific_status": "NOT_EVALUATED_TS_PMF_BARRIER",
        "attack_bracket_seed_gate": (
            "PASS_CONSTRAINED_ATTACK_BRACKET_SEED"
            if bracket else "NOT_EVALUATED_ATTACK_BRACKET_NOT_REACHED"
        ),
        "force_constant_kcal_mol_A2": manifest["protocol"]["distance_force_kcal_mol_A2"],
        "start_geometry": start_geometry,
        "endpoint_response": endpoint_response,
        "windows": windows,
        "interpretation": (
            "Corrected constrained attack scouting. A passing seed is not a TS, "
            "committor, PMF, activation barrier or mechanism proof."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "seed_gate": result["attack_bracket_seed_gate"]}))
    raise SystemExit(0 if technical else 1)


if __name__ == "__main__":
    main()
