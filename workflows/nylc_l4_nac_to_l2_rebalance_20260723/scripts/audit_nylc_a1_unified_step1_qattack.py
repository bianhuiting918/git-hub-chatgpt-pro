#!/usr/bin/env python3
"""Audit unified-core Step1 attack-only constrained minimization windows."""
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
    left, right = structure.atoms[i - 1], structure.atoms[j - 1]
    return math.sqrt(
        (left.xx - right.xx) ** 2
        + (left.xy - right.xy) ** 2
        + (left.xz - right.xz) ** 2
    )


def angle(structure, i, j, k):
    a, b, c = [structure.atoms[index - 1] for index in (i, j, k)]
    u = (a.xx - b.xx, a.xy - b.xy, a.xz - b.xz)
    v = (c.xx - b.xx, c.xy - b.xy, c.xz - b.xz)
    cosine = sum(x * y for x, y in zip(u, v)) / (
        math.sqrt(sum(x * x for x in u)) * math.sqrt(sum(x * x for x in v))
    )
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def inspect_window(prmtop, directory, target):
    output = directory / "window.out"
    restart = directory / "window.rst7"
    text = output.read_text(encoding="utf-8", errors="replace") if output.is_file() else ""
    hard = {
        name: len(re.findall(pattern, text, re.I))
        for name, pattern in HARD_PATTERNS.items()
    }
    complete = "FINAL RESULTS" in text and bool(re.search(r"Run\s+done", text))
    scc = len(re.findall(r"Convergence could not be achieved", text, re.I))
    vlimit = len(re.findall(r"vlimit\s+exceeded", text, re.I))
    bond_overflow = len(re.findall(r"BOND\s*=\s*\*+", text, re.I))
    nsteps = re.findall(
        r"^\s*(\d+)\s+(-?\d+\.\d+E[+-]\d+)\s+(\d+\.\d+E[+-]\d+)\s+(\d+\.\d+E[+-]\d+)",
        text,
        re.M,
    )
    geometry = {}
    if restart.is_file() and restart.stat().st_size:
        structure = pmd.load_file(str(prmtop), xyz=str(restart))
        geometry = {
            "og1_c12_A": distance(structure, THR267_OG1, L2_C12),
            "o2_c12_og1_deg": angle(structure, L2_O2, L2_C12, THR267_OG1),
            "c12_o2_A": distance(structure, L2_C12, L2_O2),
            "c12_n3_A": distance(structure, L2_C12, L2_N3),
            "nalpha_h_A": [distance(structure, THR267_N, h) for h in N_ALPHA_H],
        }
    geometry_plausible = bool(geometry) and (
        1.1 <= geometry["c12_o2_A"] <= 1.8
        and 1.2 <= geometry["c12_n3_A"] <= 2.2
        and 1.3 <= geometry["og1_c12_A"] <= 5.0
        and 80.0 <= geometry["o2_c12_og1_deg"] <= 130.0
        and all(0.7 <= value <= 2.5 for value in geometry["nalpha_h_A"])
    )
    clean = (
        complete
        and scc == 0
        and vlimit == 0
        and bond_overflow == 0
        and sum(hard.values()) == 0
        and geometry_plausible
    )
    return {
        "target_A": target,
        "path": str(directory),
        "complete": complete,
        "scc_warnings": scc,
        "vlimit_warnings": vlimit,
        "bond_overflow": bond_overflow,
        "hard_error_hits": hard,
        "geometry_plausible": geometry_plausible,
        "geometry": geometry,
        "last_minimization_record": (
            {
                "nstep": int(nsteps[-1][0]),
                "energy_kcal_mol": float(nsteps[-1][1]),
                "rms": float(nsteps[-1][2]),
                "gmax": float(nsteps[-1][3]),
            }
            if nsteps
            else None
        ),
        "technical_pass": clean,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--scan-root", type=pathlib.Path, required=True)
    parser.add_argument("--source-prmtop", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    scan_root = args.scan_root.resolve()
    manifest = json.loads((scan_root / "SCAN_MANIFEST.json").read_text(encoding="utf-8"))
    if manifest.get("status") != "READY_A1_UNIFIED_STEP1_QATTACK_SCAN":
        raise SystemExit("scan manifest is not READY")
    windows = [
        inspect_window(args.source_prmtop, scan_root / item["name"], item["target_A"])
        for item in manifest["windows"]
    ]
    technical = len(windows) == 9 and all(window["technical_pass"] for window in windows)
    endpoint = windows[-1]["geometry"] if windows else {}
    bracket = technical and bool(endpoint) and (
        endpoint["og1_c12_A"] <= 1.90
        and 95.0 <= endpoint["o2_c12_og1_deg"] <= 115.0
        and endpoint["c12_o2_A"] <= 1.70
        and endpoint["c12_n3_A"] <= 1.80
    )
    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_SCAN"
            if technical
            else "FAIL_TECHNICAL_A1_UNIFIED_STEP1_QATTACK_SCAN"
        ),
        "scientific_status": "NOT_EVALUATED_TS_PMF_BARRIER",
        "attack_bracket_seed_gate": (
            "PASS_CONSTRAINED_ATTACK_BRACKET_SEED"
            if bracket
            else "NOT_EVALUATED_ATTACK_BRACKET_NOT_REACHED"
        ),
        "windows": windows,
        "interpretation": (
            "Attack-only constrained minimization scouting. A bracket seed is "
            "not a TS, committor, PMF, activation barrier or mechanism proof."
        ),
    }
    args.output.write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": result["status"], "seed_gate": result["attack_bracket_seed_gate"]}))
    raise SystemExit(0 if technical else 1)


if __name__ == "__main__":
    main()
