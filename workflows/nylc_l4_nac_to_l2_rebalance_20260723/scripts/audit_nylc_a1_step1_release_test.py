#!/usr/bin/env python3
"""Independently audit local and full release from the constrained A1 attack seed."""
import argparse
import json
import math
import pathlib
import re

import parmed as pmd

THR267_N = 8949
THR267_OG1 = 8960
N_ALPHA_H = (8950, 8951, 8961)
C12_C11 = 10286
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
    carbonyl_angles = [
        angle(structure, C12_C11, L2_C12, L2_O2),
        angle(structure, L2_O2, L2_C12, L2_N3),
        angle(structure, L2_N3, L2_C12, C12_C11),
    ]
    angle_sum = sum(carbonyl_angles)
    return {
        "og1_c12_A": distance(structure, THR267_OG1, L2_C12),
        "o2_c12_og1_deg": angle(structure, L2_O2, L2_C12, THR267_OG1),
        "c12_o2_A": distance(structure, L2_C12, L2_O2),
        "c12_n3_A": distance(structure, L2_C12, L2_N3),
        "c12_c11_A": distance(structure, L2_C12, C12_C11),
        "nalpha_h_A": [distance(structure, THR267_N, h) for h in N_ALPHA_H],
        "carbonyl_neighbor_angles_deg": carbonyl_angles,
        "carbonyl_neighbor_angle_sum_deg": angle_sum,
        "carbonyl_pyramidalization_deg": 360.0 - angle_sum,
    }


def inspect_stage(prmtop, directory):
    output = directory / "stage.out"
    restart = directory / "stage.rst7"
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
    scalar_values = [
        observed.get("og1_c12_A"),
        observed.get("o2_c12_og1_deg"),
        observed.get("c12_o2_A"),
        observed.get("c12_n3_A"),
        observed.get("c12_c11_A"),
        observed.get("carbonyl_neighbor_angle_sum_deg"),
        observed.get("carbonyl_pyramidalization_deg"),
        *observed.get("nalpha_h_A", []),
    ]
    plausible = len(scalar_values) == 10 and all(
        isinstance(value, (int, float)) and math.isfinite(value) for value in scalar_values
    ) and (
        1.1 <= observed["og1_c12_A"] <= 5.0
        and 0.9 <= observed["c12_o2_A"] <= 2.5
        and 1.0 <= observed["c12_n3_A"] <= 3.0
        and 1.0 <= observed["c12_c11_A"] <= 2.0
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
    parser.add_argument("--release-root", type=pathlib.Path, required=True)
    parser.add_argument("--source-prmtop", type=pathlib.Path, required=True)
    parser.add_argument("--start-restart", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    release_root = args.release_root.resolve()
    manifest = json.loads((release_root / "RELEASE_MANIFEST.json").read_text(encoding="utf-8"))
    if manifest.get("status") != "READY_A1_STEP1_RELEASE_TEST":
        raise SystemExit("release manifest is not READY")
    names = [stage["name"] for stage in manifest.get("stages", [])]
    if names != ["local_release", "full_release"]:
        raise SystemExit("release stages do not match the frozen two-stage protocol")

    source_geometry = geometry(args.source_prmtop, args.start_restart)
    stages = {
        name: inspect_stage(args.source_prmtop, release_root / name)
        for name in names
    }
    technical = all(stage["technical_pass"] for stage in stages.values())
    endpoint = stages["full_release"]["geometry"] if technical else {}
    if technical and endpoint["og1_c12_A"] <= 1.90 and 95.0 <= endpoint["o2_c12_og1_deg"] <= 115.0:
        release_gate = "PASS_RELEASE_RETAINED_ATTACK_CONTACT"
    elif technical and endpoint["og1_c12_A"] >= 2.20:
        release_gate = "FAIL_RELEASE_RETURNED_TOWARD_REACTANT"
    else:
        release_gate = "NOT_EVALUATED_RELEASE_AMBIGUOUS"

    result = {
        "schema_version": 1,
        "status": (
            "PASS_TECHNICAL_A1_STEP1_RELEASE_TEST"
            if technical else "FAIL_TECHNICAL_A1_STEP1_RELEASE_TEST"
        ),
        "release_gate": release_gate,
        "scientific_status": "NOT_EVALUATED_TS_PMF_BARRIER",
        "source_geometry": source_geometry,
        "stages": stages,
        "interpretation": (
            "The release gate classifies persistence of a short attack contact only. "
            "Carbonyl bond lengths and pyramidalization are reported without declaring "
            "a tetrahedral intermediate. This is not a TS, committor, PMF, barrier or mechanism proof."
        ),
    }
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "release_gate": release_gate}))
    raise SystemExit(0 if technical else 1)


if __name__ == "__main__":
    main()
