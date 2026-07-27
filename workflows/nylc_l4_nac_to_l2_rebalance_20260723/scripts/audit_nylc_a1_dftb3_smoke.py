#!/usr/bin/env python3
"""Audit A1 DFTB3 one-step and 20-step numerical entry without scientific overclaim."""
import argparse
import json
import pathlib
import re

HARD_PATTERNS = {
    "sander_bomb": r"SANDER BOMB",
    "segmentation": r"segmentation",
    "forrtl": r"forrtl",
    "nan": r"\bnan\b",
    "fatal": r"FATAL",
}


def inspect(path):
    source = pathlib.Path(path)
    text = source.read_text(encoding="utf-8", errors="replace") if source.is_file() else ""
    return {
        "path": str(source),
        "run_done": bool(re.search(r"Run\s+done", text)),
        "final_results": "FINAL RESULTS" in text,
        "scc_warnings": len(re.findall(r"Convergence could not be achieved", text, re.I)),
        "vlimit_warnings": len(re.findall(r"vlimit\s+exceeded", text, re.I)),
        "hard_error_hits": {
            name: len(re.findall(pattern, text, re.I))
            for name, pattern in HARD_PATTERNS.items()
        },
        "nquant_values": [int(value) for value in re.findall(r"nquant\s*=\s*(\d+)", text, re.I)],
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--one-step", required=True)
    parser.add_argument("--segment", required=True)
    parser.add_argument("--audit", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    input_audit = json.loads(pathlib.Path(args.audit).read_text(encoding="utf-8"))
    if input_audit.get("status") != "READY_A1_DFTB3_NUMERICAL_PREFLIGHT":
        raise SystemExit("input QM contract is not READY")
    expected = {
        "qmcharge": 0,
        "spin": 1,
        "qm_atom_count": 94,
        "electron_count_including_link_h": 314,
        "link_atom_count": 1,
        "qm_theory": "DFTB3",
        "slater_koster_set": "3ob-3-1",
    }
    for key, value in expected.items():
        if input_audit.get(key) != value:
            raise SystemExit(f"input QM contract mismatch for {key}")
    stages = [inspect(args.one_step), inspect(args.segment)]
    clean = all(
        stage["run_done"]
        and stage["final_results"]
        and stage["scc_warnings"] == 0
        and stage["vlimit_warnings"] == 0
        and sum(stage["hard_error_hits"].values()) == 0
        and (not stage["nquant_values"] or set(stage["nquant_values"]) == {94})
        for stage in stages
    )
    result = {
        "schema_version": 1,
        "status": (
            "PASS_A1_DFTB3_NUMERICAL_PREFLIGHT"
            if clean
            else "FAIL_A1_DFTB3_NUMERICAL_PREFLIGHT"
        ),
        "technical_scope": "Amber18 DFTB3/3OB-3-1 numerical entry",
        "scientific_status": "NOT_EVALUATED_TS_PMF_BARRIER",
        "interpretation": (
            "This is a numerical preflight only; not a TS, reaction coordinate, "
            "PMF, barrier, or production QM-region validation."
        ),
        "production_qm_region_gate": (
            "Expand at least through Asp306/Asp308 before production Step1 work."
        ),
        "qm_audit": input_audit,
        "stages": stages,
    }
    pathlib.Path(args.output).write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"status": result["status"]}, sort_keys=True))
    raise SystemExit(0 if clean else 1)


if __name__ == "__main__":
    main()
