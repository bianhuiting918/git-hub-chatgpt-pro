#!/usr/bin/env python3
"""Audit two short Sander DFTB3 stages without overstating scientific meaning."""
import argparse
import json
import pathlib
import re


def inspect(path):
    text = pathlib.Path(path).read_text(errors="replace") if pathlib.Path(path).exists() else ""
    hard = ("SANDER BOMB", "segmentation", "forrtl", "NaN", "FATAL")
    return {
        "path": str(path),
        "run_done": "Run done" in text,
        "final_results": "FINAL RESULTS" in text,
        "scc_warnings": len(re.findall(r"Convergence could not be achieved", text, flags=re.I)),
        "vlimit_warnings": len(re.findall(r"vlimit\s+exceeded", text, flags=re.I)),
        "hard_error_hits": {token: len(re.findall(token, text, flags=re.I)) for token in hard},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--one-step", required=True)
    parser.add_argument("--segment", required=True)
    parser.add_argument("--audit", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    input_audit = json.loads(pathlib.Path(args.audit).read_text())
    stages = [inspect(args.one_step), inspect(args.segment)]
    clean = all(
        stage["run_done"]
        and stage["final_results"]
        and stage["scc_warnings"] == 0
        and stage["vlimit_warnings"] == 0
        and sum(stage["hard_error_hits"].values()) == 0
        for stage in stages
    )
    result = {
        "status": "PASS_DFTB3_PREFLIGHT" if clean else "FAIL_DFTB3_PREFLIGHT",
        "interpretation": "numerical preflight only; not a TS search, reaction coordinate, PMF, barrier, or Step2 endpoint",
        "qm_audit": input_audit,
        "stages": stages,
    }
    pathlib.Path(args.output).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps({"status": result["status"]}, sort_keys=True))
    raise SystemExit(0 if clean else 1)


if __name__ == "__main__":
    main()
