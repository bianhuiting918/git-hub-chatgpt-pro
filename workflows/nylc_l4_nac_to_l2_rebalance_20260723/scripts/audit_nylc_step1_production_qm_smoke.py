#!/usr/bin/env python3
"""Audit one-step DFTB3 output for core and network production QM regions."""
import argparse
import json
import pathlib
import re

BAD_PATTERNS = (
    "Convergence could not be achieved",
    "SANDER BOMB",
    "vlimit exceeded",
    "segmentation",
    "forrtl",
    "NaN",
    "FATAL",
)


def inspect(path):
    text = pathlib.Path(path).read_text(errors="replace")
    bad = {pattern: len(re.findall(re.escape(pattern), text, flags=re.IGNORECASE)) for pattern in BAD_PATTERNS}
    return {
        "final_results": "FINAL RESULTS" in text,
        "run_done": "Run done" in text,
        "bad_pattern_counts": bad,
        "passed": "FINAL RESULTS" in text and "Run done" in text and not any(bad.values()),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--region-audit", required=True)
    parser.add_argument("--core", required=True)
    parser.add_argument("--network", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    region_audit = json.loads(pathlib.Path(args.region_audit).read_text())
    results = {"core": inspect(args.core), "network": inspect(args.network)}
    passed = region_audit.get("status") == "READY_FOR_PRODUCTION_QM_REGION_SMOKE" and all(item["passed"] for item in results.values())
    result = {
        "schema_version": 1,
        "candidate_id": "nylc_C18_trueT267_freeGS",
        "status": "PASS_PRODUCTION_QM_REGION_SMOKE" if passed else "FAIL_PRODUCTION_QM_REGION_SMOKE",
        "regions": results,
        "protonation_microstate": region_audit.get("protonation_microstate"),
        "interpretation": "Numerical one-step result only; not a TS, RC, PMF, or barrier.",
        "next": "Build and compare the Asp306-protonated microstate before Step1 reaction-coordinate scans.",
    }
    pathlib.Path(args.output).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))
    return 0 if passed else 1


if __name__ == "__main__":
    raise SystemExit(main())
