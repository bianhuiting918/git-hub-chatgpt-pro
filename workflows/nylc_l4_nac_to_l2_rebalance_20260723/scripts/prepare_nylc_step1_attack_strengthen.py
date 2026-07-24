#!/usr/bin/env python3
"""Prepare one stronger attack-only Step1 stage after a safe weak-force plateau."""
import argparse
import json
import pathlib

from prepare_nylc_step1_proton_bracket import L2_C12, THR267_OG1, mdin, restraint

SPEC = {"name": "p00", "attack_A": 2.55, "force": 100.0, "maxcyc": 300, "attack_gate_A": 2.65, "hacc_gate_A": 99.0}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-manifest", type=pathlib.Path, required=True)
    parser.add_argument("--source-fail", type=pathlib.Path, required=True)
    parser.add_argument("--acceptor", choices=("OD1", "OD2"), required=True)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()
    source_manifest = json.loads(args.source_manifest.read_text())
    source_fail = json.loads(args.source_fail.read_text())
    if source_fail.get("status") != "FAIL_SCIENTIFIC_PROTON_BRACKET_STAGE":
        raise ValueError("source is not a retained scientific bracket failure")
    if source_fail.get("stage") != "p00":
        raise ValueError("source failure is not the p00 attack stage")
    if source_fail.get("acceptor_hypothesis") != args.acceptor:
        raise ValueError("source acceptor does not match requested branch")
    gates = source_fail.get("gates", {})
    if not gates.get("chemistry_valid") or not gates.get("proton_attached"):
        raise ValueError("source chemistry is unsafe for stronger attack")
    qm = source_manifest.get("qm_region", {})
    if qm.get("atom_count") != 110 or qm.get("qmcharge") != -1 or qm.get("electron_count") != 388:
        raise ValueError("production QM core changed")
    qmmask = qm.get("qmmask")
    if not qmmask:
        raise ValueError("audited qmmask is missing")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / "p00.in").write_text(mdin("p00", qmmask, SPEC["maxcyc"]))
    (output / "p00.RST").write_text(
        restraint(THR267_OG1, L2_C12, SPEC["attack_A"], SPEC["force"]) + "\n"
    )
    manifest = {
        "schema_version": 1,
        "status": "READY_FOR_STRENGTHENED_ATTACK_ONLY",
        "candidate_id": "nylc_C18_trueT267_freeGS",
        "acceptor_hypothesis": args.acceptor,
        "acceptor_atom": source_manifest["acceptor_atom"],
        "reactive_atoms": source_manifest["reactive_atoms"],
        "qm_region": qm,
        "schedule": [SPEC],
        "restraint_scope": "Only OG1-C12 attack distance is strengthened; OG1-HG1 remains observed but is never restrained.",
        "interpretation": "Single-variable constrained seed test; not a TS, committor, PMF, or barrier.",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, sort_keys=True))


if __name__ == "__main__":
    raise SystemExit(main())
