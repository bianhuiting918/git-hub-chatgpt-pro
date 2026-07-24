#!/usr/bin/env python3
"""Prepare staged OD1/OD2 proton-approach brackets from a passed attack seed."""
import argparse
import hashlib
import json
import pathlib

ATTACK_SEED_SHA256 = "d5e607dd722be4fd6e5d33016afb4975d18b200a78d1fceedf5bfb9fed1b87bf"
THR267_OG1 = 8961
THR267_HG1 = 8962
L2_C12 = 10287
ASP306_OD1 = 9572
ASP306_OD2 = 9573
MM_RESTRAINT_MASK = "!(@8949-8964,9568-9573,9587-9592,10273-10351,81218-81220)&!@H="
SCHEDULES = {
    "OD1": [
        {"name": "p00", "attack_A": 2.55, "proton_A": 2.80, "force": 50.0, "maxcyc": 150, "attack_gate_A": 2.65, "hacc_gate_A": 3.15},
        {"name": "p01", "attack_A": 2.30, "proton_A": 2.20, "force": 60.0, "maxcyc": 200, "attack_gate_A": 2.55, "hacc_gate_A": 2.55},
        {"name": "p02", "attack_A": 2.05, "proton_A": 1.60, "force": 70.0, "maxcyc": 200, "attack_gate_A": 2.30, "hacc_gate_A": 1.95},
        {"name": "p03", "attack_A": 1.75, "proton_A": 1.10, "force": 85.0, "maxcyc": 250, "attack_gate_A": 1.95, "hacc_gate_A": 1.30},
    ],
    "OD2": [
        {"name": "p00", "attack_A": 2.55, "proton_A": 4.55, "force": 45.0, "maxcyc": 150, "attack_gate_A": 2.65, "hacc_gate_A": 4.90},
        {"name": "p01", "attack_A": 2.30, "proton_A": 3.60, "force": 55.0, "maxcyc": 200, "attack_gate_A": 2.55, "hacc_gate_A": 3.95},
        {"name": "p02", "attack_A": 2.05, "proton_A": 2.20, "force": 70.0, "maxcyc": 200, "attack_gate_A": 2.30, "hacc_gate_A": 2.55},
        {"name": "p03", "attack_A": 1.75, "proton_A": 1.10, "force": 85.0, "maxcyc": 250, "attack_gate_A": 1.95, "hacc_gate_A": 1.30},
    ],
}


def sha256(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def restraint(atom1, atom2, target, force):
    return (
        f"&rst iat={atom1},{atom2}, r1={max(0.0,target-0.20):.3f}, "
        f"r2={max(0.0,target-0.03):.3f}, r3={target+0.03:.3f}, "
        f"r4={target+0.20:.3f}, rk2={force:.2f}, rk3={force:.2f}, /"
    )


def mdin(stage, qmmask, maxcyc):
    return f"""NylC Step1 {stage} staged proton bracket
&cntrl
  imin=1, ntmin=2, maxcyc={maxcyc}, ncyc={maxcyc}, dx0=0.005,
  ntb=1, cut=10.0, ntpr=5,
  ifqnt=1, nmropt=1, ntr=1,
  restraint_wt=1.0,
  restraintmask='{MM_RESTRAINT_MASK}',
/
&qmmm
  qmmask='{qmmask}',
  qmcharge=-1, spin=1,
  qm_theory='DFTB3',
  dftb_telec=200.0,
  qmshake=0,
/
&wt type='END' /
DISANG=../input/{stage}.RST
LISTOUT=POUT
"""


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-manifest", type=pathlib.Path, required=True)
    parser.add_argument("--attack-pass", type=pathlib.Path, required=True)
    parser.add_argument("--attack-seed", type=pathlib.Path, required=True)
    parser.add_argument("--acceptor", choices=("OD1", "OD2"), required=True)
    parser.add_argument("--output-dir", type=pathlib.Path, required=True)
    args = parser.parse_args()
    source_manifest = json.loads(args.source_manifest.read_text())
    attack_pass = json.loads(args.attack_pass.read_text())
    if attack_pass.get("status") != "PASS_ATTACK_PREREQUISITE":
        raise ValueError("attack seed did not pass prerequisite gate")
    if sha256(args.attack_seed) != ATTACK_SEED_SHA256:
        raise ValueError("attack seed SHA mismatch")
    qm = source_manifest.get("qm_region", {})
    if qm.get("atom_count") != 110 or qm.get("qmcharge") != -1 or qm.get("electron_count") != 388:
        raise ValueError("production core QM definition changed")
    qmmask = qm.get("qmmask")
    if not qmmask or not qmmask.startswith("@"):
        raise ValueError("missing audited qmmask")
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    acceptor_atom = ASP306_OD1 if args.acceptor == "OD1" else ASP306_OD2
    schedule = SCHEDULES[args.acceptor]
    for spec in schedule:
        name = spec["name"]
        (output / f"{name}.in").write_text(mdin(name, qmmask, spec["maxcyc"]))
        lines = [
            restraint(THR267_OG1, L2_C12, spec["attack_A"], spec["force"]),
            restraint(THR267_HG1, acceptor_atom, spec["proton_A"], spec["force"]),
        ]
        (output / f"{name}.RST").write_text("\n".join(lines) + "\n")
    manifest = {
        "schema_version": 1,
        "status": "READY_FOR_STEP1_PROTON_BRACKET",
        "candidate_id": "nylc_C18_trueT267_freeGS",
        "attack_seed_sha256": ATTACK_SEED_SHA256,
        "acceptor_hypothesis": args.acceptor,
        "acceptor_atom": acceptor_atom,
        "reactive_atoms": {"thr267_og1": THR267_OG1, "thr267_hg1": THR267_HG1, "l2_c12": L2_C12},
        "qm_region": qm,
        "schedule": schedule,
        "restraint_scope": "attack and HG1-to-acceptor distances; OG1-HG1 is observed but never restrained",
        "interpretation": "Constrained bracket only; not a TS, committor, PMF, or barrier.",
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, sort_keys=True))


if __name__ == "__main__":
    raise SystemExit(main())
