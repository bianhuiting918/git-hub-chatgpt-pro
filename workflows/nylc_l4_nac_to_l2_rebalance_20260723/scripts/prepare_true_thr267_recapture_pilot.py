#!/usr/bin/env python3
"""Prepare strict grompp inputs for a gentle true-Thr267 L4 recapture pilot."""
import argparse
import datetime as dt
import hashlib
import json
import os
import pathlib
import re
import subprocess
import sys

TRUE_THR267_OG1 = 8961
SUPERSEDED_THR262_OG1 = 8896
PROJECT_ROOT = pathlib.Path("/work/home/acshdt1dks/nylon_pa66_scnet_20260708")
TASK_ROOT = PROJECT_ROOT / "l4_nac_to_l2_rebalance_20260723"
SOURCE_ROOT = (
    PROJECT_ROOT / "apo_gate_l4_three_carbonyl_20260715"
    / "cyclic_gate_nac_20260719" / "branches" / "nylc_gyaq"
)
GMX = pathlib.Path("/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi")
CANDIDATES = {
    "nylc_C18_trueT267_recapture": {
        "branch": "C18", "carbonyl_c": 10297, "carbonyl_o": 10298,
        "carbonyl_group": "L4_C18", "oxygen_group": "L4_O_C18", "seed": 181267,
    },
    "nylc_C23_trueT267_recapture": {
        "branch": "C23", "carbonyl_c": 10303, "carbonyl_o": 10304,
        "carbonyl_group": "L4_C23", "oxygen_group": "L4_O_C23", "seed": 231267,
    },
}
PROTOCOLS = {
    "pilot1": {"dir": "recapture_pilot", "stem": "pilot", "nsteps": 50000,
               "rate": -0.004, "distance_k": 100, "angle_k": 20},
    "response2": {"dir": "recapture_response2", "stem": "response2", "nsteps": 50000,
                  "rate": -0.002, "distance_k": 500, "angle_k": 50},
    "response3": {"dir": "recapture_response3", "stem": "response3", "nsteps": 12500,
                  "rate": -0.004, "distance_k": 2000, "angle_k": 100},
    "response4": {"dir": "recapture_response4", "stem": "response4", "nsteps": 50000,
                  "rate": -0.004, "distance_k": 2000, "angle_k": 100,
                  "parent_dir": "recapture_response3", "parent_stem": "response3"},
    "response5": {"dir": "recapture_response5", "stem": "response5", "nsteps": 50000,
                  "rate": -0.004, "distance_k": 2000, "angle_k": 500,
                  "parent_dir": "recapture_response4", "parent_stem": "response4"},
    "response6": {"dir": "recapture_response6", "stem": "response6", "nsteps": 50000,
                  "rate": -0.004, "distance_k": 2000, "angle_k": 500,
                  "parent_dir": "recapture_response5", "parent_stem": "response5"},
    "response7": {"dir": "recapture_response7", "stem": "response7", "nsteps": 50000,
                  "rate": -0.004, "distance_k": 2000, "angle_k": 500,
                  "parent_dir": "recapture_response6", "parent_stem": "response6"},
    "response8": {"dir": "recapture_response8", "stem": "response8", "nsteps": 50000,
                  "rate": -0.0015, "distance_k": 2000, "angle_k": 500,
                  "parent_dir": "recapture_response7", "parent_stem": "response7"},
    "response9": {"dir": "recapture_response9", "stem": "response9", "nsteps": 50000,
                  "rate": -0.002, "distance_k": 2000, "angle_k": 500,
                  "parent_dir": "recapture_response8", "parent_stem": "response8"},
    "response10": {"dir": "recapture_response10", "stem": "response10", "nsteps": 50000,
                   "rate": 0.0, "distance_k": 5000, "angle_k": 500,
                   "target_distance": 0.35,
                   "parent_dir": "recapture_response9", "parent_stem": "response9"},
}


def sha256(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_index_group(text, group):
    pattern = re.compile(
        rf"(?ms)^\[\s*{re.escape(group)}\s*\]\s*\n(.*?)(?=^\s*\[|\Z)"
    )
    match = pattern.search(text)
    if not match:
        raise ValueError(f"missing index group {group}")
    return [int(x) for x in re.findall(r"\d+", match.group(1))]


def correct_index(text):
    observed = parse_index_group(text, "Thr_OG1")
    if observed != [SUPERSEDED_THR262_OG1]:
        raise RuntimeError(
            f"SUPERSEDED_WRONG_NUCLEOPHILE_IDENTITY guard: "
            f"expected legacy Thr_OG1 [{SUPERSEDED_THR262_OG1}], got {observed}"
        )
    pattern = re.compile(
        r"(?ms)(^\[\s*Thr_OG1\s*\]\s*\n)(.*?)(?=^\s*\[|\Z)"
    )
    return pattern.sub(rf"\g<1>{TRUE_THR267_OG1}\n\n", text, count=1)


def read_gro_identities(path):
    lines = pathlib.Path(path).read_text(encoding="utf-8").splitlines()
    natoms = int(lines[1].strip())
    if len(lines) != natoms + 3:
        raise ValueError("GRO atom-count mismatch")
    atoms = {}
    for index, line in enumerate(lines[2:2 + natoms], start=1):
        atoms[index] = {
            "resid": int(line[:5]),
            "resname": line[5:10].strip(),
            "atomname": line[10:15].strip(),
        }
    return atoms


def make_mdp(branch, seed, protocol="pilot1"):
    p = PROTOCOLS[protocol]
    carbon = f"L4_C{branch[1:]}"
    oxygen = f"L4_O_C{branch[1:]}"
    is_continuation = "parent_dir" in p
    continuation = "yes" if is_continuation else "no"
    velocity_block = "gen-vel                  = no" if is_continuation else (
        f"gen-vel                  = yes\ngen-temp                 = 300\ngen-seed                 = {seed}"
    )
    fixed_target = p.get("target_distance")
    distance_start = "no" if fixed_target is not None else "yes"
    distance_init = f"{fixed_target:.6f}" if fixed_target is not None else "0"
    return f"""integrator               = md
dt                       = 0.002
nsteps                   = {p["nsteps"]}
continuation             = {continuation}
constraints              = h-bonds
constraint-algorithm     = lincs
lincs-iter               = 1
lincs-order              = 4
cutoff-scheme            = Verlet
nstlist                  = 20
rlist                    = 1.0
coulombtype              = PME
rcoulomb                 = 1.0
vdwtype                  = Cut-off
rvdw                     = 1.0
DispCorr                 = EnerPres
tcoupl                   = V-rescale
tc-grps                  = System
tau-t                    = 1.0
ref-t                    = 300
pcoupl                   = no
pbc                      = xyz
nstxout-compressed       = 500
nstenergy                = 500
nstlog                   = 500
{velocity_block}
comm-mode                = Linear
nstcomm                  = 100
pull                     = yes
pull-ngroups             = 3
pull-ncoords             = 2
pull-group1-name         = Thr_OG1
pull-group2-name         = {carbon}
pull-group3-name         = {oxygen}

pull-coord1-type         = umbrella
pull-coord1-geometry     = distance
pull-coord1-groups       = 1 2
pull-coord1-dim          = Y Y Y
pull-coord1-start        = {distance_start}
pull-coord1-init         = {distance_init}
pull-coord1-rate         = {p["rate"]:.6f}
pull-coord1-k            = {p["distance_k"]}

pull-coord2-type         = umbrella
pull-coord2-geometry     = angle
pull-coord2-groups       = 2 1 2 3
pull-coord2-start        = no
pull-coord2-init         = 105.000000
pull-coord2-rate         = 0
pull-coord2-k            = {p["angle_k"]}

pull-nstxout             = 500
pull-nstfout             = 500
"""


def append_history(candidate, state, detail):
    stamp = dt.datetime.now(dt.timezone.utc).isoformat()
    with (TASK_ROOT / "run_history.tsv").open("a", encoding="utf-8") as handle:
        handle.write(f"{stamp}\ttrue_thr267_recapture_preflight\t{state}\t"
                     f"candidate={candidate};{detail}\n")
    record = {
        "timestamp_utc": stamp,
        "event": "true_thr267_recapture_preflight",
        "state": state,
        "candidate": candidate,
        "detail": detail,
    }
    with (TASK_ROOT / "run_history.jsonl").open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def prepare_one(name, protocol="pilot1"):
    cfg = CANDIDATES[name]
    p = PROTOCOLS[protocol]
    seed = cfg["seed"] + {"pilot1": 0, "response2": 1, "response3": 2, "response4": 3, "response5": 4, "response6": 5, "response7": 6, "response8": 7, "response9": 8, "response10": 9}[protocol]
    candidate_root = TASK_ROOT / "candidates" / name
    source_gro = candidate_root / "source.gro"
    source_manifest = candidate_root / "source_manifest.json"
    if not source_gro.exists() or not source_manifest.exists():
        raise FileNotFoundError(f"verified source missing for {name}")
    manifest = json.loads(source_manifest.read_text(encoding="utf-8"))
    if manifest["status"] != "VERIFIED_TRUE_THR267_RECAPTURE_START_NOT_NAC":
        raise RuntimeError(f"source status is not verified for {name}")

    checkpoint = None
    input_coord = source_gro
    if "parent_dir" in p:
        parent = candidate_root / p["parent_dir"]
        input_coord = parent / f"{p['parent_stem']}.gro"
        checkpoint = parent / f"{p['parent_stem']}.cpt"
        parent_audit = json.loads(
            (parent / f"{p['parent_stem']}_audit.json").read_text()
        )
        if not parent_audit["finished_mdrun"] or any(
            parent_audit["numerical_issue_counts"].values()
        ):
            raise RuntimeError(f"parent stage is not numerically valid for {name}")
        if parent_audit["geometry"]["distance_response_nm"] <= 0:
            raise RuntimeError(f"parent stage did not move toward Thr267 for {name}")
        if not input_coord.exists() or not checkpoint.exists():
            raise FileNotFoundError(f"parent coordinate/checkpoint missing for {name}")

    source_branch = SOURCE_ROOT / cfg["branch"]
    original_index = source_branch / "cycle.ndx"
    topology = source_branch / "topol.top"
    work = candidate_root / p["dir"]
    stem = p["stem"]
    work.mkdir(parents=True, exist_ok=True)
    outputs = [
        work / "corrected_true_thr267.ndx", work / f"{stem}.mdp",
        work / f"{stem}.tpr", work / f"{stem}.processed.mdp",
        work / "grompp.stdout", work / "grompp.stderr",
        work / "preflight.json",
    ]
    if any(path.exists() for path in outputs):
        raise FileExistsError(f"refusing to overwrite recapture preflight for {name}")

    atoms = read_gro_identities(source_gro)
    true_identity = atoms[TRUE_THR267_OG1]
    wrong_identity = atoms[SUPERSEDED_THR262_OG1]
    if true_identity != {"resid": 267, "resname": "THR", "atomname": "OG1"}:
        raise RuntimeError(f"true Thr267 identity mismatch: {true_identity}")
    if wrong_identity != {"resid": 262, "resname": "THR", "atomname": "OG1"}:
        raise RuntimeError(f"legacy Thr262 identity mismatch: {wrong_identity}")

    original_text = original_index.read_text(encoding="utf-8")
    gate_atoms = parse_index_group(original_text, "Gate")
    gate_residues = sorted({atoms[index]["resid"] for index in gate_atoms})
    if gate_residues != [261, 262, 263, 264, 265, 266]:
        raise RuntimeError(f"gate identity mismatch: residues {gate_residues}")
    if TRUE_THR267_OG1 in gate_atoms:
        raise RuntimeError("Thr267 must not be in Gate group")
    corrected_text = correct_index(original_text)
    if parse_index_group(corrected_text, "Thr_OG1") != [TRUE_THR267_OG1]:
        raise RuntimeError("corrected Thr_OG1 group failed")
    (work / "corrected_true_thr267.ndx").write_text(corrected_text, encoding="utf-8")
    (work / f"{stem}.mdp").write_text(
        make_mdp(cfg["branch"], seed, protocol=protocol), encoding="utf-8"
    )

    command = [
        str(GMX), "grompp", "-f", str(work / f"{stem}.mdp"),
        "-c", str(input_coord), "-p", str(topology),
        "-n", str(work / "corrected_true_thr267.ndx"),
        "-o", str(work / f"{stem}.tpr"),
        "-po", str(work / f"{stem}.processed.mdp"), "-maxwarn", "0",
    ]
    if checkpoint is not None:
        command.extend(["-t", str(checkpoint)])
    completed = subprocess.run(
        command, cwd=source_branch, capture_output=True, text=True,
        env={**os.environ, "GMX_MAXBACKUP": "-1"}, check=False,
    )
    (work / "grompp.stdout").write_text(completed.stdout, encoding="utf-8")
    (work / "grompp.stderr").write_text(completed.stderr, encoding="utf-8")
    status = "PASS_GROMPP_TRUE_THR267" if completed.returncode == 0 else "FAIL_GROMPP"
    audit = {
        "schema_version": 1,
        "candidate": name,
        "status": status,
        "grompp_returncode": completed.returncode,
        "maxwarn": 0,
        "source_gro_sha256": sha256(source_gro),
        "input_coordinate": str(input_coord),
        "input_coordinate_sha256": sha256(input_coord),
        "parent_checkpoint": str(checkpoint) if checkpoint else None,
        "topology_sha256": sha256(topology),
        "corrected_index_sha256": sha256(work / "corrected_true_thr267.ndx"),
        "true_thr267_og1_global_index": TRUE_THR267_OG1,
        "superseded_thr262_og1_global_index": SUPERSEDED_THR262_OG1,
        "gate_residues": gate_residues,
        "thr267_excluded_from_gate": True,
        "protocol": {
            "duration_ps": p["nsteps"] * 0.002,
            "protocol": protocol,
            "distance_reference_rate_nm_per_ps": p["rate"],
            "distance_k_kj_mol_nm2": p["distance_k"],
            "angle_reference_deg": 105.0,
            "angle_k_kj_mol_rad2": p["angle_k"],
            "gate_restraint": False,
            "fresh_velocity_seed": None if checkpoint else seed,
        },
        "scientific_gate": "NOT_EVALUATED_APPROACH_PILOT_NOT_AN_UNRESTRAINED_NAC",
    }
    (work / "preflight.json").write_text(
        json.dumps(audit, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    append_history(name, status, f"grompp_rc={completed.returncode};maxwarn=0")
    if completed.returncode != 0:
        raise RuntimeError(f"grompp failed for {name}")
    return audit


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", choices=sorted(CANDIDATES))
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--protocol", choices=sorted(PROTOCOLS), default="pilot1")
    args = parser.parse_args()
    if args.all == bool(args.candidate):
        parser.error("choose exactly one of --all or --candidate")
    names = sorted(CANDIDATES) if args.all else [args.candidate]
    failures = []
    for name in names:
        try:
            print(json.dumps(prepare_one(name, protocol=args.protocol), sort_keys=True))
        except Exception as exc:
            append_history(name, "FAIL", str(exc))
            print(f"{name}: {exc}", file=sys.stderr)
            failures.append(name)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
