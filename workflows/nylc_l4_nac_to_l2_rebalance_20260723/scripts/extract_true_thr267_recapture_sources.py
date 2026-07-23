#!/usr/bin/env python3
"""Extract corrected NylC L4 recapture starting frames from original XTC files.

These are NOT NAC frames. They are the closest fully unrestrained legacy frames
to the true catalytic Thr267 OG1 and are eligible only for a bounded L4
recapture/release experiment. Legacy atom 8896 is Thr262 OG1 and is explicitly
rejected as SUPERSEDED_WRONG_NUCLEOPHILE_IDENTITY.
"""
import argparse
import datetime as dt
import hashlib
import json
import math
import os
import pathlib
import re
import subprocess
import sys

TRUE_THR267_OG1 = 8961
SUPERSEDED_THR262_OG1 = 8896
PROJECT_ROOT = pathlib.Path(
    "/work/home/acshdt1dks/nylon_pa66_scnet_20260708"
)
TASK_ROOT = PROJECT_ROOT / "l4_nac_to_l2_rebalance_20260723"
SOURCE_ROOT = (
    PROJECT_ROOT
    / "apo_gate_l4_three_carbonyl_20260715"
    / "cyclic_gate_nac_20260719"
    / "branches"
    / "nylc_gyaq"
)
GMX = pathlib.Path(
    "/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi"
)
CANDIDATES = {
    "nylc_C18_trueT267_recapture": {
        "branch": "C18",
        "segment": "segment_000122_FREE_MONITOR",
        "local_time_ps": 16.0,
        "cumulative_time_ps": 11956.0,
        "carbonyl_c": 10297,
        "carbonyl_c_atomname": "C18",
        "carbonyl_o": 10298,
        "carbonyl_o_atomname": "O2",
        "expected_distance_nm": 1.451,
        "expected_angle_deg": 108.009,
        "selection_basis": "closest fully unrestrained true-Thr267 frame; angle compatible",
    },
    "nylc_C23_trueT267_recapture": {
        "branch": "C23",
        "segment": "segment_000320_FREE_MONITOR",
        "local_time_ps": 98.0,
        "cumulative_time_ps": 31838.0,
        "carbonyl_c": 10303,
        "carbonyl_c_atomname": "C23",
        "carbonyl_o": 10304,
        "carbonyl_o_atomname": "O3",
        "expected_distance_nm": 1.310,
        "expected_angle_deg": 97.289,
        "selection_basis": "closest angle-compatible fully unrestrained true-Thr267 frame",
    },
}


def read_gro(path):
    lines = pathlib.Path(path).read_text(encoding="utf-8").splitlines()
    if len(lines) < 4:
        raise ValueError("truncated GRO")
    natoms = int(lines[1].strip())
    if len(lines) != natoms + 3:
        raise ValueError(f"GRO count mismatch: declared {natoms}, lines {len(lines)}")
    atoms = {}
    for global_index, line in enumerate(lines[2 : 2 + natoms], start=1):
        atoms[global_index] = {
            "resid": int(line[0:5]),
            "resname": line[5:10].strip(),
            "atomname": line[10:15].strip(),
            "serial": int(line[15:20]),
            "xyz_nm": tuple(float(line[i : i + 8]) for i in (20, 28, 36)),
        }
    box_values = tuple(float(x) for x in lines[-1].split())
    return atoms, box_values[:3], lines[0]


def minimum_image(vector, box):
    if len(box) != 3 or any(x <= 0 for x in box):
        return vector
    return tuple(v - round(v / length) * length for v, length in zip(vector, box))


def subtract(a, b):
    return tuple(x - y for x, y in zip(a, b))


def norm(v):
    return math.sqrt(sum(x * x for x in v))


def angle_deg(v1, v2):
    denominator = norm(v1) * norm(v2)
    if denominator == 0:
        raise ValueError("zero-length vector in reaction angle")
    cosine = max(-1.0, min(1.0, sum(a * b for a, b in zip(v1, v2)) / denominator))
    return math.degrees(math.acos(cosine))


def validate_identity(atoms, global_index, expected):
    atom = atoms[global_index]
    observed = (atom["resid"], atom["resname"], atom["atomname"])
    if observed != expected:
        raise RuntimeError(
            f"identity mismatch at global index {global_index}: {observed} != {expected}"
        )
    return atom


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def append_history(candidate, status, summary):
    timestamp = dt.datetime.now(dt.timezone.utc).isoformat()
    tsv = TASK_ROOT / "run_history.tsv"
    jsonl = TASK_ROOT / "run_history.jsonl"
    tsv_line = (
        f"{timestamp}\textract_true_thr267_source\t{candidate}\t{status}\t{summary}\n"
    )
    with tsv.open("a", encoding="utf-8") as handle:
        handle.write(tsv_line)
    record = {
        "timestamp_utc": timestamp,
        "stage": "extract_true_thr267_source",
        "candidate": candidate,
        "status": status,
        "summary": summary,
    }
    with jsonl.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True) + "\n")


def extract_one(name):
    cfg = CANDIDATES[name]
    source_dir = SOURCE_ROOT / cfg["branch"] / "segments" / cfg["segment"]
    output_dir = TASK_ROOT / "candidates" / name
    output_dir.mkdir(parents=True, exist_ok=True)
    tmp_gro = output_dir / "source.tmp.gro"
    final_gro = output_dir / "source.gro"
    manifest_path = output_dir / "source_manifest.json"

    if final_gro.exists() or manifest_path.exists():
        raise FileExistsError(
            f"refusing to overwrite existing extracted source for {name}"
        )
    tmp_gro.unlink(missing_ok=True)

    command = [
        str(GMX),
        "trjconv",
        "-s",
        str(source_dir / "run.tpr"),
        "-f",
        str(source_dir / "run.xtc"),
        "-dump",
        f"{cfg['local_time_ps']:.3f}",
        "-o",
        str(tmp_gro),
    ]
    completed = subprocess.run(
        command,
        input="0\n",
        text=True,
        capture_output=True,
        check=False,
        env={**os.environ, "GMX_MAXBACKUP": "-1"},
    )
    (output_dir / "trjconv.log").write_text(
        completed.stdout + "\n--- STDERR ---\n" + completed.stderr,
        encoding="utf-8",
    )
    if completed.returncode != 0:
        append_history(name, "TECHNICAL_FAIL", f"trjconv_rc={completed.returncode}")
        raise RuntimeError(f"trjconv failed for {name}: rc={completed.returncode}")

    atoms, box, title = read_gro(tmp_gro)
    true_thr = validate_identity(atoms, TRUE_THR267_OG1, (267, "THR", "OG1"))
    wrong_thr = validate_identity(
        atoms, SUPERSEDED_THR262_OG1, (262, "THR", "OG1")
    )
    carbonyl_c = atoms[cfg["carbonyl_c"]]
    carbonyl_o = atoms[cfg["carbonyl_o"]]
    if (carbonyl_c["resname"], carbonyl_c["atomname"]) != ("UNL", cfg["carbonyl_c_atomname"]):
        raise RuntimeError(f"reaction carbonyl C identity mismatch: {carbonyl_c}")
    if (carbonyl_o["resname"], carbonyl_o["atomname"]) != ("UNL", cfg["carbonyl_o_atomname"]):
        raise RuntimeError(f"reaction carbonyl O identity mismatch: {carbonyl_o}")

    match = re.search(r"\bt=\s*([0-9.]+)", title)
    if not match:
        raise RuntimeError(f"no time in GRO title: {title!r}")
    extracted_time = float(match.group(1))
    if abs(extracted_time - cfg["local_time_ps"]) > 0.01:
        raise RuntimeError(
            f"wrong extracted time {extracted_time} != {cfg['local_time_ps']}"
        )

    c_to_thr = minimum_image(
        subtract(true_thr["xyz_nm"], carbonyl_c["xyz_nm"]), box
    )
    c_to_o = minimum_image(
        subtract(carbonyl_o["xyz_nm"], carbonyl_c["xyz_nm"]), box
    )
    distance = norm(c_to_thr)
    angle = angle_deg(c_to_thr, c_to_o)
    if abs(distance - cfg["expected_distance_nm"]) > 0.005:
        raise RuntimeError(
            f"distance audit mismatch {distance:.6f} vs {cfg['expected_distance_nm']:.6f}"
        )
    if abs(angle - cfg["expected_angle_deg"]) > 0.5:
        raise RuntimeError(
            f"angle audit mismatch {angle:.3f} vs {cfg['expected_angle_deg']:.3f}"
        )

    digest = sha256(tmp_gro)
    manifest = {
        "schema_version": 1,
        "candidate": name,
        "status": "VERIFIED_TRUE_THR267_RECAPTURE_START_NOT_NAC",
        "scientific_gate": "NOT_EVALUATED_RECAPTURE_AND_RELEASE_REQUIRED",
        "source": {
            "branch": cfg["branch"],
            "segment": cfg["segment"],
            "local_time_ps": extracted_time,
            "cumulative_time_ps": cfg["cumulative_time_ps"],
            "xtc": str(source_dir / "run.xtc"),
            "tpr": str(source_dir / "run.tpr"),
            "selection_basis": cfg["selection_basis"],
        },
        "identity": {
            "true_catalytic_atom": {
                "global_index": TRUE_THR267_OG1,
                "residue_number": true_thr["resid"],
                "residue_name": true_thr["resname"],
                "atom_name": true_thr["atomname"],
            },
            "superseded_wrong_atom": {
                "global_index": SUPERSEDED_THR262_OG1,
                "residue_number": wrong_thr["resid"],
                "residue_name": wrong_thr["resname"],
                "atom_name": wrong_thr["atomname"],
                "status": "SUPERSEDED_WRONG_NUCLEOPHILE_IDENTITY",
            },
            "reaction_carbonyl_c_global_index": cfg["carbonyl_c"],
            "reaction_carbonyl_o_global_index": cfg["carbonyl_o"],
        },
        "geometry": {
            "true_thr267_og1_to_carbonyl_c_nm": round(distance, 6),
            "true_thr267_attack_angle_deg": round(angle, 3),
            "nac_distance_cutoff_nm": 0.35,
            "nac_angle_range_deg": [95.0, 115.0],
            "joint_nac": bool(distance <= 0.35 and 95.0 <= angle <= 115.0),
        },
        "gate_definition": {
            "residues": [261, 262, 263, 264, 265, 266],
            "thr267_excluded": True,
        },
        "output": {
            "gro": str(final_gro),
            "sha256": digest,
        },
    }
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    os.replace(tmp_gro, final_gro)
    append_history(
        name,
        "VERIFIED",
        f"sha256={digest};distance_nm={distance:.6f};angle_deg={angle:.3f}",
    )
    return manifest


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", choices=sorted(CANDIDATES))
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()
    if args.all == bool(args.candidate):
        parser.error("choose exactly one of --all or --candidate")
    names = sorted(CANDIDATES) if args.all else [args.candidate]
    failures = []
    for name in names:
        try:
            result = extract_one(name)
            print(json.dumps(result, sort_keys=True))
        except Exception as exc:
            failures.append((name, str(exc)))
            print(f"{name}: {exc}", file=sys.stderr)
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
