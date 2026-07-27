#!/usr/bin/env python3
"""Scan the corrected nine-replica A1 ensemble for NAC-conditioned direct PT2 geometry."""
import argparse
import csv
import hashlib
import json
import math
import pathlib

import MDAnalysis as mda
import numpy as np
from MDAnalysis.lib.distances import minimize_vectors

EXPECTED_REPLICA_DENOMINATOR = 9
N_ALPHA = 8949
THR267_OG1 = 8960
TRANSFERRED_HG1 = 8961
L2_C12 = 10287
L2_O2 = 10288
L2_N3 = 10289
NAC_DISTANCE_MAX_NM = 0.35
NAC_ANGLE_MIN_DEG = 95.0
NAC_ANGLE_MAX_DEG = 115.0
PT2_DONOR_ACCEPTOR_MAX_NM = 0.35
PT2_H_ACCEPTOR_MAX_NM = 0.25
PT2_ANGLE_MIN_DEG = 135.0
AUDIT_ARRAY_ID = "61976600"


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def pbc_vector(first, second, box):
    return minimize_vectors(np.asarray(second - first, dtype=float)[None, :], box)[0]


def pbc_distance_nm(first, second, box):
    return float(np.linalg.norm(pbc_vector(first, second, box)) / 10.0)


def pbc_angle_deg(first, center, third, box):
    u = pbc_vector(center, first, box)
    v = pbc_vector(center, third, box)
    denominator = float(np.linalg.norm(u) * np.linalg.norm(v))
    cosine = float(np.dot(u, v) / denominator)
    return float(math.degrees(math.acos(max(-1.0, min(1.0, cosine)))))


def read_json_line(path):
    lines = [line for line in path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if not lines:
        raise ValueError(f"empty JSON stdout: {path}")
    return json.loads(lines[-1])


def bonded_pairs(universe):
    return {
        tuple(sorted((bond.atoms[0].index + 1, bond.atoms[1].index + 1)))
        for bond in universe.bonds
    }


def identity(universe, index1):
    atom = universe.atoms[index1 - 1]
    return {
        "index1": index1,
        "name": atom.name,
        "resname": atom.resname,
        "resid": int(atom.resid),
    }


def scan_replica(generate_stdout):
    leaf = generate_stdout.parent
    source = read_json_line(generate_stdout)
    replica_audit_path = leaf / "replica_audit.json"
    replica_audit = json.loads(replica_audit_path.read_text(encoding="utf-8"))
    tpr = pathlib.Path(source["source"]["tpr"])
    xtc = pathlib.Path(source["source"]["xtc"])
    for path in (tpr, xtc):
        if not path.is_file() or path.stat().st_size == 0:
            raise ValueError(f"missing replica source: {path}")
    observed_hashes = {"tpr": sha256(tpr), "xtc": sha256(xtc)}
    expected_hashes = {
        "tpr": source["source"]["tpr_sha256"],
        "xtc": source["source"]["xtc_sha256"],
    }
    if observed_hashes != expected_hashes:
        raise ValueError(f"source sha256 mismatch for {leaf}")

    universe = mda.Universe(str(tpr), str(xtc))
    identities = {
        str(i): identity(universe, i)
        for i in (N_ALPHA, THR267_OG1, TRANSFERRED_HG1, L2_C12, L2_O2, L2_N3)
    }
    expected_names = {
        N_ALPHA: "N", THR267_OG1: "OG1", TRANSFERRED_HG1: "HG1",
        L2_C12: "C12", L2_O2: "O2", L2_N3: "N3",
    }
    if any(identities[str(i)]["name"] != name for i, name in expected_names.items()):
        raise ValueError(f"atom identity mismatch for {leaf}: {identities}")
    bonds = bonded_pairs(universe)
    required = {
        "nalpha_hg1": tuple(sorted((N_ALPHA, TRANSFERRED_HG1))) in bonds,
        "og1_hg1_absent": tuple(sorted((THR267_OG1, TRANSFERRED_HG1))) not in bonds,
        "c12_n3": tuple(sorted((L2_C12, L2_N3))) in bonds,
    }
    if not all(required.values()):
        raise ValueError(f"A1 bond graph mismatch for {leaf}: {required}")

    atoms = universe.atoms
    nac_frames = 0
    positive = []
    best_nac = None
    for ts in universe.trajectory:
        box = ts.dimensions
        if box is None or len(box) < 6 or np.any(np.asarray(box[:3]) <= 0):
            raise ValueError(f"invalid periodic box in {xtc} frame {ts.frame}")
        nalpha = atoms[N_ALPHA - 1].position
        og1 = atoms[THR267_OG1 - 1].position
        hg1 = atoms[TRANSFERRED_HG1 - 1].position
        c12 = atoms[L2_C12 - 1].position
        o2 = atoms[L2_O2 - 1].position
        n3 = atoms[L2_N3 - 1].position
        nac_distance = pbc_distance_nm(og1, c12, box)
        nac_angle = pbc_angle_deg(o2, c12, og1, box)
        is_nac = (
            nac_distance <= NAC_DISTANCE_MAX_NM
            and NAC_ANGLE_MIN_DEG <= nac_angle <= NAC_ANGLE_MAX_DEG
        )
        if not is_nac:
            continue
        nac_frames += 1
        donor_acceptor = pbc_distance_nm(nalpha, n3, box)
        h_acceptor = pbc_distance_nm(hg1, n3, box)
        pt2_angle = pbc_angle_deg(nalpha, hg1, n3, box)
        row = {
            "candidate_id": source["candidate_id"],
            "seed": int(source["velocity_seed"]),
            "time_ps": float(ts.time),
            "frame_index": int(ts.frame),
            "nac_distance_nm": nac_distance,
            "nac_angle_deg": nac_angle,
            "nalpha_n3_nm": donor_acceptor,
            "hg1_n3_nm": h_acceptor,
            "nalpha_hg1_n3_deg": pt2_angle,
            "tpr_sha256": observed_hashes["tpr"],
            "xtc_sha256": observed_hashes["xtc"],
        }
        score = (h_acceptor, -pt2_angle, donor_acceptor)
        if best_nac is None or score < best_nac[0]:
            best_nac = (score, row)
        if (
            donor_acceptor <= PT2_DONOR_ACCEPTOR_MAX_NM
            and h_acceptor <= PT2_H_ACCEPTOR_MAX_NM
            and pt2_angle >= PT2_ANGLE_MIN_DEG
        ):
            positive.append(row)

    expected_nac = int(replica_audit["nac"]["nac_frame_count"])
    if nac_frames != expected_nac:
        raise ValueError(
            f"recomputed NAC count mismatch for {leaf}: {nac_frames} != {expected_nac}"
        )
    return {
        "candidate_id": source["candidate_id"],
        "seed": int(source["velocity_seed"]),
        "frame_count": len(universe.trajectory),
        "nac_frame_count": nac_frames,
        "pt2_preorganized_nac_frame_count": len(positive),
        "source": {
            "generate_stdout": str(generate_stdout),
            "replica_audit": str(replica_audit_path),
            "tpr": str(tpr),
            "xtc": str(xtc),
            "tpr_sha256": observed_hashes["tpr"],
            "xtc_sha256": observed_hashes["xtc"],
        },
        "atom_identity": identities,
        "bond_contract": required,
        "best_nac_by_hg1_n3": best_nac[1] if best_nac else None,
        "positive_frames": positive,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit-root", type=pathlib.Path, required=True)
    parser.add_argument("--output", type=pathlib.Path, required=True)
    args = parser.parse_args()
    root = args.audit_root.resolve()
    generate_files = sorted(root.glob(f"attempt_{AUDIT_ARRAY_ID}_*/*/*/generate.stdout"))
    actual_replica_denominator = len(generate_files)
    if actual_replica_denominator != EXPECTED_REPLICA_DENOMINATOR:
        raise SystemExit(
            f"actual_replica_denominator={actual_replica_denominator}, expected=9"
        )

    replicas = [scan_replica(path) for path in generate_files]
    all_positive = [
        row for replica in replicas for row in replica["positive_frames"]
    ]
    positive_keys = {
        (row["candidate_id"], row["seed"]) for row in all_positive
    }
    independent_positive_replica_count = len(positive_keys)
    reproduced = independent_positive_replica_count >= 2
    scientific_gate = (
        "PASS_A1_PT2_PREORGANIZED_NAC_REPRODUCED"
        if reproduced
        else "NOT_SUPPORTED_A1_PT2_PREORGANIZED_NAC_NOT_REPRODUCED"
    )

    args.output.mkdir(parents=True, exist_ok=False)
    tsv_path = args.output / "candidate_frames.tsv"
    fields = [
        "candidate_id", "seed", "time_ps", "frame_index",
        "nac_distance_nm", "nac_angle_deg", "nalpha_n3_nm",
        "hg1_n3_nm", "nalpha_hg1_n3_deg", "tpr_sha256", "xtc_sha256",
    ]
    with tsv_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(all_positive)

    result = {
        "schema_version": 1,
        "status": "PASS_TECHNICAL_A1_PT2_NAC_ENSEMBLE_SCAN",
        "scientific_status": "NOT_EVALUATED_PROTON_TRANSFER_TS_PMF_BARRIER_MECHANISM",
        "pt2_preorganization_gate": scientific_gate,
        "actual_replica_denominator": actual_replica_denominator,
        "independent_positive_replica_count": independent_positive_replica_count,
        "positive_replica_keys": [
            {"candidate_id": candidate, "seed": seed}
            for candidate, seed in sorted(positive_keys)
        ],
        "total_frame_denominator": sum(x["frame_count"] for x in replicas),
        "total_strict_nac_frames": sum(x["nac_frame_count"] for x in replicas),
        "total_pt2_preorganized_nac_frames": len(all_positive),
        "definitions": {
            "nac": {
                "distance_max_nm": NAC_DISTANCE_MAX_NM,
                "angle_min_deg": NAC_ANGLE_MIN_DEG,
                "angle_max_deg": NAC_ANGLE_MAX_DEG,
            },
            "pt2_preorganization": {
                "nalpha_n3_max_nm": PT2_DONOR_ACCEPTOR_MAX_NM,
                "hg1_n3_max_nm": PT2_H_ACCEPTOR_MAX_NM,
                "nalpha_hg1_n3_min_deg": PT2_ANGLE_MIN_DEG,
                "scope": "fixed-topology geometry screen only",
            },
        },
        "replicas": [
            {key: value for key, value in replica.items() if key != "positive_frames"}
            for replica in replicas
        ],
        "candidate_frames_tsv": str(tsv_path),
        "decision": (
            "Eligible to design a bounded direct PT2/acylation scout from independently sampled frames."
            if reproduced
            else "Do not strongly pull direct HG1-to-N3 from this A1 ensemble; test relay geometry or an alternative microstate."
        ),
        "interpretation": (
            "NAC-conditioned fixed-topology geometric preorganization is not proton transfer, "
            "a transition state, PMF, activation barrier or mechanism proof."
        ),
    }
    result_path = args.output / "RESULT.json"
    result_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "pt2_preorganization_gate": scientific_gate,
        "actual_replica_denominator": actual_replica_denominator,
        "independent_positive_replica_count": independent_positive_replica_count,
        "total_strict_nac_frames": result["total_strict_nac_frames"],
        "total_pt2_preorganized_nac_frames": len(all_positive),
    }))


if __name__ == "__main__":
    main()
