#!/usr/bin/env python3
"""Audit two frozen A1 NAC+PT2 GRO extractions using the existing frame auditor."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import os
import pathlib
import sys
from typing import Any

import numpy as np

import audit_nylc_a1_representative_frame as base

TASK_ROOT = pathlib.Path(
    "/work/home/acshdt1dks/nylon_pa66_scnet_20260708/"
    "l4_nac_to_l2_rebalance_20260723"
)
SOURCE_BASE = TASK_ROOT / "a1_activated_nac_20260726/equilibration/nac_evt25_time1462ps"
SOURCE_NDX = (
    TASK_ROOT
    / "ensemble/candidates/nac_evt25_time1462ps/"
    "build_job_61813799_11/build/source_cycle.ndx"
)
SOURCES = {
    "seed26723_t378_f189": {
        "seed": "seed26723", "time_ps": 378.0, "frame_index": 189,
        "tpr": SOURCE_BASE / "seed26723/attempt_61970146_4_61970151/npt300free/run.tpr",
        "xtc": SOURCE_BASE / "seed26723/attempt_61970146_4_61970151/npt300free/run.xtc",
        "tpr_sha256": "c60078a92c2ace51facde4ef64e453f690177fc88b4d6363427f935944fa2e43",
        "xtc_sha256": "1a54f1b5b9f139b746c22d9e0f7e9a4a94eb8154bf2b881888986eedca933d89",
        "nac_distance_nm": 0.3031186720, "nac_angle_deg": 113.41893695,
        "nalpha_n3_nm": 0.2939767713, "hg1_n3_nm": 0.1930546162,
        "pt2_angle_deg": 166.32461719,
    },
    "seed26737_t676_f338": {
        "seed": "seed26737", "time_ps": 676.0, "frame_index": 338,
        "tpr": SOURCE_BASE / "seed26737/attempt_61970146_5_61970146/npt300free/run.tpr",
        "xtc": SOURCE_BASE / "seed26737/attempt_61970146_5_61970146/npt300free/run.xtc",
        "tpr_sha256": "dbd19a399547319d10630430ed494d33f6271bab0f30cb0de5a466c6af13ba20",
        "xtc_sha256": "fcba14da98b331368061dcd990f2467628ad77b9b7a9c4ce88090f92e0831b05",
        "nac_distance_nm": 0.3186696796, "nac_angle_deg": 108.09642650,
        "nalpha_n3_nm": 0.2960844435, "hg1_n3_nm": 0.1982272745,
        "pt2_angle_deg": 158.22789159,
    },
}

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
TIME_TOLERANCE_PS = 0.001
DISTANCE_TOLERANCE_NM = 0.002
ANGLE_TOLERANCE_DEG = 0.5
PASS_STATUS = "PASS_A1_PT2_FRAME_EXTRACTION_AUDIT"
FAIL_STATUS = "NOT_EVALUATED_A1_PT2_FRAME_EXTRACTION"
SCOPE = "preorganization_only_not_proton_transfer_ts_pmf_barrier_or_mechanism"


def sha256(path: pathlib.Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: pathlib.Path, payload: dict[str, Any]) -> None:
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def source_for(candidate: str) -> dict[str, Any]:
    if candidate not in SOURCES:
        raise ValueError(f"unknown candidate {candidate}")
    return SOURCES[candidate]


def prepare(candidate: str, candidate_dir: pathlib.Path) -> None:
    source_for(candidate)
    candidate_dir.mkdir(parents=True, exist_ok=False)


def atom_identity(atom: Any, index1: int, resname: str, name: str) -> dict[str, Any]:
    record = {
        "index1": int(atom.index) + 1,
        "resid": int(atom.resid),
        "resname": str(atom.resname),
        "name": str(atom.name),
    }
    if (record["index1"], record["resname"], record["name"]) != (
        index1, resname, name
    ):
        raise ValueError(f"atom identity mismatch: {record}")
    return record


def vector(a: Any, b: Any, box: np.ndarray) -> np.ndarray:
    return base._minimum_image_A(
        np.asarray(b.position, dtype=float) - np.asarray(a.position, dtype=float),
        box,
    )


def distance_nm(a: Any, b: Any, box: np.ndarray) -> float:
    return float(np.linalg.norm(vector(a, b, box)) / 10.0)


def angle_deg(a: Any, b: Any, c: Any, box: np.ndarray) -> float:
    left, right = vector(b, a, box), vector(b, c, box)
    value = np.dot(left, right) / (np.linalg.norm(left) * np.linalg.norm(right))
    return float(np.degrees(np.arccos(np.clip(value, -1.0, 1.0))))


def fail(candidate: str, candidate_dir: pathlib.Path, stage: str, detail: str) -> None:
    payload = {
        "schema_version": 1,
        "status": FAIL_STATUS,
        "scientific_status": "NOT_EVALUATED",
        "scientific_scope": SCOPE,
        "candidate": candidate,
        "failed_stage": stage,
        "detail": detail,
    }
    write_json(candidate_dir / "manifest.json", payload)
    write_json(candidate_dir / "NOT_EVALUATED.json", payload)


def audit(candidate: str, candidate_dir: pathlib.Path, commit: str) -> None:
    source = source_for(candidate)
    source_tmp = candidate_dir / "source.tmp.gro"
    source_gro = candidate_dir / "source.gro"
    for path in (
        source_gro, candidate_dir / "manifest.json", candidate_dir / "PASS.json",
        candidate_dir / "NOT_EVALUATED.json", candidate_dir / "SHA256.tsv",
    ):
        if path.exists():
            raise FileExistsError(f"refusing to overwrite {path}")
    if not source_tmp.is_file():
        raise FileNotFoundError(source_tmp)

    hashes = {"tpr": source["tpr_sha256"], "xtc": source["xtc_sha256"]}
    base.SELECTED_TIME_PS = float(source["time_ps"])
    base.EXPECTED_SOURCE_HASHES = hashes
    base_result = base.audit_frame(
        source["tpr"], source["xtc"], source_tmp, SOURCE_NDX,
        source["time_ps"], hashes,
    )

    import MDAnalysis as mda

    topology = mda.Universe(str(source["tpr"]), str(source["xtc"]))
    ts = topology.trajectory[int(source["frame_index"])]
    if abs(float(ts.time) - float(source["time_ps"])) > TIME_TOLERANCE_PS:
        raise ValueError(
            f"frame {source['frame_index']} time {ts.time} != {source['time_ps']}"
        )
    extracted = mda.Universe(str(source_tmp))
    atoms, box = extracted.atoms, np.asarray(extracted.dimensions, dtype=float)
    identities = {
        "nalpha": atom_identity(atoms[N_ALPHA - 1], N_ALPHA, "THR", "N"),
        "thr267_og1": atom_identity(atoms[THR267_OG1 - 1], THR267_OG1, "THR", "OG1"),
        "transferred_hg1": atom_identity(
            atoms[TRANSFERRED_HG1 - 1], TRANSFERRED_HG1, "THR", "HG1"
        ),
        "l2_c12": atom_identity(atoms[L2_C12 - 1], L2_C12, "L2", "C12"),
        "l2_o2": atom_identity(atoms[L2_O2 - 1], L2_O2, "L2", "O2"),
        "l2_n3": atom_identity(atoms[L2_N3 - 1], L2_N3, "L2", "N3"),
    }

    source_atoms = topology.atoms
    n_bonds = {int(a.index) + 1 for a in source_atoms[N_ALPHA - 1].bonded_atoms}
    og1_bonds = {int(a.index) + 1 for a in source_atoms[THR267_OG1 - 1].bonded_atoms}
    c12_bonds = {int(a.index) + 1 for a in source_atoms[L2_C12 - 1].bonded_atoms}
    bond_graph = {
        "nalpha_hg1": TRANSFERRED_HG1 in n_bonds,
        "og1_hg1_absent": TRANSFERRED_HG1 not in og1_bonds,
        "c12_n3": L2_N3 in c12_bonds,
    }
    if not all(bond_graph.values()):
        raise ValueError(f"bond graph failed: {bond_graph}")

    nalpha, hg1, n3 = (
        atoms[N_ALPHA - 1], atoms[TRANSFERRED_HG1 - 1], atoms[L2_N3 - 1]
    )
    pt2 = {
        "nalpha_n3_nm": distance_nm(nalpha, n3, box),
        "hg1_n3_nm": distance_nm(hg1, n3, box),
        "pt2_angle_deg": angle_deg(nalpha, hg1, n3, box),
    }
    geometry = {
        "nac_distance_nm": float(base_result["joint_nac"]["c_og1_distance_nm"]),
        "nac_angle_deg": float(base_result["joint_nac"]["o_c_og1_angle_deg"]),
        **pt2,
    }
    for key in geometry:
        tolerance = ANGLE_TOLERANCE_DEG if key.endswith("_deg") else DISTANCE_TOLERANCE_NM
        if abs(geometry[key] - float(source[key])) > tolerance:
            raise ValueError(f"{key} differs from frozen reference")
    if not (
        geometry["nac_distance_nm"] <= NAC_DISTANCE_MAX_NM
        and NAC_ANGLE_MIN_DEG <= geometry["nac_angle_deg"] <= NAC_ANGLE_MAX_DEG
        and geometry["nalpha_n3_nm"] <= PT2_DONOR_ACCEPTOR_MAX_NM
        and geometry["hg1_n3_nm"] <= PT2_H_ACCEPTOR_MAX_NM
        and geometry["pt2_angle_deg"] >= PT2_ANGLE_MIN_DEG
    ):
        raise ValueError(f"NAC/PT2 gate failed: {geometry}")

    minimum_contact_nm = {
        key: float(value["distance_nm"])
        for key, value in base_result["minimum_contacts"].items()
    }
    box_nm = [
        float(value) / 10.0
        for value in base_result["box"]["extracted_A_deg"][:3]
    ]
    os.replace(source_tmp, source_gro)
    gro_sha = sha256(source_gro)
    manifest = {
        "schema_version": 1,
        "status": PASS_STATUS,
        "scientific_gate": "PASS_A1_PT2_PREORGANIZED_FRAME_EXTRACTION",
        "scientific_scope": SCOPE,
        "github_commit": commit,
        "candidate": candidate,
        "seed": source["seed"],
        "source": {
            "tpr": str(source["tpr"]), "xtc": str(source["xtc"]),
            "sha256": hashes, "time_ps": source["time_ps"],
            "frame_index": source["frame_index"],
            "time_tolerance_ps": TIME_TOLERANCE_PS,
        },
        "extracted": {
            "path": str(source_gro), "sha256": gro_sha,
            "atom_count": int(base_result["counts"]["system_atoms"]),
            "box_nm": box_nm,
        },
        "atom_identity": identities,
        "bond_graph": bond_graph,
        "geometry": geometry,
        "minimum_contact_nm": minimum_contact_nm,
        "base_gates": base_result["gates"],
    }
    write_json(candidate_dir / "manifest.json", manifest)
    write_json(candidate_dir / "PASS.json", {
        "status": PASS_STATUS,
        "scientific_scope": SCOPE,
        "candidate": candidate,
        "source_gro_sha256": gro_sha,
    })
    lines = []
    for name in ("source.gro", "manifest.json", "PASS.json"):
        lines.append(f"{sha256(candidate_dir / name)}  {name}\n")
    (candidate_dir / "SHA256.tsv").write_text("".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", choices=sorted(SOURCES), required=True)
    parser.add_argument("--candidate-dir", type=pathlib.Path, required=True)
    parser.add_argument("--github-commit", default="unknown")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--prepare", action="store_true")
    mode.add_argument("--audit", action="store_true")
    mode.add_argument("--fail", action="store_true")
    parser.add_argument("--stage", default="runner")
    parser.add_argument("--detail", default="unspecified")
    args = parser.parse_args()
    try:
        if args.prepare:
            prepare(args.candidate, args.candidate_dir)
            return 0
        if args.fail:
            fail(args.candidate, args.candidate_dir, args.stage, args.detail)
            return 1
        try:
            audit(args.candidate, args.candidate_dir, args.github_commit)
            return 0
        except Exception as exc:
            fail(
                args.candidate, args.candidate_dir,
                getattr(exc, "stage", "audit"), str(exc),
            )
            raise
    except Exception as exc:
        print(f"{type(exc).__name__}: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
