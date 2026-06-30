#!/usr/bin/env python3
"""Evaluate post-sequence structure predictions before PLACER.

This gate is for ESMFold/ColabFold-style protein-only predictions. It first
aligns candidate to reference with common backbone atoms, then evaluates the
catalytic residues and fixed pocket residues in the aligned coordinate frame.
Ligand checks are optional because sequence-to-structure predictors usually do
not output the ligand.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable

import numpy as np


BACKBONE_ATOMS = {"N", "CA", "C", "O"}
CATALYTIC_RESIDUES = {("A", 95), ("A", 126), ("A", 128)}
CATALYTIC_SIDECHAIN_ATOMS = {
    ("A", 95): {"ND1", "CE1", "NE2", "CD2", "CG"},
    ("A", 126): {"OG", "CB", "CA", "C", "N"},
    ("A", 128): {"OG", "CB", "CA", "C", "N"},
}
KEY_DISTANCES = {
    "His95_ND1_to_Ser128_OG": (("ATOM", "A", 95, "ND1"), ("ATOM", "A", 128, "OG")),
    "His95_NE2_to_Ser126_OG": (("ATOM", "A", 95, "NE2"), ("ATOM", "A", 126, "OG")),
    "His95_NE2_to_Ser128_OG": (("ATOM", "A", 95, "NE2"), ("ATOM", "A", 128, "OG")),
    "Ser126_OG_to_ligand_C1": (("ATOM", "A", 126, "OG"), ("HETATM", "X", 1, "C1")),
    "Ser128_OG_to_ligand_C1": (("ATOM", "A", 128, "OG"), ("HETATM", "X", 1, "C1")),
}
DEFAULT_THRESHOLDS = {
    "global_backbone_rmsd_max_A": 2.50,
    "fixed_backbone_rmsd_max_A": 1.00,
    "catalytic_heavy_rmsd_max_A": 0.75,
    "protein_key_distance_abs_delta_max_A": 0.75,
    "ligand_heavy_rmsd_max_A": 0.75,
    "ligand_key_distance_abs_delta_max_A": 0.50,
}


def parse_residue_tokens(text: str) -> set[tuple[str, int]]:
    residues = set()
    for tok in text.replace(",", " ").split():
        residues.add((tok[0], int(tok[1:])))
    return residues


def parse_pdb(path: Path) -> list[dict]:
    atoms = []
    with path.open(errors="ignore") as handle:
        for line in handle:
            rec = line[:6].strip()
            if rec not in {"ATOM", "HETATM"}:
                continue
            try:
                atom = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21].strip() or "A"
                resseq = int(line[22:26])
                coord = np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                    dtype=float,
                )
                element = (line[76:78].strip() or atom[0]).upper()
                bfactor = float(line[60:66])
            except Exception:
                continue
            atoms.append(
                {
                    "rec": rec,
                    "atom": atom,
                    "resname": resname,
                    "chain": chain,
                    "resseq": resseq,
                    "coord": coord,
                    "element": element,
                    "heavy": not (element == "H" or atom.startswith("H")),
                    "bfactor": bfactor,
                }
            )
    return atoms


def atom_map(atoms: Iterable[dict]) -> dict[tuple[str, str, int, str], dict]:
    return {(a["rec"], a["chain"], a["resseq"], a["atom"]): a for a in atoms}


def select_atoms(
    atoms: Iterable[dict],
    residues: set[tuple[str, int]] | None = None,
    backbone_only: bool = False,
    catalytic_heavy: bool = False,
    ligand: bool = False,
) -> dict[tuple[str, int, str], dict]:
    selected = {}
    for a in atoms:
        if ligand:
            if not (
                a["rec"] == "HETATM"
                and a["resname"].lower() == "bu2"
                and a["heavy"]
            ):
                continue
        else:
            if a["rec"] != "ATOM":
                continue
            if residues is not None and (a["chain"], a["resseq"]) not in residues:
                continue
            if backbone_only and a["atom"] not in BACKBONE_ATOMS:
                continue
            if catalytic_heavy:
                allowed = CATALYTIC_SIDECHAIN_ATOMS.get((a["chain"], a["resseq"]))
                if allowed is None or a["atom"] not in allowed or not a["heavy"]:
                    continue
        selected[(a["chain"], a["resseq"], a["atom"])] = a
    return selected


def kabsch_transform(ref_xyz: np.ndarray, cand_xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    ref_centroid = ref_xyz.mean(axis=0)
    cand_centroid = cand_xyz.mean(axis=0)
    ref_centered = ref_xyz - ref_centroid
    cand_centered = cand_xyz - cand_centroid
    h = cand_centered.T @ ref_centered
    u, _s, vt = np.linalg.svd(h)
    rot = u @ vt
    if np.linalg.det(rot) < 0:
        u[:, -1] *= -1
        rot = u @ vt
    return rot, ref_centroid - cand_centroid @ rot


def apply_transform(atoms: list[dict], rot: np.ndarray, trans: np.ndarray) -> list[dict]:
    out = []
    for a in atoms:
        b = dict(a)
        b["coord"] = a["coord"] @ rot + trans
        out.append(b)
    return out


def rmsd_for(ref_sel: dict, cand_sel: dict) -> tuple[float | None, int, list[str]]:
    common = sorted(set(ref_sel) & set(cand_sel))
    missing = sorted(set(ref_sel) - set(cand_sel))
    if not common:
        return None, 0, ["|".join(map(str, x)) for x in missing[:50]]
    sq = []
    for key in common:
        delta = ref_sel[key]["coord"] - cand_sel[key]["coord"]
        sq.append(float(delta @ delta))
    return math.sqrt(sum(sq) / len(sq)), len(common), ["|".join(map(str, x)) for x in missing[:50]]


def mean_bfactor(sel: dict) -> float | None:
    vals = [a["bfactor"] for a in sel.values()]
    return float(sum(vals) / len(vals)) if vals else None


def distance(atom_a: dict, atom_b: dict) -> float:
    delta = atom_a["coord"] - atom_b["coord"]
    return float(math.sqrt(delta @ delta))


def key_distances(atoms: list[dict], include_ligand: bool) -> tuple[dict[str, float], list[str]]:
    amap = atom_map(atoms)
    out = {}
    missing = []
    for name, (left, right) in KEY_DISTANCES.items():
        if not include_ligand and "ligand" in name:
            continue
        a = amap.get(left)
        b = amap.get(right)
        if a is None or b is None:
            missing.append(name)
            continue
        out[name] = distance(a, b)
    return out, missing


def evaluate(reference: Path, candidate: Path, fixed_residues: set[tuple[str, int]], require_ligand: bool) -> dict:
    ref_atoms = parse_pdb(reference)
    cand_atoms = parse_pdb(candidate)

    ref_align = select_atoms(ref_atoms, backbone_only=True)
    cand_align = select_atoms(cand_atoms, backbone_only=True)
    common_align = sorted(set(ref_align) & set(cand_align))
    if len(common_align) < 12:
        return {
            "status": "NOT_EVALUATED",
            "evaluation_status": "FAILED_BEFORE_EVALUATION",
            "fail_reasons": ["insufficient_common_backbone_atoms_for_alignment"],
            "reference": str(reference),
            "candidate": str(candidate),
            "thresholds": DEFAULT_THRESHOLDS,
        }

    ref_xyz = np.array([ref_align[k]["coord"] for k in common_align])
    cand_xyz = np.array([cand_align[k]["coord"] for k in common_align])
    rot, trans = kabsch_transform(ref_xyz, cand_xyz)
    cand_aligned = apply_transform(cand_atoms, rot, trans)

    ref_global = select_atoms(ref_atoms, backbone_only=True)
    cand_global = select_atoms(cand_aligned, backbone_only=True)
    ref_fixed = select_atoms(ref_atoms, residues=fixed_residues, backbone_only=True)
    cand_fixed = select_atoms(cand_aligned, residues=fixed_residues, backbone_only=True)
    ref_cat = select_atoms(ref_atoms, residues=CATALYTIC_RESIDUES, catalytic_heavy=True)
    cand_cat = select_atoms(cand_aligned, residues=CATALYTIC_RESIDUES, catalytic_heavy=True)
    ref_lig = select_atoms(ref_atoms, ligand=True)
    cand_lig = select_atoms(cand_aligned, ligand=True)
    candidate_has_ligand = bool(cand_lig)

    global_rmsd, global_n, global_missing = rmsd_for(ref_global, cand_global)
    fixed_rmsd, fixed_n, fixed_missing = rmsd_for(ref_fixed, cand_fixed)
    cat_rmsd, cat_n, cat_missing = rmsd_for(ref_cat, cand_cat)
    lig_rmsd, lig_n, lig_missing = rmsd_for(ref_lig, cand_lig)

    ref_kd, ref_kd_missing = key_distances(ref_atoms, include_ligand=candidate_has_ligand or require_ligand)
    cand_kd, cand_kd_missing = key_distances(cand_aligned, include_ligand=candidate_has_ligand or require_ligand)
    kd_delta = {k: cand_kd[k] - ref_kd[k] for k in sorted(ref_kd) if k in cand_kd}
    protein_delta = [abs(v) for k, v in kd_delta.items() if "ligand" not in k]
    ligand_delta = [abs(v) for k, v in kd_delta.items() if "ligand" in k]

    reasons = []
    if global_rmsd is None or global_rmsd > DEFAULT_THRESHOLDS["global_backbone_rmsd_max_A"]:
        reasons.append("global_backbone_rmsd")
    if fixed_rmsd is None or fixed_rmsd > DEFAULT_THRESHOLDS["fixed_backbone_rmsd_max_A"] or fixed_missing:
        reasons.append("fixed_backbone_rmsd_or_missing")
    if cat_rmsd is None or cat_rmsd > DEFAULT_THRESHOLDS["catalytic_heavy_rmsd_max_A"] or cat_missing:
        reasons.append("catalytic_heavy_rmsd_or_missing")
    if protein_delta and max(protein_delta) > DEFAULT_THRESHOLDS["protein_key_distance_abs_delta_max_A"]:
        reasons.append("protein_key_distance_delta")
    if require_ligand and not candidate_has_ligand:
        reasons.append("ligand_required_but_missing")
    if candidate_has_ligand:
        if lig_rmsd is None or lig_rmsd > DEFAULT_THRESHOLDS["ligand_heavy_rmsd_max_A"] or lig_missing:
            reasons.append("ligand_heavy_rmsd_or_missing")
        if ligand_delta and max(ligand_delta) > DEFAULT_THRESHOLDS["ligand_key_distance_abs_delta_max_A"]:
            reasons.append("ligand_key_distance_delta")
        if cand_kd_missing:
            reasons.append("missing_key_distance_atoms")

    return {
        "status": "PASS" if not reasons else "FAIL",
        "evaluation_status": "EVALUATED_WITH_OUTPUT",
        "reference": str(reference),
        "candidate": str(candidate),
        "alignment": {
            "method": "kabsch",
            "scope": "all_common_backbone_atoms",
            "atom_pairs": len(common_align),
        },
        "thresholds": DEFAULT_THRESHOLDS,
        "metrics": {
            "global_backbone_rmsd_A": global_rmsd,
            "global_backbone_atom_pairs": global_n,
            "fixed_backbone_rmsd_A": fixed_rmsd,
            "fixed_backbone_atom_pairs": fixed_n,
            "catalytic_heavy_rmsd_A": cat_rmsd,
            "catalytic_heavy_atom_pairs": cat_n,
            "candidate_has_ligand": candidate_has_ligand,
            "ligand_heavy_rmsd_A": lig_rmsd,
            "ligand_heavy_atom_pairs": lig_n,
            "reference_key_distances_A": ref_kd,
            "candidate_key_distances_A": cand_kd,
            "key_distance_deltas_A": kd_delta,
            "max_abs_protein_key_distance_delta_A": max(protein_delta) if protein_delta else None,
            "max_abs_ligand_key_distance_delta_A": max(ligand_delta) if ligand_delta else None,
            "fixed_backbone_mean_bfactor_or_plddt": mean_bfactor(cand_fixed),
            "catalytic_mean_bfactor_or_plddt": mean_bfactor(cand_cat),
        },
        "missing": {
            "global_backbone_atoms": global_missing[:20],
            "fixed_backbone_atoms": fixed_missing[:20],
            "catalytic_atoms": cat_missing[:20],
            "ligand_atoms": lig_missing[:20],
            "reference_key_distances": ref_kd_missing,
            "candidate_key_distances": cand_kd_missing,
        },
        "fail_reasons": reasons,
        "note": (
            "Protein-only predictors are allowed to omit ligand at this entrance gate. "
            "Ligand/key-distance gates become mandatory at PLACER crop and completion stages."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True)
    parser.add_argument("--candidate", required=True)
    parser.add_argument("--fixed-residues", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--require-ligand", action="store_true")
    args = parser.parse_args()

    fixed_residues = parse_residue_tokens(Path(args.fixed_residues).read_text())
    result = evaluate(
        reference=Path(args.reference),
        candidate=Path(args.candidate),
        fixed_residues=fixed_residues,
        require_ligand=args.require_ligand,
    )
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(result, indent=2, sort_keys=True))
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
