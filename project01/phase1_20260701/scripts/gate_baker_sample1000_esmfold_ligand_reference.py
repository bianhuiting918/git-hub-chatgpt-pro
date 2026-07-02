#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


AA3 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}

FIXED = [13, 14, 15, 16, 17, 54, 55, 56, 72, 73, 74, 148, 149, 150]
CAT_ATOMS = {
    15: ["CB", "OG"],
    55: ["CB", "CG", "ND1", "CD2", "CE1", "NE2"],
    73: ["CB", "CG", "OD1", "OD2"],
}
BB = ["N", "CA", "C", "O"]
KEY_PROTEIN = [
    ("Ser15_OG_His55_NE2", ("ATOM", "A", 15, "OG"), ("ATOM", "A", 55, "NE2")),
    ("His55_ND1_Asp73_OD1", ("ATOM", "A", 55, "ND1"), ("ATOM", "A", 73, "OD1")),
    ("His55_ND1_Asp73_OD2", ("ATOM", "A", 55, "ND1"), ("ATOM", "A", 73, "OD2")),
]
KEY_LIGAND = [
    ("Ser15_OG_bn1_C1", ("ATOM", "A", 15, "OG"), ("HETATM", "B", 391, "C1")),
    ("Ser15_OG_bn1_O1", ("ATOM", "A", 15, "OG"), ("HETATM", "B", 391, "O1")),
    ("Oxh16_N_bn1_O1", ("ATOM", "A", 16, "N"), ("HETATM", "B", 391, "O1")),
    ("Oxh149_N_bn1_O1", ("ATOM", "A", 149, "N"), ("HETATM", "B", 391, "O1")),
]
THRESH = {
    "active_motif_ca_rmsd_max_A": 1.0,
    "active_pair_delta_max_A": 1.0,
    "strict_fixed_backbone_rmsd_max_A": 1.0,
    "strict_catalytic_heavy_rmsd_max_A": 0.75,
    "strict_protein_key_distance_delta_max_A": 0.75,
    "strict_ligand_key_distance_delta_max_A": 0.75,
    "strict_mean_motif_plddt_min": 70.0,
    "strict_ligand_clash_min_A": 1.8,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference-pdb", required=True)
    parser.add_argument("--selected-tsv", required=True)
    parser.add_argument("--esmfold-status-tsv", required=True)
    parser.add_argument("--out-tsv", required=True)
    parser.add_argument("--accepted-tsv", required=True)
    parser.add_argument("--summary-json", required=True)
    parser.add_argument("--target-per-bin", type=int, default=10)
    parser.add_argument("--motif-plddt-min", type=float, default=70.0, help="Set below 0 to disable pLDDT/B-factor gate for non-ESMFold backbone-preserved structures")
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_pdb(path: Path) -> list[dict[str, object]]:
    atoms = []
    for line in path.read_text(errors="ignore").splitlines():
        rec = line[:6].strip()
        if rec not in {"ATOM", "HETATM"} or len(line) < 54:
            continue
        try:
            atom = {
                "rec": rec,
                "atom": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21].strip() or "A",
                "resseq": int(line[22:26]),
                "coord": np.array(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                    dtype=float,
                ),
                "bfactor": float(line[60:66]) if len(line) >= 66 else math.nan,
                "line": line,
            }
        except ValueError:
            continue
        element = (line[76:78].strip() if len(line) >= 78 else "") or str(atom["atom"])[0]
        atom["element"] = element.upper()
        atom["heavy"] = not (atom["element"] == "H" or str(atom["atom"]).startswith("H"))
        atoms.append(atom)
    return atoms


def atom_map(atoms: list[dict[str, object]]) -> dict[tuple[str, str, int, str], dict[str, object]]:
    return {
        (str(a["rec"]), str(a["chain"]), int(a["resseq"]), str(a["atom"])): a
        for a in atoms
    }


def sequence_of(atoms: list[dict[str, object]]) -> str:
    seen = []
    used = set()
    for atom in atoms:
        if atom["rec"] != "ATOM" or atom["chain"] != "A":
            continue
        key = (atom["chain"], atom["resseq"])
        if key in used:
            continue
        used.add(key)
        seen.append(AA3.get(str(atom["resname"]), "X"))
    return "".join(seen)


def kabsch(ref_xyz: np.ndarray, cand_xyz: np.ndarray) -> tuple[np.ndarray, np.ndarray, float]:
    ref_center = ref_xyz.mean(axis=0)
    cand_center = cand_xyz.mean(axis=0)
    ref0 = ref_xyz - ref_center
    cand0 = cand_xyz - cand_center
    u, _s, vt = np.linalg.svd(cand0.T @ ref0)
    rot = u @ vt
    if np.linalg.det(rot) < 0:
        u[:, -1] *= -1
        rot = u @ vt
    trans = ref_center - cand_center @ rot
    aligned = cand_xyz @ rot + trans
    rmsd = float(np.sqrt(((aligned - ref_xyz) ** 2).sum() / len(ref_xyz)))
    return rot, trans, rmsd


def transform_atoms(
    atoms: list[dict[str, object]], rot: np.ndarray, trans: np.ndarray
) -> list[dict[str, object]]:
    transformed = []
    for atom in atoms:
        copy = dict(atom)
        copy["coord"] = atom["coord"] @ rot + trans
        transformed.append(copy)
    return transformed


def collect_coords(amap: dict, keys: list[tuple[str, str, int, str]]) -> tuple[np.ndarray, list[str]]:
    coords = []
    missing = []
    for key in keys:
        if key not in amap:
            missing.append(":".join(map(str, key)))
        else:
            coords.append(amap[key]["coord"])
    return np.array(coords, dtype=float), missing


def rmsd_for_keys(ref_amap: dict, cand_amap: dict, keys: list[tuple[str, str, int, str]]):
    ref = []
    cand = []
    missing = []
    for key in keys:
        if key not in ref_amap or key not in cand_amap:
            missing.append(":".join(map(str, key)))
        else:
            ref.append(ref_amap[key]["coord"])
            cand.append(cand_amap[key]["coord"])
    if missing or not ref:
        return None, missing
    ref_arr = np.array(ref, dtype=float)
    cand_arr = np.array(cand, dtype=float)
    return float(np.sqrt(((ref_arr - cand_arr) ** 2).sum() / len(ref_arr))), []


def motif_pair_deltas(ref_amap: dict, cand_amap_unaligned: dict) -> tuple[float | None, float | None, list[str]]:
    values = []
    missing = []
    for i, res_i in enumerate(FIXED):
        for res_j in FIXED[i + 1 :]:
            key_i = ("ATOM", "A", res_i, "CA")
            key_j = ("ATOM", "A", res_j, "CA")
            if (
                key_i not in ref_amap
                or key_j not in ref_amap
                or key_i not in cand_amap_unaligned
                or key_j not in cand_amap_unaligned
            ):
                missing.append(f"A{res_i}-A{res_j}")
                continue
            ref_dist = np.linalg.norm(ref_amap[key_i]["coord"] - ref_amap[key_j]["coord"])
            cand_dist = np.linalg.norm(
                cand_amap_unaligned[key_i]["coord"] - cand_amap_unaligned[key_j]["coord"]
            )
            values.append(abs(float(ref_dist - cand_dist)))
    if not values:
        return None, None, missing
    return max(values), float(np.mean(values)), missing


def distance_deltas(ref_amap: dict, cand_amap: dict, pairs: list[tuple[str, tuple, tuple]]):
    values = {}
    missing = []
    max_delta = 0.0
    for name, protein_key, partner_key in pairs:
        if protein_key not in ref_amap or partner_key not in ref_amap or protein_key not in cand_amap:
            missing.append(name)
            continue
        ref_dist = float(np.linalg.norm(ref_amap[protein_key]["coord"] - ref_amap[partner_key]["coord"]))
        cand_dist = float(
            np.linalg.norm(cand_amap[protein_key]["coord"] - ref_amap[partner_key]["coord"])
        )
        delta = abs(cand_dist - ref_dist)
        values[name] = round(delta, 4)
        max_delta = max(max_delta, delta)
    return max_delta, values, missing


def ligand_clash_min(cand_atoms_aligned: list[dict[str, object]], ref_atoms: list[dict[str, object]]):
    ligand = [
        atom
        for atom in ref_atoms
        if atom["rec"] == "HETATM" and str(atom["resname"]).lower() == "bn1" and atom["heavy"]
    ]
    protein = [
        atom
        for atom in cand_atoms_aligned
        if atom["rec"] == "ATOM"
        and atom["heavy"]
        and not (atom["chain"] == "A" and int(atom["resseq"]) in FIXED)
    ]
    min_dist = None
    for lig_atom in ligand:
        for prot_atom in protein:
            dist = float(np.linalg.norm(lig_atom["coord"] - prot_atom["coord"]))
            if min_dist is None or dist < min_dist:
                min_dist = dist
    return min_dist


def evaluate(row: dict[str, str], status: dict[str, str], ref_atoms: list[dict[str, object]], ref_amap: dict):
    sample_id = row["sample_id"]
    pred = Path(status.get("output_pdb") or status.get("pdb") or row.get("postseq_predicted_pdb", ""))
    if status.get("status") not in {"OK", "SKIP_EXISTS"} or not pred.exists() or pred.stat().st_size < 1000:
        return {
            "sample_id": sample_id,
            "bin": row.get("bin", ""),
            "identity": row.get("identity", ""),
            "mutation_count": row.get("mutation_count", ""),
            "esmfold_status": status.get("status", "NOT_EVALUATED"),
            "active_pocket_ca_gate": "NOT_EVALUATED",
            "strict_sidechain_ligand_gate": "NOT_EVALUATED",
            "fail_reasons": "missing_or_failed_esmfold",
            "predicted_pdb": str(pred),
        }

    cand_atoms = parse_pdb(pred)
    cand_amap0 = atom_map(cand_atoms)
    reasons = []
    strict_reasons = []
    active_reasons = []

    seq_match = sequence_of(cand_atoms) == row.get("sequence", "")
    if not seq_match:
        reasons.append("sequence_mismatch")

    fixed_identity_mutations = []
    for res_id in FIXED:
        ca_key = ("ATOM", "A", res_id, "CA")
        if ca_key not in ref_amap or ca_key not in cand_amap0:
            fixed_identity_mutations.append(f"A{res_id}:missing")
        elif ref_amap[ca_key]["resname"] != cand_amap0[ca_key]["resname"]:
            fixed_identity_mutations.append(
                f"A{res_id}:{ref_amap[ca_key]['resname']}>{cand_amap0[ca_key]['resname']}"
            )
    if fixed_identity_mutations:
        reasons.append("fixed_identity_mutation")

    ca_keys = [("ATOM", "A", res_id, "CA") for res_id in FIXED]
    ref_xyz, missing_ref = collect_coords(ref_amap, ca_keys)
    cand_xyz, missing_cand = collect_coords(cand_amap0, ca_keys)
    if missing_ref or missing_cand or len(ref_xyz) != len(ca_keys):
        reasons.append("missing_fixed_ca")
        motif_ca_rmsd = None
        pair_max = None
        pair_mean = None
        cand_atoms_aligned = cand_atoms
        cand_amap = cand_amap0
    else:
        rot, trans, motif_ca_rmsd = kabsch(ref_xyz, cand_xyz)
        pair_max, pair_mean, _missing_pairs = motif_pair_deltas(ref_amap, cand_amap0)
        cand_atoms_aligned = transform_atoms(cand_atoms, rot, trans)
        cand_amap = atom_map(cand_atoms_aligned)

    if motif_ca_rmsd is None or motif_ca_rmsd > THRESH["active_motif_ca_rmsd_max_A"]:
        active_reasons.append("motif_ca_rmsd_gt_1A")
    if pair_max is None or pair_max > THRESH["active_pair_delta_max_A"]:
        active_reasons.append("motif_pair_delta_gt_1A")

    bb_keys = [("ATOM", "A", res_id, atom_name) for res_id in FIXED for atom_name in BB]
    cat_keys = [
        ("ATOM", "A", res_id, atom_name)
        for res_id, atom_names in CAT_ATOMS.items()
        for atom_name in atom_names
    ]
    bb_rmsd, _bb_missing = rmsd_for_keys(ref_amap, cand_amap, bb_keys)
    cat_rmsd, _cat_missing = rmsd_for_keys(ref_amap, cand_amap, cat_keys)
    protein_max, _protein_values, protein_missing = distance_deltas(ref_amap, cand_amap, KEY_PROTEIN)
    ligand_max, _ligand_values, ligand_missing = distance_deltas(ref_amap, cand_amap, KEY_LIGAND)
    motif_plddt = [
        float(atom["bfactor"])
        for atom in cand_atoms
        if atom["rec"] == "ATOM" and atom["chain"] == "A" and int(atom["resseq"]) in FIXED
    ]
    mean_plddt = float(np.mean(motif_plddt)) if motif_plddt else None
    clash = ligand_clash_min(cand_atoms_aligned, ref_atoms)

    strict_reasons.extend(reasons)
    if bb_rmsd is None or bb_rmsd > THRESH["strict_fixed_backbone_rmsd_max_A"]:
        strict_reasons.append("fixed_backbone_rmsd_gt_1A")
    if cat_rmsd is None or cat_rmsd > THRESH["strict_catalytic_heavy_rmsd_max_A"]:
        strict_reasons.append("catalytic_heavy_rmsd_gt_0p75A")
    if protein_missing or protein_max > THRESH["strict_protein_key_distance_delta_max_A"]:
        strict_reasons.append("protein_key_distance_delta_gt_0p75A")
    if ligand_missing or ligand_max > THRESH["strict_ligand_key_distance_delta_max_A"]:
        strict_reasons.append("ligand_key_distance_delta_gt_0p75A")
    plddt_min = THRESH["strict_mean_motif_plddt_min"]
    if plddt_min >= 0 and (mean_plddt is None or mean_plddt < plddt_min):
        strict_reasons.append("motif_plddt_lt_threshold")
    if clash is None or clash < THRESH["strict_ligand_clash_min_A"]:
        strict_reasons.append("ligand_clash_lt_1p8A")

    active_pass = not active_reasons and not reasons
    strict_pass = not strict_reasons
    return {
        "sample_id": sample_id,
        "bin": row.get("bin", ""),
        "identity": row.get("identity", ""),
        "mutation_count": row.get("mutation_count", ""),
        "esmfold_status": status.get("status", ""),
        "active_pocket_ca_gate": "PASS" if active_pass else "FAIL",
        "strict_sidechain_ligand_gate": "PASS" if strict_pass else "FAIL",
        "fail_reasons": ";".join(strict_reasons or active_reasons),
        "active_fail_reasons": ";".join(active_reasons),
        "motif_ca_rmsd_A": "" if motif_ca_rmsd is None else f"{motif_ca_rmsd:.4f}",
        "motif_pair_max_delta_A": "" if pair_max is None else f"{pair_max:.4f}",
        "motif_pair_mean_delta_A": "" if pair_mean is None else f"{pair_mean:.4f}",
        "fixed_backbone_rmsd_A": "" if bb_rmsd is None else f"{bb_rmsd:.4f}",
        "catalytic_heavy_rmsd_A": "" if cat_rmsd is None else f"{cat_rmsd:.4f}",
        "protein_key_delta_max_A": f"{protein_max:.4f}",
        "ligand_key_delta_max_A": f"{ligand_max:.4f}",
        "motif_mean_plddt": "" if mean_plddt is None else f"{mean_plddt:.2f}",
        "ligand_clash_min_A": "" if clash is None else f"{clash:.4f}",
        "predicted_pdb": str(pred),
    }


def main() -> None:
    args = parse_args()
    THRESH["strict_mean_motif_plddt_min"] = args.motif_plddt_min
    reference_pdb = Path(args.reference_pdb)
    selected_rows = read_tsv(Path(args.selected_tsv))
    status_rows = read_tsv(Path(args.esmfold_status_tsv))
    status_by_id = {row["sample_id"]: row for row in status_rows}
    ref_atoms = parse_pdb(reference_pdb)
    ref_amap = atom_map(ref_atoms)

    rows = [evaluate(row, status_by_id.get(row["sample_id"], {}), ref_atoms, ref_amap) for row in selected_rows]
    fields = list(rows[0].keys()) if rows else ["sample_id"]
    out_tsv = Path(args.out_tsv)
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    with out_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    by_bin = defaultdict(lambda: {"evaluated": 0, "not_evaluated": 0, "active_pass": 0, "strict_pass": 0})
    for row in rows:
        bucket = by_bin[row.get("bin", "")]
        if row["active_pocket_ca_gate"] == "NOT_EVALUATED":
            bucket["not_evaluated"] += 1
        else:
            bucket["evaluated"] += 1
        if row["active_pocket_ca_gate"] == "PASS":
            bucket["active_pass"] += 1
        if row["strict_sidechain_ligand_gate"] == "PASS":
            bucket["strict_pass"] += 1

    accepted = []
    selected_by_bin = defaultdict(int)
    strict_rows = [row for row in rows if row["strict_sidechain_ligand_gate"] == "PASS"]
    strict_rows.sort(
        key=lambda row: (
            int(row["bin"]),
            float(row["catalytic_heavy_rmsd_A"] or 999.0),
            float(row["ligand_key_delta_max_A"] or 999.0),
        )
    )
    for row in strict_rows:
        if selected_by_bin[row["bin"]] < args.target_per_bin:
            accepted.append(row)
            selected_by_bin[row["bin"]] += 1
    for bin_name, count in selected_by_bin.items():
        by_bin[bin_name]["strict_selected"] = count

    accepted_tsv = Path(args.accepted_tsv)
    with accepted_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(accepted)

    summary = {
        "status": "DONE",
        "evaluated_universe": (
            "Rows in selected-tsv; rows with absent or failed ESMFold output are "
            "NOT_EVALUATED, not structural FAIL."
        ),
        "reference_pdb": str(reference_pdb),
        "selected_tsv": args.selected_tsv,
        "esmfold_status_tsv": args.esmfold_status_tsv,
        "out_tsv": str(out_tsv),
        "accepted_tsv": str(accepted_tsv),
        "rows_total": len(rows),
        "counts_by_bin": {key: dict(value) for key, value in sorted(by_bin.items())},
        "target_per_bin": args.target_per_bin,
        "target_complete": all(
            by_bin[str(bin_name)].get("strict_selected", 0) >= args.target_per_bin
            for bin_name in [90, 80, 70, 60, 50]
        ),
        "thresholds": THRESH,
        "gate_notes": (
            "Apo ESMFold protein is motif-aligned to the reference ligand-bound "
            "sample1000_refined_0 structure. Ligand geometry uses the reference "
            "bn1 coordinates rather than requiring HETATM ligand records in the "
            "ESMFold PDB."
        ),
    }
    summary_json = Path(args.summary_json)
    summary_json.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()



