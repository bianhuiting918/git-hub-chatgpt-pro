#!/usr/bin/env python3
"""Two-layer PLACER crop gate for Project 01 phase1.

The PLACER gate intentionally reuses the same protein/pocket criteria as the
post-sequence entrance gate, then adds ligand/reaction-pose checks. A conformer
is strict-pass only when both layers pass.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Iterable

import numpy as np


HEAVY_LIGAND_ATOMS = [
    "C1",
    "O1",
    "O2",
    "C2",
    "C3",
    "C4",
    "C5",
    "C6",
    "C7",
    "C8",
    "O3",
    "C9",
    "C11",
    "O4",
    "C12",
    "C10",
    "C13",
    "C14",
]

FIXED_RESIDUES = {95, 126, 128}
CATALYTIC_RESIDUES = {95, 126, 128}

PROTEIN_KEY_PAIRS = [
    ("His95_ND1_to_Ser128_OG", ("A", 95, "HIS", "ND1"), ("A", 128, "SER", "OG")),
    ("His95_NE2_to_Ser126_OG", ("A", 95, "HIS", "NE2"), ("A", 126, "SER", "OG")),
    ("His95_NE2_to_Ser128_OG", ("A", 95, "HIS", "NE2"), ("A", 128, "SER", "OG")),
]

LIGAND_KEY_PAIRS = [
    ("Ser128_OG_to_bu2_C1", ("A", 128, "SER", "OG"), ("X", 1, "bu2", "C1")),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", required=True)
    parser.add_argument(
        "--run-dir",
        action="append",
        nargs=2,
        metavar=("LABEL", "DIR"),
        required=True,
        help="Repeatable pair: run label and directory containing *_model.pdb files.",
    )
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--global-backbone-rmsd-max-a", type=float, default=2.50)
    parser.add_argument("--fixed-backbone-rmsd-max-a", type=float, default=1.00)
    parser.add_argument("--catalytic-heavy-rmsd-max-a", type=float, default=0.75)
    parser.add_argument("--protein-key-distance-abs-delta-max-a", type=float, default=0.75)
    parser.add_argument("--ligand-heavy-rmsd-max-a", type=float, default=0.75)
    parser.add_argument("--ser128-og-to-bu2-c1-abs-delta-max-a", type=float, default=0.50)
    parser.add_argument("--minimum-crop-pass-conformers-per-sequence", type=int, default=10)
    parser.add_argument("--expected-conformers-per-sequence", type=int, default=50)
    return parser.parse_args()


def parse_atom_line(line: str) -> dict:
    return {
        "record": line[:6].strip(),
        "atom": line[12:16].strip(),
        "residue": line[17:20].strip(),
        "chain": line[21].strip(),
        "residue_index": int(line[22:26]),
        "coord": np.array(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])],
            dtype=float,
        ),
        "element": (line[76:78].strip() if len(line) >= 78 else line[12:16].strip()[0]).strip(),
    }


def parse_models(path: Path) -> list[tuple[str, list[dict]]]:
    models: list[tuple[str, list[dict]]] = []
    current: list[dict] = []
    model_id: str | None = None
    for line in path.read_text().splitlines():
        if line.startswith("MODEL"):
            if current:
                models.append((model_id or str(len(models) + 1), current))
            current = []
            parts = line.split()
            model_id = parts[1] if len(parts) > 1 else str(len(models) + 1)
        elif line.startswith("ENDMDL"):
            models.append((model_id or str(len(models) + 1), current))
            current = []
            model_id = None
        elif line.startswith(("ATOM", "HETATM")):
            current.append(parse_atom_line(line))
    if current:
        models.append((model_id or str(len(models) + 1), current))
    return models


def atom_map(atoms: Iterable[dict]) -> dict[tuple[str, int, str, str], np.ndarray]:
    return {
        (atom["chain"], atom["residue_index"], atom["residue"], atom["atom"]): atom["coord"]
        for atom in atoms
    }


def residue_heavy_atoms(
    atoms_by_key: dict[tuple[str, int, str, str], np.ndarray], residues: set[int]
) -> dict[tuple[str, int, str, str], np.ndarray]:
    return {
        key: coord
        for key, coord in atoms_by_key.items()
        if key[0] == "A" and key[1] in residues and not key[3].startswith("H")
    }


def kabsch(mobile: np.ndarray, target: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    mobile_center = mobile.mean(axis=0)
    target_center = target.mean(axis=0)
    covariance = (mobile - mobile_center).T @ (target - target_center)
    v_mat, _, wt_mat = np.linalg.svd(covariance)
    handedness = np.sign(np.linalg.det(v_mat @ wt_mat))
    rotation = v_mat @ np.diag([1.0, 1.0, handedness]) @ wt_mat
    return rotation, mobile_center, target_center


def apply_transform(
    coord: np.ndarray, rotation: np.ndarray, mobile_center: np.ndarray, target_center: np.ndarray
) -> np.ndarray:
    return (coord - mobile_center) @ rotation + target_center


def distance(coord_a: np.ndarray, coord_b: np.ndarray) -> float:
    return float(np.linalg.norm(coord_a - coord_b))


def rmsd(coords_a: Iterable[np.ndarray], coords_b: Iterable[np.ndarray]) -> float:
    arr_a = np.array(list(coords_a))
    arr_b = np.array(list(coords_b))
    return float(np.sqrt(((arr_a - arr_b) ** 2).sum(axis=1).mean()))


def evaluate_model(
    *,
    run_label: str,
    sample_id: str,
    bin_id: str,
    model_id: str,
    source_pdb: Path,
    atoms: list[dict],
    ref_map: dict[tuple[str, int, str, str], np.ndarray],
    ref_backbone: dict[tuple[str, int, str, str], np.ndarray],
    ref_fixed_backbone: dict[tuple[str, int, str, str], np.ndarray],
    ref_catalytic: dict[tuple[str, int, str, str], np.ndarray],
    ref_protein_key_distances: dict[str, float],
    ref_ligand_key_distances: dict[str, float],
    thresholds: dict[str, float | int],
) -> dict:
    row = {
        "run_label": run_label,
        "sample_id": sample_id,
        "bin": bin_id,
        "model_id": model_id,
        "source_pdb": str(source_pdb),
    }
    candidate_map = atom_map(atoms)
    candidate_backbone = {
        key: coord
        for key, coord in candidate_map.items()
        if key[0] == "A" and key[3] in {"N", "CA", "C", "O"}
    }
    common_backbone = sorted(set(ref_backbone) & set(candidate_backbone))
    if len(common_backbone) < 12:
        row.update(
            {
                "status": "NOT_EVALUATED_INSUFFICIENT_BACKBONE",
                "inherited_postseq_protein_gate": "NOT_EVALUATED",
                "ligand_reaction_pose_gate": "NOT_EVALUATED",
                "global_backbone_atom_pairs": len(common_backbone),
                "fail_reasons": "insufficient_backbone_pairs",
            }
        )
        return row

    mobile = np.array([candidate_backbone[key] for key in common_backbone])
    target = np.array([ref_backbone[key] for key in common_backbone])
    rotation, mobile_center, target_center = kabsch(mobile, target)

    def align(coord: np.ndarray) -> np.ndarray:
        return apply_transform(coord, rotation, mobile_center, target_center)

    global_rmsd = rmsd(
        [align(candidate_backbone[key]) for key in common_backbone],
        [ref_backbone[key] for key in common_backbone],
    )

    fixed_common = sorted(set(ref_fixed_backbone) & set(candidate_backbone))
    fixed_rmsd = (
        rmsd(
            [align(candidate_backbone[key]) for key in fixed_common],
            [ref_fixed_backbone[key] for key in fixed_common],
        )
        if fixed_common
        else None
    )

    candidate_catalytic = residue_heavy_atoms(candidate_map, CATALYTIC_RESIDUES)
    catalytic_common = sorted(set(ref_catalytic) & set(candidate_catalytic))
    catalytic_rmsd = (
        rmsd(
            [align(candidate_catalytic[key]) for key in catalytic_common],
            [ref_catalytic[key] for key in catalytic_common],
        )
        if catalytic_common
        else None
    )

    protein_deltas: dict[str, float] = {}
    protein_missing: list[str] = []
    for label, atom_a, atom_b in PROTEIN_KEY_PAIRS:
        if atom_a in candidate_map and atom_b in candidate_map and label in ref_protein_key_distances:
            protein_deltas[label] = (
                distance(align(candidate_map[atom_a]), align(candidate_map[atom_b]))
                - ref_protein_key_distances[label]
            )
        else:
            protein_missing.append(label)
    max_abs_protein_delta = max((abs(value) for value in protein_deltas.values()), default=None)

    ligand_ref: list[np.ndarray] = []
    ligand_candidate: list[np.ndarray] = []
    ligand_missing: list[str] = []
    for atom_name in HEAVY_LIGAND_ATOMS:
        atom_key = ("X", 1, "bu2", atom_name)
        if atom_key not in ref_map or atom_key not in candidate_map:
            ligand_missing.append(atom_name)
        else:
            ligand_ref.append(ref_map[atom_key])
            ligand_candidate.append(align(candidate_map[atom_key]))
    ligand_rmsd = rmsd(ligand_candidate, ligand_ref) if not ligand_missing else None

    ligand_deltas: dict[str, float] = {}
    ligand_key_missing: list[str] = []
    for label, atom_a, atom_b in LIGAND_KEY_PAIRS:
        if atom_a in candidate_map and atom_b in candidate_map and label in ref_ligand_key_distances:
            ligand_deltas[label] = (
                distance(align(candidate_map[atom_a]), align(candidate_map[atom_b]))
                - ref_ligand_key_distances[label]
            )
        else:
            ligand_key_missing.append(label)

    protein_reasons: list[str] = []
    if not fixed_common:
        protein_reasons.append("fixed_backbone_atoms")
    if not catalytic_common:
        protein_reasons.append("catalytic_heavy_atoms")
    protein_reasons.extend(protein_missing)
    if global_rmsd > thresholds["global_backbone_rmsd_max_A"]:
        protein_reasons.append("global_backbone_rmsd")
    if fixed_rmsd is None or fixed_rmsd > thresholds["fixed_backbone_rmsd_max_A"]:
        protein_reasons.append("fixed_backbone_rmsd")
    if catalytic_rmsd is None or catalytic_rmsd > thresholds["catalytic_heavy_rmsd_max_A"]:
        protein_reasons.append("catalytic_heavy_rmsd")
    if (
        max_abs_protein_delta is None
        or max_abs_protein_delta > thresholds["protein_key_distance_abs_delta_max_A"]
    ):
        protein_reasons.append("protein_key_distance_delta")
    inherited_gate = "PASS" if not protein_reasons else "FAIL"

    ligand_reasons: list[str] = []
    if ligand_missing:
        ligand_reasons.append("missing_ligand_atoms")
    ligand_reasons.extend(ligand_key_missing)
    if ligand_rmsd is None or ligand_rmsd > thresholds["ligand_heavy_rmsd_max_A"]:
        ligand_reasons.append("ligand_heavy_rmsd")
    ser_delta = ligand_deltas.get("Ser128_OG_to_bu2_C1")
    if (
        ser_delta is None
        or abs(ser_delta) > thresholds["ser128_og_to_bu2_c1_abs_delta_max_A"]
    ):
        ligand_reasons.append("ser128_og_to_bu2_c1_delta")
    ligand_gate = "PASS" if not ligand_reasons else "FAIL"

    if inherited_gate == "PASS" and ligand_gate == "PASS":
        status = "CROP_STRICT_PASS"
    elif inherited_gate == "FAIL" and ligand_gate == "FAIL":
        status = "FAIL_BOTH_PROTEIN_AND_LIGAND_GATES"
    elif inherited_gate == "FAIL":
        status = "FAIL_INHERITED_POSTSEQ_PROTEIN_GATE"
    else:
        status = "FAIL_PLACER_LIGAND_ADDON_GATE"

    row.update(
        {
            "status": status,
            "inherited_postseq_protein_gate": inherited_gate,
            "ligand_reaction_pose_gate": ligand_gate,
            "global_backbone_atom_pairs": len(common_backbone),
            "global_backbone_rmsd_A": global_rmsd,
            "fixed_backbone_atom_pairs": len(fixed_common),
            "fixed_backbone_rmsd_A": fixed_rmsd if fixed_rmsd is not None else "",
            "catalytic_heavy_atom_pairs": len(catalytic_common),
            "catalytic_heavy_rmsd_A": catalytic_rmsd if catalytic_rmsd is not None else "",
            "max_abs_protein_key_delta_A": (
                max_abs_protein_delta if max_abs_protein_delta is not None else ""
            ),
            "ligand_heavy_rmsd_A": ligand_rmsd if ligand_rmsd is not None else "",
            "ser128_og_to_bu2_c1_delta_A": ser_delta if ser_delta is not None else "",
            "protein_fail_reasons": ";".join(protein_reasons),
            "ligand_fail_reasons": ";".join(ligand_reasons),
            "fail_reasons": ";".join(protein_reasons + ligand_reasons),
        }
    )
    return row


def sample_from_pdb_name(path: Path) -> str:
    return path.name.split(".holo_")[0]


def main() -> int:
    args = parse_args()
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    thresholds = {
        "global_backbone_rmsd_max_A": args.global_backbone_rmsd_max_a,
        "fixed_backbone_rmsd_max_A": args.fixed_backbone_rmsd_max_a,
        "catalytic_heavy_rmsd_max_A": args.catalytic_heavy_rmsd_max_a,
        "protein_key_distance_abs_delta_max_A": args.protein_key_distance_abs_delta_max_a,
        "ligand_heavy_rmsd_max_A": args.ligand_heavy_rmsd_max_a,
        "ser128_og_to_bu2_c1_abs_delta_max_A": args.ser128_og_to_bu2_c1_abs_delta_max_a,
        "minimum_crop_pass_conformers_per_sequence": args.minimum_crop_pass_conformers_per_sequence,
        "expected_conformers_per_sequence": args.expected_conformers_per_sequence,
    }

    ref_atoms = [
        parse_atom_line(line)
        for line in Path(args.reference).read_text().splitlines()
        if line.startswith(("ATOM", "HETATM"))
    ]
    ref_map = atom_map(ref_atoms)
    ref_backbone = {
        key: coord for key, coord in ref_map.items() if key[0] == "A" and key[3] in {"N", "CA", "C", "O"}
    }
    ref_fixed_backbone = {key: coord for key, coord in ref_backbone.items() if key[1] in FIXED_RESIDUES}
    ref_catalytic = residue_heavy_atoms(ref_map, CATALYTIC_RESIDUES)
    ref_protein_key_distances = {
        label: distance(ref_map[atom_a], ref_map[atom_b])
        for label, atom_a, atom_b in PROTEIN_KEY_PAIRS
        if atom_a in ref_map and atom_b in ref_map
    }
    ref_ligand_key_distances = {
        label: distance(ref_map[atom_a], ref_map[atom_b])
        for label, atom_a, atom_b in LIGAND_KEY_PAIRS
        if atom_a in ref_map and atom_b in ref_map
    }

    detail_fields = [
        "run_label",
        "sample_id",
        "bin",
        "model_id",
        "status",
        "inherited_postseq_protein_gate",
        "ligand_reaction_pose_gate",
        "global_backbone_atom_pairs",
        "global_backbone_rmsd_A",
        "fixed_backbone_atom_pairs",
        "fixed_backbone_rmsd_A",
        "catalytic_heavy_atom_pairs",
        "catalytic_heavy_rmsd_A",
        "max_abs_protein_key_delta_A",
        "ligand_heavy_rmsd_A",
        "ser128_og_to_bu2_c1_delta_A",
        "protein_fail_reasons",
        "ligand_fail_reasons",
        "fail_reasons",
        "source_pdb",
    ]
    all_rows: list[dict] = []
    sequence_rows: list[dict] = []

    for run_label, run_dir_raw in args.run_dir:
        run_dir = Path(run_dir_raw)
        run_rows: list[dict] = []
        for pdb_path in sorted(run_dir.glob("*_model.pdb")):
            sample_id = sample_from_pdb_name(pdb_path)
            bin_id = sample_id[3:5]
            for model_id, atoms in parse_models(pdb_path):
                run_rows.append(
                    evaluate_model(
                        run_label=run_label,
                        sample_id=sample_id,
                        bin_id=bin_id,
                        model_id=model_id,
                        source_pdb=pdb_path,
                        atoms=atoms,
                        ref_map=ref_map,
                        ref_backbone=ref_backbone,
                        ref_fixed_backbone=ref_fixed_backbone,
                        ref_catalytic=ref_catalytic,
                        ref_protein_key_distances=ref_protein_key_distances,
                        ref_ligand_key_distances=ref_ligand_key_distances,
                        thresholds=thresholds,
                    )
                )

        detail_path = out_dir / f"two_layer_placer_crop_gate_{run_label}.tsv"
        with detail_path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, detail_fields, delimiter="\t", extrasaction="ignore")
            writer.writeheader()
            writer.writerows(run_rows)
        all_rows.extend(run_rows)

        by_sample: dict[str, list[dict]] = {}
        for row in run_rows:
            by_sample.setdefault(row["sample_id"], []).append(row)
        for sample_id, rows in by_sample.items():
            crop_pass_count = sum(1 for row in rows if row["status"] == "CROP_STRICT_PASS")
            protein_pass_count = sum(
                1 for row in rows if row["inherited_postseq_protein_gate"] == "PASS"
            )
            ligand_pass_count = sum(1 for row in rows if row["ligand_reaction_pose_gate"] == "PASS")
            sequence_status = (
                "PLACER_SEQUENCE_ACCEPTED_FOR_FULL_LIGAND_COMPLETION"
                if len(rows) == args.expected_conformers_per_sequence
                and crop_pass_count >= args.minimum_crop_pass_conformers_per_sequence
                else "DISCARD_SEQUENCE_REGENERATE_FROM_UPSTREAM"
            )
            ligand_values = [
                float(row["ligand_heavy_rmsd_A"])
                for row in rows
                if row.get("ligand_heavy_rmsd_A") != ""
            ]
            sequence_rows.append(
                {
                    "run_label": run_label,
                    "sample_id": sample_id,
                    "bin": rows[0]["bin"],
                    "evaluated_conformers": len(rows),
                    "protein_gate_pass_conformers": protein_pass_count,
                    "ligand_addon_gate_pass_conformers": ligand_pass_count,
                    "crop_strict_pass_conformers": crop_pass_count,
                    "sequence_status": sequence_status,
                    "best_ligand_heavy_rmsd_A": min(ligand_values) if ligand_values else "",
                }
            )

    sequence_fields = [
        "run_label",
        "sample_id",
        "bin",
        "evaluated_conformers",
        "protein_gate_pass_conformers",
        "ligand_addon_gate_pass_conformers",
        "crop_strict_pass_conformers",
        "sequence_status",
        "best_ligand_heavy_rmsd_A",
    ]
    sequence_path = out_dir / "two_layer_placer_crop_gate_sequence_summary.tsv"
    with sequence_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, sequence_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sequence_rows)

    conformer_status_counts: dict[str, dict[str, int]] = {}
    for row in all_rows:
        conformer_status_counts.setdefault(row["run_label"], {}).setdefault(row["status"], 0)
        conformer_status_counts[row["run_label"]][row["status"]] += 1

    sequence_status_counts: dict[str, dict[str, int]] = {}
    for row in sequence_rows:
        sequence_status_counts.setdefault(row["run_label"], {}).setdefault(row["sequence_status"], 0)
        sequence_status_counts[row["run_label"]][row["sequence_status"]] += 1

    summary = {
        "evaluated_universe": {
            "runs": [label for label, _ in args.run_dir],
            "total_conformers": len(all_rows),
            "total_sequences": len(sequence_rows),
        },
        "thresholds": thresholds,
        "conformer_status_counts": conformer_status_counts,
        "sequence_status_counts": sequence_status_counts,
        "outputs": {"sequence_summary_tsv": str(sequence_path)},
    }
    summary_path = out_dir / "two_layer_placer_crop_gate_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True))
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
