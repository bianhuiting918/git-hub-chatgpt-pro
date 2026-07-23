#!/usr/bin/env python3
"""Audit numerical stability and final geometry of a restrained recapture pilot."""
import argparse
import json
import pathlib
import re
import sys

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import extract_true_thr267_recapture_sources as geom

TASK_ROOT = pathlib.Path(
    "/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723"
)


def geometry(path, cfg):
    atoms, box, _ = geom.read_gro(path)
    true_thr = geom.validate_identity(
        atoms, geom.TRUE_THR267_OG1, (267, "THR", "OG1")
    )
    geom.validate_identity(
        atoms, geom.SUPERSEDED_THR262_OG1, (262, "THR", "OG1")
    )
    carbon = atoms[cfg["carbonyl_c"]]
    oxygen = atoms[cfg["carbonyl_o"]]
    c_to_thr = geom.minimum_image(
        geom.subtract(true_thr["xyz_nm"], carbon["xyz_nm"]), box
    )
    c_to_o = geom.minimum_image(
        geom.subtract(oxygen["xyz_nm"], carbon["xyz_nm"]), box
    )
    return geom.norm(c_to_thr), geom.angle_deg(c_to_thr, c_to_o)


def minimum_ligand_protein_distance(atoms, box):
    ligand = [i for i, atom in atoms.items() if atom["resname"] == "UNL"]
    if not ligand:
        raise ValueError("no UNL ligand atoms in GRO")
    first_ligand = min(ligand)
    protein = [i for i in atoms if i < first_ligand]
    if not protein:
        raise ValueError("no protein atoms before UNL ligand")
    best = None
    best_heavy = None
    for ligand_index in ligand:
        ligand_atom = atoms[ligand_index]
        for protein_index in protein:
            protein_atom = atoms[protein_index]
            vector = geom.minimum_image(
                geom.subtract(ligand_atom["xyz_nm"], protein_atom["xyz_nm"]), box
            )
            distance = geom.norm(vector)
            if best is None or distance < best[0]:
                best = (distance, protein_index, ligand_index)
            if (not ligand_atom["atomname"].startswith("H")
                    and not protein_atom["atomname"].startswith("H")
                    and (best_heavy is None or distance < best_heavy[0])):
                best_heavy = (distance, protein_index, ligand_index)
    distance, protein_index, ligand_index = best
    heavy_distance, heavy_protein_index, heavy_ligand_index = best_heavy
    return {
        "distance_nm": round(distance, 6),
        "protein_global_index": protein_index,
        "protein_identity": {
            key: atoms[protein_index][key]
            for key in ("resid", "resname", "atomname")
        },
        "ligand_global_index": ligand_index,
        "ligand_identity": {
            key: atoms[ligand_index][key]
            for key in ("resid", "resname", "atomname")
        },
        "minimum_heavy_atom_contact": {
            "distance_nm": round(heavy_distance, 6),
            "protein_global_index": heavy_protein_index,
            "protein_identity": {
                key: atoms[heavy_protein_index][key]
                for key in ("resid", "resname", "atomname")
            },
            "ligand_global_index": heavy_ligand_index,
            "ligand_identity": {
                key: atoms[heavy_ligand_index][key]
                for key in ("resid", "resname", "atomname")
            },
        },
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", choices=sorted(geom.CANDIDATES), required=True)
    parser.add_argument("--stage-dir", default="recapture_pilot")
    parser.add_argument("--stem", default="pilot")
    parser.add_argument("--min-response", type=float, default=0.15)
    parser.add_argument("--max-end-distance", type=float)
    parser.add_argument("--reference-gro")
    args = parser.parse_args()
    root = TASK_ROOT / "candidates" / args.candidate
    work = root / args.stage_dir
    source_manifest = json.loads((root / "source_manifest.json").read_text())
    cfg = geom.CANDIDATES[args.candidate]
    if args.reference_gro:
        start_distance, start_angle = geometry(pathlib.Path(args.reference_gro), cfg)
    else:
        start_distance = source_manifest["geometry"]["true_thr267_og1_to_carbonyl_c_nm"]
        start_angle = source_manifest["geometry"]["true_thr267_attack_angle_deg"]
    end_gro = work / f"{args.stem}.gro"
    end_distance, end_angle = geometry(end_gro, cfg)
    end_atoms, end_box, _ = geom.read_gro(end_gro)
    contact = minimum_ligand_protein_distance(end_atoms, end_box)
    pull_rows = []
    pull_path = work / f"{args.stem}_pullx.xvg"
    if pull_path.exists():
        for line in pull_path.read_text(errors="replace").splitlines():
            if line and line[0] not in "#@":
                values = [float(value) for value in line.split()]
                if len(values) >= 3:
                    pull_rows.append(values[:3])
    nac_flags = [row[1] <= 0.35 and 95.0 <= row[2] <= 115.0 for row in pull_rows]
    longest_frames = 0
    current_frames = 0
    for passed in nac_flags:
        current_frames = current_frames + 1 if passed else 0
        longest_frames = max(longest_frames, current_frames)
    sampling_ps = (
        pull_rows[1][0] - pull_rows[0][0] if len(pull_rows) > 1 else None
    )
    trajectory_nac = {
        "frame_count": len(pull_rows),
        "nac_frame_count": sum(nac_flags),
        "occupancy": sum(nac_flags) / len(pull_rows) if pull_rows else None,
        "sampling_ps": sampling_ps,
        "longest_continuous_nac_frames": longest_frames,
        "longest_continuous_nac_ps": (
            max(0, longest_frames - 1) * sampling_ps
            if sampling_ps is not None else None
        ),
        "minimum_distance_nm": min((row[1] for row in pull_rows), default=None),
    }
    log_text = (work / f"{args.stem}.log").read_text(errors="replace")
    counts = {
        "fatal": len(re.findall(r"(?i)fatal(?:\s+error)?", log_text)),
        "lincs_warning": len(re.findall(r"(?i)lincs\s+warning", log_text)),
        "settle_problem": len(re.findall(
            r"(?i)(?:cannot\s+be\s+settled|settle[^\n]*(?:warning|error|failed))",
            log_text,
        )),
        "nan": len(re.findall(r"(?i)\bnan\b", log_text)),
    }
    finished = bool(re.search(r"Finished mdrun", log_text))
    response = start_distance - end_distance
    numerical_pass = finished and not any(counts.values())
    response_pass = (
        response >= args.min_response
        and end_distance >= 0.30
        and 75.0 <= end_angle <= 135.0
        and contact["distance_nm"] >= 0.10
        and contact["minimum_heavy_atom_contact"]["distance_nm"] >= 0.18
        and (args.max_end_distance is None or end_distance <= args.max_end_distance)
    )
    if not numerical_pass:
        status = "FAIL_NUMERICAL"
    elif not response_pass:
        status = "FAIL_APPROACH_RESPONSE_NUMERICALLY_STABLE"
    else:
        status = "PASS_TECHNICAL_BOUNDED_APPROACH"
    result = {
        "schema_version": 1,
        "candidate": args.candidate,
        "status": status,
        "scientific_gate": "NOT_EVALUATED_RESTRAINED_PILOT_CANNOT_ESTABLISH_NAC",
        "finished_mdrun": finished,
        "minimum_required_response_nm": args.min_response,
        "maximum_allowed_end_distance_nm": args.max_end_distance,
        "reference_gro": args.reference_gro,
        "numerical_issue_counts": counts,
        "minimum_ligand_protein_contact": contact,
        "trajectory_nac": trajectory_nac,
        "severe_contact_threshold_nm": 0.10,
        "severe_heavy_atom_contact_threshold_nm": 0.18,
        "geometry": {
            "start_distance_nm": start_distance,
            "end_distance_nm": round(end_distance, 6),
            "distance_response_nm": round(response, 6),
            "start_angle_deg": start_angle,
            "end_angle_deg": round(end_angle, 3),
            "end_joint_nac": bool(
                end_distance <= 0.35 and 95.0 <= end_angle <= 115.0
            ),
        },
    }
    (work / f"{args.stem}_audit.json").write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, sort_keys=True))
    raise SystemExit(0 if status == "PASS_TECHNICAL_BOUNDED_APPROACH" else 2)


if __name__ == "__main__":
    main()
