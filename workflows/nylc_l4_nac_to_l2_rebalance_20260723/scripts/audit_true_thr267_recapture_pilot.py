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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--candidate", choices=sorted(geom.CANDIDATES), required=True)
    parser.add_argument("--stage-dir", default="recapture_pilot")
    parser.add_argument("--stem", default="pilot")
    args = parser.parse_args()
    root = TASK_ROOT / "candidates" / args.candidate
    work = root / args.stage_dir
    source_manifest = json.loads((root / "source_manifest.json").read_text())
    cfg = geom.CANDIDATES[args.candidate]
    start_distance = source_manifest["geometry"]["true_thr267_og1_to_carbonyl_c_nm"]
    start_angle = source_manifest["geometry"]["true_thr267_attack_angle_deg"]
    end_distance, end_angle = geometry(work / f"{args.stem}.gro", cfg)
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
        response >= 0.15
        and end_distance >= 0.30
        and 75.0 <= end_angle <= 135.0
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
        "numerical_issue_counts": counts,
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
