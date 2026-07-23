#!/usr/bin/env python3
"""Audit one fully unrestrained NylC true-Thr267 NPT pilot."""
import argparse
import json
import pathlib
import re
import statistics
import sys

SCRIPT_DIR = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(SCRIPT_DIR))
import audit_true_thr267_recapture_pilot as recapture
import extract_true_thr267_recapture_sources as geom

GATE_AXIS = (0.4904295935403325, 0.813080787402625, -0.3136533866174433)
GATE_BASELINE = -1.51109


def read_xvg(path):
    rows = []
    for line in pathlib.Path(path).read_text(errors="replace").splitlines():
        if not line.strip() or line[0] in "#@":
            continue
        rows.append([float(value) for value in line.split()])
    return rows


def numerical_counts(log_text):
    return {
        "fatal": len(re.findall(r"(?i)fatal(?:\s+error)?", log_text)),
        "lincs_warning": len(re.findall(r"(?i)lincs\s+warning", log_text)),
        "settle_problem": len(re.findall(
            r"(?i)(?:cannot\s+be\s+settled|settle[^\n]*(?:warning|error|failed))",
            log_text,
        )),
        "nan": len(re.findall(r"(?i)\bnan\b", log_text)),
    }


def summarize_thermo(rows):
    labels = ("temperature_K", "pressure_bar", "volume_nm3")
    result = {}
    for column, label in enumerate(labels, 1):
        values = [row[column] for row in rows]
        result[label] = {
            "mean": statistics.fmean(values),
            "stdev": statistics.pstdev(values),
            "min": min(values),
            "max": max(values),
        }
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run-root", required=True, type=pathlib.Path)
    parser.add_argument("--seed", required=True, type=int)
    parser.add_argument("--array-job-id", required=True)
    parser.add_argument("--stem", default="free_npt100")
    parser.add_argument("--window-type", default="fully_unrestrained_NPT_100ps")
    parser.add_argument("--output-name", default="unrestrained_pilot_audit.json")
    args = parser.parse_args()
    root = args.run_root
    series = json.loads((root / "nac_energy_series.json").read_text())
    log_text = (root / f"{args.stem}.log").read_text(errors="replace")
    counts = numerical_counts(log_text)
    finished = bool(re.search(r"Finished mdrun", log_text))
    atoms, box, _ = geom.read_gro(root / f"{args.stem}.gro")
    contact = recapture.minimum_ligand_protein_distance(atoms, box)
    distance, angle = recapture.geometry(
        root / f"{args.stem}.gro", geom.CANDIDATES["nylc_C18_trueT267_recapture"]
    )
    gate_rows = read_xvg(root / "gate_core_vector.xvg")
    gate = []
    for row in gate_rows:
        absolute = sum(row[i + 1] * GATE_AXIS[i] for i in range(3))
        gate.append(absolute - GATE_BASELINE)
    numerical_pass = finished and not any(counts.values())
    contact_pass = (
        contact["distance_nm"] >= 0.10
        and contact["minimum_heavy_atom_contact"]["distance_nm"] >= 0.18
    )
    nac_present = series["nac_frame_count"] > 0
    if not numerical_pass:
        scientific_status = "NOT_EVALUATED_NUMERICAL_FAIL"
    elif not contact_pass:
        scientific_status = "FAIL_UNRESTRAINED_PILOT_SEVERE_CONTACT"
    elif nac_present:
        scientific_status = "PASS_UNRESTRAINED_PILOT_NAC_PRESENT"
    else:
        scientific_status = "FAIL_UNRESTRAINED_PILOT_NO_NAC"
    result = {
        "schema_version": 1,
        "candidate_id": "nylc_C18_trueT267_recapture",
        "array_job_id": args.array_job_id,
        "velocity_seed": args.seed,
        "window_type": args.window_type,
        "technical_status": "PASS" if numerical_pass else "FAIL",
        "scientific_status": scientific_status,
        "nac": series,
        "end_geometry": {
            "distance_nm": round(distance, 6),
            "angle_deg": round(angle, 3),
            "joint_nac": bool(distance <= 0.35 and 95.0 <= angle <= 115.0),
        },
        "minimum_ligand_protein_contact": contact,
        "gate_definition": "NylC residues 261-266; Thr267 excluded",
        "gate_opening_nm": {
            "frame_count": len(gate),
            "mean": statistics.fmean(gate),
            "stdev": statistics.pstdev(gate),
            "min": min(gate),
            "max": max(gate),
        },
        "thermodynamics": summarize_thermo(read_xvg(root / "thermo.xvg")),
        "numerical_issue_counts": counts,
        "next_gate": (
            "compare all three seeds; extend at least 1 ns only from a fully "
            "unrestrained NAC-positive seed and choose the lower-potential NAC "
            "from the equilibrated free window"
        ),
    }
    (root / args.output_name).write_text(
        json.dumps(result, indent=2, sort_keys=True) + "\n"
    )
    print(json.dumps(result, sort_keys=True))
    raise SystemExit(0 if numerical_pass else 2)


if __name__ == "__main__":
    main()
