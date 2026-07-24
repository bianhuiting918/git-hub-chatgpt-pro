#!/usr/bin/env python3
"""Audit whether the NylC-C18 oxyanion hole is supplied by either protein copy."""
import argparse
import csv
import json
import pathlib

COPY1_TYR146_N = 2020
COPY1_TYR146_H = 2021
COPY1_ASN219_ND2 = 3104
COPY1_ASN219_HD21 = 3105
COPY1_ASN219_HD22 = 3106
L2_REACTIVE_O = 10288

DONOR_HEAVY_MAX_NM = 0.35
DONOR_H_MAX_NM = 0.25
DONOR_ANGLE_MIN_DEG = 135.0


def read_xvg(path):
    rows = []
    for line in pathlib.Path(path).read_text(errors="replace").splitlines():
        if line.strip() and line[0] not in "#@":
            rows.append([float(value) for value in line.split()])
    return rows


def as_bool(value):
    return str(value).lower() in ("true", "1", "yes")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-frames", required=True)
    parser.add_argument("--copy1-distances", required=True)
    parser.add_argument("--copy1-angles", required=True)
    parser.add_argument("--output-dir", required=True)
    args = parser.parse_args()

    with pathlib.Path(args.base_frames).open() as handle:
        base = list(csv.DictReader(handle, delimiter="\t"))
    distances = read_xvg(args.copy1_distances)
    angles = read_xvg(args.copy1_angles)
    if not (len(base) == len(distances) == len(angles)):
        raise ValueError(f"series length mismatch: base={len(base)} distance={len(distances)} angle={len(angles)}")
    if not base:
        raise ValueError("no frames")

    combo_names = (
        "copy1_tyr_copy1_asn",
        "copy1_tyr_copy2_asn",
        "copy2_tyr_copy1_asn",
        "copy2_tyr_copy2_asn",
    )
    combo_counts = {name: 0 for name in combo_names}
    augmented = []
    qualified = []
    for row, drow, arow in zip(base, distances, angles):
        time_ps = float(row["time_ps"])
        if max(abs(time_ps - drow[0]), abs(time_ps - arow[0])) > 1.0e-5:
            raise ValueError(f"time mismatch at {time_ps}: {drow[0]}, {arow[0]}")
        if len(drow) != 6 or len(arow) != 5:
            raise ValueError(f"column mismatch at {time_ps}: distance={len(drow)} angle={len(arow)}")
        tyr1 = drow[1] <= DONOR_HEAVY_MAX_NM and drow[2] <= DONOR_H_MAX_NM and arow[2] >= DONOR_ANGLE_MIN_DEG
        asn1_h21 = drow[3] <= DONOR_HEAVY_MAX_NM and drow[4] <= DONOR_H_MAX_NM and arow[3] >= DONOR_ANGLE_MIN_DEG
        asn1_h22 = drow[3] <= DONOR_HEAVY_MAX_NM and drow[5] <= DONOR_H_MAX_NM and arow[4] >= DONOR_ANGLE_MIN_DEG
        asn1 = asn1_h21 or asn1_h22
        tyr2 = as_bool(row["tyr_oxyanion_ready"])
        asn2 = as_bool(row["asn_oxyanion_ready"])
        strict = as_bool(row["strict_nac"])
        combos = {
            "copy1_tyr_copy1_asn": strict and tyr1 and asn1,
            "copy1_tyr_copy2_asn": strict and tyr1 and asn2,
            "copy2_tyr_copy1_asn": strict and tyr2 and asn1,
            "copy2_tyr_copy2_asn": strict and tyr2 and asn2,
        }
        for name, passed in combos.items():
            combo_counts[name] += int(passed)
        passed_names = [name for name, passed in combos.items() if passed]
        item = dict(row)
        item.update({
            "copy1_tyr_n_o_nm": drow[1],
            "copy1_tyr_h_o_nm": drow[2],
            "copy1_tyr_n_h_o_deg": arow[2],
            "copy1_asn_nd2_o_nm": drow[3],
            "copy1_asn_hd21_o_nm": drow[4],
            "copy1_asn_hd22_o_nm": drow[5],
            "copy1_asn21_n_h_o_deg": arow[3],
            "copy1_asn22_n_h_o_deg": arow[4],
            "copy1_tyr_ready": tyr1,
            "copy1_asn_ready": asn1,
            "preorganized_donor_combinations": ",".join(passed_names),
        })
        augmented.append(item)
        if passed_names:
            qualified.append({
                "time_ps": time_ps,
                "potential_energy_kj_mol": float(row["potential_energy_kj_mol"]),
                "distance_nm": float(row["distance_nm"]),
                "attack_angle_deg": float(row["attack_angle_deg"]),
                "donor_combinations": passed_names,
            })

    selected = min(qualified, key=lambda item: item["potential_energy_kj_mol"]) if qualified else None
    strict_rows = [row for row in augmented if as_bool(row["strict_nac"])]
    min_strict = {}
    for key in (
        "copy1_tyr_n_o_nm", "copy1_tyr_h_o_nm", "copy1_asn_nd2_o_nm",
        "copy1_asn_hd21_o_nm", "copy1_asn_hd22_o_nm",
        "tyr_n_o_nm", "tyr_h_o_nm", "asn_nd2_o_nm", "asn_hd21_o_nm", "asn_hd22_o_nm",
    ):
        min_strict[key] = min(float(row[key]) for row in strict_rows) if strict_rows else None

    status = (
        "PASS_CROSS_SUBUNIT_OXYANION_PREORGANIZATION_PRESENT"
        if qualified
        else "NOT_EVALUATED_NO_CROSS_SUBUNIT_OXYANION_PREORGANIZATION"
    )
    output = pathlib.Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=False)
    fieldnames = list(augmented[0])
    with (output / "step1_chain_assignment_frames.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(augmented)
    audit = {
        "schema_version": 1,
        "candidate_id": "nylc_C18_trueT267_freeGS",
        "status": status,
        "frame_count": len(base),
        "strict_nac_frame_count": len(strict_rows),
        "preorganized_frame_count": len(qualified),
        "donor_combination_counts": combo_counts,
        "lowest_potential_preorganized_frame": selected,
        "minimum_donor_distances_within_strict_nac_nm": min_strict,
        "copy_assignment": {
            "copy1": {
                "Tyr146_N_H": [COPY1_TYR146_N, COPY1_TYR146_H],
                "Asn219_ND2_H": [COPY1_ASN219_ND2, COPY1_ASN219_HD21, COPY1_ASN219_HD22],
            },
            "copy2": "values inherited from step1_network_job_61713551",
            "reactive_oxygen": L2_REACTIVE_O,
        },
        "thresholds": {
            "donor_heavy_max_nm": DONOR_HEAVY_MAX_NM,
            "donor_h_max_nm": DONOR_H_MAX_NM,
            "donor_angle_min_deg": DONOR_ANGLE_MIN_DEG,
        },
        "interpretation": "All Tyr146/Asn219 copy combinations were tested; zero counts require bounded oxyanion recapture before TS.",
    }
    (output / "step1_chain_assignment_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    print(json.dumps(audit, sort_keys=True))


if __name__ == "__main__":
    main()
