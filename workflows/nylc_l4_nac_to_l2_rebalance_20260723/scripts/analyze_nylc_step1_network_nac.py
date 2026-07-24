#!/usr/bin/env python3
"""Audit complete catalytic preorganization within free NylC-C18 L2 NAC frames."""
import argparse
import hashlib
import json
import math
import os
import pathlib

TYR146_N = 7156
TYR146_H = 7157
ASN219_ND2 = 8240
ASN219_HD21 = 8241
ASN219_HD22 = 8242
THR267_N = 8949
THR267_OG1 = 8961
ASP306_OD1 = 9572
ASP306_OD2 = 9573
ASP308_OD1 = 9591
ASP308_OD2 = 9592
L2_REACTIVE_C = 10287
L2_REACTIVE_O = 10288
L2_REACTIVE_N = 10289

donor_heavy_max_nm = 0.35
donor_h_max_nm = 0.25
donor_angle_min_deg = 135.0


def read_xvg(path):
    rows = []
    for line in pathlib.Path(path).read_text(errors="replace").splitlines():
        if line.strip() and line[0] not in "#@":
            rows.append([float(value) for value in line.split()])
    return rows


def sha256(path):
    digest = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_gro(path):
    lines = pathlib.Path(path).read_text().splitlines()
    count = int(lines[1])
    atoms = {}
    for global_index, line in enumerate(lines[2 : 2 + count], 1):
        atoms[global_index] = {
            "resid": int(line[:5]),
            "resname": line[5:10].strip(),
            "atomname": line[10:15].strip(),
            "xyz": tuple(float(line[start:end]) for start, end in ((20, 28), (28, 36), (36, 44))),
        }
    box = tuple(float(value) for value in lines[2 + count].split()[:3])
    return atoms, box


def vector(atoms, box, start, end):
    result = []
    for a, b, length in zip(atoms[start]["xyz"], atoms[end]["xyz"], box):
        delta = b - a
        delta -= round(delta / length) * length
        result.append(delta)
    return tuple(result)


def norm(value):
    return math.sqrt(sum(component * component for component in value))


def angle(atoms, box, left, vertex, right):
    first = vector(atoms, box, vertex, left)
    second = vector(atoms, box, vertex, right)
    cosine = sum(a * b for a, b in zip(first, second)) / (norm(first) * norm(second))
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def validate_gro(gro, audit_path):
    audit = json.loads(pathlib.Path(audit_path).read_text())
    selected = audit.get("lowest_potential_full_preorganization_frame")
    if not selected:
        raise ValueError("audit has no full-preorganization frame to validate")
    atoms, box = read_gro(gro)
    distance_nm = norm(vector(atoms, box, THR267_OG1, L2_REACTIVE_C))
    c_to_thr = vector(atoms, box, L2_REACTIVE_C, THR267_OG1)
    c_to_o = vector(atoms, box, L2_REACTIVE_C, L2_REACTIVE_O)
    attack_angle_deg = math.degrees(math.acos(max(-1.0, min(1.0, sum(a*b for a,b in zip(c_to_thr,c_to_o))/(norm(c_to_thr)*norm(c_to_o))))))
    tyr = (
        norm(vector(atoms, box, TYR146_N, L2_REACTIVE_O)) <= donor_heavy_max_nm
        and norm(vector(atoms, box, TYR146_H, L2_REACTIVE_O)) <= donor_h_max_nm
        and angle(atoms, box, TYR146_N, TYR146_H, L2_REACTIVE_O) >= donor_angle_min_deg
    )
    asn = (
        norm(vector(atoms, box, ASN219_ND2, L2_REACTIVE_O)) <= donor_heavy_max_nm
        and (
            norm(vector(atoms, box, ASN219_HD21, L2_REACTIVE_O)) <= donor_h_max_nm
            and angle(atoms, box, ASN219_ND2, ASN219_HD21, L2_REACTIVE_O) >= donor_angle_min_deg
            or norm(vector(atoms, box, ASN219_HD22, L2_REACTIVE_O)) <= donor_h_max_nm
            and angle(atoms, box, ASN219_ND2, ASN219_HD22, L2_REACTIVE_O) >= donor_angle_min_deg
        )
    )
    passed = distance_nm <= 0.35 and 95.0 <= attack_angle_deg <= 115.0 and tyr and asn
    if not passed:
        raise ValueError("source.tmp.gro failed full catalytic-preorganization validation")
    promoted = pathlib.Path(audit_path).parent / "selected_step1_network_gs.gro"
    audit["selected_gro_validation"] = {
        "status": "PASS",
        "distance_nm": distance_nm,
        "attack_angle_deg": attack_angle_deg,
        "sha256": sha256(gro),
    }
    os.replace(gro, promoted)
    pathlib.Path(audit_path).write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    return audit


def analyze(args):
    distance = read_xvg(args.distance)
    attack = read_xvg(args.attack_angle)
    potential = read_xvg(args.potential)
    network = read_xvg(args.network_distances)
    angles = read_xvg(args.network_angles)
    sizes = {len(item) for item in (distance, attack, potential, network, angles)}
    if len(sizes) != 1:
        raise ValueError(f"series length mismatch: {sizes}")
    frames = []
    for drow, arow, prow, nrow, grow in zip(distance, attack, potential, network, angles):
        times = [drow[0], arow[0], prow[0], nrow[0], grow[0]]
        if max(times) - min(times) > 1.0e-5:
            raise ValueError(f"time mismatch: {times}")
        distance_nm = drow[1]
        attack_angle_deg = arow[1]
        strict_nac = distance_nm <= 0.35 and 95.0 <= attack_angle_deg <= 115.0
        tyr_oxyanion_ready = nrow[1] <= donor_heavy_max_nm and nrow[2] <= donor_h_max_nm and grow[2] >= donor_angle_min_deg
        asn21 = nrow[3] <= donor_heavy_max_nm and nrow[4] <= donor_h_max_nm and grow[3] >= donor_angle_min_deg
        asn22 = nrow[3] <= donor_heavy_max_nm and nrow[5] <= donor_h_max_nm and grow[4] >= donor_angle_min_deg
        asn_oxyanion_ready = asn21 or asn22
        full = strict_nac and tyr_oxyanion_ready and asn_oxyanion_ready
        frames.append({
            "time_ps": drow[0], "distance_nm": distance_nm, "attack_angle_deg": attack_angle_deg,
            "potential_energy_kj_mol": prow[1], "strict_nac": strict_nac,
            "tyr_oxyanion_ready": tyr_oxyanion_ready, "asn_oxyanion_ready": asn_oxyanion_ready,
            "full_preorganization": full, "tyr_n_o_nm": nrow[1], "tyr_h_o_nm": nrow[2],
            "tyr_n_h_o_deg": grow[2], "asn_nd2_o_nm": nrow[3], "asn_hd21_o_nm": nrow[4],
            "asn_hd22_o_nm": nrow[5], "asn21_n_h_o_deg": grow[3], "asn22_n_h_o_deg": grow[4],
            "asp306_thr_n_min_nm": min(nrow[6:8]), "asp306_asp308_min_nm": min(nrow[8:12]),
            "lys189_thr_og1_nm": nrow[12], "asn219_od1_lys189_nz_nm": nrow[13],
        })
    qualified = [frame for frame in frames if frame["full_preorganization"]]
    selected = min(qualified, key=lambda frame: frame["potential_energy_kj_mol"]) if qualified else None
    status = "PASS_FULL_CATALYTIC_PREORGANIZATION_PRESENT" if selected else "NOT_EVALUATED_NO_FULL_CATALYTIC_PREORGANIZATION"
    output = pathlib.Path(args.output_dir)
    output.mkdir(parents=True, exist_ok=False)
    columns = list(frames[0])
    with (output / "step1_network_frames.tsv").open("w") as handle:
        handle.write("\t".join(columns) + "\n")
        for frame in frames:
            handle.write("\t".join(str(frame[key]) for key in columns) + "\n")
    audit = {
        "schema_version": 1, "candidate_id": "nylc_C18_trueT267_freeGS", "status": status,
        "frame_count": len(frames), "strict_nac_frame_count": sum(f["strict_nac"] for f in frames),
        "full_preorganization_frame_count": len(qualified),
        "lowest_potential_full_preorganization_frame": selected,
        "thresholds": {"distance_max_nm": 0.35, "attack_angle_deg": [95.0, 115.0],
            "donor_heavy_max_nm": donor_heavy_max_nm, "donor_h_max_nm": donor_h_max_nm,
            "donor_angle_min_deg": donor_angle_min_deg},
        "reactive_mapping": "NylC-C18 -> L2 N3-C12(O2)",
        "aspartate_interpretation": "Asp distances are diagnostic only until the Asp306 protonation microstate is audited.",
    }
    (output / "step1_network_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    print(json.dumps(audit, sort_keys=True))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--distance")
    parser.add_argument("--attack-angle")
    parser.add_argument("--potential")
    parser.add_argument("--network-distances")
    parser.add_argument("--network-angles")
    parser.add_argument("--output-dir")
    parser.add_argument("--validate-gro")
    parser.add_argument("--audit")
    args = parser.parse_args()
    if args.validate_gro:
        print(json.dumps(validate_gro(pathlib.Path(args.validate_gro), pathlib.Path(args.audit)), sort_keys=True))
    else:
        analyze(args)


if __name__ == "__main__":
    main()
