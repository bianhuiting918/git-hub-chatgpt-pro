#!/usr/bin/env python3
"""Audit a bounded NylC-C18 L2 oxyanion-recapture response pilot."""
import argparse
import hashlib
import json
import math
import pathlib
import re

TYR146_N = 7156
TYR146_H = 7157
ASN219_ND2 = 8240
ASN219_HD21 = 8241
THR267_OG1 = 8961
L2_C12 = 10287
L2_O2 = 10288
PROTEIN_LAST = 10272
L2_FIRST = 10273
L2_LAST = 10351


def read_gro(path):
    lines = pathlib.Path(path).read_text().splitlines()
    count = int(lines[1])
    if len(lines) != count + 3:
        raise ValueError("GRO atom-count mismatch")
    atoms = {}
    for index, line in enumerate(lines[2:2 + count], 1):
        atoms[index] = {
            "resid": int(line[:5]),
            "resname": line[5:10].strip(),
            "atomname": line[10:15].strip(),
            "xyz": tuple(float(line[a:b]) for a, b in ((20, 28), (28, 36), (36, 44))),
        }
    box = tuple(float(x) for x in lines[count + 2].split()[:3])
    return atoms, box


def vector(atoms, box, start, end):
    result = []
    for left, right, length in zip(atoms[start]["xyz"], atoms[end]["xyz"], box):
        delta = right - left
        delta -= round(delta / length) * length
        result.append(delta)
    return tuple(result)


def norm(v):
    return math.sqrt(sum(x * x for x in v))


def distance(atoms, box, a, b):
    return norm(vector(atoms, box, a, b))


def angle(atoms, box, left, vertex, right):
    a = vector(atoms, box, vertex, left)
    b = vector(atoms, box, vertex, right)
    cosine = sum(x * y for x, y in zip(a, b)) / (norm(a) * norm(b))
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def metrics(path):
    atoms, box = read_gro(path)
    expected = {
        TYR146_N: ("TYR", "N"), TYR146_H: ("TYR", "H"),
        ASN219_ND2: ("ASN", "ND2"), ASN219_HD21: ("ASN", "HD21"),
        THR267_OG1: ("THR", "OG1"), L2_C12: ("L2", "C12"), L2_O2: ("L2", "O2"),
    }
    for index, identity in expected.items():
        observed = (atoms[index]["resname"], atoms[index]["atomname"])
        if observed != identity:
            raise ValueError(f"identity mismatch {index}: {observed} != {identity}")
    heavy_min = min(
        distance(atoms, box, protein, ligand)
        for protein in range(1, PROTEIN_LAST + 1)
        if not atoms[protein]["atomname"].startswith("H")
        for ligand in range(L2_FIRST, L2_LAST + 1)
        if not atoms[ligand]["atomname"].startswith("H")
    )
    attack_distance = distance(atoms, box, THR267_OG1, L2_C12)
    attack_angle = angle(atoms, box, THR267_OG1, L2_C12, L2_O2)
    return {
        "thr267_og1_l2_c12_nm": attack_distance,
        "attack_angle_deg": attack_angle,
        "strict_nac": attack_distance <= 0.35 and 95.0 <= attack_angle <= 115.0,
        "tyr_n_o_nm": distance(atoms, box, TYR146_N, L2_O2),
        "tyr_h_o_nm": distance(atoms, box, TYR146_H, L2_O2),
        "tyr_n_h_o_deg": angle(atoms, box, TYR146_N, TYR146_H, L2_O2),
        "asn_nd2_o_nm": distance(atoms, box, ASN219_ND2, L2_O2),
        "asn_hd21_o_nm": distance(atoms, box, ASN219_HD21, L2_O2),
        "asn_n_h_o_deg": angle(atoms, box, ASN219_ND2, ASN219_HD21, L2_O2),
        "minimum_heavy_ligand_protein_distance_nm": heavy_min,
    }


def sha256(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def numerical(log_path):
    text = pathlib.Path(log_path).read_text(errors="replace")
    counts = {
        "fatal": len(re.findall(r"(?i)fatal(?:\s+error)?", text)),
        "lincs_warning": len(re.findall(r"(?i)lincs\s+warning", text)),
        "settle_problem": len(re.findall(r"(?i)(?:cannot\s+be\s+settled|settle[^\n]*(?:warning|error|failed))", text)),
        "nan": len(re.findall(r"(?i)\bnan\b", text)),
    }
    return counts, bool(re.search(r"Finished mdrun", text))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--start-gro", required=True)
    parser.add_argument("--end-gro")
    parser.add_argument("--log")
    parser.add_argument("--output", required=True)
    parser.add_argument("--source-only", action="store_true")
    args = parser.parse_args()
    start = metrics(args.start_gro)
    if args.source_only:
        passed = start["strict_nac"]
        result = {
            "status": "PASS_SOURCE_STRICT_NAC" if passed else "FAIL_SOURCE_NOT_STRICT_NAC",
            "source_sha256": sha256(args.start_gro),
            "source_geometry": start,
        }
    else:
        if not args.end_gro or not args.log:
            raise ValueError("--end-gro and --log are required")
        end = metrics(args.end_gro)
        counts, finished = numerical(args.log)
        numerical_ok = finished and not any(counts.values())
        response = (
            end["tyr_h_o_nm"] <= start["tyr_h_o_nm"] - 0.01
            and end["asn_hd21_o_nm"] <= start["asn_hd21_o_nm"] - 0.01
        )
        no_heavy_clash = end["minimum_heavy_ligand_protein_distance_nm"] >= 0.18
        passed = numerical_ok and response and no_heavy_clash
        result = {
            "status": "PASS_TECHNICAL_RESPONSE" if passed else "FAIL_TECHNICAL_OR_RESPONSE",
            "scientific_gate": "NOT_EVALUATED_RESTRAINED_OXYANION_RECAPTURE_PILOT",
            "source_sha256": sha256(args.start_gro),
            "end_sha256": sha256(args.end_gro),
            "source_geometry": start,
            "end_geometry": end,
            "response_nm": {
                "tyr_h_o_decrease": start["tyr_h_o_nm"] - end["tyr_h_o_nm"],
                "asn_hd21_o_decrease": start["asn_hd21_o_nm"] - end["asn_hd21_o_nm"],
            },
            "numerical_issue_counts": counts,
            "finished_mdrun": finished,
            "technical_numerical_ok": numerical_ok,
            "no_heavy_clash": no_heavy_clash,
            "restrained": True,
        }
    pathlib.Path(args.output).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
