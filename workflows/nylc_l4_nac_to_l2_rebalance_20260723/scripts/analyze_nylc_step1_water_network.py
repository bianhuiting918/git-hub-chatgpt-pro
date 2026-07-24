#!/usr/bin/env python3
"""Audit the paper-consistent catalytic water network in the 434 ps NylC GS hypothesis."""
import argparse
import hashlib
import json
import math
import pathlib

THR267_OG1 = 8961
LYS189_NZ = 7768
ASN219_OD1 = 8239
ASN219_ND2 = 8240
ASP306_OD1 = 9572
ASP306_OD2 = 9573
ASP308_OD1 = 9591
ASP308_OD2 = 9592
L2_C12 = 10287
L2_O2 = 10288
PROTEIN_LAST = 10272
L2_FIRST = 10273
L2_LAST = 10351
WATER_EDGE_MAX_NM = 0.40
EXPECTED_SOURCE_SHA256 = "182f68a41c72c490e0046323c9b8a0f8514ccf7af7dc4d578dcb87342c1c8e2d"


def sha256(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


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
    box = tuple(float(value) for value in lines[count + 2].split()[:3])
    return atoms, box


def vector(atoms, box, start, end):
    out = []
    for left, right, length in zip(atoms[start]["xyz"], atoms[end]["xyz"], box):
        delta = right - left
        delta -= round(delta / length) * length
        out.append(delta)
    return tuple(out)


def norm(value):
    return math.sqrt(sum(component * component for component in value))


def distance(atoms, box, left, right):
    return norm(vector(atoms, box, left, right))


def angle(atoms, box, left, vertex, right):
    first = vector(atoms, box, vertex, left)
    second = vector(atoms, box, vertex, right)
    cosine = sum(a * b for a, b in zip(first, second)) / (norm(first) * norm(second))
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def min_distance(atoms, box, left, right_group):
    return min(distance(atoms, box, left, atom) for atom in right_group)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gro", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()
    gro = pathlib.Path(args.gro)
    observed_sha = sha256(gro)
    if observed_sha != EXPECTED_SOURCE_SHA256:
        raise ValueError(f"source SHA mismatch: {observed_sha}")
    atoms, box = read_gro(gro)
    identities = {
        THR267_OG1: ("THR", "OG1"), LYS189_NZ: ("LYS", "NZ"),
        ASN219_OD1: ("ASN", "OD1"), ASN219_ND2: ("ASN", "ND2"),
        ASP306_OD1: ("ASP", "OD1"), ASP306_OD2: ("ASP", "OD2"),
        ASP308_OD1: ("ASP", "OD1"), ASP308_OD2: ("ASP", "OD2"),
        L2_C12: ("L2", "C12"), L2_O2: ("L2", "O2"),
    }
    for index, expected in identities.items():
        observed = (atoms[index]["resname"], atoms[index]["atomname"])
        if observed != expected:
            raise ValueError(f"identity mismatch at {index}: {observed} != {expected}")

    groups = {
        "Thr267_OG1": [THR267_OG1],
        "Lys189_NZ": [LYS189_NZ],
        "Asn219_OD1_ND2": [ASN219_OD1, ASN219_ND2],
        "Asp306_OD": [ASP306_OD1, ASP306_OD2],
        "Asp308_OD": [ASP308_OD1, ASP308_OD2],
    }
    waters = []
    for index, atom in atoms.items():
        if atom["resname"] not in ("SOL", "HOH", "WAT") or atom["atomname"] not in ("OW", "O"):
            continue
        distances = {name: min_distance(atoms, box, index, members) for name, members in groups.items()}
        if min(distances.values()) <= WATER_EDGE_MAX_NM:
            waters.append({"oxygen_atom": index, "water_resid": atom["resid"], "distances_nm": distances})

    thr267_asn219_water_bridges = [
        water for water in waters
        if water["distances_nm"]["Thr267_OG1"] <= WATER_EDGE_MAX_NM
        and water["distances_nm"]["Asn219_OD1_ND2"] <= WATER_EDGE_MAX_NM
    ]
    asp306_asp308_water_bridges = [
        water for water in waters
        if water["distances_nm"]["Asp306_OD"] <= WATER_EDGE_MAX_NM
        and water["distances_nm"]["Asp308_OD"] <= WATER_EDGE_MAX_NM
    ]
    lys189_asn219_direct_min_nm = min(
        distance(atoms, box, LYS189_NZ, ASN219_OD1),
        distance(atoms, box, LYS189_NZ, ASN219_ND2),
    )
    attack_distance = distance(atoms, box, THR267_OG1, L2_C12)
    attack_angle = angle(atoms, box, THR267_OG1, L2_C12, L2_O2)
    strict_nac = attack_distance <= 0.35 and 95.0 <= attack_angle <= 115.0
    heavy_min = min(
        distance(atoms, box, protein, ligand)
        for protein in range(1, PROTEIN_LAST + 1)
        if not atoms[protein]["atomname"].startswith("H")
        for ligand in range(L2_FIRST, L2_LAST + 1)
        if not atoms[ligand]["atomname"].startswith("H")
    )
    passed = (
        strict_nac
        and heavy_min >= 0.18
        and bool(thr267_asn219_water_bridges)
        and bool(asp306_asp308_water_bridges)
        and lys189_asn219_direct_min_nm <= 0.35
    )
    result = {
        "schema_version": 1,
        "candidate_id": "nylc_C18_trueT267_freeGS",
        "status": "PASS_STRUCTURAL_GS_HYPOTHESIS" if passed else "NOT_EVALUATED_INCOMPLETE_CATALYTIC_NETWORK",
        "source_time_ps": 434,
        "source_sha256": observed_sha,
        "potential_energy_kj_mol": -1832730.5,
        "strict_nac": {
            "passed": strict_nac,
            "thr267_og1_l2_c12_nm": attack_distance,
            "attack_angle_deg": attack_angle,
        },
        "minimum_heavy_ligand_protein_distance_nm": heavy_min,
        "water_edge_max_nm": WATER_EDGE_MAX_NM,
        "thr267_asn219_water_bridges": thr267_asn219_water_bridges,
        "asp306_asp308_water_bridges": asp306_asp308_water_bridges,
        "lys189_asn219_direct_min_nm": lys189_asn219_direct_min_nm,
        "direct_donor_gate": "SUPERSEDED_INVALID_DIRECT_DONOR_ASSUMPTION",
        "interpretation": "Paper-consistent structural GS hypothesis; water graph is not proof of a proton-transfer path or a transition state.",
        "next": "Audit Asp306/Asp308 protonation microstates and production QM-region boundaries before Step1 RC scans.",
    }
    pathlib.Path(args.output).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))


if __name__ == "__main__":
    main()
