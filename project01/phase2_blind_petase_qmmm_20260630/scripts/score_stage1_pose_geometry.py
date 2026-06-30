#!/usr/bin/env python3
"""Score blind Stage 1 PETase Michaelis-complex geometry.

Inputs are generated complex coordinates and ligand atom-label tables. The
filters are broad generic serine-hydrolase plausibility checks, not paper
reaction coordinates or paper transition-state criteria.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


DONOR_RESIDUES = {
    "ARG",
    "ASN",
    "GLN",
    "HID",
    "HIE",
    "HIP",
    "HIS",
    "LYS",
    "SER",
    "THR",
    "TRP",
    "TYR",
}


def parse_residue_id(text: str) -> tuple[str, int, str]:
    parts = text.split(":")
    if len(parts) == 2:
        chain, resseq = parts
        atom_name = ""
    elif len(parts) == 3:
        chain, resseq, atom_name = parts
    else:
        raise ValueError(f"Residue id must be CHAIN:RESSEQ or CHAIN:RESSEQ:ATOM, got {text!r}")
    return chain, int(resseq), atom_name


def parse_pdb(path: Path) -> list[dict[str, object]]:
    atoms: list[dict[str, object]] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            atom_name = line[12:16].strip()
            resname = line[17:20].strip()
            chain = line[21].strip()
            try:
                resseq = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except ValueError:
                continue
            element = line[76:78].strip() or atom_name[0]
            atoms.append(
                {
                    "record": line[:6].strip(),
                    "atom_name": atom_name,
                    "resname": resname,
                    "chain": chain,
                    "resseq": resseq,
                    "element": element.upper(),
                    "xyz": (x, y, z),
                }
            )
    return atoms


def read_label_candidate(path: Path, model_id: str, candidate_id: str) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("model_id") == model_id and row.get("candidate_id") == candidate_id:
                return row
    raise ValueError(f"No label row for model_id={model_id!r}, candidate_id={candidate_id!r} in {path}")


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def angle(a: tuple[float, float, float], b: tuple[float, float, float], c: tuple[float, float, float]) -> float:
    # angle ABC in degrees
    ba = [a[i] - b[i] for i in range(3)]
    bc = [c[i] - b[i] for i in range(3)]
    dot = sum(ba[i] * bc[i] for i in range(3))
    nba = math.sqrt(sum(x * x for x in ba))
    nbc = math.sqrt(sum(x * x for x in bc))
    if nba == 0 or nbc == 0:
        return float("nan")
    cos_value = max(-1.0, min(1.0, dot / (nba * nbc)))
    return math.degrees(math.acos(cos_value))


def dihedral(p0, p1, p2, p3) -> float:
    # Minimal vector implementation using scalar products.
    def sub(a, b):
        return [a[i] - b[i] for i in range(3)]

    def cross(a, b):
        return [
            a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0],
        ]

    def dot(a, b):
        return sum(a[i] * b[i] for i in range(3))

    def norm(v):
        n = math.sqrt(dot(v, v))
        return [x / n for x in v] if n else [float("nan")] * 3

    b0 = sub(p0, p1)
    b1 = sub(p2, p1)
    b2 = sub(p3, p2)
    b1n = norm(b1)
    v = sub(b0, [dot(b0, b1n) * x for x in b1n])
    w = sub(b2, [dot(b2, b1n) * x for x in b1n])
    x = dot(v, w)
    y = dot(cross(b1n, v), w)
    return math.degrees(math.atan2(y, x))


def find_residue_atom(atoms, chain: str, resseq: int, atom_name: str):
    matches = [a for a in atoms if a["chain"] == chain and a["resseq"] == resseq and a["atom_name"] == atom_name]
    if not matches:
        raise ValueError(f"Missing atom {chain}:{resseq}:{atom_name}")
    return matches[0]


def find_ligand_atom(atoms, atom_name: str):
    matches = [a for a in atoms if a["record"] == "HETATM" and a["atom_name"] == atom_name]
    if not matches:
        raise ValueError(f"Missing ligand HETATM named {atom_name}")
    if len(matches) > 1:
        raise ValueError(f"Ambiguous ligand atom name {atom_name}: {len(matches)} matches")
    return matches[0]


def donor_atoms(atoms):
    donors = []
    for atom in atoms:
        element = atom["element"]
        if element not in {"N", "O"}:
            continue
        atom_name = atom["atom_name"]
        resname = atom["resname"]
        if atom_name == "N" or resname in DONOR_RESIDUES:
            donors.append(atom)
    return donors


def trp_chi1_label(atoms, chain: str, resseq: int) -> str:
    try:
        n = find_residue_atom(atoms, chain, resseq, "N")["xyz"]
        ca = find_residue_atom(atoms, chain, resseq, "CA")["xyz"]
        cb = find_residue_atom(atoms, chain, resseq, "CB")["xyz"]
        cg = find_residue_atom(atoms, chain, resseq, "CG")["xyz"]
    except ValueError:
        return "not_available"
    chi = dihedral(n, ca, cb, cg)
    if -90 <= chi <= 0:
        label = "gauche_minus"
    elif 0 < chi <= 90:
        label = "gauche_plus"
    else:
        label = "trans"
    return f"{label}:{chi:.1f}"


def score_pose(args) -> dict[str, object]:
    atoms = parse_pdb(Path(args.complex_pdb))
    label = read_label_candidate(Path(args.label_tsv), args.model_id, args.candidate_id)
    ser_chain, ser_resseq, ser_atom_name = parse_residue_id(args.ser_og)
    his_chain, his_resseq, _ = parse_residue_id(args.his)
    trp_chain, trp_resseq, _ = parse_residue_id(args.trp)

    ser_og = find_residue_atom(atoms, ser_chain, ser_resseq, ser_atom_name or "OG")
    his_ns = [
        find_residue_atom(atoms, his_chain, his_resseq, name)
        for name in ("ND1", "NE2")
        if any(a["chain"] == his_chain and a["resseq"] == his_resseq and a["atom_name"] == name for a in atoms)
    ]
    if not his_ns:
        raise ValueError(f"Missing His ring nitrogens for {args.his}")

    c_atom = find_ligand_atom(atoms, label["scissile_carbonyl_c_atom_name"])
    o_atom = find_ligand_atom(atoms, label["scissile_carbonyl_o_atom_name"])
    leaving_o = find_ligand_atom(atoms, label["leaving_o_atom_name"])

    ser_to_c = distance(ser_og["xyz"], c_atom["xyz"])
    attack = angle(ser_og["xyz"], c_atom["xyz"], o_atom["xyz"])
    his_acceptor = min(distance(ser_og["xyz"], h["xyz"]) for h in his_ns)
    leaving_relay = min(distance(leaving_o["xyz"], h["xyz"]) for h in his_ns)
    oxyanion_count = sum(
        1
        for donor in donor_atoms(atoms)
        if donor is not o_atom and distance(o_atom["xyz"], donor["xyz"]) <= args.oxyanion_cutoff
    )
    trp_label = trp_chi1_label(atoms, trp_chain, trp_resseq)

    failed = []
    if ser_to_c > args.max_ser_c:
        failed.append("nucleophile_distance")
    if not (args.min_attack_angle <= attack <= args.max_attack_angle):
        failed.append("attack_angle")
    if oxyanion_count < 1:
        failed.append("oxyanion_hbond")
    if his_acceptor > args.max_his_acceptor:
        failed.append("his_accessibility")
    if leaving_relay > args.max_leaving_relay:
        failed.append("leaving_group_relay")

    return {
        "pose_id": args.pose_id,
        "template_pdb": args.template_pdb,
        "substrate_model": args.model_id,
        "generation_method": args.generation_method,
        "relaxed_structure_path": args.complex_pdb,
        "ser_og_to_c_A": f"{ser_to_c:.3f}",
        "attack_angle_deg": f"{attack:.2f}",
        "oxyanion_hbond_count": oxyanion_count,
        "his_acceptor_distance_A": f"{his_acceptor:.3f}",
        "leaving_group_relay_distance_A": f"{leaving_relay:.3f}",
        "trp_rotamer_label": trp_label,
        "pass_fail": "pass" if not failed else "fail",
        "rejection_reason": "none" if not failed else ",".join(failed),
        "next_step": "classical_md_candidate" if not failed else "record_in_rejected_pose_manifest",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--complex-pdb", required=True)
    parser.add_argument("--label-tsv", required=True)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--candidate-id", required=True)
    parser.add_argument("--pose-id", required=True)
    parser.add_argument("--template-pdb", required=True)
    parser.add_argument("--generation-method", default="pending")
    parser.add_argument("--ser-og", default="A:160:OG")
    parser.add_argument("--his", default="A:237")
    parser.add_argument("--trp", default="A:185")
    parser.add_argument("--max-ser-c", type=float, default=3.6)
    parser.add_argument("--min-attack-angle", type=float, default=80.0)
    parser.add_argument("--max-attack-angle", type=float, default=125.0)
    parser.add_argument("--oxyanion-cutoff", type=float, default=3.2)
    parser.add_argument("--max-his-acceptor", type=float, default=3.4)
    parser.add_argument("--max-leaving-relay", type=float, default=4.0)
    parser.add_argument("--out-tsv", default="-")
    args = parser.parse_args()

    row = score_pose(args)
    fieldnames = [
        "pose_id",
        "template_pdb",
        "substrate_model",
        "generation_method",
        "relaxed_structure_path",
        "ser_og_to_c_A",
        "attack_angle_deg",
        "oxyanion_hbond_count",
        "his_acceptor_distance_A",
        "leaving_group_relay_distance_A",
        "trp_rotamer_label",
        "pass_fail",
        "rejection_reason",
        "next_step",
    ]
    if args.out_tsv == "-":
        import sys

        writer = csv.DictWriter(sys.stdout, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)
    else:
        out = Path(args.out_tsv)
        out.parent.mkdir(parents=True, exist_ok=True)
        write_header = not out.exists()
        with out.open("a", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
            if write_header:
                writer.writeheader()
            writer.writerow(row)
    return 0 if row["pass_fail"] == "pass" else 1


if __name__ == "__main__":
    raise SystemExit(main())
