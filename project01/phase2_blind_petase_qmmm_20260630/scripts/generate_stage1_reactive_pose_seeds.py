#!/usr/bin/env python3
"""Generate blind reactive-geometry seed complexes for PETase Stage 1.

The output is a starting geometry for local relaxation/MD screening, not a
transition-state or paper-derived pose.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path


FIELDS = [
    "pose_id",
    "template_pdb",
    "substrate_model",
    "candidate_id",
    "complex_pdb",
    "ser_og_to_c_target_A",
    "attack_angle_target_deg",
    "oxyanion_target_distance_A",
    "leaving_his_target_distance_A",
    "status",
    "source",
    "note",
]


def dot(a, b):
    return sum(a[i] * b[i] for i in range(3))


def cross(a, b):
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def sub(a, b):
    return (a[0] - b[0], a[1] - b[1], a[2] - b[2])


def add(a, b):
    return (a[0] + b[0], a[1] + b[1], a[2] + b[2])


def scale(a, value):
    return (a[0] * value, a[1] * value, a[2] * value)


def norm(a):
    return math.sqrt(dot(a, a))


def unit(a):
    length = norm(a)
    if length == 0:
        raise ValueError("zero-length vector")
    return scale(a, 1.0 / length)


def distance(a, b):
    return norm(sub(a, b))


def angle(a, b, c):
    ba = sub(a, b)
    bc = sub(c, b)
    return math.degrees(math.acos(max(-1.0, min(1.0, dot(unit(ba), unit(bc))))))


def parse_triad_residue(text: str) -> tuple[str, str, int]:
    chain, residue = text.split(":", 1)
    resname = "".join(ch for ch in residue if ch.isalpha())
    resseq = int("".join(ch for ch in residue if ch.isdigit()))
    return chain, resname, resseq


def parse_pdb(path: Path) -> list[dict[str, object]]:
    atoms = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            try:
                atoms.append(
                    {
                        "line": line.rstrip("\n"),
                        "record": line[:6],
                        "atom_name": line[12:16].strip(),
                        "resname": line[17:20].strip(),
                        "chain": line[21].strip(),
                        "resseq": int(line[22:26]),
                        "element": (line[76:78].strip() or line[12:16].strip()[0]).upper(),
                        "xyz": (float(line[30:38]), float(line[38:46]), float(line[46:54])),
                    }
                )
            except ValueError:
                continue
    return atoms


def format_atom_line(atom: dict[str, object], xyz, serial: int | None = None, record: str | None = None) -> str:
    line = str(atom["line"])
    prefix = record if record is not None else line[:6]
    serial_text = f"{serial:5d}" if serial is not None else line[6:11]
    return f"{prefix}{serial_text}{line[11:30]}{xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{line[54:]}"


def find_atom(atoms, chain: str, resseq: int, atom_name: str):
    for atom in atoms:
        if atom["chain"] == chain and atom["resseq"] == resseq and atom["atom_name"] == atom_name:
            return atom
    raise ValueError(f"missing atom {chain}:{resseq}:{atom_name}")


def find_ligand_atom(atoms, atom_name: str):
    for atom in atoms:
        if atom["atom_name"] == atom_name:
            return atom
    raise ValueError(f"missing ligand atom {atom_name}")


def read_label(path: Path, model_id: str, candidate_id: str) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row.get("model_id") == model_id and row.get("candidate_id") == candidate_id:
                return row
    raise ValueError(f"no label for {model_id} {candidate_id}")


def donor_atoms(atoms, exclude: set[tuple[str, int, str]]) -> list[dict[str, object]]:
    donors = []
    for atom in atoms:
        key = (atom["chain"], atom["resseq"], atom["atom_name"])
        if key in exclude:
            continue
        if atom["element"] in {"N", "O"}:
            donors.append(atom)
    return donors


def make_frame(x_axis, y_hint):
    x = unit(x_axis)
    y_proj = sub(y_hint, scale(x, dot(y_hint, x)))
    if norm(y_proj) < 1e-6:
        y_proj = cross(x, (0.0, 0.0, 1.0))
        if norm(y_proj) < 1e-6:
            y_proj = cross(x, (0.0, 1.0, 0.0))
    y = unit(y_proj)
    z = unit(cross(x, y))
    return x, y, z


def coords_in_frame(point, origin, frame):
    rel = sub(point, origin)
    return (dot(rel, frame[0]), dot(rel, frame[1]), dot(rel, frame[2]))


def coords_from_frame(local, origin, frame):
    return add(origin, add(add(scale(frame[0], local[0]), scale(frame[1], local[1])), scale(frame[2], local[2])))


def transformed_ligand(ligand_atoms, label, target_c, target_o, target_l):
    c_atom = find_ligand_atom(ligand_atoms, label["scissile_carbonyl_c_atom_name"])
    o_atom = find_ligand_atom(ligand_atoms, label["scissile_carbonyl_o_atom_name"])
    l_atom = find_ligand_atom(ligand_atoms, label["leaving_o_atom_name"])
    source_origin = c_atom["xyz"]
    source_frame = make_frame(sub(o_atom["xyz"], source_origin), sub(l_atom["xyz"], source_origin))
    target_frame = make_frame(sub(target_o, target_c), sub(target_l, target_c))
    output = []
    for atom in ligand_atoms:
        local = coords_in_frame(atom["xyz"], source_origin, source_frame)
        output.append((atom, coords_from_frame(local, target_c, target_frame)))
    return output


def candidate_targets(protein_atoms, ser_atom, his_atoms, c_o_length: float, c_l_length: float):
    ser = ser_atom["xyz"]
    his_targets = [atom["xyz"] for atom in his_atoms]
    exclude = {(ser_atom["chain"], ser_atom["resseq"], ser_atom["atom_name"])}
    donor_pool = [atom for atom in donor_atoms(protein_atoms, exclude) if distance(atom["xyz"], ser) <= 6.5]
    seeds = []
    for donor in donor_pool:
        donor_xyz = donor["xyz"]
        for dist_sc in (2.9, 3.1, 3.3, 3.5):
            for hx in his_targets:
                # Put the carbon between Ser and His first; this favors leaving-group relay.
                base_dirs = [unit(sub(hx, ser)), unit(sub(donor_xyz, ser))]
                base_dirs.extend([(1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0), (0.0, 0.0, 1.0), (0.0, 0.0, -1.0)])
                for direction in base_dirs:
                    c_target = add(ser, scale(direction, dist_sc))
                    if distance(c_target, donor_xyz) < 0.2:
                        continue
                    o_target = add(c_target, scale(unit(sub(donor_xyz, c_target)), c_o_length))
                    attack = angle(ser, c_target, o_target)
                    oxy = distance(o_target, donor_xyz)
                    if not (80.0 <= attack <= 125.0 and oxy <= 3.2):
                        continue
                    if distance(c_target, hx) < 0.2:
                        continue
                    l_target = add(c_target, scale(unit(sub(hx, c_target)), c_l_length))
                    relay = distance(l_target, hx)
                    if relay > 4.0:
                        continue
                    seeds.append((relay + oxy + abs(105.0 - attack) / 20.0, c_target, o_target, l_target, attack, oxy, relay, donor))
    return sorted(seeds, key=lambda item: item[0])


def write_complex(path: Path, protein_atoms, transformed):
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = []
    serial = 1
    for atom in protein_atoms:
        lines.append(format_atom_line(atom, atom["xyz"], serial=serial, record=str(atom["record"])))
        serial += 1
    lines.append("TER")
    for atom, xyz in transformed:
        lines.append(format_atom_line(atom, xyz, serial=serial, record="HETATM"))
        serial += 1
    lines.append("END")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def generate(args: argparse.Namespace) -> list[dict[str, object]]:
    protein_atoms = parse_pdb(Path(args.protein_pdb))
    ligand_atoms = parse_pdb(Path(args.ligand_pdb))
    label = read_label(Path(args.label_tsv), args.model_id, args.candidate_id)
    ser_chain, _ser_name, ser_resseq = parse_triad_residue(args.triad_ser)
    his_chain, _his_name, his_resseq = parse_triad_residue(args.triad_his)
    ser_atom = find_atom(protein_atoms, ser_chain, ser_resseq, "OG")
    his_atoms = [find_atom(protein_atoms, his_chain, his_resseq, name) for name in ("ND1", "NE2") if any(a["chain"] == his_chain and a["resseq"] == his_resseq and a["atom_name"] == name for a in protein_atoms)]
    if not his_atoms:
        raise ValueError(f"no His nitrogens found for {args.triad_his}")
    c_atom = find_ligand_atom(ligand_atoms, label["scissile_carbonyl_c_atom_name"])
    o_atom = find_ligand_atom(ligand_atoms, label["scissile_carbonyl_o_atom_name"])
    l_atom = find_ligand_atom(ligand_atoms, label["leaving_o_atom_name"])
    targets = candidate_targets(protein_atoms, ser_atom, his_atoms, distance(c_atom["xyz"], o_atom["xyz"]), distance(c_atom["xyz"], l_atom["xyz"]))
    rows = []
    out_dir = Path(args.out_dir)
    for index, (_rank, c_target, o_target, l_target, attack, oxy, relay, donor) in enumerate(targets[: args.max_seeds], start=1):
        pose_id = f"REACTIVE_{args.template_pdb}_{args.model_id}_{args.candidate_id}_{index:03d}"
        complex_pdb = out_dir / f"{pose_id}.pdb"
        transformed = transformed_ligand(ligand_atoms, label, c_target, o_target, l_target)
        write_complex(complex_pdb, protein_atoms, transformed)
        rows.append(
            {
                "pose_id": pose_id,
                "template_pdb": args.template_pdb,
                "substrate_model": args.model_id,
                "candidate_id": args.candidate_id,
                "complex_pdb": str(complex_pdb),
                "ser_og_to_c_target_A": f"{distance(ser_atom['xyz'], c_target):.3f}",
                "attack_angle_target_deg": f"{attack:.2f}",
                "oxyanion_target_distance_A": f"{oxy:.3f}",
                "leaving_his_target_distance_A": f"{relay:.3f}",
                "status": "ready_for_geometry_scoring",
                "source": "blind_reactive_geometry_seed",
                "note": f"donor={donor['chain']}:{donor['resname']}{donor['resseq']}:{donor['atom_name']};not_relaxed",
            }
        )
    if not rows:
        rows.append(
            {
                "pose_id": "pending",
                "template_pdb": args.template_pdb,
                "substrate_model": args.model_id,
                "candidate_id": args.candidate_id,
                "complex_pdb": "pending",
                "ser_og_to_c_target_A": "pending",
                "attack_angle_target_deg": "pending",
                "oxyanion_target_distance_A": "pending",
                "leaving_his_target_distance_A": "pending",
                "status": "not_generated",
                "source": "blind_reactive_geometry_seed",
                "note": "no generic geometry target satisfied",
            }
        )
    write_tsv(out_dir / "reactive_pose_seed_manifest.tsv", rows)
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--protein-pdb", required=True)
    parser.add_argument("--ligand-pdb", required=True)
    parser.add_argument("--label-tsv", required=True)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--candidate-id", required=True)
    parser.add_argument("--template-pdb", required=True)
    parser.add_argument("--triad-ser", required=True)
    parser.add_argument("--triad-his", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--max-seeds", type=int, default=10)
    args = parser.parse_args()
    rows = generate(args)
    print(Path(args.out_dir) / "reactive_pose_seed_manifest.tsv")
    return 0 if any(row["status"] == "ready_for_geometry_scoring" for row in rows) else 2


if __name__ == "__main__":
    raise SystemExit(main())
