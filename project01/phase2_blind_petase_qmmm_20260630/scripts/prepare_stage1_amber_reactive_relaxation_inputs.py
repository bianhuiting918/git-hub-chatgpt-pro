#!/usr/bin/env python3
"""Prepare Amber MM restrained-relaxation inputs for docking-derived PETase poses.

This stage sits before QM/MM. It starts from no-clash docking poses and creates
classical Amber inputs that bias the ligand toward generic serine-hydrolase
Michaelis geometry without using paper transition-state coordinates.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


FIELDS = [
    "relaxation_job_id",
    "pose_id",
    "template_pdb",
    "substrate_model",
    "generation_method",
    "source_complex_pdb",
    "complex_for_amber_pdb_path",
    "relaxation_job_dir",
    "restraint_path",
    "mm_min_mdin_path",
    "runner_path",
    "status",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def sanitize_id(raw: str) -> str:
    return "".join(char if char.isalnum() or char in {"_", "-"} else "_" for char in raw)


def parse_atom_selector(text: str) -> tuple[str, int, str]:
    chain, resseq, atom_name = text.split(":", 2)
    return chain, int(resseq), atom_name


def renumber_pdb_for_amber(source: Path, destination: Path) -> None:
    lines = []
    serial = 1
    with source.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")):
                lines.append(f"{line[:6]}{serial:5d}{line[11:].rstrip()}")
                serial += 1
            elif line.startswith(("TER", "END")):
                lines.append(line.rstrip())
    if not lines or lines[-1] != "END":
        lines.append("END")
    destination.write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_pdb(path: Path) -> list[dict[str, object]]:
    atoms = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            try:
                atoms.append(
                    {
                        "serial": int(line[6:11]),
                        "record": line[:6].strip(),
                        "atom_name": line[12:16].strip(),
                        "resname": line[17:20].strip(),
                        "chain": line[21].strip(),
                        "resseq": int(line[22:26]),
                    }
                )
            except ValueError:
                continue
    return atoms


def find_protein_atom_serial(atoms: list[dict[str, object]], selector: str) -> int:
    chain, resseq, atom_name = parse_atom_selector(selector)
    for atom in atoms:
        if (
            atom["record"] == "ATOM"
            and atom["chain"] == chain
            and atom["resseq"] == resseq
            and atom["atom_name"] == atom_name
        ):
            return int(atom["serial"])
    raise ValueError(f"missing protein atom {selector}")


def find_ligand_atom_serial(atoms: list[dict[str, object]], atom_name: str) -> int:
    for atom in atoms:
        if atom["record"] == "HETATM" and atom["resname"] == "LIG" and atom["atom_name"] == atom_name:
            return int(atom["serial"])
    raise ValueError(f"missing ligand atom {atom_name}")


def read_label(path: Path, model_id: str, candidate_id: str) -> dict[str, str]:
    for row in read_tsv(path):
        if row.get("model_id") == model_id and row.get("candidate_id") == candidate_id:
            return row
    raise ValueError(f"no atom label for {model_id} {candidate_id}")


def restraint_distance(i: int, j: int, target: float, force: float, comment: str) -> str:
    lower_1 = max(0.0, target - 0.2)
    lower_2 = max(0.0, target - 0.1)
    upper_2 = target + 0.1
    upper_1 = target + 0.2
    return f"""# {comment}
&rst
  iat={i},{j},
  r1={lower_1:.2f}, r2={lower_2:.2f}, r3={upper_2:.2f}, r4={upper_1:.2f},
  rk2={force:.1f}, rk3={force:.1f},
/
"""


def restraint_angle(i: int, j: int, k: int, target_deg: float, force: float, comment: str) -> str:
    return f"""# {comment}
&rst
  iat={i},{j},{k},
  r1={target_deg - 15.0:.2f}, r2={target_deg - 5.0:.2f}, r3={target_deg + 5.0:.2f}, r4={target_deg + 15.0:.2f},
  rk2={force:.1f}, rk3={force:.1f},
/
"""


def restraints(
    ser_og: int,
    carbonyl_c: int,
    carbonyl_o: int,
    leaving_o: int,
    his_acceptor: int,
    oxyanion_donor: int,
    ser_c_target: float,
    leaving_his_target: float,
    oxyanion_target: float,
    angle_target: float,
    distance_force: float,
    angle_force: float,
) -> str:
    return "\n".join(
        [
            restraint_distance(ser_og, carbonyl_c, ser_c_target, distance_force, "Ser OG to scissile carbonyl carbon"),
            restraint_angle(ser_og, carbonyl_c, carbonyl_o, angle_target, angle_force, "Ser OG - carbonyl C - carbonyl O attack angle"),
            restraint_distance(leaving_o, his_acceptor, leaving_his_target, distance_force, "Leaving oxygen to catalytic His acceptor"),
            restraint_distance(carbonyl_o, oxyanion_donor, oxyanion_target, distance_force, "Carbonyl oxygen to oxyanion-hole donor"),
        ]
    )


def mm_min_mdin(maxcyc: int) -> str:
    return f"""Stage1 reactive-pose restrained MM minimization
&cntrl
  imin=1, maxcyc={maxcyc}, ncyc={max(1, maxcyc // 2)},
  ntb=1, cut=10.0,
  ntpr=10,
  nmropt=1,
/
DISANG=reactive_relaxation_restraints.RST
LISTOUT=POUT
"""


def runner() -> str:
    return """#!/usr/bin/env bash
set -euo pipefail

SANDER=${SANDER:-sander}
PRMTOP=${PRMTOP:?set PRMTOP to the Amber topology file for this docking complex}
INPCRD=${INPCRD:?set INPCRD to the Amber coordinate file for this docking complex}

"${SANDER}" -O -i 00_restrained_mm_min.in -o 00_restrained_mm_min.out \\
  -p "${PRMTOP}" -c "${INPCRD}" -r 00_restrained_mm_min.rst7 -ref "${INPCRD}"
"""


def readme(row: dict[str, str]) -> str:
    return "\n".join(
        [
            f"# Stage1 Amber Reactive Relaxation {row.get('pose_id', '')}",
            "",
            "Purpose: relax a no-clash docking pose toward generic Michaelis geometry before any QM/MM scan.",
            "",
            "This is a classical restrained MM preparation step, not a transition-state search.",
            "It uses generic serine-hydrolase restraints only: Ser attack distance, attack angle, oxyanion contact, and His relay contact.",
            "Build Amber topology and coordinates from `complex_for_amber.pdb`; the restraint atom indices are based on that renumbered PDB order.",
            "",
            "Boundary: no paper TS coordinates, paper reaction-coordinate terms, umbrella windows, barriers, or mechanism conclusions are used.",
        ]
    ) + "\n"


def selected_rows(rows: list[dict[str, str]], min_heavy_contact: float, max_poses: int) -> list[dict[str, str]]:
    accepted = []
    for row in rows:
        value = row.get("min_ligand_protein_heavy_contact_A", "")
        if value and value.lower() != "inf" and float(value) < min_heavy_contact:
            continue
        if row.get("complex_pdb"):
            accepted.append(row)
        if len(accepted) >= max_poses:
            break
    return accepted


def write_job(
    row: dict[str, str],
    label: dict[str, str],
    out_dir: Path,
    ser_og_selector: str,
    his_acceptor_selector: str,
    oxyanion_donor_selector: str,
    maxcyc: int,
    ser_c_target: float,
    leaving_his_target: float,
    oxyanion_target: float,
    angle_target: float,
    distance_force: float,
    angle_force: float,
) -> dict[str, str]:
    pose_id = row["pose_id"]
    source_complex = Path(row["complex_pdb"])
    job_id = f"AMBER_RELAX_{sanitize_id(pose_id)}"
    job_dir = out_dir / job_id
    job_dir.mkdir(parents=True, exist_ok=True)
    complex_for_amber = job_dir / "complex_for_amber.pdb"
    renumber_pdb_for_amber(source_complex, complex_for_amber)

    atoms = parse_pdb(complex_for_amber)
    ser_og = find_protein_atom_serial(atoms, ser_og_selector)
    his_acceptor = find_protein_atom_serial(atoms, his_acceptor_selector)
    oxyanion_donor = find_protein_atom_serial(atoms, oxyanion_donor_selector)
    carbonyl_c = find_ligand_atom_serial(atoms, label["scissile_carbonyl_c_atom_name"])
    carbonyl_o = find_ligand_atom_serial(atoms, label["scissile_carbonyl_o_atom_name"])
    leaving_o = find_ligand_atom_serial(atoms, label["leaving_o_atom_name"])

    restraint_path = job_dir / "reactive_relaxation_restraints.RST"
    restraint_path.write_text(
        restraints(
            ser_og,
            carbonyl_c,
            carbonyl_o,
            leaving_o,
            his_acceptor,
            oxyanion_donor,
            ser_c_target,
            leaving_his_target,
            oxyanion_target,
            angle_target,
            distance_force,
            angle_force,
        ),
        encoding="utf-8",
    )
    mdin_path = job_dir / "00_restrained_mm_min.in"
    mdin_path.write_text(mm_min_mdin(maxcyc), encoding="utf-8")
    runner_path = job_dir / "run_reactive_mm_relaxation.sh"
    runner_path.write_text(runner(), encoding="utf-8")
    runner_path.chmod(0o755)
    (job_dir / "00_README.md").write_text(readme(row), encoding="utf-8")

    return {
        "relaxation_job_id": job_id,
        "pose_id": pose_id,
        "template_pdb": row.get("template_pdb", ""),
        "substrate_model": row.get("substrate_model", ""),
        "generation_method": row.get("generation_method", ""),
        "source_complex_pdb": str(source_complex),
        "complex_for_amber_pdb_path": str(complex_for_amber),
        "relaxation_job_dir": str(job_dir),
        "restraint_path": str(restraint_path),
        "mm_min_mdin_path": str(mdin_path),
        "runner_path": str(runner_path),
        "status": "reactive_relaxation_inputs_ready_requires_amber_topology",
    }


def generate(args: argparse.Namespace) -> Path:
    label = read_label(Path(args.label_tsv), args.model_id, args.candidate_id)
    rows = [
        write_job(
            row,
            label,
            Path(args.out_dir),
            args.ser_og,
            args.his_acceptor,
            args.oxyanion_donor,
            args.maxcyc,
            args.ser_c_target,
            args.leaving_his_target,
            args.oxyanion_target,
            args.angle_target,
            args.distance_force,
            args.angle_force,
        )
        for row in selected_rows(read_tsv(Path(args.pose_score_summary)), args.min_ligand_protein_heavy_contact, args.max_poses)
    ]
    manifest = Path(args.out_dir) / "amber_reactive_relaxation_manifest.tsv"
    write_tsv(manifest, rows, FIELDS)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pose-score-summary", required=True)
    parser.add_argument("--label-tsv", required=True)
    parser.add_argument("--model-id", required=True)
    parser.add_argument("--candidate-id", required=True)
    parser.add_argument("--ser-og", required=True)
    parser.add_argument("--his-acceptor", required=True)
    parser.add_argument("--oxyanion-donor", required=True)
    parser.add_argument("--out-dir", required=True)
    parser.add_argument("--max-poses", type=int, default=5)
    parser.add_argument("--maxcyc", type=int, default=500)
    parser.add_argument("--min-ligand-protein-heavy-contact", type=float, default=1.2)
    parser.add_argument("--ser-c-target", type=float, default=2.8)
    parser.add_argument("--leaving-his-target", type=float, default=2.8)
    parser.add_argument("--oxyanion-target", type=float, default=2.8)
    parser.add_argument("--angle-target", type=float, default=105.0)
    parser.add_argument("--distance-force", type=float, default=50.0)
    parser.add_argument("--angle-force", type=float, default=10.0)
    args = parser.parse_args()
    manifest = generate(args)
    print(manifest)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
