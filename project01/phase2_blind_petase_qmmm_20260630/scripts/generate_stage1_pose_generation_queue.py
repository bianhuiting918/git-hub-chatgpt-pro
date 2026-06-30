#!/usr/bin/env python3
"""Generate blind PETase Stage 1 docking/pose-generation queue files."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path


QUEUE_FIELDS = [
    "pose_job_id",
    "pdb",
    "substrate_model",
    "scissile_candidate_id",
    "protein_path",
    "ligand_path",
    "atom_label_path",
    "triad_ser",
    "triad_his",
    "triad_asp",
    "box_center_x",
    "box_center_y",
    "box_center_z",
    "box_size_x",
    "box_size_y",
    "box_size_z",
    "docking_config_path",
    "status",
    "source",
    "note",
]


def read_tsv(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=QUEUE_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def parse_triad_residue(text: str) -> tuple[str, str, int]:
    chain, residue = text.split(":", 1)
    resname = "".join(ch for ch in residue if ch.isalpha())
    resseq = int("".join(ch for ch in residue if ch.isdigit()))
    return chain, resname, resseq


def parse_atom_line(line: str) -> dict[str, object] | None:
    if not line.startswith(("ATOM  ", "HETATM")):
        return None
    try:
        return {
            "atom_name": line[12:16].strip(),
            "resname": line[17:20].strip(),
            "chain": line[21].strip(),
            "resseq": int(line[22:26]),
            "x": float(line[30:38]),
            "y": float(line[38:46]),
            "z": float(line[46:54]),
        }
    except ValueError:
        parts = line.split()
        if len(parts) < 9:
            return None
        return {
            "atom_name": parts[2],
            "resname": parts[3],
            "chain": parts[4],
            "resseq": int(parts[5]),
            "x": float(parts[6]),
            "y": float(parts[7]),
            "z": float(parts[8]),
        }


def find_ser_og(path: Path, ser_text: str) -> tuple[float, float, float]:
    chain, resname, resseq = parse_triad_residue(ser_text)
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            atom = parse_atom_line(line)
            if not atom:
                continue
            if atom["chain"] == chain and atom["resname"] == resname and atom["resseq"] == resseq and atom["atom_name"] == "OG":
                return float(atom["x"]), float(atom["y"]), float(atom["z"])
    raise ValueError(f"Missing Ser OG for {ser_text} in {path}")


def ready_prepared_structures(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    output = []
    for row in rows:
        prepared = row.get("prepared_path", "")
        if not prepared or prepared == "pending":
            continue
        if row.get("protonation_status", "") not in {"reviewed_ph7", "assigned", "not_assigned"}:
            continue
        output.append(row)
    return output


def ready_ligands(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    output = []
    for row in rows:
        if row.get("status") != "built_needs_scissile_candidate_selection":
            continue
        if row.get("pdb_path") in {"", "pending"} or row.get("atom_label_path") in {"", "pending"}:
            continue
        output.append(row)
    return output


def label_candidates(path: Path, model_id: str) -> list[dict[str, str]]:
    rows = []
    for row in read_tsv(path):
        if row.get("model_id") != model_id:
            continue
        if not row.get("candidate_id"):
            continue
        rows.append(row)
    return rows


def write_config(path: Path, protein: Path, ligand: Path, center: tuple[float, float, float], box_size: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "\n".join(
            [
                f"receptor = {protein}",
                f"ligand = {ligand}",
                f"center_x = {center[0]:.3f}",
                f"center_y = {center[1]:.3f}",
                f"center_z = {center[2]:.3f}",
                f"size_x = {box_size:.3f}",
                f"size_y = {box_size:.3f}",
                f"size_z = {box_size:.3f}",
                "exhaustiveness = 16",
                "num_modes = 20",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def pending_row(note: str) -> dict[str, str]:
    return {
        "pose_job_id": "pending",
        "pdb": "pending",
        "substrate_model": "pending",
        "scissile_candidate_id": "pending",
        "protein_path": "pending",
        "ligand_path": "pending",
        "atom_label_path": "pending",
        "triad_ser": "pending",
        "triad_his": "pending",
        "triad_asp": "pending",
        "box_center_x": "pending",
        "box_center_y": "pending",
        "box_center_z": "pending",
        "box_size_x": "pending",
        "box_size_y": "pending",
        "box_size_z": "pending",
        "docking_config_path": "pending",
        "status": "not_ready",
        "source": "blind_workflow_pending",
        "note": note,
    }


def generate(prepared_manifest: Path, triad_manifest: Path, ligand_manifest: Path, out_root: Path, box_size: float) -> list[Path]:
    prepared_rows = ready_prepared_structures(read_tsv(prepared_manifest))
    triads = {row["pdb"]: row for row in read_tsv(triad_manifest) if row.get("pdb")}
    ligands = ready_ligands(read_tsv(ligand_manifest))
    output_rows: list[dict[str, str]] = []
    config_dir = out_root / "01_system_setup" / "docking_inputs"

    if not prepared_rows:
        output_rows.append(pending_row("no_ready_prepared_structure"))
    elif not ligands:
        output_rows.append(pending_row("no_ready_ligand_build_manifest"))
    else:
        for prep in prepared_rows:
            pdb_id = prep["pdb"]
            triad = triads.get(pdb_id)
            if not triad:
                continue
            protein = Path(prep["prepared_path"])
            center = find_ser_og(protein, triad["ser"])
            for ligand in ligands:
                label_path = Path(ligand["atom_label_path"])
                candidates = label_candidates(label_path, ligand["model_id"])
                for candidate in candidates:
                    candidate_id = candidate["candidate_id"]
                    job_id = f"POSE_{pdb_id}_{ligand['model_id']}_{candidate_id}"
                    config_path = config_dir / f"{job_id}.config"
                    ligand_path = Path(ligand["pdb_path"])
                    write_config(config_path, protein, ligand_path, center, box_size)
                    output_rows.append(
                        {
                            "pose_job_id": job_id,
                            "pdb": pdb_id,
                            "substrate_model": ligand["model_id"],
                            "scissile_candidate_id": candidate_id,
                            "protein_path": str(protein),
                            "ligand_path": str(ligand_path),
                            "atom_label_path": str(label_path),
                            "triad_ser": triad["ser"],
                            "triad_his": triad["his"],
                            "triad_asp": triad["asp"],
                            "box_center_x": f"{center[0]:.3f}",
                            "box_center_y": f"{center[1]:.3f}",
                            "box_center_z": f"{center[2]:.3f}",
                            "box_size_x": f"{box_size:.3f}",
                            "box_size_y": f"{box_size:.3f}",
                            "box_size_z": f"{box_size:.3f}",
                            "docking_config_path": str(config_path),
                            "status": "ready_for_docking",
                            "source": "blind_structure_geometry",
                            "note": "box_center_from_catalytic_ser_og",
                        }
                    )
    if not output_rows:
        output_rows.append(pending_row("no_matching_prepared_structure_triad_ligand"))

    queue_path = out_root / "01_system_setup" / "pose_generation_queue.tsv"
    write_tsv(queue_path, output_rows)
    return [queue_path]


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared-manifest", default="work/blind_work/01_system_setup/prepared_structure_manifest.tsv")
    parser.add_argument("--triad-manifest", default="work/blind_work/01_system_setup/ser_his_asp_triad_candidates.tsv")
    parser.add_argument("--ligand-build-manifest", default="work/blind_work/01_system_setup/ligand_build/ligand_build_manifest.tsv")
    parser.add_argument("--out-root", default="work/blind_work")
    parser.add_argument("--box-size", type=float, default=22.0)
    args = parser.parse_args()
    generated = generate(Path(args.prepared_manifest), Path(args.triad_manifest), Path(args.ligand_build_manifest), Path(args.out_root), args.box_size)
    for path in generated:
        print(path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())