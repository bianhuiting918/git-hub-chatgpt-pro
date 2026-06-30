#!/usr/bin/env python3
"""Build blind Stage 1 PETase ligand conformers and atom-label manifests.

This script intentionally uses only substrate chemistry supplied in a TSV file.
It does not read paper trajectories, TS structures, paper reaction coordinates,
or paper-derived CV definitions.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def read_smiles_table(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def stable_atom_name(atom) -> str:
    # Four characters keeps the names usable in PDB-style downstream tools.
    return f"{atom.GetSymbol().upper()}{atom.GetIdx() + 1:03d}"[-4:]


def write_tsv(path: Path, fieldnames: list[str], rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def import_rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except Exception as exc:  # pragma: no cover - depends on remote environment
        return None, None, exc
    return Chem, AllChem, None


def assign_pdb_atom_names(Chem, mol, residue_name: str = "LIG") -> None:
    for atom in mol.GetAtoms():
        name = stable_atom_name(atom)
        atom.SetProp("stage1_atom_name", name)
        info = Chem.AtomPDBResidueInfo()
        info.SetName(name.rjust(4))
        info.SetResidueName(residue_name[:3].upper())
        info.SetResidueNumber(1)
        info.SetSerialNumber(atom.GetIdx() + 1)
        atom.SetMonomerInfo(info)


def enumerate_ester_candidates(Chem, mol, model_id: str) -> list[dict[str, object]]:
    patt = Chem.MolFromSmarts("[CX3](=[OX1])-[OX2]")
    matches = mol.GetSubstructMatches(patt)
    rows: list[dict[str, object]] = []
    seen: set[tuple[int, int, int]] = set()
    for match in matches:
        carbonyl_c, carbonyl_o, single_o = match
        key = (carbonyl_c, carbonyl_o, single_o)
        if key in seen:
            continue
        seen.add(key)
        c_atom = mol.GetAtomWithIdx(carbonyl_c)
        o_atom = mol.GetAtomWithIdx(carbonyl_o)
        lg_atom = mol.GetAtomWithIdx(single_o)
        rows.append(
            {
                "model_id": model_id,
                "candidate_id": f"E{len(rows) + 1:02d}",
                "scissile_carbonyl_c_rdkit_idx0": carbonyl_c,
                "scissile_carbonyl_c_atom_name": c_atom.GetProp("stage1_atom_name"),
                "scissile_carbonyl_o_rdkit_idx0": carbonyl_o,
                "scissile_carbonyl_o_atom_name": o_atom.GetProp("stage1_atom_name"),
                "leaving_o_rdkit_idx0": single_o,
                "leaving_o_atom_name": lg_atom.GetProp("stage1_atom_name"),
                "selection_status": "candidate_not_selected",
                "selection_basis": "enumerated_from_substrate_chemistry_not_paper_result",
            }
        )
    return rows


def build_conformers(Chem, AllChem, row: dict[str, str], out_dir: Path, nconf: int, seed: int) -> dict[str, object]:
    model_id = row["model_id"]
    smiles = row.get("smiles", "").strip()
    if not smiles:
        return {
            "model_id": model_id,
            "status": "skipped_no_free_ligand_smiles",
            "formal_charge": row.get("formal_charge", ""),
            "conformer_count": 0,
            "sdf_path": "",
            "pdb_path": "",
            "atom_label_path": "",
            "sha256_sdf": "",
            "sha256_pdb": "",
            "note": "covalent precursor must be generated from protein complex",
        }

    mol0 = Chem.MolFromSmiles(smiles)
    if mol0 is None:
        return {
            "model_id": model_id,
            "status": "failed_invalid_smiles",
            "formal_charge": row.get("formal_charge", ""),
            "conformer_count": 0,
            "sdf_path": "",
            "pdb_path": "",
            "atom_label_path": "",
            "sha256_sdf": "",
            "sha256_pdb": "",
            "note": smiles,
        }

    expected_charge = row.get("formal_charge", "")
    actual_charge = Chem.GetFormalCharge(mol0)
    if expected_charge not in ("", "pka_dependent") and int(expected_charge) != actual_charge:
        return {
            "model_id": model_id,
            "status": "failed_charge_mismatch",
            "formal_charge": expected_charge,
            "conformer_count": 0,
            "sdf_path": "",
            "pdb_path": "",
            "atom_label_path": "",
            "sha256_sdf": "",
            "sha256_pdb": "",
            "note": f"rdkit_charge={actual_charge}",
        }

    mol = Chem.AddHs(mol0)
    assign_pdb_atom_names(Chem, mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = seed
    params.pruneRmsThresh = 0.25
    conf_ids = list(AllChem.EmbedMultipleConfs(mol, numConfs=nconf, params=params))
    if not conf_ids:
        return {
            "model_id": model_id,
            "status": "failed_no_conformers",
            "formal_charge": actual_charge,
            "conformer_count": 0,
            "sdf_path": "",
            "pdb_path": "",
            "atom_label_path": "",
            "sha256_sdf": "",
            "sha256_pdb": "",
            "note": "EmbedMultipleConfs returned zero conformers",
        }

    if AllChem.MMFFHasAllMoleculeParams(mol):
        for cid in conf_ids:
            AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=500)
        force_field = "MMFF94"
    else:
        for cid in conf_ids:
            AllChem.UFFOptimizeMolecule(mol, confId=cid, maxIters=500)
        force_field = "UFF"

    lig_dir = out_dir / model_id
    lig_dir.mkdir(parents=True, exist_ok=True)
    sdf_path = lig_dir / f"{model_id}_conformers.sdf"
    pdb_path = lig_dir / f"{model_id}_conf001_named_atoms.pdb"
    labels_path = out_dir / "qm_atom_labels" / f"{model_id}_atoms.tsv"

    writer = Chem.SDWriter(str(sdf_path))
    for serial, cid in enumerate(conf_ids, start=1):
        mol.SetProp("_Name", f"{model_id}_conf{serial:03d}")
        mol.SetProp("model_id", model_id)
        mol.SetProp("source_smiles", smiles)
        mol.SetProp("stage1_force_field", force_field)
        mol.SetProp("stage1_conformer_id", f"conf{serial:03d}")
        writer.write(mol, confId=cid)
    writer.close()
    Chem.MolToPDBFile(mol, str(pdb_path), confId=conf_ids[0])

    ester_rows = enumerate_ester_candidates(Chem, mol, model_id)
    write_tsv(
        labels_path,
        [
            "model_id",
            "candidate_id",
            "scissile_carbonyl_c_rdkit_idx0",
            "scissile_carbonyl_c_atom_name",
            "scissile_carbonyl_o_rdkit_idx0",
            "scissile_carbonyl_o_atom_name",
            "leaving_o_rdkit_idx0",
            "leaving_o_atom_name",
            "selection_status",
            "selection_basis",
        ],
        ester_rows,
    )

    return {
        "model_id": model_id,
        "status": "built_needs_scissile_candidate_selection",
        "formal_charge": actual_charge,
        "conformer_count": len(conf_ids),
        "sdf_path": str(sdf_path),
        "pdb_path": str(pdb_path),
        "atom_label_path": str(labels_path),
        "sha256_sdf": sha256_file(sdf_path),
        "sha256_pdb": sha256_file(pdb_path),
        "note": f"force_field={force_field};ester_candidate_count={len(ester_rows)}",
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smiles-table", default="project01/phase2_blind_petase_qmmm_20260630/blind_work/01_system_setup/ligand_smiles.tsv")
    parser.add_argument("--out-dir", default="work/blind_work/01_system_setup/ligand_build")
    parser.add_argument("--conformers", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20260630)
    args = parser.parse_args()

    smiles_table = Path(args.smiles_table)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    manifest_path = out_dir / "ligand_build_manifest.tsv"

    Chem, AllChem, import_error = import_rdkit()
    if import_error is not None:
        rows = [
            {
                "model_id": "ALL",
                "status": "blocked_missing_rdkit",
                "formal_charge": "",
                "conformer_count": 0,
                "sdf_path": "",
                "pdb_path": "",
                "atom_label_path": "",
                "sha256_sdf": "",
                "sha256_pdb": "",
                "note": repr(import_error),
            }
        ]
        write_tsv(
            manifest_path,
            ["model_id", "status", "formal_charge", "conformer_count", "sdf_path", "pdb_path", "atom_label_path", "sha256_sdf", "sha256_pdb", "note"],
            rows,
        )
        print(f"RDKit unavailable; wrote {manifest_path}", file=sys.stderr)
        return 2

    if not smiles_table.exists():
        print(f"Missing smiles table: {smiles_table}", file=sys.stderr)
        return 2

    rows = [
        build_conformers(Chem, AllChem, row, out_dir, args.conformers, args.seed)
        for row in read_smiles_table(smiles_table)
    ]
    write_tsv(
        manifest_path,
        ["model_id", "status", "formal_charge", "conformer_count", "sdf_path", "pdb_path", "atom_label_path", "sha256_sdf", "sha256_pdb", "note"],
        rows,
    )
    print(f"Wrote {manifest_path}")
    return 0 if all(str(row["status"]).startswith(("built", "skipped")) for row in rows) else 1


if __name__ == "__main__":
    raise SystemExit(main())
