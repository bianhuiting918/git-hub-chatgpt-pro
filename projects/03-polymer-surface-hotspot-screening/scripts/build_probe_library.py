#!/usr/bin/env python3
"""Validate the frozen probe manifest and build a deterministic 3D SDF library."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import tempfile
from datetime import datetime, timezone
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem, Lipinski, rdMolDescriptors


REQUIRED_COLUMNS = {
    "probe_id",
    "material_family",
    "role",
    "preferred_name",
    "canonical_smiles",
    "formal_charge",
    "heavy_atom_count",
    "total_atom_count",
    "primary_or_control",
    "stage",
    "rationale",
    "source",
    "status",
}
DERIVED_COLUMNS = ["molecular_formula", "rotatable_bonds", "validation_status"]
OUTPUT_NAMES = (
    "probes.sdf",
    "validated_probes.tsv",
    "PROBES_PASS.json",
    "PROBE_SHA256SUMS",
)
EMBED_SEED = 0xC0FFEE


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_manifest(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        missing = REQUIRED_COLUMNS.difference(fieldnames)
        if missing:
            raise ValueError(f"manifest is missing required columns: {sorted(missing)}")
        rows = list(reader)
    if not rows:
        raise ValueError("probe manifest is empty")
    probe_ids = [row["probe_id"] for row in rows]
    if any(not probe_id for probe_id in probe_ids):
        raise ValueError("probe_id cannot be empty")
    if len(probe_ids) != len(set(probe_ids)):
        raise ValueError("probe_id values must be unique")
    return fieldnames, rows


def validate_and_embed(row: dict[str, str]) -> tuple[Chem.Mol, dict[str, str]]:
    probe_id = row["probe_id"]
    if not row["source"].strip():
        raise ValueError(f"{probe_id}: source is required")
    mol = Chem.MolFromSmiles(row["canonical_smiles"])
    if mol is None:
        raise ValueError(f"{probe_id}: RDKit could not parse canonical_smiles")

    canonical = Chem.MolToSmiles(mol, canonical=True)
    if canonical != row["canonical_smiles"]:
        raise ValueError(
            f"{probe_id}: declared canonical_smiles={row['canonical_smiles']!r} "
            f"but RDKit canonicalized to {canonical!r}"
        )

    observed = {
        "formal_charge": Chem.GetFormalCharge(mol),
        "heavy_atom_count": mol.GetNumHeavyAtoms(),
        "total_atom_count": Chem.AddHs(mol).GetNumAtoms(),
    }
    for field, value in observed.items():
        if int(row[field]) != value:
            raise ValueError(
                f"{probe_id}: declared {field}={row[field]} but RDKit computed {value}"
            )

    mol_h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = EMBED_SEED
    if AllChem.EmbedMolecule(mol_h, params) != 0:
        raise ValueError(f"{probe_id}: deterministic ETKDG embedding failed")
    if AllChem.MMFFHasAllMoleculeParams(mol_h):
        AllChem.MMFFOptimizeMolecule(mol_h, maxIters=500)
    else:
        AllChem.UFFOptimizeMolecule(mol_h, maxIters=500)

    derived = {
        "molecular_formula": rdMolDescriptors.CalcMolFormula(mol),
        "rotatable_bonds": str(Lipinski.NumRotatableBonds(mol)),
        "validation_status": "PASS",
    }
    for key, value in row.items():
        mol_h.SetProp(key, value)
    for key, value in derived.items():
        mol_h.SetProp(key, value)
    mol_h.SetProp("_Name", probe_id)
    return mol_h, derived


def write_library(
    manifest: Path,
    output_dir: Path,
    fieldnames: list[str],
    rows: list[dict[str, str]],
) -> None:
    active = [row for row in rows if row["status"] != "deferred"]
    deferred = [row for row in rows if row["status"] == "deferred"]
    if not active:
        raise ValueError("probe manifest has no active probes")
    for row in deferred:
        if not row["source"].strip():
            raise ValueError(f"{row['probe_id']}: deferred source is required")

    output_dir.mkdir(parents=True, exist_ok=True)
    gate_path = output_dir / "PROBES_PASS.json"
    if gate_path.exists():
        print(f"existing probe PASS retained without overwrite: {gate_path}")
        return
    preexisting = [name for name in OUTPUT_NAMES if (output_dir / name).exists()]
    if preexisting:
        raise FileExistsError(
            f"refusing to overwrite incomplete probe outputs: {preexisting}"
        )

    with tempfile.TemporaryDirectory(prefix=".probe-build-", dir=output_dir) as temp:
        temp_dir = Path(temp)
        sdf_path = temp_dir / "probes.sdf"
        table_path = temp_dir / "validated_probes.tsv"
        temporary_rows: list[dict[str, str]] = []

        writer = Chem.SDWriter(str(sdf_path))
        try:
            for row in active:
                mol, derived = validate_and_embed(row)
                writer.write(mol)
                temporary_rows.append({**row, **derived})
        finally:
            writer.close()

        with table_path.open("w", encoding="utf-8", newline="") as handle:
            table_writer = csv.DictWriter(
                handle,
                fieldnames=fieldnames + DERIVED_COLUMNS,
                delimiter="\t",
                lineterminator="\n",
                extrasaction="raise",
            )
            table_writer.writeheader()
            table_writer.writerows(temporary_rows)

        gate = {
            "status": "PROBES_PASS",
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "scientific_scope": "probe_library_only",
            "manifest_sha256": sha256(manifest),
            "sdf_sha256": sha256(sdf_path),
            "validated_table_sha256": sha256(table_path),
            "n_active": len(active),
            "n_deferred": len(deferred),
            "n_pet_active": sum(
                row["material_family"] == "PET" for row in active
            ),
            "n_nylon_active": sum(
                row["material_family"] != "PET" for row in active
            ),
            "rdkit_version": Chem.rdBase.rdkitVersion,
            "embedding_method": "ETKDGv3_then_MMFF_or_UFF",
            "embedding_seed": EMBED_SEED,
        }
        temporary_gate = temp_dir / "PROBES_PASS.json"
        temporary_gate.write_text(
            json.dumps(gate, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )

        for name in ("probes.sdf", "validated_probes.tsv", "PROBES_PASS.json"):
            os.replace(temp_dir / name, output_dir / name)

    checksum_lines = []
    for name in ("probes.sdf", "validated_probes.tsv", "PROBES_PASS.json"):
        checksum_lines.append(f"{sha256(output_dir / name)}  {name}")
    (output_dir / "PROBE_SHA256SUMS").write_text(
        "\n".join(checksum_lines) + "\n",
        encoding="utf-8",
    )
    print(f"PROBES_PASS.json written: {gate_path}")


def main() -> int:
    args = parse_args()
    fieldnames, rows = load_manifest(args.manifest)
    write_library(args.manifest, args.output_dir, fieldnames, rows)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
