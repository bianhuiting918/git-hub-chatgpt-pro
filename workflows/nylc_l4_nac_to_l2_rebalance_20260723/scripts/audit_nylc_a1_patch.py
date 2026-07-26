#!/usr/bin/env python3
"""Independently audit the generated NylC A1 local screening patch."""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Any, Dict, Mapping, Optional, Sequence, Tuple

PASS_STATUS = "PASS_A1_SCREENING_PATCH"
REQUIRED_ATOMS = {
    "N", "H1", "H2", "HG1", "CA", "HA", "CB", "HB", "CG2",
    "HG21", "HG22", "HG23", "OG1", "C", "O",
    "NCAP", "HCAP", "CCAP", "HC1", "HC2", "HC3",
}
REQUIRED_PROVENANCE = {
    "model_definition", "net_charge_e", "charge_method", "atom_types",
    "bonded_parameters", "tool_versions", "input_sha256", "output_sha256",
    "validation_status", "parmchk2_unresolved_count", "tleap_load_status",
    "gromacs_grompp_status", "amber_energy_kcal_mol", "gromacs_energy_kj_mol",
}


class PatchAuditError(ValueError):
    """Raised when a generated patch is unsuitable even for MM screening."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise PatchAuditError(message)


def parse_mol2(text: str) -> Tuple[Dict[int, Dict[str, Any]], set]:
    section = ""
    atoms: Dict[int, Dict[str, Any]] = {}
    bonds = set()
    for raw in text.splitlines():
        line = raw.strip()
        if line.startswith("@<TRIPOS>"):
            section = line[len("@<TRIPOS>") :].lower()
            continue
        if not line:
            continue
        fields = line.split()
        if section == "atom":
            _require(len(fields) >= 9, "malformed MOL2 atom row")
            atom_id = int(fields[0])
            atoms[atom_id] = {
                "name": fields[1],
                "type": fields[5],
                "charge": float(fields[8]),
            }
        elif section == "bond":
            _require(len(fields) >= 4, "malformed MOL2 bond row")
            bonds.add(frozenset((int(fields[1]), int(fields[2]))))
    _require(bool(atoms), "MOL2 has no atoms")
    return atoms, bonds


def audit_patch(
    mol2_text: str,
    frcmod_text: str,
    provenance: Mapping[str, Any],
) -> Dict[str, Any]:
    atoms, bond_ids = parse_mol2(mol2_text)
    names = {row["name"]: atom_id for atom_id, row in atoms.items()}
    _require(set(names) == REQUIRED_ATOMS, "A1 capped-model atom set differs from authority")

    charge = sum(float(row["charge"]) for row in atoms.values())
    _require(abs(charge) <= 1.0e-5, f"A1 capped-model charge is {charge:.8f}, expected zero")
    _require(
        all(row["type"].lower() not in {"du", "xx", "x"} for row in atoms.values()),
        "untyped atom remains in A1 MOL2",
    )

    def has_bond(left: str, right: str) -> bool:
        return frozenset((names[left], names[right])) in bond_ids

    for hydrogen in ("H1", "H2", "HG1"):
        _require(has_bond("N", hydrogen), f"N-{hydrogen} bond is missing")
    _require(has_bond("N", "CA"), "N-CA bond is missing")
    _require(not has_bond("OG1", "HG1"), "forbidden OG1-HG1 bond remains")
    _require(has_bond("OG1", "CB"), "OG1-CB bond is missing")
    og1_neighbors = [
        row["name"]
        for atom_id, row in atoms.items()
        if frozenset((names["OG1"], atom_id)) in bond_ids
    ]
    _require(og1_neighbors == ["CB"] or set(og1_neighbors) == {"CB"}, "OG1 valence is not alkoxide-like")

    lowered = frcmod_text.lower()
    unresolved_markers = ("attn, need revision", "unresolved", "missing parameter")
    _require(
        not any(marker in lowered for marker in unresolved_markers),
        "unresolved parmchk2 term is present",
    )

    missing = REQUIRED_PROVENANCE.difference(provenance)
    _require(not missing, "missing provenance fields: " + ",".join(sorted(missing)))
    _require(provenance["validation_status"] == "PASS", "parameter validation did not pass")
    _require(int(provenance["net_charge_e"]) == 0, "provenance net charge is not zero")
    _require(int(provenance["parmchk2_unresolved_count"]) == 0, "unresolved parmchk2 count is nonzero")
    _require(provenance["tleap_load_status"] == "PASS", "tleap load did not pass")
    _require(provenance["gromacs_grompp_status"] == "PASS", "GROMACS grompp did not pass")
    for key in ("input_sha256", "output_sha256"):
        _require(bool(re.fullmatch(r"[0-9a-f]{64}", str(provenance[key]))), f"invalid {key}")
    for key in ("amber_energy_kcal_mol", "gromacs_energy_kj_mol"):
        try:
            value = float(provenance[key])
        except (TypeError, ValueError) as exc:
            raise PatchAuditError(f"{key} is not finite") from exc
        _require(math.isfinite(value), f"{key} is not finite")

    return {
        "status": PASS_STATUS,
        "atom_count": len(atoms),
        "bond_count": len(bond_ids),
        "charge_e": charge,
        "nalpha_hydrogens": ["H1", "H2", "HG1"],
        "ogamma_hydrogens": [],
        "charge_method": provenance["charge_method"],
        "parameter_scope": "classical A1 NAC stability screening",
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--mol2", type=Path, required=True)
    parser.add_argument("--frcmod", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        result = audit_patch(
            args.mol2.read_text(encoding="utf-8"),
            args.frcmod.read_text(encoding="utf-8"),
            json.loads(args.provenance.read_text(encoding="utf-8")),
        )
    except (OSError, json.JSONDecodeError, PatchAuditError) as exc:
        print(f"FAIL_A1_PARAMETER_AUDIT: {exc}", file=sys.stderr)
        return 2
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
