#!/usr/bin/env python3
"""Audit the immutable NylC A1 authority and its local parameter provenance."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, Mapping, Optional

PASS_STATUS = "PASS_A1_SCREENING_PATCH"
PENDING_STATUS = "NOT_EVALUATED_PENDING_GENERATION"
SHA256_RE = re.compile(r"^[0-9a-f]{64}$")
EXPECTED_PARENT_IDS = [
    "nac_evt18_time1206ps",
    "nac_evt25_time1462ps",
    "nac_evt08_time1086ps",
]
EXPECTED_TARGET = {
    "Thr267_Nalpha": "NH3+",
    "Thr267_Ogamma": "O-",
    "Asp306": "ASH",
    "Asp308": "ASP",
}
ALLOWED_LOCAL_ATOMS = {
    "Thr267:N",
    "Thr267:H1",
    "Thr267:H2",
    "Thr267:HG1",
    "Thr267:CA",
    "Thr267:CB",
    "Thr267:OG1",
    "Thr267:CG2",
    "Thr267:C",
    "Thr267:O",
}
REQUIRED_PROVENANCE = {
    "model_definition",
    "net_charge_e",
    "charge_method",
    "atom_types",
    "bonded_parameters",
    "tool_versions",
    "input_sha256",
    "output_sha256",
    "validation_status",
}


class AuditError(ValueError):
    """Raised when the authority or patch provenance violates the contract."""


def _require(condition: bool, message: str) -> None:
    if not condition:
        raise AuditError(message)


def _same_bonds(value: Any, expected: Iterable[Iterable[str]]) -> bool:
    return {frozenset(row) for row in value or []} == {
        frozenset(row) for row in expected
    }


def required_patch_fixture() -> Dict[str, Any]:
    """Return a complete synthetic patch record for contract tests only."""
    return {
        "provenance": {
            "model_definition": "N-terminal Thr(A1)-N-methylamide capped model",
            "net_charge_e": 0,
            "charge_method": "AM1-BCC",
            "atom_types": "AmberTools antechamber plus audited crosswalk",
            "bonded_parameters": "parmchk2 with zero unresolved terms",
            "tool_versions": {"ambertools": "2018"},
            "input_sha256": "0" * 64,
            "output_sha256": "1" * 64,
            "validation_status": "PASS",
        },
        "residue_charge_e": 0.0,
        "bonded_state": {
            "N_hydrogens": ["H1", "H2", "HG1"],
            "OG1_hydrogens": [],
        },
        "modified_atoms": sorted(ALLOWED_LOCAL_ATOMS),
    }


def audit_authority(
    authority: Mapping[str, Any],
    *,
    require_ready: bool = False,
) -> Dict[str, Any]:
    """Validate authority data and optionally require an MM-ready patch."""

    _require(authority.get("schema_version") == 1, "unsupported schema_version")
    audit = authority.get("parent_audit", {})
    _require(audit.get("job_id") == 61841413, "wrong parent audit job")
    _require(
        audit.get("scientific_status") == "PASS_UNRESTRAINED_M1_ENSEMBLE_NAC",
        "parent ensemble lacks scientific NAC PASS",
    )

    parents = authority.get("parents", [])
    _require(
        [row.get("id") for row in parents] == EXPECTED_PARENT_IDS,
        "parent set/order differs from frozen top3",
    )
    for row in parents:
        _require(Path(str(row.get("source_path", ""))).is_absolute(), "parent source path is not absolute")
        _require(bool(SHA256_RE.fullmatch(str(row.get("coordinate_sha256", "")))), "invalid parent SHA256")
        _require(float(row.get("distance_nm", 9.0)) <= 0.35, "parent distance exceeds NAC gate")
        _require(95.0 <= float(row.get("angle_deg", 0.0)) <= 115.0, "parent angle exceeds NAC gate")

    _require(authority.get("target_microstate") == EXPECTED_TARGET, "wrong target microstate")
    delta = authority.get("topology_delta", {})
    _require(delta.get("transferred_atom") == "Thr267:HG1", "wrong transferred atom")
    _require(
        _same_bonds(delta.get("removed_bonds"), [["Thr267:OG1", "Thr267:HG1"]]),
        "wrong removed bond",
    )
    _require(
        _same_bonds(delta.get("added_bonds"), [["Thr267:N", "Thr267:HG1"]]),
        "wrong added bond",
    )

    invariants = authority.get("invariants", {})
    _require(invariants.get("atom_count_delta") == 0, "atom count must be unchanged")
    _require(float(invariants.get("system_charge_delta_e", 999)) == 0.0, "system charge must be unchanged")
    _require(
        set(invariants.get("unchanged_microstates", [])) == {"Asp306:ASH", "Asp308:ASP"},
        "Asp306/Asp308 microstate invariant is missing",
    )
    _require(authority.get("force_field", {}).get("name") == "amber99sb-ildn", "wrong force field")

    parameter_source = authority.get("parameter_source", {})
    status = parameter_source.get("status")
    _require(
        status in {PENDING_STATUS, PASS_STATUS, "FAIL_A1_PARAMETER_AUDIT"},
        "unknown parameter-source status",
    )

    if require_ready:
        _require(status == PASS_STATUS, "parameter source is not ready for classical A1 MM")
        provenance = parameter_source.get("provenance")
        _require(isinstance(provenance, Mapping), "missing provenance")
        missing = REQUIRED_PROVENANCE.difference(provenance)
        _require(not missing, "missing provenance fields: " + ",".join(sorted(missing)))
        _require(provenance.get("validation_status") == "PASS", "parameter validation did not pass")
        _require(float(parameter_source.get("residue_charge_e", 999.0)) == 0.0, "A1 residue charge is not zero")

        bonded = parameter_source.get("bonded_state", {})
        _require(
            set(bonded.get("N_hydrogens", [])) == {"H1", "H2", "HG1"},
            "Nalpha is not NH3+ in the audited patch",
        )
        _require(not bonded.get("OG1_hydrogens", []), "OG1 remains protonated in the audited patch")

        modified_atoms = set(parameter_source.get("modified_atoms", []))
        nonlocal_atoms = modified_atoms.difference(ALLOWED_LOCAL_ATOMS)
        _require(not nonlocal_atoms, "nonlocal atoms modified: " + ",".join(sorted(nonlocal_atoms)))

    return {
        "technical_status": "PASS",
        "parameter_status": status,
        "classical_mm_ready": status == PASS_STATUS,
        "parent_count": len(parents),
        "atom_count_delta": invariants["atom_count_delta"],
        "system_charge_delta_e": invariants["system_charge_delta_e"],
    }


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main(argv: Optional[Iterable[str]] = None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("authority", type=Path)
    parser.add_argument("--require-ready", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)

    try:
        data = json.loads(args.authority.read_text(encoding="utf-8"))
        result = audit_authority(data, require_ready=args.require_ready)
        result["authority_path"] = str(args.authority.resolve())
        result["authority_sha256"] = sha256_file(args.authority)
    except (OSError, json.JSONDecodeError, AuditError) as exc:
        print(f"FAIL: {exc}", file=sys.stderr)
        return 2

    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
