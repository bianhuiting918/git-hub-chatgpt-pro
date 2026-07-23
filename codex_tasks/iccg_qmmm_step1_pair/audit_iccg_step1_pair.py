#!/usr/bin/env python3
"""Hard-gate audits for the ICCG Step1 LG1/LG2 paired pilot."""
from __future__ import annotations
import argparse, json, math
from pathlib import Path
from typing import Sequence
from build_iccg_step1_pair import Atom, ligand_atoms, protein_atoms

VDW = {"H": 1.20, "C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80}

def _dist(a: Atom, b: Atom) -> float:
    return math.sqrt(sum((a.xyz[i] - b.xyz[i]) ** 2 for i in range(3)))

def _pair(i: int, j: int) -> tuple[int, int]:
    return (i, j) if i < j else (j, i)

def topology_pairing_pass(state1_order: Sequence[str], state2_order: Sequence[str]) -> bool:
    return list(state1_order) == list(state2_order) and len(state1_order) > 0

def geometry_gate(atoms: Sequence[Atom], bonds: set[tuple[int, int]] | None = None, link_pairs: set[tuple[int, int]] | None = None) -> dict[str, object]:
    """Audit all protein-ligand nonbonded heavy-atom pairs before PASS/FAIL."""
    bonds = {_pair(*p) for p in (bonds or set())}
    link_pairs = {_pair(*p) for p in (link_pairs or set())}
    proteins = [a for a in protein_atoms(atoms) if a.element.upper() != "H"]
    ligands = [a for a in ligand_atoms(atoms) if a.element.upper() != "H"]
    if not proteins or not ligands:
        return {"pass": False, "reason": "NOT_SUBMITTED_MISSING_PROTEIN_OR_LIGAND"}
    max_overlap = -float("inf")
    min_nonbonded = float("inf")
    worst: dict[str, object] | None = None
    for a in proteins:
        if any(not math.isfinite(c) for c in a.xyz):
            return {"pass": False, "reason": "NOT_SUBMITTED_NAN_INF"}
        for b in ligands:
            if any(not math.isfinite(c) for c in b.xyz):
                return {"pass": False, "reason": "NOT_SUBMITTED_NAN_INF"}
            if _pair(a.index, b.index) in bonds or _pair(a.index, b.index) in link_pairs:
                continue
            d = _dist(a, b)
            overlap = VDW.get(a.element.upper(), 1.7) + VDW.get(b.element.upper(), 1.7) - d
            if d < min_nonbonded:
                min_nonbonded = d
            if overlap > max_overlap:
                max_overlap = overlap
                worst = {"protein_atom": a.name, "protein_residue": a.resid, "protein_element": a.element, "ligand_atom": b.name, "ligand_resname": b.resname, "ligand_element": b.element, "distance_A": d, "overlap_A": overlap, "atom_indices": [a.index, b.index]}
    passed = min_nonbonded >= 1.20 and max_overlap <= 0.80
    return {"pass": passed, "reason": "PASS" if passed else "FAIL_GEOMETRY_CLASH_NOT_LABEL", "min_nonbonded_A": min_nonbonded, "max_vdw_overlap_A": max_overlap, "worst_pair": worst}

def can_submit(preflight_report: Path) -> bool:
    """Return Stage-A structural readiness, not topology/submission readiness."""
    report = json.loads(preflight_report.read_text())
    gates = report.get("gates")
    if not isinstance(gates, list) or not gates:
        return False
    by_name = {g.get("name"): g for g in gates if isinstance(g, dict)}
    required = {"active_iccg_258_ser165_og", "lg1_54_32_two_rings", "lg2_54_32_two_rings", "lg1_lg2_atom_name_order", "paired_protein_coordinates_identical", "LG1_protein_ligand_geometry", "LG2_protein_ligand_geometry", "ile243_structural_mapping"}
    return all(by_name.get(name, {}).get("pass") is True for name in required)

def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("preflight_report", type=Path)
    p.add_argument("--submit-check", action="store_true")
    args = p.parse_args(argv)
    if args.submit_check:
        ok = can_submit(args.preflight_report)
        print("PASS" if ok else "NOT_SUBMITTED_HARD_GATE")
        return 0 if ok else 2
    print(json.dumps(json.loads(args.preflight_report.read_text()), indent=2, sort_keys=True))
    return 0
if __name__ == "__main__":
    raise SystemExit(main())
