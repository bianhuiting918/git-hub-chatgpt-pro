#!/usr/bin/env python3
"""Build the capped NylC Thr267 A1 model used for local parameter generation."""

from __future__ import annotations

import math
from collections import namedtuple
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

Vector = Tuple[float, float, float]
Atom = namedtuple(
    "Atom",
    "name element sybyl_type formal_charge coordinate_nm source_atom_id source_role",
)
Bond = namedtuple("Bond", "atom1 atom2 order")


class ModelError(ValueError):
    """Raised when the parent graph cannot define the required A1 model."""


class Model:
    def __init__(self, atoms: Sequence[Atom], bonds: Sequence[Bond]):
        self.atoms = tuple(atoms)
        self.bonds = tuple(bonds)
        self._atoms = {atom.name: atom for atom in self.atoms}
        if len(self._atoms) != len(self.atoms):
            raise ModelError("model atom names are not unique")

    def atom(self, name: str) -> Atom:
        return self._atoms[name]

    def has_bond(self, atom1: str, atom2: str) -> bool:
        key = frozenset((atom1, atom2))
        return any(frozenset((bond.atom1, bond.atom2)) == key for bond in self.bonds)

    @property
    def atom_map(self) -> List[Dict[str, Any]]:
        return [
            {
                "model_atom": atom.name,
                "source_atom_id": atom.source_atom_id,
                "source_role": atom.source_role,
            }
            for atom in self.atoms
        ]


def _sub(a: Vector, b: Vector) -> Vector:
    return tuple(x - y for x, y in zip(a, b))  # type: ignore[return-value]


def _add(a: Vector, b: Vector) -> Vector:
    return tuple(x + y for x, y in zip(a, b))  # type: ignore[return-value]


def _scale(v: Vector, factor: float) -> Vector:
    return tuple(factor * x for x in v)  # type: ignore[return-value]


def _dot(a: Vector, b: Vector) -> float:
    return sum(x * y for x, y in zip(a, b))


def _cross(a: Vector, b: Vector) -> Vector:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(v: Vector) -> float:
    return math.sqrt(_dot(v, v))


def _unit(v: Vector) -> Vector:
    length = _norm(v)
    if length < 1.0e-10:
        raise ModelError("cannot normalize a zero-length geometry vector")
    return _scale(v, 1.0 / length)


def _missing_tetrahedral_direction(center: Vector, neighbors: Iterable[Vector]) -> Vector:
    vectors = [_unit(_sub(point, center)) for point in neighbors]
    summed = (
        sum(v[0] for v in vectors),
        sum(v[1] for v in vectors),
        sum(v[2] for v in vectors),
    )
    if _norm(summed) < 1.0e-6 and len(vectors) >= 2:
        summed = _cross(vectors[0], vectors[1])
    return _unit(_scale(summed, -1.0))


def _cap_hydrogen_directions(
    center: Vector,
    bonded_n: Vector,
    existing_h: Vector,
) -> Tuple[Vector, Vector]:
    toward_n = _unit(_sub(bonded_n, center))
    toward_h = _unit(_sub(existing_h, center))
    normal = _unit(_cross(toward_n, toward_h))
    bisector = _unit(_scale(_add(toward_n, toward_h), -1.0))
    return (
        _unit(_add(_scale(bisector, 0.55), _scale(normal, 0.835))),
        _unit(_add(_scale(bisector, 0.55), _scale(normal, -0.835))),
    )


def _record_coordinate(record: Mapping[str, Any]) -> Vector:
    xyz = tuple(float(value) for value in record["coordinate_nm"])
    if len(xyz) != 3:
        raise ModelError("coordinate must have three components")
    return xyz  # type: ignore[return-value]


def build_a1_capped_model(
    thr_atoms: Mapping[str, Mapping[str, Any]],
    next_atoms: Mapping[str, Mapping[str, Any]],
    source_bonds: Iterable[frozenset],
) -> Model:
    """Transfer HG1 from Oγ to Nα and cap the following peptide N as NME."""

    required_thr = {
        "N", "H1", "H2", "CA", "HA", "CB", "HB", "CG2",
        "HG21", "HG22", "HG23", "OG1", "HG1", "C", "O",
    }
    required_next = {"N", "H", "CA", "HA"}
    missing_thr = required_thr.difference(thr_atoms)
    missing_next = required_next.difference(next_atoms)
    if missing_thr:
        raise ModelError("missing Thr267 atoms: " + ",".join(sorted(missing_thr)))
    if missing_next:
        raise ModelError("missing residue-268 cap atoms: " + ",".join(sorted(missing_next)))

    bonds = set(source_bonds)
    required_parent_bonds = {
        frozenset(("OG1", "HG1")),
        frozenset(("N", "H1")),
        frozenset(("N", "H2")),
        frozenset(("N", "CA")),
        frozenset(("C", "next:N")),
    }
    missing_bonds = required_parent_bonds.difference(bonds)
    if frozenset(("OG1", "HG1")) in missing_bonds:
        raise ModelError("parent graph lacks the required OG1-HG1 bond")
    if missing_bonds:
        rendered = ["-".join(sorted(row)) for row in missing_bonds]
        raise ModelError("parent graph lacks required bonds: " + ",".join(sorted(rendered)))

    n_xyz = _record_coordinate(thr_atoms["N"])
    hg_direction = _missing_tetrahedral_direction(
        n_xyz,
        [
            _record_coordinate(thr_atoms["CA"]),
            _record_coordinate(thr_atoms["H1"]),
            _record_coordinate(thr_atoms["H2"]),
        ],
    )
    hg_xyz = _add(n_xyz, _scale(hg_direction, 0.101))

    element_and_type = {
        "N": ("N", "N.4", 1),
        "H1": ("H", "H", 0),
        "H2": ("H", "H", 0),
        "HG1": ("H", "H", 0),
        "CA": ("C", "C.3", 0),
        "HA": ("H", "H", 0),
        "CB": ("C", "C.3", 0),
        "HB": ("H", "H", 0),
        "CG2": ("C", "C.3", 0),
        "HG21": ("H", "H", 0),
        "HG22": ("H", "H", 0),
        "HG23": ("H", "H", 0),
        "OG1": ("O", "O.3", -1),
        "C": ("C", "C.2", 0),
        "O": ("O", "O.2", 0),
    }
    atom_order = [
        "N", "H1", "H2", "HG1", "CA", "HA", "CB", "HB",
        "CG2", "HG21", "HG22", "HG23", "OG1", "C", "O",
    ]
    model_atoms: List[Atom] = []
    for name in atom_order:
        element, atom_type, charge = element_and_type[name]
        xyz = hg_xyz if name == "HG1" else _record_coordinate(thr_atoms[name])
        model_atoms.append(
            Atom(
                name,
                element,
                atom_type,
                charge,
                xyz,
                int(thr_atoms[name]["source_atom_id"]),
                "Thr267",
            )
        )

    cap_n = _record_coordinate(next_atoms["N"])
    cap_h = _record_coordinate(next_atoms["H"])
    cap_c = _record_coordinate(next_atoms["CA"])
    cap_h1 = _record_coordinate(next_atoms["HA"])
    d2, d3 = _cap_hydrogen_directions(cap_c, cap_n, cap_h1)
    cap_h2 = _add(cap_c, _scale(d2, 0.109))
    cap_h3 = _add(cap_c, _scale(d3, 0.109))
    model_atoms.extend(
        [
            Atom("CAP_N", "N", "N.am", 0, cap_n, int(next_atoms["N"]["source_atom_id"]), "Cap268"),
            Atom("CAP_H", "H", "H", 0, cap_h, int(next_atoms["H"]["source_atom_id"]), "Cap268"),
            Atom("CAP_C", "C", "C.3", 0, cap_c, int(next_atoms["CA"]["source_atom_id"]), "Cap268"),
            Atom("CAP_H1", "H", "H", 0, cap_h1, int(next_atoms["HA"]["source_atom_id"]), "Cap268"),
            Atom("CAP_H2", "H", "H", 0, cap_h2, None, "GeneratedCap"),
            Atom("CAP_H3", "H", "H", 0, cap_h3, None, "GeneratedCap"),
        ]
    )

    bond_rows = [
        ("N", "H1", "1"), ("N", "H2", "1"), ("N", "HG1", "1"), ("N", "CA", "1"),
        ("CA", "HA", "1"), ("CA", "CB", "1"), ("CA", "C", "1"),
        ("CB", "HB", "1"), ("CB", "CG2", "1"), ("CB", "OG1", "1"),
        ("CG2", "HG21", "1"), ("CG2", "HG22", "1"), ("CG2", "HG23", "1"),
        ("C", "O", "2"), ("C", "CAP_N", "am"),
        ("CAP_N", "CAP_H", "1"), ("CAP_N", "CAP_C", "1"),
        ("CAP_C", "CAP_H1", "1"), ("CAP_C", "CAP_H2", "1"), ("CAP_C", "CAP_H3", "1"),
    ]
    return Model(model_atoms, [Bond(*row) for row in bond_rows])


def render_mol2(model: Model) -> str:
    lines = [
        "@<TRIPOS>MOLECULE",
        "NTA1_CAP",
        f"{len(model.atoms)} {len(model.bonds)} 1 0 0",
        "SMALL",
        "USER_CHARGES",
        "",
        "@<TRIPOS>ATOM",
    ]
    atom_numbers = {atom.name: i for i, atom in enumerate(model.atoms, 1)}
    for i, atom in enumerate(model.atoms, 1):
        x, y, z = (_scale(atom.coordinate_nm, 10.0))
        lines.append(
            f"{i:7d} {atom.name:<8s} {x:10.4f} {y:10.4f} {z:10.4f} "
            f"{atom.sybyl_type:<6s} 1 NTA1 {float(atom.formal_charge):10.6f}"
        )
    lines.append("@<TRIPOS>BOND")
    for i, bond in enumerate(model.bonds, 1):
        lines.append(
            f"{i:6d} {atom_numbers[bond.atom1]:5d} {atom_numbers[bond.atom2]:5d} {bond.order}"
        )
    lines.extend(["@<TRIPOS>SUBSTRUCTURE", "     1 NTA1        1 TEMP              0 ****  ****    0 ROOT"])
    return "\n".join(lines) + "\n"
