#!/usr/bin/env python3
"""Bonded-graph PA66-L4 to audited PA66-L2 mapping.

The graph matcher never relies on source atom names.  The only fixed target
fact is that atom 33 of the audited L2 topology is the newly generated second
terminal carboxylate oxygen.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Dict, Iterable, List, Set, Tuple


@dataclass(frozen=True)
class Atom:
    atom_id: int
    atom_type: str
    atom_name: str
    charge: float
    mass: float
    element: str


@dataclass
class Topology:
    molecule_type: str
    atoms: Dict[int, Atom]
    neighbors: Dict[int, Set[int]]
    bond_lengths: Dict[frozenset, float]


@dataclass(frozen=True)
class CutEdge:
    target_endpoint: int
    source_endpoint: int
    removed_source_neighbor: int


@dataclass
class MappingResult:
    target_to_source: Dict[int, int]
    cut_edges: List[CutEdge]
    target_reactive_triplet: Tuple[int, int, int]


def _element(atom_type: str, mass: float) -> str:
    token = atom_type.lower()
    if token.startswith("h") or mass < 2.0:
        return "H"
    if token.startswith("n"):
        return "N"
    if token.startswith("o"):
        return "O"
    if token.startswith("c"):
        return "C"
    raise ValueError(f"Unsupported atom type {atom_type!r} with mass {mass}")


def parse_itp(path: Path) -> Topology:
    section = ""
    molecule_type = ""
    atoms: Dict[int, Atom] = {}
    bonds: List[Tuple[int, int, float]] = []
    for raw in path.read_text(encoding="utf-8").splitlines():
        body = raw.split(";", 1)[0].strip()
        if not body:
            continue
        if body.startswith("[") and body.endswith("]"):
            section = body[1:-1].strip().lower()
            continue
        fields = body.split()
        if section == "moleculetype" and not molecule_type:
            molecule_type = fields[0]
        elif section == "atoms":
            if len(fields) < 8:
                raise ValueError(f"Malformed atoms row in {path}: {raw}")
            atom_id = int(fields[0])
            atom_type = fields[1]
            atom_name = fields[4]
            charge = float(fields[6])
            mass = float(fields[7])
            atoms[atom_id] = Atom(
                atom_id=atom_id,
                atom_type=atom_type,
                atom_name=atom_name,
                charge=charge,
                mass=mass,
                element=_element(atom_type, mass),
            )
        elif section == "bonds":
            if len(fields) < 2:
                raise ValueError(f"Malformed bonds row in {path}: {raw}")
            if len(fields) < 4:\n                raise ValueError(f"Bond lacks equilibrium length in {path}: {raw}")\n            bonds.append((int(fields[0]), int(fields[1]), float(fields[3])))
    if not molecule_type or not atoms or not bonds:
        raise ValueError(f"Incomplete ITP topology: {path}")
    neighbors = {atom_id: set() for atom_id in atoms}
    bond_lengths: Dict[frozenset, float] = {}\n    for left, right, length in bonds:\n        if left not in atoms or right not in atoms:
            raise ValueError(f"Bond references missing atom: {left}-{right}")
        neighbors[left].add(right)
        neighbors[right].add(left)\n        bond_lengths[frozenset((left, right))] = length\n    return Topology(\n        molecule_type=molecule_type, atoms=atoms, neighbors=neighbors,\n        bond_lengths=bond_lengths,\n    )


def validate_reactive_triplet(
    topology: Topology, carbonyl_c: int, carbonyl_o: int, amide_n: int
) -> None:
    expected = (
        topology.atoms[carbonyl_c].element,
        topology.atoms[carbonyl_o].element,
        topology.atoms[amide_n].element,
    )
    if expected != ("C", "O", "N"):
        raise ValueError(f"Reactive triplet elements are {expected}, expected C/O/N")
    if carbonyl_o not in topology.neighbors[carbonyl_c]:
        raise ValueError("Reactive carbonyl C and O are not bonded")
    if amide_n not in topology.neighbors[carbonyl_c]:
        raise ValueError("Reactive carbonyl C and amide N are not bonded")
    if topology.atoms[carbonyl_c].atom_type.lower() != "c":
        raise ValueError("Reactive carbon is not a carbonyl atom type")


def _heavy_neighbors(topology: Topology, atom_id: int) -> Set[int]:
    return {
        neighbor
        for neighbor in topology.neighbors[atom_id]
        if topology.atoms[neighbor].element != "H"
    }


def _compatible(source: Atom, target: Atom, target_id: int) -> bool:
    if source.element != target.element:
        return False
    if target_id == 1:
        return source.element == "N"
    return source.atom_type.lower() == target.atom_type.lower()


def _target_amides(target: Topology, inherited: Set[int]) -> Iterable[Tuple[int, int, int]]:
    for carbon_id in sorted(inherited):
        atom = target.atoms[carbon_id]
        if atom.element != "C" or atom.atom_type.lower() != "c":
            continue
        oxygens = [
            n for n in target.neighbors[carbon_id]
            if n in inherited and target.atoms[n].element == "O"
        ]
        nitrogens = [
            n for n in target.neighbors[carbon_id]
            if n in inherited and target.atoms[n].element == "N"
        ]
        if len(oxygens) == 1 and len(nitrogens) == 1:
            yield carbon_id, oxygens[0], nitrogens[0]


def _search_one(
    source: Topology,
    target: Topology,
    inherited: Set[int],
    anchor: Dict[int, int],
) -> List[Dict[int, int]]:
    target_heavy_neighbors = {
        atom_id: _heavy_neighbors(target, atom_id) & inherited for atom_id in inherited
    }
    source_heavy_neighbors = {
        atom_id: _heavy_neighbors(source, atom_id)
        for atom_id, atom in source.atoms.items()
        if atom.element != "H"
    }
    results: List[Dict[int, int]] = []

    def degree_allowed(target_id: int, source_id: int) -> bool:
        target_degree = len(target_heavy_neighbors[target_id])
        source_degree = len(source_heavy_neighbors[source_id])
        if target_id in {1, 31}:
            return source_degree == target_degree + 1
        return source_degree == target_degree

    def consistent(target_id: int, source_id: int, mapping: Dict[int, int]) -> bool:
        if not _compatible(source.atoms[source_id], target.atoms[target_id], target_id):
            return False
        if not degree_allowed(target_id, source_id):
            return False
        for mapped_target, mapped_source in mapping.items():
            target_edge = mapped_target in target_heavy_neighbors[target_id]
            source_edge = mapped_source in source_heavy_neighbors[source_id]
            if target_edge != source_edge:
                return False
        return True

    if len(set(anchor.values())) != len(anchor):
        return results
    for target_id, source_id in anchor.items():
        if not consistent(target_id, source_id, {
            t: s for t, s in anchor.items() if t != target_id
        }):
            return results

    def recurse(mapping: Dict[int, int]) -> None:
        if len(mapping) == len(inherited):
            results.append(dict(mapping))
            return
        used = set(mapping.values())
        frontier = [
            target_id
            for target_id in inherited - mapping.keys()
            if target_heavy_neighbors[target_id] & mapping.keys()
        ]
        if not frontier:
            return

        ranked = []
        for target_id in frontier:
            mapped_neighbors = target_heavy_neighbors[target_id] & mapping.keys()
            candidate_sets = [
                source_heavy_neighbors[mapping[neighbor]] for neighbor in mapped_neighbors
            ]
            candidates = set.intersection(*candidate_sets) - used
            candidates = {
                source_id
                for source_id in candidates
                if consistent(target_id, source_id, mapping)
            }
            ranked.append((len(candidates), target_id, sorted(candidates)))
        _, target_id, candidates = min(ranked)
        for source_id in candidates:
            mapping[target_id] = source_id
            recurse(mapping)
            del mapping[target_id]

    recurse(dict(anchor))
    return results


def discover_heavy_mapping(
    source: Topology,
    target: Topology,
    source_reactive_triplet: Tuple[int, int, int],
) -> MappingResult:
    validate_reactive_triplet(source, *source_reactive_triplet)
    if 33 not in target.atoms or target.atoms[33].element != "O":
        raise ValueError("Audited L2 atom 33 must be the generated terminal O5")
    inherited = {
        atom_id for atom_id, atom in target.atoms.items()
        if atom.element != "H" and atom_id != 33
    }
    source_c, source_o, source_n = source_reactive_triplet
    candidates = []
    target_preference = {25: 0, 15: 1, 9: 2}
    for target_triplet in _target_amides(target, inherited):
        target_c, target_o, target_n = target_triplet
        anchor = {
            target_c: source_c,
            target_o: source_o,
            target_n: source_n,
        }
        for mapping in _search_one(source, target, inherited, anchor):
            mapped_sources = set(mapping.values())
            cuts: List[CutEdge] = []
            valid = True
            for endpoint in (1, 31):
                source_endpoint = mapping[endpoint]
                extras = sorted(_heavy_neighbors(source, source_endpoint) - mapped_sources)
                if len(extras) != 1:
                    valid = False
                    break
                cuts.append(CutEdge(endpoint, source_endpoint, extras[0]))
            if not valid:
                continue
            score = (
                target_preference.get(target_c, 99),
                tuple(mapping[target_id] for target_id in sorted(mapping)),
            )
            candidates.append((score, mapping, cuts, target_triplet))
    if not candidates:
        raise ValueError(
            f"No two-cut audited L2 subgraph contains reactive triplet "
            f"{source_reactive_triplet}"
        )
    _, mapping, cuts, target_triplet = min(candidates, key=lambda item: item[0])
    return MappingResult(
        target_to_source=dict(sorted(mapping.items())),
        cut_edges=cuts,
        target_reactive_triplet=target_triplet,
    )



Vector = Tuple[float, float, float]


def _vadd(left: Vector, right: Vector) -> Vector:
    return tuple(a + b for a, b in zip(left, right))  # type: ignore[return-value]


def _vsub(left: Vector, right: Vector) -> Vector:
    return tuple(a - b for a, b in zip(left, right))  # type: ignore[return-value]


def _vscale(vector: Vector, factor: float) -> Vector:
    return tuple(value * factor for value in vector)  # type: ignore[return-value]


def _dot(left: Vector, right: Vector) -> float:
    return sum(a * b for a, b in zip(left, right))


def _cross(left: Vector, right: Vector) -> Vector:
    return (
        left[1] * right[2] - left[2] * right[1],
        left[2] * right[0] - left[0] * right[2],
        left[0] * right[1] - left[1] * right[0],
    )


def _norm(vector: Vector) -> float:
    return math.sqrt(_dot(vector, vector))


def _unit(vector: Vector) -> Vector:
    length = _norm(vector)
    if length < 1e-10:
        raise ValueError("Cannot normalize a zero-length vector")
    return _vscale(vector, 1.0 / length)


def build_atom_mapping(
    source: Topology, target: Topology, heavy_mapping: MappingResult
) -> Dict[int, int | None]:
    """Map equivalent atoms by bonded parent; reserve rebuilt endpoint atoms."""
    mapping: Dict[int, int | None] = dict(heavy_mapping.target_to_source)
    mapping[33] = None
    for target_parent, source_parent in sorted(heavy_mapping.target_to_source.items()):
        target_hydrogens = sorted(
            neighbor
            for neighbor in target.neighbors[target_parent]
            if target.atoms[neighbor].element == "H"
        )
        if not target_hydrogens:
            continue
        if target_parent == 1:
            for target_h in target_hydrogens:
                mapping[target_h] = None
            continue
        source_hydrogens = sorted(
            neighbor
            for neighbor in source.neighbors[source_parent]
            if source.atoms[neighbor].element == "H"
        )
        if len(target_hydrogens) != len(source_hydrogens):
            raise ValueError(
                f"Hydrogen-count mismatch for target parent {target_parent} "
                f"and source parent {source_parent}: "
                f"{len(target_hydrogens)} != {len(source_hydrogens)}"
            )
        for target_h, source_h in zip(target_hydrogens, source_hydrogens):
            mapping[target_h] = source_h
    if set(mapping) != set(target.atoms):
        missing = sorted(set(target.atoms) - set(mapping))
        raise ValueError(f"Incomplete atom mapping; missing targets {missing}")
    mapped_sources = [source_id for source_id in mapping.values() if source_id is not None]
    if len(mapped_sources) != len(set(mapped_sources)):
        raise ValueError("Source atom is reused in atom mapping")
    return dict(sorted(mapping.items()))


def build_l2_coordinates(
    source: Topology,
    target: Topology,
    heavy_mapping: MappingResult,
    source_coordinates: Dict[int, Vector],
) -> Tuple[Dict[int, Vector], Dict[int, int | None]]:
    mapping = build_atom_mapping(source, target, heavy_mapping)
    coordinates: Dict[int, Vector] = {
        target_id: source_coordinates[source_id]
        for target_id, source_id in mapping.items()
        if source_id is not None
    }

    cuts = {cut.target_endpoint: cut for cut in heavy_mapping.cut_edges}
    if set(cuts) != {1, 31}:
        raise ValueError("Expected one source cut at each L2 endpoint")

    # Generate terminal O5 along the source C--removed-amide-N direction.
    carboxyl_c = coordinates[31]
    removed_n = source_coordinates[cuts[31].removed_source_neighbor]
    c_o_length = target.bond_lengths[frozenset((31, 33))]
    coordinates[33] = _vadd(
        carboxyl_c, _vscale(_unit(_vsub(removed_n, carboxyl_c)), c_o_length)
    )

    # Rebuild three ammonium hydrogens as an ideal tetrahedron about retained N--C.
    terminal_n = coordinates[1]
    retained_c = coordinates[2]
    axis = _unit(_vsub(retained_c, terminal_n))
    source_n = heavy_mapping.target_to_source[1]
    source_hydrogens = sorted(
        neighbor
        for neighbor in source.neighbors[source_n]
        if source.atoms[neighbor].element == "H"
    )
    if source_hydrogens:
        reference = _vsub(source_coordinates[source_hydrogens[0]], terminal_n)
    else:
        reference = _vsub(
            source_coordinates[cuts[1].removed_source_neighbor], terminal_n
        )
    projected = _vsub(reference, _vscale(axis, _dot(reference, axis)))
    if _norm(projected) < 1e-8:
        fallback = (1.0, 0.0, 0.0) if abs(axis[0]) < 0.8 else (0.0, 1.0, 0.0)
        projected = _cross(axis, fallback)
    basis_u = _unit(projected)
    basis_v = _unit(_cross(axis, basis_u))
    cosine = -1.0 / 3.0
    sine = math.sqrt(8.0 / 9.0)
    for index, target_h in enumerate((34, 35, 36)):
        phi = index * 2.0 * math.pi / 3.0
        radial = _vadd(
            _vscale(basis_u, math.cos(phi)),
            _vscale(basis_v, math.sin(phi)),
        )
        direction = _vadd(_vscale(axis, cosine), _vscale(radial, sine))
        n_h_length = target.bond_lengths[frozenset((1, target_h))]
        coordinates[target_h] = _vadd(
            terminal_n, _vscale(_unit(direction), n_h_length)
        )
    if set(coordinates) != set(target.atoms):
        missing = sorted(set(target.atoms) - set(coordinates))
        raise ValueError(f"Incomplete target coordinates; missing {missing}")
    return dict(sorted(coordinates.items())), mapping
