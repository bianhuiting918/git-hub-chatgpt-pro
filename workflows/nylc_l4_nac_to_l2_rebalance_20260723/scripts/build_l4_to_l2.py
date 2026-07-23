#!/usr/bin/env python3
"""Bonded-graph PA66-L4 to audited PA66-L2 mapping.

The graph matcher never relies on source atom names.  The only fixed target
fact is that atom 33 of the audited L2 topology is the newly generated second
terminal carboxylate oxygen.
"""

from __future__ import annotations

from dataclasses import dataclass
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
    bonds: List[Tuple[int, int]] = []
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
            bonds.append((int(fields[0]), int(fields[1])))
    if not molecule_type or not atoms or not bonds:
        raise ValueError(f"Incomplete ITP topology: {path}")
    neighbors = {atom_id: set() for atom_id in atoms}
    for left, right in bonds:
        if left not in atoms or right not in atoms:
            raise ValueError(f"Bond references missing atom: {left}-{right}")
        neighbors[left].add(right)
        neighbors[right].add(left)
    return Topology(molecule_type=molecule_type, atoms=atoms, neighbors=neighbors)


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
