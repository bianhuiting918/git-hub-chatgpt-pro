#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


ROOT = Path("/data/bht/project01_baker_serhyd_routeB_20260701")
DEFAULT_INPUT = ROOT / "inputs/baker_serhyd/design_pipeline/01_diffusion/inputs/theozyme.pdb"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-pdb", default=str(DEFAULT_INPUT))
    parser.add_argument("--ligand", default="bn1")
    parser.add_argument("--pocket-cutoff", type=float, default=4.0)
    parser.add_argument("--out-json", default=str(ROOT / "manifests/pocket4_first_layered_route_manifest.json"))
    parser.add_argument("--out-tsv", default=str(ROOT / "manifests/pocket4_first_layered_fixed_residues.tsv"))
    parser.add_argument("--contig-tsv", default=str(ROOT / "manifests/pocket4_first_layered_contigs.tsv"))
    return parser.parse_args()


def xyz(line: str) -> tuple[float, float, float]:
    return float(line[30:38]), float(line[38:46]), float(line[46:54])


def distance(a: tuple[float, float, float], b: tuple[float, float, float]) -> float:
    return math.sqrt(sum((a[i] - b[i]) ** 2 for i in range(3)))


def read_structure(path: Path, ligand: str) -> tuple[list[dict[str, object]], list[dict[str, object]]]:
    ligand_l = ligand.lower()
    protein: list[dict[str, object]] = []
    lig_atoms: list[dict[str, object]] = []
    seen_res: set[tuple[str, int]] = set()
    for line in path.read_text(errors="ignore").splitlines():
        if len(line) < 54:
            continue
        rec = line[:6].strip()
        atom = line[12:16].strip()
        resname = line[17:20].strip()
        chain = line[21].strip() or "A"
        if atom.startswith("H"):
            continue
        try:
            resid = int(line[22:26])
            coord = xyz(line)
        except ValueError:
            continue
        if rec == "HETATM" and resname.lower() == ligand_l:
            lig_atoms.append({"chain": chain, "resid": resid, "atom": atom, "coord": coord})
        elif rec == "ATOM":
            key = (chain, resid)
            protein.append(
                {
                    "chain": chain,
                    "resid": resid,
                    "atom": atom,
                    "resname": resname,
                    "coord": coord,
                    "first_atom_for_residue": key not in seen_res,
                }
            )
            seen_res.add(key)
    return protein, lig_atoms


def pocket_residues(
    protein: list[dict[str, object]],
    lig_atoms: list[dict[str, object]],
    cutoff: float,
) -> dict[tuple[str, int], dict[str, object]]:
    if not lig_atoms:
        raise ValueError("No ligand atoms found")
    pocket: dict[tuple[str, int], dict[str, object]] = {}
    for atom in protein:
        coord = atom["coord"]
        assert isinstance(coord, tuple)
        nearest = min(distance(coord, lig["coord"]) for lig in lig_atoms)  # type: ignore[arg-type]
        if nearest <= cutoff:
            key = (str(atom["chain"]), int(atom["resid"]))
            previous = pocket.get(key)
            if previous is None or nearest < float(previous["min_heavy_distance_A"]):
                pocket[key] = {
                    "chain": key[0],
                    "resid": key[1],
                    "resname": atom["resname"],
                    "min_heavy_distance_A": nearest,
                }
    return pocket


def grouped_segments(keys: set[tuple[str, int]]) -> list[tuple[str, int, int]]:
    segments: list[tuple[str, int, int]] = []
    for chain in sorted({chain for chain, _resid in keys}):
        vals = sorted(resid for c, resid in keys if c == chain)
        start = prev = None
        for resid in vals:
            if start is None:
                start = prev = resid
            elif prev is not None and resid == prev + 1:
                prev = resid
            else:
                assert start is not None and prev is not None
                segments.append((chain, start, prev))
                start = prev = resid
        if start is not None and prev is not None:
            segments.append((chain, start, prev))
    return segments


def generate_contig(segments: list[tuple[str, int, int]], flank: int, gap_cap: int | None = None) -> str:
    tokens: list[str] = []
    if flank > 0:
        tokens.append(str(flank))
    prev_chain: str | None = None
    prev_end: int | None = None
    for chain, start, end in segments:
        if prev_chain is not None:
            if chain == prev_chain and prev_end is not None:
                gap = max(0, start - prev_end - 1)
            else:
                gap = flank
            if gap_cap is not None:
                gap = min(gap, gap_cap)
            if gap > 0:
                tokens.append(str(gap))
        tokens.append(f"{chain}{start}-{end}")
        prev_chain = chain
        prev_end = end
    if flank > 0:
        tokens.append(str(flank))
    return ",".join(tokens)


def output_map_from_contig(contig: str) -> list[dict[str, object]]:
    out_pos = 1
    mapping: list[dict[str, object]] = []
    for token in [part.strip() for part in contig.split(",") if part.strip()]:
        if token.isdigit():
            out_pos += int(token)
            continue
        chain = token[0]
        start_s, end_s = token[1:].split("-", 1)
        start, end = int(start_s), int(end_s)
        step = 1 if end >= start else -1
        for resid in range(start, end + step, step):
            mapping.append({"reference_chain": chain, "reference_resid": resid, "output_chain": "A", "output_resid": out_pos})
            out_pos += 1
    return mapping


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    args = parse_args()
    input_pdb = Path(args.input_pdb)
    out_json = Path(args.out_json)
    out_tsv = Path(args.out_tsv)
    contig_tsv = Path(args.contig_tsv)

    protein, lig_atoms = read_structure(input_pdb, args.ligand)
    pocket = pocket_residues(protein, lig_atoms, args.pocket_cutoff)
    fixed_keys = set(pocket)
    segments = grouped_segments(fixed_keys)

    contig_specs = [
        {"contig_set": "pocket4_compact", "flank": 6, "gap_cap": 18},
        {"contig_set": "pocket4_medium", "flank": 8, "gap_cap": 28},
        {"contig_set": "pocket4_expanded", "flank": 10, "gap_cap": None},
    ]
    contig_rows: list[dict[str, object]] = []
    for spec in contig_specs:
        contig = generate_contig(segments, flank=int(spec["flank"]), gap_cap=spec["gap_cap"])
        mapping = output_map_from_contig(contig)
        contig_rows.append(
            {
                "contig_set": spec["contig_set"],
                "contig": contig,
                "fixed_residue_count": len(fixed_keys),
                "fixed_segment_count": len(segments),
                "output_fixed_residue_count": len(mapping),
                "output_length": max(int(row["output_resid"]) for row in mapping) if mapping else 0,
                "flank": spec["flank"],
                "gap_cap": "" if spec["gap_cap"] is None else spec["gap_cap"],
            }
        )

    fixed_rows = [
        {
            "chain": chain,
            "resid": resid,
            "resname": pocket[(chain, resid)]["resname"],
            "min_heavy_distance_A": f"{float(pocket[(chain, resid)]['min_heavy_distance_A']):.3f}",
        }
        for chain, resid in sorted(fixed_keys, key=lambda item: (item[0], item[1]))
    ]
    write_tsv(out_tsv, fixed_rows, ["chain", "resid", "resname", "min_heavy_distance_A"])
    write_tsv(contig_tsv, contig_rows, ["contig_set", "contig", "fixed_residue_count", "fixed_segment_count", "output_fixed_residue_count", "output_length", "flank", "gap_cap"])

    manifest = {
        "status": "DONE",
        "route": "pocket4_first_layered_ca_rfdiffusion",
        "input_pdb": str(input_pdb),
        "ligand": args.ligand,
        "ligand_heavy_atom_count": len(lig_atoms),
        "pocket_cutoff_A": args.pocket_cutoff,
        "fixed_definition": "Protein residues with any heavy atom <= 4 A from ligand heavy atoms; this is the route-defining fixed pocket.",
        "fixed_residue_count": len(fixed_rows),
        "fixed_segments": [{"chain": chain, "start": start, "end": end} for chain, start, end in segments],
        "fixed_residues_tsv": str(out_tsv),
        "contigs_tsv": str(contig_tsv),
        "contigs": contig_rows,
        "mapping_by_contig": {row["contig_set"]: output_map_from_contig(str(row["contig"])) for row in contig_rows},
        "obsolete_route_note": "Earlier motif-only compact outputs are diagnostic only and are not pocket4-first route evidence.",
    }
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(manifest, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
