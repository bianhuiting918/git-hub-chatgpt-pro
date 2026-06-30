#!/usr/bin/env python3
import argparse
import json
import re
from collections import defaultdict
from pathlib import Path


def parse_pdb(path: Path):
    residues = defaultdict(set)
    ligand_counts = defaultdict(int)
    atom_counts = 0
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith(("ATOM", "HETATM")):
            continue
        atom_counts += 1
        atom = line[12:16].strip()
        resn = line[17:20].strip()
        chain = line[21].strip() or "_"
        try:
            resi = int(line[22:26])
        except ValueError:
            continue
        residues[(chain, resi, resn)].add(atom)
        if line.startswith("HETATM"):
            ligand_counts[resn] += 1
    return residues, dict(ligand_counts), atom_counts


def parse_contig(contig: str):
    refs = []
    for token in contig.split(","):
        token = token.strip().strip("'\"")
        match = re.fullmatch(r"([A-Za-z])(\d+)-(\d+)", token)
        if not match:
            continue
        chain, start, end = match.group(1), int(match.group(2)), int(match.group(3))
        refs.extend((chain, i) for i in range(start, end + 1))
    return refs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default="/data/bht/project01_baker_serhyd_routeB_20260701")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    root = Path(args.root)
    base = root / "inputs" / "baker_serhyd"
    theozyme_pdb = base / "design_pipeline" / "01_diffusion" / "inputs" / "theozyme.pdb"
    tasks_json = base / "design_pipeline" / "01_diffusion" / "inputs" / "tasks.json"
    theozyme_cst = base / "design_pipeline" / "03_design" / "inputs" / "theozyme.cst"
    bn1_params = base / "design_pipeline" / "03_design" / "inputs" / "bn1.params"
    simple_pdb = base / "motif_gen" / "02_diffusion" / "inputs" / "simple_theozyme.pdb"
    mu1_params = base / "motif_gen" / "01_sampling_his_stub" / "inputs" / "mu1.params"
    mu1_cst = base / "motif_gen" / "01_sampling_his_stub" / "inputs" / "1LNS_mu1.cst"

    required = [theozyme_pdb, tasks_json, theozyme_cst, bn1_params, simple_pdb, mu1_params, mu1_cst]
    missing = [str(p) for p in required if not p.is_file() or p.stat().st_size == 0]

    residues, ligand_counts, atom_counts = parse_pdb(theozyme_pdb) if theozyme_pdb.is_file() else ({}, {}, 0)
    tasks = json.loads(tasks_json.read_text()) if tasks_json.is_file() else []
    contig = tasks[0].get("contigmap", {}).get("contigs", [""])[0] if tasks else ""
    contig_refs = parse_contig(contig)
    residue_index = {(c, i): resn for (c, i, resn) in residues}
    missing_contig_refs = [f"{c}{i}" for c, i in contig_refs if (c, i) not in residue_index]
    contig_residues = [{"chain": c, "resi": i, "resn": residue_index.get((c, i), "MISSING")} for c, i in contig_refs]

    cst_text = theozyme_cst.read_text(errors="replace") if theozyme_cst.is_file() else ""
    params_text = bn1_params.read_text(errors="replace") if bn1_params.is_file() else ""
    constraint_checks = {
        "has_ser_attack": "serine attack" in cst_text,
        "has_oxyanion_hole": "oxyanion hole" in cst_text,
        "has_ser_his": "ser-his" in cst_text,
        "has_his_asp": "His-Asp" in cst_text,
        "bn1_params_has_c1_o1_o2": all(f"ATOM  {x}" in params_text for x in ["C1", "O1", "O2"]),
    }

    status = "PASS" if not missing and not missing_contig_refs and ligand_counts.get("bn1", 0) > 0 and all(constraint_checks.values()) else "FAIL"
    out = {
        "status": status,
        "route_label": "baker_theozyme_new_backbone",
        "root": str(root),
        "checked_files": {p.name: str(p) for p in required},
        "missing_or_empty_files": missing,
        "theozyme_atom_records": atom_counts,
        "ligand_counts": ligand_counts,
        "tasks_contig": contig,
        "contig_reference_count": len(contig_refs),
        "missing_contig_refs": missing_contig_refs,
        "contig_residues": contig_residues,
        "constraint_checks": constraint_checks,
        "evaluated_universe_note": "This validates inputs only; no new backbone candidate is evaluated until RFAA writes output PDBs.",
    }

    out_path = Path(args.out) if args.out else root / "manifests" / "baker_theozyme_input_preflight.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2) + "\n")
    print(out_path)
    if status != "PASS":
        raise SystemExit(1)


if __name__ == "__main__":
    main()
