#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
import time
from pathlib import Path

AA1_TO_3 = {
    "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP", "C": "CYS",
    "Q": "GLN", "E": "GLU", "G": "GLY", "H": "HIS", "I": "ILE",
    "L": "LEU", "K": "LYS", "M": "MET", "F": "PHE", "P": "PRO",
    "S": "SER", "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL",
}
BACKBONE = {"N", "CA", "C", "O", "OXT"}
FIXED = {13, 14, 15, 16, 17, 54, 55, 56, 72, 73, 74, 148, 149, 150}


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--reference-pdb", required=True)
    p.add_argument("--selected-tsv", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--manifest-tsv", required=True)
    p.add_argument("--status-tsv", required=True)
    p.add_argument("--summary-json", required=True)
    p.add_argument("--bins", nargs="+", default=["90", "80", "70", "60", "50"])
    p.add_argument("--per-bin", type=int, default=2)
    p.add_argument("--python", default=sys.executable)
    p.add_argument("--minimize", action="store_true")
    return p.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def select_rows(rows: list[dict[str, str]], bins: list[str], per_bin: int) -> list[dict[str, str]]:
    counts = {b: 0 for b in bins}
    picked = []
    for row in rows:
        b = row["bin"]
        if b in counts and counts[b] < per_bin:
            picked.append(row)
            counts[b] += 1
        if all(counts[b] >= per_bin for b in counts):
            break
    return picked


def thread_backbone(reference: Path, sequence: str, threaded_pdb: Path) -> None:
    lines = []
    for line in reference.read_text(errors="ignore").splitlines():
        rec = line[:6].strip()
        if rec == "ATOM":
            if len(line) < 54:
                continue
            atom = line[12:16].strip()
            chain = line[21].strip() or "A"
            try:
                resid = int(line[22:26])
            except ValueError:
                continue
            if chain != "A" or resid < 1 or resid > len(sequence):
                continue
            if atom not in BACKBONE and resid not in FIXED:
                continue
            aa = sequence[resid - 1]
            resname = AA1_TO_3.get(aa)
            if not resname:
                raise ValueError(f"unsupported residue {aa!r} at A{resid}")
            lines.append(line[:17] + f"{resname:>3}" + line[20:])
        elif rec in {"TER", "END"}:
            continue
    lines.append("TER")
    threaded_pdb.write_text("\n".join(lines) + "\n")


def ligand_lines(reference: Path) -> list[str]:
    out = []
    for line in reference.read_text(errors="ignore").splitlines():
        if line.startswith("HETATM") and line[17:20].strip().lower() == "bn1":
            out.append(line)
    return out


def append_ligand(protein_pdb: Path, ligand: list[str], final_pdb: Path) -> None:
    body = []
    for line in protein_pdb.read_text(errors="ignore").splitlines():
        if line.startswith("END"):
            continue
        body.append(line)
    body.extend(ligand)
    body.append("END")
    final_pdb.write_text("\n".join(body) + "\n")


def run_pdbfixer(python: str, threaded: Path, fixed: Path, minimize: bool) -> None:
    code = r'''
from pathlib import Path
import sys
from pdbfixer import PDBFixer
from openmm.app import PDBFile
inp, outp, do_min = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3] == "1"
fixer = PDBFixer(filename=str(inp))
fixer.findMissingResidues()
fixer.missingResidues = {}
fixer.findMissingAtoms()
fixer.addMissingAtoms()
fixer.addMissingHydrogens(7.0)
if do_min:
    import openmm as mm
    from openmm import unit
    from openmm.app import ForceField, Modeller, NoCutoff, Simulation
    forcefield = ForceField("amber14-all.xml")
    modeller = Modeller(fixer.topology, fixer.positions)
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=NoCutoff, constraints=None)
    restraint = mm.CustomExternalForce("0.5*k*((x-x0)^2+(y-y0)^2+(z-z0)^2)")
    restraint.addGlobalParameter("k", 1000.0)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")
    fixed_res = {13,14,15,16,17,54,55,56,72,73,74,148,149,150}
    backbone = {"N", "CA", "C", "O", "OXT"}
    for atom in modeller.topology.atoms():
        try:
            resid = int(atom.residue.id)
        except Exception:
            continue
        if atom.name in backbone or resid in fixed_res:
            pos = modeller.positions[atom.index].value_in_unit(unit.nanometer)
            restraint.addParticle(atom.index, pos)
    system.addForce(restraint)
    integrator = mm.LangevinIntegrator(300*unit.kelvin, 1/unit.picosecond, 0.002*unit.picoseconds)
    simulation = Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    simulation.minimizeEnergy(maxIterations=100)
    state = simulation.context.getState(getPositions=True)
    with outp.open("w") as handle:
        PDBFile.writeFile(modeller.topology, state.getPositions(), handle)
else:
    with outp.open("w") as handle:
        PDBFile.writeFile(fixer.topology, fixer.positions, handle)
'''
    subprocess.run([python, "-c", code, str(threaded), str(fixed), "1" if minimize else "0"], check=True)


def main() -> None:
    args = parse_args()
    ref = Path(args.reference_pdb)
    selected = read_tsv(Path(args.selected_tsv))
    picked = select_rows(selected, args.bins, args.per_bin)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    lig = ligand_lines(ref)
    status_rows = []
    manifest_rows = []
    started = time.strftime("%F %T")
    for row in picked:
        sid = row["sample_id"]
        threaded = outdir / f"{sid}.threaded_backbone.pdb"
        fixed = outdir / f"{sid}.pdbfixer{'_min' if args.minimize else ''}.protein.pdb"
        final = outdir / f"{sid}.pdbfixer{'_min' if args.minimize else ''}.bn1.pdb"
        t0 = time.time()
        status = "OK"
        msg = ""
        try:
            thread_backbone(ref, row["sequence"], threaded)
            run_pdbfixer(args.python, threaded, fixed, args.minimize)
            append_ligand(fixed, lig, final)
        except Exception as exc:
            status = "FAIL"
            msg = repr(exc)
        elapsed = round(time.time() - t0, 2)
        out_row = dict(row)
        out_row["postseq_predicted_pdb"] = str(final)
        manifest_rows.append(out_row)
        status_rows.append({
            "sample_id": sid,
            "bin": row["bin"],
            "identity": row["identity"],
            "mutation_count": row["mutation_count"],
            "sequence_length": row["sequence_length"],
            "output_pdb": str(final),
            "status": status,
            "bytes": final.stat().st_size if final.exists() else 0,
            "elapsed_s": elapsed,
            "message": msg,
            "timestamp": time.strftime("%F %T"),
        })
        print(sid, row["bin"], status, elapsed, msg, flush=True)
    manifest_fields = list(picked[0].keys()) if picked else []
    if "postseq_predicted_pdb" not in manifest_fields:
        manifest_fields.append("postseq_predicted_pdb")
    with Path(args.manifest_tsv).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, manifest_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(manifest_rows)
    status_fields = ["sample_id","bin","identity","mutation_count","sequence_length","output_pdb","status","bytes","elapsed_s","message","timestamp"]
    with Path(args.status_tsv).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, status_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(status_rows)
    counts = {}
    bins = {}
    for row in status_rows:
        counts[row["status"]] = counts.get(row["status"], 0) + 1
        bins.setdefault(row["bin"], {})
        bins[row["bin"]][row["status"]] = bins[row["bin"]].get(row["status"], 0) + 1
    summary = {
        "started": started,
        "finished": time.strftime("%F %T"),
        "reference_pdb": str(ref),
        "selected_tsv": args.selected_tsv,
        "outdir": str(outdir),
        "manifest_tsv": args.manifest_tsv,
        "status_tsv": args.status_tsv,
        "rows": len(status_rows),
        "counts": counts,
        "counts_by_bin": bins,
        "method": "thread candidate sequence onto sample1000_refined backbone; keep backbone and fixed motif atoms; PDBFixer rebuilds missing sidechains; optional OpenMM restrained minimization; append reference bn1 ligand",
        "minimize": bool(args.minimize),
    }
    Path(args.summary_json).write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
