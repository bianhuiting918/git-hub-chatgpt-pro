#!/usr/bin/env python3
"""Export and audit an isolated active beta-chain ASH306 pdb2gmx probe."""
import argparse
import hashlib
import json
import pathlib
import re

CHAIN_H_FIRST = 8949
CHAIN_H_LAST = 10272
EXPECTED_CHAIN_ATOMS = 1325
EXPECTED_SOURCE_SHA256 = "182f68a41c72c490e0046323c9b8a0f8514ccf7af7dc4d578dcb87342c1c8e2d"


def sha256(path):
    h = hashlib.sha256()
    with pathlib.Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_gro(path):
    lines = pathlib.Path(path).read_text().splitlines()
    count = int(lines[1])
    atoms = []
    for index, line in enumerate(lines[2:2 + count], 1):
        atoms.append({
            "index": index,
            "resid": int(line[:5]),
            "resname": line[5:10].strip(),
            "atomname": line[10:15].strip(),
            "xyz_nm": tuple(float(line[a:b]) for a, b in ((20, 28), (28, 36), (36, 44))),
        })
    return atoms


def element(atomname):
    letters = "".join(character for character in atomname if character.isalpha()).upper()
    if not letters or letters[0] not in "HCNOS":
        raise ValueError(f"unsupported element for {atomname}")
    return letters[0]


def export_pdb(gro, pdb):
    if sha256(gro) != EXPECTED_SOURCE_SHA256:
        raise ValueError("immutable 434 ps source SHA mismatch")
    atoms = read_gro(gro)
    selected = atoms[CHAIN_H_FIRST - 1:CHAIN_H_LAST]
    if len(selected) != CHAIN_H_LAST - CHAIN_H_FIRST + 1:
        raise ValueError("active chain range mismatch")
    lines = []
    renamed = 0
    chain_id = "H"
    for serial, atom in enumerate(selected, 1):
        resid, resname = atom["resid"], atom["resname"]
        if resid == 306 and resname == "ASP":
            resname = "ASH"
            renamed += 1
        x, y, z = (value * 10.0 for value in atom["xyz_nm"])
        lines.append(
            f"ATOM  {serial:5d} {atom['atomname']:>4s} {resname:>3s} {chain_id}{resid:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {element(atom['atomname']):>2s}"
        )
    if renamed != 12:
        raise ValueError(f"expected 12 existing ASP306 atoms, renamed {renamed}")
    pathlib.Path(pdb).write_text("\n".join(lines) + "\nTER\nEND\n")


def group(atoms, resid, resname):
    return [atom for atom in atoms if atom["resid"] == resid and atom["resname"] == resname]


def audit_probe(source_gro, output_gro, itp, output):
    source = read_gro(source_gro)
    generated = read_gro(output_gro)
    ash = group(generated, 306, "ASH")
    thr = group(generated, 267, "THR")
    ash_names = {atom["atomname"] for atom in ash}
    thr_names = {atom["atomname"] for atom in thr}
    expected_thr = ("THR", {"H1", "H2", "H3", "OG1", "HG1"})
    topology_text = pathlib.Path(itp).read_text(errors="replace")
    topology_ash = bool(re.search(r";\s*residue\s+306\s+ASH\s+rtp\s+ASH\s+q\s+0\.0", topology_text))
    source_key = {(atom["resid"], atom["atomname"]): atom["xyz_nm"] for atom in source[CHAIN_H_FIRST - 1:CHAIN_H_LAST] if not atom["atomname"].startswith("H")}
    generated_key = {(atom["resid"], atom["atomname"]): atom["xyz_nm"] for atom in generated if not atom["atomname"].startswith("H")}
    max_heavy_shift = max(
        max(abs(a - b) for a, b in zip(source_key[key], generated_key[key]))
        for key in source_key
    )
    passed = (
        len(generated) == EXPECTED_CHAIN_ATOMS
        and len(ash) == 13
        and ("ASH", "HD2") == (ash[0]["resname"], "HD2") if "HD2" in ash_names else False
    )
    passed = passed and expected_thr[1].issubset(thr_names) and topology_ash and max_heavy_shift <= 0.0011
    result = {
        "schema_version": 1,
        "candidate_id": "nylc_C18_trueT267_freeGS",
        "status": "PROBE_PASS_ASH306_CHAIN" if passed else "PROBE_FAIL_ASH306_CHAIN",
        "source_sha256": sha256(source_gro),
        "generated_chain_atom_count": len(generated),
        "ash306_atom_names": sorted(ash_names),
        "thr267_atom_names": sorted(thr_names),
        "topology_has_neutral_ASH306": topology_ash,
        "max_heavy_coordinate_shift_nm": max_heavy_shift,
        "generated_gro_sha256": sha256(output_gro),
        "generated_itp_sha256": sha256(itp),
        "interpretation": "Isolated beta-chain protonation probe only; not yet a full-system topology, QM/MM microstate, TS, or PMF input.",
    }
    pathlib.Path(output).write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, sort_keys=True))
    return 0 if passed else 1


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-gro", type=pathlib.Path, required=True)
    parser.add_argument("--export-pdb", type=pathlib.Path)
    parser.add_argument("--generated-gro", type=pathlib.Path)
    parser.add_argument("--generated-itp", type=pathlib.Path)
    parser.add_argument("--audit-output", type=pathlib.Path)
    args = parser.parse_args()
    if args.export_pdb:
        export_pdb(args.source_gro, args.export_pdb)
        return 0
    if not args.generated_gro or not args.generated_itp or not args.audit_output:
        raise ValueError("audit mode requires generated GRO, ITP and output")
    return audit_probe(args.source_gro, args.generated_gro, args.generated_itp, args.audit_output)


if __name__ == "__main__":
    raise SystemExit(main())
