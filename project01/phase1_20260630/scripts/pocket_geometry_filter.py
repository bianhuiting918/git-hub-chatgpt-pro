#!/usr/bin/env python3
import argparse
import json
import math
from pathlib import Path

CATALYTIC = [('A', 95), ('A', 126), ('A', 128)]
KEY_DISTANCES = {
    'Ser128_OG_to_ligand_C1': (('ATOM', 'A', 128, 'OG'), ('HETATM', 'X', 1, 'C1')),
    'Ser126_OG_to_ligand_C1': (('ATOM', 'A', 126, 'OG'), ('HETATM', 'X', 1, 'C1')),
    'His95_NE2_to_Ser128_OG': (('ATOM', 'A', 95, 'NE2'), ('ATOM', 'A', 128, 'OG')),
    'His95_NE2_to_Ser126_OG': (('ATOM', 'A', 95, 'NE2'), ('ATOM', 'A', 126, 'OG')),
    'His95_ND1_to_Ser128_OG': (('ATOM', 'A', 95, 'ND1'), ('ATOM', 'A', 128, 'OG')),
}
THRESHOLDS = {
    'catalytic_heavy_rmsd_max_A': 0.50,
    'fixed_heavy_rmsd_max_A': 0.75,
    'ligand_heavy_rmsd_max_A': 0.50,
    'key_distance_abs_delta_max_A': 0.35,
}

def parse_residue_tokens(text):
    residues = []
    for tok in text.replace(',', ' ').split():
        chain = tok[0]
        resseq = int(tok[1:])
        residues.append((chain, resseq))
    return residues

def parse_pdb(path):
    atoms = []
    with open(path) as fh:
        for line in fh:
            rec = line[:6].strip()
            if rec not in {'ATOM', 'HETATM'}:
                continue
            try:
                atom = line[12:16].strip()
                resname = line[17:20].strip()
                chain = line[21].strip()
                resseq = int(line[22:26])
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                element = line[76:78].strip() or atom[0]
            except Exception:
                continue
            atoms.append({
                'rec': rec, 'atom': atom, 'resname': resname, 'chain': chain,
                'resseq': resseq, 'coord': (x, y, z), 'element': element.upper(),
                'heavy': not (element.upper() == 'H' or atom.startswith('H')),
            })
    return atoms

def atom_index(atoms):
    idx = {}
    for a in atoms:
        idx[(a['rec'], a['chain'], a['resseq'], a['atom'])] = a
    return idx

def dist(a, b):
    return math.sqrt(sum((a['coord'][i] - b['coord'][i]) ** 2 for i in range(3)))

def rmsd_for(ref_atoms, cand_atoms, residues=None, ligand=False):
    ref_map = {}
    for a in ref_atoms:
        if not a['heavy']:
            continue
        if ligand:
            if not (a['rec'] == 'HETATM' and a['chain'] == 'X' and a['resseq'] == 1 and a['resname'] == 'bu2'):
                continue
        elif residues is not None:
            if (a['chain'], a['resseq']) not in residues:
                continue
        else:
            continue
        key = (a['rec'], a['chain'], a['resseq'], a['resname'], a['atom'])
        ref_map[key] = a
    cand_map = {}
    for a in cand_atoms:
        if not a['heavy']:
            continue
        key = (a['rec'], a['chain'], a['resseq'], a['resname'], a['atom'])
        cand_map[key] = a
    sq = []
    missing = []
    for key, ref_atom in ref_map.items():
        cand_atom = cand_map.get(key)
        if cand_atom is None:
            missing.append('|'.join(map(str, key)))
        else:
            d = dist(ref_atom, cand_atom)
            sq.append(d * d)
    if not sq:
        return None, 0, missing
    return math.sqrt(sum(sq) / len(sq)), len(sq), missing

def key_distances(atoms):
    idx = atom_index(atoms)
    out = {}
    missing = []
    for name, (a_key, b_key) in KEY_DISTANCES.items():
        a = idx.get(a_key)
        b = idx.get(b_key)
        if a is None or b is None:
            missing.append(name)
        else:
            out[name] = dist(a, b)
    return out, missing

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--reference', required=True)
    p.add_argument('--candidate', required=True)
    p.add_argument('--fixed-residues', required=True, help='Text like: A57 A61 ...')
    p.add_argument('--out', required=True)
    args = p.parse_args()

    fixed_residues = parse_residue_tokens(Path(args.fixed_residues).read_text())
    catalytic = set(CATALYTIC)
    fixed = set(fixed_residues)
    ref_atoms = parse_pdb(args.reference)
    cand_atoms = parse_pdb(args.candidate)

    cat_rmsd, cat_n, cat_missing = rmsd_for(ref_atoms, cand_atoms, residues=catalytic)
    fixed_rmsd, fixed_n, fixed_missing = rmsd_for(ref_atoms, cand_atoms, residues=fixed)
    lig_rmsd, lig_n, lig_missing = rmsd_for(ref_atoms, cand_atoms, ligand=True)
    ref_dist, ref_missing = key_distances(ref_atoms)
    cand_dist, cand_missing = key_distances(cand_atoms)
    deltas = {k: cand_dist[k] - ref_dist[k] for k in sorted(ref_dist) if k in cand_dist}

    reasons = []
    if cat_rmsd is None or cat_rmsd > THRESHOLDS['catalytic_heavy_rmsd_max_A']:
        reasons.append('catalytic_heavy_rmsd')
    if fixed_rmsd is None or fixed_rmsd > THRESHOLDS['fixed_heavy_rmsd_max_A']:
        reasons.append('fixed_heavy_rmsd')
    if lig_rmsd is None or lig_rmsd > THRESHOLDS['ligand_heavy_rmsd_max_A']:
        reasons.append('ligand_heavy_rmsd')
    if ref_missing or cand_missing:
        reasons.append('missing_key_distance_atoms')
    max_abs_delta = max([abs(v) for v in deltas.values()] or [999.0])
    if max_abs_delta > THRESHOLDS['key_distance_abs_delta_max_A']:
        reasons.append('key_distance_delta')
    if cat_missing or fixed_missing or lig_missing:
        reasons.append('missing_matched_atoms')

    result = {
        'status': 'pass' if not reasons else 'fail',
        'reference': args.reference,
        'candidate': args.candidate,
        'thresholds': THRESHOLDS,
        'metrics': {
            'catalytic_heavy_rmsd_A': cat_rmsd,
            'catalytic_heavy_atom_pairs': cat_n,
            'fixed_heavy_rmsd_A': fixed_rmsd,
            'fixed_heavy_atom_pairs': fixed_n,
            'ligand_heavy_rmsd_A': lig_rmsd,
            'ligand_heavy_atom_pairs': lig_n,
            'reference_key_distances_A': ref_dist,
            'candidate_key_distances_A': cand_dist,
            'key_distance_deltas_A': deltas,
            'max_abs_key_distance_delta_A': max_abs_delta,
        },
        'missing': {
            'reference_key_distances': ref_missing,
            'candidate_key_distances': cand_missing,
            'catalytic_atoms': cat_missing[:20],
            'fixed_atoms': fixed_missing[:20],
            'ligand_atoms': lig_missing[:20],
        },
        'fail_reasons': reasons,
    }
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(result, indent=2, sort_keys=True))
    print(json.dumps(result, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
