#!/usr/bin/env python3
import argparse
import json
import math
from pathlib import Path

BACKBONE = {'N', 'CA', 'C', 'O'}
THRESHOLDS = {'catalytic_backbone_rmsd_max_A': 0.10, 'fixed_backbone_rmsd_max_A': 0.10, 'ligand_heavy_rmsd_max_A': 0.10}
CATALYTIC = {('A', 95), ('A', 126), ('A', 128)}

def parse_residue_tokens(text):
    return {(tok[0], int(tok[1:])) for tok in text.replace(',', ' ').split()}

def parse_pdb(path):
    atoms=[]
    with open(path) as fh:
        for line in fh:
            rec=line[:6].strip()
            if rec not in {'ATOM','HETATM'}:
                continue
            try:
                atom=line[12:16].strip(); resname=line[17:20].strip(); chain=line[21].strip(); resseq=int(line[22:26])
                coord=(float(line[30:38]), float(line[38:46]), float(line[46:54]))
                element=(line[76:78].strip() or atom[0]).upper()
            except Exception:
                continue
            atoms.append({'rec':rec,'atom':atom,'resname':resname,'chain':chain,'resseq':resseq,'coord':coord,'heavy':not(element=='H' or atom.startswith('H'))})
    return atoms

def distance(a,b):
    return math.sqrt(sum((a['coord'][i]-b['coord'][i])**2 for i in range(3)))

def rmsd(ref_atoms, cand_atoms, residues=None, backbone=False, ligand=False):
    ref={}
    for a in ref_atoms:
        if ligand:
            if not (a['rec']=='HETATM' and a['chain']=='X' and a['resseq']==1 and a['resname']=='bu2' and a['heavy']):
                continue
        else:
            if (a['chain'], a['resseq']) not in residues:
                continue
            if backbone and a['atom'] not in BACKBONE:
                continue
        ref[(a['rec'],a['chain'],a['resseq'],a['resname'],a['atom'])]=a
    cand={(a['rec'],a['chain'],a['resseq'],a['resname'],a['atom']):a for a in cand_atoms}
    sq=[]; missing=[]
    for k,a in ref.items():
        b=cand.get(k)
        if b is None:
            missing.append('|'.join(map(str,k)))
        else:
            d=distance(a,b); sq.append(d*d)
    return (math.sqrt(sum(sq)/len(sq)) if sq else None), len(sq), missing

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--reference', required=True)
    p.add_argument('--candidate', required=True)
    p.add_argument('--fixed-residues', required=True)
    p.add_argument('--out', required=True)
    args=p.parse_args()
    fixed=parse_residue_tokens(Path(args.fixed_residues).read_text())
    ref=parse_pdb(args.reference); cand=parse_pdb(args.candidate)
    cat_r, cat_n, cat_m = rmsd(ref, cand, CATALYTIC, backbone=True)
    fix_r, fix_n, fix_m = rmsd(ref, cand, fixed, backbone=True)
    lig_r, lig_n, lig_m = rmsd(ref, cand, ligand=True)
    reasons=[]
    if cat_r is None or cat_r > THRESHOLDS['catalytic_backbone_rmsd_max_A'] or cat_m: reasons.append('catalytic_backbone_rmsd_or_missing')
    if fix_r is None or fix_r > THRESHOLDS['fixed_backbone_rmsd_max_A'] or fix_m: reasons.append('fixed_backbone_rmsd_or_missing')
    if lig_r is None or lig_r > THRESHOLDS['ligand_heavy_rmsd_max_A'] or lig_m: reasons.append('ligand_heavy_rmsd_or_missing')
    result={'status':'pass' if not reasons else 'fail','reference':args.reference,'candidate':args.candidate,'thresholds':THRESHOLDS,'metrics':{'catalytic_backbone_rmsd_A':cat_r,'catalytic_backbone_atom_pairs':cat_n,'fixed_backbone_rmsd_A':fix_r,'fixed_backbone_atom_pairs':fix_n,'ligand_heavy_rmsd_A':lig_r,'ligand_heavy_atom_pairs':lig_n},'missing':{'catalytic_backbone_atoms':cat_m[:20],'fixed_backbone_atoms':fix_m[:20],'ligand_atoms':lig_m[:20]},'fail_reasons':reasons,'note':'LigandMPNN backbone outputs omit protein side chains; full side-chain key-distance filter is deferred to packed/full/PLACER structures.'}
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(result, indent=2, sort_keys=True))
    print(json.dumps(result, indent=2, sort_keys=True))
if __name__=='__main__':
    main()
