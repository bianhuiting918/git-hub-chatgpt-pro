#!/usr/bin/env python3
from __future__ import annotations
import argparse, csv, json, random
from pathlib import Path
from math import sqrt

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
ACTIVE = {40,49,56}
BINS = [90,80,70,60,50]
GROUPS = {
    'A':'AGSTV', 'G':'AGST', 'S':'STNQAG', 'T':'STNQAV',
    'V':'VILMAT', 'I':'ILVMF', 'L':'LIVMF', 'M':'MLIVF',
    'F':'FYWILM', 'Y':'YFWH', 'W':'WFY',
    'D':'DEQN', 'E':'EDQK', 'N':'NQSTDE', 'Q':'QNEDK',
    'K':'KRQE', 'R':'RKHQ', 'H':'HYRKN',
    'C':'STAG', 'P':'AGSTV'
}
FALLBACK = 'ACDEFGHIKLMNQRSTVWY'

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument('--root', default='/data/bht/project01_phase1_reset_gpu/denovo_scaffold_routeB/serine_hydrolase')
    ap.add_argument('--rfaa-pdb', default='outputs/rfaa_serhyd_motif_H95_E102_S128_bu2_batch5/sample_2.pdb')
    ap.add_argument('--manifest-in', default='manifests/routeB_serhyd_rfaabatch5_seq5_esmfold_manifest.tsv')
    ap.add_argument('--reference-sample', default='routeB_serhyd_sample_2_lmpnn_03')
    ap.add_argument('--round-id', default='round01')
    ap.add_argument('--per-bin', type=int, default=40)
    ap.add_argument('--seed', type=int, default=263001)
    ap.add_argument('--pocket-ca-cutoff', type=float, default=8.0)
    return ap.parse_args()

def read_reference_sequence(path: Path, sample_id: str) -> str:
    for row in csv.DictReader(path.open(), delimiter='\t'):
        if row['sample_id'] == sample_id:
            return row['sequence'].strip()
    raise SystemExit(f'reference sample not found: {sample_id}')

def xyz(line: str):
    try:
        return float(line[30:38]), float(line[38:46]), float(line[46:54])
    except Exception:
        parts = line.split()
        return float(parts[6]), float(parts[7]), float(parts[8])

def parse_residues_ligand_and_ca_pocket(pdb: Path, cutoff: float):
    residues = []
    ligand_atoms = []
    ca_atoms = {}
    seen = set()
    for line in pdb.open():
        rec = line[:6].strip()
        if rec not in {'ATOM','HETATM'}:
            continue
        name = line[12:16].strip()
        elem = (line[76:78].strip() or name[0]).upper()
        if elem == 'H':
            continue
        x, y, z = xyz(line)
        if rec == 'HETATM' and line[17:20].strip().lower() == 'bu2':
            ligand_atoms.append((x,y,z,name))
        elif rec == 'ATOM' and line[21].strip() == 'A':
            resid = int(line[22:26])
            key = ('A', resid, line[26].strip())
            if key not in seen:
                seen.add(key)
                residues.append((resid, AA3.get(line[17:20].strip(), 'X')))
            if name == 'CA':
                ca_atoms[resid] = (x,y,z)
    cutoff2 = cutoff * cutoff
    pocket = set()
    ca_distances = {}
    for resid, coord in ca_atoms.items():
        d2 = min((coord[0]-lx)**2 + (coord[1]-ly)**2 + (coord[2]-lz)**2 for lx,ly,lz,_ in ligand_atoms)
        ca_distances[resid] = sqrt(d2)
        if d2 <= cutoff2:
            pocket.add(resid)
    return residues, ligand_atoms, pocket, ca_distances

def choose_substitution(ref_aa: str, rng: random.Random) -> str:
    group = ''.join(a for a in GROUPS.get(ref_aa, FALLBACK) if a != ref_aa)
    if group and rng.random() < 0.85:
        return rng.choice(group)
    pool = ''.join(a for a in FALLBACK if a != ref_aa)
    return rng.choice(pool)

def target_mutations(length: int, bin_id: int) -> int:
    return int(round(length * (100 - bin_id) / 100.0))

def make_variant(ref: str, mutable_positions: list[int], nmut: int, rng: random.Random):
    selected = sorted(rng.sample(mutable_positions, nmut))
    seq = list(ref)
    muts = []
    for pos in selected:
        old = seq[pos-1]
        new = choose_substitution(old, rng)
        seq[pos-1] = new
        muts.append(f'A{pos}:{old}>{new}')
    return ''.join(seq), muts

def main() -> int:
    args = parse_args()
    root = Path(args.root)
    ref = read_reference_sequence(root / args.manifest_in, args.reference_sample)
    rfaa_pdb = root / args.rfaa_pdb
    residues, ligand_atoms, pocket, ca_distances = parse_residues_ligand_and_ca_pocket(rfaa_pdb, args.pocket_ca_cutoff)
    if len(ref) != len(residues):
        raise SystemExit(f'length mismatch ref={len(ref)} pdb_residues={len(residues)}')
    fixed = set(pocket) | set(ACTIVE)
    mutable = [i for i in range(1, len(ref)+1) if i not in fixed]
    out_dir = root / 'outputs' / f'esmfold_identity_bins_sample2_{args.round_id}'
    manifest = root / 'manifests' / f'routeB_serhyd_sample2_identity_bins_{args.round_id}_manifest.tsv'
    meta_json = root / 'manifests' / f'routeB_serhyd_sample2_identity_bins_{args.round_id}_generation_meta.json'
    out_dir.mkdir(parents=True, exist_ok=True)
    rng = random.Random(args.seed)
    fields = ['sample_id','bin','identity','mutation_count','sequence_length','sequence','postseq_predicted_pdb','reference_sample','fixed_positions','pocket_ca8_positions','active_site_positions','mutations']
    rows = []
    seen = {ref}
    skipped_bins = {}
    for bin_id in BINS:
        nmut = target_mutations(len(ref), bin_id)
        if nmut > len(mutable):
            skipped_bins[str(bin_id)] = f'target_mutations_{nmut}_gt_mutable_{len(mutable)}'
            continue
        made = 0
        attempts = 0
        while made < args.per_bin and attempts < args.per_bin * 200:
            attempts += 1
            seq, muts = make_variant(ref, mutable, nmut, rng)
            if seq in seen:
                continue
            seen.add(seq)
            made += 1
            ident = round((len(ref) - nmut) / len(ref), 6)
            sample = f'serhyd_sample2_id{bin_id}_{args.round_id}_{made:03d}'
            rows.append({
                'sample_id': sample,
                'bin': str(bin_id),
                'identity': f'{ident:.6f}',
                'mutation_count': str(nmut),
                'sequence_length': str(len(seq)),
                'sequence': seq,
                'postseq_predicted_pdb': str(out_dir / f'{sample}.pdb'),
                'reference_sample': args.reference_sample,
                'fixed_positions': ' '.join(f'A{i}' for i in sorted(fixed)),
                'pocket_ca8_positions': ' '.join(f'A{i}' for i in sorted(pocket)),
                'active_site_positions': ' '.join(f'A{i}' for i in sorted(ACTIVE)),
                'mutations': ';'.join(muts),
            })
    with manifest.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    meta = {
        'round_id': args.round_id,
        'reference_sample': args.reference_sample,
        'reference_length': len(ref),
        'reference_sequence': ref,
        'rfaa_pdb': str(rfaa_pdb),
        'ligand': 'bu2',
        'ligand_heavy_atoms': len(ligand_atoms),
        'pocket_definition_for_sequence_stage': f'active_site plus residues with CA-to-bu2 <= {args.pocket_ca_cutoff} A; RFAA sidechain/heavy atoms are diagnostic only',
        'active_site_positions': sorted(ACTIVE),
        'pocket_ca8_positions': sorted(pocket),
        'fixed_positions': sorted(fixed),
        'fixed_count': len(fixed),
        'mutable_count': len(mutable),
        'per_bin_requested': args.per_bin,
        'rows_generated': len(rows),
        'skipped_bins': skipped_bins,
        'target_mutations_by_bin': {str(b): target_mutations(len(ref), b) for b in BINS},
        'ca_distance_to_ligand_by_residue': {str(k): round(v, 4) for k,v in sorted(ca_distances.items())},
        'manifest': str(manifest),
        'output_dir': str(out_dir),
        'identity_is_relative_to': args.reference_sample,
        'hard_gate_after_prediction': 'fixed_positions_unchanged + motif_ca_rmsd_A<=1.0',
    }
    meta_json.write_text(json.dumps(meta, indent=2, sort_keys=True) + '\n')
    print(json.dumps(meta, indent=2, sort_keys=True))
    return 0

if __name__ == '__main__':
    raise SystemExit(main())
