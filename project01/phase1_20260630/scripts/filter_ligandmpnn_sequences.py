#!/usr/bin/env python3
import argparse
import json
from pathlib import Path

AA3 = {
    'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I',
    'LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'
}

def parse_residue_tokens(text):
    out = set()
    for tok in text.replace(',', ' ').split():
        out.add((tok[0], int(tok[1:])))
    return out

def pdb_sequence(path, chain='A'):
    residues = []
    seen = set()
    with open(path) as fh:
        for line in fh:
            if not line.startswith('ATOM'):
                continue
            ch = line[21].strip()
            if ch != chain:
                continue
            resseq = int(line[22:26])
            icode = line[26].strip()
            key = (ch, resseq, icode)
            if key in seen:
                continue
            seen.add(key)
            resname = line[17:20].strip()
            residues.append((ch, resseq, AA3.get(resname, 'X')))
    return residues, ''.join(a for _, _, a in residues)

def parse_fasta(path):
    records = []
    header = None
    seq = []
    with open(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header is not None:
                    records.append((header, ''.join(seq)))
                header = line[1:]
                seq = []
            else:
                seq.append(line)
    if header is not None:
        records.append((header, ''.join(seq)))
    return records

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--fasta', required=True)
    p.add_argument('--reference-pdb', required=True)
    p.add_argument('--fixed-residues', required=True)
    p.add_argument('--redesignable-residues', required=True)
    p.add_argument('--bin-label', required=True)
    p.add_argument('--mutation-min', type=int, required=True)
    p.add_argument('--mutation-max', type=int, required=True)
    p.add_argument('--out-tsv', required=True)
    p.add_argument('--out-json', required=True)
    args = p.parse_args()

    residues, ref = pdb_sequence(args.reference_pdb, 'A')
    fixed = parse_residue_tokens(Path(args.fixed_residues).read_text())
    redesignable = parse_residue_tokens(Path(args.redesignable_residues).read_text())
    records = parse_fasta(args.fasta)
    rows = []
    selected = []
    seen_seq = set()

    for idx, (header, seq) in enumerate(records):
        if len(seq) != len(ref):
            rows.append({'header': header, 'status': 'fail', 'reason': 'length_mismatch', 'length': len(seq)})
            continue
        mutations = []
        fixed_mut = []
        outside_mut = []
        for pos, ((chain, resseq, ref_aa), aa) in enumerate(zip(residues, seq), start=1):
            if aa != ref_aa:
                mutations.append(f'{chain}{resseq}:{ref_aa}>{aa}')
                if (chain, resseq) in fixed:
                    fixed_mut.append(f'{chain}{resseq}:{ref_aa}>{aa}')
                if (chain, resseq) not in redesignable:
                    outside_mut.append(f'{chain}{resseq}:{ref_aa}>{aa}')
        mut_count = len(mutations)
        identity = (len(ref) - mut_count) / len(ref)
        reasons = []
        if not (args.mutation_min <= mut_count <= args.mutation_max):
            reasons.append('mutation_count_out_of_bin')
        if fixed_mut:
            reasons.append('fixed_residue_mutation')
        if outside_mut:
            reasons.append('outside_redesignable_mutation')
        if seq in seen_seq:
            reasons.append('duplicate_sequence')
        status = 'pass' if not reasons else 'fail'
        row = {
            'bin': args.bin_label,
            'record_index': idx,
            'header': header,
            'status': status,
            'identity': round(identity, 6),
            'mutation_count': mut_count,
            'fixed_mutation_count': len(fixed_mut),
            'outside_redesignable_mutation_count': len(outside_mut),
            'sequence': seq,
            'mutations': ';'.join(mutations),
            'fail_reasons': ';'.join(reasons),
        }
        rows.append(row)
        seen_seq.add(seq)
        if status == 'pass':
            selected.append(row)

    fields = ['bin','record_index','status','identity','mutation_count','fixed_mutation_count','outside_redesignable_mutation_count','header','fail_reasons','sequence','mutations']
    Path(args.out_tsv).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out_tsv, 'w') as out:
        out.write('\t'.join(fields) + '\n')
        for r in rows:
            out.write('\t'.join(str(r.get(f, '')) for f in fields) + '\n')
    summary = {
        'bin': args.bin_label,
        'fasta': args.fasta,
        'reference_pdb': args.reference_pdb,
        'sequence_length': len(ref),
        'records_total': len(records),
        'pass_count': len(selected),
        'mutation_min': args.mutation_min,
        'mutation_max': args.mutation_max,
        'selected_headers': [r['header'] for r in selected[:10]],
        'selected_tsv': args.out_tsv,
    }
    Path(args.out_json).write_text(json.dumps(summary, indent=2, sort_keys=True))
    print(json.dumps(summary, indent=2, sort_keys=True))

if __name__ == '__main__':
    main()
