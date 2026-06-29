#!/usr/bin/env python3
import argparse
import csv
import json
import os
import sys
from pathlib import Path

import torch
from torch.nn.functional import cosine_similarity

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}

def parse_residue_tokens(path):
    return [int(tok[1:]) for tok in Path(path).read_text().replace(',', ' ').split()]

def pdb_sequence(path, chain='A'):
    residues=[]; seen=set()
    with open(path) as fh:
        for line in fh:
            if not line.startswith('ATOM') or line[21].strip()!=chain:
                continue
            resseq=int(line[22:26]); icode=line[26].strip(); key=(chain,resseq,icode)
            if key in seen:
                continue
            seen.add(key); residues.append((resseq, AA3.get(line[17:20].strip(), 'X')))
    return residues, ''.join(a for _, a in residues)

def select_positions(seq, residues, wanted):
    pos_to_i={resseq:i for i,(resseq,_) in enumerate(residues)}
    return ''.join(seq[pos_to_i[r]] for r in wanted if r in pos_to_i)

def select_structure_tokens(tokens, residues, wanted):
    pos_to_i={resseq:i for i,(resseq,_) in enumerate(residues)}
    return ''.join(tokens[pos_to_i[r]] for r in wanted if r in pos_to_i)

def to_float(x):
    return float(x.detach().cpu().item())

def main():
    p=argparse.ArgumentParser()
    p.add_argument('--selected-tsv', required=True)
    p.add_argument('--reference-pdb', required=True)
    p.add_argument('--candidate-backbone-dir', required=True)
    p.add_argument('--fixed-residues', required=True)
    p.add_argument('--shell-residues', required=True)
    p.add_argument('--out-json', required=True)
    p.add_argument('--out-tsv', required=True)
    args=p.parse_args()

    protrek_src='/Dell/Dell14/bianht/project01_tools/src/ProTrek'
    model_dir='/Dell/Dell14/bianht/project01_tools/models/ProTrek_650M'
    foldseek='/Dell/Dell14/bianht/project01_tools/foldseek/bin/foldseek'
    sys.path.insert(0, protrek_src)
    os.chdir(protrek_src)
    from utils.foldseek_util import get_struc_seq
    from model.abstract_model import AbstractModel
    AbstractModel.init_optimizers = lambda self: None
    from model.ProTrek.protrek_trimodal_model import ProTrekTrimodalModel

    residues, ref_seq = pdb_sequence(args.reference_pdb, 'A')
    pocket=parse_residue_tokens(args.fixed_residues)
    shell=parse_residue_tokens(args.shell_residues)
    layers={'global': None, 'pocket_5A': pocket, 'shell_5_12A': shell}

    ref_foldseek = get_struc_seq(foldseek, args.reference_pdb, ['A'])['A'][1].lower()
    selected=list(csv.DictReader(open(args.selected_tsv), delimiter='\t'))

    config={
        'protein_config': f'{model_dir}/esm2_t33_650M_UR50D',
        'text_config': f'{model_dir}/BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext',
        'structure_config': f'{model_dir}/foldseek_t30_150M',
        'load_protein_pretrained': False,
        'load_text_pretrained': False,
        'from_checkpoint': f'{model_dir}/ProTrek_650M.pt',
    }
    model=ProTrekTrimodalModel(**config)
    model.eval()
    model.protein_encoder.tokenizer.batch_encode_plus = model.protein_encoder.tokenizer.__call__
    model.structure_encoder.tokenizer.batch_encode_plus = model.structure_encoder.tokenizer.__call__

    ref_seq_emb={}; ref_str_emb={}
    rows=[]
    with torch.no_grad():
        for lname, wanted in layers.items():
            rseq = ref_seq if wanted is None else select_positions(ref_seq, residues, wanted)
            rstr = ref_foldseek if wanted is None else select_structure_tokens(ref_foldseek, residues, wanted)
            ref_seq_emb[lname]=model.get_protein_repr([rseq], batch_size=1)
            ref_str_emb[lname]=model.get_structure_repr([rstr], batch_size=1)
        for r in selected:
            idx=r['record_index']
            seq=r['sequence']
            pdb=str(Path(args.candidate_backbone_dir)/f'denovo_SER_hydrolase_full_input_{idx}.pdb')
            cand_foldseek=get_struc_seq(foldseek, pdb, ['A'])['A'][1].lower()
            item={'record_index': idx, 'identity': float(r['identity']), 'mutation_count': int(r['mutation_count']), 'pdb': pdb}
            for lname, wanted in layers.items():
                cseq = seq if wanted is None else select_positions(seq, residues, wanted)
                cstr = cand_foldseek if wanted is None else select_structure_tokens(cand_foldseek, residues, wanted)
                cseq_emb=model.get_protein_repr([cseq], batch_size=1)
                cstr_emb=model.get_structure_repr([cstr], batch_size=1)
                item[f'{lname}_sequence_cosine_to_ref']=to_float(cosine_similarity(cseq_emb, ref_seq_emb[lname]))
                item[f'{lname}_structure_cosine_to_ref']=to_float(cosine_similarity(cstr_emb, ref_str_emb[lname]))
                item[f'{lname}_candidate_seq_structure_cosine']=to_float(cosine_similarity(cseq_emb, cstr_emb))
                item[f'{lname}_sequence_length']=len(cseq)
                item[f'{lname}_structure_length']=len(cstr)
            rows.append(item)

    summary={'status':'pass','selected_tsv':args.selected_tsv,'reference_pdb':args.reference_pdb,'count':len(rows),'items':rows}
    Path(args.out_json).write_text(json.dumps(summary, indent=2, sort_keys=True))
    fields=list(rows[0].keys()) if rows else []
    with open(args.out_tsv,'w') as fh:
        fh.write('\t'.join(fields)+'\n')
        for row in rows:
            fh.write('\t'.join(str(row[f]) for f in fields)+'\n')
    print(json.dumps({'status':'pass','count':len(rows),'out_json':args.out_json,'out_tsv':args.out_tsv}, indent=2))

if __name__=='__main__':
    main()
