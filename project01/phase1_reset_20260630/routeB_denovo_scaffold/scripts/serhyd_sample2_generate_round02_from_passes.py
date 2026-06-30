#!/usr/bin/env python3
from __future__ import annotations
import csv, json, random, argparse, collections
from pathlib import Path

BINS=[90,80,70,60,50]
QUOTAS={90:30,80:160,70:100,60:100,50:120}
GROUPS={
 'A':'AGSTV','G':'AGST','S':'STNQAG','T':'STNQAV','V':'VILMAT','I':'ILVMF','L':'LIVMF','M':'MLIVF',
 'F':'FYWILM','Y':'YFWH','W':'WFY','D':'DEQN','E':'EDQK','N':'NQSTDE','Q':'QNEDK','K':'KRQE','R':'RKHQ','H':'HYRKN','C':'STAG','P':'AGSTV'
}
FALLBACK='ACDEFGHIKLMNQRSTVWY'

def parse_args():
 ap=argparse.ArgumentParser()
 ap.add_argument('--root',default='/data/bht/project01_phase1_reset_gpu/denovo_scaffold_routeB/serine_hydrolase')
 ap.add_argument('--round01-manifest',default='manifests/routeB_serhyd_sample2_identity_bins_round01_manifest.tsv')
 ap.add_argument('--round01-gate',default='manifests/routeB_serhyd_sample2_identity_bins_round01_active_pocket_gate.tsv')
 ap.add_argument('--round-id',default='round02')
 ap.add_argument('--seed',type=int,default=263002)
 return ap.parse_args()

def read_tsv(path):
 with Path(path).open(newline='') as f: return list(csv.DictReader(f,delimiter='\t'))

def fixed_positions(text):
 out=[]
 for tok in text.split():
  if tok.startswith('A') and tok[1:].isdigit(): out.append(int(tok[1:]))
 return sorted(set(out))

def target_mutations(length,bin_id): return int(round(length*(100-bin_id)/100.0))

def conservative(old,rng):
 pool=''.join(a for a in GROUPS.get(old,FALLBACK) if a!=old)
 if pool and rng.random()<0.85: return rng.choice(pool)
 return rng.choice(''.join(a for a in FALLBACK if a!=old))

def weighted_sample(items,weights,k,rng):
 pool=list(items); w=list(weights); out=[]
 for _ in range(k):
  total=sum(w)
  x=rng.random()*total; acc=0.0; idx=len(pool)-1
  for i,wi in enumerate(w):
   acc+=wi
   if acc>=x:
    idx=i; break
  out.append(pool.pop(idx)); w.pop(idx)
 return out

def main():
 args=parse_args(); root=Path(args.root); rng=random.Random(args.seed)
 base=read_tsv(root/args.round01_manifest); gate=read_tsv(root/args.round01_gate)
 ref=base[0]['sequence']; length=len(ref); fixed=set(fixed_positions(base[0]['fixed_positions']))
 mutable=[i for i in range(1,length+1) if i not in fixed]
 existing={r['sequence'] for r in base}
 pass_freq=collections.Counter(); fail_freq=collections.Counter(); pass_alts=collections.defaultdict(collections.Counter)
 for r in gate:
  status=r['routeB_active_pocket_only_gate']
  for m in r['mutations'].split(';'):
   if not m: continue
   pos_s,sub=m.split(':',1); pos=int(pos_s[1:]); new=sub.split('>')[1]
   if status=='PASS':
    pass_freq[pos]+=1; pass_alts[pos][new]+=1
   else:
    fail_freq[pos]+=1
 scores={}
 for pos in mutable:
  scores[pos]=(pass_freq[pos]+0.5)/(fail_freq[pos]+1.0)
 rows=[]; fields=['sample_id','bin','identity','mutation_count','sequence_length','sequence','postseq_predicted_pdb','reference_sample','fixed_positions','pocket_ca8_positions','active_site_positions','mutations','round02_position_policy']
 for bin_id in BINS:
  nmut=target_mutations(length,bin_id); made=0; attempts=0
  sorted_pos=sorted(mutable,key=lambda p:(scores[p],pass_freq[p]),reverse=True)
  pool_size=min(len(sorted_pos),max(nmut+8,int(nmut*1.45)))
  pool=sorted_pos[:pool_size]
  weights=[max(scores[p],0.01) for p in pool]
  while made<QUOTAS[bin_id] and attempts<QUOTAS[bin_id]*500:
   attempts+=1
   selected=sorted(weighted_sample(pool,weights,nmut,rng))
   seq=list(ref); muts=[]
   for pos in selected:
    old=seq[pos-1]
    if pass_alts[pos] and rng.random()<0.75:
     alts=[]; ww=[]
     for a,c in pass_alts[pos].items():
      if a!=old: alts.append(a); ww.append(c)
     new=weighted_sample(alts,ww,1,rng)[0] if alts else conservative(old,rng)
    else:
     new=conservative(old,rng)
    seq[pos-1]=new; muts.append(f'A{pos}:{old}>{new}')
   seq=''.join(seq)
   if seq in existing: continue
   existing.add(seq); made+=1
   ident=round((length-nmut)/length,6); sample=f'serhyd_sample2_id{bin_id}_{args.round_id}_{made:03d}'
   first=base[0]
   rows.append({'sample_id':sample,'bin':str(bin_id),'identity':f'{ident:.6f}','mutation_count':str(nmut),'sequence_length':str(length),'sequence':seq,'postseq_predicted_pdb':str(root/'outputs'/f'esmfold_identity_bins_sample2_{args.round_id}'/f'{sample}.pdb'),'reference_sample':first['reference_sample'],'fixed_positions':first['fixed_positions'],'pocket_ca8_positions':first.get('pocket_ca8_positions',''),'active_site_positions':first['active_site_positions'],'mutations':';'.join(muts),'round02_position_policy':f'weighted_recombine_round01_pass_fail_pool{pool_size}'})
 manifest=root/'manifests'/f'routeB_serhyd_sample2_identity_bins_{args.round_id}_manifest.tsv'
 manifest.parent.mkdir(parents=True,exist_ok=True)
 with manifest.open('w',newline='') as f:
  w=csv.DictWriter(f,fields,delimiter='\t',lineterminator='\n'); w.writeheader(); w.writerows(rows)
 meta={'round_id':args.round_id,'rows_generated':len(rows),'quotas':QUOTAS,'reference_length':length,'fixed_count':len(fixed),'mutable_count':len(mutable),'round01_pass_rows':sum(1 for r in gate if r['routeB_active_pocket_only_gate']=='PASS'),'manifest':str(manifest),'position_score_top20':[(p,round(scores[p],4),pass_freq[p],fail_freq[p]) for p in sorted(mutable,key=lambda p:scores[p],reverse=True)[:20]]}
 meta_path=root/'manifests'/f'routeB_serhyd_sample2_identity_bins_{args.round_id}_generation_meta.json'
 meta_path.write_text(json.dumps(meta,indent=2,sort_keys=True)+'\n')
 print(json.dumps(meta,indent=2,sort_keys=True))
if __name__=='__main__': main()
