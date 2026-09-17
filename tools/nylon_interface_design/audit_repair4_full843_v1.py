"""Replay four unchanged repaired PA66 poses over all 843 sites; preserve conditional denominators."""
import os
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
from pathlib import Path
import json,hashlib,datetime
import numpy as np
from scipy.spatial import cKDTree
R=Path(__file__).resolve().parent;D=R/'route56_N10_20260915_v1';B=D/'low_coverage_repair4_v1'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def angle(a,b,c):
 u=a-b;v=c-b
 return float(np.degrees(np.arccos(np.clip(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v),-1,1))))
def main():
 out=B/'full843_replay_v1.json'
 assert not out.exists()
 jobs=json.loads((D/'low_coverage_rfd12_v2/manifest.json').read_text())['jobs']
 old=json.loads((B/'peptide_donor_audit_v1.json').read_text())['reports']
 replay=json.loads((B/'independent_replay_v2.json').read_text())['reports']
 sitesfile=D/'N10_fixed_minimal_motif_material_v1/sites.json';sites=json.loads(sitesfile.read_text())
 ids=[i for i,s in enumerate(sites) if s['material']=='PA66'];assert len(ids)==843
 envfile=R.parent/'interface_surface_robustness_20260908_v1/inputs_v2/PA66_parent_atoms.npz'
 env=np.load(envfile);box=env['box_A'];px=env['xyz_A'].copy();px[:,:2]%=box[:2];el=env['elements'].astype(str)
 radii={'C':1.7,'N':1.55,'O':1.52,'S':1.8}
 trees={e:cKDTree(px[el==e],boxsize=[*box[:2],0]) for e in np.unique(el[el!='H'])}
 records=[]
 for ref in old:
  name=ref['name'];src=B/name/'after.pdb';assert sha(src)==ref['source_sha256']
  j=next(j for j in jobs if j['name']==name)
  fitfile=Path(j['source_interface'])/'fitted_fields.npz';shapefile=D/'low_coverage_rfd12_v2'/name/'shape_field.npz'
  fit=np.load(fitfile);shape=np.load(shapefile);delta=shape['centering_roundoff_A'];translation=shape['fit_to_shape_translation_A']
  atoms={};heavy=[]
  for l in src.read_text().splitlines():
   if l.startswith('ATOM'):
    key=(int(l[22:26]),l[12:16].strip());v=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]);atoms[key]=v
    if l[76:78].strip() in radii:heavy.append((key,v,radii[l[76:78].strip()]))
  q=np.array([v for _,v,_ in heavy]);rr=np.array([r for _,_,r in heavy]);qfit=q+delta-translation
  rows=[]
  for si in ids:
   st=sites[si];F=np.array(st['frame'])@fit['rot'];w=qfit@F.T+st['center'];w[:,:2]%=box[:2]
   margin=np.full(len(w),np.inf)
   for e,t in trees.items():margin=np.minimum(margin,t.query(w)[0]-rr-radii[e]+.4)
   amide=np.array(st['amide'])@fit['rot']+translation-delta;C,O=amide[:2];ds=[]
   for donor in ref['donors']:
    v=atoms[(donor['residue'],donor['atom'])]
    a=max(angle(v,atoms[(donor['residue'],h)],O) for h in donor['hydrogens'])
    dist=float(np.linalg.norm(v-O));ds.append({'distance_A':dist,'DHO_deg':a,'pass':2.7<=dist<=3.2 and a>=140})
   og=atoms[(1,'OG1')];cb=atoms[(1,'CB')]
   attack=[float(np.linalg.norm(og-C)),angle(O,C,og),angle(cb,og,C)]
   ap=3<=attack[0]<=3.4 and 95<=attack[1]<=115 and 100<=attack[2]<=130
   clear=not bool(np.any(margin< -1e-8))
   rows.append({'site_index':si,'site_id':st['site_id'],'bin':st['bin'],'clashing_atoms':int(np.sum(margin< -1e-8)),'minimum_margin_A':float(margin.min()),'material_clear':clear,'attack':attack,'attack_pass':ap,'donors':ds,'dual_clear':bool(clear and ap and all(d['pass'] for d in ds)),'single_clear':bool(clear and ap and any(d['pass'] for d in ds))})
  byid={v['site_index']:v for v in rows}
  oldrep=next(v for v in replay if v['name']==name)['states'][1]
  for v in oldrep['material_replay']:
   new=byid[v['site_index']];assert new['clashing_atoms']==v['clashing_atoms']
   assert abs(new['minimum_margin_A']-v['minimum_margin_A'])<1e-8
  for v in ref['site_rows']:
   new=byid[v['site_index']]
   assert new['dual_clear']==bool(v['attack_pass'] and v['dual_pass'] and v['material_clear'])
   assert new['single_clear']==bool(v['attack_pass'] and v['single_pass'] and v['material_clear'])
  groups=[]
  for bn in range(4):
   rrw=[v for v in rows if v['bin']==bn]
   groups.append({'bin':bn,'n':len(rrw),'material_clear':sum(v['material_clear'] for v in rrw),'dual_clear':sum(v['dual_clear'] for v in rrw),'single_clear':sum(v['single_clear'] for v in rrw)})
  rec={'name':name,'source':str(src),'sha256':sha(src),'fit_sha256':sha(fitfile),'shape_sha256':sha(shapefile),'all_sites_n':843,'conditional_n':ref['denominator'],'conditional_dual_clear':ref['attack_dual_clear_sites'],'all_dual_clear':sum(v['dual_clear'] for v in rows),'all_single_clear':sum(v['single_clear'] for v in rows),'all_material_clear':sum(v['material_clear'] for v in rows),'bins':groups,'rows':rows,'original_subset_reproduced':True}
  records.append(rec);print(json.dumps({k:v for k,v in rec.items() if k not in ['rows','bins']}),flush=True)
 report={'status':'COMPLETE_FIXED_POSE_FULL843_REPLAY_NOT_DESIGN_PASS','sites_sha256':sha(sitesfile),'parent_atoms_sha256':sha(envfile),'records':records,'method':'Same rigid pose relative to each amide frame; no new pose search. xy minimum-image vdw overlap allowance 0.4 A; unchanged attack and donor gates. Existing conditional subset replay must match exactly.','limitations':['Placeholder sequences; not independently predicted proteins','Saved explicit donor H, no H reoptimization','Does not override internal clash or peptide planarity failures','Counts are not binding affinity or activity']}
 with out.open('x') as h:json.dump(report,h,indent=2)
 with (B/'full843_replay_v1_RUNBOOK.md').open('x') as h:h.write('CPU scipy/numpy environment: ../unbiased_100chain_v1/env/bin/python audit_repair4_full843_v1.py\nInputs: unchanged after.pdb of four low_coverage_repair4_v1 cases; sites.json; PA66_parent_atoms.npz; original fitted_fields and shape transforms. Outputs: full843_replay_v1.json. Existing subset counts and margins checked exactly. Versioned exclusive output; no overwrite. All-site fixed-pose audit, not an orientation search or full design PASS.\n')
 return out
if __name__=='__main__':
 status='FAILED';error=None
 try:out=main();status='COMPLETE'
 except Exception as e:error=repr(e);raise
 finally:
  with (R/'RUN_LOG.jsonl').open('a') as h:h.write(json.dumps({'time':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script':__file__,'script_sha256':sha(Path(__file__)),'status':status,'error':error,'output':str(B/'full843_replay_v1.json'),'CPU_only':True})+'\n')
