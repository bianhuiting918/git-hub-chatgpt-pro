#!/usr/bin/env python3
"""Eight-site NylC rigid-TDD / independent peptide NH pilot; geometry proxy only."""
from pathlib import Path
import ast,json,hashlib,datetime,time
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
from scipy.stats import qmc
from scipy.spatial.distance import cdist
ROOT=Path(__file__).resolve().parent
PET=ROOT.parent.parent/'pet_dp10_400chain_direct_v1'
SOURCE=PET/'pet_donor_geometry_batch_singleNH_v1.py'
RAD={'C':1.7,'N':1.55,'O':1.52,'S':1.8}
CFG={'n_triad_trials':12000,'vdw_overlap_allowance_A':.4}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
assert sha(SOURCE)=='b796be0506a0d696ee19aad9288e2b71b6c81131c30f7e4bf2deb00c94ddb80a'
names={'norm','load_pdb','trees','clear','make_triads','sample_donor','compatibility','seed_for'}
module=ast.parse(SOURCE.read_text())
functions=[node for node in module.body if isinstance(node,ast.FunctionDef) and node.name in names]
assert {f.name for f in functions}==names
exec(compile(ast.Module(body=functions,type_ignores=[]),str(SOURCE),'exec'),globals())
templatepath=Path('/data/bht2/new_seed_structure_sequence_folddisco_20260825_v3/structure_collection_20260826_v1/structures/folddisco/PDB/PDB_3AXG_A.pdb')
rec=load_pdb(templatepath,True)
tri=[v for v in rec if int(v[0][22:26]) in (267,306,308)]
keys=[(int(l[22:26]),l[12:16].strip()) for l,_,_ in tri]
assert len(tri)==23 and len(set(keys))==23
og=keys.index((267,'OG1'));cb=keys.index((267,'CB'))
tq=np.array([q for _,q,_ in tri]);tq-=tq[og]
te=np.array([e for _,_,e in tri]);template=(tq,te,og,cb,[l for l,_,_ in tri])
aa={(int(l[22:26]),l[12:16].strip()):(l,q,e) for l,q,e in load_pdb(PET/'catalytic_motif_compare_v1/5XJH.pdb',True)}
dk=[(86,'C'),(86,'O'),(87,'N'),(87,'CA'),(87,'C'),(87,'O')]
dq=np.array([aa[k][1] for k in dk]);dq-=dq[2];de=np.array([aa[k][2] for k in dk])
z=-norm(norm(dq[0])+norm(dq[3]));x=norm(np.cross(z,[1.,0,0]) if abs(z[0])<.9 else np.cross(z,[0.,1,0]));y=np.cross(z,x)
donor_templates=[(dq,np.column_stack([x,y,z]),de,[aa[k][0] for k in dk])]
rr=np.array([RAD[e] for e in de])
def evaluate(mat,row,out):
 global rng
 rec=load_pdb(ROOT/'inputs'/f'{mat}_v10_ge90_uncapped.pdb')
 lookup={l[6:11]:q for l,q,_ in rec}
 amide=np.array([lookup[k] for k in row['serials']])
 c=amide[0];ex=norm(amide[1]-c);ey=amide[2]-c;ey=norm(ey-ex*np.dot(ex,ey));frame=np.column_stack((ex,ey,np.cross(ex,ey)))
 xyz=np.array([q for _,q,_ in rec]);el=np.array([e for _,_,e in rec])
 # Full finite environment, no local cutoff.
 env=trees(xyz,el)
 key=mat+'_'+row['site_id']
 rng=np.random.default_rng(seed_for(key+'|TDD|triad'))
 poses,stats=make_triads(amide,frame,template,env)
 q,H,*_=sample_donor(amide,frame,seed_for(key+'|NH|Sobol'))
 good=np.flatnonzero(clear(q,rr,env)>=0)
 cand=q[good];adj=compatibility(cand)
 pairs=[];accepted=[]
 for j,pose in enumerate(poses):
  valid=np.flatnonzero(clear(cand,rr,trees(pose,te))>=0) if len(cand) else np.array([],int)
  ij=np.argwhere(np.triu(adj[np.ix_(valid,valid)],1))
  accepted.append(good[valid])
  if len(ij):pairs.extend(np.column_stack((np.full(len(ij),j),good[valid[ij[:,0]]],good[valid[ij[:,1]]])).tolist())
 pairarr=np.array(pairs,dtype=int).reshape(-1,3)
 np.savez_compressed(out/(key+'.npz'),triads_A=poses,donors_A=q,donor_H_A=H,pairs=pairarr,amide_A=amide,frame=frame,origin_A=c,polymer_clear_donor_ids=good,**{f'pose_{j}_donor_ids':ids for j,ids in enumerate(accepted)})
 result=dict(material=mat,site=row['site_id'],support=row['inward_support'],side=row['nearest_side'],parent_outward_points=row['parent_outward_points'],triad_stats=stats,donor_trials=len(q),polymer_clear_donors=len(good),paired_triad_poses=len(set(pairarr[:,0].tolist())),pair_records=len(pairs),status='GEOMETRY_PROXY_ONLY')
 print(json.dumps(result),flush=True)
 return result
if __name__=='__main__':
 out=ROOT/'tdd_donor_pilot_v1';out.mkdir(exist_ok=False);rc=1
 config=dict(source_sha256=sha(SOURCE),template_sha256=sha(templatepath),core_keys=keys,donor_template='IsPETase Tyr87 peptide as generic backbone NH; not native NylC donors',triad_trials=12000,max_retained_poses=64,donor_samples=8192,attack_distance_A=[3,3.4],O_C_OG1_deg=[95,115],CB_OG1_C_deg=[100,130],N_O_A=[2.7,3.2],N_H_O_min_deg=140,NH_A=1.01,heavy_atom_overlap_A=.4,scope='reactant static geometry proxy; not validated NylC NAC; native Lys network and native donor chemistry not modeled')
 (out/'config.json').write_text(json.dumps(config,indent=2))
 try:
  results=[]
  for mat in ('PA6','PA66'):
   rows=json.loads((ROOT/'carbon_exposure_v1'/f'{mat}_sites.json').read_text())
   picked=[]
   for side in ('top','bottom'):
    for parent_visible in (True,False):
     subset=[r for r in rows if r['nearest_side']==side and (r['parent_outward_points']>0)==parent_visible]
     if subset:picked.append(sorted(subset,key=lambda r:(-r['inward_support'],-r['trimmed_outward_points'],r['site_id']))[0])
   for row in picked:results.append(evaluate(mat,row,out))
  (out/'summary.json').write_text(json.dumps(results,indent=2));rc=0
 finally:
  with (ROOT/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__).resolve()),script_sha256=sha(Path(__file__)),command='CPU Python tdd_donor_pilot_v1.py',inputs='v10 PDBs, exposure candidates, 3AXG, 5XJH, PET helper source',outputs=str(out),parameters=config,exit_code=rc))+'\n')
