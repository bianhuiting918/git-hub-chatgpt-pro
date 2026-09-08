#!/usr/bin/env python3
"""Intrinsic common-oxygen feasibility of coupled donor groups; no material yet."""
from pathlib import Path
import json,hashlib,datetime,itertools
import numpy as np
from scipy.stats import qmc
ROOT=Path(__file__).resolve().parent
RAD={'C':1.7,'N':1.55,'O':1.52,'S':1.8}
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def valid(points,xyz,d):
 D=xyz[d['atom']];r=np.linalg.norm(points-D,axis=1)
 lo,hi=(3.1,3.6) if d['element']=='S' else (2.7,3.2)
 ok=(r>=lo-1e-8)&(r<=hi+1e-8);best=np.zeros(len(points),bool)
 for h in d['hydrogens']:
  a=D-xyz[h];b=points-xyz[h];co=(b@a)/(np.linalg.norm(a)*np.linalg.norm(b,axis=1))
  best|=co<=np.cos(np.deg2rad(140))+1e-12
 return ok&best
if __name__=='__main__':
 out=ROOT/'intrinsic_dual_oxygen_v1';out.mkdir(exist_ok=False);rc=1
 try:
  info=json.loads((ROOT/'conformers_v1/summary.json').read_text());summary=[]
  for meta in info['records']:
   key=meta['key'];z=np.load(ROOT/'conformers_v1'/(key+'.npz'));donors=meta['donors']
   heavy=z['heavy_indices'];rr=np.array([RAD[str(z['elements'][i])] for i in heavy])
   allpoints=[];allconf=[];allpairs=[];counts=[]
   for ci,xyz in enumerate(z['xyz_A']):
    if z['optimization_status'][ci]!=0:counts.append(None);continue
    number=0
    for di,d in enumerate(donors):
     seed=int.from_bytes(hashlib.sha256(f'{key}|{ci}|{di}'.encode()).digest()[:4],'little')
     u=qmc.Sobol(3,scramble=True,seed=seed).random_base2(13)
     cos=2*u[:,0]-1;phi=2*np.pi*u[:,1];lo,hi=(3.1,3.6) if d['element']=='S' else (2.7,3.2)
     directions=np.column_stack((cos,np.sqrt(1-cos*cos)*np.cos(phi),np.sqrt(1-cos*cos)*np.sin(phi)))
     points=xyz[d['atom']]+(lo+(hi-lo)*u[:,2,None])*directions
     eligible=valid(points,xyz,d)
     for dj in range(di+1,len(donors)):
      ids=np.flatnonzero(eligible&valid(points,xyz,donors[dj]))
      if not len(ids):continue
      candidates=points[ids]
      margin=np.linalg.norm(candidates[:,None]-xyz[heavy][None],axis=-1)-rr[None]-1.52+.4
      accepted=candidates[(margin>=-1e-8).all(1)]
      number+=len(accepted)
      allpoints.extend(accepted.tolist());allconf.extend([ci]*len(accepted));allpairs.extend([[di,dj]]*len(accepted))
    counts.append(number)
   np.savez_compressed(out/(key+'.npz'),oxygen_A=np.array(allpoints).reshape(-1,3),conformer_indices=np.array(allconf,dtype=int),donor_pair_indices=np.array(allpairs,dtype=int).reshape(-1,2))
   r=dict(key=key,amino_acid=meta['amino_acid'],microstate=meta['microstate'],tested_conformers=sum(n is not None for n in counts),not_evaluated_conformers=sum(n is None for n in counts),conformers_with_common_oxygen=sum(n is not None and n>0 for n in counts),accepted_oxygen_pair_samples=len(allpoints),per_conformer=counts,status='COMMON_OXYGEN_FOUND' if allpoints else 'NONE_IN_THIS_FINITE_LIBRARY_NOT_ABSENCE_PROOF')
   summary.append(r);print(json.dumps({k:v for k,v in r.items() if k!='per_conformer'}),flush=True)
  report=dict(results=summary,parameters=dict(samples_per_anchor_conformer=8192,N_O_and_O_O_A=[2.7,3.2],S_O_A=[3.1,3.6],D_H_O_min_deg=140,oxygen_vdw_A=1.52,overlap_allowance_A=.4),scope='intrinsic same-residue two-distinct-donor geometry only; no PET/nylon/catalytic core; raw point counts not comparable as binding or catalytic score; finite conformer and fixed X-H orientations')
  (out/'summary.json').write_text(json.dumps(report,indent=2));rc=0
 finally:
  with (ROOT/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__).resolve()),script_sha256=sha(Path(__file__)),command='CPU Python intrinsic_dual_oxygen_v1.py',inputs='conformers_v1 converged geometries',outputs=str(out),parameters=dict(sobol_anchor_samples=8192,angle_min=140),exit_code=rc))+'\n')
