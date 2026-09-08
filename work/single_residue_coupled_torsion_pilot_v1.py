#!/usr/bin/env python3
"""Coupled single-residue torsion pilot; geometric feasibility, not energy ranking."""
from pathlib import Path
import json,hashlib,datetime,importlib.util
import numpy as np
from scipy.stats import qmc
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
R=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('intrinsic',R/'intrinsic_dual_oxygen_v1.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
OUT=R/'coupled_torsion_pilot_v1'
def torsions(m):
 result=[]
 for b in m.GetBonds():
  if b.GetBondType()!=Chem.BondType.SINGLE or b.IsInRing() or b.GetIsConjugated():continue
  j,k=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
  if min(m.GetAtomWithIdx(j).GetAtomicNum(),m.GetAtomWithIdx(k).GetAtomicNum())==1:continue
  left=[a for a in m.GetAtomWithIdx(j).GetNeighbors() if a.GetIdx()!=k]
  right=[a for a in m.GetAtomWithIdx(k).GetNeighbors() if a.GetIdx()!=j]
  if not left or not right:continue
  # Ignore terminal methyl spinning; retain heteroatom-bound hydrogen rotation.
  if not any(a.GetAtomicNum()>1 for a in left) and m.GetAtomWithIdx(j).GetAtomicNum()==6:continue
  if not any(a.GetAtomicNum()>1 for a in right) and m.GetAtomWithIdx(k).GetAtomicNum()==6:continue
  i=max(left,key=lambda a:a.GetAtomicNum()).GetIdx()
  l=max(right,key=lambda a:a.GetAtomicNum()).GetIdx()
  result.append([i,j,k,l])
 return result
if __name__=='__main__':
 OUT.mkdir(exist_ok=False);rc=1
 try:
  report=[]
  for meta in json.loads((R/'conformers_v1/summary.json').read_text())['records']:
   key=meta['key'];z=np.load(R/'conformers_v1'/(key+'.npz'))
   mols=list(Chem.SDMolSupplier(str(R/'conformers_v1'/(key+'.sdf')),removeHs=False))
   ts=torsions(mols[0]);seed=int.from_bytes(hashlib.sha256(key.encode()).digest()[:4],'little')
   u=qmc.Sobol(max(1,len(ts)),scramble=True,seed=seed).random_base2(8)
   heavy=z['heavy_indices'];rr=np.array([mod.RAD[str(z['elements'][i])] for i in heavy])
   mask=np.triu(Chem.GetDistanceMatrix(mols[0])[np.ix_(heavy,heavy)]>3,1)
   good=np.flatnonzero(z['optimization_status']==0)
   points=[];pids=[];pairs=[];geometries=[];records=[];intraclear=0;hitposes=0
   for pi,rot in enumerate(u):
    ci=int(good[pi%len(good)]);mol=Chem.Mol(mols[ci]);cf=mol.GetConformer()
    for ai,v in enumerate(z['xyz_A'][ci]):cf.SetAtomPosition(ai,v.tolist())
    for t,a in zip(ts,rot):rdMolTransforms.SetDihedralDeg(cf,*t,float(360*a-180))
    x=np.array(cf.GetPositions())
    # Rotations must preserve all covalent bond lengths.
    for b in mol.GetBonds():
     i,j=b.GetBeginAtomIdx(),b.GetEndAtomIdx()
     assert abs(np.linalg.norm(x[i]-x[j])-np.linalg.norm(z['xyz_A'][ci,i]-z['xyz_A'][ci,j]))<1e-6
    margin=np.linalg.norm(x[heavy][:,None]-x[heavy][None],axis=-1)-rr[:,None]-rr[None]+.4
    if not (margin[mask]>=-1e-8).all():continue
    intraclear+=1;gi=len(geometries);geometries.append(x)
    before=len(points)
    for di,d in enumerate(meta['donors']):
     uv=qmc.Sobol(3,scramble=True,seed=(seed+pi*17+di)%2**32).random_base2(13)
     c=2*uv[:,0]-1;p=2*np.pi*uv[:,1]
     lo,hi=(3.1,3.6) if d['element']=='S' else (2.7,3.2)
     pos=x[d['atom']]+(lo+(hi-lo)*uv[:,2,None])*np.column_stack((c,np.sqrt(1-c*c)*np.cos(p),np.sqrt(1-c*c)*np.sin(p)))
     valid=mod.valid(pos,x,d)
     for dj in range(di+1,len(meta['donors'])):
      a=pos[valid&mod.valid(pos,x,meta['donors'][dj])]
      ok=((np.linalg.norm(a[:,None]-x[heavy][None],axis=-1)-rr[None]-1.52+.4)>=-1e-8).all(1)
      a=a[ok];points.extend(a.tolist());pids.extend([gi]*len(a));pairs.extend([[di,dj]]*len(a))
    found=len(points)-before;hitposes+=int(found>0)
    records.append(dict(geometry_index=gi,trial_index=pi,original_conformer=ci,torsion_degrees=(360*rot[:len(ts)]-180).tolist(),accepted_samples=found))
   np.savez_compressed(OUT/(key+'.npz'),residue_A=np.array(geometries).reshape(-1,z['xyz_A'].shape[1],3),oxygen_A=np.array(points).reshape(-1,3),geometry_indices=np.array(pids,dtype=int),donor_pair_indices=np.array(pairs,dtype=int).reshape(-1,2))
   rec=dict(key=key,amino_acid=meta['amino_acid'],microstate=meta['microstate'],torsions_atom_indices=ts,requested_torsion_trials=256,internal_heavy_clear_poses=intraclear,common_oxygen_poses=hitposes,common_oxygen_samples=len(points),unconverged_source_excluded=int((z['optimization_status']!=0).sum()),pose_records=records)
   (OUT/(key+'.json')).write_text(json.dumps(rec,indent=2));report.append({k:v for k,v in rec.items() if k!='pose_records'})
   print(json.dumps(report[-1]),flush=True)
  # Regression: the new library must no longer lose all Thr dual-donor geometries.
  assert next(r for r in report if r['key']=='THR_16')['common_oxygen_poses']>0,'Thr torsion regression failed'
  (OUT/'summary.json').write_text(json.dumps(dict(results=report,parameters=dict(torsion_trials_per_microstate=256,oxygen_samples_per_anchor=8192,D_O_A=[2.7,3.2],S_O_A=[3.1,3.6],angle_min_deg=140,overlap_A=.4),scope='unrelaxed covalently coupled torsion pilot; not energy-ranked, not polymer-tested; finite zeros not impossibility; equal requested trial budget not convergence',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2));rc=0
 finally:
  with (R/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__).resolve()),inputs='conformers_v1',outputs=str(OUT),exit_code=rc,parameters='256 coupled torsions per microstate; 8192 oxygen samples per anchor; unchanged geometry criteria'))+'\n')
