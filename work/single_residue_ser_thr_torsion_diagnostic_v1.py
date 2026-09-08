#!/usr/bin/env python3
"""Controlled Ser/Thr torsion diagnostic; unchanged common-oxygen criteria."""
from pathlib import Path
import json,hashlib,datetime,importlib.util
import numpy as np
from scipy.stats import qmc
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
R=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('intrinsic',R/'intrinsic_dual_oxygen_v1.py')
mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
out=R/'ser_thr_torsion_diagnostic_v1';out.mkdir(exist_ok=False)
rc=1
try:
 results=[]
 for key in ['SER_15','THR_16']:
  meta=json.loads((R/'conformers_v1'/(key+'.json')).read_text())
  z=np.load(R/'conformers_v1'/(key+'.npz'))
  mols=list(Chem.SDMolSupplier(str(R/'conformers_v1'/(key+'.sdf')),removeHs=False))
  d,e=meta['donors'];n=d['atom'];o=e['atom'];h=e['hydrogens'][0];ca=meta['alpha_C']
  mol=mols[0]
  assert mol.GetAtomWithIdx(n).GetSymbol()=='N' and mol.GetAtomWithIdx(o).GetSymbol()=='O'
  cb=next(a.GetIdx() for a in mol.GetAtomWithIdx(o).GetNeighbors() if a.GetAtomicNum()==6)
  assert mol.GetBondBetweenAtoms(ca,cb) is not None
  heavy=z['heavy_indices'];rad=np.array([mod.RAD[str(z['elements'][j])] for j in heavy])
  # Independently check the angle implementation on a collinear D-H...O test.
  fake=np.array([[0.,0.,0.],[1.,0.,0.]])
  assert mod.valid(np.array([[3.,0.,0.]]),fake,dict(atom=0,element='N',hydrogens=[1]))[0]
  assert not mod.valid(np.array([[-3.,0.,0.]]),fake,dict(atom=0,element='N',hydrogens=[1]))[0]
  bonded=Chem.GetDistanceMatrix(mol)[np.ix_(heavy,heavy)]
  nonlocal_mask=np.triu(bonded>3,1)
  u=qmc.Sobol(3,scramble=True,seed=711).random_base2(13)
  c=2*u[:,0]-1;p=2*np.pi*u[:,1]
  dirs=np.column_stack((c,np.sqrt(1-c*c)*np.cos(p),np.sqrt(1-c*c)*np.sin(p)))
  rows=[];examples=[]
  for ci in [0,1]:
   original=z['xyz_A'][ci]
   for mode in ['original','H_only','chi_only','chi_and_H']:
    chis=[None] if mode in ['original','H_only'] else list(range(-180,180,30))
    hs=[None] if mode in ['original','chi_only'] else list(range(-180,180,30))
    counts=np.zeros(4,int);poses=0
    for chi in chis:
     for ht in hs:
      m=Chem.Mol(mols[ci]);conf=m.GetConformer()
      # Load full-precision original coordinates rather than rounded SDF coordinates.
      for ai,v in enumerate(original):conf.SetAtomPosition(ai,v.tolist())
      if chi is not None:rdMolTransforms.SetDihedralDeg(conf,n,ca,cb,o,float(chi))
      if ht is not None:rdMolTransforms.SetDihedralDeg(conf,ca,cb,o,h,float(ht))
      x=np.array(conf.GetPositions());poses+=1
      delta=np.linalg.norm(x[heavy][:,None]-x[heavy][None],axis=-1)-rad[:,None]-rad[None]+.4
      internal_ok=bool((delta[nonlocal_mask]>=-1e-8).all())
      pts=x[n]+(2.7+.5*u[:,2,None])*dirs
      rr=np.linalg.norm(pts-x[o],axis=1);dist=(rr>=2.7)&(rr<=3.2);counts[0]+=int(dist.sum())
      ang=dist&mod.valid(pts,x,d)&mod.valid(pts,x,e);counts[1]+=int(ang.sum())
      cand=pts[ang]
      clear=((np.linalg.norm(cand[:,None]-x[heavy][None],axis=-1)-rad[None]-1.52+.4)>=-1e-8).all(1)
      accepted=cand[clear];counts[2]+=len(accepted)
      if internal_ok:
       counts[3]+=len(accepted)
       if len(accepted) and len(examples)<8:
        examples.append(dict(conformer=ci,mode=mode,chi_deg=chi,H_torsion_deg=ht,oxygen_A=accepted[0].tolist(),residue_A=x.tolist()))
    rows.append(dict(conformer=ci,mode=mode,poses=poses,distance_samples=int(counts[0]),angle_samples=int(counts[1]),oxygen_clear_samples=int(counts[2]),internal_heavy_clear_samples=int(counts[3])))
   print(key,ci,rows[-4:],flush=True)
  rec=dict(key=key,rows=rows,examples=examples)
  results.append(rec)
 report=dict(results=results,scope='diagnostic only, two original conformers each; torsions unrelaxed; no material or activity ranking',criteria=dict(D_O_A=[2.7,3.2],D_H_O_min_deg=140,vdw_overlap_A=.4,internal_heavy_exclusion='graph distance >3'),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 (out/'summary.json').write_text(json.dumps(report,indent=2));rc=0
finally:
 with (R/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__).resolve()),outputs=str(out),exit_code=rc,parameters='two conformers Ser/Thr; chi and OH torsions 30deg; original criteria'))+'\n')
