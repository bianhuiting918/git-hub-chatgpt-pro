#!/usr/bin/env python3
"""Small equal-budget full-residue conformer library for dual-donor feasibility."""
from pathlib import Path
import json,hashlib,datetime
import numpy as np
from rdkit import Chem,rdBase
from rdkit.Chem import AllChem
ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
if __name__=='__main__':
 out=ROOT/'conformers_v1';out.mkdir(exist_ok=False);rc=1
 try:
  source=json.loads((ROOT/'chemical_eligibility.json').read_text());report=[]
  for ri,row in enumerate(source['results']):
   if row['status']!='ELIGIBLE_FOR_GEOMETRY_SCAN':continue
   key=row['amino_acid']+'_'+str(ri)
   m=Chem.MolFromSmiles(row['smiles'])
   bb=next(a.GetIdx() for a in m.GetAtoms() if a.GetAtomMapNum()==1)
   ca=next(a.GetIdx() for a in m.GetAtoms() if a.GetAtomMapNum()==2)
   cap=next(a.GetIdx() for a in m.GetAtoms() if a.GetAtomMapNum()==5)
   # Atom maps can affect RDKit CIP ranking; physical chirality is checked without maps.
   for a in m.GetAtoms():a.SetAtomMapNum(0)
   Chem.AssignStereochemistry(m,cleanIt=True,force=True)
   expected='R' if row['amino_acid']=='CYS' else 'S'
   assert m.GetAtomWithIdx(ca).GetProp('_CIPCode')==expected
   h=Chem.AddHs(m);donors=[]
   for a in h.GetAtoms():
    hs=[n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum()==1]
    if a.GetIdx()!=cap and a.GetSymbol() in ('N','O','S') and hs:
     donors.append(dict(atom=a.GetIdx(),element=a.GetSymbol(),hydrogens=hs,role='backbone' if a.GetIdx()==bb else 'sidechain'))
   assert len(donors)==row['n_distinct_donor_atoms']
   seed=int.from_bytes(hashlib.sha256(key.encode()).digest()[:4],'little')%2147483647
   ep=AllChem.ETKDGv3();ep.randomSeed=seed;ep.numThreads=1;ep.pruneRmsThresh=-1;ep.enforceChirality=True
   ids=list(AllChem.EmbedMultipleConfs(h,numConfs=32,params=ep))
   assert ids, key
   assert AllChem.MMFFHasAllMoleculeParams(h),key
   opt=AllChem.MMFFOptimizeMoleculeConfs(h,numThreads=1,maxIters=500,mmffVariant='MMFF94s')
   positions=[];convergence=[];energies=[]
   for cid,(status,energy) in zip(ids,opt):
    xyz=np.array(h.GetConformer(cid).GetPositions());assert np.isfinite(xyz).all()
    positions.append(xyz);convergence.append(status);energies.append(energy)
   heavy=[a.GetIdx() for a in h.GetAtoms() if a.GetAtomicNum()!=1]
   elements=[a.GetSymbol() for a in h.GetAtoms()]
   np.savez_compressed(out/(key+'.npz'),xyz_A=np.array(positions),heavy_indices=np.array(heavy),elements=np.array(elements),optimization_status=np.array(convergence),MMFF_energy_kcal_mol=np.array(energies))
   writer=Chem.SDWriter(str(out/(key+'.sdf')))
   for cid in ids:writer.write(h,confId=cid)
   writer.close()
   metadata=dict(key=key,amino_acid=row['amino_acid'],microstate=row['microstate'],seed=seed,requested_conformers=32,embedded=len(ids),converged=sum(s==0 for s in convergence),donors=donors,backbone_N=bb,alpha_C=ca,cap_N_excluded=cap,heavy_atoms=len(heavy),CIP_unmapped=expected,smiles=Chem.MolToSmiles(m,isomericSmiles=True),status='CONFORMERS_READY' if all(s==0 for s in convergence) else 'PARTIAL_OPTIMIZATION_CONVERGENCE')
   (out/(key+'.json')).write_text(json.dumps(metadata,indent=2));report.append(metadata)
   print(json.dumps({k:metadata[k] for k in ['key','embedded','converged','heavy_atoms','status']}),flush=True)
  (out/'summary.json').write_text(json.dumps(dict(rdkit=rdBase.rdkitVersion,source_sha256=sha(ROOT/'chemical_eligibility.json'),records=report,scope='32 equal requested ETKDGv3 conformers per eligible microstate; MMFF94s geometry cleanup only; energy is not binding or catalytic score; not conformational convergence; donor indices regenerated after canonical SMILES parsing'),indent=2));rc=0
 finally:
  with (ROOT/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__).resolve()),script_sha256=sha(Path(__file__)),command='CPU Python build_conformers_v1.py',inputs='chemical_eligibility.json',outputs=str(out),parameters=dict(requested_conformers=32,MMFF_max_iterations=500,threads=1),exit_code=rc))+'\n')
