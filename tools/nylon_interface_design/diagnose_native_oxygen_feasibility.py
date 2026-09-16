"""Diagnostic oxygen-only feasibility; immutable native N10, no material or catalytic PASS."""
from pathlib import Path
import numpy as np,json,hashlib,datetime
from scipy.stats import qmc
from route56_N10_20260915_v1.scan_N10_mature_donors_v1 import read_atoms,unit
r=Path(__file__).resolve().parent;d=r/'route56_N10_20260915_v1/mature_donor_scan_seed101_v1'
src=Path('/data/bht2/new_seed_structure_sequence_folddisco_20260825_v3/n9_n10_af_repredict_20260908_v1/outputs/nonfine_af_v1/N10_AF2fresh/seed_101/mature_relaxed_nh3.pdb')
a=read_atoms(src);x=np.array([v['xyz'] for v in a]);don=[v for v in json.loads((d/'donors.json').read_text()) if v['id'] in ['PHE135_N','ASN205_ND2']]
assert len(don)==2
reports=[]
for seed in [101,202]:
 u=qmc.Sobol(3,scramble=True,seed=seed).random_base2(18)
 z=2*u[:,0]-1;phi=2*np.pi*u[:,1];rad=2.7+.5*u[:,2]
 points=x[don[0]['atom_index']]+rad[:,None]*np.c_[np.sqrt(1-z*z)*np.cos(phi),np.sqrt(1-z*z)*np.sin(phi),z]
 distances=[];angles=[]
 for v in don:
  n=x[v['atom_index']];distances.append(np.linalg.norm(points-n,axis=1))
  angles.append(np.max(np.stack([np.degrees(np.arccos(np.clip(unit(points-x[h])@unit(n-x[h]),-1,1))) for h in v['hydrogen_indices']]),axis=0))
 D=np.array(distances).T;A=np.array(angles).T;dist=((D>=2.7)&(D<=3.2)).all(1);joint=dist&(A>=140).all(1)
 k=np.flatnonzero(dist)[np.argmax(A[dist].min(1))] if dist.any() else None
 row=dict(seed=seed,samples=len(points),both_distance=int(dist.sum()),both_distance_and_direction=int(joint.sum()),best_minimum_angle_with_both_distances=float(A[k].min()) if k is not None else None,
 best_example=dict(oxygen_A=points[k].tolist(),distances_A=D[k].tolist(),angles_deg=A[k].tolist()) if k is not None else None)
 reports.append(row)
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');out=r/'route56_N10_20260915_v1'/('native_N10_oxygen_only_'+stamp+'.json')
report=dict(status='OXYGEN_ONLY_DIAGNOSTIC_NOT_CATALYTIC_PASS',input=str(src),input_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),donors=don,reports=reports,
 limitations=['No protein motion, no Asn flip, no hydrogen reoptimization','No material collisions, attack carbon geometry, carbonyl orientation, or steric checks','Finite sampling cannot prove nonexistence','Counts refer to sampled oxygen points, not polymer sites'])
with out.open('x') as h:json.dump(report,h,indent=2)
with (r/'RUN_LOG.jsonl').open('a') as h:h.write(json.dumps(dict(time=stamp,script=str(Path(__file__)),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),output=str(out),exit_code=0))+'\n')
print(json.dumps(dict(output=str(out),reports=reports),indent=2))
