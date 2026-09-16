"""Stage-decomposed geometric diagnostic, not a catalytic or material PASS.
Run with the existing CPU Python from the N9 design project directory.
Unchanged donor/attack thresholds; no structure edits. Output JSON and RUN_LOG.
"""
from pathlib import Path
import numpy as np,json,hashlib,datetime
from scipy.stats import qmc
from route56_N10_20260915_v1.scan_N10_mature_donors_v1 import read_atoms,unit
r=Path(__file__).resolve().parent
d=r/'route56_N10_20260915_v1'
src=Path('/data/bht2/new_seed_structure_sequence_folddisco_20260825_v3/n9_n10_af_repredict_20260908_v1/outputs/nonfine_af_v1/N10_AF2fresh/seed_101/mature_relaxed_nh3.pdb')
assert hashlib.sha256(src.read_bytes()).hexdigest()=='c110940f3342b55c0ca62516df624d35d4beb4357035a17af65fa387125d9433'
a=read_atoms(src);x=np.array([v['xyz'] for v in a])
don=[v for v in json.loads((d/'mature_donor_scan_seed101_v1/donors.json').read_text()) if v['id'] in ['PHE135_N','ASN205_ND2']]
# Match Thr atoms directly from PDB labels, independent of helper dictionary keys.
raw=[l for l in src.read_text().splitlines() if l.startswith('ATOM')]
def coord(name):
 vals=[np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]) for l in raw if int(l[22:26])==258 and l[12:16].strip()==name]
 assert len(vals)==1
 return vals[0]
og=coord('OG1');cb=unit(coord('CB')-og)
def attack(points):
 v=points-og;delta=np.linalg.norm(v,axis=1);axis=unit(v);dot=axis@cb
 triangle=np.zeros(len(points),bool);full=triangle.copy()
 for radius in np.linspace(3.,3.4,401):
  cosine=(radius**2+1.23**2-delta**2)/(2*radius*1.23)
  ok=(cosine>=np.cos(np.deg2rad(115)))&(cosine<=np.cos(np.deg2rad(95)))
  t=(radius**2+delta**2-1.23**2)/(2*delta)
  h2=radius**2-t*t;ok &= h2>=0
  middle=t/radius*dot
  spread=np.sqrt(np.maximum(h2,0))/radius*np.sqrt(np.maximum(1-dot*dot,0))
  # Continuous azimuth feasibility: intersect attainable cosine interval with angle gate.
  good=ok&(middle-spread<=np.cos(np.deg2rad(100)))&(middle+spread>=np.cos(np.deg2rad(130)))
  triangle|=ok;full|=good
 return triangle,full
reports=[]
for seed in [101,202]:
 u=qmc.Sobol(3,scramble=True,seed=seed).random_base2(18)
 z=2*u[:,0]-1;p=2*np.pi*u[:,1];rad=2.7+.5*u[:,2]
 points=x[don[0]['atom_index']]+rad[:,None]*np.c_[np.sqrt(1-z*z)*np.cos(p),np.sqrt(1-z*z)*np.sin(p),z]
 D=[];A=[]
 for v in don:
  n=x[v['atom_index']];D.append(np.linalg.norm(points-n,axis=1))
  A.append(np.max(np.stack([np.degrees(np.arccos(np.clip(unit(points-x[h])@unit(n-x[h]),-1,1))) for h in v['hydrogen_indices']]),axis=0))
 D=np.array(D).T;A=np.array(A).T
 keep=((D>=2.7)&(D<=3.2)).all(1);p=points[keep];aa=A[keep]
 tri,atk=attack(p)
 masks={'distance_only':np.ones(len(p),bool),'PHE_H_direction':aa[:,0]>=140,'ASN_H_direction':aa[:,1]>=140,'both_H_directions':(aa>=140).all(1)}
 reports.append(dict(seed=seed,sampled_oxygen_points=len(points),stages={k:dict(donor_points=int(m.sum()),with_attack_triangle=int((m&tri).sum()),with_all_attack_angles=int((m&atk).sum())) for k,m in masks.items()}))
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
out=d/('native_N10_stage_decomposition_'+stamp+'.json')
report=dict(status='DIAGNOSTIC_NOT_CATALYTIC_PASS',input=str(src),input_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),donors=[v['id'] for v in don],reports=reports,protocol=dict(seeds=[101,202],oxygen_samples_per_seed=262144,oxygen_radial_sampling='uniform radius, not uniform shell volume',OG_C_grid_A=[3,3.4,401],CO_A=1.23,O_C_OG_deg=[95,115],CB_OG_C_deg=[100,130],donor_O_A=[2.7,3.2],D_H_O_min_deg=140,azimuth='analytic continuous interval'),limitations=['Counts are oxygen points, not polymer sites or probabilities','Finite oxygen and radial sampling cannot prove global nonexistence','No material or protein steric checks','Unchanged static hydrogen and sidechain positions','Carbonyl probe is preattack, not tetrahedral intermediate'])
with out.open('x') as h:json.dump(report,h,indent=2)
with (r/'RUN_LOG.jsonl').open('a') as h:h.write(json.dumps(dict(time=stamp,script=str(Path(__file__)),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input=str(src),output=str(out),exit_code=0))+'\n')
print(json.dumps(report,indent=2));print('OUTPUT',out)
