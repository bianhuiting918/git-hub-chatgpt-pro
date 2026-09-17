"""Snapshot audit of completed top19 three-seed backbone outputs.
No repairs or resubmissions. Backbone-only metrics are not catalytic PASS.
"""
import sys,json,datetime,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial.distance import cdist
from scipy.spatial import cKDTree
r=Path(__file__).resolve().parent;d=r/'route56_N10_20260915_v1'
base=d/'route6_top50_local_generation_v1';g=r.parent/'n9_core_donor_design_20260910_gpu_outputs_v1/route6_top50_local_shapes_3seeds_v2'
bb=['N','CA','C','O'];mapping=list(zip(range(1,10),[1,40,41,42,61,62,63,64,65]))
sites=json.loads((d/'N10_fixed_minimal_motif_material_v1/sites.json').read_text())
def read(p):
 out={}
 for line in p.read_text().splitlines():
  if line.startswith('ATOM  ') and line[12:16].strip() in bb:
   i=int(line[22:26]);n=line[12:16].strip();assert n not in out.setdefault(i,{})
   out[i][n]=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])
 return out
def fit(a,b):
 u,s,v=np.linalg.svd((a-a.mean(0)).T@(b-b.mean(0)));rot=u@np.diag([1,1,np.linalg.det(u@v)])@v
 shift=b.mean(0)-a.mean(0)@rot;delta=a@rot+shift-b
 return rot,shift,dict(RMSD_A=float(np.sqrt(np.mean(np.sum(delta**2,axis=1)))),max_A=float(np.linalg.norm(delta,axis=1).max()))
def tors(a,b,c,d):
 u=b-a;w=c-b;t=d-c;n=np.cross(u,w);m=np.cross(w,t)
 return float(np.degrees(np.arctan2(np.dot(np.cross(n,m),w/np.linalg.norm(w)),np.dot(n,m))))
rad={'C':1.7,'N':1.55,'O':1.52,'S':1.8};envs={}
for mat in ['PA6','PA66']:
 z=np.load(r.parent/'interface_surface_robustness_20260908_v1/inputs_v2'/(mat+'_parent_atoms.npz'))
 px=z['xyz_A'].copy();box=z['box_A'];px[:,:2]%=box[:2]
 envs[mat]=(box,{e:cKDTree(px[z['elements']==e],boxsize=[*box[:2],0]) for e in rad if np.any(z['elements']==e)})
reports=[]
for path in sorted(g.glob('*/design_*.pdb')):
 case=path.parent.name;cid=case.split('_')[1]
 j=json.loads((base/'candidate_specific_shapes_v2'/cid/case/'manifest.json').read_text())
 out=read(path);assert sorted(out)==list(range(1,127)) and all(set(v)==set(bb) for v in out.values())
 inp=read(Path(j['pdb']));a=np.array([out[k][atom] for i,k in mapping for atom in bb]);b=np.array([inp[i][atom] for i,k in mapping for atom in bb])
 rot,shift,metric=fit(a,b)
 xyz=np.array([out[i][atom] for i in range(1,127) for atom in bb]);rr=np.tile([1.55,1.7,1.7,1.52],126)
 cn=[float(np.linalg.norm(out[i]['C']-out[i+1]['N'])) for i in range(1,126)]
 omega=[180-abs(tors(out[i]['CA'],out[i]['C'],out[i+1]['N'],out[i+1]['CA'])) for i in range(1,126)]
 adjacency=[set() for _ in xyz]
 for i in range(126):
  for a,b in [(4*i,4*i+1),(4*i+1,4*i+2),(4*i+2,4*i+3)]+([(4*i+2,4*i+4)] if i<125 else []):
   adjacency[a].add(b);adjacency[b].add(a)
 margin=cdist(xyz,xyz)-rr[:,None]-rr[None,:]+.4
 for i in range(len(xyz)):
  reach={i};front={i}
  for _ in range(3):front={b for a in front for b in adjacency[a]}-reach;reach|=front
  margin[i,list(reach)]=np.inf
 original=xyz@rot+shift;hits={'PA6':[],'PA66':[]}
 for si,s in enumerate(sites):
  box,trees=envs[s['material']];world=original@np.array(s['frame']).T+s['center'];world[:,:2]%=box[:2];clear=np.full(len(world),np.inf)
  for e,t in trees.items():clear=np.minimum(clear,t.query(world,workers=1)[0]-rr-rad[e]+.4)
  if clear.min()>=-1e-8:hits[s['material']].append(si)
 badcn=[dict(residue_i=i+1,distance_A=x) for i,x in enumerate(cn) if x<1.2 or x>1.5]
 clashes=int(np.triu(margin< -1e-8,1).sum())
 reports.append(dict(case=case,path=str(path),sha256=hashlib.sha256(path.read_bytes()).hexdigest(),motif_backbone=metric,CN_bad_count=len(badcn),CN_bad_pairs=badcn,CN_range_A=[min(cn),max(cn)],omega_trans_deviation_gt30_count=sum(x>30 for x in omega),backbone_nonbonded_bad_pairs=clashes,backbone_nonbonded_min_margin_A=float(margin.min()),backbone_material_clear_counts={m:len(v) for m,v in hits.items()},backbone_material_clear_site_indices=hits,selected_sites_backbone_clear=sum(i in hits[j['material']] for i in j['selected_site_indices']),selected_sites_n=len(j['selected_site_indices']),alignment='Single rigid least-squares alignment of all nine fixed backbone residues to exact input frame; no site-specific adjustment',catalytic_sidechain_geometry='NOT_EVALUATED_BACKBONE_ONLY',donor_H_geometry='NOT_EVALUATED_BACKBONE_ONLY',full_atom_material_clear='NOT_EVALUATED_BACKBONE_ONLY'))
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
dest=base/'first19_RFD3_v2'/('backbone_audit_'+stamp+'.json')
report=dict(status='SNAPSHOT_NOT_FINAL_CATALYTIC_PASS',expected=57,audited=len(reports),pending=57-len(reports),reports=reports,limitations=['No sidechains in RFD output','Backbone material clearance is necessary not sufficient','No folding or enzyme activity inference','Final PDB may precede process exit; controller events establish technical completion'])
with dest.open('x') as f:json.dump(report,f,indent=2)
with (r/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=stamp,event='top19_backbone_audit',output=str(dest),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),audited=len(reports)))+'\n')
print(json.dumps(dict(output=str(dest),audited=len(reports),rows=[{k:v for k,v in x.items() if k in ['case','motif_backbone','CN_bad_count','omega_trans_deviation_gt30_count','backbone_nonbonded_bad_pairs','backbone_material_clear_counts']} for x in reports])))
