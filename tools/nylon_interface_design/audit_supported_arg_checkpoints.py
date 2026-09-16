"""Read-only structural audit of matched supported-motif RFD outputs.
Writes only a timestamped report and append-only run record. Does not repair,
submit, predict or classify backbone-only output as an active enzyme.
"""
from pathlib import Path
import json,datetime,hashlib
import numpy as np
from scipy.spatial.distance import cdist
r=Path.cwd().resolve();d=r/'route56_N10_20260915_v1';g=r.parent/'n9_core_donor_design_20260910_gpu_outputs_v1'
case='PA6_SC4559268_coverage10_variant1'
j=json.loads((d/'route6_supported_ARG_SC4559268_v1'/case/'manifest.json').read_text())
bb=['N','CA','C','O']
def read(path):
 res={}
 for line in path.read_text().splitlines():
  if line.startswith('ATOM  ') and line[12:16].strip() in bb:
   key=int(line[22:26]);name=line[12:16].strip()
   assert name not in res.setdefault(key,{})
   res[key][name]=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])
 return res
inp=read(Path(j['pdb']))
mapping=list(zip(range(1,10),[1,40,41,42,61,62,63,64,65]))
def fitmetric(a,b):
 u,s,v=np.linalg.svd((a-a.mean(0)).T@(b-b.mean(0)));rot=u@np.diag([1,1,np.linalg.det(u@v)])@v
 delta=(a-a.mean(0))@rot+b.mean(0)-b
 return {'RMSD_A':float(np.sqrt(np.mean(np.sum(delta*delta,axis=1)))),'max_A':float(np.linalg.norm(delta,axis=1).max())}
def tors(a,b,c,d):
 u=b-a;w=c-b;t=d-c;n=np.cross(u,w);m=np.cross(w,t)
 return float(np.degrees(np.arctan2(np.dot(np.cross(n,m),w/np.linalg.norm(w)),np.dot(n,m))))
reports=[]
for label,path in [('ActiveSite',g/'route6_supported_ARG_SC4559268_v1'/case/'design_0.pdb'),('Base',g/'route6_supported_ARG_SC4559268_base_v1/design_0.pdb')]:
 assert path.is_file(),f'NOT_READY: {path}'
 out=read(path);assert sorted(out)==list(range(1,127)) and all(set(v)==set(bb) for v in out.values())
 cn=[float(np.linalg.norm(out[i]['C']-out[i+1]['N'])) for i in range(1,126)]
 omega=[tors(out[i]['CA'],out[i]['C'],out[i+1]['N'],out[i+1]['CA']) for i in range(1,126)]
 metrics={}
 for name,pairs in [('all_fixed_backbone',mapping),('common_catalytic_support_backbone',[p for p in mapping if p[1] in [1,40,41,42,63]])]:
  a=np.array([inp[i][b] for i,k in pairs for b in bb]);v=np.array([out[k][b] for i,k in pairs for b in bb]);metrics[name]=fitmetric(a,v)
 xyz=np.array([out[i][b] for i in range(1,127) for b in bb]);rr=np.tile([1.55,1.7,1.7,1.52],126);n=len(xyz)
 adjacency=np.zeros((n,n),bool)
 for i in range(126):
  bonds=[(4*i,4*i+1),(4*i+1,4*i+2),(4*i+2,4*i+3)]+([(4*i+2,4*i+4)] if i<125 else [])
  for a,b in bonds:adjacency[a,b]=adjacency[b,a]=True
 reach=np.eye(n,dtype=bool);front=reach.copy()
 for _ in range(3):front=(front.astype(np.int16)@adjacency.astype(np.int16))>0;reach|=front
 margin=cdist(xyz,xyz)-rr[:,None]-rr[None,:]+.4;margin[reach]=np.inf
 reports.append({'checkpoint':label,'path':str(path),'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'residues':126,'atoms':504,'motif_backbone':metrics,'CN_bad_count':sum(v<1.2 or v>1.5 for v in cn),'CN_range_A':[min(cn),max(cn)],'CN_bad_pairs':[{'residue_i':i+1,'distance_A':v} for i,v in enumerate(cn) if v<1.2 or v>1.5],'support_boundaries_A':[cn[59],cn[64]],'omega_trans_deviation_gt30_count':sum(180-abs(v)>30 for v in omega),'backbone_nonbonded_bad_pairs':int(np.triu(margin< -1e-8,1).sum()),'backbone_nonbonded_min_margin_A':float(margin.min()),'all_atom_catalytic_geometry':'NOT_EVALUATED_BACKBONE_ONLY','actual_material_replay':'NOT_EVALUATED_RAW_BACKBONE'})
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
dest=d/'route6_supported_ARG_SC4559268_v1'/('matched_checkpoint_audit_'+stamp+'.json')
report={'status':'AUDIT_COMPLETE_NOT_SCIENTIFIC_PASS','source_motif_sha256':hashlib.sha256(Path(j['pdb']).read_bytes()).hexdigest(),'reports':reports}
with dest.open('x') as h:json.dump(report,h,indent=2)
with (r/'RUN_LOG.jsonl').open('a') as h:h.write(json.dumps({'time':stamp,'event':'supported_ARG_matched_checkpoint_audit','script':str(Path(__file__).resolve()),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'output':str(dest),'exit_code':0})+'\n')
print(json.dumps({'output':str(dest),'reports':[{k:v for k,v in row.items() if k!='CN_bad_pairs'} for row in reports]},indent=2))
