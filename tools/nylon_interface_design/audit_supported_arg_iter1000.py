"""Independent fixed-core repair audit. CPU; no coordinate modification."""
from pathlib import Path
import numpy as np,json,datetime,hashlib
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
r=Path.cwd().resolve();d=r/'route56_N10_20260915_v1';o=d/'route6_supported_ARG_SC4559268_v1';path=o/'fixed_core_repair_iter1000_v1/after.pdb'
a={};elements={}
for l in path.read_text().splitlines():
 if l.startswith('ATOM  '):
  i=int(l[22:26]);n=l[12:16].strip();a.setdefault(i,{})[n]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]);elements[i,n]=l[76:78].strip()
assert sorted(a)==list(range(1,127))
rad={'C':1.7,'N':1.55,'O':1.52,'S':1.8};keys=[(i,n) for i in a for n in a[i] if elements[i,n] in rad];idx={k:i for i,k in enumerate(keys)};q=np.array([a[i][n] for i,n in keys]);rr=np.array([rad[elements[k]] for k in keys]);N=len(q);adj=np.zeros((N,N),bool)
def edge(k,l):
 i,j=idx[k],idx[l];adj[i,j]=adj[j,i]=True
for i in a:
 for b,c in [('N','CA'),('CA','C'),('C','O')]:edge((i,b),(i,c))
 if i<126:edge((i,'C'),(i+1,'N'))
 if 'OXT' in a[i]:edge((i,'C'),(i,'OXT'))
 if 'CB' in a[i]:edge((i,'CA'),(i,'CB'))
 if i==1:
  for b in ['OG1','CG2']:edge((i,'CB'),(i,b))
 elif i in [40,42]:
  edge((i,'CB'),(i,'CG'))
  for b in ['OD1','OD2']:edge((i,'CG'),(i,b))
 elif i==63:
  for b,c in [('CB','CG'),('CG','CD'),('CD','NE'),('NE','CZ'),('CZ','NH1'),('CZ','NH2')]:edge((i,b),(i,c))
reach=np.eye(N,dtype=bool);front=reach.copy()
for _ in range(3):front=(front.astype(np.int16)@adj.astype(np.int16))>0;reach|=front
m=cdist(q,q)-rr[:,None]-rr[None,:]+.4;m[reach]=np.inf
pairs=np.argwhere(np.triu(m< -1e-8,1))
def angle(x,y,z):
 u=x-y;v=z-y
 return float(np.degrees(np.arccos(np.clip(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v),-1,1))))
def tors(x,y,z,w):
 u=y-x;v=z-y;t=w-z;n=np.cross(u,v);m=np.cross(v,t)
 return float(np.degrees(np.arctan2(np.dot(np.cross(n,m),v/np.linalg.norm(v)),np.dot(n,m))))
cn=[float(np.linalg.norm(a[i]['C']-a[i+1]['N'])) for i in range(1,126)];omega=[180-abs(tors(a[i]['CA'],a[i]['C'],a[i+1]['N'],a[i+1]['CA'])) for i in range(1,126)]
z=np.load(o/'PA6_SC4559268_coverage10_variant1/motif_exact.npz');rot=z['rotation_to_shape'];shift=z['translation_to_shape'];center=z['new_center_A'];chemical=(q+center-shift)@rot.T
sites=json.loads((d/'N10_fixed_minimal_motif_material_v1/sites.json').read_text());p=d/'route6_N10_single_residue_pilot_v1';panel=json.loads((p/'panel.json').read_text());k=next(i for i,v in enumerate(panel) if v['id']=='SC4559268');old=np.load(p/'coverage.npz')['coverage'][k];reports=[]
for mat in ['PA6','PA66']:
 e=np.load(r.parent/'interface_surface_robustness_20260908_v1/inputs_v2'/(mat+'_parent_atoms.npz'));xyz=e['xyz_A'].copy();box=e['box_A'];xyz[:,:2]%=box[:2];trees={el:cKDTree(xyz[e['elements']==el],boxsize=[*box[:2],0]) for el in np.unique(e['elements']) if el!='H'};rows=[]
 for sj,site in enumerate(sites):
  if site['material']!=mat:continue
  w=chemical@np.array(site['frame']).T+site['center'];w[:,:2]%=box[:2];mar=np.full(N,np.inf)
  for el,tr in trees.items():mar=np.minimum(mar,tr.query(w)[0]-rr-rad[str(el)]+.4)
  oxygen=np.array(site['amide'][1])@rot+shift-center;donors=[]
  for name,hs in [('NH1',['1HH1','2HH1']),('NH2',['1HH2','2HH2'])]:
   n=a[63][name];angles=[angle(n,a[63][h],oxygen) for h in hs];donors.append({'name':name,'N_O_A':float(np.linalg.norm(n-oxygen)),'best_N_H_O_deg':max(angles)})
  rows.append({'site_index':sj,'site_id':site['site_id'],'old_eligible':bool(old[sj]),'bin':site['bin'],'clashing_atoms':int(sum(mar< -1e-8)),'minimum_margin_A':float(mar.min()),'donors':donors})
 clear=[v for v in rows if v['clashing_atoms']==0];eligible=[v for v in clear if v['old_eligible']]
 reports.append({'material':mat,'all_sites':len(rows),'all_material_clear':len(clear),'old_eligible_denominator':sum(v['old_eligible'] for v in rows),'old_eligible_material_clear':len(eligible),'old_eligible_material_and_dual_NH_geometry':sum(all(2.7<=dd['N_O_A']<=3.2 and dd['best_N_H_O_deg']>=140 for dd in v['donors']) for v in eligible),'clear_site_details':clear})
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');dest=o/'fixed_core_repair_iter1000_v1'/('independent_audit_'+stamp+'.json')
report={'status':'INDEPENDENT_GEOMETRY_AUDIT_NOT_FINAL_SEQUENCE_PASS','source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'heavy_atoms':N,'CN_bad':sum(v<1.2 or v>1.5 for v in cn),'CN_range_A':[min(cn),max(cn)],'omega_gt30_count':sum(v>30 for v in omega),'omega_max_deviation_deg':max(omega),'self_bad_pairs':len(pairs),'worst_self_margin_A':float(m.min()),'self_pair_details':[{'atoms':[keys[i],keys[j]],'margin_A':float(m[i,j])} for i,j in pairs],'reports':reports,'limitations':['122 Gly placeholder residues; no designed sequence or independent prediction','Material-clear sites outside old eligible set are not automatically catalytic passes','Donor check is the stated N-O and N-H-O geometry screen, not an energy calculation','Complete catalytic hydrogen and attack-geometry audit remains required']}
with dest.open('x') as h:json.dump(report,h,indent=2)
with (r/'RUN_LOG.jsonl').open('a') as h:h.write(json.dumps({'time':stamp,'script':str(Path(__file__).resolve()),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'output':str(dest),'exit_code':0})+'\n')
print(json.dumps({'output':str(dest),**{k:v for k,v in report.items() if k not in ['reports','self_pair_details']},'reports':[{k:v for k,v in row.items() if k!='clear_site_details'} for row in reports]},indent=2))
