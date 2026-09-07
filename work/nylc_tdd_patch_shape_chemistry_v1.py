#!/usr/bin/env python3
"""Connected 15 A density-surface patches and nearest-VDW chemistry proxy."""
from pathlib import Path
import numpy as np,json,datetime,hashlib
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(mat,out):
 inp=ROOT/'inputs'/f'{mat}_v10_ge90_uncapped.pdb'
 lines=inp.read_text().splitlines()
 atoms={l[6:11]:l for l in lines if l.startswith(('ATOM  ','HETATM'))}
 carbonyl=set()
 for l in lines:
  if not l.startswith('CONECT'):continue
  fields=[l[i:i+5] for i in range(6,len(l),5) if l[i:i+5] in atoms]
  for b in fields[1:]:
   a=fields[0];ea=atoms[a][76:78].strip();eb=atoms[b][76:78].strip()
   if ea=='C' and eb=='O':carbonyl.add(a)
   if eb=='C' and ea=='O':carbonyl.add(b)
 classes={'polar_O':[],'polar_N':[],'carbonyl_C':[],'nonpolar_C':[]}
 for a,l in atoms.items():
  e=l[76:78].strip()
  if e=='H':continue
  key='polar_'+e if e in ('O','N') else ('carbonyl_C' if a in carbonyl else 'nonpolar_C')
  classes[key].append([float(l[k:k+8]) for k in (30,38,46)])
 radii={'polar_O':1.52,'polar_N':1.55,'carbonyl_C':1.7,'nonpolar_C':1.7}
 trees={k:cKDTree(v) for k,v in classes.items() if v}
 labels=list(trees)
 sites=json.loads((ROOT/'carbon_exposure_v1'/f'{mat}_sites.json').read_text())
 patches={};rows=[];whole={}
 for side in ('top','bottom'):
  z=np.load(ROOT/'surface_density_v2'/f'{mat}_{side}_mesh.npz')
  v=z['vertices_A'];f=z['faces'];tri=v[f];cent=tri.mean(1)
  cross=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);area=np.linalg.norm(cross,axis=1)/2
  assert np.isfinite(tri).all() and (area>0).all()
  chemistry=np.argmin(np.stack([trees[k].query(cent)[0]-radii[k] for k in labels],1),axis=1)
  # Shared-edge adjacency, not merely spatial proximity.
  edgefaces={}
  for j,face in enumerate(f):
   for a,b in [(face[0],face[1]),(face[1],face[2]),(face[2],face[0])]:edgefaces.setdefault(tuple(sorted((int(a),int(b)))),[]).append(j)
  links=[]
  for ids in edgefaces.values():
   assert len(ids)<=2
   if len(ids)==2:links.extend([(ids[0],ids[1]),(ids[1],ids[0])])
  links=np.array(links)
  adj=csr_matrix((np.ones(len(links),dtype=bool),(links[:,0],links[:,1])),shape=(len(f),len(f)))
  tree=cKDTree(cent)
  np.savez_compressed(out/f'{mat}_{side}_chemistry.npz',class_indices=chemistry,class_names=np.array(labels),face_area_A2=area)
  whole[side]=dict(area_A2=float(area.sum()),class_area_A2={k:float(area[chemistry==i].sum()) for i,k in enumerate(labels)})
  for row in sites:
   if row['nearest_side']!=side:continue
   q=np.array(row['carbonyl_xyz_A']);ids=np.array(tree.query_ball_point(q,15),dtype=int)
   assert len(ids)>0
   n,lab=connected_components(adj[ids][:,ids],directed=False)
   near=int(np.argmin(np.linalg.norm(cent[ids]-q,axis=1)));chosen=ids[lab==lab[near]]
   # Verify retained component is one connected component.
   assert connected_components(adj[chosen][:,chosen],directed=False,return_labels=False)==1
   patches[row['site_id']]=chosen
   w=area[chosen];pts=cent[chosen];mu=np.average(pts,axis=0,weights=w);delta=pts-mu
   cov=(delta*w[:,None]).T@delta/w.sum();ev,evec=np.linalg.eigh(cov);normal=evec[:,0]
   if normal[2]*(1 if side=='top' else -1)<0:normal=-normal
   heights=delta@normal
   r=dict(site_id=row['site_id'],side=side,support=row['inward_support'],parent_outward_points=row['parent_outward_points'],cut_endpoint_within_15A=row['cut_endpoint_within_15A'],triangles=len(chosen),components_before=n,area_A2=float(w.sum()),centroid_plane_rms_A=float(np.sqrt(np.average(heights**2,weights=w))),centroid_plane_height_range_A=float(np.ptp(heights)),plane_center_A=mu.tolist(),water_oriented_normal=normal.tolist(),class_area_fraction={k:float(w[chemistry[chosen]==i].sum()/w.sum()) for i,k in enumerate(labels)})
   assert abs(sum(r['class_area_fraction'].values())-1)<1e-10
   rows.append(r)
 np.savez_compressed(out/f'{mat}_patch_face_indices.npz',**patches)
 (out/f'{mat}_patches.json').write_text(json.dumps(rows,indent=2))
 return dict(material=mat,patches=len(rows),whole_surface=whole,median_patch_area_A2=float(np.median([r['area_A2'] for r in rows])),median_rms_A=float(np.median([r['centroid_plane_rms_A'] for r in rows])),scope='density boundary; centroid radius 15 A; nearest-target shared-edge component; nearest heavy-atom VDW chemistry proxy, not hydrophobic potential or atomistic SES')
if __name__=='__main__':
 out=ROOT/'patch_shape_chemistry_v1';out.mkdir(exist_ok=False);rc=1
 try:
  summary=[run(m,out) for m in ('PA6','PA66')]
  (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));rc=0
 finally:
  with (ROOT/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__).resolve()),script_sha256=sha(Path(__file__)),command='CPU Python patch_shape_chemistry_v1.py',inputs='surface_density_v2 meshes, exposure site lists, v10 PDBs',outputs=str(out),parameters=dict(radius_A=15,connectivity='shared edges',chemical_proxy='nearest heavy VDW surface'),exit_code=rc))+'\n')
