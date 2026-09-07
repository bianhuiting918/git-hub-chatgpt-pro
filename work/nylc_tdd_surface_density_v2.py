#!/usr/bin/env python3
"""Finite trimmed nylon density boundaries, using parent slab reference density."""
from pathlib import Path
import json,datetime,hashlib,sys
import numpy as np
from scipy.ndimage import gaussian_filter,label,map_coordinates
from scipy.interpolate import RegularGridInterpolator
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'surface_previews/axis_density_sixface_all_external_trim_v6'))
from density_cut_variants import parse_pdb
from terminal_density_trim_3d import center_whole_chains_in_primary_cell
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def run(mat,dp,out):
 inp=ROOT/'inputs'/f'{mat}_v10_ge90_uncapped.pdb'
 audit=json.loads((ROOT/'input_audit_v2'/f'{mat}_intact_amides.json').read_text())
 assert sha(inp)==audit['input_sha256']
 parent=ROOT.parent/f'{mat.lower()}_dp{dp}_400chain/water_slab_v1/dry_export/{mat}_DP{dp}_400chain_water_equilibrated_dry.pdb'
 assert sha(parent)==audit['parent_sha256']
 _,box,pa,_=parse_pdb(parent); box=np.array(box)
 center_whole_chains_in_primary_cell(pa,box)
 pxyz=np.array([a['xyz'] for a in pa]); pxyz[:,:2]%=box[:2]
 masses={'C':12.011,'O':15.999,'N':14.007,'H':1.008}
 pmass=np.array([masses[a['element']] for a in pa])
 lines=[l for l in inp.read_text().splitlines() if l.startswith(('ATOM  ','HETATM'))]
 xyz=np.array([[float(l[i:i+8]) for i in (30,38,46)] for l in lines])
 mass=np.array([masses[l[76:78].strip()] for l in lines])
 n=np.ceil(box).astype(int); step=box/n
 edges=[np.linspace(0,box[a],n[a]+1) for a in range(3)]
 centers=[(e[:-1]+e[1:])/2 for e in edges]
 def density(x,m,mode):
  h=np.histogramdd(x,bins=edges,weights=m)[0]
  assert abs(h.sum()-m.sum())<1e-5,'Atoms outside histogram'
  return gaussian_filter(h/np.prod(step),3/step,mode=mode)
 parentrho=density(pxyz,pmass,('wrap','wrap','constant'))
 reference=float(np.median(parentrho[:,:,abs(centers[2]-np.median(pxyz[:,2]))<10]))
 assert reference>0
 threshold=reference*.5
 rho=density(xyz,mass,'constant')
 lab,num=label(rho>=threshold); counts=np.bincount(lab.ravel()); counts[0]=0
 assert counts.max()>0
 mask=lab==counts.argmax(); valid=mask.any(2)
 side_data={}; results={}
 for side,sign in [('top',1),('bottom',-1)]:
  ix=np.where(valid,np.max(np.where(mask,np.arange(n[2]),-1),axis=2) if sign==1 else np.min(np.where(mask,np.arange(n[2]),n[2]),axis=2),0)
  adj=np.clip(ix+sign,0,n[2]-1)
  v0=np.take_along_axis(rho,ix[:,:,None],2)[:,:,0]
  v1=np.take_along_axis(rho,adj[:,:,None],2)[:,:,0]
  fraction=np.clip((threshold-v0)/np.where(abs(v1-v0)>1e-12,v1-v0,1),0,1)
  height=centers[2][ix]+fraction*(centers[2][adj]-centers[2][ix]); height[~valid]=np.nan
  xx,yy=np.meshgrid(centers[0],centers[1],indexing='ij')
  vertices=np.stack((xx,yy,height),axis=-1).reshape(-1,3)
  grid=np.arange(n[0]*n[1]).reshape(n[:2])
  good=valid[:-1,:-1]&valid[1:,:-1]&valid[1:,1:]&valid[:-1,1:]
  a,b,c,d=[v[good] for v in (grid[:-1,:-1],grid[1:,:-1],grid[1:,1:],grid[:-1,1:])]
  faces=np.concatenate((np.stack((a,b,c),1),np.stack((a,c,d),1)))
  if sign==-1:faces=faces[:,::-1]
  np.savez_compressed(out/f'{mat}_{side}_mesh.npz',vertices_A=vertices,faces=faces,height_A=height,valid=valid,step_A=step,threshold=threshold)
  interp=RegularGridInterpolator(centers[:2],height,bounds_error=False,fill_value=np.nan)
  side_data[side]=(interp,sign)
  results[side]={'faces':len(faces),'valid_xy_columns':int(valid.sum())}
 rows=[]
 disk=np.array([(x,y) for x in np.linspace(-10,10,11) for y in np.linspace(-10,10,11) if x*x+y*y<=100])
 assert len(disk)==81
 for source in audit['sites']:
  q=np.array(source['carbonyl_xyz_A']); row=dict(source)
  offsets={}
  for side,(interp,sign) in side_data.items():
   h=float(interp(q[:2][None])[0]); offsets[side]=float(sign*(q[2]-h)) if np.isfinite(h) else None
  finite=[k for k,v in offsets.items() if v is not None]
  side=min(finite,key=lambda k:abs(offsets[k])) if finite else None
  row['nearest_side']=side; row['signed_surface_offsets_A']=offsets
  row['within_3A_density_boundary']=side is not None and abs(offsets[side])<=3
  if side:
   sign=side_data[side][1]; points=np.column_stack((q[0]+disk[:,0],q[1]+disk[:,1],np.full(81,q[2]-sign*5)))
   values=map_coordinates(rho,((points/step)-.5).T,order=1,mode='constant',cval=0)
   row['inward_support']=float(np.mean(values>=threshold))
  else:row['inward_support']=None
  rows.append(row)
 summary=dict(material=mat,reference_density_Da_A3=reference,threshold_Da_A3=threshold,sigma_A=3,reference='untrimmed parent central slab median, XY periodic',object_boundary='finite v10, no periodic copies',intact_sites=len(rows),near_boundary=sum(r['within_3A_density_boundary'] for r in rows),cohorts={str(t):sum(r['within_3A_density_boundary'] and r['inward_support']>=t for r in rows) for t in (.9,.8,.7)},meshes=results,status='COARSE_DENSITY_BOUNDARY_ONLY_ATOMIC_EXPOSURE_NOT_EVALUATED')
 (out/f'{mat}_sites.json').write_text(json.dumps(rows,indent=2))
 return summary
if __name__=='__main__':
 out=ROOT/'surface_density_v2';out.mkdir(exist_ok=False)
 status='FAILED';rc=1
 try:
  summary=[run('PA6',16,out),run('PA66',8,out)]
  (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
  status='PASS_COARSE_BOUNDARY';rc=0
 finally:
  with (ROOT/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__).resolve()),script_sha256=sha(Path(__file__)),command='CPU Python surface_density_v2.py',inputs='verified v10 PDBs and original parent dry PDBs',outputs=str(out),parameters=dict(sigma_A=3,threshold_fraction=.5,boundary_proximity_A=3,support_depth_A=5,support_radius_A=10,support_points=81),exit_code=rc,result=status))+'\n')
