#!/usr/bin/env python
"""Polymer-only exposure / external-water envelope atlas. No catalyst data used."""
import argparse, hashlib, json, time, sys, itertools
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.ndimage import label, gaussian_filter, map_coordinates

R=Path(__file__).resolve().parent
RAD={'C':1.70,'N':1.55,'O':1.52,'H':1.20,'S':1.80}
PROBE=1.4
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def emit(event,**kw):
 print(json.dumps(dict(event=event,**kw)),flush=True)
def fibonacci(n):
 i=np.arange(n,dtype=float);z=1-2*(i+.5)/n;phi=i*np.pi*(3-np.sqrt(5))
 return np.column_stack([np.sqrt(1-z*z)*np.cos(phi),np.sqrt(1-z*z)*np.sin(phi),z])
class Env:
 def __init__(self,xyz,elements,box):
  self.xyz=np.array(xyz,float);self.elements=np.array(elements);self.box=np.array(box,float)
  self.z0=float(self.xyz[:,2].min()-30)
  self.trees={}
  w=self.wrap(self.xyz)
  for el in np.unique(elements):
   if el not in RAD:raise ValueError('Unsupported element '+el)
   self.trees[el]=cKDTree(w[self.elements==el],boxsize=[box[0],box[1],100000.])
 def wrap(self,p):
  q=np.array(p,float,copy=True);q[...,:2]%=self.box[:2];q[...,2]-=self.z0
  return q
 def clearance(self,p):
  w=self.wrap(p);out=np.full(len(w),np.inf)
  for el,t in self.trees.items():
   d=t.query(w,workers=2)[0]-RAD[el];np.minimum(out,d,out=out)
  return out
 def nearest(self,p):
  w=self.wrap(p);d=np.full(len(w),np.inf);e=np.full(len(w),'',dtype='U1')
  for el,t in self.trees.items():
   z=t.query(w,workers=2)[0]-RAD[el];b=z<d;d[b]=z[b];e[b]=el
  return e
def external_components(free):
 lab,n=label(free);parent=np.arange(n+1)
 def root(x):
  while parent[x]!=x:
   parent[x]=parent[parent[x]];x=parent[x]
  return x
 for a,b in [(lab[0,:,:],lab[-1,:,:]),(lab[:,0,:],lab[:,-1,:])]:
  pairs=np.unique(np.column_stack([a.ravel(),b.ravel()]),axis=0)
  for x,y in pairs:
   if x and y:
    rx,ry=root(x),root(y)
    if rx!=ry:parent[ry]=rx
 roots=np.array([root(i) for i in range(n+1)])
 seed=np.unique(np.concatenate([lab[:,:,0].ravel(),lab[:,:,-1].ravel()]))
 sr=set(roots[seed[seed>0]].tolist())
 good=np.array([v in sr for v in roots]);good[0]=False
 return good[lab]
def make_frame(normal,co,leaving):
 z=np.array(normal,float);z/=np.linalg.norm(z)
 x=co-z*np.dot(co,z)
 if np.linalg.norm(x)<1e-5:x=leaving-z*np.dot(leaving,z)
 if np.linalg.norm(x)<1e-5:
  v=np.eye(3)[np.argmin(abs(z))];x=v-z*np.dot(v,z)
 x/=np.linalg.norm(x);y=np.cross(z,x)
 return np.column_stack([x,y,z])
def exposure_bin(f):
 if f==0:return 'sampled_zero'
 if f<=.05:return 'low_0_5pct'
 if f<=.15:return 'mid_5_15pct'
 return 'high_gt15pct'
def support_bin(f):
 return '70_80' if f<.8 else ('80_90' if f<.9 else '90_100')
def build_grid(env,masses,h):
 nxy=np.ceil(env.box[:2]/h).astype(int);zlo=np.floor(env.xyz[:,2].min())-7
 nz=int(np.ceil((env.xyz[:,2].max()+7-zlo)/h))
 step=np.array([*(env.box[:2]/nxy),h]);shape=tuple(map(int,[*nxy,nz]));origin=np.array([0,0,zlo])+.5*step
 clr=np.empty(np.prod(shape),np.float32)
 for start in range(0,len(clr),100000):
  ii=np.column_stack(np.unravel_index(np.arange(start,min(start+100000,len(clr))),shape))
  clr[start:start+len(ii)]=env.clearance(origin+ii*step)
 clr=clr.reshape(shape);external=external_components(clr>=PROBE)
 xyz=env.xyz.copy();xyz[:,:2]%=env.box[:2]
 edges=[np.arange(shape[k]+1)*step[k]+origin[k]-.5*step[k] for k in range(3)]
 hist=np.histogramdd(xyz,bins=edges,weights=masses)[0].astype(np.float32)/np.prod(step)
 sigma=3/step
 grad=[]
 for k in range(3):
  order=[0,0,0];order[k]=1
  grad.append(gaussian_filter(hist,sigma,order=order,mode=['wrap','wrap','nearest'])/step[k])
 return dict(origin=origin,step=step,shape=shape,clearance=clr,external=external,gradient=grad)
def gridcoords(p,g,env):
 q=np.array(p,copy=True);q[:,:2]%=env.box[:2]
 return ((q-g['origin'])/g['step']).T
def interp(a,p,g,env,order=1):
 # XY wrapping and explicit periodic halo; Z does not wrap.
 x=np.pad(a,((1,1),(1,1),(0,0)),mode='wrap')
 c=gridcoords(p,g,env);c[:2]+=1
 return map_coordinates(x,c,order=order,mode='nearest')
def connected_shell(p,env,g):
 if len(p)==0:return np.zeros(0,bool)
 gc=gridcoords(p,g,env).T;lo=np.floor(gc).astype(int);ok=np.zeros(len(p),bool)
 for delta in itertools.product([0,1],repeat=3):
  ix=lo+delta;ix[:,:2]%=np.array(g['shape'][:2]);valid=(ix[:,2]>=0)&(ix[:,2]<g['shape'][2])&~ok
  ids=np.where(valid)[0];j=ix[ids]
  ids=ids[g['external'][j[:,0],j[:,1],j[:,2]]]
  if len(ids)==0:continue
  q=g['origin']+ix[ids]*g['step'];v=q-p[ids];v[:,:2]-=env.box[:2]*np.round(v[:,:2]/env.box[:2])
  # Eight checks along the short segment to an exterior-connected grid node.
  t=np.linspace(.125,1,8)
  seg=p[ids,None,:]+v[:,None,:]*t[None,:,None]
  clear=env.clearance(seg.reshape(-1,3)).reshape(-1,8)
  ok[ids[np.all(clear>=PROBE-1e-7,axis=1)]]=True
 return ok
def exposure(center,env,g,n=4096):
 u=fibonacci(n);p=center+(RAD['C']+PROBE)*u
 local=env.clearance(p)>=PROBE-1e-7;ext=np.zeros(n,bool)
 ext[local]=connected_shell(p[local],env,g)
 return u,local,ext
def write_dx(path,a):
 with open(path,'x') as f:
  nx,ny,nz=a.shape
  f.write(f'object 1 class gridpositions counts {nx} {ny} {nz}\norigin -15 -15 -15\ndelta 1 0 0\ndelta 0 1 0\ndelta 0 0 1\nobject 2 class gridconnections counts {nx} {ny} {nz}\nobject 3 class array type double rank 0 items {a.size} data follows\n')
  v=a.ravel()
  for i in range(0,len(v),3):f.write(' '.join(f'{x:.6g}' for x in v[i:i+3])+'\n')
  f.write('attribute "dep" string "positions"\nobject "field" class field\ncomponent "positions" value 1\ncomponent "connections" value 2\ncomponent "data" value 3\n')
def fit_groups(rows,fields,axis,out):
 d=out/'fits';d.mkdir();summ=[];rng=np.random.default_rng(20260908)
 coords=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),-1);sphere=np.linalg.norm(coords,axis=-1)<=15
 for exposure_level in ['sampled_zero','low_0_5pct','mid_5_15pct','high_gt15pct','CONNECTIVITY_UNRESOLVED']:
  for side in ['top','bottom']:
   for support in ['all','70_80','80_90','90_100']:
    ids=[i for i,r in enumerate(rows) if r['exposure_group']==exposure_level and r['side']==side and (support=='all' or r['support_group']==support) and r['normal_status']=='PASS']
    key=f'{exposure_level}__{side}__support_{support}'
    if not ids:
     summ.append(dict(group=key,n=0,status='EMPTY'));continue
    arr=fields[ids].astype(np.float32)/255;mean=arr.mean(0);sd=arr.std(0)
    band=(mean>.05)&(mean<.95)&sphere
    dist=np.mean(abs(arr[:,sphere]-mean[sphere]),axis=1)
    medoid=ids[int(np.argmin(dist))]
    keys=sorted(set(rows[i]['chain'] for i in ids));ci=None
    if len(keys)>=5:
     sums=np.stack([arr[[rows[i]['chain']==k for i in ids]].sum(0) for k in keys])
     counts=np.array([sum(rows[i]['chain']==k for i in ids) for k in keys])
     # Chain bootstrap on coarser fixed 3A diagnostic lattice, not a full-field CI.
     mask=np.zeros(mean.shape,bool);mask[::3,::3,::3]=True;mask&=band
     v=sums[:,mask];bs=[]
     for _ in range(100):
      pick=rng.integers(0,len(keys),len(keys));bs.append(v[pick].sum(0)/counts[pick].sum())
     if mask.any():
      q=np.quantile(bs,[.05,.95],axis=0);ci=float(np.mean(q[1]-q[0]))
    summary=dict(group=key,n=len(ids),chains=len(keys),status='DESCRIPTIVE_ONLY' if len(keys)<5 else 'DESCRIPTIVE_WITH_CHAIN_BOOTSTRAP',representative_site=rows[medoid]['site_id'],mean_absolute_field_deviation=float(dist.mean()),interface_band_voxels=int(band.sum()),bootstrap90_mean_width=ci,bootstrap_scope='3A diagnostic lattice; original-chain blocks; single snapshot; spatial cross-chain dependence not removed')
    np.savez_compressed(d/(key+'.npz'),occupancy_probability=mean,occupancy_sd=sd,site_indices=np.array(ids),representative_field=fields[medoid],sphere_mask=sphere)
    if support=='all':
     write_dx(d/(key+'.dx'),mean)
     write_dx(d/(key+'__representative.dx'),fields[medoid].astype(float)/255)
    summ.append(summary)
 return summ
def material(mat,sites,out,h=1,n=4096,pilot=0):
 src=R/'inputs_v2'/f'{mat}_parent_atoms.npz';z=np.load(src)
 env=Env(z['xyz_A'],z['elements'],z['box_A'])
 emit('grid_start',material=mat,h=h)
 g=build_grid(env,z['masses_Da'],h)
 emit('grid_ready',material=mat,shape=g['shape'],external_voxels=int(g['external'].sum()))
 rows=[];axis=np.arange(-15,16,dtype=float);localgrid=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),-1).reshape(-1,3)
 sphere=np.linalg.norm(localgrid,axis=1)<=15
 selected=[s for s in sites if s['surface']==mat+'_original']
 if pilot:selected=selected[:pilot]
 fields=np.empty((len(selected),31,31,31),np.uint8);checks=[]
 for j,s in enumerate(selected):
  center=np.array(s['center']);cf=np.array(s['frame']);amide=np.array(s['amide']);co=cf@amide[1];leave=cf@amide[2]
  grad=np.array([interp(a,center[None],g,env)[0] for a in g['gradient']]);gn=np.linalg.norm(grad)
  normal=-grad/gn if gn>1e-10 else np.array([0.,0.,1. if s['side']=='top' else -1.])
  outward=normal[2]*(1 if s['side']=='top' else -1)
  status='PASS' if gn>1e-10 and outward>0 else 'NORMAL_AMBIGUOUS'
  frame=make_frame(normal,co,leave)
  u,loc,ext=exposure(center,env,g,n)
  lf=float(loc.mean());ef=float(ext.mean());area=4*np.pi*3.1**2
  group=exposure_bin(ef)
  if ef==0 and lf>0:group='CONNECTIVITY_UNRESOLVED'
  elif exposure_bin(lf)!=exposure_bin(ef):group='CONNECTIVITY_UNRESOLVED'
  world=center+localgrid@frame.T
  # Interpolated global external-water occupancy, not an SES or atomistic vdW surface.
  accessible=interp(g['external'].astype(np.float32),world,g,env)
  field=1-accessible;fields[j]=np.rint(np.clip(field,0,1)*255).astype(np.uint8).reshape(31,31,31)
  clear=env.clearance(world)
  surf=(abs(clear-PROBE)<=.75)&sphere&(accessible>.1)
  chem=env.nearest(world[surf]) if surf.any() else np.array([],dtype='U1')
  band=localgrid[surf]
  opening=u[ext]@frame if ext.any() else np.empty((0,3))
  row=dict(material=mat,site_id=s['site_id'],chain=s['chain'],side=s['side'],support=s['support'],support_group=support_bin(s['support']),local_shell_fraction=lf,external_shell_fraction=ef,unresolved_fraction=lf-ef,local_SASA_A2=lf*area,external_SASA_A2=ef*area,exposure_group=group,samples=n,normal_status=status,density_gradient_norm=float(gn),normal_outward_cosine=float(outward),surface_frame=frame.tolist(),chemical_to_surface= (frame.T@cf).tolist(),center=center.tolist(),co_normal_cosine=float(np.dot(co/np.linalg.norm(co),normal)),carbon_probe_clearance_A=float(env.clearance(center[None])[0]),external_water_fraction_R15=float(accessible[sphere].mean()),local_boundary_z_quantiles_A=np.quantile(band[:,2],[.1,.5,.9]).tolist() if len(band) else None,external_shell_direction_mean_surface=opening.mean(0).tolist() if len(opening) else None,near_boundary_nearest_element_counts={el:int((chem==el).sum()) for el in ['C','N','O','H']},chemical_annotation_scope='voxel counts near probe-accessible boundary; not area fractions or hydrophobic potential')
  assert 0<=ef<=lf<=1 and np.allclose(frame.T@frame,np.eye(3),atol=1e-8)
  rows.append(row)
  if j<3:
   _,ll,ee=exposure(center,env,g,n*2)
   checks.append(dict(site_id=s['site_id'],shell_samples_refined=n*2,local_fraction_refined=float(ll.mean()),external_fraction_refined=float(ee.mean()),local_fraction_abs_change=abs(float(ll.mean())-lf),external_fraction_abs_change=abs(float(ee.mean())-ef)))
  if j%50==0:emit('sites',material=mat,done=j+1,total=len(selected),last_group=group)
 d=out/mat;d.mkdir()
 np.savez_compressed(d/'site_fields.npz',fields=fields,axis_A=axis)
 (d/'sites.json').write_text(json.dumps(rows,indent=2))
 fs=fit_groups(rows,fields,axis,d)
 summary=dict(material=mat,n=len(rows),source_sha256=sha(src),grid_spacing_A=g['step'].tolist(),grid_shape=list(g['shape']),grid_resolution_status='NOT_YET_VALIDATED_BY_REFINEMENT',sampling_checks=checks,exposure_counts={k:sum(r['exposure_group']==k for r in rows) for k in ['sampled_zero','low_0_5pct','mid_5_15pct','high_gt15pct','CONNECTIVITY_UNRESOLVED']},normal_ambiguous=sum(r['normal_status']!='PASS' for r in rows),fits=fs)
 (d/'summary.json').write_text(json.dumps(summary,indent=2))
 emit('material_complete',material=mat,n=len(rows),counts=summary['exposure_counts'])
 return summary
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--pilot',type=int,default=0);ap.add_argument('--grid',type=float,default=1);ap.add_argument('--samples',type=int,default=4096);args=ap.parse_args()
 out=R/args.output
 assert out.parent==R and not out.exists()
 out.mkdir();start=time.time()
 log=R/'RUN_LOG.jsonl'
 with log.open('a') as f:f.write(json.dumps(dict(time=time.strftime('%Y-%m-%dT%H:%M:%S%z'),event='polymer_exposure_start',command=sys.argv,source_sha256=sha(__file__),output=str(out)))+'\n')
 sites=json.load(open(R/'panel_surface_scan_v1/sites.json'))
 result=[material(mat,sites,out,args.grid,args.samples,args.pilot) for mat in ['PET','PA6','PA66']]
 doc=dict(status='TECHNICAL_COMPLETE_GRID_REFINEMENT_NOT_EVALUATED',seconds=time.time()-start,source_sha256=sha(__file__),site_source_sha256=sha(R/'panel_surface_scan_v1/sites.json'),materials=result,definitions=dict(radius_A=15,probe_radius_A=1.4,vdW_radii_A=RAD,sphere_denominator='4*pi*(C_radius+water_radius)^2; not isolated ester reference',surface='external water-center exclusion occupancy; not solvent-excluded molecular surface, not enzyme shape',normal='negative gradient of Gaussian3A mass density; carbonyl projected into tangent plane',bins='0;(0,.05];(.05,.15];(.15,1], conventions only',connectivity='6-neighbor grid flood from top/bottom, XY periodic; shell link 8 segment samples; resolution dependent',input_scope='all existing original dense-surface sites support>=.70; no catalytic filtering; density selector itself remains a cohort restriction',limitations=['single snapshot','no continuous access guarantee below grid resolution','sampled zero is below sampling resolution, not proof of zero area','no shape subclustering validated','no true electrostatic or hydrophobic potential']))
 (out/'summary.json').write_text(json.dumps(doc,indent=2))
 with (R/'RUNBOOK.md').open('a') as f:f.write('\nPolymer-only exposure atlas: '+str(out)+'. Reproduce with CPU Python polymer_exposure_v1.py --output NEW_UNIQUE_NAME --grid 1 --samples 4096. Definitions and limits in summary.json. No catalyst filters. DX maps describe water-center exclusion envelope, not atomic SES. Unresolved connectivity bins kept separate.\n')
 with log.open('a') as f:f.write(json.dumps(dict(time=time.strftime('%Y-%m-%dT%H:%M:%S%z'),event='polymer_exposure_complete',output=str(out),seconds=doc['seconds'],exit_status=0))+'\n')
 emit('complete',output=str(out),seconds=doc['seconds'])
if __name__=='__main__':main()
