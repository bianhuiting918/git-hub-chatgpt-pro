#!/usr/bin/env python3
"""Conditional spatially grouped surface fit/bootstraps; no fabricated protein surface."""
from pathlib import Path
import json,hashlib,datetime,concurrent.futures
import numpy as np
from scipy.spatial.distance import cdist
R=Path(__file__).resolve().parent;P=R/'panel_surface_scan_v1';O=R/'surface_fit_robustness_v1'
def groups(ids):
 parent=list(range(len(ids)))
 def root(i):
  while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
  return i
 seen={}
 for j,idx in enumerate(ids):
  s=SITES[idx];c=np.array(s['center']);box=BOX[s['surface']]
  if not s['surface'].endswith('trimmed'):c[:2]%=box[:2]
  for key in [('chain',s['chain']),('tile',s['side'],int(c[0]//30),int(c[1]//30))]:
   if key in seen:parent[root(j)]=root(seen[key])
   else:seen[key]=j
 return np.array([root(j) for j in range(len(ids))])
def one(task):
 name,surface,threshold,ids,kind=task
 n=len(ids);row=dict(fit_id=name,surface=surface,min_support=threshold,kind=kind,sites=n,grid_points=len(GRID))
 if not n:return dict(**row,status='NOT_EVALUATED_NO_COMPATIBLE_SITES')
 x=FIELDS[ids];charge=CHARGE[ids];occupied=x<1.3;p=occupied.mean(0);minimum=x.min(0)
 # Strict common free points: never inside the carbon-sized atom-center exclusion envelope.
 skin=(minimum>=1.3)&(minimum<=4.3)
 band=(x.min(0)<4.3)&(x.max(0)>-.5)
 means=np.full(len(GRID),np.nan,dtype=np.float32)
 if np.isfinite(charge).any():means=np.mean(charge,axis=0).astype(np.float32)
 gr=groups(ids);unique=np.unique(gr);ng=len(unique)
 row.update(independent_groups=ng,unique_chains=len(set(SITES[i]['chain'] for i in ids)),strict_common_near_interface_voxels=int(skin.sum()),strict_common_near_interface_volume_A3=float(skin.sum()),charge_status='NEAREST_ATOM_PARTIAL_CHARGE_ANNOTATION_ONLY' if np.isfinite(means).any() else 'NOT_EVALUATED_UNCAPPED_CHARGE')
 arrays=dict(occupancy_probability=p.astype(np.float32),minimum_vdw_clearance_A=minimum.astype(np.float32),strict_common_near_interface_mask=skin,mean_partial_charge_annotation_e=means,site_indices=np.array(ids))
 seed=int.from_bytes(hashlib.sha256(name.encode()).digest()[:4],'little');rng=np.random.default_rng(seed)
 if ng>=5:
  draws=rng.integers(0,ng,size=(100,ng));w=np.array([(gr[None,:]==unique[d,None]).sum(0) for d in draws],dtype=float)
  w/=w.sum(1,keepdims=True)
  bp=w@occupied.astype(float)
  low,high=np.quantile(bp,[.05,.95],axis=0)
  row.update(bootstrap_status='EVALUATED_100_CLUSTER_RESAMPLES',bootstrap_mean_absolute_occupancy_shift=float(np.mean(np.abs(bp[:,band]-p[band]))) if band.any() else None,bootstrap_mean_90CI_width=float((high-low)[band].mean()) if band.any() else None)
  arrays.update(occupancy_bootstrap_q05=low.astype(np.float32),occupancy_bootstrap_q95=high.astype(np.float32))
 else:row['bootstrap_status']='NOT_EVALUATED_LT5_INDEPENDENT_GROUPS'
 # Leave one independent group out; remove training patches within 30 A of any test center.
 cv=[];centers=np.array([SITES[i]['center'] for i in ids]);box=BOX[surface]
 for g in unique:
  test=np.flatnonzero(gr==g);train=np.flatnonzero(gr!=g)
  if not len(train):continue
  diff=centers[train,None]-centers[test][None]
  if not surface.endswith('trimmed'):diff[:,:,:2]-=np.round(diff[:,:,:2]/box[:2])*box[:2]
  dist=np.linalg.norm(diff,axis=-1);train=train[(dist>=30).all(1)]
  if len(set(gr[train]))<3:continue
  lo=x[train].min(0);allowed=(lo>=1.3)&(lo<=4.3)
  if allowed.sum()<10:continue
  # How much of the fitted near-interface allowed volume clashes with held-out polymers?
  violations=(x[test][:,allowed]<1.3).mean(1)
  cv.append(dict(test_sites=len(test),training_sites=len(train),training_groups=len(set(gr[train])),allowed_voxels=int(allowed.sum()),heldout_clash_fraction_mean=float(violations.mean()),heldout_clash_fraction_max=float(violations.max())))
 row.update(cv_folds=len(cv),cv_status='EVALUATED_SPATIALLY_PURGED_GROUP_HOLDOUT' if cv else 'NOT_EVALUATED_INSUFFICIENT_SEPARATED_GROUPS')
 if cv:row.update(heldout_clash_fraction_mean=float(np.mean([v['heldout_clash_fraction_mean'] for v in cv])),heldout_clash_fraction_worst=float(max(v['heldout_clash_fraction_max'] for v in cv)))
 row['cv_details']=cv;row['status']='FIT_COMPLETE_CONDITIONAL_GEOMETRY'
 np.savez_compressed(O/'maps'/(name+'.npz'),**arrays)
 return row
if __name__=='__main__':
 O.mkdir(exist_ok=False);(O/'maps').mkdir()
 SITES=json.loads((P/'sites.json').read_text());PANEL=json.loads((P/'panel.json').read_text());COV=np.load(P/'coverage.npz')['coverage'];GRID=np.load(P/'grid.npz')['points_A']
 BOX={}
 for sf in set(s['surface'] for s in SITES):
  mat=sf.split('_')[0];BOX[sf]=np.load(R/'inputs_v2'/(mat+'_parent_atoms.npz'))['box_A']
 FIELDS=np.full((len(SITES),len(GRID)),np.nan,np.float32);CHARGE=FIELDS.copy()
 for i in np.flatnonzero(COV.any(0)):
  z=np.load(P/'fields'/(str(i)+'.npz'));FIELDS[i]=z['polymer_vdw_clearance_A'];CHARGE[i]=z['nearest_surface_atom_partial_charge_e']
 tasks=[];scope=[]
 for sf in sorted(set(s['surface'] for s in SITES)):
  material=sf.split('_')[0];panelids=[j for j,w in enumerate(PANEL) if (w['template']=='NylC')!=(material=='PET')]
  for t in [.7,.8,.9]:
   target=np.array([s['surface']==sf and s['support']>=t for s in SITES]);suffix=sf+'_ge'+str(round(t*100))
   for j in panelids:
    ids=np.flatnonzero(target&COV[j])
    tasks.append((PANEL[j]['id']+'_'+suffix,sf,t,ids,PANEL[j]['kind']))
   for kind in ['single_residue_dual_donor','two_residue_donors','both_donor_classes']:
    pp=[j for j in panelids if kind=='both_donor_classes' or PANEL[j]['kind']==kind]
    union=target&COV[pp].any(0);intersection=target&COV[pp].all(0)
    # Sites appear once, even if supported by many selected cores.
    tasks.append(('POOLED_'+kind+'_'+suffix,sf,t,np.flatnonzero(union),'pooled_'+kind))
    scope.append(dict(surface=sf,min_support=t,kind=kind,selected_cores=len(pp),union_sites=int(union.sum()),shared_by_all_selected_cores_sites=int(intersection.sum()),interpretation='pooled common exclusion constraints; not an averaged catalytic core or proof every core fits all union sites'))
 print('FITS',len(tasks),flush=True)
 with concurrent.futures.ProcessPoolExecutor(max_workers=8) as pool:
  result=[]
  for i,row in enumerate(pool.map(one,tasks),1):
   result.append(row)
   if i%100==0:print('FIT_PROGRESS',i,len(tasks),flush=True)
 summary=dict(fits=len(result),nonempty_fits=sum(r['sites']>0 for r in result),bootstrap_evaluated=sum(r.get('bootstrap_status')=='EVALUATED_100_CLUSTER_RESAMPLES' for r in result),holdout_evaluated=sum(r.get('cv_folds',0)>0 for r in result),grid_A=1.,patch_radius_A=15.,protein_atom_center_model=dict(radius_A=1.7,overlap_A=.4,exclusion_clearance_A=1.3),grouping='connected components sharing original chain or same side/30A XY tile',holdout='leave group out plus 30A patch-overlap purge; >=3 independent training groups; >=10 fitted interface voxels',bootstrap='100 cluster resamples; >=5 independent groups',scope='conditional fit stability after historical full-data core selection; NOT independent model-selection validation; sampled atom-center admissible volume, NOT a generated protein surface; charges are nearest-atom annotations, NOT electrostatic potential; uncapped charge maps withheld',source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
 (O/'results.json').write_text(json.dumps(result,indent=2));(O/'pooled_scopes.json').write_text(json.dumps(scope,indent=2));(O/'summary.json').write_text(json.dumps(summary,indent=2))
 with (R/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__)),source_sha256=summary['source_sha256'],outputs=str(O),exit_code=0,summary=summary))+'\n')
 print(json.dumps(summary),flush=True)
