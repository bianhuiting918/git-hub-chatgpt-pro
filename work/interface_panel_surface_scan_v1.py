#!/usr/bin/env python3
"""Expanded fixed-core panel, new original-surface transfer, and cached local fields."""
from pathlib import Path
import json,importlib.util,datetime,concurrent.futures,sys
import numpy as np
from scipy.spatial import cKDTree
R=Path(__file__).resolve().parent;B=R.parent;I=R/'inputs_v2'
S=B/'single_residue_dual_donor_20260908_v1';N=B/'rectangular_400chain_v1/nylc_tdd_surface_v1';P=B/'pet_dp10_400chain_direct_v1';O=R/'panel_surface_scan_v1'
spec=importlib.util.spec_from_file_location('single',S/'surface_scan_v1.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def panel():
 out={};templates=m.templates()
 ss=json.loads((S/'fixed_coupled_transfer_v1/summary.json').read_text());tz=np.load(S/'fixed_coupled_transfer_v1/triad_coordinates.npz')
 tq={int(i):q for name in ['IsPETase','LCC','NylC'] for i,q in zip(tz[name+'_indices'],tz[name])}
 for co in ss['cohorts']:
  if co['name'] not in ['PET_ge70','PA6_ge70','PA66_ge70','NYLON_COMBINED_ge70']:continue
  selected=co['top10']+[v['best_balanced'] for v in co['per_microstate'] if 'best_balanced' in v]
  for w in selected:
   key=w['candidate_id'];c=dict(np.load(S/'surface_scan_v1/clouds'/(w['microstate']+'.npz')));j=w['donor_sample_id'];tri=tq[w['triad_index']]
   if key not in out:out[key]=dict(id=key,kind='single_residue_dual_donor',microstate=w['microstate'],template=w['template'],q=np.vstack((tri,c['q_A'][j])),elements=np.r_[templates[w['template']][1],c['elements']],donor=c['donor_A'][j],H=c['H_A'][j],donor_elements=c['donor_elements'][j],selected_by=[])
   out[key]['selected_by'].append(co['name']+' top10 or chemistry diversity')
 pm=json.loads((P/'donor_geometry_batch_singleNH_v1/joint_transfer_map_v1/motifs.json').read_text())
 pz=np.load(P/'donor_geometry_batch_singleNH_v1/joint_transfer_map_v1/joint_transfer_results.npz')
 rankings=json.loads((P/'donor_pair_existing_density_1A_weighted_v2/all473_rankings.json').read_text())
 chosen=[w['id'] for w in rankings[:10]]+[pm[i]['id'] for i in np.argsort(-pz['coverage'].sum(1),kind='stable')[:10]]
 for key in dict.fromkeys(chosen):
  j=next(i for i,v in enumerate(pm) if v['id']==key);w=pm[j];q=pz['motif_xyz_A'][j]
  out[key]=dict(id=key,kind='two_residue_donors',microstate='two_backbone_NH',template=w['template'],q=q,elements=np.array(w['elements']),donor=q[[26,32]],H=pz['motif_H_A'][j],donor_elements=np.array(['N','N']),selected_by=['PET legacy top10 and coverage top10; legacy pooled density not validation'])
 nm=json.loads((N/'fixed_core_library_v1/motifs.json').read_text());nz=np.load(N/'fixed_core_library_v1/core_coordinates.npz')
 nr=json.loads((N/'conditional_density_v1/rankings.json').read_text())
 # Generic original peptide donor element order: C,O,N,C,C,O.
 for group,ranking in nr.items():
  if not group.endswith('support0.7'):continue
  for w in ranking[:10]:
   key=w['id'];j=next(i for i,v in enumerate(nm) if v['id']==key);q=nz['cores_A'][j]
   if key not in out:out[key]=dict(id=key,kind='two_residue_donors',microstate='two_backbone_NH',template='NylC',q=q,elements=np.r_[templates['NylC'][1],np.array(['C','O','N','C','C','O']*2)],donor=q[[25,31]],H=nz['donor_H_A'][j],donor_elements=np.array(['N','N']),selected_by=[])
   out[key]['selected_by'].append(group+' top10')
 return list(out.values())
def inputs():
 oldrows,_,_=m.get_inputs();sites=[];environments={}
 for mat in ['PET','PA6','PA66']:
  z=dict(np.load(I/(mat+'_parent_atoms.npz')))
  environments[mat+'_original']=dict(xyz=z['xyz_A'],elements=z['elements'],charge=z['charges_e'],box=z['box_A'],periodic=True)
  if mat=='PET':
   for item in oldrows:
    if item['material']!='PET':continue
    row=item['row'];zz=np.load(S/'surface_scan_v1/cases'/('PET_'+row['site_id']+'.npz'))
    sites.append(dict(surface='PET_original',material='PET',site_id=row['site_id'],chain=row['segment'],side=row['side'],support=row['support_fraction'],center=zz['origin_A'],frame=zz['frame'],amide=zz['amide_A']))
  else:
   for row in json.loads((I/(mat+'_original_surface_sites.json')).read_text()):
    a=np.array(row['amide_A']);center=a[0];ex=m.norm(a[1]-center);ey=a[2]-center;ey=m.norm(ey-ex*np.dot(ex,ey));f=np.column_stack((ex,ey,np.cross(ex,ey)))
    sites.append(dict(surface=mat+'_original',material=mat,site_id=row['site_id'],chain='P'+str(row['original_chain']).zfill(3),side=row['side'],support=row['support'],center=center,frame=f,amide=(a-center)@f))
   z=dict(np.load(I/(mat+'_trimmed_atoms.npz')));environments[mat+'_trimmed']=dict(xyz=z['xyz_A'],elements=z['elements'],charge=np.full(len(z['xyz_A']),np.nan),box=z['box_A'],periodic=False)
   for item in oldrows:
    if item['material']!=mat:continue
    row=item['row'];zz=np.load(S/'surface_scan_v1/cases'/(mat+'_'+row['site_id']+'.npz'))
    sites.append(dict(surface=mat+'_trimmed',material=mat,site_id=row['site_id'],chain='P'+str(row['original_chain']).zfill(3),side=row['nearest_side'],support=row['inward_support'],center=zz['origin_A'],frame=zz['frame'],amide=zz['amide_A']))
 return sites,environments
def evaluate(i):
 s=SITES[i];e=ENV[s['surface']];q=e['xyz']-s['center']
 if e['periodic']:q[:,:2]-=np.round(q[:,:2]/e['box'][:2])*e['box'][:2]
 heavy=e['elements']!='H';local=(np.linalg.norm(q,axis=1)<=32);ids=np.flatnonzero(heavy&local);xyz=q@s['frame']
 env=m.trees(xyz[ids],e['elements'][ids]);passed=[]
 for j,w in enumerate(PANEL):
  if (w['template']=='NylC')==(s['material']=='PET'):continue
  if not m.chemical(w['donor'][None],w['H'][None],s['amide'][1],w['donor_elements'][None])[0]:continue
  if m.clear(w['q'][None],np.array([m.RAD[t] for t in w['elements']]),env)[0]>=-1e-8:passed.append(j)
 if passed:
  clearance=np.full(len(GRID),np.inf)
  for t,tree in env.items():clearance=np.minimum(clearance,tree.query(GRID)[0]-m.RAD[t])
  # Partial-charge annotation, NOT an electrostatic potential or a dielectric calculation.
  qp=np.full(len(GRID),np.nan)
  if e['periodic']:
   ai=np.flatnonzero(local);best=np.full(len(GRID),np.inf)
   for t in np.unique(e['elements'][ai]):
    ti=ai[e['elements'][ai]==t];d,ix=cKDTree(xyz[ti]).query(GRID);score=d-({'H':1.2,**m.RAD}[str(t)])
    yes=score<best;best[yes]=score[yes];qp[yes]=e['charge'][ti[ix[yes]]]
  np.savez_compressed(O/'fields'/(str(i)+'.npz'),polymer_vdw_clearance_A=clearance.astype(np.float32),nearest_surface_atom_partial_charge_e=qp.astype(np.float32))
 return i,passed
if __name__=='__main__':
 O.mkdir(exist_ok=False);(O/'fields').mkdir();PANEL=panel();SITES,ENV=inputs()
 axis=np.arange(-15,15.0001,1.);full=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),-1).reshape(-1,3);inside=np.linalg.norm(full,axis=1)<=15
 GRID=full[inside];np.savez_compressed(O/'grid.npz',points_A=GRID,axis_A=axis,sphere_mask=inside)
 assert all(np.linalg.norm(w['q'],axis=1).max()<28 for w in PANEL)
 np.savez_compressed(O/'cores.npz',**{w['id']+'_q':w['q'] for w in PANEL},**{w['id']+'_donors':w['donor'] for w in PANEL},**{w['id']+'_H':w['H'] for w in PANEL})
 (O/'panel.json').write_text(json.dumps([{k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in w.items() if k not in ['q','donor','H']} for w in PANEL],indent=2))
 (O/'sites.json').write_text(json.dumps([{k:v.tolist() if isinstance(v,np.ndarray) else v for k,v in s.items()} for s in SITES],separators=(',',':')))
 print('PANEL',len(PANEL),'SITES',len(SITES),'GRID',len(GRID),flush=True)
 cov=np.zeros((len(PANEL),len(SITES)),bool)
 with concurrent.futures.ProcessPoolExecutor(max_workers=8) as pool:
  for count,(i,yes) in enumerate(pool.map(evaluate,range(len(SITES))),1):
   cov[yes,i]=True
   if count%200==0:print('SCAN',count,len(SITES),flush=True)
 np.savez_compressed(O/'coverage.npz',coverage=cov)
 summary=dict(cores=len(PANEL),sites=len(SITES),grid_points=len(GRID),surface_counts={surface:sum(s['surface']==surface for s in SITES) for surface in ENV},cores_with_hit={surface:int(cov[:,np.array([s['surface']==surface for s in SITES])].any(1).sum()) for surface in ENV},sites_with_any_core=int(cov.any(0).sum()),scope='expanded fixed existing panel tested on newly selected original and existing trimmed sites; only passing sites have cached fields; partial charge annotation is not potential; trimmed charges NOT_EVALUATED',source_sha256=m.sha(Path(__file__)))
 (O/'summary.json').write_text(json.dumps(summary,indent=2))
 (R/'RUNBOOK.md').write_text('CPU-only existing environment. Source scripts tracked in codex/pet-surface-stats-20260907. Step 1 interface_inputs_v2.py audits charges and selects original nylon sites; v1 failed because a nylon-specific parser rejected PET, preserved. Step 2 panel_surface_scan_v1.py expands fixed candidate panels, evaluates original and trimmed sites independently, caches 1A-grid local polymer VDW-clearance fields in a 15A sphere. Partial-charge maps are nearest-surface-atom annotations, not PB/PME electrostatic potential; uncapped trimmed charge maps remain NaN. Original nylon is XY periodic. Field fitting and independent-group resampling are subsequent stages. Candidate selection used historical full-data rankings; subsequent fit validation is conditional on these selected cores, not an independent model-selection test.\n')
 with (R/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__)),source_sha256=m.sha(Path(__file__)),outputs=str(O),exit_code=0,summary=summary))+'\n')
 print(json.dumps(summary),flush=True)
