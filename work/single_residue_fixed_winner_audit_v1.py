#!/usr/bin/env python3
"""Independent full-environment replay of reported fixed-core winners and nylon parent controls."""
from pathlib import Path
import json,sys,importlib.util,datetime,concurrent.futures
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
R=Path(__file__).resolve().parent;S=R/'surface_scan_v1';F=R/'fixed_coupled_transfer_v1'
spec=importlib.util.spec_from_file_location('surface',R/'surface_scan_v1.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def clash(q,r,env):
 margin=np.full(len(q),np.inf)
 for e,tree in env.items():margin=np.minimum(margin,tree.query(q)[0]-r-m.RAD[e]+.4)
 return bool((margin>=-1e-8).all())
def check(i):
 item=ROWS[i];mat=item['material'];row=item['row'];key=mat+'_'+row['site_id']
 z=np.load(S/'cases'/(key+'.npz'));center=z['origin_A'];frame=z['frame'];a=z['amide_A']
 xyz,el=ENVS[mat];q=xyz-center
 if mat=='PET':q[:,:2]-=np.round(q[:,:2]/m.BOX[:2])*m.BOX[:2]
 env={e:cKDTree((q@frame)[el==e]) for e in np.unique(el)}
 parent=None
 if mat in PARENTS:
  px,pe=PARENTS[mat];pq=(px-center)@frame
  parent={e:cKDTree(pq[pe==e]) for e in np.unique(pe)}
 yes=[];pyes=[]
 for j,w in enumerate(WIN):
  if (w['template']=='NylC')==(mat=='PET'):continue
  d=w['donor'];h=w['H'];dist=np.linalg.norm(d-a[1],axis=1)
  cosine=np.sum((d-h)*(a[1]-h),axis=1)/(np.linalg.norm(d-h,axis=1)*np.linalg.norm(a[1]-h,axis=1))
  lo=np.where(w['donor_elements']=='S',3.1,2.7)
  chemical=bool(((dist>=lo-1e-8)&(dist<=lo+.5+1e-8)&(cosine<=np.cos(np.deg2rad(140))+1e-12)).all())
  if chemical and clash(w['q'],w['radii'],env):
   yes.append(j)
   if parent is not None and clash(w['q'],w['radii'],parent):pyes.append(j)
 return i,yes,pyes
if __name__=='__main__':
 summary=json.loads((F/'summary.json').read_text());ROWS,ENVS,LOOKUP=m.get_inputs();temps=m.templates()
 tmeta=json.loads((F/'triads.json').read_text());tz=np.load(F/'triad_coordinates.npz')
 triads={}
 for name in ['IsPETase','LCC','NylC']:
  triads.update({int(i):q for i,q in zip(tz[name+'_indices'],tz[name])})
 cloud={key:dict(np.load(S/'clouds'/(key+'.npz'))) for key in summary['microstates']}
 picked={}
 for cohort in summary['cohorts']:
  for rr in cohort['per_microstate']:
   if 'best_balanced' in rr:
    w=rr['best_balanced'];picked[w['candidate_id']]=w
 WIN=[]
 for w in picked.values():
  w=w.copy();c=cloud[w['microstate']];sid=w['donor_sample_id'];t=triads[w['triad_index']]
  te=temps[w['template']][1];r1=np.array([m.RAD[e] for e in te]);r2=np.array([m.RAD[e] for e in c['elements']])
  assert (cdist(t,c['q_A'][sid])-r1[:,None]-r2[None]+.4>=-1e-8).all()
  w.update(q=np.vstack((t,c['q_A'][sid])),radii=np.r_[r1,r2],donor=c['donor_A'][sid],H=c['H_A'][sid],donor_elements=c['donor_elements'][sid]);WIN.append(w)
 sys.path.insert(0,str(m.N.parent/'surface_previews/axis_density_sixface_all_external_trim_v6'))
 from density_cut_variants import parse_pdb
 from terminal_density_trim_3d import center_whole_chains_in_primary_cell
 PARENTS={}
 for mat,dp in [('PA6',16),('PA66',8)]:
  p=m.N.parent/f'{mat.lower()}_dp{dp}_400chain/water_slab_v1/dry_export/{mat}_DP{dp}_400chain_water_equilibrated_dry.pdb'
  _,box,atoms,_=parse_pdb(p);center_whole_chains_in_primary_cell(atoms,box)
  atoms=[a for a in atoms if a['element'] in m.RAD]
  px=np.array([a['xyz'] for a in atoms]);pe=np.array([a['element'] for a in atoms])
  maximum=float(cKDTree(px).query(ENVS[mat][0])[0].max());assert maximum<1e-5,(mat,maximum)
  PARENTS[mat]=(px,pe)
 coverage=np.zeros((len(WIN),len(ROWS)),bool);parentcoverage=coverage.copy()
 with concurrent.futures.ProcessPoolExecutor(max_workers=8) as pool:
  for i,ids,pids in pool.map(check,range(len(ROWS))):coverage[ids,i]=True;parentcoverage[pids,i]=True
 assert not np.any(parentcoverage&~coverage)
 index={w['candidate_id']:i for i,w in enumerate(WIN)};audits=[]
 for co in summary['cohorts']:
  mask=np.array([(it['material']==co['material'] or co['material']=='NYLON_COMBINED' and it['material'] in ['PA6','PA66']) and float(it['row'].get('support_fraction',it['row'].get('inward_support')))>=co['min_support'] for it in ROWS])
  for r in co['per_microstate']:
   if 'best_balanced' not in r:continue
   w=r['best_balanced'];wi=index[w['candidate_id']];actual=int(coverage[wi,mask].sum())
   assert actual==w['coverage'],(co['name'],w['candidate_id'],actual,w['coverage'])
   audits.append(dict(cohort=co['name'],microstate=w['microstate'],candidate_id=w['candidate_id'],coverage=actual,parent_control_coverage=int(parentcoverage[wi,mask].sum()) if co['material']!='PET' else None))
 np.savez_compressed(F/'winner_direct_replay.npz',coverage=coverage,parent_coverage=parentcoverage,candidate_ids=np.array(list(index)))
 out=dict(status='DIRECT_FULL_ENVIRONMENT_REPLAY_PASS',unique_winners=len(WIN),checks=len(audits),audits=audits,parent_scope='nylon finite centered original untrimmed dry parent; exact retained-coordinate correspondence checked; parent control not a periodic ingress simulation',source_sha256=m.sha(Path(__file__)))
 (F/'winner_direct_audit.json').write_text(json.dumps(out,indent=2))
 with (R/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__)),source_sha256=m.sha(Path(__file__)),outputs=str(F/'winner_direct_audit.json'),exit_code=0,summary=dict(winners=len(WIN),checks=len(audits))))+'\n')
 print(json.dumps(out),flush=True)
