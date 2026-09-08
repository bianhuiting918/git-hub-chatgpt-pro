#!/usr/bin/env python3
"""Transfer fixed triad+coupled residue motifs unchanged; count existing-neighbor incidences."""
from pathlib import Path
import json,importlib.util,datetime,concurrent.futures,time
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import csr_matrix
R=Path(__file__).resolve().parent;S=R/'surface_scan_v1';O=R/'fixed_coupled_transfer_v1'
spec=importlib.util.spec_from_file_location('stage1',R/'surface_scan_v1.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
def target(i):
 item=ROWS[i];mat=item['material'];row=item['row'];xyz,el=ENVS[mat]
 center=SITE_ARRAYS[i]['origin_A'];frame=SITE_ARRAYS[i]['frame']
 q=xyz-center
 if mat=='PET':q[:,:2]-=np.round(q[:,:2]/m.BOX[:2])*m.BOX[:2]
 keep=np.linalg.norm(q,axis=1)<=30;env=m.trees(q[keep]@frame,el[keep])
 ans=np.zeros(len(TRIADS),bool)
 for name in ['IsPETase','LCC','NylC']:
  if (name=='NylC')==(mat=='PET'):continue
  ids=np.array([j for j,t in enumerate(TRIADS) if t['template']==name],int)
  if len(ids):
   coords=np.array([TRIADS[j]['q'] for j in ids]);r=np.array([m.RAD[e] for e in TEMPLATES[name][1]])
   for b in range(0,len(ids),256):ans[ids[b:b+256]]=m.clear(coords[b:b+256],r,env)>=-1e-8
 return i,ans
if __name__=='__main__':
 assert json.loads((S/'stage1_summary.json').read_text())['completed_sites']==840
 O.mkdir(exist_ok=False);start=time.monotonic()
 ROWS,ENVS,LOOKUP=m.get_inputs();TEMPLATES=m.templates()
 SITE_ARRAYS=[];CASES=[];TRIADS=[]
 keys=sorted(p.stem for p in (S/'clouds').glob('*.npz'))
 clouds={k:dict(np.load(S/'clouds'/(k+'.npz'))) for k in keys}
 poly={k:np.zeros((len(ROWS),8192),bool) for k in keys}
 for i,item in enumerate(ROWS):
  key=item['material']+'_'+item['row']['site_id']
  rec=json.loads((S/'cases'/(key+'.json')).read_text());z=np.load(S/'cases'/(key+'.npz'))
  SITE_ARRAYS.append(dict(origin_A=z['origin_A'],frame=z['frame']))
  CASES.append(rec)
  for k in keys:poly[k][i,z[k+'_polymer_ids']]=True
  for name in rec['triad_stats']:
   for j,q in enumerate(z[name+'_triads_A']):
    # Retain all triad poses that support at least one coupled-residue candidate.
    subsets={r['microstate']:z[f"{r['microstate']}_{name}_{j}_ids"] for r in rec['results'] if r['template']==name and r['accepted_per_triad'][j]>0}
    if subsets:TRIADS.append(dict(template=name,source_index=i,source_pose=j,q=q,subsets=subsets))
 print('INPUTS',len(ROWS),'triads',len(TRIADS),'source_coupled_pairs',sum(len(s) for t in TRIADS for s in t['subsets'].values()),flush=True)
 triad_cover=np.zeros((len(ROWS),len(TRIADS)),bool)
 with concurrent.futures.ProcessPoolExecutor(max_workers=8) as pool:
  for count,(i,ans) in enumerate(pool.map(target,range(len(ROWS))),1):
   triad_cover[i]=ans
   if count%100==0:print('TRIAD_TRANSFER',count,len(ROWS),flush=True)
 assert all(triad_cover[t['source_index'],j] for j,t in enumerate(TRIADS))
 np.savez_compressed(O/'triad_coverage.npz',coverage=triad_cover)
 cohorts=[]
 for mat in ['PET','PA6','PA66','NYLON_COMBINED']:
  for threshold in [.9,.8,.7]:
   mask=np.array([(r['material']==mat or mat=='NYLON_COMBINED' and r['material'] in ['PA6','PA66']) and r['support']>=threshold for r in CASES])
   cohorts.append(dict(name=mat+'_ge'+str(round(threshold*100)),material=mat,min_support=threshold,denominator=int(mask.sum()),mask=mask))
 masks=np.column_stack([c['mask'] for c in cohorts])
 covs=[];dens=[];triadids=[];microids=[];sampleids=[]
 for ti,t in enumerate(TRIADS):
  for ki,key in enumerate(keys):
   if key not in t['subsets']:continue
   ids=t['subsets'][key];c=clouds[key];si=t['source_index']
   assert poly[key][si,ids].all()
   # Each candidate is one unchanged complete triad + unchanged complete residue.
   incidence=(triad_cover[:,ti,None]&poly[key][:,ids])
   weights=incidence.astype(np.int32).T@masks.astype(np.int32)
   assert (weights>=0).all()
   density=np.zeros((len(ids),2,len(cohorts)),np.int64)
   pair=c['pair_indices'][ids]
   for d in range(2):
    tree=cKDTree(c['donor_A'][ids,d])
    coo=tree.sparse_distance_matrix(tree,1.,output_type='coo_matrix')
    same=(pair[coo.row]==pair[coo.col]).all(1)
    adj=csr_matrix((np.ones(int(same.sum()),np.int32),(coo.row[same],coo.col[same])),shape=(len(ids),len(ids)))
    density[:,d]=adj@weights
   assert (density>=weights[:,None]).all(),'Center must contribute its own site incidences'
   covs.append(weights);dens.append(density);triadids.extend([ti]*len(ids));microids.extend([ki]*len(ids));sampleids.extend(ids.tolist())
  if (ti+1)%500==0:print('FIXED_CORE_METRICS',ti+1,len(TRIADS),flush=True)
 cov=np.concatenate(covs);den=np.concatenate(dens);triadids=np.array(triadids);microids=np.array(microids);sampleids=np.array(sampleids)
 np.savez_compressed(O/'all_fixed_candidates.npz',triad_indices=triadids,microstate_indices=microids,sample_indices=sampleids,coverage_counts=cov,neighbor_site_incidence_counts=den)
 metadata=[dict(template=t['template'],source_index=t['source_index'],source_material=ROWS[t['source_index']]['material'],source_site=ROWS[t['source_index']]['row']['site_id'],source_pose=t['source_pose']) for t in TRIADS]
 (O/'triads.json').write_text(json.dumps(metadata,separators=(',',':')))
 np.savez_compressed(O/'triad_coordinates.npz',**{n:np.array([t['q'] for t in TRIADS if t['template']==n]) for n in ['IsPETase','LCC','NylC']},**{n+'_indices':np.array([j for j,t in enumerate(TRIADS) if t['template']==n]) for n in ['IsPETase','LCC','NylC']})
 reports=[]
 for ci,cohort in enumerate(cohorts):
  robust=den[:,:,ci].min(1);mx=max(1,int(cov[:,ci].max()));md=max(1,int(robust.max()))
  score=.6*cov[:,ci]/mx+.4*robust/md
  def record(j):
   ti=int(triadids[j]);t=metadata[ti];key=keys[microids[j]];sid=int(sampleids[j])
   return dict(candidate_id='SC'+str(j+1).zfill(7),candidate_index=int(j),microstate=key,triad_index=ti,donor_sample_id=sid,template=t['template'],source_material=t['source_material'],source_site=t['source_site'],source_pose=t['source_pose'],coverage=int(cov[j,ci]),denominator=cohort['denominator'],donor_neighbor_site_incidences=den[j,:,ci].tolist(),robust_min=int(robust[j]),score=float(score[j]),donor_atom_pair=clouds[key]['pair_indices'][sid].tolist())
  eligible=np.flatnonzero(cov[:,ci]>0)
  order=eligible[np.lexsort((eligible,-robust[eligible],-cov[eligible,ci],-score[eligible]))]
  per=[]
  for ki,key in enumerate(keys):
   ix=eligible[microids[eligible]==ki]
   if not len(ix):per.append(dict(microstate=key,status='NO_HIT_IN_THIS_FIXED_LIBRARY',max_coverage=0));continue
   ix=ix[np.lexsort((ix,-robust[ix],-cov[ix,ci],-score[ix]))]
   per.append(dict(microstate=key,fixed_candidates=int((microids==ki).sum()),max_coverage=int(cov[ix,ci].max()),best_balanced=record(ix[0])))
  reports.append(dict(**{k:v for k,v in cohort.items() if k!='mask'},normalization=dict(max_coverage=mx,max_robust_min=md),top10=[record(j) for j in order[:10]],per_microstate=per))
 summary=dict(fixed_candidate_count=len(cov),triad_pose_count=len(TRIADS),microstates=keys,cohorts=reports,scope='all source-accepted triad+coupled-residue pairs in this finite cloud; unchanged coordinates transferred in C-based frame; robust density restricted to same-source accepted cloud, same triad and donor atom pair, weighted by target-site incidence; no new perturbations; not all possible conformations; no activity prediction; PET vs nylon catalytic cores not pooled',source_sha256=m.sha(Path(__file__)),elapsed_seconds=time.monotonic()-start)
 (O/'summary.json').write_text(json.dumps(summary,indent=2))
 (O/'RUNBOOK.md').write_text('Run after surface_scan_v1.py --run completes all 840 sites. Existing CPU Python, OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1, eight workers. Fixed candidates are all source-accepted coupled residue poses under a retained source triad. Rigid transfer is factorized: triad/polymer compatibility AND residue/polymer/chemical compatibility; internal triad/residue clash does not change under rigid transfer. Density uses same-source accepted existing candidates only and counts each compatible target incidence; no coordinate dedup or new perturbations. Counts are conditional on finite unrelaxed torsion library. Score = .6 normalized coverage + .4 normalized min donor density, normalized within each cohort. Cys sulfur and His microstates require separate chemical interpretation. All metrics in all_fixed_candidates.npz; row i corresponds SC(i+1).\n')
 with (R/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__)),source_sha256=m.sha(Path(__file__)),inputs=str(S),outputs=str(O),exit_code=0,summary=dict(candidates=len(cov),triads=len(TRIADS))))+'\n')
 print('FIXED_TRANSFER_COMPLETE',len(cov),flush=True)
