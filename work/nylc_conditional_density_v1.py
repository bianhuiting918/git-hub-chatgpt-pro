#!/usr/bin/env python3
"""Weighted existing-point neighborhoods conditional on the identical fixed TDD."""
import tdd_donor_pilot_v1 as p
import numpy as np,json,datetime
from collections import defaultdict
ROOT=p.ROOT
if __name__=='__main__':
 out=ROOT/'conditional_density_v1';out.mkdir(exist_ok=False);rc=1
 try:
  meta=json.loads((ROOT/'fixed_core_library_v1/motifs.json').read_text())
  poses=json.loads((ROOT/'fixed_core_library_v1/poses.json').read_text())
  sites=json.loads((ROOT/'fixed_core_transfer_v1/sites.json').read_text())
  lib=np.load(ROOT/'fixed_core_transfer_v1/coverage.npz');q=lib['cores_A'];cov=lib['coverage'];pcov=lib['parent_coverage']
  clouds=np.load(ROOT/'fixed_core_library_v1/conditional_existing_donor_clouds.npz')
  groups=defaultdict(list)
  for i,m in enumerate(meta):groups[m['pose_id']].append(i)
  triads=[];data=[]
  for pid in range(len(poses)):
   ids=groups[pid];d=clouds[f'pose_{pid}_donors_A'];h=clouds[f'pose_{pid}_H_A']
   centers=q[ids][:,[25,31]]
   membership=np.linalg.norm(d[None,None,:,2,:]-centers[:,:,None,:],axis=-1)<=1.+1e-10
   use=membership.any(axis=(0,1))
   triads.append(clouds[f'pose_{pid}_triad_A'])
   data.append((ids,d[use],h[use],membership[:,:,use]))
  triads=np.array(triads);tr=np.array([p.RAD[e] for e in p.te])
  counts=np.zeros((len(meta),len(sites),2),dtype=np.uint16)
  current=None
  for si,site in enumerate(sites):
   mat=site['material']
   if mat!=current:
    rec=p.load_pdb(ROOT/'inputs'/f'{mat}_v10_ge90_uncapped.pdb');xyz=np.array([a[1] for a in rec]);el=np.array([a[2] for a in rec])
    lookup={a['site_id']:a for a in json.loads((ROOT/'carbon_exposure_v1'/f'{mat}_sites.json').read_text())}
    atom={l[6:11]:v for l,v,e in rec};current=mat
   c=np.array(site['origin_A']);F=np.array(site['frame']);O=(atom[lookup[site['site_id']]['serials'][1]]-c)@F
   env=p.trees((xyz-c)@F,el);tpass=p.clear(triads,tr,env)>=-1e-8
   for pid in np.flatnonzero(tpass):
    ids,d,h,member=data[pid]
    if not len(d):continue
    N=d[:,2];dist=np.linalg.norm(N-O,axis=-1);v=N-h;w=O-h
    angle=np.degrees(np.arccos(np.clip(np.sum(v*w,axis=-1)/np.linalg.norm(v,axis=-1)/np.linalg.norm(w,axis=-1),-1,1)))
    ok=(dist>=2.7-1e-8)&(dist<=3.2+1e-8)&(angle>=140-1e-8)
    good=np.flatnonzero(ok)
    if len(good):ok[good]=p.clear(d[good],p.rr,env)>=-1e-8
    counts[ids,si,:]=(member&ok[None,None,:]).sum(axis=2)
   if (si+1)%25==0:print(json.dumps(dict(density_targets_done=si+1,total=len(sites))),flush=True)
  sourceindex={s['material']+'_'+s['site_id']:i for i,s in enumerate(sites)}
  assert all(np.array_equal(counts[i,sourceindex[m['source']]],m['source_density_1A']) for i,m in enumerate(meta)),'Source density replay'
  np.savez_compressed(out/'counts.npz',donor_neighbor_counts=counts)
  reports={};ranks={}
  supports=np.array([s['support'] for s in sites])
  for material in ('PA6','PA66','combined'):
   for threshold in (.9,.8,.7):
    mask=np.array([material=='combined' or s['material']==material for s in sites])&(supports>=threshold)
    C=cov[:,mask].sum(1);D=counts[:,mask].sum(1,dtype=np.int64);weak=D.min(1)
    cmax=int(C.max());dmax=int(weak.max());score=.6*C/max(cmax,1)+.4*weak/max(dmax,1)
    order=np.lexsort((np.arange(len(meta)),-score));name=material+'_support'+str(threshold)
    rows=[dict(meta[i],coverage=int(C[i]),denominator=int(mask.sum()),neighbor_counts=D[i].tolist(),weak_neighbor_count=int(weak[i]),score=float(score[i]),parent_coverage=int(pcov[i,mask].sum())) for i in order]
    ranks[name]=rows;reports[name]=dict(denominator=int(mask.sum()),max_coverage=cmax,max_weak_density=dmax,best=rows[0])
  (out/'rankings.json').write_text(json.dumps(ranks,indent=2))
  report=dict(groups=reports,definition='sum existing donor-placement/site incidences within each fixed N 1A; same fixed TDD; source accepted cloud only; center included; no coordinate dedup; single donor validity not joint neighbor-pair proof',score='0.6 normalized fixed-core site coverage + 0.4 normalized min(D1,D2), separately normalized within each cohort',scope='finite representative source libraries; no new sampling; no physical robustness or activity inference')
  (out/'summary.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True);rc=0
 finally:
  with (ROOT/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(p.Path(__file__).resolve()),script_sha256=p.sha(p.Path(__file__)),command='CPU Python conditional_density_v1.py',inputs='fixed cores, coverage, existing same-TDD donor clouds and v10 targets',outputs=str(out),parameters=dict(radius_A=1,coverage_weight=.6,density_weight=.4),exit_code=rc))+'\n')
