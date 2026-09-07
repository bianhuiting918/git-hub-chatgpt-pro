#!/usr/bin/env python3
"""Lossless compact full-site continuation of the eight-site pilot, unchanged geometry."""
import tdd_donor_pilot_v1 as p
from pathlib import Path
import numpy as np,json,datetime,time
ROOT=p.ROOT
def angle(a,b):return np.degrees(np.arccos(np.clip(np.sum(a*b,axis=-1)/np.linalg.norm(a,axis=-1)/np.linalg.norm(b,axis=-1),-1,1)))
def replay():
 report=[]
 for row in json.loads((ROOT/'tdd_donor_pilot_v1/summary.json').read_text()):
  key=row['material']+'_'+row['site'];z=np.load(ROOT/'tdd_donor_pilot_v1'/(key+'.npz'))
  pairs=z['pairs'];n=len(pairs)
  if not n:continue
  ids=np.unique(np.linspace(0,n-1,min(n,31),dtype=int))
  rec=p.load_pdb(ROOT/'inputs'/(row['material']+'_v10_ge90_uncapped.pdb'))
  env=p.trees(np.array([a[1] for a in rec]),np.array([a[2] for a in rec]))
  am=z['amide_A']; checks=[]
  for k in ids:
   j,a,b=pairs[k];tri=z['triads_A'][j];d1=z['donors_A'][a];d2=z['donors_A'][b]
   assert np.allclose(np.linalg.norm(tri[:,None]-tri[None,:],axis=-1),np.linalg.norm(p.tq[:,None]-p.tq[None,:],axis=-1),atol=1e-8)
   attack=np.linalg.norm(tri[p.og]-am[0]);ac=angle(am[1]-am[0],tri[p.og]-am[0]);ab=angle(tri[p.cb]-tri[p.og],am[0]-tri[p.og])
   assert 3-1e-8<=attack<=3.4+1e-8 and 95-1e-8<=ac<=115+1e-8 and 100-1e-8<=ab<=130+1e-8
   for d,h in [(d1,z['donor_H_A'][a]),(d2,z['donor_H_A'][b])]:
    no=np.linalg.norm(d[2]-am[1]);nh=angle(d[2]-h,am[1]-h)
    assert 2.7-1e-8<=no<=3.2+1e-8 and nh>=140-1e-8
    assert p.clear(d[None],p.rr,env)[0]>=-1e-8
    assert p.clear(d[None],p.rr,p.trees(tri,p.te))[0]>=-1e-8
   assert p.clear(tri[None],np.array([p.RAD[e] for e in p.te]),env)[0]>=-1e-8
   # Independent direct 36-distance pair calculation.
   margin=np.linalg.norm(d1[:,None]-d2[None,:],axis=-1)-p.rr[:,None]-p.rr[None,:]+.4
   assert margin.min()>=-1e-8
  report.append(dict(key=key,replayed_pairs=len(ids),status='PASS'))
 return report
def evaluate(mat,row,out,pilot_summary):
 key=mat+'_'+row['site_id']; dest=out/(key+'.npz');reportfile=out/(key+'.json')
 if reportfile.exists():return json.loads(reportfile.read_text())
 old=ROOT/'tdd_donor_pilot_v1'/(key+'.npz')
 if old.exists():
  z=np.load(old);poses=z['triads_A'];q=z['donors_A'];H=z['donor_H_A'];good=z['polymer_clear_donor_ids'];amide=z['amide_A'];frame=z['frame'];c=z['origin_A']
  accepted=[z[f'pose_{j}_donor_ids'] for j in range(len(poses))]
  stats=pilot_summary[key]['triad_stats'];reuse=True
 else:
  rec=p.load_pdb(ROOT/'inputs'/f'{mat}_v10_ge90_uncapped.pdb')
  lookup={l[6:11]:q for l,q,_ in rec};amide=np.array([lookup[k] for k in row['serials']]);c=amide[0]
  ex=p.norm(amide[1]-c);ey=amide[2]-c;ey=p.norm(ey-ex*np.dot(ex,ey));frame=np.column_stack((ex,ey,np.cross(ex,ey)))
  env=p.trees(np.array([v[1] for v in rec]),np.array([v[2] for v in rec]))
  p.rng=np.random.default_rng(p.seed_for(key+'|TDD|triad'))
  poses,stats=p.make_triads(amide,frame,p.template,env)
  q,H,*_=p.sample_donor(amide,frame,p.seed_for(key+'|NH|Sobol'))
  good=np.flatnonzero(p.clear(q,p.rr,env)>=0);cand=q[good]
  accepted=[good[np.flatnonzero(p.clear(cand,p.rr,p.trees(pose,p.te))>=0)] if len(cand) else np.array([],int) for pose in poses];reuse=False
 adj=p.compatibility(q[good]);rev={int(v):i for i,v in enumerate(good)}
 counts=[]
 for ids in accepted:
  local=np.array([rev[int(v)] for v in ids],dtype=int)
  counts.append(int(np.count_nonzero(adj[np.ix_(local,local)])//2))
 if reuse:
  assert sum(counts)==pilot_summary[key]['pair_records']
 np.savez_compressed(dest,triads_A=poses,donors_A=q,donor_H_A=H,amide_A=amide,frame=frame,origin_A=c,polymer_clear_donor_ids=good,pair_compatibility_packed=np.packbits(adj,axis=1),pair_compatibility_shape=np.array(adj.shape),pose_pair_counts=np.array(counts,dtype=np.int64),**{f'pose_{j}_donor_ids':ids for j,ids in enumerate(accepted)})
 result=dict(material=mat,site=row['site_id'],side=row['nearest_side'],support=row['inward_support'],parent_outward_points=row['parent_outward_points'],triad_stats=stats,donor_trials=len(q),polymer_clear_donors=len(good),paired_triad_poses=sum(n>0 for n in counts),pair_records=sum(counts),pilot_reused=reuse,status='GEOMETRY_PROXY_ONLY')
 reportfile.write_text(json.dumps(result));print(json.dumps(result),flush=True);return result
if __name__=='__main__':
 out=ROOT/'tdd_donor_all_v1';out.mkdir(exist_ok=True);rc=1
 start=datetime.datetime.now(datetime.timezone.utc).isoformat()
 with (ROOT/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=start,stage='full_scan_start',script=str(Path(__file__).resolve()),script_sha256=p.sha(Path(__file__)),inputs='262 candidates, unchanged pilot inputs and parameters',outputs=str(out),command='CPU Python tdd_donor_all_v1.py',exit_code=None))+'\n')
 try:
  rep=replay();(out/'pilot_replay.json').write_text(json.dumps(rep,indent=2))
  config=json.loads((ROOT/'tdd_donor_pilot_v1/config.json').read_text());config['storage']='lossless packed pair compatibility + per-pose accepted donor IDs';(out/'config.json').write_text(json.dumps(config,indent=2))
  old={r['material']+'_'+r['site']:r for r in json.loads((ROOT/'tdd_donor_pilot_v1/summary.json').read_text())}
  results=[]
  for mat in ('PA6','PA66'):
   for row in json.loads((ROOT/'carbon_exposure_v1'/f'{mat}_sites.json').read_text()):results.append(evaluate(mat,row,out,old))
  assert len(results)==262
  (out/'summary.json').write_text(json.dumps(results,indent=2));rc=0
 finally:
  with (ROOT/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),stage='full_scan_end',script=str(Path(__file__).resolve()),script_sha256=p.sha(Path(__file__)),outputs=str(out),exit_code=rc))+'\n')
