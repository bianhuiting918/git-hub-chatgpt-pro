#!/usr/bin/env python3
"""Fixed-core library from every paired TDD pose; no new sampling."""
import tdd_donor_pilot_v1 as p
import numpy as np,json,datetime
from scipy.spatial import cKDTree
ROOT=p.ROOT
if __name__=='__main__':
 out=ROOT/'fixed_core_library_v1';out.mkdir(exist_ok=False);rc=1
 try:
  summary=json.loads((ROOT/'tdd_donor_all_v1/summary.json').read_text())
  meta=[];coords=[];hydrogens=[];posemeta=[];clouds={}
  for row in summary:
   if not row['pair_records']:continue
   key=row['material']+'_'+row['site'];z=np.load(ROOT/'tdd_donor_all_v1'/(key+'.npz'))
   good=z['polymer_clear_donor_ids'];adj=np.unpackbits(z['pair_compatibility_packed'],axis=1,count=len(good)).astype(bool)
   assert np.array_equal(adj,adj.T) and not adj.diagonal().any()
   rev={int(v):i for i,v in enumerate(good)}
   frame=z['frame'];origin=z['origin_A']
   for j,count in enumerate(z['pose_pair_counts']):
    if count==0:continue
    ids=z[f'pose_{j}_donor_ids'];local=np.array([rev[int(i)] for i in ids])
    a=adj[np.ix_(local,local)];assert np.count_nonzero(a)//2==count
    # Only donor placements compatible with this identical TDD pose contribute.
    nq=z['donors_A'][ids,2]
    density=np.array(cKDTree(nq).query_ball_point(nq,1.,return_length=True))
    # Lexicographic maximum: weaker donor density, sum density, then smaller sample IDs.
    first=None;best=None;bestscore=None
    for i in range(len(ids)):
     js=np.flatnonzero(a[i,i+1:])+i+1
     if not len(js):continue
     if first is None:first=(i,int(js[0]))
     score1=np.minimum(density[i],density[js]);score2=density[i]+density[js]
     k=int(np.lexsort((ids[js],-score2,-score1))[0]);b=int(js[k])
     score=(int(score1[k]),int(score2[k]),-int(ids[i]),-int(ids[b]))
     if bestscore is None or score>bestscore:bestscore=score;best=(i,b)
    assert first is not None and best is not None
    poseid=len(posemeta);tri=(z['triads_A'][j]-origin)@frame
    posemeta.append(dict(pose_id=poseid,source=key,source_material=row['material'],source_site=row['site'],source_pose=j,source_support=row['support'],accepted_donors=len(ids),compatible_pairs=int(count)))
    clouds[f'pose_{poseid}_donors_A']=(z['donors_A'][ids]-origin)@frame
    clouds[f'pose_{poseid}_H_A']=(z['donor_H_A'][ids]-origin)@frame
    clouds[f'pose_{poseid}_sample_ids']=ids
    clouds[f'pose_{poseid}_triad_A']=tri
    selections=[('first',first)]
    if best!=first:selections.append(('dense',best))
    for kind,(aidx,bidx) in selections:
     donorids=ids[[aidx,bidx]]
     core=np.concatenate((z['triads_A'][j],z['donors_A'][donorids].reshape(12,3)))
     core=(core-origin)@frame
     H=(z['donor_H_A'][donorids]-origin)@frame
     assert np.allclose(np.linalg.norm(core[:23,None]-core[None,:23],axis=-1),np.linalg.norm(p.tq[:,None]-p.tq[None,:],axis=-1),atol=1e-8)
     idx=len(meta)
     meta.append(dict(id=f'N{idx+1:04d}',pose_id=poseid,source=key,source_site=row['site'],source_material=row['material'],source_pose=j,selection=kind,donor_sample_ids=donorids.tolist(),source_density_1A=[int(density[aidx]),int(density[bidx])],source_pair_count=int(count),source_support=row['support']))
     coords.append(core);hydrogens.append(H)
  np.savez_compressed(out/'core_coordinates.npz',cores_A=np.array(coords),donor_H_A=np.array(hydrogens))
  np.savez_compressed(out/'conditional_existing_donor_clouds.npz',**clouds)
  (out/'motifs.json').write_text(json.dumps(meta,indent=2));(out/'poses.json').write_text(json.dumps(posemeta,indent=2))
  report=dict(paired_TDD_poses=len(posemeta),fixed_cores=len(meta),source_sites=len(set(m['source'] for m in meta)),selection='first plus max(min per-donor 1A count), tie sum then sample IDs; same pair stored once',neighborhood='existing accepted donor placements under identical source TDD; center included; no coordinate dedup',universe='representative library, not exhaustive 70M pairs',status='LIBRARY_BUILT_TRANSFER_NOT_EVALUATED')
  assert len(posemeta)==1556
  (out/'summary.json').write_text(json.dumps(report,indent=2));print(json.dumps(report));rc=0
 finally:
  with (ROOT/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(p.Path(__file__).resolve()),script_sha256=p.sha(p.Path(__file__)),command='CPU Python fixed_core_library_v1.py',inputs='tdd_donor_all_v1',outputs=str(out),parameters=dict(radius_A=1,selection='first and dense pair per identical TDD'),exit_code=rc))+'\n')
