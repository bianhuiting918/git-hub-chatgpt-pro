#!/usr/bin/env python3
"""Final structural/statistical artifact audit, not scientific activity validation."""
import tdd_donor_pilot_v1 as p
import numpy as np,json,datetime
ROOT=p.ROOT
def read(rel):return json.loads((ROOT/rel).read_text())
if __name__=='__main__':
 out=ROOT/'FINAL_AUDIT.json';assert not out.exists()
 a=read('input_audit_v2/summary.json');scan=read('tdd_donor_all_v1/summary.json');meta=read('fixed_core_library_v1/motifs.json');sites=read('fixed_core_transfer_v1/sites.json')
 z=np.load(ROOT/'fixed_core_transfer_v1/coverage.npz');c=z['coverage'];pc=z['parent_coverage'];q=z['cores_A'];H=z['donor_H_A']
 d=np.load(ROOT/'conditional_density_v1/counts.npz')['donor_neighbor_counts']
 assert len(scan)==262 and len(meta)==3021 and c.shape==(3021,262) and d.shape==(3021,262,2)
 assert np.isfinite(q).all() and np.isfinite(H).all() and not np.any(pc&~c)
 assert sum(c.sum(1))==7654
 source={s['material']+'_'+s['site_id']:j for j,s in enumerate(sites)}
 assert all(c[i,source[m['source']]] for i,m in enumerate(meta))
 assert all(np.array_equal(d[i,source[m['source']]],m['source_density_1A']) for i,m in enumerate(meta))
 rank=read('conditional_density_v1/rankings.json');report=read('conditional_density_v1/summary.json')
 for material in ('PA6','PA66','combined'):
  for t in (.9,.8,.7):
   key=material+'_support'+str(t);mask=np.array([(material=='combined' or s['material']==material) and s['support']>=t for s in sites])
   C=c[:,mask].sum(1);D=d[:,mask].sum(1,dtype=np.int64);weak=D.min(1);S=.6*C/max(C.max(),1)+.4*weak/max(weak.max(),1)
   rr={r['id']:r for r in rank[key]}
   for i,m in enumerate(meta):
    r=rr[m['id']]
    assert r['coverage']==C[i] and r['neighbor_counts']==D[i].tolist() and abs(r['score']-S[i])<1e-12
   assert abs(report['groups'][key]['best']['score']-S.max())<1e-12
 for mat,n in [('PA6',126),('PA66',136)]:
  assert len(read('patch_shape_chemistry_v1/'+mat+'_patches.json'))==n
  assert len(np.load(ROOT/'patch_shape_chemistry_v1'/f'{mat}_patch_face_indices.npz').files)==n
  aa=next(x for x in a if x['material']==mat)
  assert p.sha(ROOT/'inputs'/f'{mat}_v10_ge90_uncapped.pdb')==aa['input_sha256']
 result=dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),status='TECHNICAL_AND_STATISTICAL_ARTIFACT_AUDIT_PASS',input_atoms={x['material']:x['atoms'] for x in a},intact_amides={x['material']:x['amide_states']['retained'] for x in a},screened_sites=262,paired_sites={mat:sum(r['material']==mat and r['pair_records']>0 for r in scan) for mat in ['PA6','PA66']},fixed_cores=3021,coverage_edges=int(c.sum()),parent_coverage_edges=int(pc.sum()),ranking_groups=9,source_replay='all 3021 cores and their density counts PASS',unresolved_scientific_scope=['reactant geometric proxy, not validated NylC NAC or activity','independent backbone NH donors are not native NylC donor chemistry','uncapped trimmed finite material, all represented fixed cores fail parent collision control','representative 3021 cores, not exhaustive pair universe','density mesh chemistry is a proxy, not atomistic SES or calibrated hydrophobic potential'])
 out.write_text(json.dumps(result,indent=2))
 with (ROOT/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=result['time'],script=str(p.Path(__file__).resolve()),script_sha256=p.sha(p.Path(__file__)),command='CPU Python final_audit_v1.py',inputs='input, scan, library, coverage, density and patch artifacts',outputs=str(out),exit_code=0,result=result['status']))+'\n')
 print(json.dumps(result),flush=True)
