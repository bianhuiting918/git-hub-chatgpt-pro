"""Convert eight binned dual-donor fits to same-protocol closed CA shapes; CPU only."""
import os
for k in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import hashlib,sys
r=Path(__file__).resolve().parent;p=r/'prepare_route6_complete_closed_CA_nodes_v2.py'
assert hashlib.sha256(p.read_bytes()).hexdigest()=='94b02857914140fc7270da148bc6ab365f4e2b0ef9c9d5604aa70d4ca4e43402'
s=p.read_text()
changes=[
("route6_complete_closed_CA_nodes_v2","route6_binned_lowcov_closed_CA_v1"),
(" for mat,cid in [('PA6','N2627'),('PA66','N0147')]:"," jobs=json.loads((R/'route56_N10_20260915_v1/route6_two_donor_binned_lowcov_v1/manifest.json').read_text())['jobs']\n for job in jobs:\n  mat=job['material'];cid=job['core'];bn=job['exposure_bin']"),
("fit=P/f'route6_complete_{mat}_{cid}_whole_site_v1'","fit=P/job['output']"),
("for percent in [90,70,50]:","for percent in [10,20,30]:"),
("rec=dict(material=mat,core=cid,","rec=dict(exposure_bin=bn,multi_site=(len(ix)>=2),all_exposure_group_denominator=job['all_bin_n'],material=mat,core=cid,"),
("f'{mat}_{cid}_","f'{mat}_{cid}_bin{bn}_")]
for a,b in changes:
 assert a in s,a
 s=s.replace(a,b)

s=s.replace('route6_binned_lowcov_closed_CA_v1','route6_binned_lowcov_closed_CA_v2')
bad=';assert np.all(cnt==2)'
assert s.count(bad)==1
s=s.replace(bad,"\n    if not np.all(cnt==2):\n     rec.update(status='NONMANIFOLD_MESH_REJECTED',bad_edge_count=int(np.sum(cnt!=2)));reports.append(rec);print(json.dumps(rec),flush=True);continue")
start=s.index('  for j,si in enumerate(ids):')
end=s.index('  fields=np.array(fields);',start)
body=s[start:end]
cache="""  cached=R/'route56_N10_20260915_v1/route6_binned_lowcov_closed_CA_v1'/f'{mat}_{cid}_bin{bn}_fields.npz'
  if cached.exists():
   cz=np.load(cached)
   assert np.array_equal(cz['site_indices'],ids)
   assert np.array_equal(cz['core_q_A'],q-center)
   assert np.array_equal(cz['chemical_to_fit_rotation'],rot)
   assert all(np.array_equal(cz['axis_'+k+'_A'],a) for k,a in zip('xyz',axes))
   fields=cz['clearance_minus_CA_allowance']
   print('REUSED_FIELD',str(cached),hashlib.sha256(cached.read_bytes()).hexdigest(),flush=True)
  else:
"""
s=s[:start]+cache+'\n'.join(' '+line if line else '' for line in body.split('\n'))+s[end:]

compile(s,str(p),'exec')
if '--check' in sys.argv:
 print('COMPILED: same nodes/grid/radii/components protocol, 8 original-bin jobs, 10/20/30 tiers; NO GPU')
 raise SystemExit(0)
exec(compile(s,str(p),'exec'),{'__file__':str(Path(__file__).resolve()),'__name__':'__main__'})
