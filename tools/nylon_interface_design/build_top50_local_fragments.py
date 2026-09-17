"""CPU local-fragment construction, NOT RFD generation.
Reuse checked optimizer and thresholds; test three explicit N10 five-residue templates.
All50 donor candidates retained, including cysteine with its own atom topology.
"""
import os
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
from pathlib import Path
import json,hashlib,datetime,numpy as np
from scipy.spatial import cKDTree
r=Path(__file__).resolve().parent;os.chdir(r)
base=r/'route56_N10_20260915_v1/route6_top50_local_generation_v1'
manifest=json.loads((base/'manifest.json').read_text())
out=base/'local_fragments_three_templates_v1';out.mkdir(exist_ok=False)
source=(r/'prepare_supported_arg_diagnostic_v1.py').read_text()
prefix=source.split("cov=np.load(p/'coverage.npz')")[0]
assert "fixed=[8,9,10,11]" in prefix and "max_nfev=150" in prefix
(out/'optimizer_source_snapshot.txt').write_text(prefix)
sites=json.loads((r/'route56_N10_20260915_v1/N10_fixed_minimal_motif_material_v1/sites.json').read_text())
rad={'C':1.7,'N':1.55,'O':1.52,'S':1.8};trees={}
for mat in ['PA6','PA66']:
 z=np.load(r.parent/'interface_surface_robustness_20260908_v1/inputs_v2'/(mat+'_parent_atoms.npz'))
 xyz=z['xyz_A'].copy();box=z['box_A'];xyz[:,:2]%=box[:2]
 trees[mat]=(box,{str(el):cKDTree(xyz[z['elements']==el],boxsize=[*box[:2],0]) for el in np.unique(z['elements']) if el!='H'})
rows=[]
for v in manifest['rows']:
 for start in [24,68,170]:
  tag=v['id']+'_N10_'+str(start)+'_'+str(start+4)
  row=dict(id=v['id'],microstate=v['microstate'],template_start=start,status='NOT_EVALUATED',total_old_covered=v['total_covered'])
  try:
   code=prefix.replace("=='SC4559268'","=="+repr(v['id'])).replace("170<=int(l[22:26])<=174",str(start)+"<=int(l[22:26])<="+str(start+4)).replace("range(170,175)",f"range({start},{start+5})")
   if v['microstate']=='CYS_4':
    code=code.replace("'ARG_1.npz'","'CYS_4.npz'").replace("[27,28,29,30,31,39,40,41,42]","[27,28,29,30,31,34,35,36,37]").replace("q[32:39]","q[32:34]").replace("z['r'][32:39]","z['r'][32:34]").replace("np.full((27,27),1000)","np.full((22,22),1000)").replace("range(27):","range(22):").replace("[(9,20),(20,21),(21,22),(22,23),(23,24),(24,25),(24,26)]","[(9,20),(20,21)]")
   ns={};exec(compile(code,tag,'exec'),ns)
   P=ns['P'];q=ns['q'];fragment=ns['fragment'];rr=ns['rr']
   cn=[float(np.linalg.norm(P[4*j+2]-P[4*j+4])) for j in range(4)]
   omega=[180-abs(np.degrees(ns['tors'](P,*b))) for b in ns['planes'][:4]]
   passed=bool(ns['sol'].success and ns['m'].min()>=-1e-8 and ns['cm'].min()>=-1e-8 and max(omega)<=30 and all(1.2<=x<=1.5 for x in cn) and np.array_equal(P[ns['fixed']],ns['before'][ns['fixed']]))
   clear=[]
   if passed:
    for sj in v['site_indices']:
     site=sites[sj];box,tt=trees[site['material']]
     w=fragment@np.array(site['frame']).T+site['center'];w[:,:2]%=box[:2];margin=np.full(len(w),np.inf)
     for el,tr in tt.items():margin=np.minimum(margin,tr.query(w)[0]-rr-rad[el]+.4)
     if margin.min()>=-1e-8:clear.append(sj)
   dest=out/(tag+'.npz')
   np.savez_compressed(dest,core_original=q,fragment=fragment,backbone=P,radii=rr,clear_site_indices=np.array(clear,dtype=int))
   row.update(status='LOCAL_FRAGMENT_GEOMETRY_PASS_NOT_RFD' if passed else 'LOCAL_FRAGMENT_GEOMETRY_FAIL',optimizer_success=bool(ns['sol'].success),self_margin_A=float(ns['m'].min()),core_margin_A=float(ns['cm'].min()),CN_A=cn,omega_max_deviation_deg=float(max(omega)),fixed_backbone_delta_A=float(np.linalg.norm(P[ns['fixed']]-ns['before'][ns['fixed']],axis=1).max()),material_counts={mat:sum(sites[j]['material']==mat for j in clear) for mat in ['PA6','PA66']},clear_site_indices=clear,output=str(dest),sha256=hashlib.sha256(dest.read_bytes()).hexdigest())
  except Exception as e:row.update(status='NOT_EVALUATED_EXCEPTION',error=repr(e))
  rows.append(row)
  with (out/'results.jsonl').open('a') as h:h.write(json.dumps(row)+'\n')
 print(json.dumps(dict(candidate=v['id'],finished_templates=3,local_geometry_pass=sum(x['id']==v['id'] and x['status']=='LOCAL_FRAGMENT_GEOMETRY_PASS_NOT_RFD' for x in rows))),flush=True)
report=dict(status='LOCAL_CONSTRUCTION_AUDIT_NOT_RFD',candidates=50,templates=[24,68,170],rows=rows,source_optimizer_sha256=hashlib.sha256(prefix.encode()).hexdigest(),limitations=['Template grafting and optimization, not diffusion-generated backbones','Only previously eligible material sites replayed','Support sidechains omitted; full fold and sequence not evaluated','No catalytic or activity PASS implied'])
(out/'summary.json').write_text(json.dumps(report,indent=2))
(out/'RUNBOOK.md').write_text('Run saved build_top50_local_fragments_v1.py using existing project CPU Python. Same source optimizer,150 evaluations and geometric gates as prior single-candidate preparation. Three explicit template windows24-28,68-72,170-174. Cys has separate atom indices/topology. All150 results including exceptions retained. RFD remains not submitted.\n')
with (r/'RUN_LOG.jsonl').open('a') as h:h.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__)),output=str(out),exit_code=0,records=len(rows),rfd_submitted=0))+'\n')
print('COMPLETE',len(rows))
