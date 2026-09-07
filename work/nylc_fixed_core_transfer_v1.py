#!/usr/bin/env python3
"""Transfer entire fixed 35-heavy-atom cores to all 262 target amide frames."""
import tdd_donor_pilot_v1 as p
import numpy as np,json,datetime,sys
ROOT=p.ROOT
sys.path.insert(0,str(ROOT.parent/'surface_previews/axis_density_sixface_all_external_trim_v6'))
from density_cut_variants import parse_pdb
from terminal_density_trim_3d import center_whole_chains_in_primary_cell
def chemistry(q,H,O):
 N=q[:,[25,31]] # 23 triad atoms, peptide donor N at local index 2
 d=np.linalg.norm(N-O,axis=-1)
 a=N-H;b=O-H
 ang=np.degrees(np.arccos(np.clip(np.sum(a*b,axis=-1)/np.linalg.norm(a,axis=-1)/np.linalg.norm(b,axis=-1),-1,1)))
 return ((d>=2.7-1e-8)&(d<=3.2+1e-8)&(ang>=140-1e-8)).all(1)
if __name__=='__main__':
 out=ROOT/'fixed_core_transfer_v1';out.mkdir(exist_ok=False);rc=1
 try:
  meta=json.loads((ROOT/'fixed_core_library_v1/motifs.json').read_text());lib=np.load(ROOT/'fixed_core_library_v1/core_coordinates.npz');q=lib['cores_A'];H=lib['donor_H_A']
  radii=np.array([p.RAD[e] for e in list(p.te)+list(p.de)*2])
  cov=np.zeros((len(q),262),bool);parentcov=cov.copy();chem=cov.copy();sites=[];col=0
  for mat,dp in [('PA6',16),('PA66',8)]:
   rec=p.load_pdb(ROOT/'inputs'/f'{mat}_v10_ge90_uncapped.pdb');xyz=np.array([a[1] for a in rec]);el=np.array([a[2] for a in rec]);lookup={l[6:11]:v for l,v,e in rec}
   pp=ROOT.parent/f'{mat.lower()}_dp{dp}_400chain/water_slab_v1/dry_export/{mat}_DP{dp}_400chain_water_equilibrated_dry.pdb'
   _,box,pa,_=parse_pdb(pp);center_whole_chains_in_primary_cell(pa,box)
   px=np.array([a['xyz'] for a in pa if a['element']!='H']);pe=np.array([a['element'] for a in pa if a['element']!='H'])
   for row in json.loads((ROOT/'carbon_exposure_v1'/f'{mat}_sites.json').read_text()):
    am=np.array([lookup[k] for k in row['serials']]);c=am[0];ex=p.norm(am[1]-c);ey=am[2]-c;ey=p.norm(ey-ex*np.dot(ex,ey));F=np.column_stack((ex,ey,np.cross(ex,ey)))
    O=(am[1]-c)@F;ok=chemistry(q,H,O);chem[:,col]=ok;ids=np.flatnonzero(ok)
    if len(ids):
     env=p.trees((xyz-c)@F,el);cov[ids,col]=p.clear(q[ids],radii,env)>=-1e-8
     accepted=np.flatnonzero(cov[:,col])
     if len(accepted):parentcov[accepted,col]=p.clear(q[accepted],radii,p.trees((px-c)@F,pe))>=-1e-8
    sites.append(dict(material=mat,site_id=row['site_id'],support=row['inward_support'],side=row['nearest_side'],origin_A=c.tolist(),frame=F.tolist()))
    col+=1
    if col%25==0:print(json.dumps(dict(targets_done=col,total=262)),flush=True)
  assert col==262
  source_lookup={s['material']+'_'+s['site_id']:i for i,s in enumerate(sites)}
  assert all(cov[i,source_lookup[m['source']]] for i,m in enumerate(meta)),'Source replay failure'
  assert not np.any(parentcov&~cov)
  np.savez_compressed(out/'coverage.npz',coverage=cov,parent_coverage=parentcov,chemical_pass=chem,cores_A=q,donor_H_A=H)
  (out/'sites.json').write_text(json.dumps(sites,indent=2))
  ranking=[]
  m6=np.array([s['material']=='PA6' for s in sites]);m66=~m6
  for i,m in enumerate(meta):
   ranking.append(dict(m,coverage_PA6=int(cov[i,m6].sum()),coverage_PA66=int(cov[i,m66].sum()),coverage_total=int(cov[i].sum()),parent_coverage_PA6=int(parentcov[i,m6].sum()),parent_coverage_PA66=int(parentcov[i,m66].sum())))
  ranking.sort(key=lambda r:(-r['coverage_total'],r['id']))
  (out/'rankings.json').write_text(json.dumps(ranking,indent=2))
  report=dict(fixed_cores=len(meta),target_PA6=126,target_PA66=136,coverage_edges=int(cov.sum()),best=ranking[:5],status='FIXED_GEOMETRY_TRANSFER_COMPLETE_CONDITIONAL_DENSITY_PENDING',parent_control='finite centered original dry parent, no periodic images; not activity or ingress')
  (out/'summary.json').write_text(json.dumps(report,indent=2));print(json.dumps(report),flush=True);rc=0
 finally:
  with (ROOT/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(p.Path(__file__).resolve()),script_sha256=p.sha(p.Path(__file__)),command='CPU Python fixed_core_transfer_v1.py',inputs='fixed_core_library_v1 and all 262 existing target frames',outputs=str(out),parameters=dict(core_atoms=35,internal_geometry='unchanged',overlap_A=.4),exit_code=rc))+'\n')
