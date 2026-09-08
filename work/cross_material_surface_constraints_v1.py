#!/usr/bin/env python3
"""Material-balanced common substrate-frame constraints; no averaged enzyme core."""
from pathlib import Path
import json,datetime,hashlib
import numpy as np
R=Path(__file__).resolve().parent;P=R/'panel_surface_scan_v1';O=R/'cross_material_constraints_v1'
O.mkdir(exist_ok=False)
sites=json.loads((P/'sites.json').read_text());cov=np.load(P/'coverage.npz')['coverage'];summary=[]
for label,surfaces in [('original',['PET_original','PA6_original','PA66_original']),('trimmed_nylon',['PA6_trimmed','PA66_trimmed'])]:
 for t in [.7,.8,.9]:
  groups={sf:np.array([i for i,s in enumerate(sites) if s['surface']==sf and s['support']>=t and cov[:,i].any()],int) for sf in surfaces}
  if any(len(ids)==0 for ids in groups.values()):
   summary.append(dict(object=label,min_support=t,status='NOT_EVALUATED_EMPTY_MATERIAL'));continue
  xx={sf:np.array([np.load(P/'fields'/(str(i)+'.npz'))['polymer_vdw_clearance_A'] for i in ids]) for sf,ids in groups.items()}
  probs=np.array([(x<1.3).mean(0) for x in xx.values()]);p=probs.mean(0)
  minima=np.array([x.min(0) for x in xx.values()]);minimum=minima.min(0);skin=(minimum>=1.3)&(minimum<=4.3)
  disagreement=probs.max(0)-probs.min(0)
  tests=[]
  for sf,x in xx.items():
   train_min=np.min([v.min(0) for k,v in xx.items() if k!=sf],axis=0);allowed=(train_min>=1.3)&(train_min<=4.3)
   tests.append(dict(heldout_material=sf,allowed_voxels=int(allowed.sum()),status='EVALUATED_DESCRIPTIVE_MATERIAL_HOLDOUT' if allowed.sum()>=10 else 'NOT_EVALUATED_LT10_ALLOWED_VOXELS',heldout_clash_fraction_mean=float((x[:,allowed]<1.3).mean()) if allowed.sum()>=10 else None))
  tag=label+'_ge'+str(round(t*100))
  np.savez_compressed(O/(tag+'.npz'),material_balanced_occupancy_probability=p.astype(np.float32),material_occupancy_range=disagreement.astype(np.float32),minimum_clearance_A=minimum.astype(np.float32),strict_common_near_interface_mask=skin,**{sf+'_site_indices':ids for sf,ids in groups.items()})
  summary.append(dict(object=label,min_support=t,sites_by_material={sf:len(ids) for sf,ids in groups.items()},strict_common_near_interface_voxels=int(skin.sum()),material_holdouts=tests,status='COMMON_SUBSTRATE_FRAME_CONSTRAINTS_ONLY',scope='carbonyl C/O/leaving-group frame; equal material weights, unique sites; NOT an averaged enzyme core, NOT proof one protein realizes this envelope; alternative core poses not superposed into a protein-coordinate frame; no static charge inversion'))
(O/'summary.json').write_text(json.dumps(summary,indent=2))
with (R/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__)),sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),outputs=str(O),exit_code=0))+'\n')
print(json.dumps(summary),flush=True)
