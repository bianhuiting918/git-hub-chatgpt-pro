#!/usr/bin/env python3
"""Independent direct summation of reported 1A neighbor-by-site counts."""
from pathlib import Path
import json,datetime
import numpy as np
R=Path(__file__).resolve().parent;S=R/'surface_scan_v1';F=R/'fixed_coupled_transfer_v1'
report=json.loads((F/'summary.json').read_text());triads=json.loads((F/'triads.json').read_text())
import importlib.util
spec=importlib.util.spec_from_file_location('s',R/'surface_scan_v1.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
rows,_,_=m.get_inputs();tcover=np.load(F/'triad_coverage.npz')['coverage']
keys=report['microstates'];cloud={k:dict(np.load(S/'clouds'/(k+'.npz'))) for k in keys}
polys=[];casepaths=[]
for item in rows:
 key=item['material']+'_'+item['row']['site_id'];path=S/'cases'/(key+'.npz');casepaths.append(path)
 z=np.load(path);polys.append({k:set(z[k+'_polymer_ids'].tolist()) for k in keys})
checks=[]
for co in report['cohorts']:
 targetids=[i for i,it in enumerate(rows) if (it['material']==co['material'] or co['material']=='NYLON_COMBINED' and it['material'] in ['PA6','PA66']) and float(it['row'].get('support_fraction',it['row'].get('inward_support')))>=co['min_support']]
 for rr in co['per_microstate']:
  if 'best_balanced' not in rr:continue
  w=rr['best_balanced'];k=w['microstate'];ti=w['triad_index'];tr=triads[ti];c=cloud[k];sid=w['donor_sample_id']
  z=np.load(casepaths[tr['source_index']]);ids=z[f"{k}_{tr['template']}_{tr['source_pose']}_ids"]
  same=(c['pair_indices'][ids]==c['pair_indices'][sid]).all(1);actual=[]
  for d in range(2):
   neighbors=set(ids[same&(np.linalg.norm(c['donor_A'][ids,d]-c['donor_A'][sid,d],axis=1)<=1.)].tolist())
   value=sum(len(neighbors.intersection(polys[i][k])) for i in targetids if tcover[i,ti])
   actual.append(value)
  assert actual==w['donor_neighbor_site_incidences'],(co['name'],w['candidate_id'],actual,w['donor_neighbor_site_incidences'])
  checks.append(dict(cohort=co['name'],candidate_id=w['candidate_id'],counts=actual))
out=dict(status='DIRECT_NEIGHBOR_SITE_SUM_PASS',checks=len(checks),results=checks,source_sha256=m.sha(Path(__file__)))
(F/'neighbor_count_audit.json').write_text(json.dumps(out,indent=2))
with (R/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__)),source_sha256=m.sha(Path(__file__)),outputs=str(F/'neighbor_count_audit.json'),exit_code=0,checks=len(checks)))+'\n')
print(out['status'],len(checks),flush=True)
