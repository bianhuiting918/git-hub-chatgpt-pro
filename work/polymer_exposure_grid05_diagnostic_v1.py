import json,time
from pathlib import Path
import numpy as np
import polymer_exposure_v2 as m
r=Path(__file__).resolve().parent
out=r/'polymer_exposure_grid05_diagnostic_v1';out.mkdir()
prior=json.load(open(r/'polymer_exposure_audit_v1/summary.json'));result=[]
for rec in prior['materials']:
 mat=rec['material'];a=np.load(r/'inputs_v2'/f'{mat}_parent_atoms.npz')
 env=m.Env(a['xyz_A'],a['elements'],a['box_A'])
 m.emit('grid05_start',material=mat)
 g=m.build_grid(env,a['masses_Da'],.5)
 original={x['site_id']:x for x in json.load(open(r/'polymer_exposure_full_v2'/mat/'sites.json'))}
 pick=sorted(rec['refinement'],key=lambda x:-x['grid_change'])[:4]
 data=[]
 for item in pick:
  row=original[item['site_id']];_,loc,ext=m.exposure(np.array(row['center']),env,g,8192)
  data.append(dict(site_id=item['site_id'],fraction_grid1_samples4096=row['external_shell_fraction'],fraction_grid075_samples8192=item['grid075_samples8192_exposure'],fraction_grid05_samples8192=float(ext.mean()),local_fraction_samples8192=float(loc.mean()),change075_to05=abs(float(ext.mean())-item['grid075_samples8192_exposure'])))
 v=dict(material=mat,n=len(data),sites=data);result.append(v)
 (out/(mat+'.json')).write_text(json.dumps(v,indent=2));m.emit('grid05_complete',**v)
doc=dict(status='TARGETED_DIAGNOSTIC_ONLY_NOT_FULL_ATLAS_REPLACEMENT',grid_A=.5,samples=8192,materials=result)
(out/'summary.json').write_text(json.dumps(doc,indent=2))
with (r/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(event='grid05_diagnostic_complete',output=str(out),exit_status=0))+'\n')
