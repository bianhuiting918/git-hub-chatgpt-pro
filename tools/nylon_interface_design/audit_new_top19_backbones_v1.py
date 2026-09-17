"""Incremental same-protocol audit; reuse only matching PDB SHA256 records."""
from pathlib import Path
import json,hashlib
r=Path(__file__).resolve().parent;b=r/'route56_N10_20260915_v1/route6_top50_local_generation_v1'
source=r/'audit_top19_backbones_v1.py'
assert hashlib.sha256(source.read_bytes()).hexdigest()=='95b560c8509d978722b400ea6e36f3bd21f853144a51b8900470c969dab86c23'
jobs=json.loads((b/'recovered6_RFD_v1/manifest.json').read_text())['jobs']
def load_recovered_manifest(case):
 job=next(j for j in jobs if j['name']==case);p=Path(job['input_manifest'])
 assert hashlib.sha256(p.read_bytes()).hexdigest()==job['manifest_sha256']
 j=json.loads(p.read_text())
 if 'source' in j:j=dict(j['source'],**j)
 return j
for folder,gfolder,total in [('first19_RFD3_v2','route6_top50_local_shapes_3seeds_v2',39),('recovered6_RFD_v1','route6_recovered6_three_seeds_v1',18)]:
 files=sorted((r.parent/'n9_core_donor_design_20260910_gpu_outputs_v1'/gfolder).glob('*/design_*.pdb'))
 previous={}
 for p in sorted((b/folder).glob('backbone_audit_*.json')):
  for row in json.loads(p.read_text())['reports']:previous[row['path']]=row
 cached=[];pending=set()
 for p in files:
  row=previous.get(str(p))
  if row and row['sha256']==hashlib.sha256(p.read_bytes()).hexdigest():cached.append(row)
  else:pending.add(str(p))
 print(json.dumps(dict(batch=folder,cached=len(cached),new=len(pending),planned=total)),flush=True)
 if not pending:continue
 s=source.read_text()
 old="for path in sorted(g.glob('*/design_*.pdb')):"
 assert s.count(old)==1
 s=s.replace(old,"for path in sorted(g.glob('*/design_*.pdb')):\n if str(path) not in pending_paths:continue")
 s=s.replace('reports=[]','reports=list(cached_reports)')
 s=s.replace('expected=57,audited=len(reports),pending=57-len(reports)',f'expected={total},audited=len(reports),pending={total}-len(reports)')
 if folder=='recovered6_RFD_v1':
  for old,new in [('route6_top50_local_shapes_3seeds_v2',gfolder),("j=json.loads((base/'candidate_specific_shapes_v2'/cid/case/'manifest.json').read_text())","j=load_recovered_manifest(case)"),("base/'first19_RFD3_v2'","base/'recovered6_RFD_v1'"),("'top19_backbone_audit'","'recovered6_backbone_audit'")]:
   assert s.count(old)==1;s=s.replace(old,new)
 exec(compile(s,str(source),'exec'),{'__file__':str(Path(__file__).resolve()),'__name__':'__main__','pending_paths':pending,'cached_reports':cached,'load_recovered_manifest':load_recovered_manifest})
