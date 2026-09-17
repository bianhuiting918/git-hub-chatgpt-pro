"""Identical geometric metrics for the recovered six candidates, separate report."""
from pathlib import Path
import json,hashlib
r=Path(__file__).resolve().parent;b=r/'route56_N10_20260915_v1/route6_top50_local_generation_v1'
source=r/'audit_top19_backbones_v1.py';s=source.read_text()
jobs=json.loads((b/'recovered6_RFD_v1/manifest.json').read_text())['jobs']
def load_recovered_manifest(case):
 job=next(j for j in jobs if j['name']==case);p=Path(job['input_manifest'])
 assert hashlib.sha256(p.read_bytes()).hexdigest()==job['manifest_sha256']
 j=json.loads(p.read_text())
 if 'source' in j:j=dict(j['source'],**j)
 return j
changes=[
("route6_top50_local_shapes_3seeds_v2","route6_recovered6_three_seeds_v1"),
("j=json.loads((base/'candidate_specific_shapes_v2'/cid/case/'manifest.json').read_text())","j=load_recovered_manifest(case)"),
("base/'first19_RFD3_v2'","base/'recovered6_RFD_v1'"),
("expected=57,audited=len(reports),pending=57-len(reports)","expected=18,audited=len(reports),pending=18-len(reports)"),
("'top19_backbone_audit'","'recovered6_backbone_audit'")
]
for a,c in changes:
 assert s.count(a)==1,(a,s.count(a));s=s.replace(a,c)
p=b/'recovered6_RFD_v1/auditor_provenance.json'
if not p.exists():
 with p.open('x') as f:json.dump(dict(source=str(source),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),changes=changes,input_manifests=[j['input_manifest'] for j in jobs],scope='Identical geometry gates; alternative subset/repair provenance retained in input manifests'),f,indent=2)
else:assert json.loads(p.read_text())['source_sha256']==hashlib.sha256(source.read_bytes()).hexdigest()
exec(compile(s,str(source),'exec'),{'__file__':str(Path(__file__).resolve()),'__name__':'__main__','load_recovered_manifest':load_recovered_manifest})
