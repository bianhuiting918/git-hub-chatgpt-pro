"""Reuse the same audited backbone metrics for T100/T200, separate report."""
from pathlib import Path
import hashlib,json
r=Path(__file__).resolve().parent
source=r/'audit_top19_backbones_v1.py'
s=source.read_text()
changes=[
("route6_top50_local_shapes_3seeds_v2","route6_step_count_diagnostic_v1"),
("case=path.parent.name;cid=case.split('_')[1]","case=path.parent.name;cid='SC5150928'"),
("base/'candidate_specific_shapes_v2'/cid/case/'manifest.json'","base/'candidate_specific_shapes_v2'/cid/'PA6_SC5150928_coverage10_variant1'/'manifest.json'"),
("base/'first19_RFD3_v2'","base/'step_count_diagnostic_v1'"),
("expected=57,audited=len(reports),pending=57-len(reports)","expected=2,audited=len(reports),pending=2-len(reports)"),
("'top19_backbone_audit'","'step_count_backbone_audit'")
]
for a,b in changes:
 assert s.count(a)==1,(a,s.count(a))
 s=s.replace(a,b)
p=r/'route56_N10_20260915_v1/route6_top50_local_generation_v1/step_count_diagnostic_v1/auditor_provenance.json'
if not p.exists():
 with p.open('x') as f:json.dump(dict(source=str(source),sha256=hashlib.sha256(source.read_bytes()).hexdigest(),changes=changes,scope='Only input paths and expected count change; all metrics and thresholds identical'),f,indent=2)
else:assert json.loads(p.read_text())['sha256']==hashlib.sha256(source.read_bytes()).hexdigest()
exec(compile(s,str(source),'exec'),{'__file__':str(Path(__file__).resolve()),'__name__':'__main__'})
