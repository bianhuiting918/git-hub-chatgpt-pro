"""Reuse same zero-level precision diagnostic for SC5011126; CPU only."""
from pathlib import Path
import json,hashlib
r=Path(__file__).resolve().parent
source=r/'repair_zero_level_mesh_float64_v2.py';s=source.read_text()
rows=[json.loads(x) for x in (r/'route56_N10_20260915_v1/route6_top50_local_generation_v1/local_fragments_three_templates_v1/results.jsonl').read_text().splitlines()]
row=next(x for x in rows if x['id']=='SC5011126' and x['status']=='LOCAL_FRAGMENT_GEOMETRY_PASS_NOT_RFD')
assert row['template_start']==68
for a,b in [
("cid='SC4559268';name='PA6_SC4559268_coverage10_variant1'","cid='SC5011126';name='PA6_SC5011126_coverage10_variant1'"),
("SC4559268_N10_170_174.npz","SC5011126_N10_68_72.npz")]:
 assert s.count(a)==1;s=s.replace(a,b)
exec(compile(s,str(source),'exec'),{'__file__':str(Path(__file__).resolve()),'__name__':'__main__'})
