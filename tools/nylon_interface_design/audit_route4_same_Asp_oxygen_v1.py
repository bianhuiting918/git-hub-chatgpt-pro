"""Versioned same-Asp-oxygen correction for route4 and matched references."""
from pathlib import Path
import json,hashlib,datetime,collections
import numpy as np
R=Path(__file__).resolve().parent
src=R/'route4_expanded_matched_snapshot_20260916T011223Z.json'
old=json.loads(src.read_text());assert len(old['records'])==735
records=[]
for r in old['records']:
 p=Path(r['path']);raw=p.read_bytes()
 assert hashlib.sha256(raw).hexdigest()==r['source_sha256']
 t,d1,d2=r['core'];a={}
 for l in raw.decode().splitlines():
  if not l.startswith('ATOM') or int(l[22:26]) not in [t,d1,d2]:continue
  k=(int(l[22:26]),l[12:16].strip())
  assert k not in a
  a[k]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
 o=min(['OD1','OD2'],key=lambda o:np.linalg.norm(a[t,'N']-a[d1,o]))
 v=a[d1,o];o2=min(['OD1','OD2'],key=lambda z:np.linalg.norm(v-a[d2,z]))
 ds=[float(np.linalg.norm(a[t,'N']-v)),float(np.linalg.norm(a[t,'N']-a[t,'OG1'])),float(np.linalg.norm(v-a[d2,o2]))]
 gap=[max(lo-d,0,d-hi) for d,(lo,hi) in zip(ds,[(2.5,3.5),(2.5,3.3),(3.3,4.7)])]
 contacts=[c for c in r['N90_or_258_H_D1_contacts'] if c['acceptor']==o]
 direct=any(2.5<=c['DO']<=3.5 and c['HO']<=2.6 and c['DHO_deg']>=130 for c in contacts)
 records.append(dict(name=r['name'],family=r['family'],stage=r['stage'],seed=r['seed'],path=str(p),source_sha256=r['source_sha256'],core=r['core'],
                     selected_Asp1_oxygen=o,selected_Asp2_oxygen=o2,distances_A=ds,distance_violation_A=gap,distance_gap_sum_A=sum(gap),
                     distance_gate=not any(gap),direction_diagnostic=direct,distance_and_direction_diagnostic=not any(gap) and direct,
                     selected_oxygen_contacts=contacts,legacy_distances_A=r['distances_A'],legacy_distance_gate=r['distance_gate'],
                     technical_status=r.get('technical_status'),donor_material_joint='NOT_RECOMPUTED_NO_COORDINATE_CHANGE'))
summary=[]
for family,stage in sorted({(r['family'],r['stage']) for r in records}):
 rr=[r for r in records if r['family']==family and r['stage']==stage]
 summary.append(dict(family=family,stage=stage,n=len(rr),distance_pass=sum(r['distance_gate'] for r in rr),
                     joint_pass=sum(r['distance_and_direction_diagnostic'] for r in rr)))
rank=sorted([r for r in records if r['family']=='route4_expanded' and r['stage']=='after'],
            key=lambda r:(r['distance_gap_sum_A'],-max(c['DHO_deg'] for c in r['selected_oxygen_contacts']),r['name'],r['seed']))
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
out=R/('route4_same_Asp_oxygen_audit_'+stamp+'.json')
with out.open('x') as f:
 json.dump(dict(status='CORRECTED_STATIC_GEOMETRY_NOT_ACTIVITY',source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),
                definition='Same Asp1 oxygen nearest Thr N reused in Asp1-Asp2 distance; nearest Asp2 oxygen; unchanged intervals.',
                records=records,summary=summary,closest_after=[dict(name=r['name'],seed=r['seed'],gap_A=r['distance_gap_sum_A']) for r in rank[:10]],
                limitations=['Legacy files preserved','Contacts reuse stage-specific legacy hydrogen geometry','Not full donor, steric, folding or activity validation',
                             'N10 two mature seeds have identical coordinates and are not independent conformations']),f,indent=2)
with (R/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=stamp,script=str(Path(__file__)),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),output=str(out),exit_code=0))+'\n')
print(json.dumps(dict(output=str(out),summary=summary)))
