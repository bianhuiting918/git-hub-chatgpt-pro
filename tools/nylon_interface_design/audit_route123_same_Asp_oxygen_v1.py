"""Re-audit route1/2 using one consistent Asp1 oxygen; preserve legacy results."""
from pathlib import Path
import json,hashlib,datetime,collections
import numpy as np
R=Path(__file__).resolve().parent
src=R/'route123_canonical_joint_audit_20260915T232917Z.json'
old=json.loads(src.read_text())
assert len(old['records'])==501
records=[]
for r in old['records']:
 p=Path(r['prediction']);raw=p.read_bytes()
 assert hashlib.sha256(raw).hexdigest()==r['prediction_sha256']
 a={};chains=set()
 for l in raw.decode().splitlines():
  if not l.startswith('ATOM'):continue
  ri=int(l[22:26]);at=l[12:16].strip()
  if ri not in [90,129,131]:continue
  key=(ri,at)
  assert key not in a,('Duplicate core atom',p,key)
  a[key]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
  chains.add(l[21])
 assert len(chains)==1
 o=min(['OD1','OD2'],key=lambda k:np.linalg.norm(a[90,'N']-a[129,k]))
 v=a[129,o]
 d2=min(['OD1','OD2'],key=lambda k:np.linalg.norm(v-a[131,k]))
 distances=[float(np.linalg.norm(a[90,'N']-v)),float(np.linalg.norm(a[90,'N']-a[90,'OG1'])),float(np.linalg.norm(v-a[131,d2]))]
 intervals=[(2.5,3.5),(2.5,3.3),(3.3,4.7)]
 violations=[max(lo-d,0,d-hi) for d,(lo,hi) in zip(distances,intervals)]
 contacts=[c for c in r.get('N90_H_D129_contacts',[]) if c['acceptor']==o]
 angle=max((c['DHO_deg'] for c in contacts),default=None)
 direction=any(2.5<=c['DO']<=3.5 and c['HO']<=2.6 and c['DHO_deg']>=130 for c in contacts)
 records.append(dict(id=r['id'],route=r['route'],prediction=str(p),prediction_sha256=r['prediction_sha256'],
                     selected_Asp1_oxygen=o,selected_Asp2_oxygen=d2,distances_A=distances,
                     distance_violation_A=violations,total_distance_violation_A=sum(violations),
                     distance_gate=not any(violations),selected_oxygen_best_NH_angle_deg=angle,
                     direction_gate=direction if contacts else None,joint_distance_direction=not any(violations) and direction,
                     mean_CA_pLDDT=r['mean_CA_pLDDT'],legacy_distances_A=r['distances_A'],legacy_distance_gate=r['distance_gate'],
                     changed_DD_distance=abs(distances[2]-r['distances_A'][2])>1e-4,
                     stage='uncleaved_AF2_precursor',material_reevaluation='NOT_RECOMPUTED_NO_COORDINATE_CHANGE',
                     donor_reevaluation='NOT_RECOMPUTED_NO_COORDINATE_CHANGE'))
summaries=[]
for route in ['route1','route2','reference']:
 rows=[r for r in records if r['route']==route]
 ranked=sorted(rows,key=lambda r:(r['total_distance_violation_A'],-r['mean_CA_pLDDT'],r['id']))
 summaries.append(dict(route=route,n=len(rows),distance_pass=sum(r['distance_gate'] for r in rows),
                       joint_pass=sum(r['joint_distance_direction'] for r in rows),
                       changed_DD_distance=sum(r['changed_DD_distance'] for r in rows),
                       distance_pass_ids=[r['id'] for r in ranked if r['distance_gate']],ranking=[r['id'] for r in ranked]))
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
out=R/('route123_same_Asp_oxygen_audit_'+stamp+'.json')
result=dict(status='CORRECTED_DEFINITION_PRECURSOR_AUDIT_NOT_ACTIVITY',
            source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),
            definition='Asp1 oxygen nearest Thr90 N is reused for its distance to Asp2; choose nearest Asp2 oxygen. All numerical intervals unchanged.',
            records=records,summaries=summaries,
            limitations=['Legacy JSON and structures unchanged','Direction reuses legacy added-H contacts, not newly optimized hydrogens',
                         'Mature comparisons require separate stage-specific records','Distance pass alone is not catalytic success'])
with out.open('x') as f:json.dump(result,f,indent=2)
with (R/'RUN_LOG.jsonl').open('a') as f:
 f.write(json.dumps(dict(time=stamp,script=str(Path(__file__)),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                        output=str(out),exit_code=0,status=result['status']))+'\n')
print(json.dumps(dict(output=str(out),summaries=[{k:v for k,v in r.items() if k!='ranking'} for r in summaries])))
