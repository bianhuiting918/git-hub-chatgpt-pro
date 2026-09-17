"""Strict geometry snapshot of the five route2 matched-maturation diagnostics."""
from pathlib import Path
import json,hashlib,datetime,math
import numpy as np
R=Path(__file__).resolve().parent
B=R/'matched_cleavage_route2_sameO_five_v1'
manifest=json.loads((B/'manifest.json').read_text())
def angle(x,y,z):
 u=x-y;v=z-y
 return float(np.degrees(np.arccos(np.clip(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v),-1,1))))
def dihedral(p):
 b0=-(p[1]-p[0]);b1=p[2]-p[1];b2=p[3]-p[2];b1/=np.linalg.norm(b1)
 v=b0-np.dot(b0,b1)*b1;w=b2-np.dot(b2,b1)*b1
 return float(np.degrees(np.arctan2(np.dot(np.cross(b1,v),w),np.dot(v,w))))
rows=[];pending=[]
for name in manifest['ids']:
 for seed in manifest['seeds']:
  ap=B/'outputs'/name/('seed_'+str(seed))/'audit.json'
  if not ap.exists():pending.append([name,seed]);continue
  audit=json.loads(ap.read_text())
  if 'pdb_path' not in audit:
   rows.append(dict(name=name,seed=seed,status='NOT_EVALUATED_NO_PDB',technical_status=audit['status']));continue
  src=Path(audit['pdb_path']);raw=src.read_bytes()
  assert hashlib.sha256(raw).hexdigest()==audit['pdb_sha256']
  a={};order=[];elements={}
  for l in raw.decode().splitlines():
   if not l.startswith('ATOM'):continue
   k=(l[21],int(l[22:26]));at=l[12:16].strip()
   if k not in a:a[k]={};order.append(k)
   assert at not in a[k]
   a[k][at]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
   elements[(k,at)]=l[76:78].strip()
  assert len({k[1] for k in order})==len(order)
  res={k[1]:v for k,v in a.items()}
  N=res[90]['N'];OG=res[90]['OG1'];o=min(['OD1','OD2'],key=lambda k:np.linalg.norm(N-res[129][k]));O=res[129][o]
  o2=min(['OD1','OD2'],key=lambda k:np.linalg.norm(O-res[131][k]))
  ds=[float(np.linalg.norm(N-O)),float(np.linalg.norm(N-OG)),float(np.linalg.norm(O-res[131][o2]))]
  gaps=[max(lo-d,0,d-hi) for d,(lo,hi) in zip(ds,[(2.5,3.5),(2.5,3.3),(3.3,4.7)])]
  contacts=[]
  for key,h in res[90].items():
   if key not in ['1H','2H','3H','H','H1','H2','H3']:continue
   contacts.append(dict(H=key,NO_A=ds[0],HO_A=float(np.linalg.norm(h-O)),NHO_deg=angle(N,h,O)))
  direction=any(2.5<=c['NO_A']<=3.5 and c['HO_A']<=2.6 and c['NHO_deg']>=130 for c in contacts)
  bonds=[]
  for i,j in zip(order,order[1:]):
   if i[0]!=j[0] or j[1]!=i[1]+1:continue
   d=float(np.linalg.norm(a[i]['C']-a[j]['N']))
   omega=dihedral(np.array([a[i]['CA'],a[i]['C'],a[j]['N'],a[j]['CA']]))
   bonds.append(dict(i=list(i),j=list(j),CN_A=d,omega_deg=omega,trans_deviation_deg=abs(180-abs(omega))))
  bad=[b for b in bonds if not 1.2<=b['CN_A']<=1.5]
  obad=[b for b in bonds if b['trans_deviation_deg']>30]
  coordinate_record=[(list(k),at,xyz.tolist()) for k in order for at,xyz in sorted(a[k].items())]
  ch=hashlib.sha256(json.dumps(coordinate_record,separators=(',',':')).encode()).hexdigest()
  rows.append(dict(name=name,seed=seed,status='STRICT_STATIC_GEOMETRY_AUDIT_NOT_ACTIVITY',source=str(src),source_sha256=audit['pdb_sha256'],coordinate_sha256=ch,
                   technical_status=audit['status'],selected_Asp1_oxygen=o,selected_Asp2_oxygen=o2,distances_A=ds,distance_gap_sum_A=sum(gaps),
                   distance_pass=not any(gaps),NH_contacts=contacts,direction_pass=direction,
                   TDD_distance_direction_pass=not any(gaps) and direction,CN_bad_count=len(bad),CN_bad_pairs=bad,omega_bad_count=len(obad),omega_bad_pairs=obad,
                   full_internal_clash_audit='NOT_EVALUATED_BY_THIS_SCRIPT',donor_material_joint='NOT_EVALUATED_BY_THIS_SCRIPT'))
duplicates=[]
for name in manifest['ids']:
 rs=[r for r in rows if r['name']==name and 'coordinate_sha256' in r]
 if len(rs)==2:duplicates.append(dict(name=name,identical_coordinates=rs[0]['coordinate_sha256']==rs[1]['coordinate_sha256']))
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
out=B/('strict_geometry_'+stamp+'.json')
result=dict(status='COMPLETE_GEOMETRY_AUDIT' if len(rows)==10 and not pending else 'SNAPSHOT_GEOMETRY_AUDIT',expected=10,records=rows,pending=pending,coordinate_replication=duplicates,
            limitations=['Technical relaxed-state result, not experimental maturation or MD ensemble','Same-coordinate seeds are not independent conformations','Strict distance/direction/CN/omega only; no enzyme activity inference'])
with out.open('x') as f:json.dump(result,f,indent=2)
with (R/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=stamp,script=str(Path(__file__)),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),output=str(out),exit_code=0,status=result['status']))+'\n')
print(json.dumps(dict(output=str(out),status=result['status'],audited=len(rows),TDD_pass=sum(r.get('TDD_distance_direction_pass',False) for r in rows),replication=duplicates,
                     rows=[{k:r.get(k) for k in ['name','seed','distances_A','distance_gap_sum_A','CN_bad_count','omega_bad_count']} for r in rows])))
