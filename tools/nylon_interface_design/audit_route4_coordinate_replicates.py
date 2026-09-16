"""Read-only source audit of exact coordinate duplicates; append versioned report.
Run with existing CPU Python from N9 design project directory. No recomputation.
"""
from pathlib import Path
from collections import defaultdict,Counter
import json,hashlib,datetime
r=Path(__file__).resolve().parent
src=r/'route4_expanded_matched_snapshot_20260916T011223Z.json'
j=json.loads(src.read_text());rows=[];groups=defaultdict(list)
for v in j['records']:
 p=Path(v['path']);row={k:v.get(k) for k in ['name','family','stage','seed','path','distance_gate','direction_diagnostic']}
 if not p.exists():row['status']='NOT_EVALUATED_MISSING';rows.append(row);continue
 data=p.read_bytes();row['sha256']=hashlib.sha256(data).hexdigest()
 if row['sha256']!=v['source_sha256']:
  row['status']='NOT_EVALUATED_SOURCE_CHANGED';rows.append(row);continue
 lines=[s for s in data.decode().splitlines() if s.startswith(('ATOM  ','HETATM'))]
 # Exact serialized atom identities and coordinates, excluding scores, headers, B and occupancy.
 payload='\n'.join(s[:6]+s[12:27]+s[30:54] for s in lines).encode()
 row.update(status='AUDITED',atom_count=len(lines),coordinate_sha256=hashlib.sha256(payload).hexdigest())
 rows.append(row);groups[(v['family'],v['name'],v['stage'])].append(row)
 pairs=[]
for k,vs in groups.items():
 if len(vs)>1:
  pairs.append(dict(family=k[0],name=k[1],stage=k[2],records=len(vs),seeds=[v['seed'] for v in vs],unique_coordinates=len(set(v['coordinate_sha256'] for v in vs)),coordinate_hashes=[v['coordinate_sha256'] for v in vs]))
summary=[]
for family,stage in sorted(set((v['family'],v['stage']) for v in rows)):
 vs=[v for v in rows if v['family']==family and v['stage']==stage];good=[v for v in vs if v['status']=='AUDITED']
 ps=[p for p in pairs if p['family']==family and p['stage']==stage]
 summary.append(dict(family=family,stage=stage,records=len(vs),audited=len(good),unique_coordinates=len(set(v['coordinate_sha256'] for v in good)),multi_seed_groups=len(ps),identical_coordinate_groups=sum(p['unique_coordinates']==1 for p in ps),distance_pass_records=sum(bool(v['distance_gate']) for v in good),joint_distance_direction_records=sum(bool(v['distance_gate']) and bool(v['direction_diagnostic']) for v in good)))
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ');out=r/('route4_coordinate_replication_audit_'+stamp+'.json')
report=dict(status='COORDINATE_REPLICATION_AUDIT_NOT_ACTIVITY',input=str(src),input_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),summary=summary,groups=pairs,records=rows,limitations=['Exact PDB coordinate equality at stored precision, not RMSD clustering','Different coordinates do not establish statistically independent sampling','Existing geometric flags inherited only after source hash validation','No new molecular calculations and no changed acceptance gates'])
with out.open('x') as h:json.dump(report,h,indent=2)
with (r/'RUN_LOG.jsonl').open('a') as h:h.write(json.dumps(dict(time=stamp,script=str(Path(__file__)),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),input=str(src),output=str(out),exit_code=0))+'\n')
print(json.dumps(dict(output=str(out),summary=summary,status_counts=dict(Counter(v['status'] for v in rows))),indent=2))
