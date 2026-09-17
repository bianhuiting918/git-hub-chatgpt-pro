"""Complete five missing route2 matched-maturation diagnostics; CPU, unchanged protocol."""
import os
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
from pathlib import Path
import json,hashlib,sys,subprocess,datetime,concurrent.futures,importlib.util
R=Path(__file__).resolve().parent
B=R/'matched_cleavage_route2_sameO_five_v1'
M=Path('/data/bht2/new_seed_structure_sequence_folddisco_20260825_v3/n9_n10_af_repredict_20260908_v1/scripts/local_relax_nonfine_matched_v1.py')
MSHA='a1cd4de44eabfaad566f0e94c12e3fee9c17f54657ce218a4d8b9ac6f724a5d0'
PY='/data/bht2/new_seed_structure_sequence_folddisco_20260825_v3/nylon_maturation_relax_20260906_v1/env/bin/python'
assert hashlib.sha256(M.read_bytes()).hexdigest()==MSHA
def log(**kw):
 kw['time']=datetime.datetime.now(datetime.timezone.utc).isoformat()
 for p in [R/'RUN_LOG.jsonl',B/'events.jsonl']:
  with p.open('a') as f:f.write(json.dumps(kw)+'\n')
if len(sys.argv)>1 and sys.argv[1]=='worker':
 name,seed=sys.argv[2],int(sys.argv[3])
 row=json.loads((B/'inputs'/name/'preparation.json').read_text())
 assert hashlib.sha256(Path(row['source']).read_bytes()).hexdigest()==row['source_sha256']
 spec=importlib.util.spec_from_file_location('matched',M);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 result=m.run(row,seed,B/'outputs')
 raise SystemExit(0 if result['status']=='LOCAL_RELAX_TECHNICAL_PASS' else 2)
def cleave(p):
 lines=[];switched=False
 for l in p.read_text().splitlines():
  if not l.startswith('ATOM') or l[76:78].strip()=='H':continue
  n=int(l[22:26])
  if n>=90 and not switched:lines.append('TER');switched=True
  lines.append(l[:21]+('A' if n<=89 else 'B')+l[22:])
 return '\n'.join(lines+['TER','END'])+'\n'
old=R/'matched_cleavage_e07_12_20260916_v1/inputs/N9_route2_e07_12'
oldrow=json.loads((old/'preparation.json').read_text())
assert cleave(Path(oldrow['source'])).splitlines()==(old/'cleaved_unit_raw.pdb').read_text().splitlines()
src=R/'route123_same_Asp_oxygen_audit_20260917T062707422842Z.json'
data=json.loads(src.read_text())
rows=[r for r in data['records'] if r['route']=='route2' and r['distance_gate']]
assert {r['id'] for r in rows}=={'N9_route2_e08_13','N9_route2_e17_05','N9_route2_e04_01','N9_route2_e09_03','N9_route2_e02_07'}
# Validate all input identities before creating any run directory.
for r in rows:
 p=Path(r['prediction']);assert hashlib.sha256(p.read_bytes()).hexdigest()==r['prediction_sha256']
 residues={}
 for l in p.read_text().splitlines():
  if l.startswith('ATOM'):residues[int(l[22:26])]=l[17:20]
 assert sorted(residues)==list(range(1,len(residues)+1))
 assert [residues[i] for i in [90,129,131]]==['THR','ASP','ASP']
 r['nres']=len(residues)
if '--check' in sys.argv:
 print(json.dumps(dict(status='READ_ONLY_INPUT_CHECK_PASS',candidates=len(rows),workers=2,seeds=[101,202],module_sha256=MSHA)))
 raise SystemExit(0)
B.mkdir(exist_ok=False)
for r in rows:
 d=B/'inputs'/r['id'];d.mkdir(parents=True)
 (d/'cleaved_unit_raw.pdb').write_text(cleave(Path(r['prediction'])))
 row=dict(site_id=r['id'],seed=r['id'],nres=r['nres'],core_pose_indices='90;129;131',cut_after_pose_index=89,prepared_dir=str(d),source=r['prediction'],source_sha256=r['prediction_sha256'],status='MATCHED_MATURATION_DIAGNOSTIC_NOT_ACTIVITY')
 (d/'preparation.json').write_text(json.dumps(row,indent=2))
(B/'manifest.json').write_text(json.dumps(dict(source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),module=str(M),module_sha256=MSHA,ids=[r['id'] for r in rows],seeds=[101,202],workers=2,script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2))
(B/'RUNBOOK.md').write_text('Run this script with existing CPU Python. Five missing route2 candidates, two seeds each, two CPU workers one thread each. Identical N9_N10_MATCHED_NONFINE_V1 parameters; no catalytic restraints or sequence design. New directory only. Original stage outputs retained. Technical PASS is not strict geometry PASS: separately audit same-Asp-oxygen distances, H directions, CN 1.2-1.5 A, omega and sterics. No GPU, no automatic retries.\n')
def run(job):
 name,seed=job;cmd=[PY,str(Path(__file__)), 'worker',name,str(seed)]
 with (B/(name+'_'+str(seed)+'.stdout')).open('x') as f:
  p=subprocess.Popen(cmd,cwd=R,stdout=f,stderr=subprocess.STDOUT)
  log(step='route2_sameO_maturation_start',name=name,seed=seed,pid=p.pid,command=cmd)
  code=p.wait()
 rec=dict(step='route2_sameO_maturation_finished',name=name,seed=seed,exit_code=code)
 log(**rec);return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
 results=list(pool.map(run,[(r['id'],seed) for r in rows for seed in [101,202]]))
(B/'completion.json').write_text(json.dumps(results,indent=2))
log(step='route2_sameO_maturation_batch_complete_REQUIRES_GEOMETRY_AUDIT',completed=len(results),nonzero=sum(r['exit_code']!=0 for r in results))
raise SystemExit(0 if all(r['exit_code']==0 for r in results) else 2)
