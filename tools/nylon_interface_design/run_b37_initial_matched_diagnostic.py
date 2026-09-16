"""Missing original B37_1 diagnostic using existing matched CPU protocol.
Run with existing CPU Python from project root; two worker processes, one thread each.
No sequence changes or catalytic restraints. New directory only; no automatic retries.
"""
import os
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
from pathlib import Path
import json,hashlib,sys,subprocess,datetime,concurrent.futures,importlib.util
r=Path(__file__).resolve().parent
base=r/'matched_cleavage_original_B37_1_v1'
module=Path('/data/bht2/new_seed_structure_sequence_folddisco_20260825_v3/n9_n10_af_repredict_20260908_v1/scripts/local_relax_nonfine_matched_v1.py')
python='/data/bht2/new_seed_structure_sequence_folddisco_20260825_v3/nylon_maturation_relax_20260906_v1/env/bin/python'
def log(rec):
 rec['time']=datetime.datetime.now(datetime.timezone.utc).isoformat()
 with (r/'RUN_LOG.jsonl').open('a') as h:h.write(json.dumps(rec)+'\n')
if len(sys.argv)>1 and sys.argv[1]=='worker':
 spec=importlib.util.spec_from_file_location('matched',module);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m)
 row=json.loads((base/'inputs/N9_R4_B37_1/preparation.json').read_text())
 result=m.run(row,int(sys.argv[2]),base/'outputs')
 raise SystemExit(0 if result['status']=='LOCAL_RELAX_TECHNICAL_PASS' else 2)
def cleave(p):
 rows=[];switched=False
 for l in p.read_text().splitlines():
  if not l.startswith('ATOM'):continue
  if l[76:78].strip()=='H':continue
  n=int(l[22:26])
  if n>=90 and not switched:rows.append('TER');switched=True
  rows.append(l[:21]+('A' if n<=89 else 'B')+l[22:])
 return '\n'.join(rows+['TER','END'])+'\n'
# Reproduce the existing validated preparation before constructing new input.
old=r/'matched_cleavage_route4_nearpair_20260916_v1/inputs/N9_R4_B37_8'
oldrow=json.loads((old/'preparation.json').read_text())
assert cleave(Path(oldrow['source'])).splitlines()==(old/'cleaved_unit_raw.pdb').read_text().splitlines()
audit=json.loads((r/'route4_fixed_donor_pair_20260916T181221Z.json').read_text())
v=next(v for v in audit['records'] if v['family']=='initial130' and v['name']=='N9_R4_B37_1')
src=Path(v['path']);assert hashlib.sha256(src.read_bytes()).hexdigest()==v['sha256']
base.mkdir(exist_ok=False);prepared=base/'inputs/N9_R4_B37_1';prepared.mkdir(parents=True)
(prepared/'cleaved_unit_raw.pdb').write_text(cleave(src))
row=dict(site_id='N9_R4_B37_1',seed='N9_R4_B37_1',nres=169,core_pose_indices='90;129;131',cut_after_pose_index=89,prepared_dir=str(prepared),source=str(src),source_sha256=v['sha256'],status='DIAGNOSTIC_NOT_EXPERIMENTAL_MATURATION')
(prepared/'preparation.json').write_text(json.dumps(row,indent=2))
(base/'RUNBOOK.md').write_text('CPU diagnostic. Existing N9_N10_MATCHED_NONFINE_V1 protocol unchanged. Two seeds 101,202; one thread each. Source and module hashes in manifest. No catalytic geometry restraints. Strict catalytic, donor and CN audit required separately. Do not rerun into this directory.\n')
(base/'manifest.json').write_text(json.dumps(dict(input=row,module=str(module),module_sha256=hashlib.sha256(module.read_bytes()).hexdigest(),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),seeds=[101,202]),indent=2))
def run(seed):
 cmd=[python,str(Path(__file__)), 'worker',str(seed)]
 with (base/f'seed_{seed}.stdout').open('x') as h:
  p=subprocess.Popen(cmd,cwd=r,stdout=h,stderr=subprocess.STDOUT)
  log(dict(step='original_B37_matched_start',seed=seed,pid=p.pid,command=cmd,output=str(base)))
  code=p.wait()
 rec=dict(step='original_B37_matched_finished',seed=seed,exit_code=code,output=str(base));log(rec);return rec
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:results=list(pool.map(run,[101,202]))
(base/'completion.json').write_text(json.dumps(results,indent=2));print(json.dumps(results))
raise SystemExit(0 if all(v['exit_code']==0 for v in results) else 2)
