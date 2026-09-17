"""Six previously unsubmitted donor candidates, three seeds each.
Starts only after the existing two-job step diagnostic completes successfully.
"""
from pathlib import Path
import json,hashlib,datetime,subprocess,time
r=Path(__file__).resolve().parent;b=r/'route56_N10_20260915_v1/route6_top50_local_generation_v1'
o=b/'recovered6_RFD_v1';o.mkdir(exist_ok=False)
spec=[('SC4559268',1,True),('SC4934131',3,False),('SC4677763',3,False),('SC4402181',2,False),('SC4614163',3,False),('SC5011126',1,True)]
jobs=[]
for cid,var,repair in spec:
 original=json.loads((b/'candidate_specific_shapes_v2'/cid/f'PA6_{cid}_coverage10_variant1/manifest.json').read_text())
 assert original['status']=='NONMANIFOLD_MESH_REJECTED'
 name=f'PA6_{cid}_coverage10_variant{var}'
 p=b/'zero_level_mesh_repair_v2'/name/'manifest.json' if repair else b/'candidate_specific_shapes_v2'/cid/name/'manifest.json'
 j=json.loads(p.read_text())
 assert j['status']==('NUMERICAL_MESH_REPAIR_PREPARED_NOT_SUBMITTED' if repair else 'PREPARED_NOT_SUBMITTED')
 for key in ['pdb','shape']:assert hashlib.sha256(Path(j[key]).read_bytes()).hexdigest()==j[key+'_sha256']
 out=r.parent/'n9_core_donor_design_20260910_gpu_outputs_v1/route6_recovered6_three_seeds_v1'/name
 cmd=[('inference.output_prefix='+str(out/'design')) if a.startswith('inference.output_prefix=') else 'inference.num_designs=3' if a.startswith('inference.num_designs=') else a for a in j['command']]
 assert cmd[0]=='/usr/local/bin/gpurun' and 'CUDA_VISIBLE_DEVICES=1' in cmd and 'diffuser.T=50' in cmd
 jobs.append(dict(name=name,core=cid,variant=var,numerical_mesh_repair=repair,input_manifest=str(p),manifest_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),command=cmd,output=str(out),comparison_note='Alternative subset or numerical repair, not identical interface to failed variant1'))
with (o/'manifest.json').open('x') as f:json.dump(dict(status='WAITING_EXISTING_STEP_DIAGNOSTIC',jobs=jobs,expected_backbones=18,concurrency=1),f,indent=2)
def log(stage,**kw):
 x=dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),step=stage,**kw)
 for p in [o/'events.jsonl',r/'RUN_LOG.jsonl']:
  with p.open('a') as f:f.write(json.dumps(x)+'\n')
(o/'RUNBOOK.md').write_text('Six donors omitted by the original batch because variant1 had invalid mesh. Four use already valid alternate PA6 10% subsets; two use strictly checked numerical mesh repairs. Three seeds each, T50 unchanged. Start only after step-count diagnostic completed with success; no automatic retries. At most one child of this controller, plus existing main batch. Do not modify main batch. Preserve nonzero exits and outputs. Compare interface differences explicitly; no catalytic PASS implied.\n')
log('recovered6_waiting_for_step_diagnostic',jobs=6)
start=time.monotonic()
while True:
 events=[json.loads(x) for x in (b/'step_count_diagnostic_v1/events.jsonl').read_text().splitlines()]
 if any(x['step']=='step_diagnostic_batch_complete' for x in events):break
 if any(x['step']=='step_diagnostic_exit' and x['exit_code']!=0 for x in events):raise RuntimeError('Step diagnostic failed; recovery generation not launched')
 proc=Path('/proc/1504746/cmdline')
 if not proc.exists() or b'run_steps100_200_diagnostic_v1.py' not in proc.read_bytes():raise RuntimeError('Step diagnostic controller absent without successful completion')
 if time.monotonic()-start>3600:raise RuntimeError('Wait limit reached, no submission')
 time.sleep(30)
for j in jobs:
 assert hashlib.sha256(Path(j['input_manifest']).read_bytes()).hexdigest()==j['manifest_sha256']
 out=Path(j['output']);out.mkdir(parents=True,exist_ok=False,mode=0o2770);out.chmod(0o2770)
 with (o/(j['name']+'.stdout')).open('x') as f:
  p=subprocess.Popen(j['command'],cwd='/data/bht2/shape_conditioned_rfdiffusion_20260905_v1/source/RFdiffusion',stdout=f,stderr=subprocess.STDOUT)
  log('recovered_candidate_start',name=j['name'],pid=p.pid,command=j['command']);rc=p.wait()
 log('recovered_candidate_exit',name=j['name'],exit_code=rc)
 if rc:raise SystemExit(rc)
 for seed in range(3):
  p=out/f'design_{seed}.pdb';assert p.exists() and p.stat().st_size>1000
  log('recovered_backbone_technical_complete_NOT_GEOMETRY_PASS',name=j['name'],seed=seed,path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest())
log('recovered6_batch_complete_NOT_GEOMETRY_PASS',backbones=18)
