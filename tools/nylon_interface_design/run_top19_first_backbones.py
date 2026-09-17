"""One first backbone per distinct donor candidate, sequential A100 via gpurun.
Wait only for the active CPU preparation to produce each input. No retry.
Other interface variants remain prepared but not submitted.
"""
import os,sys,time,json,hashlib,subprocess,datetime
from pathlib import Path
r=Path(__file__).resolve().parent
base=r/'route56_N10_20260915_v1/route6_top50_local_generation_v1'
src=base/'local_fragments_three_templates_v1/results.jsonl'
rows=[json.loads(x) for x in src.read_text().splitlines()]
good=[x for x in rows if x['status']=='LOCAL_FRAGMENT_GEOMETRY_PASS_NOT_RFD']
assert len(good)==19 and len({x['id'] for x in good})==19
o=base/'first19_RFD_v1';o.mkdir(exist_ok=False)
shape=base/'candidate_specific_shapes_v1'
def log(stage,**kw):
 x=dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),step=stage,**kw)
 print(json.dumps(x),flush=True)
 for p in [o/'events.jsonl',r/'RUN_LOG.jsonl']:
  with p.open('a') as f:f.write(json.dumps(x)+'\n')
with (o/'manifest.json').open('x') as f:json.dump(dict(candidates=[x['id'] for x in good],designs_per_candidate=1,material='PA6',conditional_fraction=.1,variant=1,concurrency=1,no_retry=True,other_interfaces='PREPARATION_ONLY_NOT_SUBMITTED',purpose='Cross-candidate backbone feasibility screen; not complete route6 coverage'),f,indent=2)
(o/'RUNBOOK.md').write_text('Run from remote N9/N10 project with existing CPU Python. This controller starts only gpurun GPU commands, one at a time, from hash-verified candidate-specific manifests. It waits at most 60 minutes for a CPU-prepared input while confirming the preparation process remains present. Every nonzero GPU command exit stops the batch; no automatic retries. Outputs must not preexist. Preserve failed outputs. Generation technical completion is not geometry PASS. All 19 distinct candidates are prioritized before expanding shape variants.\n')
completed=[]
for row in good:
 cid=row['id'];name=f'PA6_{cid}_coverage10_variant1';mf=shape/cid/name/'manifest.json'
 start=time.monotonic()
 while not mf.exists():
  p=Path('/proc/1441232/cmdline')
  if not p.exists() or b'prepare_top50_fragment_shapes_v1.py' not in p.read_bytes():raise RuntimeError('Preparation stopped before input '+name)
  if time.monotonic()-start>3600:raise RuntimeError('Preparation wait limit '+name)
  time.sleep(30)
 # A published JSON may still be finishing its write; no GPU command has started.
 time.sleep(1)
 j=json.loads(mf.read_text())
 if j['status']!='PREPARED_NOT_SUBMITTED':
  log('candidate_input_not_ready',name=name,status=j['status']);continue
 for key in ['pdb','shape']:assert hashlib.sha256(Path(j[key]).read_bytes()).hexdigest()==j[key+'_sha256']
 cmd=j['command'];assert cmd[0]=='/usr/local/bin/gpurun' and 'CUDA_VISIBLE_DEVICES=1' in cmd
 assert any('ActiveSite_ckpt.pt' in x for x in cmd)
 out=Path(j['output']);out.mkdir(parents=True,exist_ok=False,mode=0o2770);out.chmod(0o2770)
 assert not list(out.iterdir())
 with (o/(name+'.stdout')).open('x') as f:
  p=subprocess.Popen(cmd,cwd='/data/bht2/shape_conditioned_rfdiffusion_20260905_v1',stdout=f,stderr=subprocess.STDOUT)
  log('candidate_rfd_started',name=name,pid=p.pid,command=cmd,manifest_sha256=hashlib.sha256(mf.read_bytes()).hexdigest())
  rc=p.wait()
 log('candidate_rfd_exit',name=name,exit_code=rc)
 if rc:raise SystemExit(rc)
 pdb=out/'design_0.pdb'
 assert pdb.exists() and pdb.stat().st_size>1000
 completed.append(dict(name=name,path=str(pdb),sha256=hashlib.sha256(pdb.read_bytes()).hexdigest()))
 log('candidate_rfd_technical_complete_NOT_GEOMETRY_PASS',**completed[-1])
with (o/'completion.json').open('x') as f:json.dump(dict(status='RFD_TECHNICAL_COMPLETE_REQUIRES_GEOMETRY_AUDIT',expected=19,completed=completed),f,indent=2)
log('first19_batch_finished',completed=len(completed),expected=19)
