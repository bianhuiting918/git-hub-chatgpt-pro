"""Run the prepared zero-guide diagnostic after both existing RFD batches exit."""
from pathlib import Path
import json,hashlib,time,datetime,subprocess,sys
R=Path(__file__).resolve().parent
B=R/'route56_N10_20260915_v1/route6_top50_local_generation_v1'
O=B/'SC5150928_zero_guide_diagnostic_v1'
mf=O/'manifest.json'
assert hashlib.sha256(mf.read_bytes()).hexdigest()=='8553facc530ca633c39d739cead29de8e807d5341acacfaae36c49b0f1d97145'
m=json.loads(mf.read_text());cmd=m['command']
assert cmd[0]=='/usr/local/bin/gpurun' and 'potentials.guide_scale=0' in cmd and 'inference.num_designs=3' in cmd
def log(step,**kw):
 rec=dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),step=step,**kw)
 for p in [O/'events.jsonl',R/'RUN_LOG.jsonl']:
  with p.open('a') as f:f.write(json.dumps(rec)+'\n')
def state(folder,terminal,pid,script):
 es=[json.loads(l) for l in (B/folder/'events.jsonl').read_text().splitlines()]
 bad=[e for e in es if e.get('step') in ['candidate_rfd_exit','recovered_candidate_exit'] and e.get('exit_code',0)!=0]
 if bad:raise RuntimeError('Upstream failure; no retry: '+str(bad[-1]))
 if any(e.get('step')==terminal for e in es):
  return True
 p=Path('/proc')/str(pid)/'cmdline'
 if not p.exists():raise RuntimeError('Upstream controller absent without completion: '+script)
 if script.encode() not in p.read_bytes():raise RuntimeError('Upstream PID identity mismatch')
 return False
def ready():
 a=state('first19_RFD3_v2','first19_batch_finished',1471149,'run_top19_three_backbones_v2.py')
 b=state('recovered6_RFD_v1','recovered6_batch_complete_NOT_GEOMETRY_PASS',1618827,'run_six_recovered_candidates_v1.py')
 return a and b
if '--check' in sys.argv:
 print(json.dumps(dict(status='READ_ONLY_QUEUE_CHECK_PASS',both_RFD_complete=ready(),gpu_submitted=False)));raise SystemExit(0)
with (O/'controller.claim').open('x') as f:f.write(str(__import__('os').getpid()))
log('zero_guide_waiting_for_both_RFD_batches')
start=time.monotonic()
while not ready():
 if time.monotonic()-start>12*3600:raise RuntimeError('Upstream wait timeout; diagnostic not submitted')
 time.sleep(30)
# Recheck immutable inputs and all baseline structures immediately before launch.
for key in ['input_pdb','shape']:
 assert hashlib.sha256(Path(m[key]).read_bytes()).hexdigest()==m[key+'_sha256']
for b in m['baseline']:assert hashlib.sha256(Path(b['path']).read_bytes()).hexdigest()==b['sha256']
# Check outputs, not stale expected=57 metadata in the main controller.
main=json.loads((B/'first19_RFD3_v2/completion.json').read_text())
assert len(main['completed'])==39
out=Path(next(a.split('=',1)[1] for a in cmd if a.startswith('inference.output_prefix='))).parent
out.mkdir(parents=True,exist_ok=False,mode=0o2770);out.chmod(0o2770)
with (O/'generation.stdout').open('x') as f:
 p=subprocess.Popen(cmd,cwd=m['cwd'],stdout=f,stderr=subprocess.STDOUT)
 log('zero_guide_started',pid=p.pid,command=cmd,output=str(out))
 code=p.wait()
log('zero_guide_exit',exit_code=code)
if code:raise SystemExit(code)
records=[]
for seed in range(3):
 p=out/('design_'+str(seed)+'.pdb');assert p.is_file() and p.stat().st_size>1000
 records.append(dict(seed=seed,path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
with (O/'completion.json').open('x') as f:json.dump(dict(status='TECHNICAL_COMPLETE_NOT_GEOMETRY_PASS',records=records),f,indent=2)
log('zero_guide_complete_REQUIRES_GEOMETRY_AUDIT',backbones=3)
