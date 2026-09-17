"""Reuse pinned exploratory AF2 protocol after two MPNN jobs succeed; no retries."""
from pathlib import Path
import json,hashlib,datetime,time,sys
r=Path(__file__).resolve().parent
b=r/'route56_N10_20260915_v1/route6_top50_local_generation_v1/two_repaired_exploratory_mpnn20_v1'
source=r/'predict_supported_arg_core4_control_v1.py'
assert hashlib.sha256(source.read_bytes()).hexdigest()=='908c0d486d3f24d29ff402551882b2d11e885f60813812f00aa237757ffd7f75'
assert hashlib.sha256((r/'route6_lowcov_af40_v1/config.json').read_bytes()).hexdigest()=='22a55d5cedd27bdfe9fdc9965c6292a98869b3908029a32ca954043c892a0521'
manifest=json.loads((b/'manifest.json').read_text())
assert len(manifest['jobs'])==2
programs=[]
for j in manifest['jobs']:
 name=j['name'];cid=name.split('_')[0]
 text=source.read_text()
 old="outroot/'supported_ARG_SC4559268_mpnn20_core4_v1/seqs/supported_ARG_SC4559268_mpnn20_core4_v1.fa'"
 new="Path("+repr(str(Path(j['output'])/'seqs'/(name+'.fa')))+")"
 assert text.count(old)==1
 text=text.replace(old,new)
 replacements={
 'R6_SUP_CORE4_':'R6_REP_'+cid+'_',
 'route6_supported_arg_core4_af20_v1':'route6_repaired_'+cid+'_af20_v1',
 'route6_supported_arg_core4_msa20_v1':'route6_repaired_'+cid+'_msa20_v1',
 'af2_route6_supported_arg_core4_20_v1':'af2_route6_repaired_'+cid+'_20_v1',
 'supported_arg_core4_af20_':'repaired_'+cid+'_af20_'}
 for old,new in replacements.items():
  assert old in text
  text=text.replace(old,new)
 compile(text,str(source), 'exec')
 programs.append((name,text))
if '--check' in sys.argv:
 print(json.dumps(dict(status='READ_ONLY_CONFIGURATION_CHECK',targets=[n for n,t in programs],sequences=40,model='alphafold2_ptm',models=1,seeds=1,recycles=3,gpu='existing_H100_config',no_retry=True)))
 raise SystemExit(0)
with (b/'prediction_controller.claim').open('x') as f:f.write(str(__file__))
def log(step,**kw):
 row=dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),step=step,**kw)
 for p in [b/'prediction_events.jsonl',r/'RUN_LOG.jsonl']:
  with p.open('a') as f:f.write(json.dumps(row)+'\n')
with (b/'RUNBOOK.md').open('a') as f:
 f.write('\nIndependent prediction: run predict_two_repaired_after_mpnn_v1.py with project CPU Python. --check is read-only. Waits for both MPNN jobs to succeed, validates unique 20 sequences per target with fixed TDDR, then reuses SHA-pinned existing MSA and AF2 protocol, serially on existing H100 configuration via gpurun. All outputs are new per-candidate directories. Stops without retry on failure or missing precursor process. Full fold/mature catalytic and material audits remain required.\n')
log('repaired_af_waiting_for_mpnn',script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
start=time.monotonic()
while True:
 events=[json.loads(x) for x in (b/'events.jsonl').read_text().splitlines()]
 if any(e['step']=='exploratory_mpnn_batch_complete_NOT_PREDICTED' for e in events):break
 if any(e['step']=='exploratory_mpnn_exit' and e['exit_code']!=0 for e in events):
  log('repaired_af_NOT_SUBMITTED_mpnn_failed');raise SystemExit(2)
 proc=Path('/proc/1777986/cmdline')
 if not proc.exists() or b'run_two_exploratory_mpnn_after_recovery_v1.py' not in proc.read_bytes():
  log('repaired_af_NOT_SUBMITTED_predecessor_absent');raise SystemExit(2)
 if time.monotonic()-start>46800:
  log('repaired_af_NOT_SUBMITTED_wait_expired');raise SystemExit(2)
 time.sleep(30)
try:
 for name,text in programs:
  log('repaired_af_target_start',name=name,derived_script_sha256=hashlib.sha256(text.encode()).hexdigest())
  exec(compile(text,str(source),'exec'),{'__file__':str(source),'__name__':'__main__'})
  log('repaired_af_target_technical_complete_NOT_GEOMETRY_PASS',name=name)
except BaseException as exc:
 log('repaired_af_stopped_no_retry',error_type=type(exc).__name__,error=str(exc));raise
log('repaired_af_batch_technical_complete_NOT_GEOMETRY_PASS',targets=2,sequences=40)
