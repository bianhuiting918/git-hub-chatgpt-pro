"""Resume only verified ready MSAs; no automatic retries or prediction parameter changes."""
from pathlib import Path
import json,hashlib,shutil,subprocess,datetime
r=Path(__file__).resolve().parent;g=r/'route6_supported_arg_core4_af20_v1'
c=json.loads((g/'config.json').read_text());qs=json.loads((g/'sequence_mapping.json').read_text());msa=Path(c['msa'][9]);out=Path(c['af2'][9])
assert not list(out.glob('*.pdb')) and not list(out.glob('*.done.txt'))
ready=[];missing=[]
for n,s in qs.items():
 p=msa/(n+'.a3m')
 if not p.exists():missing.append(n);continue
 q='';started=False
 for l in p.read_text().splitlines():
  if l.startswith('#'):continue
  if l.startswith('>'):
   if started:break
   started=True
  elif started:q+=l.strip()
 assert q==s,('QUERY_MISMATCH',n)
 ready.append((n,p))
assert len(ready)==17 and set(missing)=={'R6_SUP_CORE4_01','R6_SUP_CORE4_02','R6_SUP_CORE4_17'}
inp=g/'msa_inputs_ready17_v1';inp.mkdir(exist_ok=False)
for n,p in ready:shutil.copy2(p,inp/p.name)
cmd=c['af2'].copy();cmd[8]=str(inp);assert cmd[0]=='/usr/local/bin/gpurun'
manifest=dict(status='READY_SUBSET_NOT_COMPLETE',total_expected=20,ready=[n for n,p in ready],missing=missing,command=cmd,msa_hashes={n:hashlib.sha256(p.read_bytes()).hexdigest() for n,p in ready},source_config_sha256=hashlib.sha256((g/'config.json').read_bytes()).hexdigest())
with (g/'resume_ready17_manifest_v1.json').open('x') as h:json.dump(manifest,h,indent=2)
with (g/'RUNBOOK.md').open('a') as h:h.write('\nRecovery: original MSA command exited zero but 01,02,17 lack a3m after DNS failures. resume_core4_ready17_v1.py validates query identity and predicts remaining 17 with unchanged parameters. Missing three remain NOT_EVALUATED; no denominator change.\n')
def log(stage,**kw):
 rec=dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),step=stage,**kw)
 for p in [g/'events.jsonl',r/'RUN_LOG.jsonl']:
  with p.open('a') as h:h.write(json.dumps(rec)+'\n')
log('core4_ready17_af2_start',**manifest)
with (g/'af2_ready17_v1.stdout').open('x') as h:
 p=subprocess.Popen(cmd,cwd=r,stdout=h,stderr=subprocess.STDOUT)
 log('core4_ready17_af2_process',pid=p.pid)
 rc=p.wait()
log('core4_ready17_af2_exit',exit_code=rc)
if rc:raise SystemExit(rc)
for n,p in ready:assert (out/(n+'.done.txt')).exists() and len(list(out.glob(n+'_unrelaxed*.pdb')))==1
log('core4_ready17_technical_complete',completed=17,total_expected=20,missing=missing)
