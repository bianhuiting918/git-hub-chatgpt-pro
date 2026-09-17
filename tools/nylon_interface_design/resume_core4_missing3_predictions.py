"""Resume only the three verified recovered queries; unchanged AF parameters."""
from pathlib import Path
import json,hashlib,subprocess,datetime,shutil
r=Path(__file__).resolve().parent;g=r/'route6_supported_arg_core4_af20_v1'
c=json.loads((g/'config.json').read_text());qs=json.loads((g/'sequence_mapping.json').read_text())
rec=json.loads((g/'missing3_msa_recovery_v1/completion.json').read_text())
assert rec['status']=='MSA_VERIFIED_NOT_PREDICTED' and not rec['missing']
names={x['name'] for x in rec['ready']};assert names=={'R6_SUP_CORE4_01','R6_SUP_CORE4_02','R6_SUP_CORE4_17'}
out=Path(c['af2'][9]);events=[json.loads(x) for x in (g/'events.jsonl').read_text().splitlines()]
assert any(x.get('step')=='core4_ready17_technical_complete' for x in events)
assert len(list(out.glob('*.done.txt')))==17
inp=g/'msa_inputs_recovered3_v1';inp.mkdir(exist_ok=False)
for x in rec['ready']:
 n=x['name'];p=Path(x['path']);assert hashlib.sha256(p.read_bytes()).hexdigest()==x['sha256']
 assert not (out/(n+'.done.txt')).exists() and not list(out.glob(n+'_unrelaxed*.pdb'))
 q='';started=False
 for line in p.read_text().splitlines():
  if line.startswith('#'):continue
  if line.startswith('>'):
   if started:break
   started=True
  elif started:q+=line.strip()
 assert q==qs[n]
 shutil.copy2(p,inp/p.name)
cmd=c['af2'].copy();cmd[8]=str(inp);assert cmd[0]=='/usr/local/bin/gpurun'
with (g/'resume_recovered3_manifest_v1.json').open('x') as h:json.dump(dict(command=cmd,queries=sorted(names),source=rec,config_sha256=hashlib.sha256((g/'config.json').read_bytes()).hexdigest()),h,indent=2)
def log(stage,**kw):
 x=dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),step=stage,**kw)
 for p in [g/'events.jsonl',r/'RUN_LOG.jsonl']:
  with p.open('a') as h:h.write(json.dumps(x)+'\n')
with (g/'RUNBOOK.md').open('a') as h:h.write('\nresume_core4_missing3_predictions_v1.py validates recovered 01/02/17 query identities and hashes, preserves completed seventeen, reuses AF command parameters and records exit. No automatic retry.\n')
with (g/'af2_recovered3_v1.stdout').open('x') as h:
 p=subprocess.Popen(cmd,cwd=r,stdout=h,stderr=subprocess.STDOUT);log('core4_recovered3_af2_start',pid=p.pid,command=cmd);rc=p.wait()
log('core4_recovered3_af2_exit',exit_code=rc)
if rc:raise SystemExit(rc)
for n in names:assert (out/(n+'.done.txt')).exists() and len(list(out.glob(n+'_unrelaxed*.pdb')))==1
log('core4_all20_technical_complete',completed=20)
subprocess.run([str(r.parent/'unbiased_100chain_v1/env/bin/python'),str(r/'audit_supported_arg_core4_predictions_v1.py')],cwd=r,check=True)
