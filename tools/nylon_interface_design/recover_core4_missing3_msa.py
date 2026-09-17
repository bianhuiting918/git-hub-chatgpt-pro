"""One bounded recovery attempt for missing MSA outputs. No GPU retry loop."""
from pathlib import Path
import json,datetime,subprocess,hashlib,socket
r=Path(__file__).resolve().parent;g=r/'route6_supported_arg_core4_af20_v1'
c=json.loads((g/'config.json').read_text());qs=json.loads((g/'sequence_mapping.json').read_text())
names=['R6_SUP_CORE4_01','R6_SUP_CORE4_02','R6_SUP_CORE4_17']
old=Path(c['msa'][9]);assert all(not (old/(n+'.a3m')).exists() for n in names)
socket.getaddrinfo('api.colabfold.com',443,type=socket.SOCK_STREAM)
base=g/'missing3_msa_recovery_v1';base.mkdir(exist_ok=False)
fa=base/'queries.fasta';fa.write_text(''.join('>'+n+'\n'+qs[n]+'\n' for n in names))
out=old.parent/'route6_supported_arg_core4_missing3_msa_v1';out.mkdir(mode=0o2770,exist_ok=False);out.chmod(0o2770)
cmd=c['msa'].copy();cmd[8]=str(fa);cmd[9]=str(out);assert cmd[0]=='/usr/local/bin/gpurun'
manifest=dict(names=names,total_expected=20,command=cmd,query_sha256=hashlib.sha256(fa.read_bytes()).hexdigest(),original_config_sha256=hashlib.sha256((g/'config.json').read_bytes()).hexdigest(),reason='Original MSA DNS failures; 17 valid MSAs retained')
(base/'manifest.json').write_text(json.dumps(manifest,indent=2))
(base/'RUNBOOK.md').write_text('Run recovery script once with project CPU Python. Same MSA parameters, only queries 01,02,17. Separate outputs preserve original failure evidence. No AF2 launched here; validate a3m query identities first. No automatic task retry.\n')
def log(step,**kw):
 rec=dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),step=step,**kw)
 for p in [base/'events.jsonl',r/'RUN_LOG.jsonl']:
  with p.open('a') as h:h.write(json.dumps(rec)+'\n')
log('missing3_msa_start',**manifest)
with (base/'msa.stdout').open('x') as h:
 p=subprocess.Popen(cmd,cwd=r,stdout=h,stderr=subprocess.STDOUT);log('missing3_msa_process',pid=p.pid);rc=p.wait()
log('missing3_msa_exit',exit_code=rc)
if rc:raise SystemExit(rc)
ready=[];missing=[]
for n in names:
 f=out/(n+'.a3m')
 if not f.exists():missing.append(n);continue
 q='';start=False
 for l in f.read_text().splitlines():
  if l.startswith('#'):continue
  if l.startswith('>'):
   if start:break
   start=True
  elif start:q+=l.strip()
 assert q==qs[n],('QUERY_MISMATCH',n)
 ready.append(dict(name=n,path=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest()))
report=dict(status='MSA_VERIFIED_NOT_PREDICTED' if not missing else 'MSA_INCOMPLETE',ready=ready,missing=missing)
(base/'completion.json').write_text(json.dumps(report,indent=2));log('missing3_msa_verified',**report)
raise SystemExit(0 if not missing else 2)
