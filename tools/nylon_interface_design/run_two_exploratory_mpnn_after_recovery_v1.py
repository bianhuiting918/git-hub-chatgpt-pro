"""Wait for existing recovery batch success, then two serial exploratory MPNN jobs."""
from pathlib import Path
import json,hashlib,subprocess,time,datetime
r=Path(__file__).resolve().parent;b=r/'route56_N10_20260915_v1/route6_top50_local_generation_v1';o=b/'two_repaired_exploratory_mpnn20_v1'
with (o/'controller.claim').open('x') as f:f.write(str(__file__))
manifest=json.loads((o/'manifest.json').read_text())
def log(step,**kw):
 row=dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),step=step,**kw)
 for p in [o/'events.jsonl',r/'RUN_LOG.jsonl']:
  with p.open('a') as f:f.write(json.dumps(row)+'\n')
log('exploratory_mpnn_waiting_for_recovered6')
start=time.monotonic()
while True:
 events=[json.loads(x) for x in (b/'recovered6_RFD_v1/events.jsonl').read_text().splitlines()]
 if any(e['step']=='recovered6_batch_complete_NOT_GEOMETRY_PASS' for e in events):break
 if any(e['step']=='recovered_candidate_exit' and e['exit_code']!=0 for e in events):
  log('exploratory_mpnn_not_submitted_predecessor_failed');raise SystemExit(2)
 proc=Path('/proc/1618827/cmdline')
 if not proc.exists() or b'run_six_recovered_candidates_v1.py' not in proc.read_bytes():
  log('exploratory_mpnn_not_submitted_predecessor_absent');raise SystemExit(2)
 if time.monotonic()-start>43200:
  log('exploratory_mpnn_wait_expired_NOT_SUBMITTED');raise SystemExit(2)
 time.sleep(30)
for j in manifest['jobs']:
 assert hashlib.sha256(Path(j['source']).read_bytes()).hexdigest()==j['source_sha256']
 assert j['command'][0]=='gpurun'
 dest=Path(j['output']);dest.mkdir(parents=True,exist_ok=False,mode=0o2770);dest.chmod(0o2770)
 with (o/(j['name']+'.stdout')).open('x') as f:
  p=subprocess.Popen(j['command'],cwd=r,stdout=f,stderr=subprocess.STDOUT);log('exploratory_mpnn_started',name=j['name'],pid=p.pid,command=j['command']);rc=p.wait()
 log('exploratory_mpnn_exit',name=j['name'],exit_code=rc)
 if rc:raise SystemExit(rc)
 fasta=dest/'seqs'/(j['name']+'.fa')
 lines=fasta.read_text().splitlines();seqs=[lines[i+1] for i,line in enumerate(lines[:-1]) if line.startswith('>') and 'sample=' in line]
 assert len(seqs)==20 and all(len(s)==126 and ''.join(s[i-1] for i in [1,40,42,63])=='TDDR' for s in seqs)
 log('exploratory_mpnn_sequences_complete_NOT_PREDICTED',name=j['name'],sequences=len(seqs),unique_sequences=len(set(seqs)),fasta=str(fasta),sha256=hashlib.sha256(fasta.read_bytes()).hexdigest())
log('exploratory_mpnn_batch_complete_NOT_PREDICTED',targets=2)
