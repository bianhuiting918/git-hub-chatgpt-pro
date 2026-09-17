"""One seed each T100/T200; unchanged motif, shape, checkpoint and other settings.
This tests the full step-count intervention including cumulative shape guidance.
No automatic retries, no overwrites, no activity inference.
"""
from pathlib import Path
import json,hashlib,subprocess,datetime
r=Path(__file__).resolve().parent
b=r/'route56_N10_20260915_v1/route6_top50_local_generation_v1'
j=json.loads((b/'candidate_specific_shapes_v2/SC5150928/PA6_SC5150928_coverage10_variant1/manifest.json').read_text())
assert j['status']=='PREPARED_NOT_SUBMITTED'
for key in ['pdb','shape']:assert hashlib.sha256(Path(j[key]).read_bytes()).hexdigest()==j[key+'_sha256']
o=b/'step_count_diagnostic_v1';o.mkdir(exist_ok=False)
cwd=Path('/data/bht2/shape_conditioned_rfdiffusion_20260905_v1/source/RFdiffusion')
assert (cwd/'scripts/run_inference.py').is_file()
def log(stage,**kw):
 x=dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),step=stage,**kw)
 for p in [o/'events.jsonl',r/'RUN_LOG.jsonl']:
  with p.open('a') as f:f.write(json.dumps(x)+'\n')
jobs=[]
for t in [100,200]:
 out=r.parent/'n9_core_donor_design_20260910_gpu_outputs_v1/route6_step_count_diagnostic_v1'/('T'+str(t))
 cmd=[('diffuser.T='+str(t)) if a=='diffuser.T=50' else ('inference.output_prefix='+str(out/'design')) if a.startswith('inference.output_prefix=') else a for a in j['command']]
 assert cmd[0]=='/usr/local/bin/gpurun' and 'CUDA_VISIBLE_DEVICES=1' in cmd and 'inference.num_designs=1' in cmd
 assert sum(a!=bb for a,bb in zip(cmd,j['command']))==2
 jobs.append(dict(T=t,output=str(out),command=cmd))
with (o/'manifest.json').open('x') as f:json.dump(dict(status='DIAGNOSTIC_NOT_PRODUCTION',jobs=jobs,source=j,seed=0,baseline='existing T50 design_0',confound='Same per-step guide scale, more cumulative guiding updates; not resolution-only causal isolation',stop='Any nonzero exit stops batch without retry',rollback='Existing T50 untouched; diagnostic outputs retained'),f,indent=2)
(o/'RUNBOOK.md').write_text('Execute with existing CPU Python; child commands use gpurun A100. Two bounded diagnostic jobs, one seed each. Only diffuser.T and output prefix differ from the same input manifest. Source code confirms beta endpoints scale by 200/T and model time feature 1-t/T. More steps include more guidance updates. Keep all original geometry gates and original outputs. Inspect C-N, clashes, motif retention and material clearance before expanding replicates. No retries on failed execution.\n')
for j in jobs:
 out=Path(j['output']);out.mkdir(parents=True,exist_ok=False,mode=0o2770);out.chmod(0o2770)
 with (o/('T'+str(j['T'])+'.stdout')).open('x') as f:
  p=subprocess.Popen(j['command'],cwd=cwd,stdout=f,stderr=subprocess.STDOUT);log('step_diagnostic_started',T=j['T'],pid=p.pid,command=j['command']);rc=p.wait()
 log('step_diagnostic_exit',T=j['T'],exit_code=rc)
 if rc:raise SystemExit(rc)
 pdb=out/'design_0.pdb';assert pdb.is_file()
 log('step_diagnostic_technical_complete_NOT_GEOMETRY_PASS',T=j['T'],pdb=str(pdb),sha256=hashlib.sha256(pdb.read_bytes()).hexdigest())
log('step_diagnostic_batch_complete',cases=2)
