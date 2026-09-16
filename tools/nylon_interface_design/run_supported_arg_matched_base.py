"""One matched Base-checkpoint diagnostic. Run from the remote project.
Keeps the supported motif, CAD shape, seed and T50 unchanged; swaps checkpoint.
Never overwrites results and never retries or submits another tier automatically.
"""
from pathlib import Path
import json,subprocess,datetime,hashlib
r=Path.cwd().resolve();d=r/'route56_N10_20260915_v1'
source=d/'route6_supported_ARG_SC4559268_v1/PA6_SC4559268_coverage10_variant1/manifest.json'
j=json.loads(source.read_text())
base=next(v for v in json.loads((d/'generator_checkpoint_controls126_v1.json').read_text()) if v['name']=='base')
checkpoint=next(a for a in base['command'] if a.startswith('inference.ckpt_override_path='))
assert Path(checkpoint.split('=',1)[1]).is_file()
for key in ['pdb','shape']:
 assert hashlib.sha256(Path(j[key]).read_bytes()).hexdigest()==j[key+'_sha256']
o=d/'route6_supported_ARG_SC4559268_base_v1';o.mkdir(exist_ok=False)
g=r.parent/'n9_core_donor_design_20260910_gpu_outputs_v1/route6_supported_ARG_SC4559268_base_v1';g.mkdir(mode=0o2770,exist_ok=False);g.chmod(0o2770)
cmd=[checkpoint if a.startswith('inference.ckpt_override_path=') else 'inference.output_prefix='+str(g/'design') if a.startswith('inference.output_prefix=') else a for a in j['command']]
assert cmd[0]=='/usr/local/bin/gpurun' and 'diffuser.T=50' in cmd
assert len([(a,b) for a,b in zip(cmd,j['command']) if a!=b])==2
assert not any('CUDA_MPS_' in a or a=='sudo' for a in cmd)
rec={'command':cmd,'source_manifest':str(source),'source_manifest_sha256':hashlib.sha256(source.read_bytes()).hexdigest(),'pdb_sha256':j['pdb_sha256'],'shape_sha256':j['shape_sha256'],'output':str(g),'scope':'one matched checkpoint diagnostic, not a production batch'}
(o/'manifest.json').write_text(json.dumps(rec,indent=2))
(o/'RUNBOOK.md').write_text('Run this script with the existing CPU Python from the remote project. It launches one gpurun A100 job, changing only the checkpoint and output path. Exclusive-create outputs; no automatic retry. Read completion.json then independently compare motif retention, all backbone connections, self-clashes and actual material access against the ActiveSite run. Exit 0 is technical completion only.\n')
def log(event,**kw):
 row={'time':datetime.datetime.now(datetime.timezone.utc).isoformat(),'event':event,'script':str(Path(__file__).resolve()),'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'output':str(g),**kw}
 for path in [o/'events.jsonl',r/'RUN_LOG.jsonl']:
  with path.open('a') as h:h.write(json.dumps(row)+'\n')
 print(json.dumps(row),flush=True)
with (o/'generation.log').open('x') as h:
 p=subprocess.Popen(cmd,cwd='/data/bht2/shape_conditioned_rfdiffusion_20260905_v1/source/RFdiffusion',stdout=h,stderr=subprocess.STDOUT)
 (o/'launch.json').write_text(json.dumps({'wrapper_pid':p.pid,'command':cmd},indent=2))
 log('supported_ARG_matched_Base_start',wrapper_pid=p.pid)
 rc=p.wait()
(o/'completion.json').write_text(json.dumps({'exit_code':rc,'output_exists':(g/'design_0.pdb').is_file(),'scientific_status':'NOT_EVALUATED'}))
log('supported_ARG_matched_Base_exit',exit_code=rc)
raise SystemExit(rc)
