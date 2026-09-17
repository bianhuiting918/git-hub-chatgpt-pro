"""Prepare a paired zero-guide diagnostic without submitting GPU work."""
from pathlib import Path
import json,hashlib,datetime
R=Path(__file__).resolve().parent
B=R/'route56_N10_20260915_v1/route6_top50_local_generation_v1'
O=B/'SC5150928_zero_guide_diagnostic_v1'
name='PA6_SC5150928_coverage10_variant1'
events=[json.loads(x) for x in (B/'first19_RFD3_v2/events.jsonl').read_text().splitlines()]
matches=[x for x in events if x.get('step')=='candidate_rfd_started' and x.get('name')==name]
assert len(matches)==1
base=matches[0]['command']
assert base[0]=='/usr/local/bin/gpurun'
for token in ['diffuser.T=50','inference.num_designs=3','inference.deterministic=true','potentials.guide_scale=2']:
 assert token in base,token
m=json.loads((B/'candidate_specific_shapes_v2/SC5150928'/name/'manifest.json').read_text())
for key in ['pdb','shape']:
 assert hashlib.sha256(Path(m[key]).read_bytes()).hexdigest()==m[key+'_sha256']
prefix=next(x.split('=',1)[1] for x in base if x.startswith('inference.output_prefix='))
baseline=[]
for seed in range(3):
 p=Path(prefix+'_'+str(seed)+'.pdb')
 assert p.is_file()
 baseline.append(dict(seed=seed,path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
newprefix=R.parent/'n9_core_donor_design_20260910_gpu_outputs_v1/route6_SC5150928_zero_guide_diagnostic_v1/design'
command=[('potentials.guide_scale=0' if x=='potentials.guide_scale=2' else 'inference.output_prefix='+str(newprefix) if x.startswith('inference.output_prefix=') else x) for x in base]
assert sum(x!=y for x,y in zip(base,command))==2
assert not newprefix.parent.exists(),'Do not reuse an existing output directory'
result=dict(status='PREPARED_NOT_SUBMITTED',purpose='Paired guide_scale 2 versus 0, same shape module and sampling retained; diagnose backbone connectivity, not production shape success',baseline=baseline,baseline_command=base,command=command,cwd='/data/bht2/shape_conditioned_rfdiffusion_20260905_v1/source/RFdiffusion',input_pdb=m['pdb'],input_pdb_sha256=m['pdb_sha256'],shape=m['shape'],shape_sha256=m['shape_sha256'],gates='Unchanged CN 1.2-1.5 A, omega, backbone clashes, fixed motif RMSD and material replay; no sidechain/folding/activity claim',limitations=['One core and three seeds only','Same seeds do not alone prove identical starting noise; compare available initial trajectory/state before claiming paired-noise equivalence','Pooled legacy shape is used solely for diagnosis; not a bin-specific production design'],stop='No automatic GPU retries; do not launch while two project GPU tasks are active')
O.mkdir(exist_ok=False)
(O/'manifest.json').write_text(json.dumps(result,indent=2))
(O/'RUNBOOK.md').write_text('Preparation only. Run this preparation script from the N9 project with the existing CPU Python. The manifest contains the exact gpurun argument array and cwd for a later explicitly scheduled diagnostic. No GPU work was launched. Preserve all baseline outputs and all gates. Evaluate all three seeds, including failures. Verify initial-state equivalence separately; fixed seed alone is insufficient. Stop on GPU failure under project MPS rules. No source or production configuration was modified.\n')
with (R/'RUN_LOG.jsonl').open('a') as f:
 f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__)),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),output=str(O),status=result['status'],exit_code=0))+'\n')
print(json.dumps(dict(output=str(O),status=result['status'],baseline_count=len(baseline),gpu_submitted=False)))
