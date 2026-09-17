"""CPU audit adapter for two repaired exploratory targets; no predictions or GPU actions."""
from pathlib import Path
import json,hashlib,sys
R=Path(__file__).resolve().parent
B=R/'route56_N10_20260915_v1/route6_top50_local_generation_v1/two_repaired_exploratory_mpnn20_v1'
source=R/'audit_supported_arg_core4_predictions_v1.py'
assert hashlib.sha256(source.read_bytes()).hexdigest()=='292f826b7540433c440c351f2559b348ecda94e1f98e347c38b723367af54a63'
jobs=json.loads((B/'manifest.json').read_text())['jobs']
assert len(jobs)==2
for job in jobs:
 cid=job['name'].split('_')[0]
 ref=Path(job['source']);assert hashlib.sha256(ref.read_bytes()).hexdigest()==job['source_sha256']
 s=source.read_text()
 replacements={
 "r/'route6_supported_arg_core4_af20_v1'":"r/"+repr('route6_repaired_'+cid+'_af20_v1'),
 "'n9_core_donor_design_20260910_gpu_outputs_v1/af2_route6_supported_arg_core4_20_v1'":repr('n9_core_donor_design_20260910_gpu_outputs_v1/af2_route6_repaired_'+cid+'_20_v1'),
 "r/'route56_N10_20260915_v1/route6_supported_ARG_SC4559268_v1/fixed_core_repair_v1/after.pdb'":"Path("+repr(str(ref))+")"}
 for a,b in replacements.items():
  assert s.count(a)==1,(cid,a)
  s=s.replace(a,b)
 compile(s,str(source),'exec')
 mapping=R/('route6_repaired_'+cid+'_af20_v1')/'sequence_mapping.json'
 if '--check' in sys.argv:
  print(json.dumps({'candidate':cid,'reference_sha256':job['source_sha256'],'derived_audit_sha256':hashlib.sha256(s.encode()).hexdigest(),'mapping_ready':mapping.exists(),'status':'CONFIGURATION_VALIDATED_NOT_PREDICTION_AUDIT'}))
  continue
 if not mapping.exists():
  print(json.dumps({'candidate':cid,'status':'NOT_EVALUATED_SEQUENCE_MAPPING_PENDING'}));continue
 queries=json.loads(mapping.read_text());assert len(queries)==20
 for name,seq in queries.items():
  assert name.startswith('R6_REP_'+cid+'_') and len(seq)==126
  assert all(seq[i-1]==a for i,a in [(1,'T'),(40,'D'),(42,'D'),(63,'R')])
 out=R.parent/('n9_core_donor_design_20260910_gpu_outputs_v1/af2_route6_repaired_'+cid+'_20_v1')
 if not any((out/(name+'.done.txt')).exists() for name in queries):
  print(json.dumps({'candidate':cid,'status':'NOT_EVALUATED_PREDICTIONS_PENDING'}));continue
 exec(compile(s,str(source),'exec'),{'__file__':str(source),'__name__':'__main__'})
