"""Eight exposure-specific two-residue donor fits at requested 10/20/30%; CPU only."""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
from pathlib import Path
import json,hashlib,sys,datetime
r=Path(__file__).resolve().parent;d=r/'route56_N10_20260915_v1'
source=r/'fit_route6_complete_donor_whole_site_v1.py'
assert hashlib.sha256(source.read_bytes()).hexdigest()=='f5fbda548821b6fbaebb729ece14903e2c50cee8ebcd093cf323f9abc5626bc1'
text=source.read_text()
assert text.count('for frac in [.9,.7,.5]:')==1
text=text.replace('for frac in [.9,.7,.5]:','for frac in [.1,.2,.3]:')
text=text.replace("rec.update(fraction=frac,required_n=k,","rec.update(multi_site=(k>=2),selected_site_count=k,fraction=frac,required_n=k,")
text=text.replace('not 90/70/50 percent of all PET sites','not 10/20/30 percent of all material sites')
text=text.replace('rerun new_exposure_whole_site_fit_20260909_v1.py','protocol wrapper fit_two_donor_bins_lowcov_v1.py; underlying source fit_route6_complete_donor_whole_site_v1.py')
compile(text,str(source),'exec')
summary_path=d/'route6_complete_donor_material_v1/summary.json'
summary=json.loads(summary_path.read_text())
jobs=[]
for x in summary['reports']:
 if x['bin']<0:continue
 top=x['top10'][0]
 jobs.append(dict(material=x['material'],core=top['id'],exposure_bin=x['bin'],original_eligible_n=top['sites'],all_bin_n=x['denominator'],output='route6_complete_'+x['material']+'_'+top['id']+'_bin'+str(x['bin'])+'_lowcov_v1'))
assert len(jobs)==8 and len({x['output'] for x in jobs})==8
if '--check' in sys.argv:
 print(json.dumps(dict(status='READ_ONLY_CHECK',jobs=jobs,fractions=[.1,.2,.3],CPU_only=True)))
 raise SystemExit(0)
out=d/'route6_two_donor_binned_lowcov_v1';out.mkdir(exist_ok=False)
manifest=dict(status='CPU_FITTING_RUNNING_NOT_RFD',jobs=jobs,fractions=[.1,.2,.3],selection='Top complete-fragment material-eligible count within each original bin; existing ID tie order; not global optimality',summary_sha256=hashlib.sha256(summary_path.read_bytes()).hexdigest(),source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),derived_script_sha256=hashlib.sha256(text.encode()).hexdigest())
with (out/'manifest.json').open('x') as f:json.dump(manifest,f,indent=2)
with (out/'RUNBOOK.md').open('x') as f:f.write('Run fit_two_donor_bins_lowcov_v1.py with existing project CPU Python. --check is read-only. Uses eight top complete-donor cores, one per material/bin. Original geometric criteria, fields and chemistry retained; only requested coverage tiers change to 10/20/30. Original results preserved. Each tier records single/multi-site distinction. No GPU jobs submitted. Stops at first error; inspect completed output and logs before any manual recovery. No folding or activity conclusion.\n')
def log(step,**kw):
 row=dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),step=step,**kw)
 print(json.dumps(row),flush=True)
 for p in [out/'events.jsonl',r/'RUN_LOG.jsonl']:
  with p.open('a') as f:f.write(json.dumps(row)+'\n')
for j in jobs:
 log('two_donor_bin_fit_start',**j)
 sys.argv=[str(source),'--material',j['material'],'--core',j['core'],'--bin',str(j['exposure_bin']),'--output',j['output']]
 try:
  exec(compile(text,str(source),'exec'),{'__file__':str(source),'__name__':'__main__'})
  a=json.loads((r.parent/'interface_surface_robustness_20260908_v1'/j['output']/'summary.json').read_text())
  assert a['denominator']==j['original_eligible_n'] and a['exposure_bin']==j['exposure_bin'] and a['exposure_group_denominator']==j['all_bin_n']
  assert [x['fraction'] for x in a['tiers']]==[.1,.2,.3]
  log('two_donor_bin_fit_technical_complete',core=j['core'],material=j['material'],bin=j['exposure_bin'],tiers=[dict(fraction=x['fraction'],required=x['required_n'],status=x['status']) for x in a['tiers']])
 except BaseException as exc:
  log('two_donor_bin_fit_failed_no_retry',core=j['core'],error_type=type(exc).__name__,error=str(exc));raise
log('two_donor_eight_bin_fit_complete_NOT_RFD')
