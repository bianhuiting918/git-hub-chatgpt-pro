"""Audit bin-specific route6 manifests and separate single-site from multisite inputs."""
from pathlib import Path
import json,hashlib,datetime,math,collections
import numpy as np
R=Path(__file__).resolve().parent
D=R/'route56_N10_20260915_v1'
B=D/'route6_top50_local_generation_v1/candidate_specific_binned_shapes_v1'
sites=json.loads((D/'N10_fixed_minimal_motif_material_v1/sites.json').read_text())
panel=json.loads((D/'route6_N10_single_residue_pilot_v1/panel.json').read_text())
coverage=np.load(D/'route6_N10_single_residue_pilot_v1/coverage.npz')['coverage']
pidx={p['id']:i for i,p in enumerate(panel)}
rows=[];errors=[];pending=[]
for p in sorted(B.glob('*/*/manifest.json')):
 try:r=json.loads(p.read_text())
 except json.JSONDecodeError:
  pending.append(str(p));continue
 ids=r['selected_site_indices'];mat=r['material'];bn=r['exposure_bin'];i=pidx[r['core']]
 den=sum(s['material']==mat and s['bin']==bn and bool(coverage[i,j]) for j,s in enumerate(sites))
 group=sum(s['material']==mat and s['bin']==bn for s in sites)
 ok=(den==r['conditional_denominator'] and group==r['all_exposure_group_denominator']
     and len(ids)==math.ceil(den*r['fraction']-1e-9) and len(set(ids))==len(ids)
     and len(ids)==r['selected_site_count'] and (len(ids)>=2)==r['multi_site']
     and all(sites[j]['material']==mat and sites[j]['bin']==bn and coverage[i,j] for j in ids))
 if not ok:errors.append([r['name'],'subset_or_denominator'])
 if r['status']=='PREPARED_NOT_SUBMITTED':
  for key in ['pdb','shape']:
   if hashlib.sha256(Path(r[key]).read_bytes()).hexdigest()!=r[key+'_sha256']:
    errors.append([r['name'],key+'_hash'])
  assert r['command'][0]=='/usr/local/bin/gpurun'
 rows.append(dict(name=r['name'],material=mat,core=r['core'],bin=bn,fraction=r['fraction'],
                  selected_n=len(ids),conditional_denominator=den,group_denominator=group,
                  multi_site=len(ids)>=2,status=r['status'],manifest=str(p)))
complete=(B/'summary.json').exists()
if complete:
 summary=json.loads((B/'summary.json').read_text())
 if summary['total']!=len(rows):errors.append(['TOTAL','summary_manifest_count_mismatch'])
 if {j['name'] for j in summary['jobs']}!={r['name'] for r in rows}:errors.append(['TOTAL','summary_names_mismatch'])
counts=collections.Counter((r['material'],r['fraction'],r['status'],r['multi_site']) for r in rows)
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%S%fZ')
result=dict(status='AUDIT_FAIL' if errors else 'COMPLETE_INPUT_BOOKKEEPING_AUDIT' if complete else 'SNAPSHOT_INPUT_BOOKKEEPING_AUDIT',
            rows=rows,total=len(rows),preparation_complete=complete,errors=errors,pending=pending,
            groups=[dict(material=k[0],fraction=k[1],status=k[2],multi_site=k[3],count=v) for k,v in sorted(counts.items())],
            ready_multisite=[r['manifest'] for r in rows if r['multi_site'] and r['status']=='PREPARED_NOT_SUBMITTED'],
            ready_single_site=[r['manifest'] for r in rows if not r['multi_site'] and r['status']=='PREPARED_NOT_SUBMITTED'],
            limitations=['Input bookkeeping and file hashes only; not new continuous mesh or atomistic steric validation',
                         'No GPU submission, backbone generation, independent folding, or activity pass',
                         'Do not count alternative subsets as independent proteins or polymer sites'])
out=B/('input_audit_'+stamp+'.json')
with out.open('x') as f:json.dump(result,f,indent=2)
with (R/'RUN_LOG.jsonl').open('a') as f:
 f.write(json.dumps(dict(time=stamp,script=str(Path(__file__)),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                        output=str(out),status=result['status'],exit_code=1 if errors else 0,total=len(rows)))+'\n')
print(json.dumps(dict(output=str(out),status=result['status'],total=len(rows),ready_multisite=len(result['ready_multisite']),
                      ready_single_site=len(result['ready_single_site']),errors=errors)))
raise SystemExit(1 if errors else 0)
