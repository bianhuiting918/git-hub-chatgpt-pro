"""Prepare all top50 donor candidates, no RFD execution or sequence design.
Rank by total existing material coverage; tie by identifier. Preserve all candidates.
"""
from pathlib import Path
import numpy as np,json,hashlib,datetime,collections
r=Path(__file__).resolve().parent;d=r/'route56_N10_20260915_v1';p=d/'route6_N10_single_residue_pilot_v1';out=d/'route6_top50_local_generation_v1'
panel=json.loads((p/'panel.json').read_text());cov=np.load(p/'coverage.npz')['coverage'];sites=json.loads((d/'N10_fixed_minimal_motif_material_v1/sites.json').read_text())
assert cov.shape==(len(panel),len(sites))==(892,1702)
order=sorted(range(len(panel)),key=lambda i:(-int(cov[i].sum()),panel[i]['id']))[:50]
assert collections.Counter(panel[i]['microstate'] for i in order)=={'ARG_1':49,'CYS_4':1}
out.mkdir(exist_ok=False);rows=[]
for rank,i in enumerate(order,1):
 v=panel[i];z=np.load(p/(v['microstate']+'.npz'));ix=np.flatnonzero(z['indices']==i);assert len(ix)==1;k=int(ix[0]);folder=out/v['id'];folder.mkdir()
 q=z['q'][k];assert np.isfinite(q).all()
 np.savez_compressed(folder/'candidate_exact.npz',q=q,r=z['r'],elements=z['elements'],donor=z['donor'][k],H=z['H'][k],donor_elements=z['donor_elements'][k],coverage=cov[i],panel_index=i)
 counts={mat:int(sum(bool(cov[i,j]) and s['material']==mat for j,s in enumerate(sites))) for mat in ['PA6','PA66']}
 row=dict(rank=rank,id=v['id'],microstate=v['microstate'],panel_index=i,total_covered=int(cov[i].sum()),material_counts=counts,site_indices=np.flatnonzero(cov[i]).tolist(),source_npz=str(p/(v['microstate']+'.npz')),source_sha256=hashlib.sha256((p/(v['microstate']+'.npz')).read_bytes()).hexdigest(),exact_input=str(folder/'candidate_exact.npz'),exact_input_sha256=hashlib.sha256((folder/'candidate_exact.npz').read_bytes()).hexdigest(),status='EXACT_INPUT_READY_RFD_NOT_SUBMITTED',limitations=['Existing small-core coverage, not full-enzyme coverage','Supported local fragment construction and peptide/steric audits still required','Cys donor chemistry retained separately; not silently substituted with Arg'])
 (folder/'manifest.json').write_text(json.dumps(row,indent=2));rows.append(row)
report=dict(status='50_EXACT_INPUTS_READY_0_NEW_RFD_SUBMITTED',total=50,rank_definition='Existing PA6+PA66 unique site counts, descending; identifier tie-break',source_panel_sha256=hashlib.sha256((p/'panel.json').read_bytes()).hexdigest(),source_coverage_sha256=hashlib.sha256((p/'coverage.npz').read_bytes()).hexdigest(),rows=rows)
(out/'manifest.json').write_text(json.dumps(report,indent=2))
(out/'RUNBOOK.md').write_text('Run prepare_top50_candidate_inputs_v1.py from project with existing CPU Python. Creates exclusive versioned exact-input packages for all top50 candidates. No GPU run is submitted. Next: build and audit supported local fragments without moving TDD/donor geometry, then RFD generation; report each candidate separately. Do not treat input preparation or template grafting as RFD generation. Preserve Cys chemistry and all failed candidates.\n')
with (r/'RUN_LOG.jsonl').open('a') as h:h.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__)),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),output=str(out),exit_code=0,prepared=50,submitted=0))+'\n')
print(json.dumps(dict(output=str(out),prepared=len(rows),submitted=0,top5=[{k:v[k] for k in ['id','total_covered','material_counts']} for v in rows[:5]]),indent=2))
