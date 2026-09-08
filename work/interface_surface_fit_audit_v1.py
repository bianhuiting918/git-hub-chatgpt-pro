#!/usr/bin/env python3
"""Replay old-coverage invariants and independently verify every fitted field."""
from pathlib import Path
import json,datetime,hashlib
import numpy as np
R=Path(__file__).resolve().parent;B=R.parent;P=R/'panel_surface_scan_v1';F=R/'surface_fit_robustness_v1'
panel=json.loads((P/'panel.json').read_text());sites=json.loads((P/'sites.json').read_text());coverage=np.load(P/'coverage.npz')['coverage'];rows=json.loads((F/'results.json').read_text());grid=np.load(P/'grid.npz')['points_A']
lookup={w['id']:j for j,w in enumerate(panel)}
oldchecks=0
root=B/'single_residue_dual_donor_20260908_v1/fixed_coupled_transfer_v1'
for co in json.loads((root/'summary.json').read_text())['cohorts']:
 mask=np.array([((co['material']=='PET' and x['surface']=='PET_original') or x['surface']==co['material']+'_trimmed' or co['material']=='NYLON_COMBINED' and x['surface'] in ['PA6_trimmed','PA66_trimmed']) and x['support']>=co['min_support'] for x in sites])
 for w in co['top10']+[v['best_balanced'] for v in co['per_microstate'] if 'best_balanced' in v]:
  if w['candidate_id'] in lookup:
   actual=int(coverage[lookup[w['candidate_id']],mask].sum());assert actual==w['coverage'],(co['name'],w['candidate_id'],actual,w['coverage']);oldchecks+=1
pet=B/'pet_dp10_400chain_direct_v1/donor_geometry_batch_singleNH_v1/joint_transfer_map_v1'
pm=json.loads((pet/'motifs.json').read_text());pz=np.load(pet/'joint_transfer_results.npz')
newsite={(x['surface'],x['site_id']):i for i,x in enumerate(sites)}
petcols=[newsite[('PET_original',str(s))] for s in pz['site_ids']]
for j,w in enumerate(pm):
 if w['id'] in lookup:assert np.array_equal(coverage[lookup[w['id']],petcols],pz['coverage'][j]);oldchecks+=1
ny=B/'rectangular_400chain_v1/nylc_tdd_surface_v1'
nm=json.loads((ny/'fixed_core_library_v1/motifs.json').read_text());nz=np.load(ny/'fixed_core_transfer_v1/coverage.npz');ns=json.loads((ny/'fixed_core_transfer_v1/sites.json').read_text())
cols=[newsite[(x['material']+'_trimmed',x['site_id'])] for x in ns]
for j,w in enumerate(nm):
 if w['id'] in lookup:assert np.array_equal(coverage[lookup[w['id']],cols],nz['coverage'][j]);oldchecks+=1
cache={i:dict(np.load(P/'fields'/(str(i)+'.npz'))) for i in np.flatnonzero(coverage.any(0))}
maps=0
for row in rows:
 if not row['sites']:continue
 z=np.load(F/'maps'/(row['fit_id']+'.npz'));ids=z['site_indices'];assert len(ids)==row['sites'] and len(set(ids.tolist()))==len(ids)
 x=np.array([cache[int(i)]['polymer_vdw_clearance_A'] for i in ids]);expected=(x<1.3).mean(0);minimum=x.min(0)
 assert np.allclose(z['occupancy_probability'],expected,atol=1e-7)
 assert np.array_equal(z['strict_common_near_interface_mask'],(minimum>=1.3)&(minimum<=4.3))
 assert int(z['strict_common_near_interface_mask'].sum())==row['strict_common_near_interface_voxels']
 if row['surface'].endswith('trimmed'):assert np.isnan(z['mean_partial_charge_annotation_e']).all()
 else:assert np.isfinite(z['mean_partial_charge_annotation_e']).all()
 if 'occupancy_bootstrap_q05' in z:
  assert (z['occupancy_bootstrap_q05']>=0).all() and (z['occupancy_bootstrap_q95']<=1).all()
  assert (z['occupancy_bootstrap_q05']<=z['occupancy_bootstrap_q95']+1e-7).all()
 for fold in row['cv_details']:assert fold['training_groups']>=3 and fold['allowed_voxels']>=10 and 0<=fold['heldout_clash_fraction_mean']<=1
 maps+=1
summary=dict(status='PASS_REPLAY_AND_MAP_AUDIT',historical_coverage_checks=oldchecks,nonempty_maps_checked=maps,withheld_trimmed_charge='PASS_NAN',scope='does not validate dielectric electrostatics or activity; spatial holdout conditional on selected fixed cores')
(F/'FINAL_AUDIT.json').write_text(json.dumps(summary,indent=2))
with (R/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__)),sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),exit_code=0,summary=summary))+'\n')
print(json.dumps(summary),flush=True)
