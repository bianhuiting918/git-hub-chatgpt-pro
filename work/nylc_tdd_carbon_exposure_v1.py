#!/usr/bin/env python3
"""Carbon solvent-center shell diagnostic. Not SES, enzyme clearance or a hard gate."""
from pathlib import Path
import sys,json,datetime,hashlib
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parent/'surface_previews/axis_density_sixface_all_external_trim_v6'))
from density_cut_variants import parse_pdb
from terminal_density_trim_3d import center_whole_chains_in_primary_cell
RAD={'C':1.7,'O':1.52,'N':1.55,'H':1.2}
PROBE=1.4;N=1024
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def shell(q,xyz,r,tree):
 points=q+3.1*DIR
 near=tree.query_ball_point(q,6.21)
 d2=((points[:,None,:]-xyz[near][None,:,:])**2).sum(2)
 return np.all(d2>=(r[near][None,:]+PROBE)**2-1e-8,axis=1)
i=np.arange(N);z=1-2*(i+.5)/N;phi=i*np.pi*(3-np.sqrt(5))
DIR=np.column_stack((np.sqrt(1-z*z)*np.cos(phi),np.sqrt(1-z*z)*np.sin(phi),z))
assert np.max(abs(np.linalg.norm(DIR,axis=1)-1))<1e-12
assert shell(np.zeros(3),np.zeros((1,3)),np.array([1.7]),cKDTree(np.zeros((1,3)))).all()
assert not shell(np.zeros(3),np.array([[0.,0.,0.]]),np.array([2.0]),cKDTree(np.zeros((1,3)))).any()
def run(mat,dp,out):
 inp=ROOT/'inputs'/f'{mat}_v10_ge90_uncapped.pdb'
 audit=json.loads((ROOT/'input_audit_v2'/f'{mat}_intact_amides.json').read_text())
 assert sha(inp)==audit['input_sha256']
 lines=[l for l in inp.read_text().splitlines() if l.startswith(('ATOM  ','HETATM'))]
 xyz=np.array([[float(l[k:k+8]) for k in (30,38,46)] for l in lines]);r=np.array([RAD[l[76:78].strip()] for l in lines])
 parent=ROOT.parent/f'{mat.lower()}_dp{dp}_400chain/water_slab_v1/dry_export/{mat}_DP{dp}_400chain_water_equilibrated_dry.pdb'
 assert sha(parent)==audit['parent_sha256']
 _,box,pa,_=parse_pdb(parent);center_whole_chains_in_primary_cell(pa,box)
 px=np.array([a['xyz'] for a in pa]);pr=np.array([RAD[a['element']] for a in pa])
 # Same finite centered-parent coordinates: no periodic copy changes in this diagnostic.
 tree=cKDTree(xyz);ptree=cKDTree(px)
 candidates=[a for a in json.loads((ROOT/'surface_density_v2'/f'{mat}_sites.json').read_text()) if a['within_3A_density_boundary'] and a['inward_support']>=.7]
 rows=[];tm=[];pm=[]
 for a in candidates:
  q=np.array(a['carbonyl_xyz_A']);t=shell(q,xyz,r,tree);p=shell(q,px,pr,ptree)
  assert not np.any(p&~t),'Parent accessibility cannot exceed subset'
  sign=1 if a['nearest_side']=='top' else -1
  outward=DIR[:,2]*sign>0
  b=dict(a,trimmed_clear_points=int(t.sum()),parent_clear_points=int(p.sum()),trimmed_outward_points=int((t&outward).sum()),parent_outward_points=int((p&outward).sum()),deletion_unmasked_points=int((t&~p).sum()),shell_area_proxy_A2=float(t.sum()/N*4*np.pi*3.1**2))
  rows.append(b);tm.append(t);pm.append(p)
 np.savez_compressed(out/f'{mat}_shell_masks.npz',directions=DIR,trimmed=np.array(tm),parent=np.array(pm),site_ids=np.array([r['site_id'] for r in rows]))
 (out/f'{mat}_sites.json').write_text(json.dumps(rows,indent=2))
 return dict(material=mat,candidates=len(rows),any_clear_trimmed=sum(r['trimmed_clear_points']>0 for r in rows),any_clear_parent=sum(r['parent_clear_points']>0 for r in rows),any_outward_trimmed=sum(r['trimmed_outward_points']>0 for r in rows),any_outward_parent=sum(r['parent_outward_points']>0 for r in rows),any_deletion_unmasked=sum(r['deletion_unmasked_points']>0 for r in rows),status='SHELL_CLEARANCE_DIAGNOSTIC_NOT_ENZYME_ACCESS',samples=N,probe_radius_A=PROBE,carbon_shell_radius_A=3.1,boundary='finite centered parent comparison; no solvent connectivity or ingress proof; not a hard exclusion')
if __name__=='__main__':
 out=ROOT/'carbon_exposure_v1';out.mkdir(exist_ok=False);rc=1
 try:
  summary=[run('PA6',16,out),run('PA66',8,out)]
  (out/'summary.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2));rc=0
 finally:
  with (ROOT/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__).resolve()),script_sha256=sha(Path(__file__)),command='CPU Python carbon_exposure_v1.py',inputs='v10 and original dry parent PDBs; surface_density_v2 sites',outputs=str(out),parameters=dict(samples=N,probe_A=PROBE,cohort_support=.7),exit_code=rc))+'\n')
