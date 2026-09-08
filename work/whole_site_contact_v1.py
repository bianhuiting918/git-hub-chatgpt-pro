#!/usr/bin/env python
"""Fixed whole-site subsets; paired material chemistry and enzyme contact faces."""
import sys,json,math,time,datetime,argparse,hashlib
from pathlib import Path
sys.dont_write_bytecode=True
R=Path(__file__).resolve().parent
sys.path.insert(0,str(R))
import numpy as np
from scipy.spatial import cKDTree
from scipy.ndimage import map_coordinates
import chemical_contact_surface_v3 as m
def target_count(n,fraction):return int(math.ceil(n*fraction-1e-10))
def whole_pass(values,tolerance=.4):return np.all(np.asarray(values)<=tolerance+1e-7,axis=1)
def subsets(fields,k,weight):
 fields=np.asarray(fields);n=len(fields);out={}
 for seed in range(n):
  chosen=[seed];mx=fields[seed].copy()
  while len(chosen)<k:
   scores=(np.maximum(np.maximum(fields,mx)+.5,0)*weight).sum(1)
   scores[chosen]=np.inf;j=int(np.argmin(scores));chosen.append(j);mx=np.maximum(mx,fields[j])
  key=tuple(sorted(chosen));score=float((np.maximum(mx+.5,0)*weight).sum())
  out[key]=score
 return [list(key) for key in sorted(out,key=lambda x:(out[x],x))]
class Environment:
 def __init__(self,mat):
  self.ce=m.ChemEnv(mat);z=self.ce.data;self.trees={}
  for e in ['C','N','O','S']:
   ix=np.where(z['elements']==e)[0]
   if len(ix):self.trees[e]=cKDTree(self.ce.env.wrap(z['xyz_A'][ix]),boxsize=[*z['box_A'][:2],100000.])
 def clear(self,points):
  w=self.ce.env.wrap(points);val=np.full(len(w),np.inf)
  for e,t in self.trees.items():val=np.minimum(val,t.query(w,workers=1)[0]-m.RAD[e])
  return val
def sample_mesh(v,f):
 edges=np.unique(np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1),axis=0)
 return np.vstack([v,v[edges].mean(1),v[f].mean(1)])
def chemistry(v,rows,frames,env):
 vals=[];weights=[];support=[]
 for s,F in zip(rows,frames):
  world=np.array(s['center'])+v@F.T;dist=env.clear(world)
  w=np.exp(-.5*(dist/1.5)**2);w[abs(dist)>4.5]=0
  vals.append(env.ce.chemistry(world));weights.append(w);support.append(abs(dist)<=1.5)
 vals=np.array(vals);w=np.array(weights);den=w.sum(0);mean=(vals*w[:,:,None]).sum(0)/np.maximum(den[:,None],1e-12)
 sd=np.sqrt(np.maximum((vals**2*w[:,:,None]).sum(0)/np.maximum(den[:,None],1e-12)-mean**2,0))
 sup=np.array(support).mean(0);valid=(sup>=.2)&(den>.01)
 return mean,sd,sup,valid
def payload_export(out,name,field,rows,frames,env,enzyme=False,q=None,rad=None):
 v,f,n=m.mesh_patch(field,enzyme=enzyme)
 if not len(f):return None,None,None
 props,sd,sup,valid=chemistry(v,rows,frames,env)
 fixed=(abs(m.core_clear(v,q,rad))<.5) if enzyme else None
 stats,pay=m.export_mesh(out,name,v,f,n,props,sd,sup,valid,enzyme=enzyme,fixedmask=fixed)
 # Independent material HBA/HBD channels, not enzyme-preference colors.
 extras=[]
 for key,col,target in [('HBA',1,np.array([.15,.55,.85])),('HBD',2,np.array([.7,.25,.7]))]:
  colors=.93+(target-.93)*props[:,col,None];colors[~valid]=.55
  extras.append(dict(pay,name=name+'_'+key,c=colors.round(4).tolist()))
 return stats,pay,extras
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--material',default='PET');ap.add_argument('--core',default='M329');ap.add_argument('--output',required=True);a=ap.parse_args()
 out=R/a.output;assert out.parent==R;out.mkdir(exist_ok=False)
 start=time.time();mat=a.material;cid=a.core
 def log(event,**kw):
  rec=dict(time=datetime.datetime.now().astimezone().isoformat(),event=event,core=cid,output=str(out),**kw)
  print(json.dumps(rec),flush=True)
  with (R/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(rec)+'\n')
 log('whole_site_start',script=__file__,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),argv=sys.argv)
 panel=json.load(open(R/'panel_surface_scan_v1/panel.json'));allrows=json.load(open(R/'panel_surface_scan_v1/sites.json'))
 p=next(x for x in panel if x['id']==cid);pi=panel.index(p)
 cov=np.load(R/'panel_surface_scan_v1/coverage.npz')['coverage']
 ids=[i for i,s in enumerate(allrows) if s['surface']==mat+'_original' and cov[pi,i]]
 rows=[allrows[i] for i in ids];N=len(ids)
 z=np.load(R/'chemical_contact_surface_full_v3'/mat/cid/'matched_material_fields.npz')
 assert ids==z['site_indices'].tolist()
 q=z['core_q_A'];don=z['donors_A'];rot=z['chemical_to_core_frame'];rad=np.array([m.RAD[e] for e in p['elements']])
 core=np.load(R/'panel_surface_scan_v1/cores.npz');assert np.allclose(q@rot.T,core[cid+'_q'])
 frames=[np.array(s['frame'])@rot for s in rows];env=Environment(mat)
 fields=[];core_clearance=[]
 for i,(s,F) in enumerate(zip(rows,frames)):
  world=np.array(s['center'])+m.GRID@F.T
  fields.append(-env.clear(world).reshape(37,37,37))
  core_clearance.append(float(np.min(env.clear(np.array(s['center'])+q@F.T)-rad)))
 fields=np.array(fields,np.float32);cc=np.array(core_clearance)
 assert np.all(cc>=-.400001),('ORIGINAL_CORE_CRITERION_MISMATCH',cc.tolist())
 np.savez_compressed(out/'material_fields_heavy_vdw.npz',fields=fields,core_q_A=q,donors_A=don,rot=rot,site_indices=ids,core_clearance_A=cc)
 log('original_core_criterion_revalidated',sites=N,min_clearance=float(cc.min()))
 inside=np.linalg.norm(m.GRID,axis=1)<=15
 # Retain local space; higher weight near the core, no planarity term.
 weight=np.exp(-np.linalg.norm(m.GRID[inside],axis=1)/10)
 short=fields.reshape(N,-1)[:,inside]
 result=dict(material=mat,core=cid,denominator=N,original_site_ids=[s['site_id'] for s in rows],core_coordinates_unchanged=True,core_clearance_A=cc.tolist(),tiers=[],search='all-seed greedy; top3 subsets per tier evaluated; no global optimality claim',collision='exact periodic heavy-atom VDW, original 0.4A overlap allowance; whole mesh vertices/edge-midpoints/centroids plus local occupied 1A volume samples and unchanged core',chemical_proxy=env.ce.audit,limitations=['one material snapshot','heuristic subset optimization','discrete mesh/volume checks, not continuous proof','backbone support and foldability NOT_EVALUATED','contact area measured but no validated energetic cutoff','coverage conditional on original core-eligible site denominator'])
 corepay=m.core_payload(q,p['elements'],don,p['template'],p['kind'])
 for frac in [.9,.7,.5]:
  k=target_count(N,frac);sets=subsets(short,k,weight);evaluated=[]
  for si,ix in enumerate(sets[:3]):
   fieldmax=fields[ix].max(0)
   pf=m.protein_field(fieldmax.ravel(),m.GRID,q,rad,.5).reshape(37,37,37)
   v,f,n=m.mesh_patch(pf,enzyme=True)
   if not len(f):continue
   samples=sample_mesh(v,f);vol=m.GRID[inside&(pf.ravel()>0)]
   va=m.vertex_areas(v,f);meshmax=[];volmax=[];contacts=[]
   for s,F in zip(rows,frames):
    origin=np.array(s['center'])
    distances=env.clear(origin+samples@F.T)
    meshmax.append(float(-distances.min()))
    volmax.append(float(-env.clear(origin+vol@F.T).min()) if len(vol) else 0)
    vd=distances[:len(v)]
    contacts.append(float(va[(vd>=-.4)&(vd<=1.5)].sum()/va.sum()))
   meshmax=np.array(meshmax);volmax=np.array(volmax)
   passed=(meshmax<=.4000001)&(volmax<=.4000001)&(cc>=-.4000001)
   # Any selected site failing disqualifies this subset; never swap IDs per spatial point.
   selected_pass=bool(passed[ix].all())
   rec=dict(selected_indices=ix,selected_site_ids=[rows[i]['site_id'] for i in ix],subset_all_pass=selected_pass,actual_coverage_n=int(passed.sum()),passed_site_ids=[rows[i]['site_id'] for i in np.flatnonzero(passed)],mesh_max_penetration_A=meshmax.tolist(),volume_max_penetration_A=volmax.tolist(),per_site_contact_area_fraction=contacts,selected_mean_contact_fraction=float(np.mean(np.array(contacts)[ix])),local_protein_volume_grid_A3=int(len(vol)),mesh_vertices=len(v),triangles=len(f))
   evaluated.append((rec,pf,fieldmax))
   log('whole_subset_evaluated',fraction=frac,subset_rank=si,required=k,coverage=int(passed.sum()),subset_pass=selected_pass,contact=rec['selected_mean_contact_fraction'])
  feasible=[x for x in evaluated if x[0]['subset_all_pass'] and x[0]['actual_coverage_n']>=k]
  candidates=feasible or evaluated
  if not candidates:
   result['tiers'].append(dict(fraction=frac,required_n=k,status='NO_MESH'));continue
  # Under fixed coverage requirement choose greatest shared contact, then retained volume.
  chosen=max(candidates,key=lambda x:(x[0]['subset_all_pass'],x[0]['selected_mean_contact_fraction'],x[0]['local_protein_volume_grid_A3']))
  rec,pf,fieldmax=chosen;ix=rec['selected_indices'];ss=[rows[i] for i in ix];ff=[frames[i] for i in ix]
  td=out/('coverage'+str(int(frac*100)));td.mkdir()
  rec.update(fraction=frac,required_n=k,status='GEOMETRIC_COVERAGE_PASS_NOT_FULL_ENZYME' if feasible else 'TARGET_NOT_MET_IN_SEARCH',candidate_subsets_searched=len(sets),candidate_subsets_validated=len(evaluated))
  meshes=[]
  st,pay,ex=payload_export(td,'enzyme_contact',pf,ss,ff,env,True,q,rad);meshes.append(pay);rec['enzyme_mesh']=st
  median=np.median(fields[ix],axis=0)
  st,pay,ex=payload_export(td,'material_consensus',median,ss,ff,env);meshes+=[pay]+ex;rec['material_consensus']=st
  medoid=min(ix,key=lambda j:float(np.mean(abs(short[j]-np.median(short[ix],axis=0)))))
  st,pay,ex=payload_export(td,'material_real_representative',fields[medoid],[rows[medoid]],[frames[medoid]],env);meshes+=[pay]+ex;rec['representative_site']=rows[medoid]['site_id']
  np.savez_compressed(td/'fitted_fields.npz',protein_field=pf,selected_material_max_field=fieldmax,material_median_field=median,selected_indices=ix)
  m.viewer(td/(mat+'_'+cid+'_coverage'+str(int(frac*100))+'_view.py'),meshes,corepay,'enzyme_contact')
  with (td/'summary.json').open('x') as f:json.dump(rec,f,indent=2)
  with (td/'search_candidates.json').open('x') as f:json.dump([x[0] for x in evaluated],f,indent=2)
  result['tiers'].append(rec)
 with (out/'summary.json').open('x') as f:json.dump(result,f,indent=2)
 lines=['# Whole-site coverage pilot: '+mat+' '+cid,'','Fixed-core denominator: '+str(N)+'. Not spatial quantiles. Material/enzymatic surfaces are paired on exactly the same selected site IDs.','',
 '|Target|Required|Actual whole-site geometry count|Selected subset all pass|Mean selected contact area fraction|Status|','|---|---:|---:|---|---:|---|']
 for t in result['tiers']:lines.append('|%d%%|%d|%d/%d|%s|%.3f|%s|'%(t['fraction']*100,t['required_n'],t['actual_coverage_n'],N,t['subset_all_pass'],t['selected_mean_contact_fraction'],t['status']))
 lines+=['','View files: each coverage folder contains one self-contained *_view.py. Default enzyme_contact. Disable it and enable material_consensus or material_real_representative to view material geometry/lipophilicity; *_HBA/HBD independently show material donor/acceptor channels. Cyan dot is carbonyl C, not exposed surface.','',
 'Geometry uses heavy-atom VDW union to match original core screening; old q-series used SES and is not numerically directly comparable. Existing 0.4A allowance unchanged. 0.5A stand-off, 1A grid, 15A local patch. No new MD or donor sampling. Shape search optimizes whole subsets; no imposed flatness or planarity objective. Global optimum not proven.','',
 'PASS is geometric only. Contact area is reported separately, not declared favorable binding. No backbone, scaffold support, folding, dynamic robustness or energetic validation. Material consensus averages can hide patch variability; actual medoid is provided. All results are conditional on this core original eligible sites, not 90/70/50 percent of all PET sites.']
 with (out/'README.md').open('x') as f:f.write('\n'.join(lines))
 log('whole_site_complete',seconds=time.time()-start,exit_status=0)
 with (R/'RUNBOOK.md').open('a') as f:f.write('\nWhole-site pilot: '+str(out)+'; rerun whole_site_contact_v1.py --material '+mat+' --core '+cid+' --output NEW_NAME with existing CPU env. Read per-tier whole-site coverage, contact areas and scope.\\n')
if __name__=='__main__':main()
