#!/usr/bin/env python3
"""Single coupled amino-acid donor surface scan; CPU, resumable, finite library."""
from pathlib import Path
import ast,json,hashlib,datetime,sys,time,concurrent.futures
import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
from scipy.stats import qmc
R=Path(__file__).resolve().parent
P=R.parent/'pet_dp10_400chain_direct_v1'
N=R.parent/'rectangular_400chain_v1/nylc_tdd_surface_v1'
OUT=R/'surface_scan_v1'
RAD={'C':1.7,'N':1.55,'O':1.52,'S':1.8}
CFG={'n_triad_trials':12000,'vdw_overlap_allowance_A':.4}
BOX=np.array([125.697,125.697,182.848])
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
source=P/'pet_donor_geometry_batch_singleNH_v1.py'
assert sha(source)=='b796be0506a0d696ee19aad9288e2b71b6c81131c30f7e4bf2deb00c94ddb80a'
names={'norm','load_pdb','trees','clear','make_triads','seed_for'}
nodes=[x for x in ast.parse(source.read_text()).body if isinstance(x,ast.FunctionDef) and x.name in names]
assert {n.name for n in nodes}==names
exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'),globals())
def chemical(q,h,o,els):
 dist=np.linalg.norm(q-o,axis=-1)
 lo=np.where(els=='S',3.1,2.7);hi=lo+.5
 a=q-h;b=o-h
 cosine=np.sum(a*b,axis=-1)/(np.linalg.norm(a,axis=-1)*np.linalg.norm(b,axis=-1))
 return ((dist>=lo-1e-8)&(dist<=hi+1e-8)&(cosine<=np.cos(np.deg2rad(140))+1e-12)).all(-1)
def build_clouds():
 folder=OUT/'clouds';folder.mkdir(parents=True,exist_ok=True)
 clouds={}
 for row in json.loads((R/'coupled_torsion_pilot_v1/summary.json').read_text())['results']:
  key=row['key']
  if row['common_oxygen_samples']==0:continue
  meta=json.loads((R/'conformers_v1'/(key+'.json')).read_text())
  raw=np.load(R/'conformers_v1'/(key+'.npz'));z=np.load(R/'coupled_torsion_pilot_v1'/(key+'.npz'))
  path=folder/(key+'.npz')
  if not path.exists():
   gids=np.unique(z['geometry_indices']);u=qmc.Sobol(2,scramble=True,seed=seed_for(key+'|surface')).random_base2(13)
   chosen=[]
   for uv in u:
    g=gids[min(int(uv[0]*len(gids)),len(gids)-1)]
    oi=np.flatnonzero(z['geometry_indices']==g)
    chosen.append(int(oi[min(int(uv[1]*len(oi)),len(oi)-1)]))
   chosen=np.array(chosen);gi=z['geometry_indices'][chosen];pi=z['donor_pair_indices'][chosen]
   x=z['residue_A'][gi];o=z['oxygen_A'][chosen]
   dn=np.zeros((len(x),2,3));hh=np.zeros_like(dn);els=np.empty((len(x),2),dtype='U1')
   for i in range(len(x)):
    for j in range(2):
     d=meta['donors'][pi[i,j]];dn[i,j]=x[i,d['atom']];els[i,j]=d['element']
     h=x[i,d['hydrogens']];a=dn[i,j]-h;b=o[i]-h
     c=np.sum(a*b,axis=1)/(np.linalg.norm(a,axis=1)*np.linalg.norm(b,axis=1))
     hh[i,j]=h[np.argmin(c)]
   rot=Rotation.random(len(x),random_state=np.random.default_rng(seed_for(key+'|rotation'))).as_matrix()
   def tr(v):return np.einsum('naj,nkj->nak',v-o[:,None],rot)+np.array([1.23,0,0])
   heavy=raw['heavy_indices'];q=tr(x[:,heavy]);dn=tr(dn);hh=tr(hh)
   assert chemical(dn,hh,np.array([1.23,0,0]),els).all()
   assert np.max(np.linalg.norm(q,axis=-1))<26,'30 A environment cutoff unsafe'
   np.savez_compressed(path,q_A=q,donor_A=dn,H_A=hh,elements=raw['elements'][heavy],donor_elements=els,pair_indices=pi,intrinsic_sample_indices=chosen)
  z=np.load(path);clouds[key]={k:z[k] for k in z.files}
 return clouds
def get_inputs():
 rows=[];envs={};lookup={}
 for mat in ['PET','PA6','PA66']:
  path=P/'pet_dp10_400chain/water_slab_v1/dry_export/PET_DP10_400chain_water_equilibrated_dry.pdb' if mat=='PET' else N/'inputs'/(mat+'_v10_ge90_uncapped.pdb')
  rec=load_pdb(path);xyz=np.array([q for _,q,_ in rec]);el=np.array([e for _,_,e in rec])
  if mat=='PET':xyz[:,:2]%=BOX[:2]
  envs[mat]=(xyz,el)
  lookup[mat]={(l[72:76].strip(),l[12:16].strip()) if mat=='PET' else l[6:11]:i for i,(l,_,_) in enumerate(rec)}
  rr=json.loads((P/'donor_geometry_batch_singleNH_v1/extension70_v1/eligible_sites.json' if mat=='PET' else N/'carbon_exposure_v1'/(mat+'_sites.json')).read_text())
  for row in rr:rows.append(dict(material=mat,row=row))
 assert [sum(r['material']==m for r in rows) for m in ['PET','PA6','PA66']]==[578,126,136]
 return rows,envs,lookup
def templates():
 out={}
 for name,pdb,res in [('IsPETase','5XJH',[160,237,206]),('LCC','4EB0',[165,242,210])]:
  rec=load_pdb(P/'catalytic_motif_compare_v1'/(pdb+'.pdb'),True)
  tri=[v for v in rec if int(v[0][22:26]) in res]
  keys=[(int(l[22:26]),l[12:16].strip()) for l,_,_ in tri]
  og=keys.index((res[0],'OG'));cb=keys.index((res[0],'CB'))
  q=np.array([v[1] for v in tri]);e=np.array([v[2] for v in tri])
  out[name]=(q-q[og],e,og,cb,[v[0] for v in tri])
 rec=load_pdb(Path('/data/bht2/new_seed_structure_sequence_folddisco_20260825_v3/structure_collection_20260826_v1/structures/folddisco/PDB/PDB_3AXG_A.pdb'),True)
 tri=[v for v in rec if int(v[0][22:26]) in [267,306,308]]
 keys=[(int(l[22:26]),l[12:16].strip()) for l,_,_ in tri];og=keys.index((267,'OG1'));cb=keys.index((267,'CB'))
 q=np.array([v[1] for v in tri]);out['NylC']=(q-q[og],np.array([v[2] for v in tri]),og,cb,[v[0] for v in tri])
 return out
def evaluate(item):
 global rng
 mat,row=item['material'],item['row'];sid=row['site_id'];key=mat+'_'+sid
 folder=OUT/('smoke' if SMOKE else 'cases');folder.mkdir(exist_ok=True)
 jf=folder/(key+'.json')
 if jf.exists():return json.loads(jf.read_text())
 start=time.monotonic();xyz,el=ENVS[mat]
 ids=[LOOKUP[mat][(row['segment'],row[k])] for k in ['carbonyl_C_name','carbonyl_O_name','ester_O_name']] if mat=='PET' else [LOOKUP[mat][s] for s in row['serials']]
 a=xyz[ids].copy();center=a[0].copy()
 if mat=='PET':a[:,:2]-=np.round((a[:,:2]-center[:2])/BOX[:2])*BOX[:2]
 ex=norm(a[1]-center);ey=a[2]-center;ey=norm(ey-ex*np.dot(ex,ey));frame=np.column_stack((ex,ey,np.cross(ex,ey)))
 q=xyz-center
 if mat=='PET':q[:,:2]-=np.round(q[:,:2]/BOX[:2])*BOX[:2]
 keep=np.linalg.norm(q,axis=1)<=30
 env=trees(q[keep]@frame,el[keep]);aa=(a-center)@frame
 poses={};stats={}
 if mat=='PET':
  # Replay original world-coordinate RNG sequence before converting to common local frame.
  worldenv=trees(q[keep]+center,el[keep])
  oldfile=next((f for f in [P/'donor_geometry_batch_singleNH_v1/extension70_v1/cases'/(sid+'.json'),P/'donor_geometry_batch_singleNH_v1/production/cases'/(sid+'.json')] if f.exists()),None)
  assert oldfile is not None
  old=json.loads(oldfile.read_text())
  for name in ['IsPETase','LCC']:
   rng=np.random.default_rng(seed_for(sid+'|'+name+'|triad'))
   pp,ss=make_triads(a,frame,TEMPLATES[name],worldenv)
   reference=next(t['triad'] for t in old['templates'] if t['template']==name)
   assert ss==reference,(key,name,ss,reference)
   poses[name]=(pp-center)@frame;stats[name]=ss
 else:
  z=np.load(N/'tdd_donor_all_v1'/(key+'.npz'));assert np.allclose(z['amide_A'],a)
  poses['NylC']=(z['triads_A']-center)@frame;stats['NylC']=dict(poses_evaluated=len(z['triads_A']),source='saved original TDD poses including no-pair poses')
 arrays={'amide_A':aa,'origin_A':center,'frame':frame};results=[]
 for name,pp in poses.items():arrays[name+'_triads_A']=pp
 for ck,c in CLOUDS.items():
  chemical_ids=np.flatnonzero(chemical(c['donor_A'],c['H_A'],aa[1],c['donor_elements']))
  good=chemical_ids[clear(c['q_A'][chemical_ids],np.array([RAD[e] for e in c['elements']]),env)>=-1e-8] if len(chemical_ids) else chemical_ids
  arrays[ck+'_polymer_ids']=good
  for name,pp in poses.items():
   counts=[];densitybest=[];representatives=[]
   for j,pose in enumerate(pp):
    ok=good[clear(c['q_A'][good],np.array([RAD[e] for e in c['elements']]),trees(pose,TEMPLATES[name][1]))>=-1e-8] if len(good) else good
    arrays[f'{ck}_{name}_{j}_ids']=ok
    counts.append(len(ok));best=[0,0];chosen=[]
    if len(ok):
     # Density uses only accepted coupled poses, identical triad and donor atom pair.
     dd=np.zeros((len(ok),2),int)
     for pair in np.unique(c['pair_indices'][ok],axis=0):
      ix=np.flatnonzero((c['pair_indices'][ok]==pair).all(1))
      for d in range(2):
       pos=c['donor_A'][ok[ix],d]
       dd[ix,d]=cKDTree(pos).query_ball_point(pos,1.,return_length=True)
     order=np.lexsort((ok,-dd.sum(1),-dd.min(1)));b=int(order[0])
     best=dd[b].tolist();chosen=sorted(set([int(ok[0]),int(ok[b])]))
     # Explicit independent angular replay on selected full-residue representatives.
     assert chemical(c['donor_A'][chosen],c['H_A'][chosen],aa[1],c['donor_elements'][chosen]).all()
    densitybest.append(best);representatives.append(chosen)
   results.append(dict(microstate=ck,template=name,chemical_candidates=len(chemical_ids),polymer_clear_candidates=len(good),accepted_per_triad=counts,best_density_per_triad=densitybest,representatives_per_triad=representatives,site_has_coupled_donor=any(counts),status='EVALUATED_FINITE_LIBRARY' if len(pp) else 'NOT_EVALUATED_NO_SAMPLED_TRIAD'))
 rec=dict(material=mat,site_id=sid,support=float(row['support_fraction'] if mat=='PET' else row['inward_support']),side=row['side'] if mat=='PET' else row['nearest_side'],triad_stats=stats,results=results,seconds=time.monotonic()-start)
 tmp=folder/(key+'.tmp.npz');np.savez_compressed(tmp,**arrays);tmp.replace(folder/(key+'.npz'))
 jf.write_text(json.dumps(rec,separators=(',',':')))
 print(json.dumps(dict(material=mat,site=sid,seconds=round(rec['seconds'],2),hit_microstates=sorted(set(r['microstate'] for r in results if r['site_has_coupled_donor'])))),flush=True)
 return rec
if __name__=='__main__':
 SMOKE='--smoke' in sys.argv;OUT.mkdir(exist_ok=True)
 config=dict(source_sha256=sha(Path(__file__)),intrinsic_sha256=sha(R/'coupled_torsion_pilot_v1/summary.json'),sites=dict(PET=578,PA6=126,PA66=136),samples_per_microstate=8192,oxygen_reference_A=[1.23,0,0],density='existing accepted coupled placements; same triad and donor atom pair; center included; no dedup',scope='finite unrelaxed intrinsic library; single static frame per material; heavy-atom geometry, not activity; per-site union is not fixed-core coverage')
 cp=OUT/'config.json'
 if cp.exists():assert json.loads(cp.read_text())==config
 else:cp.write_text(json.dumps(config,indent=2))
 (OUT/'RUNBOOK.md').write_text('CPU only, existing unbiased_100chain_v1/env/bin/python. Set OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1. Run surface_scan_v1.py --smoke then --run. Eight processes. Resume identical script/config; completed site JSON files are checkpoints. Original files preserved. PET XY minimum-image, 30A environment with donor extent <26A; original triad replay counts required. Nylon uses saved original triads and finite trimmed environment. Common C-based frame holds candidate residue geometry fixed, true oxygen geometry retested. Full coupled residues including caps are sterically tested. Density counts only pre-existing accepted poses under same triad and atom pair. No coordinates are newly perturbed for density. No per-site PDBs. The first stage site-wise union is not fixed-core coverage. Stage two must transfer selected full triad+residue motifs unchanged before claiming a robust fixed-core winner.\n')
 CLOUDS=build_clouds();ROWS,ENVS,LOOKUP=get_inputs();TEMPLATES=templates()
 # Existing analytic clash checks.
 assert clear(np.array([[[3.,0,0]]]),np.array([1.7]),trees(np.zeros((1,3)),np.array(['C'])))[0]>=-1e-12
 assert clear(np.array([[[2.9,0,0]]]),np.array([1.7]),trees(np.zeros((1,3)),np.array(['C'])))[0]<0
 if SMOKE:
  chosen=[]
  for mat in ['PET','PA6','PA66']:
   rr=[r for r in ROWS if r['material']==mat]
   if mat=='PET':chosen.append(next(r for r in rr if r['row']['site_id']=='P244_C42_O18'))
   else:
    chosen.append(next(r for r in rr if len(np.load(N/'tdd_donor_all_v1'/(mat+'_'+r['row']['site_id']+'.npz'))['triads_A'])))
  ans=[evaluate(r) for r in chosen]
  (OUT/'SMOKE_PASS.json').write_text(json.dumps(dict(source_sha256=config['source_sha256'],sites=len(ans),site_ids=[r['material']+'_'+r['site_id'] for r in ans]),indent=2))
 else:
  assert json.loads((OUT/'SMOKE_PASS.json').read_text())['source_sha256']==config['source_sha256']
  with concurrent.futures.ProcessPoolExecutor(max_workers=8) as pool:ans=list(pool.map(evaluate,ROWS))
  (OUT/'stage1_summary.json').write_text(json.dumps(dict(completed_sites=len(ans),fixed_core_ranking='NOT_EVALUATED_STAGE2_REQUIRED',results=[dict(material=m,microstate=k,denominator=sum(a['material']==m for a in ans),sites_with_any_placement=sum(a['material']==m and any(r['microstate']==k and r['site_has_coupled_donor'] for r in a['results']) for a in ans)) for m in ['PET','PA6','PA66'] for k in CLOUDS]),indent=2))
 with (R/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__)),command='--smoke' if SMOKE else '--run',source_sha256=config['source_sha256'],output=str(OUT),completed_sites=len(ans),exit_code=0))+'\n')
 print('SMOKE_PASS' if SMOKE else 'STAGE1_COMPLETE',flush=True)
