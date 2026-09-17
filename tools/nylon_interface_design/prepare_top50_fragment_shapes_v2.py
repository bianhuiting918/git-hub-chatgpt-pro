"""Prepare candidate-specific whole-site CAD shapes from audited local fragments.
CPU only. Existing whole-site subset rule, 25 A envelope, 1 A grid and
0.4 A overlap allowance retained. No RFD submission in this script.
"""
import os
for key in ['OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS']:os.environ[key]='1'
import sys,json,hashlib,datetime,math
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.ndimage import label,map_coordinates
sys.dont_write_bytecode=True
R=Path(__file__).resolve().parent;D=R/'route56_N10_20260915_v1'
P=R.parent/'interface_surface_robustness_20260908_v1'
sys.path.insert(0,str(P/'deps_chem_contact_v1'))
from skimage.measure import marching_cubes
import fit_route6_single_N10_whole_site_v1 as h
base=D/'route6_top50_local_generation_v1';src=base/'local_fragments_three_templates_v1'
O=base/'candidate_specific_shapes_v2';O.mkdir(exist_ok=False)
def dump(p,x):
 with p.open('x') as f:json.dump(x,f,indent=2)
def log(event,**kw):
 x=dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),event=event,**kw)
 print(json.dumps(x),flush=True)
 for p in [O/'events.jsonl',R/'RUN_LOG.jsonl']:
  with p.open('a') as f:f.write(json.dumps(x)+'\n')
rows=[json.loads(x) for x in (src/'results.jsonl').read_text().splitlines()]
good=[x for x in rows if x['status']=='LOCAL_FRAGMENT_GEOMETRY_PASS_NOT_RFD']
assert len(good)==19 and len({x['id'] for x in good})==19
sites=json.loads((D/'N10_fixed_minimal_motif_material_v1/sites.json').read_text())
panel=json.loads((D/'route6_N10_single_residue_pilot_v1/panel.json').read_text())
coverage=np.load(D/'route6_N10_single_residue_pilot_v1/coverage.npz')['coverage']
old=json.loads((D/'route6_supported_ARG_SC4559268_v1/manifest.json').read_text())['jobs'][0]
ax=np.arange(-32,33,dtype=float);grid=np.stack(np.meshgrid(ax,ax,ax,indexing='ij'),-1).reshape(-1,3)
rad={'C':1.7,'N':1.55,'O':1.52,'S':1.8};bb=['N','CA','C','O']
inside=np.linalg.norm(h.m.GRID,axis=1)<=15;weight=np.exp(-np.linalg.norm(h.m.GRID[inside],axis=1)/10)
dump(O/'manifest.json',dict(status='CPU_PREPARATION_RUNNING_NOT_RFD',candidates=19,materials=['PA6','PA66'],fractions=[.1,.2,.3],variants_max=3,exposure_scope='Pooled first batch; bin labels retained, bin-specific shapes not yet prepared',denominator='Original core eligible sites, NOT reduced after adding local backbone',source_sha256=hashlib.sha256((src/'results.jsonl').read_bytes()).hexdigest()))
(O/'RUNBOOK.md').write_text('Run this script using existing CPU Python from the N9/N10 project. Exclusive versioned outputs. No GPU submission. Local fragment and TDD coordinates unchanged; peripheral sidechains absent. Shape is a soft CA guide, not a full atom collision proof. Original core-eligible denominator retained; selected sites must also accommodate the added fragment. 25 A envelope, 1 A grid, 1.7 A CA radius minus 0.4 A allowance. All-seed greedy subsets and max-min Jaccard diversity reuse prior protocol. Pooled exposure cohort first; bin-specific jobs remain pending. Each prepared command uses gpurun and unchanged ActiveSite parameters.\n')
def diverse(sets):
 result=[sets[0]] if sets else [];remaining=sets[1:]
 while remaining and len(result)<3:
  k=max(range(len(remaining)),key=lambda i:(min(1-len(set(remaining[i])&set(y))/len(set(remaining[i])|set(y)) for y in result),-i))
  result.append(remaining.pop(k))
 return result
envs={}
for mat in ['PA6','PA66']:
 z=np.load(P/'inputs_v2'/(mat+'_parent_atoms.npz'));px=z['xyz_A'].copy();box=z['box_A'];px[:,:2]%=box[:2]
 envs[mat]=(box,{e:cKDTree(px[z['elements']==e],boxsize=[*box[:2],0]) for e in rad if np.any(z['elements']==e)})
def clear(points,mat):
 box,trees=envs[mat];w=points.copy();w[:,:2]%=box[:2];v=np.full(len(w),np.inf)
 for e,t in trees.items():v=np.minimum(v,t.query(w,workers=1)[0]-rad[e])
 return v
jobs=[]
for row in good:
 cid=row['id'];assert row['microstate']=='ARG_1'
 z=np.load(row['output']);q=z['core_original'];b=z['backbone'];fragment=z['fragment']
 assert hashlib.sha256(Path(row['output']).read_bytes()).hexdigest()==row['sha256']
 blocks=[('THR',['N','CA','C','O','CB','OG1','CG2'],q[:7]),('ASP',bb+['CB','CG','OD1','OD2'],q[7:15]),('GLY',bb,q[23:27]),('ASP',bb+['CB','CG','OD1','OD2'],q[15:23])]
 for j in range(5):blocks.append(('ARG',bb+['CB','CG','CD','NE','CZ','NH1','NH2'],np.vstack([b[8:12],q[32:39]])) if j==2 else ('GLY',bb,b[4*j:4*j+4]))
 lines=[];serial=0;ca=[]
 for ri,(rn,names,coords) in enumerate(blocks,1):
  for name,xyz in zip(names,coords):
   serial+=1;line=f'ATOM  {serial:5d} {name:^4s} {rn:3s} A{ri:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1:6.2f}{0:6.2f}          {name[0]:>2s}';lines.append(line)
   if name=='CA':ca.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
 ca=np.array(ca);center=ca.mean(0);dst=O/cid;dst.mkdir();pdb=dst/'motif.pdb';pdb.write_text('\n'.join(lines)+'\nEND\n')
 pi=next(i for i,p in enumerate(panel) if p['id']==cid)
 for mat in ['PA6','PA66']:
  ids=[int(i) for i in z['clear_site_indices'] if sites[int(i)]['material']==mat]
  original=[i for i,s in enumerate(sites) if s['material']==mat and coverage[pi,i]]
  assert set(ids)<=set(original)
  fields=[];local=[]
  for si in ids:
   s=sites[si];F=np.array(s['frame']);origin=np.array(s['center'])
   margin=clear(fragment@F.T+origin,mat)-z['radii']+.4
   assert margin.min()>=-1e-8
   fields.append((clear((grid+center)@F.T+origin,mat)-1.3).astype(np.float32))
   local.append(-clear(h.m.GRID[inside]@F.T+origin,mat))
  fields=np.array(fields);local=np.array(local)
  log('candidate_material_fields_ready',core=cid,material=mat,original_denominator=len(original),fragment_clear=len(ids))
  for fraction in [.1,.2,.3]:
   k=h.target_count(len(original),fraction)
   if k>len(ids):log('insufficient_fragment_clear',core=cid,material=mat,fraction=fraction);continue
   sets=diverse(h.subsets(local,k,weight))
   for vi,ix in enumerate(sets,1):
    name=f'{mat}_{cid}_coverage{round(fraction*100)}_variant{vi}';td=dst/name;td.mkdir()
    field=np.minimum(25-np.linalg.norm(grid,axis=1),fields[ix].min(0)).reshape(65,65,65)
    points=(ca-center+32).T;values=map_coordinates(field,points,order=1);labs,n=label(field>0);cl=labs[tuple(np.rint(points).astype(int))]
    rec=dict(name=name,core=cid,material=mat,fraction=fraction,variant=vi,conditional_denominator=len(original),fragment_clear_n=len(ids),selected_site_indices=[ids[i] for i in ix],selected_bins=[sites[ids[i]]['bin'] for i in ix],fixed_CA_field_A=values.tolist(),status='CORE_COMPONENT_OR_FIELD_FAILED')
    if np.all(values>0) and np.all(cl>0) and len(set(cl))==1:
     field[(labs!=cl[0])&(field>0)]=-1
     v,f,_,_=marching_cubes(field.astype(np.float32),0,spacing=(1,1,1),allow_degenerate=False);v-=32
     edges=np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1);_,counts=np.unique(edges,axis=0,return_counts=True);if not np.all(counts==2):
      rec.update(status='NONMANIFOLD_MESH_REJECTED',bad_edge_count=int(np.sum(counts!=2)));dump(td/'manifest.json',rec);jobs.append(rec);log('shape_case_rejected',name=name,status=rec['status']);continue
     obj=td/'shape.obj'
     with obj.open('x') as out:
      for a in v:out.write('v %.6f %.6f %.6f\n'%tuple(a))
      for a in f+1:out.write('f %d %d %d\n'%tuple(a))
     np.savez_compressed(td/'shape_field.npz',field=field,axis_A=ax,center_A=center,selected_site_indices=rec['selected_site_indices'])
     output=R.parent/'n9_core_donor_design_20260910_gpu_outputs_v1/route6_top50_local_shapes_v1'/name
     cmd=[a.replace(old['pdb'],str(pdb)).replace(old['shape'],str(obj)).replace(old['output'],str(output)) for a in old['command']]
     assert cmd[0]=='/usr/local/bin/gpurun'
     rec.update(status='PREPARED_NOT_SUBMITTED',command=cmd,pdb=str(pdb),shape=str(obj),output=str(output),pdb_sha256=hashlib.sha256(pdb.read_bytes()).hexdigest(),shape_sha256=hashlib.sha256(obj.read_bytes()).hexdigest(),fixed_output_residues=old['fixed_output_residues'],catalytic_output_residues=old['catalytic_output_residues'],support_identity_note=old['support_identity_note'])
    dump(td/'manifest.json',rec);jobs.append(rec)
    log('shape_case_prepared',name=name,status=rec['status'])
  del fields,local
dump(O/'summary.json',dict(status='CPU_PREPARATION_COMPLETE_NOT_RFD',jobs=jobs,prepared=sum(x['status']=='PREPARED_NOT_SUBMITTED' for x in jobs),total=len(jobs)))
log('candidate_shape_preparation_complete',cases=len(jobs),prepared=sum(x['status']=='PREPARED_NOT_SUBMITTED' for x in jobs))
