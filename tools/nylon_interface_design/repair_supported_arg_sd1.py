import os
os.environ.update(OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1')
from pathlib import Path
import numpy as np,json,datetime,hashlib
from scipy.spatial import cKDTree
import pyrosetta as py
r=Path(__file__).resolve().parent;d=r/'route56_N10_20260915_v1';p=r.parent/'interface_surface_robustness_20260908_v1';out=Path('/data/bht2/polymer_material_reference_20260827/simulation_slabs/n9_core_donor_design_20260910_v1/route56_N10_20260915_v1/route6_supported_ARG_SC4559268_v1/fixed_core_repair_sd1_v1');out.mkdir(exist_ok=False)
py.init('-mute all -constant_seed -jran 101 -ignore_unrecognized_res false');R=py.rosetta
pose=py.pose_from_pdb('/data/bht2/polymer_material_reference_20260827/simulation_slabs/n9_core_donor_design_20260910_gpu_outputs_v1/route6_supported_ARG_SC4559268_v1/PA6_SC4559268_coverage10_variant1/design_0.pdb');assert pose.size()==126 and pose.num_chains()==1 and pose.residue(1).is_lower_terminus()
refatoms={}
for l in Path('/data/bht2/polymer_material_reference_20260827/simulation_slabs/n9_core_donor_design_20260910_v1/route56_N10_20260915_v1/route6_supported_ARG_SC4559268_v1/PA6_SC4559268_coverage10_variant1/motif.pdb').read_text().splitlines():
 if l.startswith('ATOM'):refatoms.setdefault(int(l[22:26]),{})[l[12:16].strip()]=np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])
delta=np.mean([v['CA'] for v in refatoms.values()],axis=0);mapping=dict(zip([1,40,41,42,61,62,63,64,65],range(1,10)));fixed=[]
for i,j in mapping.items():
 for name,q in refatoms[j].items():
  aid=R.core.id.AtomID(pose.residue(i).atom_index(name),i);pose.set_xyz(aid,R.numeric.xyzVector_double_t(*map(float,q-delta)));fixed.append(aid)
fixedkeys={(a.rsd(),a.atomno()) for a in fixed};ref=pose.clone();mm=R.core.kinematics.MoveMap();mm.set_bb(False);mm.set_chi(False);mm.set_jump(False);heavy=[];start=[];radii=[];rad={'C':1.7,'N':1.55,'O':1.52,'S':1.8}
for i in range(1,pose.size()+1):
 for a in range(1,pose.residue(i).natoms()+1):
  aid=R.core.id.AtomID(a,i);mm.set_atom(aid,(i,a) not in fixedkeys)
  if a<=pose.residue(i).nheavyatoms():
   name=pose.residue(i).atom_name(a).strip();heavy.append(aid);v=pose.xyz(aid);start.append([v.x,v.y,v.z]);radii.append(rad[name[0]])
start=np.array(start);radii=np.array(radii);cm=R.core.optimization.CartesianMinimizerMap();cm.setup(pose,mm);assert all(not cm.atom_is_moving(a) for a in fixed)
z=np.load(d/'route6_supported_ARG_SC4559268_v1/PA6_SC4559268_coverage10_variant1/motif_exact.npz');center=-z['translation_to_shape'];site=json.loads((d/'N10_fixed_minimal_motif_material_v1/sites.json').read_text())[257];F=np.array(site['frame'])@z['rotation_to_shape'];assert np.allclose(delta,z['new_center_A'],atol=1e-8);data=np.load(p/'inputs_v2/PA6_parent_atoms.npz');els=data['elements'].astype(str);box=data['box_A'];trees={};material={}
for el in np.unique(els[els!='H']):
 v=data['xyz_A'][els==el].copy();v[:,:2]%=box[:2];material[el]=v;trees[el]=cKDTree(v,boxsize=[box[0],box[1],0])
def world(q):
 w=(q+center+delta)@F.T+site['center'];w[:,:2]%=box[:2];return w
anchor=R.core.id.AtomID(pose.residue(1).atom_index('CA'),1);nconstraints=0
for k,(aid,q) in enumerate(zip(heavy,start)):
 if (aid.rsd(),aid.atomno()) in fixedkeys:continue
 pose.add_constraint(R.core.scoring.constraints.CoordinateConstraint(aid,anchor,R.numeric.xyzVector_double_t(*map(float,q)),R.core.scoring.func.HarmonicFunc(0,1.0)))
 w=world(q[None,:])[0]
 for el,t in trees.items():
  for idx in t.query_ball_point(w,8.):
   dx=material[el][idx]-w;dx[:2]-=np.round(dx[:2]/box[:2])*box[:2];target=q+dx@F;lower=float(radii[k]+rad[el]-.4);fun=R.core.scoring.constraints.BoundFunc(lower,100000.,.2,'material_exclusion')
   pose.add_constraint(R.core.scoring.constraints.CoordinateConstraint(aid,anchor,R.numeric.xyzVector_double_t(*map(float,target)),fun));nconstraints+=1
sf=py.create_score_function('ref2015_cart');sf.set_weight(R.core.scoring.cart_bonded,.5);sf.set_weight(R.core.scoring.pro_close,0);sf.set_weight(R.core.scoring.coordinate_constraint,1)
def vec(i,n):
 v=pose.residue(i).xyz(n);return np.array([v.x,v.y,v.z])
def angle(a,b,c):
 v=a-b;w=c-b;return float(np.degrees(np.arccos(np.clip(np.dot(v,w)/np.linalg.norm(v)/np.linalg.norm(w),-1,1))))
def audit():
 q=np.array([[pose.xyz(a).x,pose.xyz(a).y,pose.xyz(a).z] for a in heavy]);w=world(q);mar=np.full(len(q),np.inf)
 for el,t in trees.items():mar=np.minimum(mar,t.query(w)[0]-radii-rad[el]+.4)
 cn=[float(np.linalg.norm(vec(i,'C')-vec(i+1,'N'))) for i in range(1,126)];angles=[dict(pair=[i,i+1],CA_C_N_deg=angle(vec(i,'CA'),vec(i,'C'),vec(i+1,'N')),C_N_CA_deg=angle(vec(i,'C'),vec(i+1,'N'),vec(i+1,'CA'))) for i in range(1,126)]
 return dict(CN_min_A=min(cn),CN_max_A=max(cn),bad_CN=sum(v<1.2 or v>1.5 for v in cn),peptide_angles=angles,fixed_max_shift_A=max(pose.xyz(a).distance(ref.xyz(a)) for a in fixed),material_clashing_heavy_atoms=int((mar< -1e-8).sum()),minimum_material_margin_A=float(mar.min()),max_heavy_displacement_A=float(np.linalg.norm(q-start,axis=1).max()),energy_REU=float(sf(pose)),thr_N_terminal=pose.residue(1).is_lower_terminus())
report=dict(status='RUNNING_DIAGNOSTIC',before=audit(),material_site_index=257,material_site=site['site_id'],material_constraints=nconstraints,max_iterations=500,restraint_sd_A=1.0,exclusion_sd_A=.2,material_neighbor_radius_A=8,fixed_heavy_atoms=len(fixed),limitations=['Placeholder sequence, not complete designed protein','Fixed motif heavy atoms exact; hydrogens allowed to optimize','Soft material constraints require independent final replay','No folding or catalytic PASS'])
(out/'before.json').write_text(json.dumps(report,indent=2));pose.dump_pdb(str(out/'before.pdb'));print('START',flush=True)
mini=R.protocols.minimization_packing.MinMover();mini.movemap(mm);mini.score_function(sf);mini.min_type('lbfgs_armijo_nonmonotone');mini.tolerance(1e-6);mini.max_iter(500);mini.cartesian(True);mini.apply(pose)
report['after']=audit();report['status']='TECHNICAL_COMPLETE_REQUIRES_GEOMETRY_AUDIT';pose.dump_pdb(str(out/'after.pdb'));(out/'audit.json').write_text(json.dumps(report,indent=2))
with (r/'RUN_LOG.jsonl').open('a') as h:h.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__)),sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),exit_status=0,output=str(out)))+'\n')
print(json.dumps({label:{k:v for k,v in report[label].items() if k!='peptide_angles'} for label in ['before','after']}),flush=True)
