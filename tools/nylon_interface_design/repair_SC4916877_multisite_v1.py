"""Versioned six-site CPU repair adapter; preserves prior minimization settings."""
from pathlib import Path
import hashlib
r=Path(__file__).resolve().parent
p=r/'repair_supported_arg_same_protocol_v1.py'
s=p.read_text()
assert hashlib.sha256(p.read_bytes()).hexdigest()=='ef7eed4578fdd8a6f1ff2066e5a37d372fc4f5fd099c328744aca8a689a1803e'
oldroot='route6_supported_ARG_SC4559268_v1'
s=s.replace(oldroot+'/fixed_core_repair_v1','route6_top50_local_generation_v1/SC4916877_seed2_multisite_repair_v1')
s=s.replace('n9_core_donor_design_20260910_gpu_outputs_v1/'+oldroot+'/PA6_SC4559268_coverage10_variant1/design_0.pdb','n9_core_donor_design_20260910_gpu_outputs_v1/route6_top50_local_shapes_3seeds_v2/PA6_SC4916877_coverage10_variant1/design_2.pdb')
s=s.replace(oldroot+'/PA6_SC4559268_coverage10_variant1/motif.pdb','route6_top50_local_generation_v1/candidate_specific_shapes_v2/SC4916877/motif.pdb')
needle='for i,j in mapping.items():'
assert s.count(needle)==1
s=s.replace(needle,'''# Align all output atoms once before restoring exact motif.
aa=[];bb=[]
for i,j in mapping.items():
 for name in ['N','CA','C','O']:
  v=pose.residue(i).xyz(name);aa.append([v.x,v.y,v.z]);bb.append(refatoms[j][name]-delta)
aa=np.array(aa);bb=np.array(bb);u,sv,vh=np.linalg.svd((aa-aa.mean(0)).T@(bb-bb.mean(0)));rot=u@np.diag([1,1,np.linalg.det(u@vh)])@vh;shift=bb.mean(0)-aa.mean(0)@rot
for ii in range(1,pose.size()+1):
 for an in range(1,pose.residue(ii).natoms()+1):
  aid=R.core.id.AtomID(an,ii);v=pose.xyz(aid);q=np.array([v.x,v.y,v.z])@rot+shift;pose.set_xyz(aid,R.numeric.xyzVector_double_t(*map(float,q)))
for i,j in mapping.items():''')
start=s.index("z=np.load(d/")
end=s.index("anchor=R.core.id.AtomID",start)
s=s[:start]+'''j=json.loads((d/'route6_top50_local_generation_v1/candidate_specific_shapes_v2/SC4916877/PA6_SC4916877_coverage10_variant1/manifest.json').read_text())
assert j['selected_site_indices']==[204,358,460,543,549,648]
assert hashlib.sha256(Path(j['pdb']).read_bytes()).hexdigest()==j['pdb_sha256']
sites=json.loads((d/'N10_fixed_minimal_motif_material_v1/sites.json').read_text())
data=np.load(p/'inputs_v2/PA6_parent_atoms.npz');els=data['elements'].astype(str);box=data['box_A'];trees={};material={}
for el in np.unique(els[els!='H']):
 v=data['xyz_A'][els==el].copy();v[:,:2]%=box[:2];material[el]=v;trees[el]=cKDTree(v,boxsize=[box[0],box[1],0])
def world(q,site):
 w=(q+delta)@np.array(site['frame']).T+site['center'];w[:,:2]%=box[:2];return w
''' +s[end:]
start=s.index(" w=world(q[None,:])[0]")
end=s.index("sf=py.create_score_function",start)
s=s[:start]+''' for si in j['selected_site_indices']:
  site=sites[si];F=np.array(site['frame']);w=world(q[None,:],site)[0]
  for el,t in trees.items():
   for idx in t.query_ball_point(w,8.):
    dx=material[el][idx]-w;dx[:2]-=np.round(dx[:2]/box[:2])*box[:2];target=q+dx@F;lower=float(radii[k]+rad[el]-.4);fun=R.core.scoring.constraints.BoundFunc(lower,100000.,.2,'material_exclusion')
    pose.add_constraint(R.core.scoring.constraints.CoordinateConstraint(aid,anchor,R.numeric.xyzVector_double_t(*map(float,target)),fun));nconstraints+=1
''' +s[end:]
start=s.index(" q=np.array([[pose.xyz(a).x",s.index('def audit():'))
end=s.index(" cn=[",start)
s=s[:start]+''' q=np.array([[pose.xyz(a).x,pose.xyz(a).y,pose.xyz(a).z] for a in heavy]);mar=np.full(len(q),np.inf);per_site=[]
 for si in j['selected_site_indices']:
  w=world(q,sites[si]);local=np.full(len(q),np.inf)
  for el,t in trees.items():local=np.minimum(local,t.query(w)[0]-radii-rad[el]+.4)
  mar=np.minimum(mar,local);per_site.append(dict(site_index=si,min_margin_A=float(local.min()),clashing_atoms=int((local< -1e-8).sum())))
''' +s[end:]
s=s.replace("return dict(CN_min_A=", "return dict(selected_site_audit=per_site,CN_min_A=")
s=s.replace("material_site_index=257,material_site=site['site_id']", "selected_site_indices=j['selected_site_indices'],conditional_denominator=j['conditional_denominator']")
s=s.replace("print('START',flush=True)", """(out/'RUNBOOK.md').write_text('Single CPU diagnostic, same ref2015_cart and 500-iteration settings as baseline. All nine motif heavy-atom groups fixed after one rigid alignment. Six original PA6 site constraints simultaneously applied with 8 A neighbor radius. New version only, no retries. Coordinate SD 0.5 A, exclusion SD 0.2 A. Material and geometry independent audit required; placeholder sequence is not designed enzyme.\\n')
print('START',flush=True)""")
compile(s,str(p),'exec')
exec(compile(s,str(p),'exec'),dict(__file__=str(Path(__file__).resolve()),__name__='__main__'))
