"""CPU-only zero-level degeneracy diagnostic, one previously rejected shape."""
import sys,json,hashlib,datetime
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.ndimage import label,map_coordinates
r=Path(__file__).resolve().parent;d=r/'route56_N10_20260915_v1';b=d/'route6_top50_local_generation_v1'
sys.path.insert(0,str(r.parent/'interface_surface_robustness_20260908_v1/deps_chem_contact_v1'))
from skimage.measure import marching_cubes
cid='SC4559268';name='PA6_SC4559268_coverage10_variant1';src=b/'candidate_specific_shapes_v2'/cid
j=json.loads((src/name/'manifest.json').read_text());assert j['status']=='NONMANIFOLD_MESH_REJECTED'
o=b/'zero_level_mesh_repair_v2'/name;o.mkdir(parents=True,exist_ok=False)
sites=json.loads((d/'N10_fixed_minimal_motif_material_v1/sites.json').read_text())
ca=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in (src/'motif.pdb').read_text().splitlines() if l.startswith('ATOM  ') and l[12:16].strip()=='CA'])
center=ca.mean(0);ax=np.arange(-32,33,dtype=float);grid=np.stack(np.meshgrid(ax,ax,ax,indexing='ij'),-1).reshape(-1,3)
z=np.load(r.parent/'interface_surface_robustness_20260908_v1/inputs_v2/PA6_parent_atoms.npz');px=z['xyz_A'].copy();box=z['box_A'];px[:,:2]%=box[:2]
rad={'C':1.7,'N':1.55,'O':1.52,'S':1.8};trees={e:cKDTree(px[z['elements']==e],boxsize=[*box[:2],0]) for e in rad if np.any(z['elements']==e)}
def clear(w):
 w=w.copy();w[:,:2]%=box[:2];v=np.full(len(w),np.inf)
 for e,t in trees.items():v=np.minimum(v,t.query(w)[0]-rad[e])
 return v
field=25-np.linalg.norm(grid,axis=1)
for i in j['selected_site_indices']:
 s=sites[i];field=np.minimum(field,(clear((grid+center)@np.array(s['frame']).T+s['center'])-1.3).astype(np.float32))
field=field.reshape(65,65,65);labs,n=label(field>0);points=(ca-center+32).T;cl=labs[tuple(np.rint(points).astype(int))]
assert np.all(cl>0) and len(set(cl))==1
field[(labs!=cl[0])&(field>0)]=-1
meshes={};stats=[]
for lev in [0.,1e-6]:
 v,f,_,_=marching_cubes(field.astype(np.float32),lev,allow_degenerate=False);v=v.astype(np.float64);v-=32
 edges=np.sort(np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]),axis=1);_,cnt=np.unique(edges,axis=0,return_counts=True)
 stats.append(dict(level=lev,bad_edges=int((cnt!=2).sum()),vertices=len(v),faces=len(f)));meshes[lev]=(v,f)
assert stats[0]['bad_edges']>0 and stats[1]['bad_edges']==0
v0,_=meshes[0.];v,f=meshes[1e-6]
haus=max(float(cKDTree(v0).query(v)[0].max()),float(cKDTree(v).query(v0)[0].max()))
fixed_values=map_coordinates(field,points,order=1);assert np.all(fixed_values>1e-6)
frag=np.load(b/'local_fragments_three_templates_v1/SC4559268_N10_170_174.npz')
core=np.load(d/'route6_N10_single_residue_pilot_v1/ARG_1.npz');q=frag['core_original'];rr=core['r']
checks=[]
for i in j['selected_site_indices']:
 s=sites[i];F=np.array(s['frame']);origin=np.array(s['center'])
 cm=float((clear(q@F.T+origin)-rr+.4).min())
 fm=float((clear(frag['fragment']@F.T+origin)-frag['radii']+.4).min())
 assert cm>=-1e-8 and fm>=-1e-8
 checks.append(dict(site_index=i,core_margin_A=cm,fragment_margin_A=fm,CA_mesh_min_margin_A=float((clear((v+center)@F.T+origin)-1.3).min())))
obj=o/'shape.obj'
with obj.open('x') as h:
 for a in v:h.write('v %.9f %.9f %.9f\n'%tuple(a))
 for a in f+1:h.write('f %d %d %d\n'%tuple(a))
serialized=np.array([[float(t) for t in line.split()[1:]] for line in obj.read_text().splitlines() if line.startswith('v ')])
assert len(np.unique(serialized,axis=0))==len(serialized)
np.savez_compressed(o/'shape_field.npz',field=field,axis_A=ax,center_A=center,level=1e-6)
old=json.loads((d/'route6_supported_ARG_SC4559268_v1/manifest.json').read_text())['jobs'][0]
out=r.parent/'n9_core_donor_design_20260910_gpu_outputs_v1/route6_zero_level_mesh_repair_v2'/name
cmd=[a.replace(old['pdb'],str(src/'motif.pdb')).replace(old['shape'],str(obj)).replace(old['output'],str(out)).replace('inference.num_designs=1','inference.num_designs=3') for a in old['command']]
assert cmd[0]=='/usr/local/bin/gpurun'
report=dict(status='NUMERICAL_MESH_REPAIR_PREPARED_NOT_SUBMITTED',source=j,mesh_stats=stats,nearest_vertex_symmetric_max_A=haus,positive_grid_class_changes=int(((field>0)!=(field>1e-6)).sum()),fixed_CA_field_A=fixed_values.tolist(),exact_material_checks=checks,pdb=str(src/'motif.pdb'),shape=str(obj),pdb_sha256=hashlib.sha256((src/'motif.pdb').read_bytes()).hexdigest(),shape_sha256=hashlib.sha256(obj.read_bytes()).hexdigest(),command=cmd,output=str(out),limitations=['Nearest vertex distance is not continuous surface Hausdorff distance','CA envelope interpolation error retained and reported; not full atom collision proof','No core, material subset, denominator or collision gate altered','Finite sampling and sequence/foldability not evaluated'])
with (o/'manifest.json').open('x') as h:json.dump(report,h,indent=2)
(o/'RUNBOOK.md').write_text('Recompute the identical selected-site field. Compare marching cubes level zero and positive 1e-6. Strict closure and unique serialized vertices required. Exact unchanged core/fragment material collision checks required. Keep original rejected mesh evidence. No GPU submission by this script. Potential prepared generation command uses three seeds; needs scheduler coordination before launch.\n')
with (r/'RUN_LOG.jsonl').open('a') as h:h.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),event='zero_level_mesh_repair',output=str(o),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),status=report['status']))+'\n')
print(json.dumps({k:v for k,v in report.items() if k not in ['source','command','exact_material_checks']}))
