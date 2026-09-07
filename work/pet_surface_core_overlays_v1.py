"""Align existing dense-PET water-facing surface patches to six fixed catalytic cores.
No new enzyme pose, MD, density parameter change, or atomistic surface calculation.
"""
from pathlib import Path
import json,gzip,hashlib,datetime,shutil
import numpy as np
from scipy.spatial.distance import cdist
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
P=Path("/data/bht2/polymer_material_reference_20260827/simulation_slabs/pet_dp10_400chain_direct_v1")
out=P/"pet_surface_core_overlays_v1";out.mkdir(exist_ok=False)
B=P/"donor_geometry_batch_singleNH_v1/joint_transfer_map_v1"
E=P/"donor_geometry_batch_singleNH_v1/extension70_v1/eligible_sites.json"
S=P/"surface_chemistry_shape_20260907_v1/bulk_density_preview_v1/bulk_density_preview.json.gz"
D=P/"pet_dp10_400chain/water_slab_v1/dry_export/PET_DP10_400chain_water_equilibrated_dry.pdb"
CORE=P/"motif_top6_vs_native_lcc_v1"
paths=[B/"joint_transfer_results.npz",B/"motifs.json",E,S,D]+list(CORE.glob("*_triad_aligned.pdb"))
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}
assert hashes[str(D)]=="e8bbd2c435a4b84f0333c595add80b9e7b0918471f70867d30378bb71eab2318"
z=np.load(B/"joint_transfer_results.npz");meta=json.loads((B/"motifs.json").read_text());rows=json.loads(E.read_text())
assert np.array_equal(z["site_ids"],[r["site_id"] for r in rows])
raw=json.load(gzip.open(S,"rt"));box=np.array([125.697,125.697,182.848]);n=np.ceil(box[:2]).astype(int);step=box[:2]/n
assert abs(raw["threshold_Da_A3"]-.36123088845440243)<1e-10
# Re-use saved boundary heights, adding only omitted periodic seam cells.
meshes={}
for side in ["top","bottom"]:
 old=np.array(raw["sides"][side]);h=np.full(tuple(n),np.nan)
 pts=old.reshape(-1,3);ix=np.rint(pts[:,:2]/step-.5).astype(int)
 h[ix[:,0],ix[:,1]]=pts[:,2]
 assert np.isfinite(h).all()
 assert np.max(abs(h[ix[:,0],ix[:,1]]-pts[:,2]))<1e-9
 tris=[]
 for i in range(n[0]):
  for j in range(n[1]):
   corners=[]
   for a,b in [(i,j),(i+1,j),(i+1,j+1),(i,j+1)]:
    corners.append([(a+.5)*step[0],(b+.5)*step[1],h[a%n[0],b%n[1]]])
   q=np.array(corners)
   tris.extend([q[[0,1,2]],q[[0,2,3]]] if side=="top" else [q[[0,2,1]],q[[0,3,2]]])
 meshes[side]=np.array(tris)
 assert len(meshes[side])==31752
atomlines=[l for l in D.read_text().splitlines() if l.startswith(("ATOM  ","HETATM"))]
coords=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in atomlines]);coords[:,:2]%=box[:2]
lookup={(l[72:76].strip(),l[12:16].strip()):j for j,l in enumerate(atomlines)}
assert len(coords)==92400
def norm(x):return x/np.linalg.norm(x,axis=-1,keepdims=True)
def frame(r):
 ids=[lookup[r["segment"],r[k]] for k in ["carbonyl_C_name","carbonyl_O_name","ester_O_name"]]
 e=coords[ids].copy();c=e[0].copy();e[:,:2]-=np.round((e[:,:2]-c[:2])/box[:2])*box[:2]
 x=norm(e[1]-c);y=e[2]-c;y=norm(y-x*np.dot(x,y));f=np.stack([x,y,np.cross(x,y)],axis=1)
 assert abs(np.linalg.det(f)-1)<1e-10
 return e,f
def component(tri):
 _,inv=np.unique(np.round(tri.reshape(-1,3),5),axis=0,return_inverse=True);vf=inv.reshape(-1,3)
 edges=np.sort(np.concatenate([vf[:,[0,1]],vf[:,[1,2]],vf[:,[2,0]]]),axis=1);owners=np.tile(np.arange(len(tri)),3)
 order=np.lexsort((edges[:,1],edges[:,0]));same=np.all(edges[order][1:]==edges[order][:-1],axis=1)
 a=owners[order][:-1][same];b=owners[order][1:][same]
 graph=coo_matrix((np.ones(len(a)*2),(np.r_[a,b],np.r_[b,a])),shape=(len(tri),len(tri))).tocsr()
 nc,lab=connected_components(graph,directed=False)
 nearest=int(np.argmin(np.linalg.norm(tri.mean(1),axis=1)))
 keep=lab==lab[nearest]
 assert connected_components(graph[keep][:,keep],directed=False,return_labels=False)==1
 return keep,int(nc)
names=["M449","M346","M218","M318","M329","M450"]
indices=[next(i for i,m in enumerate(meta) if m["id"]==name) for name in names]
used=np.flatnonzero(z["coverage"][indices].any(axis=0))
assert len(used)==58 and int(z["coverage"][indices].sum())==153
patches={};audits=[]
for j in used:
 r=rows[j];e,F=frame(r);c=e[0]
 tri=meshes[r["side"]].copy()
 tri[:,:,:2]-=np.round((tri.mean(1)[:,:2]-c[:2])/box[:2])[:,None,:]*box[:2]
 local=tri-c
 keep=np.linalg.norm(local.mean(1),axis=1)<=15
 local=local[keep];which,ncomp=component(local);local=local[which]
 q=local@F
 v=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);normals=norm(v)
 near=int(np.argmin(np.linalg.norm(q.mean(1),axis=1)))
 le=(e-c)@F
 assert np.linalg.norm(q.mean(1),axis=1).max()<=15+1e-8
 assert np.linalg.norm(q@F.T+c-tri[keep][which],axis=-1).max()<1e-8
 outward=np.array([0,0,1 if r["side"]=="top" else -1])@F
 assert np.all(normals@outward>0)
 patches[r["site_id"]]={"side":r["side"],"support":r["support_fraction"],"triangles_ester_frame_A":np.round(q,5).tolist(),"ester_A":le.tolist(),"water_arrow_start_A":q[near].mean(0).tolist(),"water_arrow_unit":normals[near].tolist()}
 audits.append({"site":r["site_id"],"triangles":len(q),"components_before":ncomp,"components_kept":1,"radius_centroid_max_A":float(np.linalg.norm(q.mean(1),axis=1).max()),"nearest_bulk_surface_to_C_A":float(np.linalg.norm(q.mean(1),axis=1).min())})
models=[]
for i,name in zip(indices,names):
 path=CORE/(name+"_triad_aligned.pdb")
 aligned=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in path.read_text().splitlines() if l.startswith("ATOM") and l[76:78].strip()!="H"])
 q=z["motif_xyz_A"][i];ac=q.mean(0);bc=aligned.mean(0);u,_,vt=np.linalg.svd((q-ac).T@(aligned-bc));R=u@vt
 if np.linalg.det(R)<0:u[:,-1]*=-1;R=u@vt
 t=bc-ac@R
 err=float(np.linalg.norm(q@R+t-aligned,axis=1).max())
 assert err<.001 and abs(np.linalg.det(R)-1)<1e-10
 sites=[rows[j]["site_id"] for j in np.flatnonzero(z["coverage"][i])]
 sites.sort(key=lambda sid:(-patches[sid]["support"],sid))
 models.append({"id":name,"triad_source":meta[i]["template"],"sites":sites,"R":R.tolist(),"t":t.tolist(),"alignment_recovery_max_error_A":err,"cohort_counts":{str(s):int((z["coverage"][i]&(z["site_support"]>=s)).sum()) for s in [.9,.8,.7]}})
 shutil.copy2(path,out/(name+"_core.pdb"))
 # Independent target-ester geometry check against the old model coordinates.
 for sid in sites:
  j=next(k for k,r in enumerate(rows) if r["site_id"]==sid);e,F=frame(rows[j])
  le=(e-e[0])@F
  og=q[meta[i]["functional_indices"]["og"]]
  assert 3-1e-8<=np.linalg.norm(og)<=3.4+1e-8
  assert np.max(abs(cdist((q@R+t),le@R+t)-cdist(q,le)))<1e-8
payload={"models":models,"patches":patches,"radius_A":15,"surface_kind":"existing sigma=3 A half-density dense-PET top/bottom water-facing envelope, NOT atomistic SES","threshold_Da_A3":raw["threshold_Da_A3"]}
with gzip.open(out/"overlay_data.json.gz","wt") as f:json.dump(payload,f,separators=(",",":"))
summary={"status":"EXISTING_SURFACES_REGISTERED","candidate_count":6,"core_site_pairings":153,"unique_sites":58,"counts":{m["id"]:len(m["sites"]) for m in models},"surface_kind":payload["surface_kind"],"patch":"same water-facing side only, triangle centroids within 15 A of carbonyl C; retain shared-edge component nearest C; periodic XY seam cells reconstructed from unchanged saved heights","registration":"original site ester frame, then same proper rigid transform as the extracted triad-aligned core; no independent surface alignment","interpretation":"Descriptive overlay of matched surface environments; not an atomic exclusion surface, unique protein shape, hydrophobicity field, or validated RFD constraint.","source_sha256":hashes,"patch_audit":audits,"models":models}
(out/"summary.json").write_text(json.dumps(summary,indent=2))
(out/"RUNBOOK.md").write_text("# PET surfaces relative to fixed catalytic cores\n\nCPU extraction and coordinate transformations only. No new docking, density parameters, MD, or sampling. Inputs are the existing 473x578 coverage matrix and existing sigma 3 A half-density top/bottom boundary. The boundary removes atomistic chain texture by coarse graining; it is not proof of dynamic free-chain classification. Same-side, circular 15 A centroid patches; keep nearest shared-edge connected component. Add periodic seam cells from saved boundary heights rather than leaving artificial XY gaps.\n\nEvery site is registered by C origin, C->O X axis, ester-plane Y axis, right-handed Z, then moved with the candidate core into the previous triad-aligned frame. Do not independently align PET patches or average away occlusion. Six core-site counts are 23,20,24,26,29,31; shared sites are not independent replicates. Surface colors identify sites, not chemistry. Water arrows are local outward normals to the density envelope. Output is a shape-distribution preview, not a ready-made RFD shape or atomic collision guarantee.\n")
for p,h in hashes.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
(out/"RUN_LOG.jsonl").write_text(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=__file__,status=summary["status"],inputs_unchanged=True))+"\n")
(out/"SHA256.json").write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()},indent=2))
print(json.dumps({"output":str(out),"counts":summary["counts"],"unique_sites":58,"pairings":153,"triangle_count_unique":sum(a["triangles"] for a in audits),"file_bytes":(out/"overlay_data.json.gz").stat().st_size,"status":summary["status"]}),flush=True)
