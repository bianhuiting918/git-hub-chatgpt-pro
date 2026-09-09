#!/usr/bin/env python3
"""Tangential spherical probe: exact segment collisions, finite straight-ray accessibility."""
import sys,json,time,datetime,argparse,math
from pathlib import Path
import numpy as np
from scipy.spatial import ConvexHull
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
R=Path(__file__).resolve().parent;sys.path.insert(0,str(R))
import chemical_contact_surface_v3 as chem
old=chem.old
RADII=np.array([.5,1,1.5,2,3,4,5,6,8,10.])
def score(r,a):
 return float(np.trapz(np.asarray(a,float),np.log(r))/(np.log(r[-1])-np.log(r[0])))
def segment_hits(env,mid,u,length,r):
 hits=np.zeros(len(mid),bool)
 for el,tree in env.trees.items():
  search=r+old.RAD[el]+np.asarray(length)/2
  assert np.max(search)<min(env.box[:2])/2
  w=env.wrap(mid);neighbors=tree.query_ball_point(w,search,workers=2)
  sizes=np.array([len(q) for q in neighbors])
  if not sizes.sum():continue
  ids=np.repeat(np.arange(len(mid)),sizes);ix=np.concatenate([np.asarray(q,int) for q in neighbors if len(q)])
  d=tree.data[ix]-w[ids];d[:,:2]-=env.box[:2]*np.round(d[:,:2]/env.box[:2])
  t=np.clip(np.sum(d*u[ids],axis=1),-length[ids]/2,length[ids]/2)
  dsq=np.sum((d-t[:,None]*u[ids])**2,axis=1)
  hits[np.unique(ids[dsq<(r+old.RAD[el]-1e-7)**2])]=True
 return hits
def trace(env,center,u,r,max_path=120.,step=2.):
 u=np.asarray(u,float);assert np.allclose(np.linalg.norm(u,axis=1),1,atol=1e-7)
 origin=np.asarray(center)+(1.7+r)*u
 endpoint=env.clearance(origin)>=r-1e-7
 reached=np.zeros(len(u),bool);blocked=~endpoint;unresolved=np.zeros(len(u),bool)
 radmax=max(old.RAD[e] for e in env.trees)
 high=env.xyz[:,2].max()+radmax+r+.01;low=env.xyz[:,2].min()-radmax-r-.01
 uz=u[:,2];distance_to_exit=np.full(len(u),np.inf)
 pos=uz>1e-12;neg=uz<-1e-12
 distance_to_exit[pos]=np.maximum((high-origin[pos,2])/uz[pos],0)
 distance_to_exit[neg]=np.maximum((low-origin[neg,2])/uz[neg],0)
 traveled=np.zeros(len(u));active=endpoint.copy()
 for _ in range(int(math.ceil(max_path/step))+2):
  done=active&(traveled>=distance_to_exit-1e-9);reached[done]=True;active[done]=False
  cap=active&(traveled>=max_path-1e-9);unresolved[cap]=True;active[cap]=False
  ids=np.flatnonzero(active)
  if not len(ids):break
  length=np.minimum(step,np.minimum(distance_to_exit[ids]-traveled[ids],max_path-traveled[ids]))
  mid=origin[ids]+(traveled[ids]+length/2)[:,None]*u[ids]
  collision=segment_hits(env,mid,u[ids],length,r)
  blocked[ids[collision]]=True;active[ids[collision]]=False
  traveled[ids]+=length
 unresolved[active]=True
 assert np.all(reached<=endpoint) and np.all(reached.astype(int)+blocked+unresolved==1)
 return dict(endpoint=endpoint,reached=reached,blocked=blocked,unresolved=unresolved,traveled_A=traveled)
def graph(directions):
 f=ConvexHull(directions).simplices;i=np.r_[f[:,0],f[:,1],f[:,2]];j=np.r_[f[:,1],f[:,2],f[:,0]]
 return coo_matrix((np.ones(len(i)),(i,j)),shape=(len(directions),len(directions))).tocsr()
def opening(mask,u,g):
 ix=np.flatnonzero(mask)
 if not len(ix):return dict(component_count=0,largest_fraction=0.,largest_solid_angle_sr=0.,mean_direction=None)
 n,lab=connected_components(g[ix][:,ix],directed=False);counts=np.bincount(lab);best=ix[lab==counts.argmax()]
 mean=u[best].mean(0);norm=np.linalg.norm(mean)
 return dict(component_count=int(n),largest_fraction=len(best)/len(u),largest_solid_angle_sr=4*np.pi*len(best)/len(u),mean_direction=(mean/norm).tolist() if norm>1e-9 else None)
def catalogue(mat,ce):
 mol=chem.Chem.SDMolSupplier(str(chem.sdf_for(mat)),removeHs=False)[0];n=mol.GetNumAtoms();motifs=[]
 for atom in mol.GetAtoms():
  if atom.GetSymbol()!="C":continue
  oo=[b.GetOtherAtom(atom).GetIdx() for b in atom.GetBonds() if b.GetBondTypeAsDouble()==2 and b.GetOtherAtom(atom).GetSymbol()=="O"]
  leaves=[b.GetOtherAtom(atom) for b in atom.GetBonds() if b.GetBondTypeAsDouble()==1 and b.GetOtherAtom(atom).GetSymbol()==("O" if mat=="PET" else "N")]
  for l in leaves:
   if oo and any(v.GetSymbol()=="C" and v.GetIdx()!=atom.GetIdx() for v in l.GetNeighbors()):motifs.append([atom.GetIdx(),oo[0],l.GetIdx()])
 rows=[];xyz=ce.env.xyz
 for chain in range(400):
  for motif in motifs:
   ids=np.array(motif)+chain*n;a=xyz[ids].copy();a[:,:2]-=ce.env.box[:2]*np.round((a[:,:2]-a[0,:2])/ce.env.box[:2])
   ex=a[1]-a[0];ex/=np.linalg.norm(ex);ey=a[2]-a[0];ey-=ex*np.dot(ex,ey);ey/=np.linalg.norm(ey)
   f=np.column_stack([ex,ey,np.cross(ex,ey)])
   rows.append(dict(site_id=f"P{chain+1:03d}_{ce.data['atom_names'][ids[0]]}_{ce.data['atom_names'][ids[1]]}",chain=chain+1,parent_indices=ids.tolist(),center=xyz[ids[0]].tolist(),frame=f.tolist()))
 return rows
def sample_sites(ce,rows,k):
 # All parent bond candidates participate; no density/support selection.
 u=old.fibonacci(64);fractions=[]
 for start in range(0,len(rows),128):
  ss=rows[start:start+128];c=np.array([s["center"] for s in ss]);f=np.array([s["frame"] for s in ss])
  p=c[:,None,:]+2.2*np.einsum("nj,bkj->bnk",u,f)
  fractions.extend((ce.env.clearance(p.reshape(-1,3)).reshape(len(ss),64)>=.5-1e-7).mean(1).tolist())
 a=np.array(fractions);groups=[np.where(a==0)[0],np.where((a>0)&(a<=.05))[0],np.where((a>.05)&(a<=.15))[0],np.where(a>.15)[0]]
 chosen=[]
 for ids in groups:
  if len(ids):chosen.extend(ids[np.linspace(0,len(ids)-1,min(k//4,len(ids)),dtype=int)].tolist())
 if len(chosen)<k:
  left=np.array([i for i in range(len(rows)) if i not in chosen]);chosen.extend(left[np.linspace(0,len(left)-1,min(k-len(chosen),len(left)),dtype=int)].tolist())
 return chosen,a
def material(mat,out,n=1024,pilot=12):
 ce=chem.ChemEnv(mat);rows=catalogue(mat,ce);d=out/mat;d.mkdir()
 old.emit("catalogue_ready",material=mat,bonds=len(rows))
 selected,pre=sample_sites(ce,rows,pilot)
 np.savez_compressed(d/"all_parent_endpoint_prescreen.npz",endpoint_fraction=pre,selected_indices=selected)
 with (d/"all_parent_bonds.json").open("x") as f:json.dump(rows,f,separators=(",",":"))
 with (d/"selected_sites.json").open("x") as f:json.dump([dict(rows[i],prescreen_endpoint_fraction=pre[i]) for i in selected],f,indent=2)
 ul=old.fibonacci(n);g=graph(ul);results=[]
 legacy={s["site_id"] for s in json.load(open(R/"panel_surface_scan_v1/sites.json")) if s["surface"]==mat+"_original"}
 for ii,i in enumerate(selected):
  s=rows[i];u=ul@np.array(s["frame"]).T;rr=[];masks=[];endmasks=[];unresolved=[]
  for r in RADII:
   a=trace(ce.env,np.array(s["center"]),u,r)
   masks.append(a["reached"]);endmasks.append(a["endpoint"]);unresolved.append(a["unresolved"])
   rr.append(dict(radius_A=float(r),endpoint_n=int(a["endpoint"].sum()),accessible_n=int(a["reached"].sum()),unresolved_n=int(a["unresolved"].sum()),A=float(a["reached"].mean()),endpoint_fraction=float(a["endpoint"].mean()),solid_angle_sr=float(4*np.pi*a["reached"].mean()),**opening(a["reached"],ul,g)))
  masks=np.array(masks);endmasks=np.array(endmasks)
  assert not np.any(masks[1:]&~masks[:-1]),"Fixed-direction tangential nesting violation"
  assert not np.any(endmasks[1:]&~endmasks[:-1])
  np.savez_compressed(d/(s["site_id"]+"_directions.npz"),directions_chemical=ul,frame=np.array(s["frame"]),radii_A=RADII,accessible=masks,endpoint=endmasks,unresolved=unresolved)
  ar=np.array([x["A"] for x in rr]);passed=RADII[ar>0]
  rec=dict(material=mat,site_id=s["site_id"],parent_indices=s["parent_indices"],directions=n,in_legacy_ge70_cohort=s["site_id"] in legacy,prescreen_endpoint_fraction=float(pre[i]),S_log=score(RADII,ar),S_linear=float(np.trapz(ar,RADII)/(RADII[-1]-RADII[0])),sampled_Rmax_A=float(passed.max()) if len(passed) else None,Rmax_right_censored=bool(len(passed) and passed.max()==RADII[-1]),radii=rr,topology_scope="MODEL_PARENT_BOND; LOOSE_TAIL_AND_CAP_CLASSIFICATION_NOT_EVALUATED",sampling_status="NO_ACCESS_DETECTED_AT_SAMPLED_DIRECTIONS" if not len(passed) else "SAMPLED_GEOMETRIC_ACCESS")
  results.append(rec);old.emit("probe_site",material=mat,done=ii+1,total=len(selected),site=s["site_id"],S_log=rec["S_log"],Rmax=rec["sampled_Rmax_A"])
 # Refine three deterministic nonzero-access cases spanning score range; use zero if none.
 order=np.argsort([v["S_log"] for v in results],kind="stable");nonzero=[j for j in order if results[j]["S_log"]>0]
 refs=np.array(nonzero if nonzero else list(order));pick=refs[np.linspace(0,len(refs)-1,min(3,len(refs)),dtype=int)]
 refinement=[]
 for j in pick:
  s=rows[selected[j]];u2=old.fibonacci(2*n)@np.array(s["frame"]).T
  arr=np.array([trace(ce.env,np.array(s["center"]),u2,r)["reached"].mean() for r in RADII])
  refinement.append(dict(site_id=s["site_id"],directions=2*n,A=arr.tolist(),max_absolute_A_change=float(np.max(abs(arr-np.array([x["A"] for x in results[j]["radii"]])))),S_log=score(RADII,arr)))
 doc=dict(material=mat,model_parent_bond_count=len(rows),pilot_count=len(selected),selected_original_indices=selected,source_atom_sha256=old.sha(R/"inputs_v2"/(mat+"_parent_atoms.npz")),chemical_mapping=ce.audit,results=results,refinement=refinement)
 with (d/"summary.json").open("x") as f:json.dump(doc,f,indent=2)
 return doc
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--output",required=True);ap.add_argument("--pilot",type=int,default=12);ap.add_argument("--directions",type=int,default=1024);args=ap.parse_args()
 out=R/args.output;assert out.parent==R;out.mkdir(exist_ok=False);t=time.time()
 with (R/"RUN_LOG.jsonl").open("a") as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),event="tangent_probe_start",command=sys.argv,source_sha256=old.sha(__file__),output=str(out)))+"\n")
 mats=[material(mat,out,args.directions,args.pilot) for mat in ["PET","PA6","PA66"]]
 doc=dict(status="PILOT_COMPUTED_REQUIRES_AUDIT",materials=mats,seconds=time.time()-t,definitions=dict(endpoint="carbonyl C center+(1.7+r)*u; target atoms retained",path="exact collision test of straight swept spheres; <=2A segments, max120A; external plane beyond ALL parent vdW atoms; XY periodic; unresolved never pass",score="S_log=integral A(r) dlog(r)/log(10/0.5); empirical geometric score0..1, no density term",direction_scope="4pi full sphere; C=O x, leaving-atom plane y, cross-product z; not restricted to attack angles",opening="connected sampled spherical Delaunay graph, not continuous cone proof",limitations=["single fixed snapshot","finite angular and radius sampling","straight paths only","large isotropic sphere is not an enzyme core","parent model caps/loose tails not filtered in pilot","not full-population result; 64-direction endpoint-only screen is not final exposure","no new triad-donor scan"]),radii_A=RADII.tolist())
 with (out/"summary.json").open("x") as f:json.dump(doc,f,indent=2)
 with (R/"RUN_LOG.jsonl").open("a") as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),event="tangent_probe_complete",output=str(out),exit_status=0,seconds=time.time()-t))+"\n")
 with (R/"RUNBOOK.md").open("a") as f:f.write("\nTangent probe pilot: "+str(out)+". CPU existing Python tangent_probe_v1.py --output NEW_UNIQUE_NAME --pilot 12 --directions 1024. All-atom exact segment collision, r0.5..10A, 120A path cap; preserve unresolved/censored flags and no density in score. Read summary definitions before any biological interpretation.\n")
 old.emit("tangent_probe_complete",seconds=time.time()-t,output=str(out))
if __name__=="__main__":main()
