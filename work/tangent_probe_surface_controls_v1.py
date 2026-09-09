import sys,pathlib,json,datetime,time,numpy as np
R=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(R))
import tangent_probe_v1 as p
OUT=R/"tangent_probe_surface_controls_20260909_v1";OUT.mkdir(exist_ok=False)
start=time.time();docs=[]
for mat in ["PET","PA6","PA66"]:
 ce=p.chem.ChemEnv(mat);cat=p.catalogue(mat,ce)
 bycarbon={s["site_id"].rsplit("_",1)[0]:s for s in cat}
 oldrows=json.load(open(R/"polymer_exposure_full_v2"/mat/"sites.json"))
 exposed=sorted([s for s in oldrows if s["external_shell_fraction"]>0],key=lambda s:(-s["external_shell_fraction"],s["site_id"]))
 chosen=[exposed[i] for i in np.linspace(0,len(exposed)-1,min(12,len(exposed)),dtype=int)]
 d=OUT/mat;d.mkdir();ul=p.old.fibonacci(1024);graph=p.graph(ul);records=[]
 for j,old in enumerate(chosen):
  s=bycarbon[old["site_id"].rsplit("_",1)[0]]
  delta=np.array(s["center"])-old["center"];delta[:2]-=ce.env.box[:2]*np.round(delta[:2]/ce.env.box[:2])
  assert np.linalg.norm(delta)<1e-5
  u=ul@np.array(s["frame"]).T;rr=[];ma=[];ep=[];un=[]
  for r in p.RADII:
   result=p.trace(ce.env,np.array(s["center"]),u,r)
   ma.append(result["reached"]);ep.append(result["endpoint"]);un.append(result["unresolved"])
   rr.append(dict(radius_A=float(r),A=float(result["reached"].mean()),accessible_n=int(result["reached"].sum()),endpoint_n=int(result["endpoint"].sum()),unresolved_n=int(result["unresolved"].sum()),**p.opening(result["reached"],ul,graph)))
  ma=np.array(ma);ep=np.array(ep)
  assert not np.any(ma[1:]&~ma[:-1]);assert not np.any(ep[1:]&~ep[:-1])
  ar=np.array([x["A"] for x in rr]);rs=p.RADII[ar>0]
  record=dict(material=mat,site_id=s["site_id"],legacy_site_id=old["site_id"],parent_indices=s["parent_indices"],carbon_identity_verified=True,old_water_exposure_for_control_selection_only=old["external_shell_fraction"],score_uses_water_or_density=False,S_log=p.score(p.RADII,ar),sampled_Rmax_A=float(rs.max()) if len(rs) else None,Rmax_right_censored=bool(len(rs) and rs.max()==10),directions=1024,radii=rr)
  np.savez_compressed(d/(s["site_id"]+"_directions.npz"),directions_chemical=ul,frame=s["frame"],radii_A=p.RADII,accessible=ma,endpoint=ep,unresolved=un)
  records.append(record);p.old.emit("surface_control",material=mat,done=j+1,total=len(chosen),site=s["site_id"],legacy_id=old["site_id"],score=record["S_log"],Rmax=record["sampled_Rmax_A"])
 # Refine the most open control, with the same 120A straight path constraint.
 best=max(records,key=lambda x:x["S_log"]);s=bycarbon[best["site_id"].rsplit("_",1)[0]];u=p.old.fibonacci(2048)@np.array(s["frame"]).T
 arr=np.array([p.trace(ce.env,np.array(s["center"]),u,r)["reached"].mean() for r in p.RADII])
 doc=dict(material=mat,control_count=len(records),selection="12 evenly spaced ranks among previously externally exposed sites; legacy support-selected positive controls, not unbiased population",results=records,refinement=dict(site_id=s["site_id"],directions=2048,A=arr.tolist(),S_log=p.score(p.RADII,arr),max_absolute_A_change=float(np.max(abs(arr-np.array([v["A"] for v in best["radii"]]))))))
 with (d/"summary.json").open("x") as f:json.dump(doc,f,indent=2)
 docs.append(doc)
with (OUT/"summary.json").open("x") as f:json.dump(dict(status="SURFACE_CONTROLS_COMPLETE_REQUIRES_AUDIT",materials=docs,seconds=time.time()-start),f,indent=2)
with (R/"RUN_LOG.jsonl").open("a") as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(pathlib.Path(__file__)),source_sha256=p.old.sha(__file__),output=str(OUT),exit_status=0,seconds=time.time()-start))+"\n")
print("SURFACE_CONTROLS_COMPLETE",flush=True)
