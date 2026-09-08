#!/usr/bin/env python3
"""CPU pilot: measured multiscale material shape matching; not enzyme design."""
import sys,json,math,argparse,time,datetime,hashlib
from pathlib import Path
import numpy as np
from scipy.spatial import cKDTree
from scipy.ndimage import gaussian_filter,map_coordinates
R=Path(__file__).resolve().parent
sys.path.insert(0,str(R))
import chemical_contact_surface_v3 as m
AX=np.arange(-22.,23.)
GRID=np.stack(np.meshgrid(AX,AX,AX,indexing="ij"),-1).reshape(-1,3)
def wquantile(x,w,q):
 ix=np.argsort(x,kind="stable");x=np.asarray(x)[ix];w=np.asarray(w)[ix]
 return float(x[min(np.searchsorted(np.cumsum(w),q*w.sum(),side="left"),len(x)-1)])
def mesh_rep(v,f,refine=1):
 if len(f)==0:return None
 tri=np.asarray(v)[f];a=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
 if refine>1:
  # Split each triangle into four, preserving its exact area and surface.
  for _ in range(refine-1):
   p,q,r=tri[:,0],tri[:,1],tri[:,2];pq=(p+q)/2;qr=(q+r)/2;rp=(r+p)/2
   tri=np.concatenate([np.stack([p,pq,rp],1),np.stack([pq,q,qr],1),np.stack([rp,qr,r],1),np.stack([pq,qr,rp],1)])
  a=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)/2
 centers=tri.mean(1)
 samples=np.concatenate([tri[:,0],tri[:,1],tri[:,2],(tri[:,0]+tri[:,1])/2,(tri[:,1]+tri[:,2])/2,(tri[:,2]+tri[:,0])/2,centers])
 return dict(centers=centers,weights=a,tree=cKDTree(samples),area=float(a.sum()))
def distance(a,b):
 if a is None or b is None:return float("inf"),float("inf")
 da=b["tree"].query(a["centers"])[0];db=a["tree"].query(b["centers"])[0]
 return max(wquantile(da,a["weights"],.95),wquantile(db,b["weights"],.95)),float(max(da.max(),db.max()))
def required_tolerance(errors,fraction):
 return float(np.sort(errors)[int(math.ceil(len(errors)*fraction))-1])
def jsonwrite(p,d):
 with p.open("x") as f:json.dump(d,f,indent=2,allow_nan=False)
def finite(x):
 if isinstance(x,(float,np.floating)) and not np.isfinite(x):return None
 if isinstance(x,dict):return {str(k):finite(v) for k,v in x.items()}
 if isinstance(x,(list,tuple)):return [finite(v) for v in x]
 if isinstance(x,np.ndarray):return finite(x.tolist())
 if isinstance(x,np.integer):return int(x)
 return x
def density_grid(ce,g):
 xyz=ce.data["xyz_A"].copy();xyz[:,:2]%=ce.env.box[:2]
 edges=[g["origin"][k]-.5*g["step"][k]+np.arange(g["shape"][k]+1)*g["step"][k] for k in range(3)]
 hist=np.histogramdd(xyz,bins=edges,weights=ce.data["masses_Da"])[0]/np.prod(g["step"])
 rho=gaussian_filter(hist,3/g["step"],mode=("wrap","wrap","nearest"))
 z=g["origin"][2]+np.arange(g["shape"][2])*g["step"][2]
 central=abs(z-np.median(xyz[:,2]))<10
 threshold=.5*float(np.median(rho[:,:,central]))
 assert threshold>0
 return rho,threshold
DISK=np.array([(x,y,-5.) for x in np.linspace(-10,10,11) for y in np.linspace(-10,10,11) if x*x+y*y<=100])
def support_at(v,s,g,ce,rho,threshold):
 frame=np.array(s["surface_frame"])
 points=np.array(s["center"])+(v[:,None,:]+DISK[None,:,:])@frame.T
 coords=m.old.gridcoords(points.reshape(-1,3),g,ce.env)
 val=map_coordinates(rho,coords,order=1,mode="grid-wrap").reshape(len(v),-1)
 return (val>=threshold).mean(1)
def connected(v,f,n):
 if not len(f):return np.empty((0,3)),np.empty((0,3),int),np.empty((0,3))
 from scipy.sparse import coo_matrix
 from scipy.sparse.csgraph import connected_components
 ii=np.r_[f[:,0],f[:,1],f[:,2]];jj=np.r_[f[:,1],f[:,2],f[:,0]]
 _,lab=connected_components(coo_matrix((np.ones(len(ii)),(ii,jj)),shape=(len(v),len(v))),directed=False)
 choices=[]
 for k in np.unique(lab[f[:,0]]):
  ff=f[lab[f[:,0]]==k];ids=np.unique(ff)
  area=np.linalg.norm(np.cross(v[ff[:,1]]-v[ff[:,0]],v[ff[:,2]]-v[ff[:,0]]),axis=1).sum()/2
  choices.append(((np.linalg.norm(v[ids],axis=1).min(),-area),ff))
 f=min(choices,key=lambda a:a[0])[1];ids=np.unique(f);lookup=np.full(len(v),-1);lookup[ids]=np.arange(len(ids))
 return v[ids],lookup[f],n[ids]
def make_mesh(field,radius,validity=None):
 v,f,n=m.mesh_patch(field,radius=radius)
 if validity is not None and len(f):f=f[np.all(validity(v)[f],axis=1)];v,f,n=connected(v,f,n)
 return v,f,n
def main():
 ap=argparse.ArgumentParser();ap.add_argument("--material",default="PET",choices=["PET","PA6","PA66"]);ap.add_argument("--pilot",type=int,default=8);ap.add_argument("--output",required=True)
 a=ap.parse_args();out=R/a.output;assert out.parent==R;out.mkdir(exist_ok=False)
 started=time.time();log=dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__)),source_sha256=m.old.sha(__file__),command=sys.argv,output=str(out))
 with (R/"RUN_LOG.jsonl").open("a") as f:f.write(json.dumps(dict(log,event="multiscale_start"))+"\n")
 rows=json.load(open(R/"polymer_exposure_full_v2"/a.material/"sites.json"))
 panelrows=json.load(open(R/"panel_surface_scan_v1/sites.json"));cov=np.load(R/"panel_surface_scan_v1/coverage.npz")["coverage"]
 hits={s["site_id"] for i,s in enumerate(panelrows) if s["surface"]==a.material+"_original" and cov[:,i].any()}
 alln=len(rows)
 if a.pilot:
  # Deterministic, evenly spread sites in each stratum, not favorable-site picking.
  eligible=[s for s in rows if s["site_id"] in hits];other=[s for s in rows if s["site_id"] not in hits]
  rows=[ss[i] for ss,k in [(eligible,a.pilot//2),(other,a.pilot-a.pilot//2)] for i in np.linspace(0,len(ss)-1,min(k,len(ss)),dtype=int)]
 assert len({s["site_id"] for s in rows})==len(rows)
 jsonwrite(out/"sites.json",rows)
 ce=m.ChemEnv(a.material)
 m.old.emit("multiscale_grid_start",material=a.material)
 g=m.ses_grid(ce)
 rho,threshold=density_grid(ce,g)
 m.AX=AX;m.GRID=GRID
 fields=[];meshes={r:[] for r in [10,15,20]};reps={r:[] for r in meshes};audits=[]
 for i,s in enumerate(rows):
  world=np.array(s["center"])+GRID@np.array(s["surface_frame"]).T
  field=m.interp(g["phi"],world,g,ce.env).reshape((len(AX),)*3);fields.append(field.astype(np.float16))
  rec=dict(site_id=s["site_id"],normal_status=s["normal_status"],radii={})
  for radius in meshes:
   v,f,n=m.mesh_patch(field,radius=radius);before=float(m.vertex_areas(v,f).sum()) if len(f) else 0.
   if s["normal_status"]!="PASS":v,f,n=np.empty((0,3)),np.empty((0,3),int),np.empty((0,3))
   elif len(f):
    sup=support_at(v,s,g,ce,rho,threshold)
    f=f[np.all(sup[f]>=.7,axis=1)];v,f,n=connected(v,f,n)
   meshes[radius].append((v,f,n));reps[radius].append(mesh_rep(v,f))
   area=float(m.vertex_areas(v,f).sum()) if len(f) else 0.
   rec["radii"][radius]=dict(raw_area_A2=before,dense_supported_area_A2=area,retained_fraction=area/before if before else None,status="MESH_AVAILABLE" if len(f) else "NOT_EVALUATED_NO_DENSE_PATCH")
  audits.append(rec);m.old.emit("multiscale_site",material=a.material,done=i+1,total=len(rows))
 fields=np.asarray(fields)
 np.savez_compressed(out/"fields_R22.npz",axis_A=AX,fields=fields)
 jsonwrite(out/"mesh_scope_audit.json",finite(audits))
 groups={"all":list(range(len(rows))),"current_core_accessible":[i for i,s in enumerate(rows) if s["site_id"] in hits]}
 results=[];payloads=[]
 for group,ids in groups.items():
  if not ids:continue
  # Pilot candidate bank: each observed real field and unflattened pointwise median.
  bank=[("real_"+rows[i]["site_id"],fields[i].astype(float),[i]) for i in ids if rows[i]["normal_status"]=="PASS"]
  bank.append(("median_all",np.median(fields[ids].astype(float),axis=0),ids))
  candidates=[]
  for name,field,src in bank:
   def valid(v):
    val=np.mean([support_at(v,rows[i],g,ce,rho,threshold)>=.7 for i in src],axis=0)
    return val>=.5
   cm={r:make_mesh(field,r,valid) for r in meshes}
   cr={r:mesh_rep(*cm[r][:2]) for r in meshes}
   errors={r:[distance(cr[r],reps[r][i]) for i in ids] for r in meshes}
   candidates.append(dict(name=name,field=field,mesh=cm,errors=errors,source=src))
  for fraction in [.9,.7,.5]:
   best=min(candidates,key=lambda c:(required_tolerance([e[0] for e in c["errors"][20]],fraction),np.mean([e[0] for e in c["errors"][20]]),c["name"]))
   err20=np.array(best["errors"][20])[:,0];k=math.ceil(len(ids)*fraction)
   chosen_local=np.argsort(err20,kind="stable")[:k];chosen=[ids[j] for j in chosen_local]
   for radius in meshes:
    v,f,n=best["mesh"][radius];errors=np.array(best["errors"][radius]);tol=float(errors[chosen_local,0].max())
    actual=np.where(errors[:,0]<=tol+1e-9)[0] if np.isfinite(tol) else np.array([],int)
    key=f"{a.material}_{group}_coverage{int(fraction*100)}_R{radius}"
    d=out/key;d.mkdir()
    rec=dict(material=a.material,cohort=group,radius_A=radius,target_fraction=fraction,denominator=len(ids),required_n=k,candidate=best["name"],chosen_site_ids=[rows[i]["site_id"] for i in chosen],tolerance_for_same_subset_A=tol,actual_n_at_tolerance=len(actual),actual_site_ids=[rows[ids[j]]["site_id"] for j in actual],coverage_at_fixed_A={str(t):int((errors[:,0]<=t).sum()) for t in [1,2,3]},per_site_errors=[dict(site_id=rows[i]["site_id"],D95_A=float(errors[j,0]),sampled_max_A=float(errors[j,1])) for j,i in enumerate(ids)],status="PILOT_DESCRIPTIVE_IN_SAMPLE" if a.pilot else "DESCRIPTIVE_IN_SAMPLE",candidate_bank_size=len(bank),mesh_status="AVAILABLE" if len(f) else "NOT_EVALUATED")
    if len(f) and len(actual):
     actualrows=[rows[ids[j]] for j in actual]
     props,sd,sup,valid=m.chem_on_mesh(v,actualrows,lambda s:np.array(s["surface_frame"]),ce,g)
     stats,pay=m.export_mesh(d,key,v,f,n,props,sd,sup,valid);rec["mesh_stats"]=stats
     channels=[pay]
     for channel,idx,target in [("HBA",1,np.array([.15,.55,.85])),("HBD",2,np.array([.7,.25,.7]))]:
      col=.93+(target-.93)*props[:,idx,None];col[~valid]=.55
      channels.append(dict(pay,name=key+"_"+channel,c=col.round(4).tolist()))
     m.viewer(d/(key+"_view.py"),channels,default=key)
     payloads.append(dict(key=key,meshes=channels))
    jsonwrite(d/"summary.json",finite(rec));results.append(rec)
    m.old.emit("multiscale_fit",key=key,tolerance_A=finite(tol),actual=len(actual),denominator=len(ids))
 # Quadrature refinement on first two real patches, never interpreted as grid convergence.
 refine={}
 for radius in meshes:
  if len(rows)>1:
   x,y=meshes[radius][:2]
   coarse=distance(mesh_rep(*x[:2]),mesh_rep(*y[:2]))[0]
   fine=distance(mesh_rep(*x[:2],refine=2),mesh_rep(*y[:2],refine=2))[0]
   refine[radius]=dict(coarse_D95_A=coarse,refined_D95_A=fine,absolute_change_A=abs(coarse-fine))
 doc=dict(status="PILOT_COMPLETE_REQUIRES_SCIENTIFIC_REVIEW" if a.pilot else "COMPUTED_REQUIRES_AUDIT",material=a.material,pilot=bool(a.pilot),parent_denominator=alln,evaluated_sites=len(rows),radii_A=[10,15,20],chemical_typing=ce.audit,definitions=dict(geometry="original parent exterior-connected approximate water SES; 0.5 A global grid;1 A local mesh",bulk_filter="per-vertex >=70% density support on 81-point R10 disk 5A inward along fixed local outward normal; Gaussian3A rho, threshold half median central slab density; nearest carbonyl connected component; proxy, not chain topology removal",frame="density outward normal + projected C=O; same material frame across scales; not core-specific rotation",D95="maximum of two area-weighted 95th-percentile nearest sampled-mesh distances; triangle centroids weighted by area; vertices/midpoints/centroids target samples",coverage="fixed best-candidate matched subset selected at R20, same IDs for R10/R15; measured tolerance is max D95 in that subset; actual matches at tolerance may include extra sites",candidate_bank="real site fields + pointwise median, choose lowest required tolerance at R20; finite-library in-sample result, not global optimum",chemistry="post-hoc Crippen fragment index/HBA/HBD over actually matched sites; not electrostatic potential or binding energy"),quadrature_refinement=refine,results=results,seconds=time.time()-started,source_sha256=m.old.sha(__file__),limitations=["No new MD or donor sampling","Current selected 101-core panel, not full historical library","Original frame uncertainty preserved","No enzyme shape compatibility inferred from material matching","Dense support is approximate; visually review tails/sidewall scope","Local grid convergence NOT_EVALUATED","Independent validation NOT_EVALUATED"])
 jsonwrite(out/"summary.json",finite(doc));jsonwrite(out/"viewer_payloads.json",payloads)
 with (R/"RUN_LOG.jsonl").open("a") as f:f.write(json.dumps(dict(log,event="multiscale_complete",exit_status=0,seconds=time.time()-started))+"\n")
 with (R/"RUNBOOK.md").open("a") as f:f.write("\nMultiscale measured matching pilot: "+str(out)+". Reproduce with CPU Python multiscale_material_match_v1.py --material "+a.material+" --pilot "+str(a.pilot)+" --output NEW_UNIQUE_NAME. Existing 0.5A SES, new ±22A local field, dense-support-filtered mesh; D95 not enzyme collision; see summary definitions.\n")
 m.old.emit("multiscale_complete",output=str(out),seconds=time.time()-started)
if __name__=="__main__":main()
