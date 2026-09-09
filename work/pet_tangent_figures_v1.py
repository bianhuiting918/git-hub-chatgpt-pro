#!/usr/bin/env python3
"""PET-only postprocessing; immutable tangent scan; no catalyst scan."""
import json,time,datetime,sys
from pathlib import Path
import numpy as np
from scipy.ndimage import map_coordinates
from scipy.sparse.csgraph import connected_components
import tangent_probe_v1 as t
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
SRC=R/"tangent_probe_all_parent_20260909_v1"/"PET"
OUT=R/"pet_tangent_figures_20260909_v1"
BOUNDS=[.0025,.01,.02,np.inf]
LABELS=["0 < S <= 0.0025","0.0025 < S <= 0.01","0.01 < S <= 0.02","S > 0.02"]
COLORS=["#516b99","#229aab","#bd8b2a","#bd5346"]
def bin_index(s):
 if s<=0:return None
 return int(np.searchsorted(BOUNDS,s,side="left"))
def frame_from_opening(z,co):
 z=np.array(z,float);z/=np.linalg.norm(z)
 x=np.array(co,float)-z*np.dot(co,z)
 if np.linalg.norm(x)<1e-8:
  v=np.eye(3)[np.argmin(abs(z))];x=v-z*np.dot(v,z)
 x/=np.linalg.norm(x)
 return np.column_stack([x,np.cross(z,x),z])
def cap_kind(mol,ci,li):
 c=mol.GetAtomWithIdx(int(ci));le=mol.GetAtomWithIdx(int(li))
 if not any(a.GetSymbol()=="C" and a.GetIsAromatic() for a in c.GetNeighbors()):return "ACYL_CAP"
 for a in le.GetNeighbors():
  if a.GetIdx()!=ci and a.GetSymbol()=="C" and sum(x.GetAtomicNum()>1 for x in a.GetNeighbors())==1:return "ALKYL_CAP"
 return "BACKBONE"
def make_grid(env,h):
 nxy=np.ceil(env.box[:2]/h).astype(int);lo=np.floor(env.xyz[:,2].min())-5
 nz=int(np.ceil((env.xyz[:,2].max()+5-lo)/h))
 shape=tuple(int(x) for x in [*nxy,nz]);step=np.array([*(env.box[:2]/nxy),h]);origin=np.array([0.,0.,lo])+.5*step
 free=np.empty(np.prod(shape),bool)
 for start in range(0,len(free),100000):
  idx=np.column_stack(np.unravel_index(np.arange(start,min(start+100000,len(free))),shape))
  free[start:start+len(idx)]=env.clearance(origin+idx*step)>=.5-1e-7
 ext=t.old.external_components(free.reshape(shape))
 return dict(origin=origin,step=step,shape=shape,external=ext)
def interp_ext(grid,p,env):
 q=np.array(p,float);q[:,:2]%=env.box[:2]
 c=((q-grid["origin"])/grid["step"]).T
 outside=(c[2]<0)|(c[2]>grid["shape"][2]-1)
 c[2]=np.clip(c[2],0,grid["shape"][2]-1)
 v=map_coordinates(grid["external"].astype(np.float32),c,order=1,mode="grid-wrap")
 v[outside]=1.
 return v
def read_mask(s):
 return np.load(SRC/(s["site_id"]+"_directions.npz"))
def main():
 start=time.time();OUT.mkdir(exist_ok=False)
 summary=json.load(open(SRC/"summary.json"));rows=json.load(open(SRC/"all_parent_bonds.json"))
 selected=json.load(open(SRC/"selected_sites.json"));results=summary["results"]
 assert len(rows)==len(selected)==len(results)==8400
 byid={x["site_id"]:x for x in rows};assert len(byid)==8400
 assert set(byid)=={x["site_id"] for x in results}=={x["site_id"] for x in selected}
 ce=t.chem.ChemEnv("PET");assert t.old.sha(R/"inputs_v2/PET_parent_atoms.npz")==summary["source_atom_sha256"]
 mol=t.chem.Chem.SDMolSupplier(str(t.chem.sdf_for("PET")),removeHs=False)[0];n=mol.GetNumAtoms()
 ul=t.old.fibonacci(1024);g=t.graph(ul);assignments=[];positive=[];unresolved_sites=0;uncertain=[]
 for ix,s in enumerate(results):
  a=read_mask(s);mask=a["accessible"];end=a["endpoint"];un=a["unresolved"]
  assert mask.shape==end.shape==un.shape==(10,1024)
  assert np.allclose(a["directions_chemical"],ul) and np.array_equal(a["radii_A"],t.RADII)
  assert not np.any(mask&~end) and not np.any(mask&un)
  assert not np.any(mask[1:]&~mask[:-1]) and not np.any(end[1:]&~end[:-1])
  assert np.allclose(a["frame"],byid[s["site_id"]]["frame"])
  ar=mask.mean(1);assert abs(t.score(t.RADII,ar)-s["S_log"])<1e-12
  for j,r in enumerate(s["radii"]):
   assert int(mask[j].sum())==r["accessible_n"] and int(end[j].sum())==r["endpoint_n"] and int(un[j].sum())==r["unresolved_n"]
  unresolved_sites+=int(un.any())
  upper=t.score(t.RADII,(mask|un).mean(1))
  ci,oi,li=s["parent_indices"];kind=cap_kind(mol,ci%n,li%n)
  bi=bin_index(s["S_log"]);boundary=bi!=bin_index(upper)
  if boundary:uncertain.append(s["site_id"])
  rec=dict(site_id=s["site_id"],chain=byid[s["site_id"]]["chain"],kind=kind,S_log=s["S_log"],S_upper_unresolved=upper,bin=bi,bin_connectivity_uncertain=boundary,has_unresolved=bool(un.any()),sampled_Rmax_A=s["sampled_Rmax_A"],A=ar.tolist(),loose_tail_status="NOT_EVALUATED")
  assignments.append(rec)
  if bi is not None and kind=="BACKBONE" and not boundary:
   ids=np.flatnonzero(mask[0]);nc,lab=connected_components(g[ids][:,ids],directed=False);counts=np.bincount(lab);best=ids[lab==counts.argmax()]
   z=ul[best].mean(0);z/=np.linalg.norm(z);z=ul[best[np.argmax(ul[best]@z)]]
   chemical=np.array(byid[s["site_id"]]["frame"])
   F=frame_from_opening(chemical@z,chemical[:,0])
   assert np.allclose(F.T@F,np.eye(3)) and np.linalg.det(F)>.99999
   rec["surface_frame"]=F.tolist();rec["center"]=byid[s["site_id"]]["center"];positive.append(rec)
  if ix%2000==0:t.old.emit("audit",done=ix+1,total=len(results))
 # Verify unchanged numerical results against all previous PET pilot/control records.
 comparisons=0
 for folder in ["tangent_probe_pilot_20260909_v1","tangent_probe_surface_controls_20260909_v1"]:
  p=R/folder/"PET"/"summary.json"
  for s in json.load(open(p))["results"]:
   prior=np.load(p.parent/(s["site_id"]+"_directions.npz"));now=read_mask(s)
   for k in ["accessible","endpoint","unresolved"]:assert np.array_equal(prior[k],now[k])
   comparisons+=1
 t.old.emit("audit_pass",records=8400,prior_comparisons=comparisons,fit_sites=len(positive))
 axis=np.arange(-15,15.001,.5);xx,zz=np.meshgrid(axis,axis,indexing="ij")
 plane0=np.stack([xx,np.zeros_like(xx),zz],-1);plane1=np.stack([np.zeros_like(xx),xx,zz],-1)
 local=np.stack([plane0,plane1]).reshape(-1,3);circle=xx**2+zz**2<=225
 t.old.emit("surface_grid_start",spacing=.5,probe=.5)
 grid=make_grid(ce.env,.5)
 # Cache float once, to avoid repeated allocation of the full global grid.
 grid["external"]=grid["external"].astype(np.float32)
 t.old.emit("surface_grid_ready",shape=grid["shape"])
 fields=np.empty((len(positive),2,len(axis),len(axis)),np.float32)
 for i,s in enumerate(positive):
  world=np.array(s["center"])+local@np.array(s["surface_frame"]).T
  fields[i]=(1-interp_ext(grid,world,ce.env)).reshape(2,len(axis),len(axis))
  if i%200==0:t.old.emit("section_fields",done=i+1,total=len(positive))
 # Independent coarser geometry diagnostic (does not recompute probe scores).
 diagnostic_ids=[]
 for b in range(4):
  ids=[i for i,s in enumerate(positive) if s["bin"]==b]
  diagnostic_ids.extend(ids[j] for j in np.linspace(0,len(ids)-1,min(2,len(ids)),dtype=int))
 coarse=make_grid(ce.env,1.);coarse["external"]=coarse["external"].astype(np.float32);grid_checks=[]
 for i in diagnostic_ids:
  s=positive[i];world=np.array(s["center"])+local@np.array(s["surface_frame"]).T
  p=(1-interp_ext(coarse,world,ce.env)).reshape(fields[i].shape)
  grid_checks.append(dict(site_id=s["site_id"],mean_absolute_occupancy_difference=float(abs(p[:,circle]-fields[i][:,circle]).mean()),binary_disagreement_fraction=float(((p[:,circle]>=.5)!=(fields[i][:,circle]>=.5)).mean())))
 plt.rcParams.update({"font.size":10,"svg.fonttype":"none","axes.spines.top":False,"axes.spines.right":False})
 fig,axs=plt.subplots(2,4,figsize=(14,7.5),sharex=True,sharey=True,constrained_layout=True);groups=[]
 for b in range(4):
  ids=[i for i,s in enumerate(positive) if s["bin"]==b];arr=fields[ids];mean=arr.mean(0)
  distance=np.mean(abs(arr[:,:,circle]-mean[:,circle]),axis=(1,2));rep=ids[int(np.argmin(distance))]
  av=np.array([positive[i]["A"] for i in ids]);rm=np.array([positive[i]["sampled_Rmax_A"] for i in ids])
  group=dict(bin=b,label=LABELS[b],n=len(ids),chains=len({positive[i]["chain"] for i in ids}),representative=positive[rep]["site_id"],A_median=np.median(av,axis=0).tolist(),A_q25=np.quantile(av,.25,axis=0).tolist(),A_q75=np.quantile(av,.75,axis=0).tolist(),Rmax_median=float(np.median(rm)),Rmax_ge2=int((rm>=2).sum()),Rmax_ge4=int((rm>=4).sum()),Rmax_ge8=int((rm>=8).sum()),Rmax_censored10=int((rm>=10).sum()),field_heterogeneity=float(distance.mean()))
  groups.append(group)
  np.savez_compressed(OUT/f"group_{b}_sections.npz",axis_A=axis,occupancy_probability=mean,site_ids=np.array([positive[i]["site_id"] for i in ids]),representative_field=fields[rep],circle_mask=circle)
  for k in range(2):
   ax=axs[k,b];p=mean[k].copy();p[~circle]=np.nan
   ax.contourf(xx,zz,p,levels=[.25,.75,1.001],colors=["#f0e8d5","#b9a784"],alpha=.65)
   for lev,ls,lw in [(.25,":",1),(.5,"-",2),(.75,"--",1)]:
    if np.nanmin(p)<lev<np.nanmax(p):ax.contour(xx,zz,p,levels=[lev],colors=[COLORS[b]],linestyles=[ls],linewidths=[lw])
   pr=fields[rep,k].copy();pr[~circle]=np.nan
   ax.contour(xx,zz,pr,levels=[.5],colors=["#737373"],linewidths=.65,alpha=.65)
   ax.plot(0,0,"ko",ms=4);ax.annotate("C",(0,0),xytext=(4,-12),textcoords="offset points",fontsize=8)
   ax.set_aspect("equal");ax.set(xlim=(-15,15),ylim=(-15,15),xticks=[-15,0,15],yticks=[-15,0,15])
   ax.axvline(0,color=".85",lw=.5,zorder=0);ax.axhline(0,color=".85",lw=.5,zorder=0)
   if k==0:ax.set_title(LABELS[b]+f"\nn = {len(ids)}",fontsize=11)
   ax.set_xlabel(("x" if k==0 else "y")+" (A)")
   if b==0:ax.set_ylabel("z: aligned accessible direction (A)")
 fig.suptitle("PET backbone | carbonyl-centered interface by geometric exposure\nTop: x-z section; bottom: y-z section. Original slab top/bottom pooled.",fontsize=14)
 fig.supxlabel("Colored contours: blocked-space probability 0.25 / 0.50 / 0.75; gray: actual representative.\n0.5 A probe-center exclusion envelope, not atomic SES. C center = origin; patch radius 15 A. Loose tails not independently removed.",fontsize=10)
 for ext in ["png","svg"]:fig.savefig(OUT/("pet_exposure_interfaces."+ext),dpi=160)
 plt.close(fig)
 fig,axs=plt.subplots(1,3,figsize=(14,4.7),constrained_layout=True)
 for b,gr in enumerate(groups):
  color=COLORS[b];axs[0].plot(t.RADII,100*np.array(gr["A_median"]),"-o",color=color,ms=3,label=f'{LABELS[b]} (n={gr["n"]})')
  axs[0].fill_between(t.RADII,100*np.array(gr["A_q25"]),100*np.array(gr["A_q75"]),color=color,alpha=.12)
  ids=[i for i,s in enumerate(positive) if s["bin"]==b]
  a=np.array([positive[i]["sampled_Rmax_A"] for i in ids])
  axs[1].plot(t.RADII,[100*np.mean(a>=r) for r in t.RADII],"-o",color=color,ms=3)
 axs[0].set(xlabel="Sphere radius (A)",ylabel="Accessible directions (% of full sphere)",title="A(r): median and site interquartile range",xlim=(.5,10),ylim=(0,None))
 axs[0].legend(fontsize=8,frameon=False)
 axs[1].set(xlabel="Sphere radius (A)",ylabel="Sites with >=1 accessible direction (%)",title="How many sites accommodate this size?",xlim=(.5,10),ylim=(0,101))
 kinds=["BACKBONE","ACYL_CAP","ALKYL_CAP"];tot=[sum(s["kind"]==k for s in assignments) for k in kinds];pos=[sum(s["kind"]==k and s["S_log"]>0 for s in assignments) for k in kinds]
 axs[2].bar(range(3),tot,color="#d7d7d7",label="No access detected");axs[2].bar(range(3),pos,color="#2692a5",label="Access detected")
 for j,(n,p) in enumerate(zip(tot,pos)):axs[2].text(j,n+120,f"{p}/{n}",ha="center",fontsize=10)
 axs[2].set(xticks=range(3),xticklabels=["Backbone","Acyl cap","Alkyl cap"],ylabel="Model bond count",title="Denominator audit: accessible / total",ylim=(0,max(tot)*1.18))
 axs[2].legend(fontsize=8,frameon=False)
 fig.suptitle("PET geometric accessibility | 8400 model ester bonds audited",fontsize=14)
 fig.supxlabel("S is a log-radius averaged geometric score, not support %, water exposure %, or activity.\nGrouping by S makes separation of A(r) expected, not independent biological validation. Sampling at 10 A is right-censored.",fontsize=10)
 for ext in ["png","svg"]:fig.savefig(OUT/("pet_exposure_statistics."+ext),dpi=160)
 plt.close(fig)
 for gr in groups:
  r=next(s for s in positive if s["site_id"]==gr["representative"]);gr["representative_S_log"]=r["S_log"]
 audit=dict(status="PET_FULL_NUMERICAL_AUDIT_PASS_DESCRIPTIVE_FIGURES",source=str(SRC),source_summary_sha256=t.old.sha(SRC/"summary.json"),atom_sha256=summary["source_atom_sha256"],script_sha256=t.old.sha(__file__),model_bonds=8400,accessible_model_bonds=sum(s["S_log"]>0 for s in assignments),no_access_detected=sum(s["S_log"]==0 for s in assignments),topology_counts={k:dict(total=n,accessible=p) for k,n,p in zip(kinds,tot,pos)},unresolved_sites=unresolved_sites,bin_connectivity_uncertain=uncertain,fit_sites=len(positive),groups=groups,prior_mask_comparisons=comparisons,refinement=summary["refinement"],surface_grid=dict(spacing_A=grid["step"].tolist(),probe_A=.5,shape=grid["shape"],connectivity="6-neighbor with XY periodic boundary and exterior top/bottom seeds"),surface_grid_checks=grid_checks,limitations=["single snapshot","loose tails not independently labeled or removed","caps excluded from main fits but retained in all counts","zero access is not proof of burial or no curved path","aligned z maximizes consistency of observed openings, not material normal","surface envelope uses curved grid connectivity; score uses straight continuous sphere paths","contours are per-voxel population frequencies, not joint whole-site coverage","finite-grid morphology is descriptive, not atomistic SES or enzyme shape","strata are reporting conventions, not validated biological thresholds"],seconds=time.time()-start)
 (OUT/"summary.json").write_text(json.dumps(audit,indent=2))
 (OUT/"site_assignments.json").write_text(json.dumps(assignments,separators=(",",":")))
 report=["# PET: full-catalogue geometric exposure and aligned interface sections","","Numerical audit complete. Figures are descriptive material geometry, not enzyme activity.","",f"All {len(results)} model bonds audited; {audit['accessible_model_bonds']} have sampled external access. Main fits use {len(positive)} backbone sites. Caps are separate; loose-tail classification remains NOT_EVALUATED.","","## Figures","","![Interfaces](pet_exposure_interfaces.png)","","![Statistics](pet_exposure_statistics.png)","","## Groups","","| S_log group | Sites | Median sampled Rmax (A) | Rmax >=4 A | Representative |","|---|---:|---:|---:|---|"]
 for gr in groups:report.append(f"| {gr['label']} | {gr['n']} | {gr['Rmax_median']} | {gr['Rmax_ge4']} | {gr['representative']} |")
 report+=["","## Interpretation","","Main shape panels pool slab top/bottom, translate carbonyl C to the origin and rotate each site's largest r=0.5 A accessible opening to +z. x is projected C=O. No mirror reflection. Two perpendicular cuts are shown to expose anisotropy.","","Shape contour values are blocked-space frequency 25/50/75% across each fixed exposure group. They are NOT exposure levels and NOT a guarantee that one enzyme fits that percentage of sites. Gray contour is one actual representative (minimum mean section-field deviation).","", "S_log = normalized integral of A(r) over log radius, r=0.5..10 A. The score, masks and direction counts are unchanged from the production scan. Display geometry is a 0.5 A probe-center exclusion envelope sampled at 0.5 A, NOT water exposure or atomic SES. The grid uses external connected void, so its curved-path envelope is not a substitute for the strict straight-path score.","","Small values of S do not directly mean small water-exposed surface area. No-access sampled sites are reported separately because an opening-based alignment cannot be defined for them without introducing another criterion.","","## Limitations",""]+["- "+x for x in audit["limitations"]]
 report+=["","Raw per-site masks remain on the server. GitHub publishes figures, summaries, source and assignments; no claim of a full NPZ backup.","","Remote output: "+str(OUT)]
 (OUT/"README.md").write_text("\n".join(report)+"\n")
 with (R/"RUN_LOG.jsonl").open("a") as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),event="pet_tangent_figures_complete",source_sha256=t.old.sha(__file__),output=str(OUT),exit_status=0,seconds=audit["seconds"]))+"\n")
 with (R/"RUNBOOK.md").open("a") as f:f.write("\nPET tangent full-audit/figures: existing CPU Python pet_tangent_figures_v1.py; output "+str(OUT)+". Use NEW OUT path to rerun; no edits to active all-material scan. Definitions/denominators/loose-tail limitations in README and summary.\n")
 t.old.emit("pet_figures_complete",output=str(OUT),fit_sites=len(positive),groups=groups,seconds=audit["seconds"])
if __name__=="__main__":main()
