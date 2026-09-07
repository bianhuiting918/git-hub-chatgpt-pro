"""All 473 saved joint candidates: descriptive coordinate distributions, no resampling."""
from pathlib import Path
import json,hashlib,datetime
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image
P=Path("/data/bht2/polymer_material_reference_20260827/simulation_slabs/pet_dp10_400chain_direct_v1")
R=P/"donor_geometry_batch_singleNH_v1/joint_transfer_map_v1"
O=R/"all473_distribution_v1"
O.mkdir(exist_ok=False)
inputs=["motifs.json","joint_transfer_results.npz","summary.json"]
before={f:hashlib.sha256((R/f).read_bytes()).hexdigest() for f in inputs}
z=np.load(R/"joint_transfer_results.npz")
M=z["motif_xyz_A"];H=z["motif_H_A"]
motifs=json.loads((R/"motifs.json").read_text())
assert M.shape==(473,36,3) and H.shape==(473,2,3) and len(motifs)==473
keys=["Ser OG","His NE2","Asp O (OD1)","NH donor N (pooled)"]
colors=["#239566","#4268b7","#a951a1","#d69224"]
# Keep the same explicitly labeled Asp atom as prior figures. This is OD1,
# not a claim that OD1 is always the His-contacting carboxylate oxygen.
idx=[]
for m in motifs:
    a=m["functional_indices"]; od=a["asp_o"][0]
    assert m["atom_names"][od]=="OD1",m["id"]
    idx.append([a["og"],a["his_ne2"],od])
idx=np.array(idx);ii=np.arange(473)
G=[M[ii,idx[:,j]] for j in range(3)]+[M[:,[26,32]].reshape(-1,3)]
assert [len(q) for q in G]==[473,473,473,946]
assert all(np.isfinite(q).all() for q in G)
def stats(q):
    return {a:{"min":float(np.min(q[:,j])),"q05":float(np.quantile(q[:,j],.05)),"median":float(np.median(q[:,j])),"q95":float(np.quantile(q[:,j],.95)),"max":float(np.max(q[:,j]))} for j,a in enumerate("XYZ")}
ranges={k:stats(q) for k,q in zip(keys,G)}
templates=np.array([m["template"] for m in motifs])
bytemplate={}
for t in np.unique(templates):
    mask=templates==t
    groups=[q[mask] for q in G[:3]]+[M[mask][:,[26,32]].reshape(-1,3)]
    bytemplate[t]={"motifs":int(mask.sum()),"ranges_A":{k:stats(q) for k,q in zip(keys,groups)}}
E=np.array([m["source_ester_local_A"] for m in motifs])
ref=np.median(E,axis=0)
maxcount=z["coverage"].sum(1)
report={"status":"TECHNICAL_COMPLETE_DESCRIPTIVE_ONLY","candidate_count":473,"unique_source_sites":len(set(m["source_site"] for m in motifs)),
 "template_counts":{t:int((templates==t).sum()) for t in np.unique(templates)},
 "point_counts":dict(zip(keys,[len(q) for q in G])),
 "coordinate_frame":"Carbonyl C at origin; +X C->O; +Y towards ester oxygen perpendicular to X; +Z cross(X,Y). Units angstrom. Not lab-frame top/bottom.",
 "weights":"Equal weight per candidate; pooled donors contribute two points per candidate. No coverage weighting or candidate subsampling.",
 "intervals":"Min-max and marginal coordinate q05-q95. Not confidence intervals, not a joint 90% volume, not a continuous feasible region.",
 "reference_ester":"Median of saved source-ester coordinates, solely for visual orientation.",
 "Asp_marker":"OD1 retained consistently with previous figures; not asserted to be the closest oxygen to His.",
 "selection_bias":"One FIRST-SUCCESS exemplar per successful source site/template. This finite library is not all accepted geometries or equilibrium sampling.",
 "ranges_A":ranges,"by_template":bytemplate}
plt.rcParams.update({"font.size":11,"figure.facecolor":"white","axes.spines.top":False,"axes.spines.right":False})
def reference(ax,d):
    ax.plot(ref[[1,0,2],d[0]],ref[[1,0,2],d[1]],lw=2,color="#444444",zorder=10)
    ax.scatter(ref[:,d[0]],ref[:,d[1]],c=["#222222","#ce303f","#ce303f"],s=[38,25,22],zorder=11)
    ax.axhline(0,color="#aaaaaa",lw=.6);ax.axvline(0,color="#aaaaaa",lw=.6)
def axes_style(ax,d,limit):
    ax.set_xlim(-limit,limit);ax.set_ylim(-limit,limit);ax.set_aspect("equal")
    ax.set_xlabel("XYZ"[d[0]]+" (angstrom)");ax.set_ylabel("XYZ"[d[1]]+" (angstrom)")
    ax.grid(alpha=.12,lw=.5)
allq=np.concatenate(G);limit=np.ceil(np.abs(allq).max())+1
dims=[(0,1),(0,2),(1,2)]
fig,axs=plt.subplots(1,3,figsize=(14.5,5.4))
for ax,d in zip(axs,dims):
    for j in [2,1,0,3]:
        q=G[j];ax.scatter(q[:,d[0]],q[:,d[1]],s=10,c=colors[j],alpha=.32,edgecolors="none",rasterized=True)
    reference(ax,d);axes_style(ax,d,limit)
    ax.set_title("XYZ"[d[0]]+"-"+"XYZ"[d[1]]+" projection")
handles=[Line2D([],[],lw=0,marker="o",color=c,markersize=7,label=k+"  n="+str(len(q))) for k,q,c in zip(keys,G,colors)]
fig.suptitle("All 473 joint candidates: where the functional atoms lie\nEvery saved candidate included; two equivalent NH donors pooled (946 N positions)",fontsize=14)
fig.legend(handles=handles,loc="lower center",ncol=2,frameon=False,fontsize=10)
fig.subplots_adjust(top=.78,bottom=.22,left=.06,right=.985,wspace=.25)
fig.savefig(O/"all473_positions.png",dpi=180);plt.close(fig)
# Separate atom panels avoid overplotting. Color is cross-site coverage, not density.
fig,axs=plt.subplots(4,3,figsize=(13.5,15),layout="constrained")
norm=matplotlib.colors.Normalize(vmin=0,vmax=int(maxcount.max()))
for j,(name,q) in enumerate(zip(keys,G)):
    cc=np.repeat(maxcount,2) if j==3 else maxcount
    lim=np.ceil(np.abs(q).max())+1
    order=np.argsort(cc)
    for ax,d in zip(axs[j],dims):
        im=ax.scatter(q[order,d[0]],q[order,d[1]],c=cc[order],cmap="viridis",norm=norm,s=13,alpha=.8,edgecolors="none")
        reference(ax,d);axes_style(ax,d,lim)
        ax.set_title(name+" | "+"XYZ"[d[0]]+"-"+"XYZ"[d[1]]+" | n="+str(len(q)),fontsize=11)
fig.suptitle("All-candidate positions, separated by functional atom\nColor = sites fitted by the parent joint motif (out of 578); not point density\nAxis limits differ between atom rows; units are angstrom",fontsize=13)
fig.colorbar(im,ax=axs.ravel().tolist(),label="Passing sites / 578",shrink=.6,pad=.015)
fig.savefig(O/"all473_atom_panels.png",dpi=170);plt.close(fig)
# Marginal intervals. Thin line = complete range; thick = q05-q95; dot=median.
fig,axs=plt.subplots(1,3,figsize=(14.5,4.8))
for ax,dim in zip(axs,"XYZ"):
    for j,(name,c) in enumerate(zip(keys,colors)):
        v=ranges[name][dim]
        ax.plot([v["min"],v["max"]],[j,j],c=c,lw=1.2)
        ax.plot([v["q05"],v["q95"]],[j,j],c=c,lw=7,solid_capstyle="butt")
        ax.scatter(v["median"],j,s=45,c=c,edgecolor="white",lw=.7,zorder=4)
        ax.text((v["q05"]+v["q95"])/2,j-.19,f'{v["q05"]:.2f} to {v["q95"]:.2f}',ha="center",va="bottom",fontsize=10)
    ax.set_ylim(3.6,-.7);ax.set_xlim(-limit,limit);ax.axvline(0,color="#aaaaaa",lw=.7)
    ax.set_yticks(range(4),keys if ax is axs[0] else [""]*4);ax.set_xlabel(dim+" (angstrom)");ax.grid(axis="x",alpha=.2)
fig.suptitle("Coordinate intervals across ALL candidates (not design tolerances)\nThin: min-max | Thick: 5th-95th percentile | Dot: median",fontsize=14)
fig.subplots_adjust(top=.74,bottom=.15,left=.16,right=.985,wspace=.16)
fig.savefig(O/"all473_coordinate_intervals.png",dpi=180);plt.close(fig)
# Pair-correlation descriptor, invariant to swapping donor labels.
N=M[:,[26,32]];v=N-E[:,1,None,:]
a=np.linalg.norm(v,axis=2)
nn=np.degrees(np.arccos(np.clip(np.sum(v[:,0]*v[:,1],axis=1)/(a[:,0]*a[:,1]),-1,1)))
report["donor_pair_N_O_N_deg"]={k:float(x) for k,x in zip(["min","q05","median","q95","max"],np.quantile(nn,[0,.05,.5,.95,1]))}
report["donor_pair_N_N_A"]={k:float(x) for k,x in zip(["min","q05","median","q95","max"],np.quantile(np.linalg.norm(N[:,0]-N[:,1],axis=1),[0,.05,.5,.95,1]))}
(O/"statistics.json").write_text(json.dumps(report,indent=2))
(O/"RUNBOOK.md").write_text("# All-473 distribution plot\n\nSaved joint_transfer_results.npz + motifs.json -> all points, atom-separated projections, marginal coordinate intervals and descriptive JSON. No geometry resampling or changed thresholds. CPU runtime: simulation_slabs/unbiased_100chain_v1/env/bin/python. Run saved script once; new output directory must not exist.\n\n473 exemplars are equal-weight candidates, not 473 independent PET sites or all accepted geometries. 235 IsPETase + 238 LCC. Both identical donor Ns are pooled. Per-coordinate 5-95% bounds are not a joint 90% volume. The source library was selected by first success, so point density is not binding probability. Full coordinates and pair identity remain in the parent NPZ.\n")
for f,h in before.items():assert hashlib.sha256((R/f).read_bytes()).hexdigest()==h
for file in O.glob("*.png"):
    with Image.open(file) as im:im.verify()
assert np.array_equal(np.sort(M[:,[26,32]].reshape(-1,3),axis=0),np.sort(M[:,[32,26]].reshape(-1,3),axis=0))
for group,q in zip(keys,G):
    for j,dim in enumerate("XYZ"):
        v=ranges[group][dim]
        assert v["min"]<=v["q05"]<=v["median"]<=v["q95"]<=v["max"]
        assert abs(v["min"]-q[:,j].min())<1e-12
manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.iterdir() if p.is_file()}
(O/"SHA256.json").write_text(json.dumps(manifest,indent=2))
with (R/"RUN_LOG.jsonl").open("a") as f:
    f.write(json.dumps({"time":datetime.datetime.now(datetime.timezone.utc).isoformat(),"event":"ALL473_DESCRIPTIVE_DISTRIBUTIONS","script":str(Path(__file__)),"output":str(O),"input_sha256":before,"source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"status":"COMPLETE","parent_unchanged":True})+"\n")
print("ALL473_DISTRIBUTIONS_VERIFIED",json.dumps(report),flush=True)
