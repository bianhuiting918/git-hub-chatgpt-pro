"""Readable PET joint motif figures; replot saved statistics only."""
from pathlib import Path
import json,hashlib,datetime
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.spatial.distance import cdist
P=Path("/data/bht2/polymer_material_reference_20260827/simulation_slabs/pet_dp10_400chain_direct_v1")
R=P/"donor_geometry_batch_singleNH_v1/joint_transfer_map_v1"
O=R/"figures_v2";O.mkdir(exist_ok=False)
summary=json.loads((R/"summary.json").read_text());motifs=json.loads((R/"motifs.json").read_text())
z=np.load(R/"joint_transfer_results.npz");M=z["motif_xyz_A"];H=z["motif_H_A"]
colors=["#239566","#4268b7","#a951a1","#dc9426","#069da7"]
names=["Ser OG","His NE2","Asp O","NH-a","NH-b"]
plt.rcParams.update({"font.size":11,"axes.spines.top":False,"axes.spines.right":False,"figure.facecolor":"white"})
def functional(i):
    idx=motifs[i]["functional_indices"]
    return [idx["og"],idx["his_ne2"],idx["asp_o"][0],26,32]
def draw(ax,i,dims,cloud=None):
    m=motifs[i];q=M[i];h=H[i];ester=np.array(m["source_ester_local_A"])
    dims=list(dims);ids=functional(i)
    res=np.array(m["residues"])
    if cloud is not None:
        for j in cloud:
            pts=M[j][functional(j)]
            for n in range(5):
                ax.scatter(*pts[n,dims],s=32,facecolor=colors[n],edgecolor="none",alpha=.24,zorder=2)
    for rn,(r,c) in enumerate(zip(["SER","HIS","ASP","NH1","NH2"],colors)):
        sub=q[res==r];d=cdist(sub,sub)
        for a,b in np.argwhere(np.triu((d>.8)&(d<1.85),1)):
            ax.plot(sub[[a,b],dims[0]],sub[[a,b],dims[1]],color=c,alpha=.35,lw=1.5)
    ax.plot(ester[[1,0,2],dims[0]],ester[[1,0,2],dims[1]],color="#666666",lw=2.5,zorder=5)
    for j,c in enumerate(["#272727","#d63742","#d63742"]):
        ax.scatter(*ester[j,dims],s=90 if j==0 else 65,c=c,zorder=6)
    # Direct ester labels with distinct leader offsets, avoid nearby catalytic markers.
    offsets=[(-20,-23),(4,17),(-47,12)] if dims[1]==2 else [(-18,-23),(1,16),(-48,12)]
    for j,(text,off) in enumerate(zip(["C","O","O-ester"],offsets)):
        ax.annotate(text,ester[j,dims],xytext=off,textcoords="offset points",fontsize=10,
                    arrowprops=dict(arrowstyle="-",lw=.6,color="#666666"))
    markers=["o","s","D","^","v"]
    for n,a in enumerate(ids):
        ax.scatter(*q[a,dims],s=110,c=colors[n],marker=markers[n],edgecolor="white",lw=.7,zorder=7)
    for k,a in enumerate([26,32]):
        seg=np.vstack([q[a],h[k],ester[1]])
        ax.plot(seg[:,dims[0]],seg[:,dims[1]],ls="--",lw=1.3,color=colors[3+k])
        ax.annotate("",h[k,dims],q[a,dims],arrowprops=dict(arrowstyle="->",color=colors[3+k],lw=1))
    attack=np.vstack([ester[0],q[ids[0]]])
    ax.plot(attack[:,dims[0]],attack[:,dims[1]],ls=":",lw=1.5,color=colors[0])
    ax.axhline(0,color="#999999",lw=.5);ax.axvline(0,color="#999999",lw=.5)
    data=M[cloud].reshape(-1,3) if cloud is not None else q
    lo=np.minimum(data.min(0),ester.min(0))-1.5;hi=np.maximum(data.max(0),ester.max(0))+1.5
    ax.set_xlim(lo[dims[0]],hi[dims[0]]);ax.set_ylim(lo[dims[1]],hi[dims[1]])
    ax.set_aspect("equal",adjustable="box")
    ax.set_xlabel("X: carbonyl C -> O (angstrom)")
    ax.set_ylabel(("Y: towards ester O" if dims[1]==1 else "Z: out of ester plane")+" (angstrom)")
handles=[Line2D([],[],color=c,marker=mk,lw=0,markersize=9,label=n) for c,mk,n in zip(colors,["o","s","D","^","v"],names)]
handles += [Line2D([],[],color="#666666",ls="--",label="NH -> carbonyl O"),
            Line2D([],[],color=colors[0],ls=":",label="Ser OG -> carbonyl C")]
best=summary["top10"]["support70_all"][0];i=best["index"]
fig,axs=plt.subplots(1,2,figsize=(12.5,6))
for ax,dims,title in zip(axs,[(0,1),(0,2)],["View in ester plane (X-Y)","Side view (X-Z)"]):
    draw(ax,i,dims);ax.set_title(title)
fig.suptitle(f"Overall widest observed coverage: {best['id']} | {best['template']} | {best['count']}/{best['total']} sites (5.4%)\nOne intact, rigid triad + two-NH arrangement; no pose optimization",fontsize=14)
fig.legend(handles=handles,loc="lower center",ncol=4,frameon=False,fontsize=10)
fig.subplots_adjust(top=.79,bottom=.19,wspace=.27,left=.08,right=.97)
fig.savefig(O/"overall_joint_geometry.png",dpi=180);plt.close(fig)
# All cohorts, compact functional-only labels through a single legend.
fig,axs=plt.subplots(3,2,figsize=(11,14))
selected=[]
for row,key in enumerate(["support90","support80","support70_all"]):
    v=summary["top10"][key][0];j=v["index"];selected.append(j)
    for ax,dims in zip(axs[row],[(0,1),(0,2)]):
        draw(ax,j,dims)
        ax.set_title(f"Support >= {['90%','80%','70% / ALL'][row]} | {v['id']} {v['template']}\n{v['count']}/{v['total']} sites ({v['fraction']*100:.1f}%)",fontsize=12)
# identical scales across cohorts within each projection
lo=M[selected].min((0,1))-1.5;hi=M[selected].max((0,1))+1.5
for row in range(3):
    for col in range(2):
        axs[row,col].set_xlim(lo[0],hi[0]);axs[row,col].set_ylim(lo[col+1],hi[col+1])
fig.suptitle("Best-observed joint placements by PET support condition\nLeft: X-Y | Right: X-Z | same coordinates and scales",fontsize=14)
fig.legend(handles=handles,loc="lower center",ncol=4,frameon=False,fontsize=10)
fig.subplots_adjust(top=.91,bottom=.085,hspace=.47,wspace=.35,left=.1,right=.97)
fig.savefig(O/"cohort_joint_geometry.png",dpi=170);plt.close(fig)
cloud=[v["index"] for v in summary["top10"]["support70_all"]]
fig,axs=plt.subplots(1,2,figsize=(12.5,7))
for ax,dims,title in zip(axs,[(0,1),(0,2)],["X-Y: in-plane positions","X-Z: out-of-plane positions"]):
    draw(ax,i,dims,cloud);ax.set_title(title)
fig.suptitle("Spatial spread of the 10 highest-coverage joint motifs\nSmall translucent points: top 10 | large shaped markers and sticks: M450\nPositions are correlated: do not mix donors from different motifs",fontsize=13)
fig.legend(handles=handles,loc="lower center",ncol=4,frameon=False,fontsize=10)
fig.subplots_adjust(top=.77,bottom=.16,wspace=.3,left=.09,right=.97)
fig.savefig(O/"top10_joint_spread.png",dpi=180);plt.close(fig)
# Coordinates and actual geometric descriptors for the three selected exemplars.
desc=[]
for j in selected:
    q=M[j];h=H[j];ester=np.array(motifs[j]["source_ester_local_A"]);ii=functional(j)
    v=q[ii[0]];theta=float(np.degrees(np.arccos(v[0]/np.linalg.norm(v))));phi=float(np.degrees(np.arctan2(v[2],v[1])))
    u=q[[26,32]]-ester[1];nn=float(np.degrees(np.arccos(np.clip(np.dot(u[0],u[1])/np.linalg.norm(u[0])/np.linalg.norm(u[1]),-1,1))))
    desc.append({"id":motifs[j]["id"],"Ser_OG_C_A":float(np.linalg.norm(v)),"O_C_Ser_OG_deg":theta,"Ser_attack_azimuth_deg":phi,"N_O_A":np.linalg.norm(u,axis=1).tolist(),"N_O_N_deg":nn,"functional_positions_A":{n:q[k].tolist() for n,k in zip(names,ii)}})
(O/"descriptors.json").write_text(json.dumps(desc,indent=2))
manifest={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in O.iterdir() if p.is_file()}
(O/"SHA256.json").write_text(json.dumps(manifest,indent=2))
with (R/"RUN_LOG.jsonl").open("a") as f:
    f.write(json.dumps({"time":datetime.datetime.now(datetime.timezone.utc).isoformat(),"event":"REPLOT_ONLY","script":str(Path(__file__)),"new_output":str(O),"geometry_unchanged":True,"status":"COMPLETE"})+"\n")
print("REPLOT_COMPLETE",json.dumps(desc))
