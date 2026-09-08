import pathlib,json,numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
R=pathlib.Path(__file__).resolve().parent
OUT=R/"multiscale_material_match_pilot_report_v1";OUT.mkdir(exist_ok=True)
roots={mat:R/("multiscale_material_match_pilot_v1" if mat=="PET" else "multiscale_material_match_"+mat+"_pilot_v1") for mat in ["PET","PA6","PA66"]}
summary=[]
for mat,root in roots.items():
 if not (root/"summary.json").exists():continue
 d=json.load(open(root/"summary.json"));summary.append(d)
 for group in ["all","current_core_accessible"]:
  fig=plt.figure(figsize=(12,10))
  for row,channel in enumerate(["lipo","HBA","HBD"]):
   for col,radius in enumerate([10,15,20]):
    rec=next(x for x in d["results"] if x["cohort"]==group and x["radius_A"]==radius and x["target_fraction"]==.7)
    key=f"{mat}_{group}_coverage70_R{radius}"
    p=root/key/(key+".npz")
    ax=fig.add_subplot(3,3,row*3+col+1,projection="3d")
    if p.exists():
     z=np.load(p);v=z["vertices_A"];f=z["faces"];colors=z["colors_rgb"]
     if channel!="lipo":
      val=z["material_"+channel+"_weight"]
      target=np.array([.15,.55,.85] if channel=="HBA" else [.7,.25,.7])
      colors=.93+(target-.93)*val[:,None];colors[~z["chemical_supported"]]=.55
     poly=Poly3DCollection(v[f],facecolors=colors[f].mean(1),edgecolors="none",alpha=.95,rasterized=True)
     ax.add_collection3d(poly)
    ax.scatter([0],[0],[0],color="cyan",edgecolor="black",s=12,linewidth=.4)
    ax.set(xlim=(-21,21),ylim=(-21,21),zlim=(-21,21),xlabel="x (A)",ylabel="y (A)",zlabel="z (A)")
    ax.set_xticks([-20,0,20]);ax.set_yticks([-20,0,20]);ax.set_zticks([-20,0,20])
    ax.view_init(elev=30,azim=-60);ax.set_box_aspect((1,1,1))
    tol=rec["tolerance_for_same_subset_A"]
    ax.set_title(f"R{radius} A | {channel}\n{rec['actual_n_at_tolerance']}/{rec['denominator']} at D95 <= {tol:.2f} A",fontsize=10)
  fig.suptitle(f"PILOT ONLY: {mat}, {group}, target >=70%\nLarge tolerance is NOT good shape matching",fontsize=13)
  fig.text(.5,.015,"Lipo: brown positive / blue negative fragment index; HBA blue; HBD purple; gray low support.\nChemical annotation is descriptive, not energy. Cyan point = carbonyl C.",ha="center",fontsize=9)
  fig.subplots_adjust(top=.88,bottom=.08,hspace=.25,wspace=.12)
  for ext in ["png","svg"]:fig.savefig(OUT/(mat+"_"+group+"_pilot_surfaces."+ext),dpi=110)
  plt.close(fig)
fig,axes=plt.subplots(2,3,figsize=(12,7),sharex=True,sharey=True)
for col,d in enumerate(summary):
 for row,group in enumerate(["all","current_core_accessible"]):
  ax=axes[row,col]
  for frac,color in [(.9,"#b54835"),(.7,"#2878ac"),(.5,"#588846")]:
   vals=[next(x for x in d["results"] if x["cohort"]==group and x["radius_A"]==r and x["target_fraction"]==frac)["tolerance_for_same_subset_A"] for r in [10,15,20]]
   ax.plot([10,15,20],vals,"o-",color=color,label=f"target {int(frac*100)}%")
  ax.axhspan(0,3,color="gray",alpha=.10);ax.set_title(d["material"]+" | "+group,fontsize=10)
  ax.set_xlabel("Patch radius (A)");ax.set_ylabel("Required D95 tolerance (A)")
  ax.set_xticks([10,15,20]);ax.set_ylim(bottom=0);ax.grid(alpha=.2);ax.legend(fontsize=8)
fig.suptitle("Pilot: actual shape matching requires a stated tolerance\n8 sites/material; accessible subset 4 sites/material; NOT full-population results")
fig.tight_layout(rect=(0,0,1,.91))
for ext in ["png","svg"]:fig.savefig(OUT/("pilot_required_tolerance."+ext),dpi=130)
plt.close(fig)
print("PILOT_FIGURES",len(summary),str(OUT))
