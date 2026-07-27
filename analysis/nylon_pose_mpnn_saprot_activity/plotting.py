import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa:F401
import plotly.express as px
from .schema import ENDPOINT_LABELS,PRIMARY_PREDICTORS

def _md5_set_sha(md5s):
    return hashlib.sha256("\n".join(sorted(map(str,md5s))).encode()).hexdigest()

def _regression_band(x,y,grid,seed=20260727,n_boot=500):
    if len(x)<3 or np.std(x)==0:
        return np.full_like(grid,np.nan),np.full_like(grid,np.nan),np.full_like(grid,np.nan)
    coef=np.polyfit(x,y,1); mid=np.polyval(coef,grid)
    rng=np.random.default_rng(seed); preds=[]
    for _ in range(n_boot):
        idx=rng.integers(0,len(x),len(x))
        if np.std(x[idx])==0: continue
        preds.append(np.polyval(np.polyfit(x[idx],y[idx],1),grid))
    if not preds: return mid,np.full_like(grid,np.nan),np.full_like(grid,np.nan)
    lo,hi=np.percentile(np.asarray(preds),[2.5,97.5],axis=0)
    return mid,lo,hi

def render_endpoint_figures(frame,endpoint,out_dir):
    out=Path(out_dir); out.mkdir(parents=True,exist_ok=True)
    cols=["enzyme","sequence_md5"]+PRIMARY_PREDICTORS+[endpoint]
    d=frame[cols].dropna().copy()
    label=ENDPOINT_LABELS[endpoint]
    hover={"sequence_md5":True,"enzyme":True}
    fig=px.scatter_3d(d,x="pose_strict_nac",y="mpnn_vanilla_mean",
        z="saprot_full_mean",color=endpoint,hover_name="enzyme",hover_data=hover,
        labels={"pose_strict_nac":"Strict NAC pose-retention fraction",
                "mpnn_vanilla_mean":"ProteinMPNN vanilla v_48_020 mean log-likelihood",
                "saprot_full_mean":"SaProt full-protein mean log-likelihood",
                endpoint:label},color_continuous_scale="Viridis")
    fig.update_traces(marker={"size":5,"opacity":0.85})
    html=out/f"{endpoint}_3d.html"
    fig.write_html(html,include_plotlyjs=True,full_html=True)
    f=plt.figure(figsize=(8,6)); ax=f.add_subplot(111,projection="3d")
    sc=ax.scatter(d.pose_strict_nac,d.mpnn_vanilla_mean,d.saprot_full_mean,
                  c=d[endpoint],cmap="viridis",s=38,edgecolor="k",linewidth=.25)
    ax.set_xlabel("Strict NAC retention"); ax.set_ylabel("ProteinMPNN mean LL")
    ax.set_zlabel("SaProt mean LL"); f.colorbar(sc,ax=ax,pad=.12,label=label)
    f.tight_layout(); f.savefig(out/f"{endpoint}_3d.png",dpi=220); plt.close(f)
    f,axes=plt.subplots(1,3,figsize=(15,4.4))
    xlabels=["Strict NAC pose-retention fraction",
             "ProteinMPNN vanilla v_48_020 mean log-likelihood",
             "SaProt full-protein mean log-likelihood"]
    for i,(ax,pred,xlabel) in enumerate(zip(axes,PRIMARY_PREDICTORS,xlabels)):
        x=d[pred].to_numpy(float); y=d[endpoint].to_numpy(float)
        ax.scatter(x,y,c=d[endpoint],cmap="viridis",edgecolor="k",linewidth=.25,s=38)
        if len(x) and np.ptp(x)>0:
            grid=np.linspace(x.min(),x.max(),100)
            mid,lo,hi=_regression_band(x,y,grid,seed=20260727+i)
            ax.plot(grid,mid,color="#b2182b",lw=1.5)
            ax.fill_between(grid,lo,hi,color="#b2182b",alpha=.18)
        ax.set_xlabel(xlabel); ax.set_ylabel(label)
    f.tight_layout(); f.savefig(out/f"{endpoint}_2d_panels.png",dpi=220); plt.close(f)
    audit={"endpoint":endpoint,"n_points":int(len(d)),
           "sequence_md5_sha256":_md5_set_sha(d.sequence_md5),
           "y_axis_label":label,
           "files":[html.name,f"{endpoint}_3d.png",f"{endpoint}_2d_panels.png"]}
    (out/f"{endpoint}_figure_audit.json").write_text(json.dumps(audit,indent=2)+"\n")
    return audit
