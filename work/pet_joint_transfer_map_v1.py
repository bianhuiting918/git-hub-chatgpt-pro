"""Registered PET triad + two-NH motif transfer screen; CPU, geometry only.
This ranks an existing FIRST-SUCCESS candidate library, not global optima.
No per-site structures; no changes to parent protocol/results.
"""
from pathlib import Path
import json, hashlib, datetime, time, sys, csv
import numpy as np
from scipy.spatial.distance import cdist
P=Path("/data/bht2/polymer_material_reference_20260827/simulation_slabs/pet_dp10_400chain_direct_v1")
SRC=P/"pet_donor_geometry_batch_singleNH_v1.py"
EXPECTED="b796be0506a0d696ee19aad9288e2b71b6c81131c30f7e4bf2deb00c94ddb80a"
assert hashlib.sha256(SRC.read_bytes()).hexdigest()==EXPECTED
g={"__name__":"parent_functions"}
s=SRC.read_text()
exec(compile(s[:s.index("rows=collect()")],str(SRC),"exec"),g)
root=P/"donor_geometry_batch_singleNH_v1"/"joint_transfer_map_v1"
root.mkdir(exist_ok=False)
t0=time.monotonic()
def log(event,**kw):
    with (root/"RUN_LOG.jsonl").open("a") as f:
        f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),event=event,**kw))+"\n")
log("START",source=str(Path(__file__)),parent_sha256=EXPECTED)
base=P/"donor_geometry_batch_singleNH_v1"/"extension70_v1"
rows=json.loads((base/"eligible_sites.json").read_text())
assert len(rows)==578 and len(set(r["site_id"] for r in rows))==578
results=[json.loads((base/"cases"/(r["site_id"]+".json")).read_text()) for r in rows]
parent_hashes={str(base/"cases"/(r["site_id"]+".json")):hashlib.sha256((base/"cases"/(r["site_id"]+".json")).read_bytes()).hexdigest() for r in rows}
np.seterr(invalid="raise")
def angle(a,b):
    return np.degrees(np.arccos(np.clip(np.sum(a*b,axis=-1)/np.linalg.norm(a,axis=-1)/np.linalg.norm(b,axis=-1),-1,1)))
def local_env(r):
    ester,frame,ids=g["site_frame"](r)
    q=g["xyz"].copy(); q[:,:2]-=np.round((q[:,:2]-ester[0,:2])/g["BOX"][:2])*g["BOX"][:2]
    mask=np.linalg.norm(q-ester[0],axis=1)<=30
    return ester,frame,g["trees"](q[mask],g["elems"][mask])
motifs=[]; coords=[]; hydrogens=[]; radii=[]; indices=[]; self_errors=[]
for k,(r,result) in enumerate(zip(rows,results)):
    hits=[t for t in result["templates"] if t["example_audit"] is not None]
    if not hits:continue
    ester,frame,env=local_env(r)
    q,H,*_=g["sample_donor"](ester,frame,g["seed_for"](r["site_id"]+"|singleNH|Sobol"))
    for tr in hits:
        name=tr["template"]; ex=tr["example_audit"]; template=g["templates"][name]
        g["rng"]=np.random.default_rng(g["seed_for"](r["site_id"]+"|"+name+"|triad"))
        poses,stats=g["make_triads"](ester,frame,template,env)
        assert stats==tr["triad"]
        pose=poses[ex["pose"]]; pair=[ex["sample_i"],ex["sample_j"]]
        m=np.concatenate([pose,q[pair[0]],q[pair[1]]])
        h=H[pair]
        els=list(template[1])+list(g["donor_templates"][0][2])*2
        rr=np.array([g["RAD"][e] for e in els])
        assert g["clear"](m[None],rr,env)[0]>=-1e-8
        dr=rr[24:30]
        assert (cdist(m[24:30],m[30:36])-dr[:,None]-dr[None,:]+.4).min()>=-1e-8
        for frag in [m[24:30],m[30:36]]:
            assert g["clear"](frag[None],dr,g["trees"](pose,template[1]))[0]>=-1e-8
        lc=(m-ester[0])@frame; lh=(h-ester[0])@frame
        atomlines=template[4]
        names=[l[12:16].strip() for l in atomlines]; res=[l[17:20] for l in atomlines]
        idx={"og":template[2],"cb":template[3],
             "his_ne2":next(i for i in range(24) if names[i]=="NE2" and res[i]=="HIS"),
             "asp_o": [i for i in range(24) if res[i]=="ASP" and names[i] in ["OD1","OD2"]]}
        motifs.append({"id":"M%03d"%(len(motifs)+1),"template":name,"source_site":r["site_id"],"source_support":r["support_fraction"],"source_pose":ex["pose"],"source_donor_samples":pair,"elements":els,"atom_names":names+["Cprev","Oprev","N","CA","C","O"]*2,"residues":res+["NH1"]*6+["NH2"]*6,"functional_indices":idx,"source_ester_local_A":((ester-ester[0])@frame).tolist()})
        coords.append(lc);hydrogens.append(lh);radii.append(rr);indices.append(idx)
    if k%40==0:print("RECONSTRUCT",k,"motifs",len(motifs),flush=True)
M=np.array(coords); HH=np.array(hydrogens); RR=np.array(radii)
assert len(motifs)==473,(len(motifs),"expected 235 IsPETase + 238 LCC")
og=np.array([x["og"] for x in indices]);cb=np.array([x["cb"] for x in indices]);ii=np.arange(len(M))
N=M[:,[26,32]]; OG=M[ii,og];CB=M[ii,cb]
# The full motif is rigid. Carbonyl C at origin; carbonyl O lies on +X.
assert np.max(np.linalg.norm(M,axis=-1))+3.6<30
counts={}; okmat=np.zeros((len(M),len(rows)),bool); geommat=np.zeros_like(okmat)
margins=np.full(okmat.shape,np.nan)
chemmax=[]; chemmin=[]
for k,r in enumerate(rows):
    ester,frame,env=local_env(r); targetO=(ester[1]-ester[0])@frame
    da=np.linalg.norm(OG,axis=1); aa=angle(OG,np.broadcast_to(targetO,OG.shape))
    sa=angle(CB-OG,-OG)
    no=np.linalg.norm(N-targetO,axis=-1); ha=angle(N-HH,targetO-HH)
    valid=(da>=3-1e-8)&(da<=3.4+1e-8)&(aa>=95-1e-8)&(aa<=115+1e-8)&(sa>=100-1e-8)&(sa<=130+1e-8)&np.all((no>=2.7-1e-8)&(no<=3.2+1e-8)&(ha>=140-1e-8),axis=1)
    geommat[:,k]=valid
    which=np.flatnonzero(valid)
    w=M[which]@frame.T+ester[0]
    v=g["clear"](w,RR[which],env)
    margins[which,k]=v;okmat[which,k]=v>=-1e-8
    if k%80==0:print("TRANSFER",k,"/",len(rows),flush=True)
for i,m in enumerate(motifs):
    k=next(j for j,r in enumerate(rows) if r["site_id"]==m["source_site"])
    assert okmat[i,k],("source replay failed",m["id"],r["site_id"])
# Swap NH fragments must not change any chemical or collision decision.
swap=np.r_[np.arange(24),np.arange(30,36),np.arange(24,30)]
assert np.allclose(np.sort(RR,axis=1),np.sort(RR[:,swap],axis=1))
support=np.array([r["support_fraction"] for r in rows])
cohorts={"support90":support>=.9,"support80":support>=.8,"support70_all":support>=.7,
         "only80to90":(support>=.8)&(support<.9),"only70to80":(support>=.7)&(support<.8)}
den={k:int(v.sum()) for k,v in cohorts.items()}
assert list(den.values())==[344,470,578,126,108]
rankings={}
for key,mask in cohorts.items():
    c=okmat[:,mask].sum(axis=1)
    order=sorted(range(len(M)),key=lambda i:(-int(c[i]),-float(np.nanmedian(margins[i,mask][okmat[i,mask]])) if c[i] else 0,motifs[i]["id"]))
    rankings[key]=[dict(id=motifs[i]["id"],template=motifs[i]["template"],count=int(c[i]),total=int(mask.sum()),fraction=float(c[i]/mask.sum()),index=int(i)) for i in order]
    counts[key]=c
    with (root/(key+"_ranking.csv")).open("w") as f:
        w=csv.DictWriter(f,fieldnames=list(rankings[key][0]));w.writeheader();w.writerows(rankings[key])
# Independent full-polymer cdist replay on 6 positives and 6 negatives of overall best.
best=rankings["support70_all"][0]["index"]
audit=[]
for target in list(np.flatnonzero(okmat[best])[:6])+list(np.flatnonzero(~okmat[best])[:6]):
    r=rows[target];ester,frame,_=g["site_frame"](r);q=g["xyz"].copy()
    q[:,:2]-=np.round((q[:,:2]-ester[0,:2])/g["BOX"][:2])*g["BOX"][:2]
    fullrad=np.array([g["RAD"][e] for e in g["elems"]])
    candidate=M[best]@frame.T+ester[0]
    margin=float((cdist(candidate,q)-RR[best,:,None]-fullrad[None,:]+.4).min())
    exact=bool(geommat[best,target] and margin>=-1e-8)
    assert exact==bool(okmat[best,target])
    audit.append({"site":r["site_id"],"passed":exact,"full_PET_min_margin_A":margin})
for p,h in parent_hashes.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
np.savez_compressed(root/"joint_transfer_results.npz",motif_xyz_A=M,motif_H_A=HH,coverage=okmat,chemical_pass=geommat,min_steric_margin_A=margins,site_support=support,site_ids=np.array([r["site_id"] for r in rows]))
(root/"motifs.json").write_text(json.dumps(motifs))
summary={"status":"TECHNICAL_COMPLETE_GEOMETRY_ONLY","candidate_count":len(motifs),"unique_target_sites":len(rows),"denominators":den,"top10":{k:v[:10] for k,v in rankings.items()},"candidate_library":"one previously audited first-success pair per site/template; not exhaustive; adaptive source sample budget preserved","registration":"C origin; X C->O; Y towards ester-O projected perpendicular X; Z right-hand normal; no rigid-body pose optimization in transfer","interpretation":"coverage in this finite candidate library at fixed ester registration; not energy, activity, binding affinity or global optimum","source_replay_pass":len(motifs),"independent_full_PET_replays":audit,"parent_cases_hashes_unchanged":len(parent_hashes),"elapsed_s":time.monotonic()-t0}
(root/"summary.json").write_text(json.dumps(summary,indent=2))
(root/"config.json").write_text(json.dumps({"parent":g["BATCH_CFG"],"parent_sha256":EXPECTED,"candidate_count":len(motifs),"registration":summary["registration"],"source_script_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest()},indent=2))
(root/"RUNBOOK.md").write_text("# PET joint motif transfer map v1\n\nCPU only. Run the saved script with OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 and the simulation_slabs/unbiased_100chain_v1/env/bin/python runtime. Output directory must not already exist. Parent data are read-only.\n\nOne previously found joint motif per successful site/template is replayed exactly, then transferred rigidly in the local ester coordinate system to all 578 sites. No pose optimization or global optimum claim. Source library has first-hit and finite/adaptive-sampling selection bias. Donors are exchangeable; only one ordering evaluated. The 70% set is the complete deduplicated set, not an independent fourth cohort.\n\nAll chemical and heavy-atom steric criteria match the parent protocol. This is a reactant-state geometric screen; no TS, energies, full-backbone closure, dynamics or entry trajectory evaluated. Thresholds were not relaxed.\n")
# Publication-style static figures: numbers and actual coordinates only.
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
plt.rcParams.update({"font.size":10,"axes.spines.top":False,"axes.spines.right":False,"figure.facecolor":"white"})
topids=[]
for key in ["support90","support80","support70_all"]:
    for v in rankings[key][:4]:
        if v["index"] not in topids:topids.append(v["index"])
mat=np.array([[counts[key][i]/den[key] for key in ["support90","support80","support70_all"]] for i in topids])
fig,ax=plt.subplots(figsize=(8.6,max(4.2,len(topids)*.47+1.6)),layout="constrained")
im=ax.imshow(mat,cmap="YlGnBu",vmin=0,vmax=max(.01,float(mat.max())),aspect="auto")
for a,i in enumerate(topids):
    for b,key in enumerate(["support90","support80","support70_all"]):
        ax.text(b,a,f"{counts[key][i]}/{den[key]}\n{mat[a,b]*100:.1f}%",ha="center",va="center",fontsize=10,color="white" if mat[a,b]>mat.max()*.55 else "#172c3e")
ax.set_yticks(range(len(topids)),[motifs[i]["id"]+" | "+motifs[i]["template"] for i in topids])
ax.set_xticks(range(3),["Support >=90%","Support >=80%","Support >=70%\nALL unique sites"])
ax.set_title("Which fixed triad + two-NH arrangement fits more PET sites?\nTop-four union from each cohort; finite 473-motif library",pad=14)
fig.colorbar(im,ax=ax,label="Fraction of sites passing geometry")
fig.savefig(root/"coverage_ranking.png",dpi=180);plt.close(fig)
# Coordinate maps: three best motifs, two orthogonal views; identical axes.
wins=[rankings[k][0]["index"] for k in ["support90","support80","support70_all"]]
colors={"SER":"#24a06b","HIS":"#4967bc","ASP":"#b05ea8","NH1":"#e89826","NH2":"#cf5446"}
fig,axs=plt.subplots(3,2,figsize=(11,13),layout="constrained")
lo=np.min(M[wins],axis=(0,1))-1.0;hi=np.max(M[wins],axis=(0,1))+1.0
for rownum,(i,key) in enumerate(zip(wins,["support90","support80","support70_all"])):
    m=motifs[i];xyz=M[i];h=HH[i];ester=np.array(m["source_ester_local_A"]);idx=indices[i]
    res=np.array(m["residues"])
    groups=[np.flatnonzero(res==r) for r in ["SER","HIS","ASP","NH1","NH2"]]
    for col,dims in enumerate([(0,1),(0,2)]):
        ax=axs[rownum,col]
        for r,ids in zip(["SER","HIS","ASP","NH1","NH2"],groups):
            q=xyz[ids];dd=cdist(q,q)
            for a,b in np.argwhere(np.triu((dd>.8)&(dd<1.85),1)):
                ax.plot(q[[a,b],dims[0]],q[[a,b],dims[1]],color=colors[r],alpha=.48,lw=1.8)
            ax.scatter(q[:,dims[0]],q[:,dims[1]],s=20,color=colors[r],alpha=.5,zorder=3)
        ax.plot(ester[[1,0,2],dims[0]],ester[[1,0,2],dims[1]],color="#555555",lw=3,zorder=5)
        ax.scatter(ester[:,dims[0]],ester[:,dims[1]],s=[90,75,55],c=["#343434","#d52a38","#d52a38"],zorder=6)
        for a,label in enumerate(["C","O","O-ester"]):
            ax.annotate(label,(ester[a,dims[0]],ester[a,dims[1]]),xytext=(4,-13 if a==0 else 5),textcoords="offset points",fontsize=10)
        for a,label,clr in [(idx["og"],"Ser OG",colors["SER"]),(idx["his_ne2"],"His NE2",colors["HIS"]),(idx["asp_o"][0],"Asp O",colors["ASP"]),(26,"NH-a",colors["NH1"]),(32,"NH-b",colors["NH2"])]:
            ax.scatter(xyz[a,dims[0]],xyz[a,dims[1]],s=70,color=clr,edgecolor="white",zorder=7)
            ax.annotate(label,(xyz[a,dims[0]],xyz[a,dims[1]]),xytext=(6,5),textcoords="offset points",fontsize=9,color=clr)
        for d,a in enumerate([26,32]):
            seg=np.stack([xyz[a],h[d],ester[1]])
            ax.plot(seg[:,dims[0]],seg[:,dims[1]],ls="--",lw=1,color=colors["NH1" if d==0 else "NH2"])
        attack=np.stack([ester[0],xyz[idx["og"]]])
        ax.plot(attack[:,dims[0]],attack[:,dims[1]],ls=":",color=colors["SER"],lw=1.6)
        ax.axhline(0,color="#aaaaaa",lw=.5);ax.axvline(0,color="#aaaaaa",lw=.5)
        ax.set_xlim(lo[dims[0]],hi[dims[0]]);ax.set_ylim(lo[dims[1]],hi[dims[1]])
        ax.set_aspect("equal",adjustable="box");ax.set_xlabel("X: carbonyl C -> O (angstrom)")
        ax.set_ylabel(("Y: towards ester O" if col==0 else "Z: out of ester plane")+" (angstrom)")
        ax.set_title(f"{key.replace('support','Support ').replace('_all',' / ALL')} | {m['id']} {m['template']}\n{counts[key][i]}/{den[key]} sites ({100*counts[key][i]/den[key]:.1f}%) | "+("ester-plane view" if col==0 else "side view"),fontsize=11)
fig.suptitle("Actual joint arrangements with the widest observed coverage\nSame ester-based coordinates and axis limits; not independent atom averages",fontsize=14)
fig.savefig(root/"joint_position_maps.png",dpi=180);plt.close(fig)
# A compact overall winner in two views, reuse figure via per-row creation would be redundant.
# Full top-10 coordinate table allows downstream PyMOL/de novo motif use without PDB flood.
bestset=sorted(set(v["index"] for key in ["support90","support80","support70_all"] for v in rankings[key][:10]))
with (root/"top_motif_coordinates.csv").open("w") as f:
    w=csv.writer(f);w.writerow(["motif","template","residue","atom","x_A","y_A","z_A"])
    for i in bestset:
        for j,q in enumerate(M[i]):
            w.writerow([motifs[i]["id"],motifs[i]["template"],motifs[i]["residues"][j],motifs[i]["atom_names"][j],*q])
hashes={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir() if p.is_file() and p.name!="RUN_LOG.jsonl"}
(root/"SHA256.json").write_text(json.dumps(hashes,indent=2))
log("COMPLETE",seconds=time.monotonic()-t0,candidates=len(motifs),target_sites=len(rows),status=summary["status"])
print("JOINT_TRANSFER_COMPLETE",json.dumps(summary),flush=True)
