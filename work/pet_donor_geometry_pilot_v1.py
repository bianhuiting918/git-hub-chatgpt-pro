"""PET catalytic geometry pilot. CPU only. No activity/energy inference."""
from pathlib import Path
import json, csv, hashlib, datetime, sys
import numpy as np
from scipy.ndimage import gaussian_filter, label, map_coordinates
from scipy.spatial import cKDTree
from scipy.spatial.transform import Rotation
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

P=Path("/data/bht2/polymer_material_reference_20260827/simulation_slabs/pet_dp10_400chain_direct_v1")
OUT=P/"donor_geometry_pilot_v1"
BOX=np.array([125.697,125.697,182.848])
RAD={"C":1.70,"O":1.52,"N":1.55,"S":1.80}
CFG={"seed":20260907,"n_triad_trials":12000,"n_donor_trials_per_type":16000,
 "max_conditioned_poses":64,"max_donors_per_type_for_pair_test":64,
 "attack_distance_A":[3.0,3.4],"attack_angle_Ocarbon_C_SerOG_deg":[95,115],
 "SerCB_OG_C_deg":[100,130],"donor_N_O_A":[2.7,3.2],"N_H_O_min_deg":140,
 "NH_A":1.01,"vdw_overlap_allowance_A":0.4,
 "scope":"single-frame reactant-state geometry; isolated peptide donors; no connectivity requirement; no electronic energies, TS, ingress paths or dynamic chain classification",
 "sampling":"uniform attack cos(angle), azimuth, distance and SO3; uniform donor solid angle and radius, H cone 40deg and peptide roll; NOT equilibrium probabilities"}
rng=np.random.default_rng(CFG["seed"])

def norm(x): return x/np.linalg.norm(x,axis=-1,keepdims=True)
def load_pdb(path,protein=False):
    records=[]
    for l in path.read_text().splitlines():
        if not l.startswith(("ATOM  ","HETATM")):continue
        if protein and (not l.startswith("ATOM  ") or l[21]!="A" or l[16] not in (" ","A")):continue
        e=l[76:78].strip()
        if e not in RAD:continue
        records.append((l,np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]),e))
    return records

rawpath=P/"pet_dp10_400chain/water_slab_v1/dry_export/PET_DP10_400chain_water_equilibrated_dry.pdb"
raw=load_pdb(rawpath)
xyz=np.array([a[1] for a in raw]); elems=np.array([a[2] for a in raw]); xyz[:,:2]%=BOX[:2]
lookup={(l[72:76].strip(),l[12:16].strip()):i for i,(l,_,_) in enumerate(raw)}
# Same mass-density bulk support selection as preview, now retaining all positive carbon areas.
alllines=[l for l in rawpath.read_text().splitlines() if l.startswith(("ATOM  ","HETATM"))]
ax=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in alllines]);ax[:,:2]%=BOX[:2]
mass=np.array([{"C":12.011,"O":15.999,"H":1.008}[l[76:78].strip()] for l in alllines])
n=np.ceil(BOX).astype(int);step=BOX/n; edges=[np.linspace(0,BOX[j],n[j]+1) for j in range(3)]
rho=gaussian_filter(np.histogramdd(ax,bins=edges,weights=mass)[0]/np.prod(step),3/step,mode=("wrap","wrap","constant"))
zc=(edges[2][1:]+edges[2][:-1])/2
threshold=float(np.median(rho[:,:,abs(zc-np.median(ax[:,2]))<10])*.5)
lab,nlab=label(rho>=threshold);sizes=np.bincount(lab.ravel());sizes[0]=0
dense=lab==sizes.argmax();height=zc[np.max(np.where(dense,np.arange(n[2]),0),axis=2)]
disk=np.array([(x,y,-5) for x in np.arange(-10,10.1,2) for y in np.arange(-10,10.1,2) if x*x+y*y<=100])
def dens(q):return map_coordinates(rho,(q/step-.5).T,order=1,mode="nearest")
eligible=[]
with (P/"surface_chemistry_shape_20260907_v1/results_exterior_v3_uncertainty/esters_top_15A.csv").open() as f:
    for row in csv.DictReader(f):
        c=np.array([float(row["center_"+a+"_A"]) for a in "xyz"])
        if row["shape_scope_status"]!="COMPLETE_FOR_CLASSIFIED_MESH" or float(row["carbonyl_C_surface_area_A2"])<=0 or not all(18<c[j]<BOX[j]-18 for j in (0,1)):continue
        h=float(map_coordinates(height,(c[:2]/step[:2]-.5)[:,None],order=1)[0])
        support=float(np.mean(dens(c+disk)>=threshold))
        if abs(c[2]-h)<=3 and support>=.9:
            row.update(density_boundary_offset_A=float(c[2]-h),inward_support=support)
            eligible.append(row)
eligible.sort(key=lambda r:float(r["carbonyl_C_surface_area_A2"]))
assert len(eligible)>=6
selected=[eligible[round(q*(len(eligible)-1))] for q in [.05,.2,.4,.6,.8,.95]]
assert len({r["site_id"] for r in selected})==6
for rank,row in enumerate(selected):row["exposure_rank_group"]=["low","low","middle","middle","high","high"][rank]

def trees(q,e):return {t:cKDTree(q[e==t]) for t in np.unique(e)}
def clear(q,radii,env):
    shape=q.shape[:-2];flat=q.reshape(-1,3); r=np.broadcast_to(radii,q.shape[:-1]).ravel()
    margin=np.full(len(flat),np.inf)
    for t,tree in env.items():
        d=tree.query(flat,workers=1)[0];margin=np.minimum(margin,d-r-RAD[t]+CFG["vdw_overlap_allowance_A"])
    return margin.reshape(*shape,-1).min(-1)

# Exact analytic collision tests including the declared overlap allowance.
assert clear(np.array([[[3.,0,0]]]),np.array([1.7]),trees(np.zeros((1,3)),np.array(["C"])))[0]>=-1e-12
assert clear(np.array([[[2.9,0,0]]]),np.array([1.7]),trees(np.zeros((1,3)),np.array(["C"])))[0]<0

templates={}
donor_templates=[]
for name,pdb,res,donors in [("IsPETase","5XJH",[160,237,206],[87,161]),("LCC","4EB0",[165,242,210],[95,166])]:
    rec=load_pdb(P/"catalytic_motif_compare_v1"/(pdb+".pdb"),True)
    aa={(int(l[22:26]),l[12:16].strip()):(l,q,e) for l,q,e in rec}
    tri=[(k,v) for k,v in aa.items() if k[0] in res]
    # Fit/placement uses whole three residues, not isolated functional points.
    keys=[k for k,v in tri]; tq=np.array([v[1] for k,v in tri]); te=np.array([v[2] for k,v in tri])
    og=keys.index((res[0],"OG"));cb=keys.index((res[0],"CB"))
    tq=tq-tq[og]
    templates[name]=(tq,te,og,cb,[v[0] for k,v in tri])
    if name=="IsPETase":
        for d in donors:
            keys_d=[(d-1,"C"),(d-1,"O"),(d,"N"),(d,"CA"),(d,"C"),(d,"O")]
            dq=np.array([aa[k][1] for k in keys_d]);de=np.array([aa[k][2] for k in keys_d]);dq-=dq[2]
            hv=-norm(norm(dq[0])+norm(dq[3]))
            z=hv; x=norm(np.cross(z,[1.,0,0]) if abs(z[0])<.9 else np.cross(z,[0.,1,0]));y=np.cross(z,x)
            donor_templates.append((dq,np.column_stack([x,y,z]),de,[aa[k][0] for k in keys_d]))
assert all(len(t[0])==24 for t in templates.values())

def site_frame(row):
    ids=[lookup[(row["segment"],row[k])] for k in ["carbonyl_C_name","carbonyl_O_name","ester_O_name"]]
    ester=xyz[ids].copy(); c=ester[0].copy()
    ester[:,:2]-=np.round((ester[:,:2]-c[:2])/BOX[:2])*BOX[:2]
    ex=norm(ester[1]-c);ey=ester[2]-c;ey=norm(ey-ex*np.dot(ex,ey));ez=np.cross(ex,ey)
    return ester,np.column_stack([ex,ey,ez]),ids

def make_triads(ester,frame,template,env):
    tq,te,og,cb,_=template;N=CFG["n_triad_trials"]
    ct=rng.uniform(np.cos(np.deg2rad(115)),np.cos(np.deg2rad(95)),N);phi=rng.uniform(0,2*np.pi,N)
    direction=np.column_stack([ct,np.sqrt(1-ct*ct)*np.cos(phi),np.sqrt(1-ct*ct)*np.sin(phi)])@frame.T
    dist=rng.uniform(3,3.4,N);pos=ester[0]+dist[:,None]*direction
    rot=Rotation.random(N,random_state=rng).as_matrix()
    q=np.einsum("aj,nkj->nak",tq,rot)+pos[:,None,:]
    cosine=np.sum(norm(q[:,cb]-q[:,og])*(-direction),axis=1)
    angle=np.degrees(np.arccos(np.clip(cosine,-1,1)));geometry=(angle>=100)&(angle<=130)
    cand=q[geometry];margin=clear(cand,np.array([RAD[e] for e in te]),env)
    good=cand[margin>=0]
    keep=np.linspace(0,len(good)-1,min(len(good),64),dtype=int) if len(good) else np.array([],int)
    picked=good[keep]
    if len(picked):
        assert np.allclose(np.linalg.norm(picked[0,:,None]-picked[0,None,:],axis=-1),np.linalg.norm(tq[:,None]-tq[None,:],axis=-1),atol=1e-8)
    return picked,{"trials":N,"attack_and_Ser_angle_pass":int(geometry.sum()),"polymer_clear":len(good),"poses_evaluated":len(picked)}

def donor_library(ester,frame,template,env):
    dq,F,de,_=template;N=CFG["n_donor_trials_per_type"]
    u=norm(rng.normal(size=(N,3)));rr=rng.uniform(2.7,3.2,N);Npos=ester[1]+rr[:,None]*u
    aim=-u;tmp=np.tile([1.,0,0],(N,1));tmp[np.abs(aim[:,0])>.9]=[0,1,0]
    x=norm(np.cross(aim,tmp));y=np.cross(aim,x)
    ct=rng.uniform(np.cos(np.deg2rad(40)),1,N);phi=rng.uniform(0,2*np.pi,N)
    hdir=ct[:,None]*aim+np.sqrt(1-ct*ct)[:,None]*(np.cos(phi)[:,None]*x+np.sin(phi)[:,None]*y)
    tmp=np.tile([1.,0,0],(N,1));tmp[abs(hdir[:,0])>.9]=[0,1,0]
    bx=norm(np.cross(hdir,tmp));by=np.cross(hdir,bx);roll=rng.uniform(0,2*np.pi,N)
    dx=np.cos(roll)[:,None]*bx+np.sin(roll)[:,None]*by;dy=np.cross(hdir,dx)
    D=np.stack([dx,dy,hdir],axis=2);rot=np.einsum("nij,kj->nik",D,F)
    q=np.einsum("aj,nkj->nak",dq,rot)+Npos[:,None]
    H=Npos+1.01*hdir
    angle=np.degrees(np.arccos(np.clip(np.sum(norm(Npos-H)*norm(ester[1]-H),axis=1),-1,1)))
    hb=angle>=140
    geo_ids=np.where(hb)[0];margin=clear(q[hb],np.array([RAD[e] for e in de]),env);ids=geo_ids[margin>=0]
    local_dir=u@frame
    bins_phi=np.linspace(-180,180,25);bins_z=np.linspace(-1,1,13)
    ph=np.degrees(np.arctan2(local_dir[:,2],local_dir[:,1]))
    hist=np.histogram2d(ph[hb],local_dir[hb,0],bins=[bins_phi,bins_z])[0]
    return {"q":q[ids],"H":H[ids],"direction":local_dir[ids],"r":rr[ids],"angle":angle[ids],"roll_deg":np.degrees(roll[ids]),"hist_den":hist,"geom_trials":int(hb.sum()),"trials":N,"polymer_clear":len(ids)}

def write_pdb(path,records,coords,chain=None):
    lines=[]
    for i,(rec,q) in enumerate(zip(records,coords)):
        if i and rec[22:26]!=records[i-1][22:26] and int(rec[22:26])>int(records[i-1][22:26])+1:lines.append("TER")
        text=rec[:30]+"".join(f"{v:8.3f}" for v in q)+rec[54:]
        if chain:text=text[:21]+chain+text[22:]
        lines.append(text)
    path.write_text("\n".join(lines)+"\nTER\nEND\n")

def run_case(row,name,control=False):
    ester,frame,ids=site_frame(row); center=ester[0]
    image_xyz=xyz.copy();image_xyz[:,:2]-=np.round((image_xyz[:,:2]-center[:2])/BOX[:2])*BOX[:2]
    if control:qenv=image_xyz[ids];eenv=elems[ids]
    else:qenv=image_xyz;eenv=elems
    env=trees(qenv,eenv);t=templates[name]
    poses,stats=make_triads(ester,frame,t,env)
    ds=[donor_library(ester,frame,d,env) for d in donor_templates]
    ma=np.zeros((len(poses),len(ds[0]["q"])),bool);mb=np.zeros((len(poses),len(ds[1]["q"])),bool)
    pose_stats=[];pair_samples=[];best=None;sampled_pairs_total=0;feasible_pairs_total=0
    for k,pose in enumerate(poses):
        tenv=trees(pose,t[1])
        for z,m in enumerate([ma,mb]):
            if len(ds[z]["q"]):m[k]=clear(ds[z]["q"],np.array([RAD[e] for e in donor_templates[z][2]]),tenv)>=0
        pools=[np.flatnonzero(ma[k]),np.flatnonzero(mb[k])]
        ix=[rng.choice(pool,min(len(pool),64),replace=False) for pool in pools]
        tested=len(ix[0])*len(ix[1]);count=0
        if tested:
            a,b=ds[0]["q"][ix[0]],ds[1]["q"][ix[1]]
            delta=a[:,None,:,None,:]-b[None,:,None,:,:]
            d=np.sqrt(np.sum(delta*delta,axis=-1))
            rad_a=np.array([RAD[e] for e in donor_templates[0][2]]);rad_b=np.array([RAD[e] for e in donor_templates[1][2]])
            margin=(d-rad_a[None,None,:,None]-rad_b[None,None,None,:]+.4).min(axis=(2,3))
            valid=np.argwhere(margin>=0);count=len(valid)
            for ia,ib in valid[:20]:pair_samples.append([k,int(ix[0][ia]),int(ix[1][ib])])
            if len(valid) and best is None:
                ia,ib=valid[0];best=(k,int(ix[0][ia]),int(ix[1][ib]))
        sampled_pairs_total+=tested;feasible_pairs_total+=count
        pose_stats.append({"pose":k,"A_feasible":len(pools[0]),"B_feasible":len(pools[1]),"A_pair_subsample":len(ix[0]),"B_pair_subsample":len(ix[1]),"pairs_tested":tested,"pairs_feasible":count})
    folder=OUT/(row["site_id"]+"_"+name+("_ester3_control" if control else ""));folder.mkdir()
    arrays={"triads_xyz":poses,"ester_xyz":ester,"ester_frame_columns":frame,"mask_A":ma,"mask_B":mb,"pair_examples":np.array(pair_samples,dtype=int).reshape(-1,3)}
    for z,tag in enumerate(["A","B"]):
        for key in ["q","H","direction","r","angle","roll_deg","hist_den"]:arrays[tag+"_"+key]=ds[z][key]
    np.savez_compressed(folder/"conditional_space.npz",**arrays)
    res={"site_id":row["site_id"],"template":name,"control":control,"carbon_area_A2":float(row["carbonyl_C_surface_area_A2"]),"exposure_rank_group":row["exposure_rank_group"],"roughness_RMS_A":float(row["roughness_plane_normal_rms_A"]),"bulk_support":row["inward_support"],"triad":stats,"donor_A":{k:ds[0][k] for k in ["trials","geom_trials","polymer_clear"]},"donor_B":{k:ds[1][k] for k in ["trials","geom_trials","polymer_clear"]},"poses_with_independent_A_and_B":sum(x["A_feasible"]>0 and x["B_feasible"]>0 for x in pose_stats),"poses_with_sampled_compatible_pair":sum(x["pairs_feasible"]>0 for x in pose_stats),"sampled_pairs_tested":sampled_pairs_total,"sampled_pairs_feasible":feasible_pairs_total,"status":"FOUND_GEOMETRIC_COMBINATIONS" if best else "NOT_FOUND_IN_FINITE_SAMPLE","per_pose":pose_stats}
    if best:
        k,ia,ib=best
        # Independent replay of the selected pair checks geometry and every environment.
        for j,idx in enumerate([ia,ib]):
            q=ds[j]["q"][idx:idx+1];hh=ds[j]["H"][idx]
            assert clear(q,np.array([RAD[e] for e in donor_templates[j][2]]),env)[0]>=-1e-8
            assert clear(q,np.array([RAD[e] for e in donor_templates[j][2]]),trees(poses[k],t[1]))[0]>=-1e-8
            assert 2.7-1e-8<=np.linalg.norm(q[0,2]-ester[1])<=3.2+1e-8
            assert ds[j]["angle"][idx]>=140
        write_pdb(folder/"triad.pdb",t[4],poses[k])
        write_pdb(folder/"donor_A.pdb",donor_templates[0][3],ds[0]["q"][ia],chain="D")
        write_pdb(folder/"donor_B.pdb",donor_templates[1][3],ds[1]["q"][ib],chain="E")
        near=np.linalg.norm(image_xyz-center,axis=1)<=20
        write_pdb(folder/"PET_context_20A.pdb",[raw[z][0] for z in np.where(near)[0]],image_xyz[near])
        write_pdb(folder/"target_ester.pdb",[raw[z][0] for z in ids],ester)
        res["example_indices"]={"pose":k,"donor_A":ia,"donor_B":ib}
    (folder/"summary.json").write_text(json.dumps(res,indent=2))
    print(json.dumps({k:v for k,v in res.items() if k!="per_pose"}),flush=True)
    return res,arrays

OUT.mkdir(exist_ok=False)
(OUT/"config.json").write_text(json.dumps(CFG,indent=2))
(OUT/"selected_sites.json").write_text(json.dumps({"eligible_count":len(eligible),"selection":"top surface, positive C area, <=3 A from dense envelope, >=90% support 5 A inward, no XY edges; six exposure quantiles","sites":selected},indent=2))
results=[];plotdata=[]
for row in selected:
    for name in templates:
        res,arr=run_case(row,name);results.append(res)
        if name=="IsPETase":plotdata.append((res,arr))
res,_=run_case(selected[0],"IsPETase",control=True);results.append(res)
compact=[{k:v for k,v in r.items() if k!="per_pose"} for r in results]
(OUT/"summary.json").write_text(json.dumps(compact,indent=2))
fig,axes=plt.subplots(2,3,figsize=(14,8),constrained_layout=True)
for ax,(res,arr) in zip(axes.ravel(),plotdata):
    if len(arr["mask_A"]):
        d=arr["A_direction"];ph=np.degrees(np.arctan2(d[:,2],d[:,1]));counts=arr["mask_A"].sum(axis=0)
        hist=np.histogram2d(ph,d[:,0],bins=[np.linspace(-180,180,25),np.linspace(-1,1,13)],weights=counts)[0]
        denom=arr["A_hist_den"]*len(arr["mask_A"])
        rate=np.divide(hist,denom,out=np.zeros_like(hist),where=denom>0)
        im=ax.imshow(rate.T,origin="lower",extent=[-180,180,-1,1],aspect="auto",vmin=0,vmax=1,cmap="viridis")
    else:
        im=ax.imshow(np.zeros((12,24)),origin="lower",extent=[-180,180,-1,1],aspect="auto",vmin=0,vmax=1,cmap="viridis")
        ax.text(0,0,"NO TRIAD POSE FOUND",ha="center",color="white")
    ax.set_title(res["site_id"]+" | C area %.2f A2"%res["carbon_area_A2"]);ax.set_xlabel("Azimuth around C=O axis (deg)");ax.set_ylabel("cos(polar angle to C->O)")
fig.colorbar(im,ax=axes.ravel().tolist(),label="Feasible sampled A fragments / tested (HB-valid fragments x retained poses)")
fig.suptitle("Reactant-state geometric pilot: independent donor A | IsPETase triad\nSampling acceptance, NOT binding probability or catalytic activity")
fig.savefig(OUT/"donor_direction_maps.png",dpi=140);plt.close(fig)
(OUT/"RUNBOOK.md").write_text("CPU-only pilot. Run the immutable source with existing unbiased_100chain_v1/env/bin/python -B. See config.json for all search windows and steric overlap allowance. The polymer is all original heavy atoms, with minimum-image XY coordinates about each site. No chain atoms are deleted for collision checks. Bulk support is a density heuristic. Donors are six-heavy-atom peptide fragments plus an idealized virtual NH; they are independent, have no imposed linkage, and are checked for simultaneous non-overlap. Terminal peptide connectivity, full sidechains, backbone completion, full-enzyme access, water/TS states, electrostatics, dynamics and catalysis remain NOT_EVALUATED. Candidate direction maps are proposal-dependent finite-sampling acceptance fractions. Triad poses are limited to 64 and donor pair tests to 64x64 per pose. No-solution means NOT_FOUND_IN_FINITE_SAMPLE, not proven impossible. Distances/angles are initial declared search windows, not fitted enzyme activity thresholds. Control retains only target ester C/O/O steric obstacles, not a chemically complete molecule.\n")
manifest={str(f.relative_to(OUT)):hashlib.sha256(f.read_bytes()).hexdigest() for f in OUT.rglob("*") if f.is_file()}
manifest["SOURCE_SHA256"]=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
manifest["INPUT_PDB_SHA256"]=hashlib.sha256(rawpath.read_bytes()).hexdigest()
(OUT/"SHA256.json").write_text(json.dumps(manifest,indent=2))
with (OUT/"RUN_LOG.jsonl").open("a") as f:f.write(json.dumps({"time":datetime.datetime.now(datetime.timezone.utc).isoformat(),"status":"TECHNICAL_COMPLETE","cases":len(results),"config":CFG})+"\n")
print("PILOT_COMPLETE",len(results),flush=True)
