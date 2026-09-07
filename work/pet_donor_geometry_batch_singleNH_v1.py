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
OUT=P/"donor_geometry_batch_singleNH_v1"
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


from scipy.stats import qmc
from scipy.spatial.distance import cdist
import time
BATCH_CFG={"version":"singleNH_v1","triad_trials_per_template":12000,"max_triad_poses":64,
 "donor_template":"IsPETase Tyr87 peptide: C86 O86 N87 CA87 C87 O87; virtual NH=1.01 A",
 "donor_sampler":"scrambled Sobol 6D, common nested prefixes; direct valid NHO construction",
 "donor_stages":[2048,4096,8192],"budget_rule":"all sites with sampled triad poses evaluated at 2048 and 4096; extend to 8192 if either template with poses has no pair",
 "support_thresholds":[0.9,0.8],"other_geometry":CFG,
 "collision_environment":"original heavy atoms, minimum-image XY; exact safe 30 A neighborhood after template-size bound check",
 "output":"statistics and angular histograms only; no per-site PDB",
 "interpretation":"sampling-dependent geometry, not activity, TS stabilization, complete protein realizability, or absence proof"}
assert max(np.linalg.norm(t[0],axis=1).max()+3.4 for t in templates.values())+3.6<30
assert np.linalg.norm(donor_templates[0][0],axis=1).max()+3.2+1.5+3.6<30
def seed_for(text):return int.from_bytes(hashlib.sha256(text.encode()).digest()[:4],"little")

def collect():
    byside={}
    for side in ("top","bottom"):
        hh=zc[np.max(np.where(dense,np.arange(n[2]),0),axis=2)] if side=="top" else zc[np.min(np.where(dense,np.arange(n[2]),n[2]-1),axis=2)]
        dd=np.array([(x,y,-5 if side=="top" else 5) for x in np.arange(-10,10.1,2) for y in np.arange(-10,10.1,2) if x*x+y*y<=100])
        got=[]
        with (P/"surface_chemistry_shape_20260907_v1/results_exterior_v3_uncertainty"/("esters_"+side+"_15A.csv")).open() as f:
            for row in csv.DictReader(f):
                if row["shape_scope_status"]!="COMPLETE_FOR_CLASSIFIED_MESH" or float(row["carbonyl_C_surface_area_A2"])<=0:continue
                c=np.array([float(row["center_"+a+"_A"]) for a in "xyz"])
                h=float(map_coordinates(hh,(c[:2]/step[:2]-.5)[:,None],order=1,mode="grid-wrap")[0])
                support=float(np.mean(map_coordinates(rho,((c+dd)/step-.5).T,order=1,mode="grid-wrap")>=threshold))
                if abs(c[2]-h)<=3 and support>=.8:
                    row.update(side=side,support_fraction=support,boundary_offset_A=float(c[2]-h))
                    got.append(row)
        byside[side]=got
    rows=byside["top"]+byside["bottom"]
    assert len(rows)==len({r["site_id"] for r in rows}),"Unexpected cross-side duplicate: adjudication required"
    return sorted(rows,key=lambda r:(r["side"],r["site_id"]))

def sample_donor(ester,frame,seed):
    u=qmc.Sobol(d=6,scramble=True,seed=seed).random_base2(13)
    co=2*u[:,0]-1;ph=2*np.pi*u[:,1]
    direction=np.column_stack([co,np.sqrt(1-co*co)*np.cos(ph),np.sqrt(1-co*co)*np.sin(ph)])@frame.T
    rr=2.7+.5*u[:,2];Npos=ester[1]+rr[:,None]*direction;aim=-direction
    gmin=np.deg2rad(140);bmax=np.pi-gmin-np.arcsin(1.01/rr*np.sin(gmin))
    beta=np.arccos(1-u[:,3]*(1-np.cos(bmax)))
    tmp=np.tile([1.,0,0],(len(u),1));tmp[np.abs(aim[:,0])>.9]=[0,1,0]
    bx=norm(np.cross(aim,tmp));by=np.cross(aim,bx);phi=2*np.pi*u[:,4]
    hdir=np.cos(beta)[:,None]*aim+np.sin(beta)[:,None]*(np.cos(phi)[:,None]*bx+np.sin(phi)[:,None]*by)
    tmp=np.tile([1.,0,0],(len(u),1));tmp[np.abs(hdir[:,0])>.9]=[0,1,0]
    bx=norm(np.cross(hdir,tmp));by=np.cross(hdir,bx);roll=2*np.pi*u[:,5]
    dx=np.cos(roll)[:,None]*bx+np.sin(roll)[:,None]*by;dy=np.cross(hdir,dx)
    D=np.stack([dx,dy,hdir],axis=2);dq,F,de,_=donor_templates[0];rot=np.einsum("nij,kj->nik",D,F)
    q=np.einsum("aj,nkj->nak",dq,rot)+Npos[:,None];H=Npos+1.01*hdir
    angles=np.degrees(np.arccos(np.clip(np.sum(norm(Npos-H)*norm(ester[1]-H),axis=1),-1,1)))
    assert angles.min()>=140-1e-8 and angles.max()<=180
    assert np.allclose(np.linalg.norm(q[:,2]-ester[1],axis=1),rr)
    assert np.allclose(np.linalg.norm(q[0,:,None]-q[0,None,:],axis=-1),np.linalg.norm(dq[:,None]-dq[None,:],axis=-1))
    return q,H,angles,(direction@frame),(hdir@frame),np.degrees(roll)

def compatibility(q):
    N=len(q);adj=np.zeros((N,N),bool)
    if N<2:return adj
    rr=np.array([RAD[e] for e in donor_templates[0][2]])
    # Necessary N-N separation eliminates pairs before testing all 36 atom pairs.
    ij=np.argwhere(np.triu(cdist(q[:,2],q[:,2])>=2*RAD["N"]-.4,k=1))
    for offset in range(0,len(ij),2048):
        z=ij[offset:offset+2048]
        d=np.linalg.norm(q[z[:,0],:,None,:]-q[z[:,1],None,:,:],axis=-1)
        ok=(d-rr[None,:,None]-rr[None,None,:]+.4).min(axis=(1,2))>=0
        a,b=z[ok].T;adj[a,b]=True;adj[b,a]=True
    assert not adj.diagonal().any() and np.array_equal(adj,adj.T)
    return adj

def hist(v,weights=None):
    return np.histogram2d(np.degrees(np.arctan2(v[:,2],v[:,1])),v[:,0],
      bins=[np.linspace(-180,180,25),np.linspace(-1,1,13)],weights=weights)[0]

def evaluate(row,root):
    global rng
    key=row["site_id"];path=root/"cases"/(key+".json")
    if path.exists():return json.loads(path.read_text())
    tic=time.monotonic();ester,frame,ids=site_frame(row);center=ester[0]
    allq=xyz.copy();allq[:,:2]-=np.round((allq[:,:2]-center[:2])/BOX[:2])*BOX[:2]
    local=np.linalg.norm(allq-center,axis=1)<=30;env=trees(allq[local],elems[local])
    rr=np.array([RAD[e] for e in donor_templates[0][2]])
    tposes={};tstats={}
    for name in templates:
        rng=np.random.default_rng(seed_for(key+"|"+name+"|triad"))
        tposes[name],tstats[name]=make_triads(ester,frame,templates[name],env)
    records=[];arrays={};examples={};stages=[];stage_results={}
    if sum(len(v) for v in tposes.values()):
        q,H,angles,directions,hdir,roll=sample_donor(ester,frame,seed_for(key+"|singleNH|Sobol"))
        accepted=[];last=0
        for budget in [2048,4096,8192]:
            clear_ids=np.arange(last,budget)[clear(q[last:budget],rr,env)>=0];accepted.extend(clear_ids.tolist());last=budget
            inds=np.array(accepted,dtype=int);cand=q[inds]
            masks={};union=np.zeros(len(inds),bool)
            for name,poses in tposes.items():
                mm=np.zeros((len(poses),len(inds)),bool)
                for j,pose in enumerate(poses):
                    if len(inds):mm[j]=clear(cand,rr,trees(pose,templates[name][1]))>=0
                masks[name]=mm
                if len(mm):union|=mm.any(axis=0)
            uids=np.flatnonzero(union);adj=compatibility(cand[uids])
            this={}
            for name,mm in masks.items():
                counts=[];donors=[];pairs=0;posehit=0;participant=np.zeros(len(inds),bool)
                for j in range(len(mm)):
                    pos=np.flatnonzero(mm[j,uids]);sub=adj[np.ix_(pos,pos)];ij=np.argwhere(np.triu(sub,k=1))
                    cc=len(ij);counts.append(cc);donors.append(int(mm[j].sum()));pairs+=cc
                    if cc:
                        posehit+=1;participants=uids[pos[np.unique(ij)]];participant[participants]=True
                        if name not in examples:
                            aa,bb=uids[pos[ij[0]]];examples[name]={"pose":j,"sample_i":int(inds[aa]),"sample_j":int(inds[bb])}
                this[name]={"budget":budget,"polymer_clear_donors":len(inds),"poses_with_pair":posehit,"pose_pair_combinations":pairs,"donors_per_pose":donors,"unordered_pairs_per_pose":counts,"distinct_participating_donor_poses":int(participant.sum())}
                if len(mm):
                    arrays[name+"_N_position_hist"]=hist(directions[inds],mm.sum(0))
                    arrays[name+"_NH_direction_hist"]=hist(hdir[inds],mm.sum(0))
                    arrays[name+"_roll_hist"]=np.histogram(roll[inds],bins=np.linspace(0,360,25),weights=mm.sum(0))[0]
                    arrays[name+"_pair_participant_N_hist"]=hist(directions[inds][participant])
                    arrays[name+"_tested_position_hist"]=hist(directions[:budget])*len(mm)
            stages.append({"budget":budget,"templates":this});stage_results=this
            if budget==4096 and all(not len(tposes[name]) or this[name]["poses_with_pair"]>0 for name in templates):break
        # Replay first discovered pair for each template, independent of adjacency construction.
        for name,ex in examples.items():
            aq=q[ex["sample_i"]];bq=q[ex["sample_j"]];pose=tposes[name][ex["pose"]]
            margin=float((cdist(aq,bq)-rr[:,None]-rr[None,:]+.4).min())
            assert margin>=-1e-8
            for frag in [aq,bq]:
                assert clear(frag[None],rr,env)[0]>=-1e-8
                assert clear(frag[None],rr,trees(pose,templates[name][1]))[0]>=-1e-8
            ex.update(pair_clearance_A=margin,N_O_A=[float(np.linalg.norm(aq[2]-ester[1])),float(np.linalg.norm(bq[2]-ester[1]))],N_H_O_deg=[float(angles[ex["sample_i"]]),float(angles[ex["sample_j"]])])
    for name in templates:
        final=stage_results.get(name,{"budget":0,"polymer_clear_donors":0,"poses_with_pair":0,"pose_pair_combinations":0,"donors_per_pose":[],"unordered_pairs_per_pose":[],"distinct_participating_donor_poses":0})
        records.append({"template":name,"triad":tstats[name],"donor":final,"status":"FOUND_GEOMETRIC_COMBINATIONS" if final["poses_with_pair"] else "NOT_FOUND_IN_FINITE_SAMPLE","donor_stage_status":"EVALUATED" if len(tposes[name]) else "NOT_EVALUATED_NO_SAMPLED_TRIAD_POSE","example_audit":examples.get(name)})
    result={"site_id":key,"side":row["side"],"support_fraction":row["support_fraction"],"boundary_offset_A":row["boundary_offset_A"],"carbon_area_A2":float(row["carbonyl_C_surface_area_A2"]),"roughness_RMS_A":float(row["roughness_plane_normal_rms_A"]),"templates":records,"sampling_stages":stages,"elapsed_s":time.monotonic()-tic}
    if arrays:np.savez_compressed(root/"cases"/(key+"_angular_histograms.npz"),**arrays)
    tmp=path.with_suffix(".tmp");tmp.write_text(json.dumps(result,separators=(",",":")));tmp.replace(path)
    with (root/"RUN_LOG.jsonl").open("a") as f:f.write(json.dumps({"time":datetime.datetime.now(datetime.timezone.utc).isoformat(),"site":key,"status":"COMPLETE","seconds":result["elapsed_s"]})+"\n")
    return result

def summarize(rows,root,label_name):
    results=[json.loads((root/"cases"/(r["site_id"]+".json")).read_text()) for r in rows]
    summary={"cohort":label_name,"sites":len(results),"top":sum(r["side"]=="top" for r in results),"bottom":sum(r["side"]=="bottom" for r in results),
      "sites_with_any_template_pair":sum(any(t["donor"]["poses_with_pair"] for t in r["templates"]) for r in results),
      "sites_with_any_sampled_triad":sum(any(t["triad"]["polymer_clear"] for t in r["templates"]) for r in results),
      "templates":{},"total_case_seconds":sum(r["elapsed_s"] for r in results),"status":"TECHNICAL_COMPLETE_GEOMETRY_ONLY"}
    for name in templates:
        ts=[next(t for t in r["templates"] if t["template"]==name) for r in results]
        summary["templates"][name]={"site_count":len(ts),"sites_with_pair":sum(t["donor"]["poses_with_pair"]>0 for t in ts),"sites_with_triad":sum(t["triad"]["polymer_clear"]>0 for t in ts),"sampled_triad_poses":sum(t["triad"]["poses_evaluated"] for t in ts),"poses_with_pair":sum(t["donor"]["poses_with_pair"] for t in ts)}
    (root/(label_name+"_summary.json")).write_text(json.dumps(summary,indent=2))
    with (root/(label_name+"_sites.tsv")).open("w",newline="") as f:
        fields=["site_id","side","support_fraction","carbon_area_A2","roughness_RMS_A","template","triad_clear","triad_poses_tested","donor_budget","poses_with_pair","pose_pair_combinations","status","donor_stage_status"]
        w=csv.DictWriter(f,fieldnames=fields,delimiter="\t");w.writeheader()
        for r in results:
            for t in r["templates"]:
                w.writerow({"site_id":r["site_id"],"side":r["side"],"support_fraction":r["support_fraction"],"carbon_area_A2":r["carbon_area_A2"],"roughness_RMS_A":r["roughness_RMS_A"],"template":t["template"],"triad_clear":t["triad"]["polymer_clear"],"triad_poses_tested":t["triad"]["poses_evaluated"],"donor_budget":t["donor"]["budget"],"poses_with_pair":t["donor"]["poses_with_pair"],"pose_pair_combinations":t["donor"]["pose_pair_combinations"],"status":t["status"],"donor_stage_status":t["donor_stage_status"]})
    print("COHORT_COMPLETE",json.dumps(summary),flush=True)
    return summary

rows=collect()
smoke="--smoke" in sys.argv
root=OUT/("smoke" if smoke else "production");root.mkdir(parents=True,exist_ok=True);(root/"cases").mkdir(exist_ok=True)
sourcehash=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
config=dict(BATCH_CFG,source_sha256=sourcehash,input_sha256=hashlib.sha256(rawpath.read_bytes()).hexdigest())
configpath=root/"config.json"
if configpath.exists():assert json.loads(configpath.read_text())==config,"Resume config mismatch"
else:configpath.write_text(json.dumps(config,indent=2))
(root/"eligible_sites.json").write_text(json.dumps(rows,indent=2))
selection={"support90":sum(r["support_fraction"]>=.9 for r in rows),"support80":len(rows),"new_at80":sum(r["support_fraction"]<.9 for r in rows)}
print("COUNTS",json.dumps(selection),flush=True)
(root/"RUNBOOK.md").write_text("CPU-only, existing NumPy/SciPy environment. Run python -B this_source.py --smoke then python -B this_source.py --run. Config and source SHA must match on resume. Results are checkpointed after each site. 90% cohort completes before 80% additions; overlap reuses identical per-site results. 81-point inward XY disk, 10 A radius, 5 A inward Z, >=half central density after sigma=3 A smoothing; <=3 A from envelope, positive carbon area, complete classified surface, periodic XY included. This is a bulk-support heuristic, not a validated dynamic tail classifier. Fixed source-derived NH peptide template, two independent instances with no linkage constraint. No PDB generation. Every unordered surviving pair is checked at each evaluated budget; counts are pose-dependent, not biochemical probabilities. All unchanged chemistry/steric windows are in config. Atom typing uses heavy atoms and allows 0.4 A VDW overlap. Sidechain completion, full-protein access, TS/water states, electronic energies and catalysis remain NOT_EVALUATED. Sobol position and orientation density is the declared proposal distribution; adaptive budgets must be retained when comparing counts.\n")
if smoke:
    low=min((r for r in rows if r["side"]=="top" and r["support_fraction"]>=.9),key=lambda r:float(r["carbonyl_C_surface_area_A2"]))
    chosen=[low,next(r for r in rows if r["site_id"]=="P244_C42_O18"),next(r for r in rows if r["side"]=="bottom" and .8<=r["support_fraction"]<.9)]
    for r in chosen:
        z=evaluate(r,root);print("SMOKE_CASE",r["site_id"],[(t["template"],t["triad"]["polymer_clear"],t["donor"]["poses_with_pair"],t["donor"]["budget"]) for t in z["templates"]],flush=True)
    summarize(chosen,root,"smoke")
    print("SMOKE_PASS",flush=True)
else:
    smokeconfig=OUT/"smoke/config.json"
    assert smokeconfig.exists() and json.loads(smokeconfig.read_text())==config,"Smoke run required with same config"
    for lab,minsupport in [("support90",.9),("support80",.8)]:
        cohort=[r for r in rows if r["support_fraction"]>=minsupport]
        for i,row in enumerate(cohort):
            evaluate(row,root)
            if (i+1)%10==0 or i+1==len(cohort):print("PROGRESS",lab,i+1,len(cohort),flush=True)
        summarize(cohort,root,lab)
    hashes={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in root.rglob("*") if f.is_file()}
    (root/"SHA256.json").write_text(json.dumps(hashes,indent=2))
    print("BATCH_COMPLETE",flush=True)
