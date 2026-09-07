"""Native LCC/IsPETase triad + OWN backbone donors. CPU diagnostic only."""
from pathlib import Path
import json,hashlib,datetime,time
import numpy as np
from scipy.stats import qmc
from scipy.spatial.transform import Rotation
from scipy.spatial.distance import cdist
P=Path("/data/bht2/polymer_material_reference_20260827/simulation_slabs/pet_dp10_400chain_direct_v1")
OUT=P/"native_motif_surface_compare_v1";OUT.mkdir(exist_ok=False)
start=time.monotonic()
def log(event,**kw):
    with (OUT/"RUN_LOG.jsonl").open("a") as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),event=event,**kw))+"\n")
source=P/"pet_donor_geometry_batch_singleNH_v1.py"
assert hashlib.sha256(source.read_bytes()).hexdigest()=="b796be0506a0d696ee19aad9288e2b71b6c81131c30f7e4bf2deb00c94ddb80a"
g={"__name__":"parent_functions"};s=source.read_text();exec(compile(s[:s.index("rows=collect()")],str(source),"exec"),g)
norm=g["norm"]
def angle(a,b):return np.degrees(np.arccos(np.clip(np.sum(norm(a)*norm(b),axis=-1),-1,1)))
config={"seeds":[20260907,20260908],"trials_per_seed":1048576,"block":65536,"parent_geometry":g["CFG"],
 "method":"Two common scrambled Sobol 6D pose pools. Native motifs rigid. Distances/angles unchanged. Uniform attack cos(theta), azimuth, distance and Shoemake SO3.",
 "scope":"34 unique heavy atoms: full triad residues + OWN two backbone peptide fragments, with shared atoms deduplicated. Donor side chains absent consistently with previous fragment protocol. Virtual H for NH angle only.",
 "comparison":"Native LCC versus native IsPETase; not a rank insertion into the earlier 473 fixed placements.",
 "status":"DIAGNOSTIC_GEOMETRY_ONLY"}
log("START",config=config,script=str(Path(__file__)))
B=P/"donor_geometry_batch_singleNH_v1/extension70_v1"
rows=json.loads((B/"eligible_sites.json").read_text())
assert len(rows)==578 and len(set(r["site_id"] for r in rows))==578
esters=[];frames=[]
for r in rows:
    e,f,_=g["site_frame"](r);esters.append(e);frames.append(f)
E=np.array(esters);F=np.array(frames);localE=np.einsum("naj,njk->nak",E-E[:,:1],F)
lengths=localE[:,1,0];Lmin=float(lengths.min());Lmax=float(lengths.max())
assert np.abs(localE[:,1,1:]).max()<1e-10
config["target_carbonyl_lengths_A"]=[Lmin,Lmax]
specs=[("IsPETase","5XJH",[160,237,206],[87,161],"a0410e1407aa9ccf5fa16396a3fe329c46504fc4d6182fb438607ec23c697680"),("LCC","4EB0",[165,242,210],[95,166],"952288d4a49893ee259638bc72e17868b6903f0027c4bd4d207f76b9ddf85718")]
models={}
for name,pdb,res,ds,sha in specs:
    path=P/"catalytic_motif_compare_v1"/(pdb+".pdb");assert hashlib.sha256(path.read_bytes()).hexdigest()==sha
    rec=g["load_pdb"](path,True);aa={(int(l[22:26]),l[12:16].strip()):(l,q,e) for l,q,e in rec}
    keys=[k for k in aa if k[0] in res]
    for d in ds:
        for k in [(d-1,"C"),(d-1,"O"),(d,"N"),(d,"CA"),(d,"C"),(d,"O")]:
            if k not in keys:keys.append(k)
    q=np.array([aa[k][1] for k in keys]);og=keys.index((res[0],"OG"));cb=keys.index((res[0],"CB"));origin=q[og].copy();q-=origin
    nidx=[keys.index((d,"N")) for d in ds];H=[]
    for d in ds:
        n=aa[d,"N"][1]
        h=n-1.01*norm(norm(aa[d-1,"C"][1]-n)+norm(aa[d,"CA"][1]-n))
        H.append(h-origin)
    elements=np.array([aa[k][2] for k in keys]);radii=np.array([g["RAD"][e] for e in elements])
    assert len(q)==34
    models[name]={"q":q,"H":np.array(H),"keys":keys,"og":og,"cb":cb,"nidx":nidx,"elements":elements,"radii":radii,"res":res,"donors":ds,"pdb":pdb}
# Save an exact triad-aligned native overlay, with no fictitious substrate.
ref=models["IsPETase"];mov=models["LCC"]
def fit_ids(m):
    a,b,c=m["res"]
    return [m["keys"].index(k) for k in [(a,"CB"),(a,"OG"),(b,"CG"),(b,"ND1"),(b,"CD2"),(b,"CE1"),(b,"NE2"),(c,"CG"),(c,"OD1"),(c,"OD2")]]
A=mov["q"][fit_ids(mov)];T=ref["q"][fit_ids(ref)];ac=A.mean(0);tc=T.mean(0)
u,_,vt=np.linalg.svd((A-ac).T@(T-tc));rot=u@vt
if np.linalg.det(rot)<0:u[:,-1]*=-1;rot=u@vt
rmsd=float(np.sqrt(np.mean(np.sum(((A-ac)@rot+tc-T)**2,axis=1))))
assert abs(rmsd-.2503020615823227)<1e-8
scenes={"overlay_rmsd_A":rmsd,"native":[],"surface_examples":[]}
for name,m in models.items():
    q=m["q"].copy();h=m["H"].copy()
    if name=="LCC":q=(q-ac)@rot+tc;h=(h-ac)@rot+tc
    bonds=np.argwhere(np.triu((cdist(q,q)>.9)&(cdist(q,q)<1.85),1)).tolist()
    scenes["native"].append({"name":name,"pdb":m["pdb"],"q":q.tolist(),"H":h.tolist(),"labels":[str(a)+":"+b for a,b in m["keys"]],"elements":m["elements"].tolist(),"bonds":bonds,"og":m["og"],"nidx":m["nidx"],"donor_residues":m["donors"],"roles":[("Ser" if a==m["res"][0] else "His" if a==m["res"][1] else "Asp" if a==m["res"][2] else "donor") for a,b in m["keys"]]})
(OUT/"native_3d_data.json").write_text(json.dumps(scenes))
# Persist/replay the earlier small geometry-only preflight exactly.
preflight={}
for name,m in models.items():
    rng=np.random.default_rng(20260907);n=12000
    ct=rng.uniform(np.cos(np.deg2rad(115)),np.cos(np.deg2rad(95)),n);ph=rng.uniform(0,2*np.pi,n)
    dr=np.column_stack([ct,np.sqrt(1-ct*ct)*np.cos(ph),np.sqrt(1-ct*ct)*np.sin(ph)])
    dist=rng.uniform(3,3.4,n);pos=dist[:,None]*dr;rt=Rotation.random(n,random_state=rng).as_matrix()
    N=np.einsum("aj,nkj->nak",m["q"][m["nidx"]],rt)+pos[:,None]
    H=np.einsum("aj,nkj->nak",m["H"],rt)+pos[:,None]
    CB=np.einsum("j,nkj->nk",m["q"][m["cb"]],rt)+pos
    sa=angle(CB-pos,-pos);ok=(sa>=100)&(sa<=130);O=np.array([1.22,0,0])
    no=np.linalg.norm(N-O,axis=-1);ha=angle(N-H,O-H);each=(no>=2.7)&(no<=3.2)&(ha>=140)
    preflight[name]={"trials":n,"Ser_geometry":int(ok.sum()),"donor1":int((ok&each[:,0]).sum()),"donor2":int((ok&each[:,1]).sum()),"both":int((ok&each.all(1)).sum())}
(OUT/"preflight12000.json").write_text(json.dumps(preflight,indent=2))
# Master pools: screen with exact necessary distance bounds over observed C=O lengths.
pools={}
for name,m in models.items():
    savedq=[];savedh=[];seedids=[];trialids=[];diag=[]
    for seed in config["seeds"]:
        sampler=qmc.Sobol(6,scramble=True,seed=seed);ng=0;nd=0
        for b in range(config["trials_per_seed"]//config["block"]):
            v=sampler.random(config["block"])
            ct=np.cos(np.deg2rad(115))+v[:,1]*(np.cos(np.deg2rad(95))-np.cos(np.deg2rad(115)));ph=2*np.pi*v[:,2]
            dr=np.column_stack([ct,np.sqrt(1-ct*ct)*np.cos(ph),np.sqrt(1-ct*ct)*np.sin(ph)])
            pos=(3+.4*v[:,0,None])*dr
            quat=np.column_stack([np.sqrt(1-v[:,3])*np.sin(2*np.pi*v[:,4]),np.sqrt(1-v[:,3])*np.cos(2*np.pi*v[:,4]),np.sqrt(v[:,3])*np.sin(2*np.pi*v[:,5]),np.sqrt(v[:,3])*np.cos(2*np.pi*v[:,5])])
            rt=Rotation.from_quat(quat).as_matrix()
            CB=np.einsum("j,nkj->nk",m["q"][m["cb"]],rt)+pos;sa=angle(CB-pos,-pos);ix=np.flatnonzero((sa>=100)&(sa<=130));ng+=len(ix)
            rt=rt[ix];pos=pos[ix];N=np.einsum("aj,nkj->nak",m["q"][m["nidx"]],rt)+pos[:,None]
            near=np.clip(N[:,:,0],Lmin,Lmax);dmin=np.sqrt((N[:,:,0]-near)**2+(N[:,:,1:]**2).sum(-1))
            dmax=np.maximum(np.linalg.norm(N-np.array([Lmin,0,0]),axis=-1),np.linalg.norm(N-np.array([Lmax,0,0]),axis=-1))
            keep=np.all((dmin<=3.2+1e-8)&(dmax>=2.7-1e-8),axis=1);nd+=int(keep.sum())
            if keep.any():
                qq=np.einsum("aj,nkj->nak",m["q"],rt[keep])+pos[keep,None]
                hh=np.einsum("aj,nkj->nak",m["H"],rt[keep])+pos[keep,None]
                savedq.append(qq);savedh.append(hh);seedids.extend([seed]*len(qq));trialids.extend((b*config["block"]+ix[keep]).tolist())
        diag.append({"seed":seed,"trials":config["trials_per_seed"],"Ser_geometry":ng,"necessary_distance_pass":nd})
    Q=np.concatenate(savedq) if savedq else np.empty((0,34,3));HH=np.concatenate(savedh) if savedh else np.empty((0,2,3))
    pools[name]={"Q":Q,"H":HH,"seed":np.array(seedids),"trial":np.array(trialids),"diagnostics":diag}
    np.savez_compressed(OUT/(name+"_necessary_pose_pool.npz"),Q=Q,H=HH,seed=np.array(seedids),trial=np.array(trialids))
    print("POOL",name,len(Q),diag,flush=True)
records=[];examples={};audit=[]
for k,row in enumerate(rows):
    e=E[k];f=F[k];le=localE[k];allq=g["xyz"].copy();allq[:,:2]-=np.round((allq[:,:2]-e[0,:2])/g["BOX"][:2])*g["BOX"][:2]
    local=np.linalg.norm(allq-e[0],axis=1)<=30;env=g["trees"](allq[local],g["elems"][local]);isolated=g["trees"](le,np.array(["C","O","O"]))
    record={"site_id":row["site_id"],"support_fraction":row["support_fraction"],"side":row["side"],"templates":{}}
    for name,m in models.items():
        pool=pools[name];Q=pool["Q"];HH=pool["H"];N=Q[:,m["nidx"]];no=np.linalg.norm(N-le[1],axis=-1);ha=angle(N-HH,le[1]-HH)
        both=np.all((no>=2.7-1e-8)&(no<=3.2+1e-8)&(ha>=140-1e-8),axis=1);ix=np.flatnonzero(both)
        ids=ix[g["clear"](Q[ix],m["radii"],isolated)>=-1e-8] if len(ix) else np.array([],int)
        w=Q[ids]@f.T+e[0]
        margin=g["clear"](w,m["radii"],env) if len(ids) else np.array([])
        ok=ids[margin>=-1e-8]
        rep={str(seed):int(np.sum(pool["seed"][ok]==seed)) for seed in config["seeds"]}
        record["templates"][name]={"chemical_pass":len(ix),"isolated_ester_pass":len(ids),"PET_pass":len(ok),"per_seed_pass":rep,"status":"FOUND_NATIVE_GEOMETRY" if len(ok) else "NOT_FOUND_IN_FINITE_SAMPLE"}
        if len(ok) and (name not in examples or len(ok)>examples[name]["pose_count"]):
            which=int(ids[np.argmax(margin)]);q=Q[which];hh=HH[which]
            examples[name]={"name":name,"site_id":row["site_id"],"pose_count":len(ok),"seed":int(pool["seed"][which]),"trial":int(pool["trial"][which]),"q":q.tolist(),"H":hh.tolist(),"ester":le.tolist(),"min_margin_A":float(margin.max()),"N_O_A":no[which].tolist(),"N_H_O_deg":ha[which].tolist(),"pool_index":which}
        if len(ok) and sum(x["name"]==name for x in audit)<3:
            which=int(ok[0]);q=Q[which]@f.T+e[0]
            rr=np.array([g["RAD"][x] for x in g["elems"]]);exact=float((cdist(q,allq)-m["radii"][:,None]-rr[None,:]+.4).min())
            assert exact>=-1e-8
            audit.append({"name":name,"site":row["site_id"],"exact_full_PET_margin_A":exact})
    records.append(record)
    if k%100==0:print("SITE",k,flush=True)
support=np.array([r["support_fraction"] for r in rows]);summary={"config":config,"preflight12000":preflight,"sites":578,"cohorts":{},"audit":audit,"status":"TECHNICAL_COMPLETE_GEOMETRY_ONLY"}
for key,mask in [("support90",support>=.9),("support80",support>=.8),("support70_all",support>=.7)]:
    summary["cohorts"][key]={"sites":int(mask.sum()),"templates":{}}
    for name in models:
        rr=[r["templates"][name] for r,use in zip(records,mask) if use]
        summary["cohorts"][key]["templates"][name]={"sites_chemical":sum(r["chemical_pass"]>0 for r in rr),"sites_isolated_ester":sum(r["isolated_ester_pass"]>0 for r in rr),"sites_PET":sum(r["PET_pass"]>0 for r in rr),"PET_pose_total":sum(r["PET_pass"] for r in rr),"sites_PET_by_seed":{str(seed):sum(r["per_seed_pass"][str(seed)]>0 for r in rr) for seed in config["seeds"]}}
summary["pool_diagnostics"]={name:p["diagnostics"] for name,p in pools.items()}
summary["elapsed_s"]=time.monotonic()-start
scenes["surface_examples"]=list(examples.values());scenes["summary"]=summary
(OUT/"native_3d_data.json").write_text(json.dumps(scenes))
(OUT/"summary.json").write_text(json.dumps(summary,indent=2));(OUT/"sites.json").write_text(json.dumps(records))
(OUT/"config.json").write_text(json.dumps(config,indent=2))
(OUT/"RUNBOOK.md").write_text("# Native motif surface comparison v1\n\nUse simulation_slabs/unbiased_100chain_v1/env/bin/python with OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1. Run the saved script once; output must not exist. CPU only.\n\nTwo native structures: WT IsPETase 5XJH A and LCC 4EB0 A. Full catalytic triad plus each enzyme's own two backbone donor fragments, with native relative coordinates fixed. No donor side chains, full enzyme, substrate relaxation, energies or TS. Initial 12000 geometry-only trials at C=O 1.22 A were sparse, so two equally sized Sobol 2^20 diagnostic pools were used without relaxing any chemical/steric thresholds. All 578 PET sites evaluated exactly in their ester coordinate systems. The no-bulk isolated control retains three target ester atoms. This ranks two native motifs under finite sampling, not complete enzymes and not comparable to ranks of the prior 473 fixed poses.\n\nVisual native overlay is a 10-atom triad fit with no substrate present. Surface scenes, if any, are independently geometry-validated examples and preserve native donor layout.\n")
hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file() and p.name!="RUN_LOG.jsonl"}
(OUT/"SHA256.json").write_text(json.dumps(hashes,indent=2))
log("COMPLETE",status=summary["status"],seconds=summary["elapsed_s"])
print("NATIVE_COMPARISON_COMPLETE",json.dumps(summary),flush=True)
