"""PET donor-pair translation robustness, CPU only; existing candidate centers.
No geometric threshold relaxation, energy estimate, or global optimum claim.
"""
from pathlib import Path
import argparse,json,hashlib,time,datetime,sys
import numpy as np
from scipy.stats import qmc
from scipy.spatial.distance import cdist
ap=argparse.ArgumentParser()
ap.add_argument("--power",type=int,default=10)
ap.add_argument("--limit",type=int,default=473)
ap.add_argument("--tag",required=True)
args=ap.parse_args()
assert 4<=args.power<=13 and 1<=args.limit<=473
P=Path("/data/bht2/polymer_material_reference_20260827/simulation_slabs/pet_dp10_400chain_direct_v1")
source=P/"pet_donor_geometry_batch_singleNH_v1.py"
assert hashlib.sha256(source.read_bytes()).hexdigest()=="b796be0506a0d696ee19aad9288e2b71b6c81131c30f7e4bf2deb00c94ddb80a"
g={"__name__":"parent_functions"}
s=source.read_text();exec(compile(s[:s.index("rows=collect()")],str(source),"exec"),g)
base=P/"donor_geometry_batch_singleNH_v1/joint_transfer_map_v1"
files=[base/"joint_transfer_results.npz",base/"motifs.json",P/"donor_geometry_batch_singleNH_v1/extension70_v1/eligible_sites.json",source]
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
z=np.load(files[0]);meta=json.loads(files[1].read_text());rows=json.loads(files[2].read_text())
M=z["motif_xyz_A"];H=z["motif_H_A"];old=z["coverage"];support=z["site_support"]
assert len(rows)==578 and len(meta)==473 and np.array_equal(z["site_ids"],[r["site_id"] for r in rows])
RR=np.array([[g["RAD"][e] for e in m["elements"]] for m in meta]);dr=RR[0,24:30]
assert np.all(RR[:,24:30]==dr) and np.all(RR[:,30:36]==dr)
root=P/"donor_pair_robustness_v1"/args.tag;root.mkdir(parents=True,exist_ok=False)
t0=time.monotonic()
def log(event,**kw):
 with (root/"RUN_LOG.jsonl").open("a") as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),event=event,**kw))+"\n")
log("START",argv=sys.argv,input_sha256=hashes)
config={"candidate_centers":args.limit,"target_sites":578,"radius_A":1.0,"replicate_seeds":[2026090711,2026090712],"samples_per_replicate":2**args.power,"sampling":"scrambled Sobol 6D mapped to independent uniform-volume 3D unit balls; translate each entire 6-heavy-atom peptide and its virtual H; no rotations; triad fixed","thresholds":g["CFG"],"primary_score":"sum over all sites of empirical joint-valid perturbation fraction = mean sites covered by one common perturbed pair; not union of independently optimized sites","reported_scopes":["donor_only_no_triad","joint_with_own_fixed_triad"],"uncertainty":"two scrambled sequences and nested prefixes, not equilibrium probabilities or guaranteed full-ball robustness","input_sha256":hashes}
(root/"config.json").write_text(json.dumps(config,indent=2))
(root/"RUNBOOK.md").write_text("# PET donor-pair robustness v1\n\nRun on CPU from project directory using existing unbiased_100chain_v1/env/bin/python -B pet_donor_pair_robustness_v1.py --power "+str(args.power)+" --limit "+str(args.limit)+" --tag "+args.tag+". Output tag must be new. OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1. Parent files read-only.\n\nIndependent uniform 1 A ball translations of two intact donor fragments, fixed orientations and triad. Center candidates are the previous finite first-success library, not native donor layouts. Both NH geometries and heavy-atom clashes must pass simultaneously. Donor-only and joint scopes are kept separate. Primary robust score is expected site coverage under specified uniform perturbations, not physical thermal probability. No assertion of stability, full protein realizability, activity, or global optimization.\n")
def angle(a,b):
 return np.degrees(np.arccos(np.clip(np.sum(a*b,axis=-1)/np.linalg.norm(a,axis=-1)/np.linalg.norm(b,axis=-1),-1,1)))
def geometry(q,h,o):
 n=q[:,:,2];no=np.linalg.norm(n-o,axis=-1);ha=angle(n-h,o-h)
 return np.all((no>=2.7-1e-8)&(no<=3.2+1e-8)&(ha>=140-1e-8),axis=1)
def pairclear(q):
 d=np.linalg.norm(q[:,0,:,None,:]-q[:,1,None,:,:],axis=-1)
 return (d-dr[None,:,None]-dr[None,None,:]+.4).min(axis=(1,2))>=-1e-8
n=2**args.power
offsets=[]
for seed in config["replicate_seeds"]:
 u=qmc.Sobol(6,scramble=True,seed=seed).random_base2(args.power).reshape(n,2,3)
 cost=2*u[:,:,0]-1;phi=2*np.pi*u[:,:,1];radius=u[:,:,2]**(1/3)
 offsets.append(radius[:,:,None]*np.stack([cost,np.sqrt(1-cost*cost)*np.cos(phi),np.sqrt(1-cost*cost)*np.sin(phi)],axis=-1))
offsets=np.concatenate(offsets)
assert np.linalg.norm(offsets,axis=-1).max()<=1 and np.allclose(np.mean(np.linalg.norm(offsets,axis=-1)**3,axis=0),.5,atol=.04)
# Prepare actual PET neighborhoods in each target ester coordinate frame.
envs=[];oxygen=[];esters=[];triad_ok=np.zeros((args.limit,578),bool)
for j,r in enumerate(rows):
 ester,frame,ids=g["site_frame"](r);q=g["xyz"].copy()
 q[:,:2]-=np.round((q[:,:2]-ester[0,:2])/g["BOX"][:2])*g["BOX"][:2]
 local=np.linalg.norm(q-ester[0],axis=1)<=30
 env=g["trees"]((q[local]-ester[0])@frame,g["elems"][local]);envs.append(env)
 o=(ester[1]-ester[0])@frame;oxygen.append(o);esters.append((ester-ester[0])@frame)
 for i in range(args.limit):
  idx=meta[i]["functional_indices"];og=M[i,idx["og"]];cb=M[i,idx["cb"]]
  da=np.linalg.norm(og);aa=angle(og,o);sa=angle(cb-og,-og)
  triad_ok[i,j]=3-1e-8<=da<=3.4+1e-8 and 95-1e-8<=aa<=115+1e-8 and 100-1e-8<=sa<=130+1e-8
 which=np.flatnonzero(triad_ok[:,j])
 if len(which):triad_ok[which,j]&=g["clear"](M[which,:24],RR[which,:24],env)>=-1e-8
 if j%150==0:print("ENV",j,flush=True)
oxygen=np.array(oxygen)
assert np.linalg.norm(M,axis=-1).max()+1+3.6<30
counts=np.zeros((2,args.limit,2,578),np.uint16) # scope,candidate,replicate,site
prefix_counts=np.zeros((2,args.limit,2,578),np.uint16)
nominal=np.zeros((2,args.limit,578),bool)
sample_coverage=np.zeros((2,args.limit,2*n),np.uint16)
best_samples=[];audit=[];clouds={}
for i in range(args.limit):
 q0=M[i,24:].reshape(1,2,6,3);h0=H[i][None]
 q=q0+offsets[:,:,None,:];h=h0+offsets
 assert np.allclose(np.linalg.norm(q[:,:,2]-h,axis=-1),1.01,atol=1e-7)
 pair=pairclear(q);pc0=bool(pairclear(q0)[0]);assert pc0
 tenv=g["trees"](M[i,:24],np.array(meta[i]["elements"][:24]))
 tc=g["clear"](q.reshape(-1,12,3),RR[i,24:],tenv)>=-1e-8
 tc0=bool(g["clear"](q0.reshape(1,12,3),RR[i,24:],tenv)[0]>=-1e-8);assert tc0
 for j,env in enumerate(envs):
  o=oxygen[j];which=np.flatnonzero(pair&geometry(q,h,o))
  passed=np.zeros(2*n,bool)
  if len(which):passed[which]=g["clear"](q[which].reshape(-1,12,3),RR[i,24:],env)>=-1e-8
  joint=passed&tc&triad_ok[i,j]
  for scope,ok in enumerate([passed,joint]):
   byrep=ok.reshape(2,n);counts[scope,i,:,j]=byrep.sum(axis=1)
   prefix_counts[scope,i,:,j]=byrep[:,:n//2].sum(axis=1)
   sample_coverage[scope,i]+=ok.astype(np.uint16)
  nom=bool(geometry(q0,h0,o)[0] and pc0)
  if nom:nom=bool(g["clear"](q0.reshape(1,12,3),RR[i,24:],env)[0]>=-1e-8)
  nominal[0,i,j]=nom;nominal[1,i,j]=nom and tc0 and triad_ok[i,j]
 assert np.array_equal(nominal[1,i],old[i]),("ZERO_SHIFT_REPLAY_FAILED",meta[i]["id"])
 assert np.all(counts[1,i]<=counts[0,i])
 # Store only compact aggregate state; no per-site PDB flood.
 for scope in range(2):
  t=int(np.argmax(sample_coverage[scope,i]))
  best_samples.append(dict(scope=scope,candidate=meta[i]["id"],trial=t,coverage=int(sample_coverage[scope,i,t]),offsets_A=offsets[t].tolist()))
 # Deterministic independent full-PET verification of nominal + perturbed positive and negative.
 if i in [0,args.limit-1]:
  for shift in [-1,int(np.argmax(sample_coverage[1,i])),0]:
   cq=q0[0] if shift<0 else q[shift];ch=h0[0] if shift<0 else h[shift]
   for j in [0,int(np.argmax(counts[1,i].sum(axis=0)))]:
    ester,frame,_=g["site_frame"](rows[j]);full=g["xyz"].copy()
    full[:,:2]-=np.round((full[:,:2]-ester[0,:2])/g["BOX"][:2])*g["BOX"][:2]
    localq=(full-ester[0])@frame
    pr=np.array([g["RAD"][e] for e in g["elems"]])
    margin=float((cdist(cq.reshape(12,3),localq)-RR[i,24:,None]-pr[None,:]+.4).min())
    tmargin=float((cdist(cq.reshape(12,3),M[i,:24])-RR[i,24:,None]-RR[i,None,:24]+.4).min())
    gm=bool(geometry(cq[None],ch[None],oxygen[j])[0]);pp=bool(pairclear(cq[None])[0])
    exact=gm and pp and margin>=-1e-8
    pred=gm and pp and bool(g["clear"](cq.reshape(1,12,3),RR[i,24:],envs[j])[0]>=-1e-8)
    assert exact==pred
    exact_joint=exact and tmargin>=-1e-8 and triad_ok[i,j]
    pred_joint=pred and bool(g["clear"](cq.reshape(1,12,3),RR[i,24:],tenv)[0]>=-1e-8) and triad_ok[i,j]
    assert exact_joint==pred_joint
    audit.append(dict(candidate=meta[i]["id"],trial=shift,site=rows[j]["site_id"],donor_pass=bool(exact),joint_pass=bool(exact_joint),PET_margin_A=margin))
 if i%10==0:print("CANDIDATE",i+1,"of",args.limit,"seconds",round(time.monotonic()-t0,1),flush=True)
cohorts={"support90":support>=.9,"support80":support>=.8,"support70_all":support>=.7}
rankings={}
for scope,name in enumerate(config["reported_scopes"]):
 rankings[name]={}
 for cohort,mask in cohorts.items():
  rr=[]
  for i in range(args.limit):
   fr=counts[scope,i][:,mask]
   perrep=fr.sum(axis=1)/n
   p=fr.sum(axis=0)/(2*n)
   pref=prefix_counts[scope,i][:,mask].sum()/(n)
   rr.append(dict(id=meta[i]["id"],index=i,triad_source=meta[i]["template"],center_sites=int(nominal[scope,i,mask].sum()),total_sites=int(mask.sum()),mean_covered_sites=float(p.sum()),replicate_scores=perrep.tolist(),half_prefix_score=float(pref),sites_at_least_10pct=int((p>=.1).sum()),sites_at_least_50pct=int((p>=.5).sum()),sites_any_sample=int((p>0).sum()),center_N_A=M[i,[26,32]].tolist(),NH_unit_vectors=((H[i]-M[i,[26,32]])/1.01).tolist()))
  rr.sort(key=lambda d:(-d["mean_covered_sites"],-d["center_sites"],d["id"]))
  for k,r in enumerate(rr):r["robust_rank"]=k+1
  center_order=sorted(rr,key=lambda d:(-d["center_sites"],-d["mean_covered_sites"],d["id"]))
  rankings[name][cohort]={"robust":rr,"center":center_order}
np.savez_compressed(root/"robustness.npz",counts=counts,prefix_counts=prefix_counts,nominal=nominal,sample_coverage=sample_coverage,offsets_A=offsets,triad_ok=triad_ok,site_ids=z["site_ids"],support=support)
(root/"rankings.json").write_text(json.dumps(rankings))
(root/"best_sampled_offsets.json").write_text(json.dumps(best_samples))
summary={"status":"TECHNICAL_COMPLETE_GEOMETRY_ONLY","candidate_centers":args.limit,"target_sites":578,"samples_per_center":2*n,"radius_A":1,"zero_shift_replays_passed":args.limit,"independent_full_PET_audit":audit,"elapsed_s":time.monotonic()-t0,"top5":{scope:{cohort:{"robust":v["robust"][:5],"center":v["center"][:5]} for cohort,v in cohorts_.items()} for scope,cohorts_ in rankings.items()},"interpretation":config}
(root/"summary.json").write_text(json.dumps(summary,indent=2))
for p,hx in hashes.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==hx
(root/"SHA256.json").write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir() if p.is_file() and p.name!="RUN_LOG.jsonl"},indent=2))
log("COMPLETE",elapsed_s=time.monotonic()-t0,zero_shift_replays=args.limit,input_unchanged=True)
print("COMPLETE",json.dumps({"output":str(root),"elapsed_s":time.monotonic()-t0,"top":{k:v["support70_all"]["robust"][:3] for k,v in rankings.items()}},default=float),flush=True)
