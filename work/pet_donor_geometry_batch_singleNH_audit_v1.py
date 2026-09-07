from pathlib import Path
import json,hashlib,csv,datetime
import numpy as np
from scipy.spatial.distance import cdist
P=Path("/data/bht2/polymer_material_reference_20260827/simulation_slabs/pet_dp10_400chain_direct_v1")
source=P/"pet_donor_geometry_batch_singleNH_v1.py"
text=source.read_text()
exec(compile(text[:text.index("rows=collect()")],str(source),"exec"))
root=OUT/"production"
assert (root/"support80_summary.json").exists()
manifest=json.loads((root/"SHA256.json").read_text())
for name,sha in manifest.items():assert hashlib.sha256((root/name).read_bytes()).hexdigest()==sha,name
rows=json.loads((root/"eligible_sites.json").read_text());byid={r["site_id"]:r for r in rows}
cases={p.stem:json.loads(p.read_text()) for p in (root/"cases").glob("*.json")}
assert len(cases)==len(rows)==470
counts={}
for tag,threshold_support in [("support90",.9),("support80",.8)]:
 group=[r for r in cases.values() if r["support_fraction"]>=threshold_support]
 summary=json.loads((root/(tag+"_summary.json")).read_text())
 found=sum(any(t["donor"]["poses_with_pair"]>0 for t in r["templates"]) for r in group)
 triad=sum(any(t["triad"]["polymer_clear"]>0 for t in r["templates"]) for r in group)
 assert summary["sites"]==len(group) and summary["sites_with_any_template_pair"]==found and summary["sites_with_any_sampled_triad"]==triad
 table=list(csv.DictReader((root/(tag+"_sites.tsv")).open(),delimiter="\t"))
 assert len(table)==2*len(group) and len({(x["site_id"],x["template"]) for x in table})==len(table)
 counts[tag]={"sites":len(group),"found_pair":found,"sampled_triad_but_no_pair":triad-found,"no_sampled_triad":len(group)-triad,"top":sum(r["side"]=="top" for r in group),"bottom":sum(r["side"]=="bottom" for r in group)}
 for r in group:
  for t in r["templates"]:
   f=t["donor"]
   assert len(f["unordered_pairs_per_pose"])==t["triad"]["poses_evaluated"]
   assert sum(x>0 for x in f["unordered_pairs_per_pose"])==f["poses_with_pair"]
   assert sum(f["unordered_pairs_per_pose"])==f["pose_pair_combinations"]
# Smoke reproducibility: identical seed/config gives identical numerical outcomes.
for f in (OUT/"smoke/cases").glob("*.json"):
 a=json.loads(f.read_text());b=cases[a["site_id"]]
 a.pop("elapsed_s");b=dict(b);b.pop("elapsed_s");assert a==b
# Independent distance replay on one success from each cohort.
replays=[]
for lower,upper in [(.9,1.01),(.8,.9)]:
 r=next(r for r in cases.values() if lower<=r["support_fraction"]<upper and any(t["example_audit"] for t in r["templates"]))
 row=byid[r["site_id"]];ester,frame,ids=site_frame(row);center=ester[0]
 fullq=xyz.copy();fullq[:,:2]-=np.round((fullq[:,:2]-center[:2])/BOX[:2])*BOX[:2]
 local=np.linalg.norm(fullq-center,axis=1)<=30;env=trees(fullq[local],elems[local])
 q,H,angles,directions,hdir,roll=sample_donor(ester,frame,seed_for(r["site_id"]+"|singleNH|Sobol"))
 t=next(t for t in r["templates"] if t["example_audit"]);name=t["template"];ex=t["example_audit"]
 rng=np.random.default_rng(seed_for(r["site_id"]+"|"+name+"|triad"))
 poses,stats=make_triads(ester,frame,templates[name],env);assert stats==t["triad"]
 pose=poses[ex["pose"]];frags=[q[ex["sample_i"]],q[ex["sample_j"]]]
 pe=np.array([RAD[e] for e in elems]);rr=np.array([RAD[e] for e in donor_templates[0][2]])
 tr=np.array([RAD[e] for e in templates[name][1]])
 margins={"triad_PET":float((cdist(pose,fullq)-tr[:,None]-pe[None,:]+.4).min())}
 for j,frag in enumerate(frags):
  margins["donor"+str(j)+"_PET"]=float((cdist(frag,fullq)-rr[:,None]-pe[None,:]+.4).min())
  margins["donor"+str(j)+"_triad"]=float((cdist(frag,pose)-rr[:,None]-tr[None,:]+.4).min())
 margins["donor_donor"]=float((cdist(frags[0],frags[1])-rr[:,None]-rr[None,:]+.4).min())
 assert min(margins.values())>=-1e-8
 replays.append({"site_id":r["site_id"],"template":name,"margins_with_0.4A_allowance":margins})
assert not list(root.rglob("*.pdb"))
result={"status":"AUDIT_PASS","hashes_checked":len(manifest),"cohorts":counts,"new_at80":{k:counts["support80"][k]-counts["support90"][k] for k in ["sites","found_pair","sampled_triad_but_no_pair","no_sampled_triad"]},"smoke_reproduction":"PASS","independent_full_PET_distance_replays":replays,"per_site_PDB_count":0}
(root/"independent_audit.json").write_text(json.dumps(result,indent=2))
print(json.dumps(result),flush=True)
