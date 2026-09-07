"""Extend frozen single-NH protocol to support >=70%; compute only new sites."""
from pathlib import Path
import json,hashlib,shutil,datetime
import numpy as np
from scipy.spatial.distance import cdist
P=Path("/data/bht2/polymer_material_reference_20260827/simulation_slabs/pet_dp10_400chain_direct_v1")
base_source=P/"pet_donor_geometry_batch_singleNH_v1.py"
source_text=base_source.read_text()
assert hashlib.sha256(base_source.read_bytes()).hexdigest()=="b796be0506a0d696ee19aad9288e2b71b6c81131c30f7e4bf2deb00c94ddb80a"
exec(compile(source_text[:source_text.index("rows=collect()")],str(base_source),"exec"))
parent=OUT/"production"
old_config=json.loads((parent/"config.json").read_text())
assert old_config["input_sha256"]==hashlib.sha256(rawpath.read_bytes()).hexdigest()
manifest=json.loads((parent/"SHA256.json").read_text())
for f,h in manifest.items():assert hashlib.sha256((parent/f).read_bytes()).hexdigest()==h
# Only selection threshold changes. All sampling/evaluation functions are untouched.
start=source_text.index("def collect():");end=source_text.index("\ndef sample_donor(",start)
collect_source=source_text[start:end]
assert collect_source.count("support>=.8")==1
exec(compile(collect_source.replace("support>=.8","support>=.7"),"<collect_threshold70>","exec"))
rows70=collect();oldrows=json.loads((parent/"eligible_sites.json").read_text())
oldids={r["site_id"] for r in oldrows};ids70={r["site_id"] for r in rows70}
assert len(oldids)==470 and len(ids70)==578 and oldids<=ids70
newrows=[r for r in rows70 if r["site_id"] not in oldids]
assert len(newrows)==108 and all(.7<=r["support_fraction"]<.8 for r in newrows)
root=OUT/"extension70_v1";root.mkdir(exist_ok=True);(root/"cases").mkdir(exist_ok=True)
cfg={"parent_config":old_config,"support_threshold":.7,"extension_source_sha256":hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),"new_sites":108,"reused_sites":470}
if (root/"config.json").exists():assert json.loads((root/"config.json").read_text())==cfg
else:(root/"config.json").write_text(json.dumps(cfg,indent=2))
(root/"eligible_sites.json").write_text(json.dumps(rows70,indent=2))
(root/"reused_results.json").write_text(json.dumps({"parent":str(parent),"site_ids":sorted(oldids),"angular_histograms":"For reused sites, read the parent cases directory; unchanged and not duplicated."},indent=2))
for sid in oldids:
 src=parent/"cases"/(sid+".json");dst=root/"cases"/(sid+".json")
 if dst.exists():assert dst.read_bytes()==src.read_bytes()
 else:shutil.copy2(src,dst)
(root/"RUNBOOK.md").write_text("Run existing CPU environment python -B this extension script. Original source is SHA-pinned and loaded only up to definitions, without executing its batch driver. Only collect support>=.8 changes to >=.7. Previous 470 results are verified and reused byte-for-byte; only 108 new sites are evaluated. All geometry windows, random seeds, single NH fragment, 2048/4096/8192 Sobol schedule and 0.4 A overlap allowance remain unchanged. No per-site PDBs. Per-site histograms for reused sites remain in production/cases; new histograms are here in cases. Results are finite-sample geometric candidates, not catalyst activity or an impossibility proof. Density support is a heuristic, not a validated dynamic tail classification.\n")
print("START_EXTENSION",len(rows70),len(newrows),flush=True)
for i,row in enumerate(newrows):
 evaluate(row,root)
 if (i+1)%10==0 or i+1==len(newrows):print("PROGRESS_70_NEW",i+1,len(newrows),flush=True)
summary=summarize(rows70,root,"support70")
newsummary=summarize(newrows,root,"new70_only")
# Recompute status categories and validate all per-pose denominators.
cases=[json.loads((root/"cases"/(r["site_id"]+".json")).read_text()) for r in rows70]
def counts(group):
 tr=sum(any(t["triad"]["polymer_clear"]>0 for t in r["templates"]) for r in group)
 pair=sum(any(t["donor"]["poses_with_pair"]>0 for t in r["templates"]) for r in group)
 return {"sites":len(group),"found_pair":pair,"triad_but_no_pair":tr-pair,"no_sampled_triad":len(group)-tr}
for r in cases:
 for t in r["templates"]:
  d=t["donor"]
  assert len(d["unordered_pairs_per_pose"])==t["triad"]["poses_evaluated"]
  assert sum(d["unordered_pairs_per_pose"])==d["pose_pair_combinations"]
  assert sum(v>0 for v in d["unordered_pairs_per_pose"])==d["poses_with_pair"]
for sid in oldids:assert (root/"cases"/(sid+".json")).read_bytes()==(parent/"cases"/(sid+".json")).read_bytes()
for f,h in manifest.items():assert hashlib.sha256((parent/f).read_bytes()).hexdigest()==h
assert not list(root.rglob("*.pdb"))
newcases=[r for r in cases if r["site_id"] not in oldids]
# Independently replay one new successful pair against ALL polymer heavy atoms.
chosen=next((r for r in newcases if any(t["example_audit"] for t in r["templates"])),None)
replay=None
if chosen:
 row=next(r for r in newrows if r["site_id"]==chosen["site_id"])
 ester,frame,ids=site_frame(row);center=ester[0]
 fullq=xyz.copy();fullq[:,:2]-=np.round((fullq[:,:2]-center[:2])/BOX[:2])*BOX[:2]
 local=np.linalg.norm(fullq-center,axis=1)<=30;env=trees(fullq[local],elems[local])
 q,H,angles,directions,hdir,roll=sample_donor(ester,frame,seed_for(row["site_id"]+"|singleNH|Sobol"))
 t=next(t for t in chosen["templates"] if t["example_audit"]);name=t["template"];ex=t["example_audit"]
 rng=np.random.default_rng(seed_for(row["site_id"]+"|"+name+"|triad"))
 poses,stats=make_triads(ester,frame,templates[name],env);assert stats==t["triad"]
 pose=poses[ex["pose"]];a=q[ex["sample_i"]];b=q[ex["sample_j"]]
 pr=np.array([RAD[e] for e in elems]);rr=np.array([RAD[e] for e in donor_templates[0][2]]);tr=np.array([RAD[e] for e in templates[name][1]])
 def margin(x,xr,y,yr):return float((cdist(x,y)-xr[:,None]-yr[None,:]+.4).min())
 margins={"triad_PET":margin(pose,tr,fullq,pr),"A_PET":margin(a,rr,fullq,pr),"B_PET":margin(b,rr,fullq,pr),"A_triad":margin(a,rr,pose,tr),"B_triad":margin(b,rr,pose,tr),"A_B":margin(a,rr,b,rr)}
 assert min(margins.values())>=-1e-8
 replay={"site_id":row["site_id"],"template":name,"margins_with_0.4A_allowance":margins}
audit={"status":"AUDIT_PASS","support70":counts(cases),"new70_only":counts(newcases),"reused_unchanged":470,"parent_hashes_unchanged":len(manifest),"independent_full_PET_replay":replay,"per_site_PDB_count":0}
(root/"independent_audit.json").write_text(json.dumps(audit,indent=2))
with (root/"RUN_LOG.jsonl").open("a") as f:f.write(json.dumps({"time":datetime.datetime.now(datetime.timezone.utc).isoformat(),"status":"COMPLETE","audit":audit})+"\n")
hashes={str(f.relative_to(root)):hashlib.sha256(f.read_bytes()).hexdigest() for f in root.rglob("*") if f.is_file() and f.name!="SHA256.json"}
(root/"SHA256.json").write_text(json.dumps(hashes,indent=2))
print("EXTENSION70_COMPLETE",json.dumps(audit),flush=True)
