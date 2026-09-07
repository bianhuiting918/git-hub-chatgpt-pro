"""Statistics only: existing accepted NH positions. No new pose generation or collision tests."""
from pathlib import Path
import json,hashlib,datetime,numpy as np
from scipy.spatial.distance import cdist
P=Path("/data/bht2/polymer_material_reference_20260827/simulation_slabs/pet_dp10_400chain_direct_v1")
base=P/"donor_geometry_batch_singleNH_v1/joint_transfer_map_v1"
source_files=[base/"joint_transfer_results.npz",base/"motifs.json"]
input_hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in source_files}
z=np.load(base/"joint_transfer_results.npz");meta=json.loads((base/"motifs.json").read_text())
Q=z["motif_xyz_A"][:,[26,32],:];coverage=z["coverage"];C=coverage.sum(axis=1)
assert Q.shape==(473,2,3) and coverage.shape==(473,578) and C.max()==31 and C.min()>0
# Count every existing site-candidate-donor incidence, including coincident locations.
points=Q.reshape(-1,3)
weights=np.repeat(C,2)
dist=cdist(Q.reshape(-1,3),points).reshape(473,2,-1)
near=dist<=1.0+1e-12
D=near@weights
density=D.min(axis=1)
score=.6*C/C.max()+.4*density/density.max()
# Optional audit: are nearby donors also present together in an existing pair?
same=np.linalg.norm(Q[:,None,0]-Q[None,:,0],axis=-1)<=1+1e-12
same&=np.linalg.norm(Q[:,None,1]-Q[None,:,1],axis=-1)<=1+1e-12
swap=np.linalg.norm(Q[:,None,0]-Q[None,:,1],axis=-1)<=1+1e-12
swap&=np.linalg.norm(Q[:,None,1]-Q[None,:,0],axis=-1)<=1+1e-12
pair_near=same|swap
def record(i):
 return {"id":meta[i]["id"],"triad_source":meta[i]["template"],"coverage":int(C[i]),"N1_accepted_point_records_1A":int(D[i,0]),"N2_accepted_point_records_1A":int(D[i,1]),"density_min":int(density[i]),"score":float(score[i]),"N_A":Q[i].tolist(),"existing_joint_neighbor_candidates_including_self":int(pair_near[i].sum())}
order=sorted(range(473),key=lambda i:(-score[i],-C[i],meta[i]["id"]))
ranks={i:k+1 for k,i in enumerate(order)}
records=[dict(rank=ranks[i],**record(i)) for i in order]
# Repeated colors at identical positions count once per accepted site-candidate incidence.
assert int(weights.sum())==6526
assert np.all(D>=C[:,None])
# Swapping donor labels leaves the chosen density and ranking score unchanged.
assert np.array_equal(D[:,::-1].min(axis=1),density)
winner=order[0];brute=[sum(int(w) for p,w in zip(points,weights) if np.linalg.norm(p-n)<=1+1e-12) for n in Q[winner]]
assert brute==D[winner].tolist()
# Exact location cutoff sensitivity check, not another model calculation.
boundary=np.abs(dist-1).min()
out={"status":"EXISTING_DATA_STATISTICS_ONLY","candidate_count":473,"ester_sites":578,"existing_accepted_site_candidate_edges":int(coverage.sum()),"existing_N_records":946,"accepted_N_point_records":int(weights.sum()),"radius_A":1,"coverage_normalizer":int(C.max()),"density_normalizer":int(density.max()),"weights":{"coverage":.6,"existing_point_density":.4},"density_definition":"min of accepted site-candidate-N incidence counts within 1 A of each reference N; includes center; coincident coordinates retain multiplicity; no deduplication","scope":"same donor chemical fragment; full original library pooled; NOT conditioned on an identical triad pose; no energy, perturbation sampling, new geometry tests, or proven robustness","top10":records[:10],"old_coverage_winner":dict(rank=ranks[449],**record(449)),"cutoff_nearest_gap_A":float(boundary)}
root=P/"donor_pair_existing_density_1A_weighted_v2";root.mkdir(exist_ok=False)
out["input_sha256"]=input_hashes
(root/"summary.json").write_text(json.dumps(out,indent=2))
(root/"all473_rankings.json").write_text(json.dumps(records,indent=2))
(root/"RUNBOOK.md").write_text("# Existing accepted PET donor point density, multiplicity-weighted v2\n\nNo new poses, perturbations or geometry validation are calculated. Run the saved pet_donor_existing_density_1A_weighted_v2.py on CPU using simulation_slabs/unbiased_100chain_v1/env/bin/python -B; output folder must not exist.\n\nEvery original accepted site-candidate-donor incidence is counted; same position fitting n esters contributes n. The 473 by 578 coverage matrix has 3263 accepted incidences and 6526 donor records. The reference center is included. Radius 1 A about each NH N; density = smaller of the two record counts. Score = 0.6 C/max(C) + 0.4 D/max(D). No coordinate or ester-site deduplication is used in the density. The exact center coverage counts distinct ester sites.\n\nDensity pools existing accepted motifs across different triad poses; it is not conditional compatibility with an identical triad or simultaneous pair validity. It is a positional-density proxy, not a perturbation success probability, energy, or activity. Candidate enzyme names identify triad source only, not native full motifs. Old deduplicated-density source and independent-perturbation pilot are not used.\n")
for p,h in input_hashes.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
(root/"RUN_LOG.jsonl").write_text(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=__file__,status=out["status"],weights=out["weights"],radius_A=1,input_unchanged=True,ranking_count=len(records)))+"\n")
(root/"SHA256.json").write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir() if p.is_file()},indent=2))
print(json.dumps({"output":str(root),"top":records[:6],"source_inputs_unchanged":True,"candidate_count":len(records)}),flush=True)
