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
# Round only for identity deduplication; distances use the original first coordinate.
_,first,inverse=np.unique(np.round(Q.reshape(-1,3),6),axis=0,return_index=True,return_inverse=True)
points=Q.reshape(-1,3)[first];dist=cdist(Q.reshape(-1,3),points).reshape(473,2,-1)
near=(dist<=1.0+1e-12)&(dist>1e-6)
D=near.sum(axis=2);density=D.min(axis=1);score=.6*C/C.max()+.4*density/density.max()
# Optional audit: are nearby donors also present together in an existing pair?
same=np.linalg.norm(Q[:,None,0]-Q[None,:,0],axis=-1)<=1+1e-12
same&=np.linalg.norm(Q[:,None,1]-Q[None,:,1],axis=-1)<=1+1e-12
swap=np.linalg.norm(Q[:,None,0]-Q[None,:,1],axis=-1)<=1+1e-12
swap&=np.linalg.norm(Q[:,None,1]-Q[None,:,0],axis=-1)<=1+1e-12
pair_near=same|swap
def record(i):
 return {"id":meta[i]["id"],"triad_source":meta[i]["template"],"coverage":int(C[i]),"N1_other_positions_1A":int(D[i,0]),"N2_other_positions_1A":int(D[i,1]),"density_min":int(density[i]),"score":float(score[i]),"N_A":Q[i].tolist(),"existing_joint_neighbor_candidates_including_self":int(pair_near[i].sum())}
order=sorted(range(473),key=lambda i:(-score[i],-C[i],meta[i]["id"]))
ranks={i:k+1 for k,i in enumerate(order)}
records=[dict(rank=ranks[i],**record(i)) for i in order]
# A position accepted at multiple ester sites is one point, not multiple density votes.
assert len(points)<=946
# Swapping donor labels leaves the chosen density and ranking score unchanged.
assert np.array_equal(D[:,::-1].min(axis=1),density)
winner=order[0];brute=[sum(1 for p in points if 1e-6<np.linalg.norm(p-n)<=1+1e-12) for n in Q[winner]]
assert brute==D[winner].tolist()
# Exact location cutoff sensitivity check, not another model calculation.
boundary=np.abs(dist-1).min()
out={"status":"EXISTING_DATA_STATISTICS_ONLY","candidate_count":473,"ester_sites":578,"existing_accepted_site_candidate_edges":int(coverage.sum()),"existing_N_records":946,"unique_N_positions":len(points),"radius_A":1,"coverage_normalizer":int(C.max()),"density_normalizer":int(density.max()),"weights":{"coverage":.6,"existing_point_density":.4},"density_definition":"min of counts of distinct OTHER accepted N positions within 1 A of each fixed donor N; own position excluded; union over the existing accepted library; repeated site colors deduplicated","scope":"same donor chemical fragment; full original library pooled; NOT conditioned on an identical triad pose; no energy, perturbation sampling, new geometry tests, or proven robustness","top10":records[:10],"old_coverage_winner":dict(rank=ranks[449],**record(449)),"all_rankings":records,"cutoff_nearest_gap_A":float(boundary),"winner_neighbor_position_indices":[np.flatnonzero(near[winner,d]).tolist() for d in range(2)],"unique_positions_A":points.tolist()}
root=P/"donor_pair_existing_density_1A_v1";root.mkdir(exist_ok=False)
out["input_sha256"]=input_hashes
(root/"ranking_and_density.json").write_text(json.dumps(out,indent=2))
(root/"RUNBOOK.md").write_text("# Existing PET donor cloud density, no new geometry calculation\n\nRun the saved pet_donor_existing_density_1A_v1.py using the existing CPU Python environment. Output directory must not exist. Parent data remain read-only.\n\n473 fixed motifs, 578 target esters. Count distinct OTHER accepted N positions within 1 A of each reference donor. Deduplicate repeated positions at 1e-6 A; distance tests use original coordinates. Density = smaller of the two counts. Score = 0.6 * center coverage / maximum center coverage + 0.4 * density / maximum density. Neighbor counts pool the existing accepted library across triad poses; they do not establish compatibility with one identical fixed triad pose, nor simultaneous pair feasibility. An optional paired-neighbor count is also retained. No new perturbations, collision tests, energy, or activity inference.\n")
for p,h in input_hashes.items():assert hashlib.sha256(Path(p).read_bytes()).hexdigest()==h
(root/"RUN_LOG.jsonl").write_text(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=__file__,status=out["status"],weights=out["weights"],radius_A=1,input_unchanged=True))+"\n")
(root/"SHA256.json").write_text(json.dumps({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in root.iterdir() if p.is_file()},indent=2))
print(json.dumps({k:v for k,v in out.items() if k not in ["all_rankings","unique_positions_A","winner_neighbor_position_indices"]}),flush=True)
