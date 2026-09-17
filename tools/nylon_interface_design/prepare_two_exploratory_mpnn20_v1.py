"""Prepare exploratory MPNN inputs only. Does not launch GPU."""
from pathlib import Path
import json,hashlib,datetime
r=Path(__file__).resolve().parent
b=r/'route56_N10_20260915_v1/route6_top50_local_generation_v1'
out=b/'two_repaired_exploratory_mpnn20_v1'
assert not out.exists()
fixed=[1,40,42,63]
records=[];jobs=[]
for cid in ['SC4916877','SC4242958']:
 src=b/(cid+'_seed2_multisite_repair_v3_rebuildH')/'after.pdb'
 report=json.loads(sorted(src.parent.glob('idealH_independent_audit_*.json'))[-1].read_text())
 assert hashlib.sha256(src.read_bytes()).hexdigest()==report['source_sha256']
 assert report['CN_bad']==0 and report['omega_gt30_count']==0
 a={};rn={}
 for line in src.read_text().splitlines():
  if line.startswith('ATOM'):
   assert line[21]=='A'
   i=int(line[22:26]);rn[i]=line[17:20];a.setdefault(i,{})[line[12:16].strip()]=[float(line[k:k+8]) for k in [30,38,46]]
 assert sorted(a)==list(range(1,127))
 codes={'GLY':'G','THR':'T','ASP':'D','ARG':'R'}
 seq=''.join(codes[rn[i]] for i in range(1,127))
 assert ''.join(seq[i-1] for i in fixed)=='TDDR'
 name=cid+'_seed2_repaired_exploratory'
 record=dict(name=name,num_of_chains=1,seq=seq,seq_chain_A=seq,coords_chain_A={n+'_chain_A':[a[i][n] for i in range(1,127)] for n in ['N','CA','C','O']})
 records.append(record)
 dest=r.parent/'n9_core_donor_design_20260910_gpu_outputs_v1/two_repaired_exploratory_mpnn20_v1'/name
 inp=out/name
 cmd=['gpurun','env','CUDA_VISIBLE_DEVICES=1','PATH=/data/bht2/SS-RFD3/.venv/bin:/usr/local/bin:/usr/bin:/bin','python',str(r/'software/ProteinMPNN_api_pinned/protein_mpnn_run.py'),'--path_to_model_weights',str(r/'software/ProteinMPNN_api_pinned/vanilla_model_weights'),'--model_name','v_48_020','--jsonl_path',str(inp/'parsed.jsonl'),'--chain_id_jsonl',str(inp/'assigned.jsonl'),'--fixed_positions_jsonl',str(inp/'fixed.jsonl'),'--out_folder',str(dest),'--num_seq_per_target','20','--batch_size','1','--sampling_temp','0.1','--seed','101','--backbone_noise','0.0']
 jobs.append(dict(name=name,input=str(inp),output=str(dest),source=str(src),source_sha256=report['source_sha256'],command=cmd,known_generic_radius_contacts=report['self_bad_pairs'],status='EXPLORATORY_READY_NOT_SUBMITTED'))
out.mkdir()
for record,j in zip(records,jobs):
 inp=Path(j['input']);inp.mkdir()
 for file,obj in [('parsed.jsonl',record),('assigned.jsonl',{j['name']:[['A'],[]]}),('fixed.jsonl',{j['name']:{'A':fixed}})]:
  with (inp/file).open('x') as f:f.write(json.dumps(obj)+'\n')
with (out/'manifest.json').open('x') as f:json.dump(dict(status='EXPLORATORY_READY_NOT_SUBMITTED',fixed_positions=fixed,planned_sequences_per_target=20,jobs=jobs,limitations=['Internal contacts not resolved','Free Cartesian hydrogen optimization invalidates donor-H inference','Idealized-H geometry only; not stable catalytic core evidence','Independent folding, mature geometry, material replay and complete sequence validation required']),f,indent=2)
(out/'RUNBOOK.md').write_text('Preparation only; no GPU launched. Use manifest commands through gpurun after an existing task slot is available. Create each output directory with project group access, never modify GPU infrastructure. Preserve these exploratory labels. Validate 20 unique sequences per target, lengths and fixed TDDR identities; then independent structure prediction and full geometry/material audit. No automatic retries on GPU failure. Stop on nonzero command exit. No existing files overwritten.\n')
with (r/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),step='two_exploratory_mpnn_inputs_prepared_NOT_SUBMITTED',output=str(out),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),targets=2,planned_sequences=40))+'\n')
print(json.dumps(dict(output=str(out),targets=2,planned_sequences=40,status='EXPLORATORY_READY_NOT_SUBMITTED')))
