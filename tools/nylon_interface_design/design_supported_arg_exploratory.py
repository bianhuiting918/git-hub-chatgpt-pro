"""Exploratory ProteinMPNN designs; source backbone is NOT a geometry PASS."""
from pathlib import Path
import subprocess,json,datetime,hashlib,os
r=Path(__file__).resolve().parent
name='supported_ARG_SC4559268_mpnn20_v1'
inp=r/'inputs'/name;out=r.parent/'n9_core_donor_design_20260910_gpu_outputs_v1'/name
src=r/'route56_N10_20260915_v1/route6_supported_ARG_SC4559268_v1/fixed_core_repair_v1/after.pdb'
assert hashlib.sha256(src.read_bytes()).hexdigest()=='bd6c4d376b4a950a3f1a22457dbb54dac40d162823d6015b5a2baad8d0e99ab54'
a={};rn={}
for l in src.read_text().splitlines():
 if l.startswith('ATOM  '):
  assert l[21]=='A'
  i=int(l[22:26]);rn[i]=l[17:20];a.setdefault(i,{})[l[12:16].strip()]=[float(l[30:38]),float(l[38:46]),float(l[46:54])]
assert sorted(a)==list(range(1,127))
codes={'GLY':'G','THR':'T','ASP':'D','ARG':'R'};seq=''.join(codes[rn[i]] for i in range(1,127))
fixed=[1,40,41,42,61,62,63,64,65]
assert ''.join(seq[i-1] for i in fixed)=='TDGDGGRGG'
inp.mkdir(exist_ok=False);out.mkdir(mode=0o2770,exist_ok=False);out.chmod(0o2770)
record=dict(name=name,num_of_chains=1,seq=seq,seq_chain_A=seq,coords_chain_A={n+'_chain_A':[a[i][n] for i in range(1,127)] for n in ['N','CA','C','O']})
for file,obj in [('parsed.jsonl',record),('assigned.jsonl',{name:[['A'],[]]}),('fixed.jsonl',{name:{'A':fixed}})]:
 with (inp/file).open('x') as h:h.write(json.dumps(obj)+'\n')
cmd=['gpurun','env','CUDA_VISIBLE_DEVICES=1','PATH=/data/bht2/SS-RFD3/.venv/bin:/usr/local/bin:/usr/bin:/bin','python',str(r/'software/ProteinMPNN_api_pinned/protein_mpnn_run.py'),'--path_to_model_weights',str(r/'software/ProteinMPNN_api_pinned/vanilla_model_weights'),'--model_name','v_48_020','--jsonl_path',str(inp/'parsed.jsonl'),'--chain_id_jsonl',str(inp/'assigned.jsonl'),'--fixed_positions_jsonl',str(inp/'fixed.jsonl'),'--out_folder',str(out),'--num_seq_per_target','20','--batch_size','1','--sampling_temp','0.1','--seed','101','--backbone_noise','0.0']
manifest={'status':'EXPLORATORY_NOT_GEOMETRY_PASS','source':str(src),'source_sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'fixed':fixed,'command':cmd,'known_defects':{'raw_self_close_pairs':9,'fixed_core_pairs':2,'omega_gt30':1},'output':str(out)}
with (inp/'manifest.json').open('x') as h:json.dump(manifest,h,indent=2)
with (inp/'RUNBOOK.md').open('x') as h:h.write('Exploratory fixed-core design, not a accepted backbone. Run this saved script once from the project root using the existing CPU Python; it calls ProteinMPNN via gpurun. Twenty sequences, seed101, temperature0.1, backbone_noise0. No automatic retry. Independent folding and joint catalytic/material audit required. Preserve baseline and defects.\n')
with (r/'RUN_LOG.jsonl').open('a') as h:h.write(json.dumps({'time':datetime.datetime.now(datetime.timezone.utc).isoformat(),'step':'supported_ARG_mpnn_start',**manifest})+'\n')
p=subprocess.run(cmd,cwd=r)
with (r/'RUN_LOG.jsonl').open('a') as h:h.write(json.dumps({'time':datetime.datetime.now(datetime.timezone.utc).isoformat(),'step':'supported_ARG_mpnn_exit','exit_code':p.returncode,'output':str(out)})+'\n')
raise SystemExit(p.returncode)
