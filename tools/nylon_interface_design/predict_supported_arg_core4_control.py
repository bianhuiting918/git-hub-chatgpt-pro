"""Independent prediction of exploratory supported-Arg sequences using existing matched settings."""
from pathlib import Path
import subprocess,json,datetime,hashlib,shutil
r=Path(__file__).resolve().parent;outroot=r.parent/'n9_core_donor_design_20260910_gpu_outputs_v1'
src=outroot/'supported_ARG_SC4559268_mpnn20_core4_v1/seqs/supported_ARG_SC4559268_mpnn20_core4_v1.fa'
records=[]
for line in src.read_text().splitlines():
 if line.startswith('>'):records.append([line,''])
 elif line.strip():records[-1][1]+=line.strip()
seqs=[x[1] for x in records[1:]]
assert len(seqs)==len(set(seqs))==20
fixed=[1,40,42,63]
assert all(len(s)==126 and ''.join(s[i-1] for i in fixed)=='TDDR' for s in seqs)
queries={f'R6_SUP_CORE4_{i:02d}':s for i,s in enumerate(seqs,1)}
g=r/'route6_supported_arg_core4_af20_v1';g.mkdir(exist_ok=False)
fa=g/'queries.fasta';fa.write_text(''.join('>'+n+'\n'+s+'\n' for n,s in queries.items()))
(g/'sequence_mapping.json').write_text(json.dumps(queries,indent=2))
c=json.loads((r/'route6_lowcov_af40_v1/config.json').read_text());msa=c['msa'].copy();af=c['af2'].copy()
msadir=outroot/'route6_supported_arg_core4_msa20_v1';afdir=outroot/'af2_route6_supported_arg_core4_20_v1';inp=g/'msa_inputs';inp.mkdir()
for p in [msadir,afdir]:p.mkdir(mode=0o2770,exist_ok=False);p.chmod(0o2770)
msa[8]=str(fa);msa[9]=str(msadir);af[8]=str(inp);af[9]=str(afdir)
assert msa[0]==af[0]=='/usr/local/bin/gpurun'
(g/'config.json').write_text(json.dumps(dict(msa=msa,af2=af,source=str(src),source_sha256=hashlib.sha256(src.read_bytes()).hexdigest(),status='EXPLORATORY_NOT_GEOMETRY_PASS'),indent=2))
(g/'RUNBOOK.md').write_text('Run this saved script once with the project CPU Python. It uses existing matched MSA and AF2 commands via gpurun, writes versioned directories, and stops on any failure without retry. MSA query identity must match. Source backbone has known defects and this batch is exploratory. Prediction completion is not catalytic acceptance. Do not overwrite outputs.\n')
def log(stage,**kw):
 rec=dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),step='supported_arg_core4_af20_'+stage,**kw)
 for p in [g/'events.jsonl',r/'RUN_LOG.jsonl']:
  with p.open('a') as h:h.write(json.dumps(rec)+'\n')
def run(stage,args):
 log(stage+'_start',command=args)
 with (g/(stage+'.stdout')).open('x') as h:rc=subprocess.call(args,cwd=r,stdout=h,stderr=subprocess.STDOUT)
 log(stage+'_exit',exit_code=rc)
 if rc:raise SystemExit(rc)
run('msa',msa)
for n,seq in queries.items():
 f=msadir/(n+'.a3m');q='';started=False
 for line in f.read_text().splitlines():
  if line.startswith('#'):continue
  if line.startswith('>'):
   if started:break
   started=True
  elif started:q+=line.strip()
 assert q==seq,('MSA_QUERY_MISMATCH',n)
 shutil.copy2(f,inp/f.name)
run('af2',af)
for n in queries:assert (afdir/(n+'.done.txt')).exists() and len(list(afdir.glob(n+'_unrelaxed*.pdb')))==1
log('technical_complete',structures=20,status='PREDICTIONS_COMPLETE_NOT_GEOMETRY_PASS')
