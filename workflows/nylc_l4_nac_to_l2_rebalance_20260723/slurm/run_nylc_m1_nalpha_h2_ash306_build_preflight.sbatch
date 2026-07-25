#!/bin/bash
#SBATCH -J nylc_m1_build
#SBATCH -p xahcnormal
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH --mem-per-cpu=2500M
#SBATCH -t 01:00:00
#SBATCH -o /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_m1_build_%j.out
#SBATCH -e /work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723/logs/nylc_m1_build_%j.err
set -euo pipefail
TASK=/work/home/acshdt1dks/nylon_pa66_scnet_20260708/l4_nac_to_l2_rebalance_20260723
FLOW="$TASK/repo/workflows/nylc_l4_nac_to_l2_rebalance_20260723"
CAND="$TASK/candidates/nylc_C18_trueT267_freeGS"
SRC="$CAND/runs/mm_proton_geometry_ext19ns_job_61780674"
OLD_BUILD="$CAND/build"
PY=/work/home/acshdt1dks/opt/interface-stability-tools/envs/interface/bin/python
GMX=/public/software/apps/Gromacs-DCU2/2022.1/mpi/bin/gmx_mpi
job="${SLURM_JOB_ID:?}"
OUT="$CAND/microstates/M1_nalphaH2_Asp306H/build_job_${job}"
mkdir -p "$TASK/logs" "$(dirname "$OUT")"
test ! -e "$OUT"
mkdir -p "$OUT/selection" "$OUT/chain_gen" "$OUT/build"
record() {
  local state=$1 detail=$2 stamp
  stamp="$(date -Is)"
  ( flock -x 9; printf '%s\tnylc_m1_nalphaH2_Asp306H_build\t%s\tcandidate=nylc_C18_trueT267_freeGS;job=%s;%s\n' "$stamp" "$state" "$job" "$detail" >&9 ) 9>>"$TASK/run_history.tsv"
  ( flock -x 8; "$PY" - "$stamp" "$state" "$detail" "$job" >&8 <<'PY'
import json,sys
stamp,state,detail,job=sys.argv[1:]
print(json.dumps({"timestamp":stamp,"event":"nylc_m1_nalphaH2_Asp306H_build","state":state,"detail":detail,"candidate_id":"nylc_C18_trueT267_freeGS","slurm_job_id":job},sort_keys=True))
PY
  ) 8>>"$TASK/run_history.jsonl"
}
fail() { local rc=$?; trap - ERR; set +e; record FAIL_TECHNICAL "exit_code=$rc;out=$OUT"; exit "$rc"; }
trap fail ERR
record START "source_job=61780674_cancelled;microstate=NalphaH2_ThrOH_Asp306H_Asp308minus;selection=lowest_potential_fully_unrestrained_NAC;scope=classical_MM_proxy"
test -s "$SRC/run.tpr"
test -s "$SRC/run.part0002.xtc"
test -s "$SRC/run.part0002.edr"
test -s "$OLD_BUILD/topol.top"
sha256sum "$SRC/run.tpr" "$SRC/run.part0002.xtc" "$SRC/run.part0002.edr" > "$OUT/source_sha256.tsv"
set +u
source /work/home/acshdt1dks/opt/gromacs-fastest/env.sh
set -u
export GMX_MAXBACKUP=-1
"$GMX" distance -f "$SRC/run.part0002.xtc" -s "$SRC/run.tpr" \
  -select 'atomnr 8961 plus atomnr 10287' -oall "$OUT/selection/nac_distance.xvg" \
  >"$OUT/selection/distance.stdout" 2>"$OUT/selection/distance.stderr"
printf '[ nac_angle ]\n8961 10287 10288\n' > "$OUT/selection/nac_angle.ndx"
printf '0\n' | "$GMX" angle -f "$SRC/run.part0002.xtc" -n "$OUT/selection/nac_angle.ndx" \
  -ov "$OUT/selection/nac_angle.xvg" -type angle \
  >"$OUT/selection/angle.stdout" 2>"$OUT/selection/angle.stderr"
printf 'Potential\n0\n' | "$GMX" energy -f "$SRC/run.part0002.edr" \
  -o "$OUT/selection/potential_energy.xvg" \
  >"$OUT/selection/energy.stdout" 2>"$OUT/selection/energy.stderr"
"$PY" "$FLOW/scripts/analyze_nac_series.py" \
  --distance-xvg "$OUT/selection/nac_distance.xvg" \
  --angle-xvg "$OUT/selection/nac_angle.xvg" \
  --potential-xvg "$OUT/selection/potential_energy.xvg" \
  --distance-max-nm 0.35 --angle-min-deg 95 --angle-max-deg 115 \
  --output "$OUT/selection/nac_energy_series.json"
selected_time="$($PY - "$OUT/selection/nac_energy_series.json" <<'PY'
import json,sys
x=json.load(open(sys.argv[1]))
r=x.get('lowest_potential_nac_frame')
if r is None: raise SystemExit('no fully unrestrained NAC frame in cancelled M0 trajectory')
print(r['time_ps'])
PY
)"
printf '0\n' | "$GMX" trjconv -f "$SRC/run.part0002.xtc" -s "$SRC/run.tpr" \
  -dump "$selected_time" -o "$OUT/selection/source_selected_NAC.gro" \
  >"$OUT/selection/trjconv.stdout" 2>"$OUT/selection/trjconv.stderr"
"$PY" - "$FLOW" "$OUT/selection/source_selected_NAC.gro" "$OUT/chain_gen/chainH_ASH306.pdb" "$OUT/chain_gen/export_audit.json" <<'PY'
import json,pathlib,sys
flow,source,pdb,audit=map(pathlib.Path,sys.argv[1:])
sys.path.insert(0,str(flow/'scripts'))
from prepare_nylc_nalpha_h2_ash306 import export_chain_pdb
text,result=export_chain_pdb(source.read_text(),8949,10272,ash_resid=306)
assert result['chain_atom_count']==1324 and result['renamed_ash306_atom_count']==12
pdb.write_text(text); audit.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
PY
(
  cd "$OUT/chain_gen"
  "$GMX" pdb2gmx -f chainH_ASH306.pdb -o chainH_ASH306.gro -p chainH.top -i posre_chainH.itp \
    -ff amber99sb-ildn -water tip3p -ignh >pdb2gmx.stdout 2>pdb2gmx.stderr
)
cp "$OLD_BUILD/topol.top" "$OLD_BUILD/PA66_L2_GMX.itp" "$OLD_BUILD/source_cycle.ndx" "$OUT/build/"
cp "$OLD_BUILD"/topol_Protein_chain_A.itp "$OLD_BUILD"/topol_Protein_chain_D.itp "$OLD_BUILD"/topol_Protein_chain_E.itp "$OUT/build/"
cp "$OLD_BUILD"/posre_Protein_chain_A*.itp "$OLD_BUILD"/posre_Protein_chain_D*.itp "$OLD_BUILD"/posre_Protein_chain_E*.itp "$OUT/build/"
cp "$OLD_BUILD"/posre_l2_*.itp "$OUT/build/"
"$PY" - "$FLOW" "$OLD_BUILD/topol_Protein_chain_H.itp" "$OUT/chain_gen/chainH.top" "$OUT/chain_gen/chainH_ASH306.gro" "$OUT/selection/source_selected_NAC.gro" "$OUT/build" "$OUT/build_audit.json" "$selected_time" <<'PY'
import hashlib,json,math,pathlib,re,sys
flow,old_itp,standalone,chain_gro,full_gro,build,audit_path=map(pathlib.Path,sys.argv[1:8]); selected_time=float(sys.argv[8])
sys.path.insert(0,str(flow/'scripts'))
from prepare_nylc_nalpha_h2_ash306 import drop_gro_atom,extract_molecule_itp,replace_gro_slice,rewrite_chain_itp

def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def parse_itp(text):
    section=''; out=[]
    for line in text.splitlines():
        m=re.match(r'\s*\[\s*([^]]+?)\s*\]',line)
        if m: section=m.group(1).strip().lower(); continue
        if section=='atoms' and line.strip() and not line.lstrip().startswith((';','#')):
            t=line.split(';',1)[0].split()
            if len(t)>=8 and t[0].isdigit(): out.append({'nr':int(t[0]),'type':t[1],'resid':int(t[2]),'resname':t[3],'name':t[4],'charge':float(t[6])})
    return out
def gro_atoms(text):
    lines=text.splitlines(); n=int(lines[1]); assert len(lines)==n+3
    return [{'nr':i+1,'resid':int(x[:5]),'resname':x[5:10].strip(),'name':x[10:15].strip(),'xyz':tuple(float(x[a:b]) for a,b in ((20,28),(28,36),(36,44)))} for i,x in enumerate(lines[2:2+n])]
def heavy(name):
    letters=''.join(x for x in name if x.isalpha()).upper(); return bool(letters) and letters[0]!='H'
rewritten,a_itp=rewrite_chain_itp(standalone.read_text())
itp=extract_molecule_itp(rewritten).replace('posre_chainH.itp','posre_Protein_chain_H.itp')
if 'POSRES_WEAK' not in itp:
    itp += '\n#ifdef POSRES_WEAK\n#include "posre_Protein_chain_H_weak.itp"\n#endif\n'
chain_text,a_gro=drop_gro_atom(chain_gro.read_text())
full_text,a_merge=replace_gro_slice(full_gro.read_text(),chain_text,8949,10272)
new_atoms=parse_itp(itp); old_atoms=parse_itp(old_itp.read_text()); chain_atoms=gro_atoms(chain_text)
assert len(new_atoms)==len(chain_atoms)==1324
assert [(x['nr'],x['resid'],x['resname'],x['name']) for x in new_atoms]==[(x['nr'],x['resid'],x['resname'],x['name']) for x in chain_atoms]
old_charge=sum(x['charge'] for x in old_atoms); new_charge=sum(x['charge'] for x in new_atoms)
assert abs(old_charge-new_charge)<1e-5 and abs(new_charge+4.0)<1e-5
old_other={(x['resid'],x['name']):(x['type'],x['charge']) for x in old_atoms if x['resid'] not in (267,306)}
new_other={(x['resid'],x['name']):(x['type'],x['charge']) for x in new_atoms if x['resid'] not in (267,306)}
assert old_other==new_other
thr=[x for x in new_atoms if x['resid']==267]; ash=[x for x in new_atoms if x['resid']==306]; asp308=[x for x in new_atoms if x['resid']==308]
assert {x['name'] for x in thr}.issuperset({'N','H1','H2','OG1','HG1'}) and 'H3' not in {x['name'] for x in thr}
assert {x['name'] for x in ash}.issuperset({'OD1','OD2','HD2'}) and {x['resname'] for x in ash}=={'ASH'}
assert {x['name'] for x in asp308}.issuperset({'OD1','OD2'}) and 'HD2' not in {x['name'] for x in asp308}
source_atoms=gro_atoms(full_gro.read_text()); merged_atoms=gro_atoms(full_text)
assert len(source_atoms)==len(merged_atoms)
for a,b in zip(source_atoms[:8948],merged_atoms[:8948]): assert (a['resid'],a['resname'],a['name'],a['xyz'])==(b['resid'],b['resname'],b['name'],b['xyz'])
for a,b in zip(source_atoms[10272:],merged_atoms[10272:]): assert (a['resid'],a['resname'],a['name'],a['xyz'])==(b['resid'],b['resname'],b['name'],b['xyz'])
source_key={(x['resid'],x['name']):x['xyz'] for x in source_atoms[8948:10272] if heavy(x['name'])}
new_key={(x['resid'],x['name']):x['xyz'] for x in chain_atoms if heavy(x['name'])}
ordinary=set(source_key)&set(new_key)-{(355,'OC1'),(355,'OC2')}
assert ordinary and all(max(abs(a-b) for a,b in zip(source_key[k],new_key[k]))<=0.0011 for k in ordinary)
full_lookup={(x['resid'],x['resname'],x['name']):x for x in merged_atoms}
def local_global(resid,name):
    x=next(x for x in new_atoms if x['resid']==resid and x['name']==name); return 8948+x['nr']
mapping={'thr267_N':local_global(267,'N'),'thr267_H1':local_global(267,'H1'),'thr267_H2':local_global(267,'H2'),'thr267_OG1':local_global(267,'OG1'),'thr267_HG1':local_global(267,'HG1'),'asp306_OD1':local_global(306,'OD1'),'asp306_OD2':local_global(306,'OD2'),'asp306_HD2':local_global(306,'HD2'),'asp308_OD1':local_global(308,'OD1'),'asp308_OD2':local_global(308,'OD2'),'l2_C':10287,'l2_O':10288,'l2_N':10289}
assert merged_atoms[10286]['name']=='C12' and merged_atoms[10287]['name']=='O2' and merged_atoms[10288]['name']=='N3'
def dist(i,j): return math.dist(merged_atoms[i-1]['xyz'],merged_atoms[j-1]['xyz'])
def angle(i,j,k):
    a=merged_atoms[i-1]['xyz']; b=merged_atoms[j-1]['xyz']; c=merged_atoms[k-1]['xyz']; u=[a[q]-b[q] for q in range(3)]; v=[c[q]-b[q] for q in range(3)]; co=sum(u[q]*v[q] for q in range(3))/(math.sqrt(sum(z*z for z in u))*math.sqrt(sum(z*z for z in v))); return math.degrees(math.acos(max(-1,min(1,co))))
attack_d=dist(mapping['thr267_OG1'],mapping['l2_C']); attack_a=angle(mapping['thr267_OG1'],mapping['l2_C'],mapping['l2_O'])
assert attack_d<=0.35 and 95<=attack_a<=115
(build/'topol_Protein_chain_H.itp').write_text(itp)
(build/'system_M1.gro').write_text(full_text)
for fname,k in [('posre_Protein_chain_H.itp',1000),('posre_Protein_chain_H_weak.itp',100)]:
    lines=['; generated for M1 NalphaH2/Asp306H chain H','[ position_restraints ]','; atom type fx fy fz']
    lines += [f"{x['nr']:6d} 1 {k} {k} {k}" for x in new_atoms if heavy(x['name'])]
    (build/fname).write_text('\n'.join(lines)+'\n')
result={'schema_version':1,'status':'PASS_M1_BUILD_AUDIT_PRE_GROMPP','candidate_id':'nylc_C18_trueT267_freeGS','source_job_id':'61780674','source_selected_time_ps':selected_time,'microstate':{'Thr267_Nalpha':'NH2_neutral_proxy','Thr267_Ogamma':'OH','Asp306':'ASH_HD2_neutral','Asp308':'ASP_minus'},'forcefield':'amber99sb-ildn','parameterization_note':'Nalpha neutral amine proxy retains N3/H bonded types; N charge -0.62540 gives exact residue neutralization and N+2H group charge -0.23860 e, close to force-field neutral LYN amine group -0.26373 e; not a new RESP fit.','old_chain_charge_e':round(old_charge,6),'new_chain_charge_e':round(new_charge,6),'full_atom_count':len(merged_atoms),'active_chain_atom_count':len(new_atoms),'atom_mapping_global_index1':mapping,'reaction_geometry':{'attack_distance_nm':attack_d,'attack_angle_deg':attack_a},'itp_rewrite':a_itp,'gro_rewrite':a_gro,'merge':a_merge,'sha256':{'source_selected_gro':sha(full_gro),'system_M1_gro':sha(build/'system_M1.gro'),'chainH_itp':sha(build/'topol_Protein_chain_H.itp')}}
audit_path.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
print(json.dumps(result,sort_keys=True))
PY
(
  cd "$OUT/build"
  "$GMX" grompp -f "$FLOW/mdp/em.mdp" -c system_M1.gro -r system_M1.gro -p topol.top \
    -o preflight.tpr -po preflight.processed.mdp -maxwarn 0 >grompp.stdout 2>grompp.stderr
)
"$PY" - "$OUT/build/system_M1.gro" "$OUT/build/chain_rest.ndx" <<'PY'
import pathlib,sys
p=pathlib.Path(sys.argv[1]); n=int(p.read_text().splitlines()[1]); out=pathlib.Path(sys.argv[2])
def block(name,values):
    lines=[f'[ {name} ]']; vals=list(values)
    lines += [' '.join(str(x) for x in vals[i:i+15]) for i in range(0,len(vals),15)]; return '\n'.join(lines)
out.write_text(block('ChainH',range(8949,10273))+'\n'+block('Rest',list(range(1,8949))+list(range(10273,n+1)))+'\n')
PY
printf '0\n1\n' | "$GMX" mindist -f "$OUT/build/system_M1.gro" -s "$OUT/build/preflight.tpr" \
  -n "$OUT/build/chain_rest.ndx" -od "$OUT/build/chain_rest_mindist.xvg" \
  >"$OUT/build/mindist.stdout" 2>"$OUT/build/mindist.stderr"
"$PY" - "$OUT/selection/nac_energy_series.json" "$OUT/build_audit.json" "$OUT/build/chain_rest_mindist.xvg" "$OUT/PASS.json" "$job" <<'PY'
import json,sys
series_path,build_path,mindist_path,out_path,job=sys.argv[1:]
series=json.load(open(series_path)); build=json.load(open(build_path))
vals=[]
for line in open(mindist_path):
    if line.strip() and not line.lstrip().startswith(('#','@')): vals.append(float(line.split()[1]))
assert vals and min(vals)>=0.08
chosen=series['lowest_potential_nac_frame']; assert chosen is not None
payload={'schema_version':1,'status':'PASS_TECHNICAL_M1_BUILD_GROMPP','scientific_status':'M1_CLASSICAL_PROXY_READY_FOR_EM_NOT_MECHANISM_PASS','candidate_id':'nylc_C18_trueT267_freeGS','slurm_job_id':job,'source_selection':{'rule':'lowest instantaneous potential among fully unrestrained NAC frames','selected':chosen,'trajectory_frame_count':series['frame_count'],'nac_frame_count':series['nac_frame_count'],'nac_occupancy':series['nac_occupancy']},'build_audit':build,'minimum_chainH_rest_distance_nm':min(vals),'grompp_maxwarn':0,'next_gate':'EM then restrained low-temperature heating/release then fully unrestrained M1 NPT'}
open(out_path,'w').write(json.dumps(payload,indent=2,sort_keys=True)+'\n')
PY
record PASS_TECHNICAL "status=PASS_TECHNICAL_M1_BUILD_GROMPP;microstate=NalphaH2_Asp306H;out=$OUT/PASS.json"
trap - ERR