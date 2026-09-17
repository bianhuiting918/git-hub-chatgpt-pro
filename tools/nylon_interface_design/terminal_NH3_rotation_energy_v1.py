"""CPU-only fixed-heavy-atom NH3 orientation diagnostic; no structure repair."""
import json, math, hashlib, datetime
from pathlib import Path
import numpy as np
import pyrosetta
from pyrosetta import rosetta
R=Path(__file__).resolve().parent
D=R/'route56_N10_20260915_v1'
OUT=D/'terminal_NH3_rotation_energy_v1'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def xyz(v): return np.array([v.x,v.y,v.z],dtype=float)
def angle(a,b,c):
 u=a-b; v=c-b
 return float(np.degrees(np.arccos(np.clip(np.dot(u,v)/np.linalg.norm(u)/np.linalg.norm(v),-1,1))))
def run():
 OUT.mkdir(exist_ok=False)
 (OUT/'RUNBOOK.md').write_text("CPU PyRosetta ref2015_cart; cart_bonded=0.5, pro_close=0. No minimization or restraints. Rotate only three terminal N hydrogens around N-CA, 0..355 deg by 5 deg. Preserve all heavy atoms and other hydrogens. Score is Rosetta energy, only within-case comparisons. Original PDBs unchanged. Geometry uses same Asp1 oxygen nearest Thr N; NO 2.5-3.5 A, HO <=2.6 A, NHO >=130 deg. Not donor-pair, folding, material or activity validation. Rerun with a NEW output version; never overwrite. Command: CPU PyRosetta Python terminal_NH3_rotation_energy_v1.py\n")
 audit=json.loads((R/'route4_same_Asp_oxygen_audit_20260917T063959469270Z.json').read_text())
 cases=[]
 for name in ['N10','N9_R4_B37_1_x20_20']:
  rec=next(r for r in audit['records'] if r['name']==name and r['stage']=='after' and r['seed']==101)
  cases.append((name,Path(rec['path']),rec['source_sha256'],rec['core']))
 cases.append(('NylC',D/'NylC_hydrogens_no_relax_v1.pdb','5701314028a0df85ed9072ffd3deafdb12d1b96f4d5aa51b719fb53e863d49d2',[267,306,308]))
 pyrosetta.init('-mute all -constant_seed -jran 101 -ignore_unrecognized_res false')
 sf=pyrosetta.create_score_function('ref2015_cart')
 sf.set_weight(rosetta.core.scoring.cart_bonded,0.5)
 sf.set_weight(rosetta.core.scoring.pro_close,0.0)
 summaries=[]
 for name,path,expected,core in cases:
  assert sha(path)==expected,(name,'SOURCE_SHA_CHANGED')
  p=pyrosetta.pose_from_pdb(str(path))
  ids=[]
  for number,resname in zip(core,['THR','ASP','ASP']):
   matches=[i for i in range(1,p.total_residue()+1) if p.pdb_info().number(i)==number and p.residue(i).name3()==resname]
   assert len(matches)==1,(name,number,matches)
   ids.append(matches[0])
  ti,di,d2i=ids
  tr=p.residue(ti); ni=tr.atom_index('N')
  hs=[i for i in range(1,tr.natoms()+1) if tr.atom_is_hydrogen(i) and tr.atom_base(i)==ni]
  assert len(hs)==3,(name,hs,tr.type().name())
  # Explicit restoration of supplied coordinates prevents loader-built target H from becoming an unreported input.
  original={}
  for line in path.read_text().splitlines():
   if line.startswith(('ATOM  ','HETATM')) and int(line[22:26])==core[0] and line[17:20].strip()=='THR':
    original[line[12:16].strip()]=np.array([float(line[30:38]),float(line[38:46]),float(line[46:54])])
  for hi in hs:
   an=tr.atom_name(hi).strip()
   assert an in original,(name,'MISSING_ORIGINAL_H',an)
   p.set_xyz(rosetta.core.id.AtomID(hi,ti),rosetta.numeric.xyzVector_double_t(*original[an]))
  fixed=[(i,j,xyz(p.residue(i).xyz(j))) for i in range(1,p.total_residue()+1) for j in range(1,p.residue(i).natoms()+1) if not (i==ti and j in hs)]
  n=xyz(p.residue(ti).xyz('N')); ca=xyz(p.residue(ti).xyz('CA'))
  axis=(ca-n)/np.linalg.norm(ca-n)
  vectors=[xyz(p.residue(ti).xyz(h))-n for h in hs]
  oname=min(['OD1','OD2'],key=lambda a:np.linalg.norm(xyz(p.residue(di).xyz(a))-n))
  o=xyz(p.residue(di).xyz(oname)); no=float(np.linalg.norm(n-o))
  rows=[]
  for deg in range(0,360,5):
   q=p.clone(); t=math.radians(deg); contacts=[]
   for hi,v in zip(hs,vectors):
    w=v*math.cos(t)+np.cross(axis,v)*math.sin(t)+axis*np.dot(axis,v)*(1-math.cos(t))
    h=n+w
    q.set_xyz(rosetta.core.id.AtomID(hi,ti),rosetta.numeric.xyzVector_double_t(*h))
    ho=float(np.linalg.norm(h-o)); a=angle(n,h,o)
    contacts.append({'H':tr.atom_name(hi).strip(),'HO_A':ho,'NHO_deg':a,'pass':2.5<=no<=3.5 and ho<=2.6 and a>=130})
   q.energies().clear()
   score=float(sf(q))
   assert math.isfinite(score)
   displacement=max(float(np.linalg.norm(xyz(q.residue(i).xyz(j))-v)) for i,j,v in fixed)
   assert displacement<1e-10,(name,deg,'FIXED_ATOMS_MOVED')
   components={term:float(q.energies().total_energies()[getattr(rosetta.core.scoring,term)]) for term in ['fa_rep','fa_elec','hbond_sc','hbond_bb_sc','cart_bonded']}
   rows.append({'deg':deg,'score_REU':score,'contacts':contacts,'direction_pass':any(c['pass'] for c in contacts),'fixed_atom_max_shift_A':displacement,'unweighted_terms':components})
  minimum=min(rows,key=lambda r:r['score_REU'])
  passing=[r for r in rows if r['direction_pass']]
  bestpass=min(passing,key=lambda r:r['score_REU']) if passing else None
  for row in rows: row['delta_from_min_REU']=row['score_REU']-minimum['score_REU']
  result={'name':name,'input':str(path),'input_sha256':expected,'core':core,'selected_Asp1_O':oname,'NO_A':no,'hydrogens':[tr.atom_name(i).strip() for i in hs],'rows':rows,'baseline':rows[0],'minimum':minimum,'lowest_energy_direction_pass':bestpass,'passing_orientations':len(passing),'orientations':len(rows),'limitations':['Fixed protein and all other hydrogens; not free energy or pH sampling','No biological activity or complete catalytic pass inference','NylC default-added H and relaxed N10/candidate H have different preparation histories']}
  (OUT/(name+'.json')).write_text(json.dumps(result,indent=2))
  brief={k:result[k] for k in ['name','NO_A','passing_orientations','orientations']}
  brief.update(baseline_angle=max(c['NHO_deg'] for c in rows[0]['contacts']),minimum_deg=minimum['deg'],minimum_direction_pass=minimum['direction_pass'],minimum_angle=max(c['NHO_deg'] for c in minimum['contacts']),best_pass_delta_REU=None if bestpass is None else bestpass['delta_from_min_REU'])
  summaries.append(brief)
  print(json.dumps(brief),flush=True)
 (OUT/'summary.json').write_text(json.dumps({'status':'COMPLETE_FIXED_HEAVY_NH3_DIAGNOSTIC_NOT_CATALYTIC_PASS','cases':summaries},indent=2))
 return summaries
if __name__=='__main__':
 start=datetime.datetime.now(datetime.timezone.utc).isoformat()
 try:
  results=run(); status='COMPLETE'; error=None
 except Exception as exc:
  status='FAILED';error=repr(exc);raise
 finally:
  with (R/'RUN_LOG.jsonl').open('a') as f:
   f.write(json.dumps({'time':start,'script':__file__,'script_sha256':sha(Path(__file__)),'output':str(OUT),'status':status,'error':error,'parameters':{'rotation_step_deg':5,'heavy_atoms':'fixed','minimization':False,'GPU':False}})+'\n')
