#!/usr/bin/env python3
"""Audit v10 nylon atoms/bonds against original SDF and dry export; no pose sampling."""
from pathlib import Path
import sys,json,hashlib,datetime
import numpy as np
from scipy.spatial import cKDTree
ROOT=Path(__file__).resolve().parent
BASE=ROOT.parent
SRC=BASE/'surface_previews/axis_density_sixface_all_external_trim_v6'
sys.path.insert(0,str(SRC))
from density_cut_variants import parse_pdb
from terminal_density_trim_3d import center_whole_chains_in_primary_cell
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def atoms_edges(p):
 lines=p.read_text().splitlines()
 atoms={l[6:11]:l for l in lines if l.startswith(('ATOM  ','HETATM'))}
 edges=set()
 for l in lines:
  if l.startswith('CONECT'):
   fields=[l[i:i+5] for i in range(6,len(l),5)]
   for v in fields[1:]:
    if v.strip(): edges.add(tuple(sorted((fields[0],v))))
 assert all(a in atoms and b in atoms for a,b in edges)
 return atoms,edges
def audit(mat,dp,expected):
 inp=ROOT/'inputs'/f'{mat}_v10_ge90_uncapped.pdb'
 parent=BASE/f'{mat.lower()}_dp{dp}_400chain/water_slab_v1/dry_export/{mat}_DP{dp}_400chain_water_equilibrated_dry.pdb'
 sdf=BASE.parent/'unbiased_100chain_v1/00_matched_chains'/f'{mat}_DP{dp}.sdf'
 assert sha(inp)==expected
 _,box,pa,_=parse_pdb(parent)
 center_whole_chains_in_primary_cell(pa,box)
 kept,ke=atoms_edges(inp)
 original,oe=atoms_edges(parent)
 assert set(kept)<=set(original)
 assert ke=={e for e in oe if all(a in kept for a in e)}
 byid={a['serial']:a for a in pa}
 maxerr=0.
 for key,l in kept.items():
  a=byid[key]; assert l[12:30]==a['line'][12:30] and l[72:78]==a['line'][72:78]
  maxerr=max(maxerr,float(np.max(np.abs(np.array([float(l[i:i+8]) for i in (30,38,46)])-a['xyz']))))
 assert maxerr<0.002, maxerr
 lines=sdf.read_text().splitlines(); na=int(lines[3][:3]); nb=int(lines[3][3:6]); assert na==316
 el=[l[31:34].strip() for l in lines[4:4+na]]
 bonds=[(int(l[:3])-1,int(l[3:6])-1,int(l[6:9])) for l in lines[4+na:4+na+nb]]
 adj=[[] for _ in el]
 for a,b,o in bonds: adj[a].append((b,o)); adj[b].append((a,o))
 motifs=[]
 for c in range(na):
  if el[c]!='C': continue
  oxy=[b for b,o in adj[c] if el[b]=='O' and o==2]
  ns=[b for b,o in adj[c] if el[b]=='N' and o==1]
  if not oxy: continue
  for n in ns:
   assert len(oxy)==1
   motifs.append((c,oxy[0],n))
 assert len(motifs)==17
 # Every parent bond must agree with replicated SDF connectivity.
 expected_edges={tuple(sorted((pa[k*316+a]['serial'],pa[k*316+b]['serial']))) for k in range(400) for a,b,_ in bonds}
 assert expected_edges==oe
 cut=[e for e in oe if sum(a in kept for a in e)==1]
 cutids=sorted({a for e in cut for a in e if a in kept})
 xyz=lambda a: np.array([float(kept[a][i:i+8]) for i in (30,38,46)])
 tree=cKDTree([xyz(a) for a in cutids]) if cutids else None
 rows=[]; states={'retained':0,'one_sided_cut':0,'removed':0}
 for k in range(400):
  for j,(c,o,n) in enumerate(motifs):
   ids=[pa[k*316+a]['serial'] for a in (c,o,n)]
   cc,oo,nn=ids; exists=[a in kept for a in ids]
   if all(exists):
    states['retained']+=1
    d=float(tree.query(xyz(cc))[0]) if tree else None
    rows.append(dict(site_id=f'P{k+1:03d}_{pa[k*316+c]["line"][12:16].strip()}_{pa[k*316+o]["line"][12:16].strip()}',serials=ids,original_chain=k+1,amide_index=j+1,carbonyl_xyz_A=xyz(cc).tolist(),nearest_retained_cut_endpoint_A=d,cut_endpoint_within_15A=d is not None and d<=15))
   elif exists[0] != exists[2]: states['one_sided_cut']+=1
   else: states['removed']+=1
 assert sum(states.values())==6800
 assert states['one_sided_cut']==len(cut)
 result=dict(material=mat,input_sha256=sha(inp),parent_sha256=sha(parent),sdf_sha256=sha(sdf),atoms=len(kept),bonds=len(ke),coordinate_mapping_max_error_A=maxerr,parent_amides=6800,amide_states=states,retained_cut_endpoints=len(cutids),intact_amides_with_cut_endpoint_within_15A=sum(r['cut_endpoint_within_15A'] for r in rows),status='INPUT_CONNECTIVITY_AUDITED_SURFACE_NOT_EVALUATED',boundary='finite v10 trimmed object; original whole-chain translations reproduced; no periodic copies introduced',sites=rows)
 return result
if __name__=='__main__':
 out=ROOT/'input_audit_v2'; out.mkdir(exist_ok=False)
 try:
  results=[audit('PA6',16,'1e93b32c62f442c34531ca65767cffbeffd26375b9c0f786db8a48e4b1a630bf'),audit('PA66',8,'125963a9af3f6d1cb3b62e567cb6247bf7adde330f6489933bbd2a5c12087fe4')]
  for r in results: (out/(r['material']+'_intact_amides.json')).write_text(json.dumps(r,indent=2))
  summary=[{k:v for k,v in r.items() if k!='sites'} for r in results]
  (out/'summary.json').write_text(json.dumps(summary,indent=2))
  print(json.dumps(summary,indent=2))
  status='PASS'; exitcode=0
 except BaseException as e:
  status=repr(e); exitcode=1
  raise
 finally:
  with (ROOT/'RUN_LOG.jsonl').open('a') as f:
   f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__).resolve()),script_sha256=sha(Path(__file__)),command='CPU Python audit_inputs.py',inputs='inputs/ v10 PDBs; original dry exports; original matched-chain SDFs',outputs=str(out),parameters=dict(cut_proximity_A=15,mapping_tolerance_A=0.002),exit_code=exitcode,result=status))+'\n')
