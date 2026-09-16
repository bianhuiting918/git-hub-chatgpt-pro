"""Read-only structure evaluation with versioned reports; never predicts missing structures."""
from pathlib import Path
import numpy as np,json,hashlib,datetime
from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
r=Path(__file__).resolve().parent;o=r/'route6_supported_arg_af20_v1'
out=r.parent/'n9_core_donor_design_20260910_gpu_outputs_v1/af2_route6_supported_arg20_v1'
queries=json.loads((o/'sequence_mapping.json').read_text());parser=PDBParser(QUIET=True)
refpath=r/'route56_N10_20260915_v1/route6_supported_ARG_SC4559268_v1/fixed_core_repair_v1/after.pdb'
ref=parser.get_structure('reference',refpath)[0]['A'];rows=[];pending=[]
def fit(a,ids):
 keys=[(i,n) for i in ids for n in ['N','CA','C','O']]
 x=np.array([a[i][n].coord for i,n in keys]);y=np.array([ref[i][n].coord for i,n in keys])
 u,sv,vt=np.linalg.svd((x-x.mean(0)).T@(y-y.mean(0)));fix=np.eye(3);fix[-1,-1]=np.linalg.det(u@vt);rot=u@fix@vt
 trans=y.mean(0)-x.mean(0)@rot
 return rot,trans,float(np.sqrt(np.mean(np.sum((x@rot+trans-y)**2,axis=1))))
for name,seq in queries.items():
 if not (out/(name+'.done.txt')).exists():pending.append(name);continue
 fs=list(out.glob(name+'_unrelaxed*.pdb'));assert len(fs)==1
 f=fs[0];a=parser.get_structure(name,f)[0]['A'];assert ''.join(seq1(v.resname) for v in a)==seq
 def dist(i,n,j,m):return float(np.linalg.norm(a[i][n].coord-a[j][m].coord))
 chosen=min(['OD1','OD2'],key=lambda n:dist(1,'N',40,n))
 ds=[dist(1,'N',40,chosen),dist(1,'N',1,'OG1'),min(dist(40,chosen,42,m) for m in ['OD1','OD2'])]
 lo=np.array([2.5,2.5,3.3]);hi=np.array([3.5,3.3,4.7]);gap=np.maximum(lo-ds,0)+np.maximum(np.array(ds)-hi,0)
 rot,trans,core_rms=fit(a,[1,40,42]);_,_,all_rms=fit(a,[1,40,41,42,61,62,63,64,65])
 donors={n:float(np.linalg.norm(a[63][n].coord@rot+trans-ref[63][n].coord)) for n in ['NH1','NH2']}
 cn=[dist(i,'C',i+1,'N') for i in range(1,126)]
 rows.append(dict(name=name,source=str(f),sha256=hashlib.sha256(f.read_bytes()).hexdigest(),sequence_match=True,
 mean_CA_pLDDT=float(np.mean([v['CA'].bfactor for v in a])),TDD_backbone_RMSD_A=core_rms,all_fixed_backbone_RMSD_A=all_rms,
 TDD_distances_A=ds,TDD_Asp1_selected=chosen,TDD_distance_gate=bool(np.all(gap==0)),TDD_distance_gap_sum_A=float(gap.sum()),
 donor_displacement_after_TDD_alignment_A=donors,donor_NN_A=dist(63,'NH1',63,'NH2'),CN_range_A=[min(cn),max(cn)],
 bad_CN_count=sum(not 1.2<=v<=1.5 for v in cn),donor_H_direction='NOT_EVALUATED',material_accessibility='NOT_EVALUATED',status='NOT_EVALUATED_JOINT_GEOMETRY'))
stamp=datetime.datetime.now(datetime.timezone.utc).strftime('%Y%m%dT%H%M%SZ')
report=dict(status='PARTIAL_BASIC_GEOMETRY_AUDIT' if pending else 'COMPLETE_BASIC_GEOMETRY_AUDIT',expected=20,audited=len(rows),pending=pending,rows=rows,
reference=str(refpath),reference_sha256=hashlib.sha256(refpath.read_bytes()).hexdigest(),
limitations=['Exploratory source has known peptide and close-contact defects','One-model one-seed three-recycle predictions are not activity evidence','Arg atom-specific displacements can be affected by equivalent NH1/NH2 naming','TDD-only fit is a diagnostic alignment, not material-accessibility validation','Selected Asp oxygen is shared by the two reported interresidue TDD distances'])
dest=o/('prediction_audit_'+stamp+'.json')
with dest.open('x') as h:json.dump(report,h,indent=2)
with (r/'RUN_LOG.jsonl').open('a') as h:h.write(json.dumps(dict(time=stamp,script=str(Path(__file__)),script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),output=str(dest),audited=len(rows),pending=len(pending),exit_code=0))+'\n')
print(json.dumps(dict(output=str(dest),audited=len(rows),pending=len(pending),distance_pass=sum(x['TDD_distance_gate'] for x in rows))))
