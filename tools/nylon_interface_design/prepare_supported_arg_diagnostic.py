"""Prepare a supported ARG motif; CPU only, no GPU submission.
Run from the existing remote N9/N10 project using its CPU environment.
Creates a new output version and appends RUN_LOG.jsonl. Existing results
are never overwritten. Local geometry is not full-protein validation.
"""
import os,datetime,hashlib
import pathlib,json,numpy as np
from scipy.optimize import least_squares
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
D=pathlib.Path('route56_N10_20260915_v1');p=D/'route6_N10_single_residue_pilot_v1';panel=json.load(open(p/'panel.json'));k=next(i for i,v in enumerate(panel) if v['id']=='SC4559268');z=np.load(p/'ARG_1.npz');q=z['q'][int(np.where(z['indices']==k)[0][0])];src=pathlib.Path('/data/bht2/new_seed_structure_sequence_folddisco_20260825_v3/n9_n10_af_repredict_20260908_v1/outputs/nonfine_af_v1/N10_AF2fresh/seed_101/mature_relaxed_nh3.pdb');res={};bb=['N','CA','C','O']
for l in src.read_text().splitlines():
 if l.startswith('ATOM  ') and l[21]=='A' and 170<=int(l[22:26])<=174 and l[12:16].strip() in bb:res.setdefault(int(l[22:26]),{})[l[12:16].strip()]=[float(l[30:38]),float(l[38:46]),float(l[46:54])]
T=np.array([res[n][b] for n in range(170,175) for b in bb]);ai=[5,6,7,8,9,10,11,12,13];target=q[[27,28,29,30,31,39,40,41,42]];A=T[ai];u,s,v=np.linalg.svd((A-A.mean(0)).T@(target-target.mean(0)));R=u@np.diag([1,1,np.linalg.det(u@v)])@v;P=(T-A.mean(0))@R+target.mean(0);P[ai]=target;before=P.copy();fixed=[8,9,10,11];free=np.array([j for j in range(20) if j not in fixed]);allb=[]
for j in range(5):
 allb += [(4*j,4*j+1),(4*j+1,4*j+2),(4*j+2,4*j+3)]
 if j<4:allb.append((4*j+2,4*j+4))
nei={i:[] for i in range(20)}
for i,j in allb:nei[i].append(j);nei[j].append(i)
angles=[(i,j,k) for j in range(20) for ii,i in enumerate(nei[j]) for k in nei[j][ii+1:]]
def ang(p,i,j,k):
 a=p[i]-p[j];b=p[k]-p[j];return np.arccos(np.clip(np.dot(a,b)/np.linalg.norm(a)/np.linalg.norm(b),-1,1))
bonds=[b for b in allb if not all(t in fixed for t in b)];angles=[b for b in angles if not all(t in fixed for t in b)];bd=np.array([np.linalg.norm(T[i]-T[j]) for i,j in bonds]);at=np.array([ang(T,*b) for b in angles]);planes=[(4*j+1,4*j+2,4*j+4,4*j+5) for j in range(4)]+[(4*j+3,4*j+2,4*j+4,4*j+5) for j in range(4)]
def tors(p,i,j,k,l):
 a=p[j]-p[i];b=p[k]-p[j];c=p[l]-p[k];n=np.cross(a,b);m=np.cross(b,c);return np.arctan2(np.dot(np.cross(n,m),b/np.linalg.norm(b)),np.dot(n,m))
tt=np.array([tors(T,*b) for b in planes]);ref=before[free].copy()
def residual(x):
 p=before.copy();p[free]=x.reshape(-1,3);d=np.array([np.linalg.norm(p[i]-p[j]) for i,j in bonds]);aa=np.array([ang(p,*b) for b in angles]);dt=np.array([tors(p,*b) for b in planes])-tt;dt=np.arctan2(np.sin(dt),np.cos(dt));return np.r_[(d-bd)/.02,(aa-at)/np.deg2rad(3),dt/np.deg2rad(5),((p[free]-ref)/2).ravel()]
sol=least_squares(residual,before[free].ravel(),max_nfev=150,ftol=1e-9,xtol=1e-9,gtol=1e-9);P[free]=sol.x.reshape(-1,3)
rad={'C':1.7,'N':1.55,'O':1.52,'S':1.8};rbb=np.tile([1.55,1.7,1.7,1.52],5);side=q[32:39];rside=z['r'][32:39];fragment=np.vstack([P,side]);rr=np.r_[rbb,rside];graph=np.full((27,27),1000);np.fill_diagonal(graph,0)
for i,j in allb+[(9,20),(20,21),(21,22),(22,23),(23,24),(24,25),(24,26)]:graph[i,j]=graph[j,i]=1
for j in range(27):graph=np.minimum(graph,graph[:,j,None]+graph[None,j,:])
m=cdist(fragment,fragment)-rr[:,None]-rr[None,:]+.4;m[graph<=3]=np.inf;self_bad=np.argwhere(np.triu(m< -1e-8,1));cm=cdist(fragment,q[:27])-rr[:,None]-z['r'][None,:27]+.4
cov=np.load(p/'coverage.npz')['coverage'][k];sites=json.load(open(D/'N10_fixed_minimal_motif_material_v1/sites.json'));report=[];hitall=[]
for mat in ['PA6','PA66']:
 e=np.load(pathlib.Path('../interface_surface_robustness_20260908_v1/inputs_v2')/(mat+'_parent_atoms.npz'));px=e['xyz_A'].copy();box=e['box_A'];px[:,:2]%=box[:2];trees={el:cKDTree(px[e['elements']==el],boxsize=[*box[:2],0]) for el in np.unique(e['elements']) if el!='H'};counts={}
 for sj,site in enumerate(sites):
  if site['material']!=mat or not cov[sj]:continue
  w=fragment@np.array(site['frame']).T+site['center'];w[:,:2]%=box[:2];margin=np.full(len(w),np.inf)
  for el,tr in trees.items():margin=np.minimum(margin,tr.query(w)[0]-rr-rad[str(el)]+.4)
  if margin.min()>=-1e-8:hitall.append(sj);counts[site['bin']]=counts.get(site['bin'],0)+1
 report.append({'material':mat,'parent_eligible_sites':sum(bool(cov[j]) and s['material']==mat for j,s in enumerate(sites)),'repaired_fragment_clear':sum(counts.values()),'bin_counts':counts})
omega=[180-abs(np.degrees(tors(P,*b))) for b in planes[:4]]
print(json.dumps({'id':'SC4559268','template':'N10 A170-174','optimizer_success':bool(sol.success),'fixed_backbone_displacement':float(np.linalg.norm(P[fixed]-before[fixed],axis=1).max()),'sidechain_and_TDD_unchanged':True,'self_min_margin_A':float(m.min()),'self_bad_pairs':self_bad.tolist(),'fragment_vs_TDD_Gly298_min_margin_A':float(cm.min()),'omega_trans_deviation_max_deg':max(omega),'CN_A':[float(np.linalg.norm(P[4*j+2]-P[4*j+4])) for j in range(4)],'coverage':report,'scope':'CPU in-memory geometric feasibility diagnostic; only old eligible sites replayed; full scaffold and peripheral sidechains not evaluated'},indent=2),flush=True)

from scipy.interpolate import RegularGridInterpolator
jobs=json.load(open(D/'low_coverage_rfd12_v2/manifest.json'))['jobs']
for job in jobs:
 if 'PA6_SC4559268' not in job['name']:continue
 pp=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in pathlib.Path(job['pdb']).read_text().splitlines() if l.startswith('ATOM  ')])
 srcq=q[:23];dstq=np.vstack([pp[:15],pp[19:27]]);u,s,v=np.linalg.svd((srcq-srcq.mean(0)).T@(dstq-dstq.mean(0)));rot=u@np.diag([1,1,np.linalg.det(u@v)])@v;shift=dstq.mean(0)-srcq.mean(0)@rot;err=np.linalg.norm(srcq@rot+shift-dstq,axis=1);assert max(err)<.002
 f=np.load(pathlib.Path(job['pdb']).parent/'shape_field.npz');ca=np.vstack([q[[1,8,24,16]],P[[1,5,9,13,17]]])@rot+shift
 values=RegularGridInterpolator([f['axis_'+s+'_A'] for s in 'xyz'],f['field'],bounds_error=False,fill_value=-100)(ca)
 print(json.dumps({'interface':job['name'],'selected_count':len(job['selected_site_indices']),'selected_material_clear_after_repair':sum(j in hitall for j in job['selected_site_indices']),'rigid_core_fit_max_A':float(err.max()),'new_fixed_CA_field_A':values.tolist(),'new_fixed_CA_outside_count':int(sum(values<=0))}),flush=True)

# Serialize the checked local fragment and retain the original CAD fields.
project=pathlib.Path.cwd().resolve()
assert (project/'route56_N10_20260915_v1').is_dir()
out=project/'route56_N10_20260915_v1/route6_supported_ARG_SC4559268_v1'
out.mkdir(exist_ok=False)
records=[]
assert sol.success and m.min()>=-1e-8 and cm.min()>=-1e-8
assert max(omega)<=30 and np.array_equal(P[fixed],before[fixed])
assert all(1.2<=np.linalg.norm(P[4*j+2]-P[4*j+4])<=1.5 for j in range(4))
try:
 for job in jobs:
  if 'PA6_SC4559268' not in job['name']:continue
  assert all(j in hitall for j in job['selected_site_indices'])
  pp=np.array([[float(l[30:38]),float(l[38:46]),float(l[46:54])] for l in pathlib.Path(job['pdb']).read_text().splitlines() if l.startswith('ATOM  ')])
  srcq=q[:23];dstq=np.vstack([pp[:15],pp[19:27]])
  u,s,v=np.linalg.svd((srcq-srcq.mean(0)).T@(dstq-dstq.mean(0)));rot=u@np.diag([1,1,np.linalg.det(u@v)])@v;shift=dstq.mean(0)-srcq.mean(0)@rot
  assert np.linalg.norm(srcq@rot+shift-dstq,axis=1).max()<.002
  blocks=[('THR',['N','CA','C','O','CB','OG1','CG2'],q[:7]),('ASP',['N','CA','C','O','CB','CG','OD1','OD2'],q[7:15]),('GLY',bb,q[23:27]),('ASP',['N','CA','C','O','CB','CG','OD1','OD2'],q[15:23])]
  for j in range(5):
   if j==2:blocks.append(('ARG',bb+['CB','CG','CD','NE','CZ','NH1','NH2'],np.vstack([P[8:12],q[32:39]])))
   else:blocks.append(('GLY',bb,P[4*j:4*j+4]))
  dst=out/job['name'];dst.mkdir();lines=[];serial=0;ca=[]
  for ri,(rn,names,coords) in enumerate(blocks,1):
   for name,xyz in zip(names,coords@rot+shift):
    serial+=1;line=f'ATOM  {serial:5d} {name:^4s} {rn:3s} A{ri:4d}    {xyz[0]:8.3f}{xyz[1]:8.3f}{xyz[2]:8.3f}{1:6.2f}{0:6.2f}          {name[0]:>2s}';lines.append(line)
    if name=='CA':ca.append([float(line[30:38]),float(line[38:46]),float(line[46:54])])
  pdb=dst/'motif.pdb';pdb.write_text('\n'.join(lines)+'\nEND\n')
  f=np.load(pathlib.Path(job['pdb']).parent/'shape_field.npz')
  fieldvals=RegularGridInterpolator([f['axis_'+s+'_A'] for s in 'xyz'],f['field'],bounds_error=False,fill_value=-100)(ca)
  assert np.all(fieldvals>0)
  oldcenter=f['centering_roundoff_A'];newcenter=np.mean(ca,axis=0);obj=dst/'shape.obj';mesh=[]
  for line in pathlib.Path(job['shape']).read_text().splitlines():
   if line.startswith('v '):
    xyz=np.array([float(t) for t in line.split()[1:4]])+oldcenter-newcenter
    mesh.append('v %.6f %.6f %.6f'%tuple(xyz))
   else:mesh.append(line)
  obj.write_text('\n'.join(mesh)+'\n')
  np.savez_compressed(dst/'motif_exact.npz',core_original_A=q,fragment_original_A=fragment,rotation_to_shape=rot,translation_to_shape=shift,new_center_A=newcenter,old_center_A=oldcenter)
  output=project.parent/'n9_core_donor_design_20260910_gpu_outputs_v1/route6_supported_ARG_SC4559268_v1'/job['name']
  output.mkdir(parents=True,exist_ok=False,mode=0o2770);output.chmod(0o2770)
  cmd=[]
  for arg in job['command']:
   if arg.startswith('inference.input_pdb='):arg='inference.input_pdb='+str(pdb)
   elif arg.startswith('inference.output_prefix='):arg='inference.output_prefix='+str(output/'design')
   elif arg.startswith('contigmap.contigs='):arg='contigmap.contigs=[A1-1/38-38/A2-4/18-18/A5-9/61-61]'
   elif arg.startswith('potentials.guiding_potentials='):
    arg=arg.replace(job['shape'],str(obj))
   cmd.append(arg)
  assert cmd[0]=='/usr/local/bin/gpurun'
  rec={'name':job['name'],'status':'PREPARED_NOT_SUBMITTED','command':cmd,'output':str(output),'pdb':str(pdb),'shape':str(obj),'pdb_sha256':hashlib.sha256(pdb.read_bytes()).hexdigest(),'shape_sha256':hashlib.sha256(obj.read_bytes()).hexdigest(),'original_manifest':str(project/'route56_N10_20260915_v1/low_coverage_rfd12_v2/manifest.json'),'selected_site_indices':job['selected_site_indices'],'conditional_denominator':52,'fixed_output_residues':[1,40,41,42,61,62,63,64,65],'catalytic_output_residues':[1,40,42,63],'support_identity_note':'Gly placeholders are not required final sequence identities','new_fixed_CA_field_A':fieldvals.tolist(),'limitations':['Peripheral sidechains absent','Only original eligible sites replayed','No RFD output or independent prediction has been evaluated','Shape remains a soft generation guide']}
  (dst/'manifest.json').write_text(json.dumps(rec,indent=2));records.append(rec)
 (out/'manifest.json').write_text(json.dumps({'status':'PREPARED_NOT_SUBMITTED','jobs':records},indent=2))
 (out/'RUNBOOK.md').write_text('CPU preparation: run this downloaded script from the remote project with ../unbiased_100chain_v1/env/bin/python. Outputs are exclusive-create and cannot be overwritten. Each manifest records a gpurun command for a separate, explicitly audited generation. First test coverage10 only; inspect raw motif retention, all backbone bonds/angles, clashes and material replay before continuing. No automatic GPU retries. Preserve all failed outputs. This is a local-fragment diagnostic, not a complete enzyme.\n')
 print(json.dumps({'prepared':len(records),'output':str(out),'status':'PREPARED_NOT_SUBMITTED'}),flush=True)
finally:
 with (project/'RUN_LOG.jsonl').open('a') as log:
  log.write(json.dumps({'time':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script':str(pathlib.Path(__file__).resolve()),'script_sha256':hashlib.sha256(pathlib.Path(__file__).read_bytes()).hexdigest(),'output':str(out),'prepared_cases':len(records),'status':'PREPARED_NOT_SUBMITTED' if len(records)==3 else 'PREPARATION_INCOMPLETE'})+'\n')
