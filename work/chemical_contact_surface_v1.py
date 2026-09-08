#!/usr/bin/env python
"""Chemical polymer SES atlas and explicit fixed-core contact-face candidates."""
import sys,json,time,base64,zlib,itertools,argparse,datetime
from pathlib import Path
import numpy as np
R=Path(__file__).resolve().parent
sys.path.insert(0,str(R/'deps_chem_contact_v1'))
from skimage.measure import marching_cubes
from scipy.ndimage import distance_transform_edt,map_coordinates
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
from rdkit import Chem,RDConfig,rdBase
from rdkit.Chem import rdMolDescriptors,ChemicalFeatures,Crippen
import polymer_exposure_v2 as old
B=R.parent
AX=np.arange(-18,19,dtype=float)
GRID=np.stack(np.meshgrid(AX,AX,AX,indexing='ij'),-1).reshape(-1,3)
RAD=old.RAD
CASES={'PET':['SC1345656','SC1989827','SC3000704','M450','M329','M318'],
       'PA6':['SC4399876','SC4415918','SC4927359','N2240','N1626','N2403'],
       'PA66':['SC4874583','SC4927359','SC4309501','N0564','N2240','N2303']}
def fragment_props(mol):
 c=np.array(rdMolDescriptors._CalcCrippenContribs(mol))[:,0]
 for a in mol.GetAtoms():
  if a.GetAtomicNum()==1:
   ns=list(a.GetNeighbors());assert len(ns)==1
   c[ns[0].GetIdx()]+=c[a.GetIdx()];c[a.GetIdx()]=0
 factory=ChemicalFeatures.BuildFeatureFactory(str(Path(RDConfig.RDDataDir)/'BaseFeatures.fdef'))
 hba=np.zeros(len(c));hbd=hba.copy()
 for f in factory.GetFeaturesForMol(mol):
  if f.GetFamily()=='Acceptor':hba[list(f.GetAtomIds())]=1
  if f.GetFamily()=='Donor':hbd[list(f.GetAtomIds())]=1
 assert abs(c.sum()-Crippen.MolLogP(mol))<1e-6
 return dict(logp=c,hba=hba,hbd=hbd)
def sdf_for(mat):
 return B/'pet_dp10_400chain_direct_v1/00_chain/PET_DP10.sdf' if mat=='PET' else B/'unbiased_100chain_v1/00_matched_chains'/f'{mat}_DP{16 if mat=="PA6" else 8}.sdf'
class ChemEnv:
 def __init__(self,mat):
  z=np.load(R/'inputs_v2'/f'{mat}_parent_atoms.npz');self.data=z
  self.env=old.Env(z['xyz_A'],z['elements'],z['box_A'])
  mol=Chem.SDMolSupplier(str(sdf_for(mat)),removeHs=False)[0];n=mol.GetNumAtoms()
  assert len(z['elements'])==n*400
  assert [a.GetSymbol() for a in mol.GetAtoms()]==z['elements'][:n].tolist()
  meta=next(x for x in json.load(open(R/'inputs_v2/summary.json')) if x['material']==mat)
  section='';bonds=set()
  for line in Path(meta['itp_path']).read_text().splitlines():
   line=line.split(';')[0].strip()
   if line.startswith('['):section=line.strip('[] ').strip();continue
   if section=='bonds' and line:
    w=line.split();bonds.add(tuple(sorted((int(w[0])-1,int(w[1])-1))))
  assert bonds=={tuple(sorted((x.GetBeginAtomIdx(),x.GetEndAtomIdx()))) for x in mol.GetBonds()}
  props=fragment_props(mol);heavy=np.where(z['elements']!='H')[0]
  self.hxyz=z['xyz_A'][heavy];self.hr=np.array([RAD[e] for e in z['elements'][heavy]])
  self.props=np.column_stack([np.tile(props[k],400)[heavy] for k in ['logp','hba','hbd']])
  self.tree=cKDTree(self.env.wrap(self.hxyz),boxsize=[z['box_A'][0],z['box_A'][1],100000.])
  self.audit=dict(sdf=str(sdf_for(mat)),sdf_sha256=old.sha(sdf_for(mat)),itp_sha256=old.sha(meta['itp_path']),graph_atom_mapping='PASS',atom_count=n,logP_fragment_sum=float(props['logp'].sum()),rdkit_version=rdBase.rdkitVersion)
 def chemistry(self,p):
  d,ix=self.tree.query(self.env.wrap(p),k=8,workers=2)
  gap=d-self.hr[ix];w=np.exp(-.5*((gap-gap.min(1)[:,None])/1.5)**2);w/=w.sum(1)[:,None]
  return np.sum(self.props[ix]*w[:,:,None],axis=1)
def interp(a,p,g,env):
 q=np.array(p,copy=True);q[:,:2]%=env.box[:2];c=((q-g['origin'])/g['step']).T
 outside=(c[2]<0)|(c[2]>a.shape[2]-1)
 c[2]=np.clip(c[2],0,a.shape[2]-1)
 v=map_coordinates(a,c,order=1,mode='grid-wrap');v[outside]=-3.
 return v
def ses_grid(ce,h=.5):
 g=old.build_grid(ce.env,ce.data['masses_Da'],h)
 del g['gradient'];del g['clearance']
 pad=int(np.ceil(2.5/min(g['step'][:2])))
 ext=np.pad(g['external'],((pad,pad),(pad,pad),(0,0)),mode='wrap')
 dist=distance_transform_edt(~ext,sampling=g['step'])
 phi=(dist[pad:-pad,pad:-pad,:]-.5*min(g['step'])-1.4).astype(np.float32)
 del dist,ext
 g['phi']=phi
 return g
def mesh_patch(field,radius=15.,enzyme=False):
 field=np.asarray(field,np.float32)
 if field.min()>=0 or field.max()<=0:return np.empty((0,3)),np.empty((0,3),int),np.empty((0,3))
 v,f,_,_=marching_cubes(field,level=0,spacing=(1,1,1),allow_degenerate=False)
 v+=AX[0]
 grad=np.stack([map_coordinates(a,(v-AX[0]).T,order=1,mode='nearest') for a in np.gradient(field)],axis=1)
 normals=-grad/np.maximum(np.linalg.norm(grad,axis=1)[:,None],1e-10)
 valid=np.all(np.linalg.norm(v[f],axis=2)<=radius,axis=1)
 if enzyme:valid&=normals[f].mean(1)[:,2]<.2
 f=f[valid]
 if len(f)==0:return np.empty((0,3)),np.empty((0,3),int),np.empty((0,3))
 ii=np.r_[f[:,0],f[:,1],f[:,2]];jj=np.r_[f[:,1],f[:,2],f[:,0]]
 _,lab=connected_components(coo_matrix((np.ones(len(ii)),(ii,jj)),shape=(len(v),len(v))),directed=False)
 present=np.unique(lab[f[:,0]]);best=None
 for k in present:
  ff=f[lab[f[:,0]]==k];vid=np.unique(ff)
  area=np.linalg.norm(np.cross(v[ff[:,1]]-v[ff[:,0]],v[ff[:,2]]-v[ff[:,0]]),axis=1).sum()/2
  metric=(float(np.linalg.norm(v[vid],axis=1).min()),-float(area))
  if best is None or metric<best[0]:best=(metric,ff)
 f=best[1];used=np.unique(f);lookup=np.full(len(v),-1);lookup[used]=np.arange(len(used))
 return v[used],lookup[f],normals[used]
def core_clear(p,q,radii):
 out=np.full(len(p),np.inf)
 for a,r in zip(q,radii):out=np.minimum(out,np.linalg.norm(p-a,axis=1)-r)
 return out
def protein_field(material_phi,points,q,radii,gap):
 return np.maximum(-np.asarray(material_phi)-gap,-core_clear(points,q,radii))
def lipo_colors(v):
 t=np.clip(np.asarray(v)/.5,-1,1)
 neutral=np.array([.85,.85,.82]);polar=np.array([.12,.40,.83]);nonpolar=np.array([.78,.48,.10])
 return np.where((t>=0)[:,None],neutral+(nonpolar-neutral)*t[:,None],neutral+(polar-neutral)*(-t[:,None]))
def contact_colors(props,supported):
 out=lipo_colors(props[:,0]);accept=props[:,1]>.25;donate=props[:,2]>.25
 out[accept]=[.22,.65,.90]
 out[donate]=[.70,.35,.75]
 out[accept&donate]=[.34,.72,.50]
 out[~supported]=[.55,.55,.55]
 return out
def vertex_areas(v,f):
 area=np.linalg.norm(np.cross(v[f[:,1]]-v[f[:,0]],v[f[:,2]]-v[f[:,0]]),axis=1)/2
 va=np.zeros(len(v))
 for k in range(3):np.add.at(va,f[:,k],area/3)
 return va
def chem_on_mesh(v,sites,frame_getter,ce,g):
 total=np.zeros((len(v),3));sq=total.copy();den=np.zeros(len(v));hits=np.zeros(len(v),int)
 for s in sites:
  F=frame_getter(s);p=np.array(s['center'])+v@F.T
  phi=interp(g['phi'],p,g,ce.env)
  w=np.exp(-.5*(phi/1.5)**2);w[np.abs(phi)>4.5]=0
  val=ce.chemistry(p)
  total+=val*w[:,None];sq+=val*val*w[:,None];den+=w;hits+=abs(phi)<=1.5
 mean=total/np.maximum(den[:,None],1e-12)
 sd=np.sqrt(np.maximum(sq/np.maximum(den[:,None],1e-12)-mean*mean,0))
 support=hits/len(sites);valid=(support>=.2)&(den>.01)
 return mean,sd,support,valid
def export_mesh(folder,name,v,f,n,props,sd,support,valid,enzyme=False,fixedmask=None):
 folder.mkdir(exist_ok=True)
 colors=contact_colors(props,valid) if enzyme else lipo_colors(props[:,0])
 if not enzyme:colors[~valid]=[.55,.55,.55]
 if fixedmask is not None:colors[fixedmask]=[.92,.73,.22]
 np.savez_compressed(folder/(name+'.npz'),vertices_A=v,faces=f,normals=n,fragment_lipophilicity_index=props[:,0],material_HBA_weight=props[:,1],material_HBD_weight=props[:,2],chemical_sd=sd,surface_support_fraction=support,chemical_supported=valid,colors_rgb=colors)
 with (folder/(name+'.ply')).open('x') as z:
  z.write(f'ply\nformat ascii 1.0\nelement vertex {len(v)}\nproperty float x\nproperty float y\nproperty float z\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nelement face {len(f)}\nproperty list uchar int vertex_indices\nend_header\n')
  for a,c in zip(v,np.rint(colors*255).astype(int)):z.write(' '.join(map(str,[*np.round(a,4),*c]))+'\n')
  for a in f:z.write('3 '+' '.join(map(str,a))+'\n')
 va=vertex_areas(v,f);den=va.sum()
 stats=dict(vertices=len(v),triangles=len(f),area_A2=float(den),chemically_supported_area_fraction=float(va[valid].sum()/den),mean_fragment_lipophilicity_index=float(np.sum(props[:,0]*va)/den),positive_index_area_fraction=float(va[(props[:,0]>0)&valid].sum()/den),negative_index_area_fraction=float(va[(props[:,0]<0)&valid].sum()/den),material_HBA_weighted_area_fraction=float(np.sum(props[:,1]*va)/den),material_HBD_weighted_area_fraction=float(np.sum(props[:,2]*va)/den))
 payload=dict(name=name,v=np.round(v,4).tolist(),f=f.tolist(),n=np.round(n,5).tolist(),c=np.round(colors,4).tolist())
 return stats,payload
def viewer(path,meshes,core=None,default=None):
 payload=dict(meshes=meshes,core=core,default=default)
 blob=base64.b64encode(zlib.compress(json.dumps(payload,separators=(',',':')).encode(),9)).decode()
 script="# Self-contained PyMOL viewer; geometry/chemical preference model, not a folded protein.\n"
 script+="import json,zlib,base64\nfrom pymol import cmd\nfrom pymol.cgo import BEGIN,END,TRIANGLES,NORMAL,VERTEX,COLOR,SPHERE,CYLINDER\n"
 script+="D=json.loads(zlib.decompress(base64.b64decode("+repr(blob)+")))\n"
 script+="""
for mesh in D['meshes']:
 obj=[BEGIN,TRIANGLES]
 for face in mesh['f']:
  for i in face:obj += [COLOR]+mesh['c'][i]+[NORMAL]+mesh['n'][i]+[VERTEX]+mesh['v'][i]
 obj += [END];cmd.load_cgo(obj,mesh['name'])
 if D['default'] and mesh['name']!=D['default']:cmd.disable(mesh['name'])
cmd.pseudoatom('target_carbonyl_C',pos=[0,0,0],color='cyan',label='carbonyl C')
cmd.show('spheres','target_carbonyl_C');cmd.set('sphere_scale',.25,'target_carbonyl_C')
if D['core']:
 c=D['core'];q=c['q'];obj=[];cols={'C':[.25,.7,.25],'N':[.2,.3,.9],'O':[.9,.2,.2],'S':[.9,.8,.2]}
 for a,e in zip(q,c['elements']):obj += [COLOR]+cols.get(e,[.6,.6,.6])+[SPHERE]+a+[.22]
 for i,j in c['bonds']:
  obj += [CYLINDER]+q[i]+q[j]+[.10]+cols.get(c['elements'][i],[.6,.6,.6])+cols.get(c['elements'][j],[.6,.6,.6])
 cmd.load_cgo(obj,'fixed_catalytic_core')
 for i,a in enumerate(c['donors']):
  cmd.pseudoatom('donor_'+str(i+1),pos=a,color='yellow',label='D'+str(i+1))
  cmd.show('spheres','donor_'+str(i+1));cmd.set('sphere_scale',.3,'donor_'+str(i+1))
cmd.bg_color('white');cmd.set('two_sided_lighting',1);cmd.set('cgo_transparency',.12)
cmd.zoom('target_carbonyl_C',18)
print('Material lipophilicity: blue negative/polar-fragment index, brown positive/lipophilic. Gray insufficient support.')
print('Enzyme contact preference: cyan donor-facing material acceptor; purple acceptor-facing material donor; green mixed; gold fixed core.')
print('This is a local contact-face hypothesis, not a protein structure, free energy or activity prediction.')
"""
 with path.open('x') as f:f.write(script)
def core_payload(q,elements,donors,template,kind):
 bounds=[0,7,15,23] if template=='NylC' else [0,6,14,24]
 if kind=='two_residue_donors':bounds+=[bounds[-1]+6,bounds[-1]+12]
 else:bounds+=[len(q)]
 group=np.zeros(len(q),int)
 for j in range(len(bounds)-1):group[bounds[j]:bounds[j+1]]=j
 cov={'C':.76,'N':.71,'O':.66,'S':1.05};bonds=[]
 for i in range(len(q)):
  for j in range(i):
   dist=np.linalg.norm(q[i]-q[j])
   if group[i]==group[j] and .8<dist<cov[elements[i]]+cov[elements[j]]+.3:bonds.append([i,j])
 return dict(q=q.tolist(),elements=list(elements),donors=donors.tolist(),bonds=bonds,bond_scope='geometry inferred within known residue fragments; original core coordinates unchanged')
def material(mat,all_sites,panel,coverage,cores,out,pilot):
 ce=ChemEnv(mat);old.emit('SES_grid_start',material=mat)
 g=ses_grid(ce)
 old.emit('SES_grid_ready',material=mat,shape=g['phi'].shape)
 d=out/mat;d.mkdir();(d/'chemical_typing_audit.json').write_text(json.dumps(ce.audit,indent=2))
 oldrows=json.load(open(R/'polymer_exposure_full_v2'/mat/'sites.json'))
 bins={r['site_id']:r['fine_bin'] for r in json.load(open(R/'polymer_exposure_fine_strata_v2'/mat/'site_bin_assignments.json'))}
 items=[dict(s,bin=bins[s['site_id']]) for s in oldrows]
 if pilot:items=items[:3]
 F=np.empty((len(items),37,37,37),np.float16)
 for i,s in enumerate(items):
  world=np.array(s['center'])+GRID@np.array(s['surface_frame']).T
  F[i]=interp(g['phi'],world,g,ce.env).reshape(37,37,37)
  if i%100==0:old.emit('material_fields',material=mat,done=i+1,total=len(items))
 np.savez_compressed(d/'polymer_SES_fields.npz',fields=F,axis_A=AX)
 (d/'site_index.json').write_text(json.dumps([dict(site_id=s['site_id'],bin=s['bin'],side=s['side'],normal_status=s['normal_status']) for s in items],indent=2))
 atlas=[];ad=d/'material_atlas';ad.mkdir()
 for group,side in sorted(set((s['bin'],s['side']) for s in items)):
  if 'UNCERTAIN' in group or 'UNRESOLVED' in group:continue
  ix=[i for i,s in enumerate(items) if s['bin']==group and s['side']==side and s['normal_status']=='PASS']
  if not ix:continue
  phi=np.median(F[ix].astype(np.float32),axis=0)
  v,f,n=mesh_patch(phi)
  if not len(f):atlas.append(dict(group=group,side=side,status='NO_SURFACE'));continue
  ss=[items[i] for i in ix]
  props,sd,support,valid=chem_on_mesh(v,ss,lambda s:np.array(s['surface_frame']),ce,g)
  key=mat+'_'+group+'_'+side
  stats,pay=export_mesh(ad,key,v,f,n,props,sd,support,valid)
  viewer(ad/(key+'_view.py'),[pay])
  rec=dict(group=group,side=side,n=len(ix),status='DESCRIPTIVE_LOW_N' if len(ix)<5 else 'EMPIRICAL_CHEMICAL_SURFACE',**stats);atlas.append(rec)
 (d/'material_atlas_summary.json').write_text(json.dumps(atlas,indent=2))
 cr=[];selected=CASES[mat][:1] if pilot else CASES[mat]
 for cid in selected:
  pi=next(i for i,p in enumerate(panel) if p['id']==cid);p=panel[pi]
  ids=[i for i,s in enumerate(all_sites) if s['surface']==mat+'_original' and coverage[pi,i]]
  ss=[all_sites[i] for i in ids];assert ss
  q=cores[cid+'_q'];donors=cores[cid+'_donors'];nt=23 if p['template']=='NylC' else 24
  rot=old.make_frame(q[:nt].mean(0),donors[1]-donors[0],np.array([1.,0.,0.]))
  qc=q@rot;dc=donors@rot;radii=np.array([RAD[x] for x in p['elements']])
  fields=[]
  for s in ss:
   world=np.array(s['center'])+GRID@(np.array(s['frame'])@rot).T
   fields.append(interp(g['phi'],world,g,ce.env).reshape(37,37,37))
  fields=np.array(fields,np.float32)
  cd=d/cid;cd.mkdir()
  np.savez_compressed(cd/'matched_material_fields.npz',fields=fields,site_indices=ids,core_q_A=qc,donors_A=dc,chemical_to_core_frame=rot)
  meshes=[];rec=dict(core_id=cid,template=p['template'],donor_kind=p['kind'],coverage_n=len(ids),denominator=sum(s['surface']==mat+'_original' for s in all_sites),sites=[s['site_id'] for s in ss],interfaces=[])
  phi50=np.median(fields,axis=0);v,f,n=mesh_patch(phi50)
  if len(f):
   props,sd,sup,valid=chem_on_mesh(v,ss,lambda s:np.array(s['frame'])@rot,ce,g)
   st,pay=export_mesh(cd,'matched_material_median',v,f,n,props,sd,sup,valid);meshes.append(pay)
   rec['matched_material_mesh']=st
  for quantile in [.5,.9,1.]:
   phi=np.quantile(fields,quantile,axis=0)
   pf=protein_field(phi.ravel(),GRID,qc,radii,.5).reshape(37,37,37)
   v,f,n=mesh_patch(pf,enzyme=True)
   if not len(f):
    rec['interfaces'].append(dict(quantile=quantile,status='NO_CONTACT_PATCH'));continue
   props,sd,sup,valid=chem_on_mesh(v,ss,lambda s:np.array(s['frame'])@rot,ce,g)
   fixed=np.abs(core_clear(v,qc,radii))<.5
   st,pay=export_mesh(cd,'enzyme_contact_q'+str(int(quantile*100)),v,f,n,props,sd,sup,valid,enzyme=True,fixedmask=fixed);meshes.append(pay)
   area=vertex_areas(v,f);clashes=[];clashes04=[]
   for s in ss:
    world=np.array(s['center'])+v@(np.array(s['frame'])@rot).T
    z=interp(g['phi'],world,g,ce.env)
    clashes.append(float(area[z>0].sum()/area.sum()));clashes04.append(float(area[z>.4].sum()/area.sum()))
   inrange=np.max(abs(qc),axis=1)<17
   qfield=map_coordinates(pf,(qc[inrange]-AX[0]).T,order=1,mode='nearest')
   rec['interfaces'].append(dict(quantile=quantile,status='CANDIDATE_CONTACT_FACE_NOT_FULL_PROTEIN',core_centers_in_local_volume=bool(np.all(qfield>=-1e-4)),core_centers_checked=int(inrange.sum()),mean_material_penetrating_area_fraction=float(np.mean(clashes)),worst_material_penetrating_area_fraction=float(max(clashes)),mean_penetrating_area_fraction_gt04A=float(np.mean(clashes04)),worst_penetrating_area_fraction_gt04A=float(max(clashes04)),sites_with_zero_penetration_gt04A=sum(x==0 for x in clashes04),per_site_penetration=clashes,per_site_penetration_gt04A=clashes04,**st))
  viewer(cd/(mat+'_'+cid+'_view.py'),meshes,core_payload(qc,p['elements'],dc,p['template'],p['kind']),'enzyme_contact_q90')
  (cd/'summary.json').write_text(json.dumps(rec,indent=2));cr.append(rec)
  old.emit('core_contact_complete',material=mat,core=cid,sites=len(ids),interfaces=len(rec['interfaces']))
 result=dict(material=mat,total_sites=len(items),chemical_typing=ce.audit,material_atlas=atlas,core_conditioned=cr)
 (d/'summary.json').write_text(json.dumps(result,indent=2))
 return result
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);ap.add_argument('--pilot',action='store_true');a=ap.parse_args()
 out=R/a.output;assert out.parent==R;out.mkdir(exist_ok=False)
 start=time.time();log=R/'RUN_LOG.jsonl'
 with log.open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now().astimezone().isoformat(),event='chemical_contact_start',command=sys.argv,output=str(out),source_sha256=old.sha(__file__)))+'\n')
 panel=json.load(open(R/'panel_surface_scan_v1/panel.json'));sites=json.load(open(R/'panel_surface_scan_v1/sites.json'));cov=np.load(R/'panel_surface_scan_v1/coverage.npz')['coverage'];cores=np.load(R/'panel_surface_scan_v1/cores.npz')
 results=[material(mat,sites,panel,cov,cores,out,a.pilot) for mat in (['PET'] if a.pilot else ['PET','PA6','PA66'])]
 result=dict(status='TECHNICAL_COMPLETE_CHEMICAL_PROXY_AND_CONTACT_HYPOTHESES',seconds=time.time()-start,materials=results,source_sha256=old.sha(__file__),definitions=dict(material_surface='Approximate SES from exterior-water EDT at0.5A; local mesh1A; open R15 patch',lipophilicity='Wildman-Crippen atomic logP contributions, explicit-H contributions aggregated to heavy atoms; 8-neighbor1.5A relative-distance smoothing; dimensionless empirical fragment index, not hydration free energy',material_color='brown positive fragment index, blue negative; gray low near-surface support',enzyme_surface='local candidate complement of quantile polymer SES with0.5A stand-off, union with fixed-core vdW volumes; material-facing connected component only; NOT full protein',enzyme_chemistry='hypothesized contact preference; HBA->HBD, HBD->HBA; nonpolar->nonpolar; no charge sign inversion; gold fixed-core area',quantile='local-field quantile, NOT fraction of entire sites jointly accommodated',limitations=['single snapshot','inherited exposure bins and numerical uncertainty','approximate SES and mesh resolution','no protein backbone/foldability evidence','no energetic validation','surface chemistry is empirical lipophilicity plus donor/acceptor annotation']))
 (out/'summary.json').write_text(json.dumps(result,indent=2))
 with log.open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now().astimezone().isoformat(),event='chemical_contact_complete',output=str(out),exit_status=0,seconds=result['seconds']))+'\n')
 with (R/'RUNBOOK.md').open('a') as f:f.write('\nChemical surface/contact candidates: '+str(out)+'. Reproduce: CPU Python chemical_contact_surface_v1.py --output NEW_NAME. Dependencies isolated in deps_chem_contact_v1. Read definitions and per-core penetration audits; not a complete enzyme. Portable *_view.py files require PyMOL only.\n')
 old.emit('chemical_contact_complete',output=str(out),seconds=result['seconds'])
if __name__=='__main__':main()
