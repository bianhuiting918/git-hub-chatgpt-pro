#!/usr/bin/env python3
"""Audit polymer charges and identify original nylon bulk-water interface sites."""
from pathlib import Path
import json,re,sys,hashlib,datetime
import numpy as np
from scipy.ndimage import gaussian_filter,label,map_coordinates
from scipy.spatial import cKDTree
R=Path(__file__).resolve().parent;B=R.parent;N=B/'rectangular_400chain_v1/nylc_tdd_surface_v1'
sys.path.insert(0,str(N.parent/'surface_previews/axis_density_sixface_all_external_trim_v6'))
from density_cut_variants import parse_pdb
from terminal_density_trim_3d import center_whole_chains_in_primary_cell
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def itp_atoms(p):
 section='';out=[]
 for line in p.read_text().splitlines():
  line=line.split(';')[0].strip()
  if line.startswith('['):section=line.strip('[] ').strip();continue
  if section=='atoms' and line:
   f=line.split();out.append(dict(index=int(f[0]),name=f[4],charge=float(f[6]),mass=float(f[7]),atomtype=f[1]))
 assert [a['index'] for a in out]==list(range(1,len(out)+1))
 return out
def amides(sdf):
 lines=sdf.read_text().splitlines();na=int(lines[3][:3]);nb=int(lines[3][3:6])
 el=[x[31:34].strip() for x in lines[4:4+na]];adj=[[] for _ in el]
 for l in lines[4+na:4+na+nb]:
  a,b,o=int(l[:3])-1,int(l[3:6])-1,int(l[6:9]);adj[a].append((b,o));adj[b].append((a,o))
 motifs=[]
 for c in range(na):
  if el[c]!='C':continue
  oo=[a for a,o in adj[c] if el[a]=='O' and o==2];nn=[a for a,o in adj[c] if el[a]=='N' and o==1]
  if oo:
   assert len(oo)==1
   motifs.extend((c,oo[0],n) for n in nn)
 return motifs,na
if __name__=='__main__':
 O=R/'inputs_v1';O.mkdir(exist_ok=False);report=[];rc=1
 try:
  for mat,dp in [('PET',10),('PA6',16),('PA66',8)]:
   water=(B/'pet_dp10_400chain_direct_v1/pet_dp10_400chain/water_slab_v1' if mat=='PET' else N.parent/f'{mat.lower()}_dp{dp}_400chain/water_slab_v1')
   pdb=water/'dry_export'/f'{mat}_DP{dp}_400chain_water_equilibrated_dry.pdb'
   top=water/'system_water.top';inc=[Path(p) for p in re.findall(r'#include\s+"([^"]+)"',top.read_text()) if p.endswith('_GMX.itp')]
   assert len(inc)==1;aa=itp_atoms(inc[0])
   _,box,atoms,_=parse_pdb(pdb);box=np.array(box)
   if mat!='PET':center_whole_chains_in_primary_cell(atoms,box)
   assert len(atoms)==400*len(aa)
   for i,a in enumerate(atoms):assert a['line'][12:16].strip()==aa[i%len(aa)]['name'],(mat,i,a['line'][12:16],aa[i%len(aa)])
   q=np.tile([a['charge'] for a in aa],400);mass=np.tile([a['mass'] for a in aa],400)
   xyz=np.array([a['xyz'] for a in atoms]);els=np.array([a['element'] for a in atoms]);serial=np.array([a['serial'] for a in atoms])
   np.savez_compressed(O/(mat+'_parent_atoms.npz'),xyz_A=xyz,elements=els,charges_e=q,masses_Da=mass,serials=serial,box_A=box,atom_names=np.array([a['line'][12:16].strip() for a in atoms]))
   rec=dict(material=mat,parent_atoms=len(atoms),atoms_per_chain=len(aa),parent_total_charge_e=float(q.sum()),single_chain_charge_e=float(q[:len(aa)].sum()),parent_sha256=sha(pdb),itp_path=str(inc[0]),itp_sha256=sha(inc[0]),atom_order_name_mapping='PASS_ALL_ATOMS',charge_scope='fixed AM1-BCC force-field partial charges; not surface potential, not solvent response')
   if mat!='PET':
    trimmed=N/'inputs'/(mat+'_v10_ge90_uncapped.pdb')
    lines=[l for l in trimmed.read_text().splitlines() if l.startswith(('ATOM  ','HETATM'))]
    lookup={v:i for i,v in enumerate(serial)};ids=np.array([lookup[l[6:11]] for l in lines])
    tx=np.array([[float(l[k:k+8]) for k in (30,38,46)] for l in lines]);assert np.max(np.abs(tx-xyz[ids]))<.002
    np.savez_compressed(O/(mat+'_trimmed_atoms.npz'),xyz_A=tx,elements=els[ids],inherited_charges_e=q[ids],parent_indices=ids,serials=serial[ids],box_A=box)
    rec.update(trimmed_atoms=len(ids),trimmed_inherited_total_charge_e=float(q[ids].sum()),trimmed_charge_status='NOT_VALIDATED_UNCAPPED_CHEMISTRY_NO_ESP_INTERPRETATION')
    # Original parent boundaries, not the previous interior positions exposed by trimming.
    x=xyz.copy();x[:,:2]%=box[:2];n=np.ceil(box).astype(int);step=box/n
    edges=[np.linspace(0,box[k],n[k]+1) for k in range(3)]
    hist=np.histogramdd(x,bins=edges,weights=mass)[0]
    assert abs(hist.sum()-mass.sum())<1e-4
    rho=gaussian_filter(hist/np.prod(step),3/step,mode=('wrap','wrap','constant'))
    zc=(edges[2][1:]+edges[2][:-1])/2
    threshold=.5*float(np.median(rho[:,:,abs(zc-np.median(x[:,2]))<10]))
    labs,num=label(rho>=threshold);sz=np.bincount(labs.ravel());sz[0]=0;mask=labs==sz.argmax();valid=mask.any(2)
    heights={}
    for side,sign in [('top',1),('bottom',-1)]:
     ix=np.where(valid,np.max(np.where(mask,np.arange(n[2]),-1),axis=2) if sign==1 else np.min(np.where(mask,np.arange(n[2]),n[2]),axis=2),0)
     adj=np.clip(ix+sign,0,n[2]-1)
     v0=np.take_along_axis(rho,ix[:,:,None],2)[:,:,0];v1=np.take_along_axis(rho,adj[:,:,None],2)[:,:,0]
     frac=np.clip((threshold-v0)/np.where(abs(v1-v0)>1e-12,v1-v0,1),0,1)
     hh=zc[ix]+frac*(zc[adj]-zc[ix]);hh[~valid]=np.nan;heights[side]=hh
    motifs,na=amides(B/'unbiased_100chain_v1/00_matched_chains'/f'{mat}_DP{dp}.sdf')
    assert na==len(aa) and len(motifs)==17
    disk=np.array([(dx,dy) for dx in np.linspace(-10,10,11) for dy in np.linspace(-10,10,11) if dx*dx+dy*dy<=100])
    assert len(disk)==81
    rows=[]
    for chain in range(400):
     for mi,motif in enumerate(motifs):
      ii=np.array(motif)+chain*na;c=x[ii[0]]
      offsets={side:float(sign*(c[2]-map_coordinates(hh,(c[:2]/step[:2]-.5)[:,None],order=1,mode='grid-wrap')[0])) for (side,hh),sign in zip(heights.items(),[1,-1])}
      side=min(offsets,key=lambda side:abs(offsets[side]));sign=1 if side=='top' else -1
      if not np.isfinite(offsets[side]) or abs(offsets[side])>3:continue
      points=np.column_stack((c[0]+disk[:,0],c[1]+disk[:,1],np.full(len(disk),c[2]-sign*5)))
      assert (points[:,2]>=0).all() and (points[:,2]<box[2]).all()
      support=float((map_coordinates(rho,(points/step-.5).T,order=1,mode='grid-wrap')>=threshold).mean())
      if support<.7:continue
      a=x[ii].copy();a[:,:2]-=np.round((a[:,:2]-c[:2])/box[:2])*box[:2]
      rows.append(dict(site_id=f'P{chain+1:03d}_{aa[motif[0]]["name"]}_{aa[motif[1]]["name"]}',original_chain=chain+1,amide_index=mi+1,parent_indices=ii.tolist(),serials=serial[ii].tolist(),carbonyl_xyz_A=c.tolist(),amide_A=a.tolist(),side=side,support=support,boundary_offset_A=offsets[side]))
    (O/(mat+'_original_surface_sites.json')).write_text(json.dumps(rows,indent=2))
    rec.update(original_parent_amides=6800,original_surface_sites=len(rows),original_surface_cohorts={str(t):sum(r['support']>=t for r in rows) for t in [.9,.8,.7]},surface_status='DENSITY_SELECTED_CORE_ACCESSIBILITY_NOT_EVALUATED',surface_protocol=dict(sigma_A=3,threshold_fraction=.5,offset_A=3,support_depth_A=5,support_radius_A=10,support_points=81,XY_periodic=True))
   report.append(rec);print(json.dumps(rec),flush=True)
  (O/'summary.json').write_text(json.dumps(report,indent=2));rc=0
 finally:
  with (R/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(Path(__file__)),source_sha256=sha(Path(__file__)),output=str(O),exit_code=rc))+'\n')
