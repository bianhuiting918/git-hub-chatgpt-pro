#!/usr/bin/env python3
"""Export selected fixed cores and matched density surfaces into common TDD frame."""
import tdd_donor_pilot_v1 as p
import numpy as np,json,gzip,datetime
ROOT=p.ROOT
def pdb_line(i,name,res,chain,num,q,elem):
 return f"ATOM  {i:5d} {name:>4s} {res:>3s} {chain}{num:4d}    {q[0]:8.3f}{q[1]:8.3f}{q[2]:8.3f}  1.00  0.00          {elem:>2s}\n"
if __name__=='__main__':
 out=ROOT/'pymol_export_v1';out.mkdir(exist_ok=False);rc=1
 try:
  report=json.loads((ROOT/'conditional_density_v1/summary.json').read_text())
  ranks=json.loads((ROOT/'fixed_core_transfer_v1/rankings.json').read_text())
  meta=json.loads((ROOT/'fixed_core_library_v1/motifs.json').read_text())
  sites=json.loads((ROOT/'fixed_core_transfer_v1/sites.json').read_text())
  z=np.load(ROOT/'fixed_core_transfer_v1/coverage.npz');cores=z['cores_A'];H=z['donor_H_A'];cov=z['coverage'];pcov=z['parent_coverage']
  chosen=['N0108']
  for g in report['groups'].values():
   if g['best']['id'] not in chosen:chosen.append(g['best']['id'])
  maximum=max(r['coverage_total'] for r in ranks)
  for r in ranks:
   if r['coverage_total']==maximum and r['id'] not in chosen:chosen.append(r['id'])
  mesh={};patchids={};chem={}
  for mat in ('PA6','PA66'):
   patchids[mat]=np.load(ROOT/'patch_shape_chemistry_v1'/f'{mat}_patch_face_indices.npz')
   for side in ('top','bottom'):
    mesh[mat,side]=np.load(ROOT/'surface_density_v2'/f'{mat}_{side}_mesh.npz')
    chem[mat,side]=np.load(ROOT/'patch_shape_chemistry_v1'/f'{mat}_{side}_chemistry.npz')
  models=[];patches={}
  for mid in chosen:
   i=next(i for i,m in enumerate(meta) if m['id']==mid);q=cores[i]
   a=q[:23];b=p.tq;ac=a.mean(0);bc=b.mean(0);u,sv,vt=np.linalg.svd((a-ac).T@(b-bc))
   d=np.eye(3);d[-1,-1]=np.linalg.det(u@vt);R=u@d@vt;t=bc-ac@R
   aligned=q@R+t;alignedH=H[i]@R+t
   assert np.max(abs(aligned[:23]-b))<1e-8
   lines=['REMARK Artificial dual backbone NH design. NylC TDD rigid. Virtual H. Not native NylC donors.\n']
   for j,(v,l) in enumerate(zip(aligned[:23],p.template[4]),1):
    lines.append(pdb_line(j,l[12:16].strip(),l[17:20],'A',int(l[22:26]),v,l[76:78].strip()))
   for donor,chain in enumerate(('B','C')):
    for k,l in enumerate(p.donor_templates[0][3]):
     j=23+donor*6+k;lines.append(pdb_line(j+1,l[12:16].strip(),l[17:20],chain,int(l[22:26]),aligned[j],l[76:78].strip()))
    lines.append(pdb_line(36+donor,'HN','TYR',chain,87,alignedH[donor],'H'))
   lines.append('END\n');(out/f'{mid}_core.pdb').write_text(''.join(lines))
   matched=np.flatnonzero(cov[i]);matched=sorted(matched,key=lambda k:(-sites[k]['support'],sites[k]['material'],sites[k]['site_id']))
   entries=[]
   for k in matched:
    site=sites[k];mat=site['material'];side=site['side'];pid=mat+'_'+site['site_id']
    if pid not in patches:
     mm=mesh[mat,side];ids=patchids[mat][site['site_id']]
     tri=mm['vertices_A'][mm['faces'][ids]]
     F=np.array(site['frame']);c=np.array(site['origin_A']);tri=(tri-c)@F
     cc=chem[mat,side]
     patches[pid]=dict(material=mat,site_id=site['site_id'],side=side,support=site['support'],triangles_amide_frame_A=tri.tolist(),chemistry_indices=cc['class_indices'][ids].tolist(),chemistry_names=cc['class_names'].tolist(),water_direction_amide=(np.array([0,0,1 if side=='top' else -1])@F).tolist())
    entries.append(dict(patch_id=pid,parent_clear=bool(pcov[i,k])))
   models.append(dict(id=mid,source=meta[i],R=R.tolist(),t=t.tolist(),sites=entries,coverage_total=len(entries),coverage_PA6=sum(patches[e['patch_id']]['material']=='PA6' for e in entries),coverage_PA66=sum(patches[e['patch_id']]['material']=='PA66' for e in entries),parent_coverage=int(pcov[i].sum())))
  data=dict(models=models,patches=patches,default='N0108',scope='finite v10 coarse density surface, connected 15 A patches; core is fixed TDD plus artificial backbone NH donors',chemistry='nearest heavy atom VDW class, not hydrophobic potential',rank_groups=report['groups'])
  with gzip.open(out/'overlay_data.json.gz','wt') as f:json.dump(data,f,separators=(',',':'))
  summary=dict(selected_cores=chosen,unique_patches=len(patches),surface_pairings=sum(m['coverage_total'] for m in models),all_library_parent_coverage_edges=int(pcov.sum()),all_library_parent_covered_sites=int(pcov.any(0).sum()),models=[{k:v for k,v in m.items() if k not in ('sites','R','t','source')} for m in models])
  (out/'summary.json').write_text(json.dumps(summary,indent=2))
  hashes={f.name:p.sha(f) for f in out.iterdir() if f.is_file()}
  (out/'SHA256.json').write_text(json.dumps(hashes,indent=2));print(json.dumps(summary),flush=True);rc=0
 finally:
  with (ROOT/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(p.Path(__file__).resolve()),script_sha256=p.sha(p.Path(__file__)),command='CPU Python pymol_export_v1.py',inputs='fixed cores, conditional ranks, matched connected meshes',outputs=str(out),parameters=dict(alignment='proper rigid fit all 23 TDD atoms'),exit_code=rc))+'\n')
