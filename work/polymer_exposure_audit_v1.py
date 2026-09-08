#!/usr/bin/env python
"""Read-only primary-data audit plus versioned refinement / shape summaries."""
import json,sys,time,hashlib
from pathlib import Path
import numpy as np
import polymer_exposure_v2 as m
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score,adjusted_rand_score
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

R=Path(__file__).resolve().parent
def heights(field):
 h=np.full((31,31),np.nan)
 for ix in range(31):
  for iy in range(31):
   v=field[ix,iy]
   crossings=np.where((v[:-1]>=.5)&(v[1:]<.5))[0]
   if len(crossings):
    k=int(crossings[-1]);h[ix,iy]=-15+k+(.5-v[k])/(v[k+1]-v[k])
 return h
def shape_descriptor(field):
 h=heights(field);a=np.arange(-15,16);x,y=np.meshgrid(a,a,indexing='ij')
 mask=(x*x+y*y<=100)&np.isfinite(h)&(abs(h)<12)
 if mask.sum()<40:return dict(status='NOT_EVALUATED_INSUFFICIENT_SINGLE_VALUED_PATCH')
 X=np.column_stack([.5*x[mask]**2,x[mask]*y[mask],.5*y[mask]**2,x[mask],y[mask],np.ones(mask.sum())])
 b=np.linalg.lstsq(X,h[mask],rcond=None)[0];e=np.linalg.eigvalsh([[b[0],b[1]],[b[1],b[2]]])
 return dict(status='QUADRATIC_HEIGHT_PROXY',valid_columns=int(mask.sum()),height_at_C_axis_A=float(b[5]),tangent_slope=float(np.linalg.norm(b[3:5])),height_Hessian_eigenvalues_invA=e.tolist(),quadratic_residual_RMS_A=float(np.sqrt(np.mean((X@b-h[mask])**2))),scope='outermost water-center boundary height fit over R10; Hessian eigenvalues are not exact principal curvatures; overhangs not represented')
def group_cluster(arr):
 n=len(arr)
 if n<20:return None,dict(status='NOT_EVALUATED_N_LT20')
 f=arr[:,::3,::3,::3].reshape(n,-1).astype(float)
 valid=f.std(0)>.05;f=f[:,valid]
 if not valid.any():return None,dict(status='NO_VARIATION')
 # No feature-wise scaling: each occupancy voxel receives equal geometric weight.
 fits=[]
 for k in [2,3,4]:
  km=KMeans(n_clusters=k,random_state=20260908,n_init=10).fit(f)
  ct=np.bincount(km.labels_)
  if ct.min()<5:continue
  ss=silhouette_score(f,km.labels_)
  fits.append((ss,k,km.labels_))
 if not fits:return None,dict(status='NO_ADMISSIBLE_CLUSTERING')
 ss,k,lab=max(fits,key=lambda x:x[0])
 # Stability at a second initialization is necessary but not independent validation.
 alt=KMeans(n_clusters=k,random_state=20260909,n_init=10).fit_predict(f)
 ari=adjusted_rand_score(lab,alt)
 rec=dict(status='EXPLORATORY_ONLY',silhouette=float(ss),k=k,initialization_ARI=float(ari))
 if ss<.15 or ari<.8:
  rec['status']='NO_CLEAR_STABLE_SEPARATION';return None,rec
 return lab,rec
def main():
 src=R/'polymer_exposure_full_v2';out=R/'polymer_exposure_audit_v1';out.mkdir()
 start=time.time();summaries=[];allrows={}
 axis=np.arange(-15,16);grid=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),-1).reshape(-1,3);sphere=np.linalg.norm(grid,axis=1)<=15
 primary=json.load(open(src/'summary.json'))
 assert [x['n'] for x in primary['materials']]==[578,413,445]
 for mat in ['PET','PA6','PA66']:
  rows=json.load(open(src/mat/'sites.json'));z=np.load(src/mat/'site_fields.npz');fields=z['fields']
  assert fields.shape==(len(rows),31,31,31)
  assert len(set(r['site_id'] for r in rows))==len(rows)
  d=out/mat;d.mkdir();allrows[mat]=rows
  desc=[dict(site_id=r['site_id'],exposure_group=r['exposure_group'],**shape_descriptor(f.astype(float)/255)) for r,f in zip(rows,fields)]
  (d/'shape_descriptors.json').write_text(json.dumps(desc,indent=2))
  clusters=[]
  for side in ['top','bottom']:
   for exposure in ['sampled_zero','low_0_5pct','mid_5_15pct','high_gt15pct']:
    ids=[i for i,r in enumerate(rows) if r['side']==side and r['exposure_group']==exposure and r['normal_status']=='PASS']
    lab,rec=group_cluster(fields[ids].astype(float)/255)
    rec.update(side=side,exposure=exposure,n=len(ids))
    if lab is not None:
     rec['clusters']=[]
     for k in sorted(set(lab)):
      ii=np.array(ids)[lab==k];a=fields[ii].astype(float)/255;mean=a.mean(0)
      rep=int(ii[np.argmin(np.mean(abs(a.reshape(len(ii),-1)-mean.ravel()),axis=1))])
      key=f'{exposure}__{side}__shape{k+1}'
      np.savez_compressed(d/(key+'.npz'),occupancy_probability=mean,site_indices=ii)
      m.write_dx(d/(key+'.dx'),mean)
      rec['clusters'].append(dict(cluster=k+1,n=len(ii),representative_site=rows[rep]['site_id']))
    clusters.append(rec)
  # Stratified deterministic grid/sampling sensitivity; keep primary results unchanged.
  pick=set()
  for group in sorted(set(r['exposure_group'] for r in rows)):
   ids=[i for i,r in enumerate(rows) if r['exposure_group']==group]
   if ids:pick.update(np.array(ids)[np.linspace(0,len(ids)-1,min(6,len(ids))).astype(int)].tolist())
  positive=[i for i,r in enumerate(rows) if r['external_shell_fraction']>0]
  pick.update(sorted(positive,key=lambda i:rows[i]['external_shell_fraction'])[:3])
  pick.update(sorted(range(len(rows)),key=lambda i:-rows[i]['unresolved_fraction'])[:3])
  for t in [.05,.15]:pick.update(sorted(range(len(rows)),key=lambda i:abs(rows[i]['external_shell_fraction']-t))[:3])
  pick=sorted(pick);atoms=np.load(R/'inputs_v2'/f'{mat}_parent_atoms.npz')
  env=m.Env(atoms['xyz_A'],atoms['elements'],atoms['box_A'])
  m.emit('refine_grid_start',material=mat,sites=len(pick))
  fine=m.build_grid(env,atoms['masses_Da'],.75)
  refined=[]
  for i in pick:
   r=rows[i];center=np.array(r['center'])
   _,loc,ext=m.exposure(center,env,fine,4096)
   _,loc2,ext2=m.exposure(center,env,fine,8192)
   q=center+grid@np.array(r['surface_frame']).T
   af=1-m.interp(fine['external'].astype(np.float32),q,fine,env)
   df=abs(af-fields[i].ravel()/255)
   refined.append(dict(site_id=r['site_id'],primary_exposure=r['external_shell_fraction'],grid075_exposure=float(ext.mean()),grid075_samples8192_exposure=float(ext2.mean()),grid_change=abs(float(ext.mean())-r['external_shell_fraction']),angular_change=abs(float(ext2.mean())-float(ext.mean())),field_mean_absolute_change_R15=float(df[sphere].mean()),primary_bin=r['exposure_group'],refined_raw_bin=m.exposure_bin(float(ext2.mean()))))
  maxgrid=max(x['grid_change'] for x in refined);maxang=max(x['angular_change'] for x in refined);maxfield=max(x['field_mean_absolute_change_R15'] for x in refined)
  rec=dict(material=mat,sites=len(rows),audit_sample_n=len(pick),max_exposure_fraction_grid_change=maxgrid,max_exposure_fraction_angular_change=maxang,max_field_MAE_R15=maxfield,grid_sensitivity_status='PASS_SELECTED_CASES_ONLY' if maxgrid<=.005 and maxfield<=.05 else 'SENSITIVE_SELECTED_CASES',sampling_status='PASS_SELECTED_CASES_ONLY' if maxang<=.005 else 'SENSITIVE_SELECTED_CASES',refinement=refined,clusters=clusters,input_sha256=m.sha(src/mat/'sites.json'),shape_proxy_status='DESCRIPTIVE_NOT_CLASSIFICATION_GROUND_TRUTH')
  summaries.append(rec);(d/'audit.json').write_text(json.dumps(rec,indent=2));m.emit('audit_material_complete',material=mat,maxgrid=maxgrid,maxfield=maxfield,clusters=[x['status'] for x in clusters])
 # Plot exact same exposure axes and outward-aligned y=0 cross-sections.
 plt.rcParams.update({'font.size':10,'svg.fonttype':'none'})
 fig,axs=plt.subplots(1,3,figsize=(12,3.5),constrained_layout=True)
 for ax,mat in zip(axs,['PET','PA6','PA66']):
  rr=allrows[mat];v=np.array([r['external_shell_fraction']*100 for r in rr])
  for side,col in [('top','#2378aa'),('bottom','#d87832')]:
   a=np.sort([r['external_shell_fraction']*100 for r in rr if r['side']==side])
   ax.step(a,np.arange(1,len(a)+1)/len(a),where='post',label=side,color=col)
  ax.set(title=f'{mat}, n={len(rr)}',xlabel='Carbon external SASA / full sphere (%)',ylabel='Cumulative fraction',xlim=(0,max(20,max(r['external_shell_fraction']*100 for rows in allrows.values() for r in rows))),ylim=(0,1))
  ax.axvline(5,color='0.7',lw=.6);ax.axvline(15,color='0.7',lw=.6);ax.legend(frameon=False)
 fig.savefig(out/'exposure_distributions.svg');plt.close(fig)
 fig,axs=plt.subplots(3,4,figsize=(13,10),sharex=True,sharey=True,constrained_layout=True)
 groups=['sampled_zero','low_0_5pct','mid_5_15pct','high_gt15pct']
 xx,zz=np.meshgrid(axis,axis,indexing='ij')
 for row,mat in enumerate(['PET','PA6','PA66']):
  for col,group in enumerate(groups):
   ax=axs[row,col];parts=[]
   for side,color in [('top','#2378aa'),('bottom','#d87832')]:
    path=src/mat/'fits'/f'{group}__{side}__support_all.npz'
    if not path.exists():continue
    a=np.load(path);p=a['occupancy_probability'][:,15,:].copy()
    p[xx*xx+zz*zz>225]=np.nan
    for lev,ls,lw in [(.25,':',.8),(.5,'-',1.5),(.75,'--',.8)]:
     if np.nanmin(p)<lev<np.nanmax(p):ax.contour(xx,zz,p,levels=[lev],colors=[color],linestyles=[ls],linewidths=[lw])
    parts.append(f'{side} n={len(a["site_indices"])}')
   ax.plot(0,0,'ko',ms=3);ax.set_aspect('equal');ax.set_xlim(-15,15);ax.set_ylim(-15,15)
   ax.set_title(f'{mat}: {group}\n'+', '.join(parts),fontsize=9)
   if not parts:ax.text(0,0,'No sites',ha='center',va='bottom')
   if row==2:ax.set_xlabel('Tangential x (A)')
   if col==0:ax.set_ylabel('Outward z (A)')
 fig.suptitle('Polymer-only interface cross-sections: solid P=0.50; dotted 0.25; dashed 0.75\nBlue top, orange bottom. Origin: carbonyl C. Water-center exclusion envelope, not molecular SES.',fontsize=11)
 fig.savefig(out/'interface_cross_sections.svg');plt.close(fig)
 doc=dict(status='AUDIT_COMPLETE_WITH_EXPLICIT_LIMITS',seconds=time.time()-start,materials=summaries)
 (out/'summary.json').write_text(json.dumps(doc,indent=2))
 with (R/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(event='polymer_exposure_audit_complete',output=str(out),exit_status=0,seconds=doc['seconds']))+'\n')
 m.emit('audit_complete',seconds=doc['seconds'])
if __name__=='__main__':main()
