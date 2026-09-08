#!/usr/bin/env python
"""Supplementary fine bins. Reuses primary arrays; no new exposure evaluation."""
from pathlib import Path
import json,numpy as np
import polymer_exposure_v2 as m
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
R=Path(__file__).resolve().parent
GROUPS=['sampled_zero','positive_le1pct','gt1_le2pct','gt2_le5pct','gt5pct']
def binof(v):
 if v==0:return GROUPS[0]
 if v<=.01:return GROUPS[1]
 if v<=.02:return GROUPS[2]
 if v<=.05:return GROUPS[3]
 return GROUPS[4]
def main():
 src=R/'polymer_exposure_full_v2';out=R/'polymer_exposure_fine_strata_v1';out.mkdir();result=[]
 axis=np.arange(-15,16);coords=np.stack(np.meshgrid(axis,axis,axis,indexing='ij'),-1);sphere=np.linalg.norm(coords,axis=-1)<=15
 plt.rcParams.update({'font.size':9,'svg.fonttype':'none'})
 fig,axs=plt.subplots(3,5,figsize=(15,9),sharex=True,sharey=True,constrained_layout=True)
 xx,zz=np.meshgrid(axis,axis,indexing='ij')
 for ri,mat in enumerate(['PET','PA6','PA66']):
  rows=json.load(open(src/mat/'sites.json'));fields=np.load(src/mat/'site_fields.npz')['fields'];d=out/mat;d.mkdir()
  bins=[binof(r['external_shell_fraction']) if binof(r['external_shell_fraction'])==binof(r['local_shell_fraction']) else 'CONNECTIVITY_UNRESOLVED' for r in rows]
  rec=dict(material=mat,n=len(rows),counts={g:bins.count(g) for g in GROUPS+['CONNECTIVITY_UNRESOLVED']},groups=[])
  for gi,group in enumerate(GROUPS):
   ax=axs[ri,gi];parts=[]
   for side,color in [('top','#2378aa'),('bottom','#d87832')]:
    ids=[i for i,r in enumerate(rows) if bins[i]==group and r['side']==side and r['normal_status']=='PASS']
    key=group+'__'+side
    if not ids:continue
    arr=fields[ids].astype(np.float32)/255;mean=arr.mean(0);sd=arr.std(0)
    dev=np.mean(abs(arr[:,sphere]-mean[sphere]),axis=1);rep=ids[int(np.argmin(dev))]
    np.savez_compressed(d/(key+'.npz'),occupancy_probability=mean,occupancy_sd=sd,site_indices=ids,sphere_mask=sphere)
    m.write_dx(d/(key+'.dx'),mean)
    # Do not hide group size or merge top/bottom.
    p=mean[:,15,:].copy();p[xx*xx+zz*zz>225]=np.nan
    for lev,ls,lw in [(.25,':',.7),(.5,'-',1.5),(.75,'--',.7)]:
     if np.nanmin(p)<lev<np.nanmax(p):ax.contour(xx,zz,p,levels=[lev],colors=[color],linestyles=[ls],linewidths=[lw])
    rec['groups'].append(dict(group=group,side=side,n=len(ids),chains=len(set(rows[i]['chain'] for i in ids)),representative_site=rows[rep]['site_id'],mean_absolute_field_deviation=float(dev.mean()),status='DESCRIPTIVE_ONLY'))
    parts.append(f'{side} {len(ids)}')
   ax.plot(0,0,'ko',ms=3);ax.set_aspect('equal');ax.set(xlim=(-15,15),ylim=(-15,15))
   ax.set_title(mat+' '+['sampled 0','0-1%','1-2%','2-5%','>5%'][gi]+'\n'+', '.join(parts),fontsize=9)
   if not parts:ax.text(0,1,'No resolved sites',ha='center',fontsize=8)
   if ri==2:ax.set_xlabel('Tangential x (A)')
   if gi==0:ax.set_ylabel('Outward z (A)')
  (d/'site_bin_assignments.json').write_text(json.dumps([dict(site_id=r['site_id'],fine_bin=b) for r,b in zip(rows,bins)],indent=2))
  result.append(rec)
 fig.suptitle('Polymer-only interface by carbon exposure: y=0 sections, radius 15 A\nBlue: top; orange: bottom. Solid P=0.50; dotted 0.25; dashed 0.75. Carbonyl C at origin.\nExternal water-center exclusion envelope, not atomic SES. Percent denominator: full inflated C sphere.',fontsize=11)
 fig.savefig(out/'fine_interface_sections.svg');fig.savefig(out/'fine_interface_sections.png',dpi=130);plt.close(fig)
 summary=dict(status='TECHNICAL_COMPLETE_DESCRIPTIVE_FINE_STRATA',primary=str(src),source_sha256=m.sha(__file__),bin_scope='supplementary reporting bins; same values and denominator; uncertainty evaluated against each finer boundary',materials=result)
 (out/'summary.json').write_text(json.dumps(summary,indent=2))
 with (R/'RUN_LOG.jsonl').open('a') as f:f.write(json.dumps(dict(event='polymer_exposure_fine_strata_complete',output=str(out),exit_status=0))+'\n')
 print(json.dumps(summary),flush=True)
if __name__=='__main__':main()
