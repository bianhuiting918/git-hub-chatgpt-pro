import sys,json,ast,base64,zlib
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
R=Path(__file__).resolve().parent;sys.path.insert(0,str(R))
import whole_site_contact_v1 as w
P=R/'whole_site_M329_pilot_v1';O=P/'water_preview';O.mkdir(exist_ok=False)
fig=plt.figure(figsize=(16,12))
titles=['Actual water-facing PET: lipo','Selected PET consensus: lipo','Selected PET consensus: HBA','Candidate enzyme contact face']
summary=json.load(open(P/'summary.json'));q=np.load(P/'material_fields_heavy_vdw.npz')['core_q_A']
for row,pct in enumerate([90,70,50]):
 d=P/('coverage'+str(pct));rec=json.load(open(d/'summary.json'))
 for col in range(4):
  key='material_water_real' if col==0 else 'enzyme_contact' if col==3 else 'material_water_consensus'
  z=np.load((d if col==3 else d/'water_surface_views')/(key+'.npz'));v=z['vertices_A'];f=z['faces'];colors=z['colors_rgb'].copy()
  if col==2:
   values=z['material_HBA_weight'];colors=.93+(np.array([.15,.55,.85])-.93)*values[:,None];colors[~z['chemical_supported']]=.55
  ax=fig.add_subplot(3,4,row*4+col+1,projection='3d')
  ax.add_collection3d(Poly3DCollection(v[f],facecolors=colors[f].mean(1),edgecolors='none',alpha=.85,rasterized=True))
  if col==3:ax.scatter(*q.T,c='darkgreen',s=6,depthshade=False)
  ax.scatter([0],[0],[0],c='cyan',s=25,edgecolors='black',depthshade=False)
  ax.set_xlim(-15,15);ax.set_ylim(-15,15);ax.set_zlim(-15,15);ax.set_box_aspect((1,1,1));ax.view_init(25,-65)
  ax.set_xticks([-10,0,10]);ax.set_yticks([-10,0,10]);ax.set_zticks([-10,0,10]);ax.tick_params(labelsize=7,pad=0)
  ax.set_xlabel('x (A)',fontsize=8,labelpad=-1);ax.set_ylabel('y (A)',fontsize=8,labelpad=-1);ax.set_zlabel('z (A)',fontsize=8,labelpad=-1)
  ax.set_title(str(pct)+'% ('+str(rec['actual_coverage_n'])+'/29)\n'+titles[col],fontsize=10)
fig.suptitle('PET M329 | fixed whole-site subsets, paired material chemistry and enzyme faces',fontsize=15,y=.985)
fig.text(.5,.015,'Each row uses one fixed PET site subset. Lipo: brown positive / blue negative empirical fragment index. HBA: white 0 to blue 1.\nEnzyme: cyan donor preference; gold fixed-core area; gray low chemical support. Geometry-only coverage, not binding energy or foldability.',ha='center',fontsize=10)
fig.subplots_adjust(left=.01,right=.98,top=.93,bottom=.09,wspace=.10,hspace=.30)
fig.savefig(O/'M329_paired_water_material_enzyme.png',dpi=120);fig.savefig(O/'M329_paired_water_material_enzyme.svg',dpi=100)
print('PREVIEW_READY')
