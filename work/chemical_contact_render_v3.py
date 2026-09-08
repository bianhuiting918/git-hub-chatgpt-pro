#!/usr/bin/env python
"""Audit and render chemical-contact v3 outputs without new donor sampling."""
import sys,json,ast,base64,zlib,hashlib,types
from pathlib import Path
sys.dont_write_bytecode=True
R=Path(__file__).resolve().parent
sys.path.insert(0,str(R))
import chemical_contact_surface_v3 as m
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
P=R/'chemical_contact_surface_full_v3'
O=P/'publication'
summary=json.load(open(P/'summary.json'))
def data(path):
 return {k:v for k,v in np.load(path).items()}
def payload(z,name,channel='lipo'):
 c=z['colors_rgb'].copy()
 if channel in ['hba','hbd']:
  v=z['material_HBA_weight' if channel=='hba' else 'material_HBD_weight']
  target=np.array([.15,.55,.85]) if channel=='hba' else np.array([.70,.25,.70])
  c=np.ones((len(v),3))*.93+(target-.93)*v[:,None]
 elif channel=='support':
  v=z['surface_support_fraction'];c=plt.get_cmap('viridis')(v)[:,:3]
 if channel!='support':c[~z['chemical_supported']]=.55
 return dict(name=name,v=z['vertices_A'].round(4).tolist(),f=z['faces'].tolist(),n=z['normals'].round(5).tolist(),c=c.round(4).tolist())
def decode(path):
 node=ast.parse(path.read_text())
 assign=next(n for n in node.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='D' for t in n.targets))
 blob=ast.literal_eval(assign.value.args[0].args[0].args[0])
 return json.loads(zlib.decompress(base64.b64decode(blob)))
def axes(ax,title):
 ax.set_title(title,fontsize=11);ax.set_xlim(-15,15);ax.set_ylim(-15,15);ax.set_zlim(-15,15)
 ax.set_box_aspect((1,1,1));ax.view_init(elev=32,azim=-65)
 ax.set_xticks([-10,0,10]);ax.set_yticks([-10,0,10]);ax.set_zticks([-10,0,10])
 ax.tick_params(labelsize=7,pad=0);ax.set_xlabel('x (A)',fontsize=8,labelpad=-1);ax.set_ylabel('y (A)',fontsize=8,labelpad=-1);ax.set_zlabel('z (A)',fontsize=8,labelpad=-1)
 ax.scatter([0],[0],[0],s=28,c='cyan',edgecolors='black',linewidths=.4,depthshade=False)
def surface(ax,z,colors=None,alpha=1):
 v=z['vertices_A'];f=z['faces'];c=z['colors_rgb'] if colors is None else colors
 obj=Poly3DCollection(v[f],facecolors=c[f].mean(1),edgecolors='none',alpha=alpha,rasterized=True)
 ax.add_collection3d(obj)
def drawcore(ax,c):
 q=np.array(c['q']);e=c['elements'];colors={'C':'#238b45','O':'#e3342f','N':'#2446ce','S':'#e3c62d'}
 for i,j in c['bonds']:ax.plot(*q[[i,j]].T,color='#343434',lw=1.7)
 ax.scatter(*q.T,c=[colors.get(x,'gray') for x in e],s=12,depthshade=False)
 d=np.array(c['donors']);ax.scatter(*d.T,c='yellow',edgecolors='black',s=35,depthshade=False)
 for i,a in enumerate(d):ax.text(*a,'D'+str(i+1),fontsize=8)

selected={'PET':['SC1345656','M329'],'PA6':['SC4415918','N2240'],'PA66':['SC4874583','N0564']}
fig=plt.figure(figsize=(11,12))
for row,mat in enumerate(selected):
 z=data(P/mat/'material_atlas'/(mat+'_gt2_le5pct_top.npz'))
 for col,ch in enumerate(['lipo','hba','hbd']):
  ax=fig.add_subplot(3,3,row*3+col+1,projection='3d')
  pay=payload(z,'tmp',ch);surface(ax,z,np.array(pay['c']));axes(ax,mat+' | '+{'lipo':'fragment lipophilicity','hba':'material H-bond acceptor','hbd':'material H-bond donor'}[ch])
fig.suptitle('Material-only consensus surfaces | carbonyl exposure 2-5%, top | radius 15 A',fontsize=13,y=.98)
fig.text(.5,.015,'Lipo: blue negative / brown positive (scale -0.5 to +0.5). HBA/HBD: white 0 to blue/purple 1.\nGray = insufficient near-surface support. Cyan dot = carbonyl C, not its exposed surface.',ha='center',fontsize=9)
fig.subplots_adjust(left=.02,right=.97,top=.92,bottom=.09,wspace=.02,hspace=.30)
fig.savefig(O/'material_chemistry.png',dpi=130);fig.savefig(O/'material_chemistry.svg',dpi=110);plt.close(fig)
fig=plt.figure(figsize=(11,13))
for row,(mat,ids) in enumerate(selected.items()):
 for col,cid in enumerate(ids):
  z=data(P/mat/cid/'enzyme_contact_q90.npz');d=decode(P/mat/cid/(mat+'_'+cid+'_view.py'))
  ax=fig.add_subplot(3,2,row*2+col+1,projection='3d');surface(ax,z,alpha=.65);drawcore(ax,d['core']);axes(ax,mat+' '+cid+' | '+('one-residue dual donor' if col==0 else 'two-residue donors'))
fig.suptitle('Fixed catalytic cores + candidate enzyme contact faces (q90)',fontsize=14,y=.98)
fig.text(.5,.015,'Cyan face: donor preference; purple: acceptor preference; green: mixed.\nGold: fixed-core surface; gray: low support. Green/red/blue atoms = C/O/N of unchanged core; yellow points = donor N. Not a folded protein.',ha='center',fontsize=9)
fig.subplots_adjust(left=.03,right=.97,top=.93,bottom=.085,wspace=.02,hspace=.30)
fig.savefig(O/'core_contact_faces.png',dpi=130);fig.savefig(O/'core_contact_faces.svg',dpi=110);plt.close(fig)
