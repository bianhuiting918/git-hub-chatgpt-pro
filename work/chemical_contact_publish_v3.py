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
O=P/'publication';O.mkdir(exist_ok=False)
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
audit={'mesh_count':0,'original_coverage_and_core_coordinates':'PASS','pymol_gui':'NOT_EVALUATED'}
for f in P.rglob('*.npz'):
 z=np.load(f)
 if 'vertices_A' not in z:continue
 v=z['vertices_A'];t=z['faces'];assert np.isfinite(v).all() and np.linalg.norm(v,axis=1).max()<=15.00001
 assert t.min()>=0 and t.max()<len(v)
 edges=np.concatenate([t[:,[0,1]],t[:,[1,2]],t[:,[2,0]]])
 graph=coo_matrix((np.ones(len(edges)),(edges[:,0],edges[:,1])),shape=(len(v),len(v)))
 assert connected_components(graph,directed=False)[0]==1
 _,counts=np.unique(np.sort(edges,axis=1),axis=0,return_counts=True);assert (counts==1).any()
 for key in ['material_HBA_weight','material_HBD_weight','surface_support_fraction','colors_rgb']:
  assert np.isfinite(z[key]).all() and z[key].min()>=0 and z[key].max()<=1.00001
 audit['mesh_count']+=1
panel=json.load(open(R/'panel_surface_scan_v1/panel.json'));sites=json.load(open(R/'panel_surface_scan_v1/sites.json'))
cov=np.load(R/'panel_surface_scan_v1/coverage.npz')['coverage'];cores=np.load(R/'panel_surface_scan_v1/cores.npz')
for mat in summary['materials']:
 for c in mat['core_conditioned']:
  cid=c['core_id'];pi=next(i for i,x in enumerate(panel) if x['id']==cid)
  ids=[i for i,s in enumerate(sites) if s['surface']==mat['material']+'_original' and cov[pi,i]]
  z=np.load(P/mat['material']/cid/'matched_material_fields.npz')
  assert ids==z['site_indices'].tolist() and len(ids)==c['coverage_n']
  rot=z['chemical_to_core_frame']
  assert np.allclose(rot.T@rot,np.eye(3),atol=1e-7)
  assert np.allclose(z['core_q_A']@rot.T,cores[cid+'_q'],atol=1e-7)
  assert np.allclose(z['donors_A']@rot.T,cores[cid+'_donors'],atol=1e-7)
# Check generated viewer payload and execute against a recording PyMOL API stub.
calls=[]
class Cmd:
 def __getattr__(self,name):
  def f(*args,**kwargs):
   if name=='load_cgo':assert len(args[0])>0
   calls.append(name)
  return f
pm=types.ModuleType('pymol');pm.cmd=Cmd();cg=types.ModuleType('pymol.cgo')
for i,k in enumerate('BEGIN END TRIANGLES NORMAL VERTEX COLOR SPHERE CYLINDER'.split()):setattr(cg,k,i+1)
sys.modules['pymol']=pm;sys.modules['pymol.cgo']=cg
for f in P.rglob('*_view.py'):
 d=decode(f)
 for a in d['meshes']:assert len(a['v'])==len(a['c'])==len(a['n']) and len(a['f'])>0
 exec(compile(f.read_text(),str(f),'exec'),{'print':lambda *x:None})
audit['viewer_payload_and_mock_execution_count']=len(list(P.rglob('*_view.py')))
audit['mesh_finite_connected_open_R15']='PASS'
# Selected portable viewers: all geometries included, separate chemistry/support channels.
selected={'PET':['SC1345656','M329'],'PA6':['SC4415918','N2240'],'PA66':['SC4874583','N0564']}
for mat,ids in selected.items():
 f=P/mat/'material_atlas'/(mat+'_gt2_le5pct_top.npz');z=data(f)
 meshes=[payload(z,mat+'_'+ch,ch) for ch in ['lipo','hba','hbd','support']]
 m.viewer(O/(mat+'_material_view.py'),meshes,default=mat+'_lipo')
 for cid in ids:
  source=P/mat/cid;d=decode(source/(mat+'_'+cid+'_view.py'))
  for q in [50,90,100]:
   z=data(source/('enzyme_contact_q'+str(q)+'.npz'))
   d['meshes'].append(payload(z,'support_q'+str(q),'support'))
  m.viewer(O/(mat+'_'+cid+'_view.py'),d['meshes'],d['core'],'enzyme_contact_q90')
fig=plt.figure(figsize=(11,10))
for row,mat in enumerate(selected):
 z=data(P/mat/'material_atlas'/(mat+'_gt2_le5pct_top.npz'))
 for col,ch in enumerate(['lipo','hba','hbd']):
  ax=fig.add_subplot(3,3,row*3+col+1,projection='3d')
  pay=payload(z,'tmp',ch);surface(ax,z,np.array(pay['c']));axes(ax,mat+' | '+{'lipo':'fragment lipophilicity','hba':'material H-bond acceptor','hbd':'material H-bond donor'}[ch])
fig.suptitle('Material-only consensus surfaces | carbonyl exposure 2-5%, top | radius 15 A',fontsize=13,y=.98)
fig.text(.5,.015,'Lipo: blue negative / brown positive (scale -0.5 to +0.5). HBA/HBD: white 0 to blue/purple 1.\nGray = insufficient near-surface support. Cyan dot = carbonyl C, not its exposed surface.',ha='center',fontsize=9)
fig.subplots_adjust(left=.02,right=.97,top=.92,bottom=.09,wspace=.02,hspace=.12)
fig.savefig(O/'material_chemistry.png',dpi=130);fig.savefig(O/'material_chemistry.svg',dpi=110);plt.close(fig)
fig=plt.figure(figsize=(10,11))
for row,(mat,ids) in enumerate(selected.items()):
 for col,cid in enumerate(ids):
  z=data(P/mat/cid/'enzyme_contact_q90.npz');d=decode(P/mat/cid/(mat+'_'+cid+'_view.py'))
  ax=fig.add_subplot(3,2,row*2+col+1,projection='3d');surface(ax,z,alpha=.65);drawcore(ax,d['core']);axes(ax,mat+' '+cid+' | '+('one-residue dual donor' if col==0 else 'two-residue donors'))
fig.suptitle('Fixed catalytic cores + candidate enzyme contact faces (q90)',fontsize=14,y=.98)
fig.text(.5,.015,'Cyan face: donor preference; purple: acceptor preference; green: mixed; gold: fixed-core surface; gray: low support.\nGreen/red/blue atoms = C/O/N of unchanged core; yellow points = donor N. Not a folded protein.',ha='center',fontsize=9)
fig.subplots_adjust(left=.03,right=.97,top=.93,bottom=.085,wspace=.02,hspace=.10)
fig.savefig(O/'core_contact_faces.png',dpi=130);fig.savefig(O/'core_contact_faces.svg',dpi=110);plt.close(fig)
with (O/'audit.json').open('x') as f:json.dump(audit,f,indent=2)
# Preserve complete detailed summary remotely; publish compact, traceable report.
lines=['# PET / PA6 / PA66 材料化学表面与固定核心接触面：v3','',
'状态：计算与网格审计完成；化学指标为经验代理，酶面为局部接触面假设；未验证完整蛋白、折叠、结合能或催化活性。','',
'## 本轮实际补齐','',
'- 材料本身：1436个位点，按原暴露分层、上下侧分开，排除68个数值/连通性不确定与2个法向不明确位点后，1366个位点构成29组材料共识面。没有用催化核心筛选材料总体。',
'- 固定核心：18个材料/核心组合，每个计算q50、q90、q100，共54个候选酶接触面；另有18个适配材料中位面。所有均为羰基碳15 A范围内、最近的连通开放网格，无方盒或球壳封口。',
'- 原始核心与供体坐标、适配位点名单均逐项核对，未增加供体采样或改变排名分母。','',
'## 查看','',
'![材料化学共识面](material_chemistry.svg)','',
'![固定核心接触面](core_contact_faces.svg)','',
'材料文件是2-5%暴露组的上表面代表，其他暴露组在远程完整输出中。下载下列 .py 到任意位置，在PyMOL用 File > Run 打开。无需另外下载原PDB或安装RDKit。当前仅完成语法、载荷和模拟API执行检查，实际PyMOL GUI为NOT_EVALUATED。一次打开一个文件，避免同名对象覆盖。','']
for mat,ids in selected.items():
 lines.append('- '+mat+': [材料表面]('+mat+'_material_view.py)；'+'；'.join('['+cid+']('+mat+'_'+cid+'_view.py)' for cid in ids))
lines+=['','材料对象可分别启用 *_lipo、*_hba、*_hbd、*_support。核心文件默认显示enzyme_contact_q90，可禁用后启用enzyme_contact_q50/q100；matched_material_median是材料面，support_q90是接触支持比例。','',
'## 全部18套结果','',
'下表“无穿入”仅指网格顶点按面积加权的材料穿入深度>0.4 A为零，不是所有原子无碰撞、可折叠或有活性。q90不是90%位点整面适配。','',
'|材料|核心|原核心适配数|q90平均穿入面积%|q90无穿入位点数|q100平均穿入面积%|q100无穿入位点数|q100化学支持面积%|',
'|---|---|---:|---:|---:|---:|---:|---:|']
for mat in summary['materials']:
 for c in mat['core_conditioned']:
  a,b=c['interfaces'][1:]
  lines.append('|%s|%s|%d/%d|%.3f|%d|%.3f|%d|%.1f|'%(mat['material'],c['core_id'],c['coverage_n'],c['denominator'],100*a['mean_penetrating_area_fraction_gt04A'],a['sites_with_zero_penetration_gt04A'],100*b['mean_penetrating_area_fraction_gt04A'],b['sites_with_zero_penetration_gt04A'],100*b['chemically_supported_area_fraction']))
lines+=['','## 科学解释与边界','',
'1. 亲疏水不是元素色：由完整SDF化学图和ITP逐键核对后计算Wildman-Crippen原子片段logP贡献，显式H归到相邻重原子；8近邻、1.5 A平滑。它不是水合自由能、表面真实logP或静电势。正负符号不能直接定义严格亲水/疏水面积。',
'2. 平均会稀释斑块：按羰基和法向对齐后，远端芳环/亚甲基/极性基团方向可能不一致。均值较均匀不证明每个真实位点均匀；NPZ保存标准差及接近材料面的位点支持比例。HBA/HBD应独立查看，不能被片段logP代替。',
'3. 酶化学面只表示接触偏好：面对材料HBA可考虑酶HBD，反之亦然；不意味着必须形成该氢键，不代表电荷反号，更没有评估去溶剂化。',
'4. 几何面由材料局部场分位数、0.5 A退让距离与固定核心vdW体积并集合成。它是局部接触边界，不是完整酶外形；只验证核心原子中心在局部体积内，未验证骨架连通、核心所有原子表面包埋、反应通道及动态进攻。',
'5. q50偏贴合但碰撞较多；q100更保守却可能失去共同接触支持。不能直接将覆盖最多的核心宣布为最佳完整酶。应在这些固定条件内，以整面碰撞与接触保留共同评价。',
'6. 单一模拟快照、继承原暴露定义与不确定标记。0.5 A全局网格、1 A局部网格；孤立碳测试的平移网格半径最大误差约0.24 A，不能声称亚网格精度。',
'7. 当前输出覆盖各暴露水平的材料独立共识面，以及各固定核心的全部匹配位点共识面；尚未把每套核心进一步按暴露层细分重拟合，也没有跨材料唯一通用酶面结论。','',
'## 版本更正','',
'- v1存在额外半网格偏移，保留为诊断产物。',
'- v2去掉偏移，但水侧距离截为-1.4 A，导致化学支持错误虚高至100%；v2化学支持统计作废。',
'- v3保留零等值面附近定义，增加水侧距离延拓并扩大XY周期边距至20 A。合成平面测试确认远端水距递减、材料侧场不变。',
'- 之前几何界面报告的“完成”不包含材料亲疏水图或酶接触面，本报告才补充这两类明确输出。历史结果不删除，不混作已验证酶结构。','',
'## 可复现性与存储','',
'远程完整输出：'+str(P),'',
'重跑：在项目根目录用既有CPU环境运行 chemical_contact_surface_v3.py --output NEW_VERSION_NAME；随后运行 chemical_contact_publish_v3.py（其输入固定为full_v3，输出必须不存在）。',
'脚本、运行日志和RUNBOOK保留在远程项目。GitHub本轮发布代表性查看脚本、审计、图和报告；全部NPZ/PLY未上传，不是全量备份。','',
'输入与计算source SHA256：'+summary['source_sha256'],'',
'## 方法来源','',
'- [RDKit Wildman-Crippen](https://www.rdkit.org/docs/source/rdkit.Chem.Crippen.html)',
'- [Euclidean distance transform分子表面方法](https://pmc.ncbi.nlm.nih.gov/articles/PMC2779860/)，仅作为方法依据，不替代本实现离散误差审计。',
'- [scikit-image marching cubes](https://scikit-image.org/docs/stable/api/skimage.measure.html#skimage.measure.marching_cubes)','']
with (O/'README.md').open('x') as f:f.write('\n'.join(lines))
manifest=[dict(file=str(f.relative_to(P)),bytes=f.stat().st_size,sha256=hashlib.sha256(f.read_bytes()).hexdigest()) for f in P.rglob('*') if f.is_file()]
with (O/'manifest.json').open('x') as f:json.dump(manifest,f,indent=2)
print(json.dumps(audit));print('PUBLICATION_READY',str(O))
