import sys,json,ast,zlib,base64,hashlib
from pathlib import Path
import numpy as np
from scipy.ndimage import map_coordinates
R=Path(__file__).resolve().parent;sys.path.insert(0,str(R))
import whole_site_contact_v1 as w
m=w.m
OUT=R/'whole_site_panel_v1/publication';OUT.mkdir(exist_ok=False)
allrows=json.load(open(R/'panel_surface_scan_v1/sites.json'))
summary=[]
def decode(path):
 tree=ast.parse(path.read_text());node=next(n for n in tree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='D' for t in n.targets))
 return json.loads(zlib.decompress(base64.b64decode(ast.literal_eval(node.value.args[0].args[0].args[0]))))
for mat,cids in m.CASES.items():
 env=w.Environment(mat)
 for cid in cids:
  root=R/('whole_site_M329_pilot_v1' if mat=='PET' and cid=='M329' else 'whole_site_'+mat+'_'+cid+'_v1')
  rec=json.load(open(root/'summary.json'));old=np.load(R/'chemical_contact_surface_full_v3'/mat/cid/'matched_material_fields.npz')
  fields=old['fields'];rows=[allrows[i] for i in old['site_indices']];rot=old['chemical_to_core_frame']
  for t in rec['tiers']:
   pct=int(t['fraction']*100);d=root/('coverage'+str(pct));sub=d/'water_surface_views';sub.mkdir(exist_ok=False)
   ix=t['selected_indices'];rep=next(i for i,s in enumerate(rows) if s['site_id']==t['representative_site'])
   payload=[]
   for name,inds,field in [('material_water_consensus',ix,np.median(fields[ix],axis=0)),('material_water_real',[rep],fields[rep])]:
    v,f,n=m.mesh_patch(field);vals=[];weights=[];hits=[]
    for i in inds:
     ds=map_coordinates(fields[i],(v-m.AX[0]).T,order=1)
     wt=np.exp(-.5*(ds/1.5)**2);wt[abs(ds)>4.5]=0
     world=np.array(rows[i]['center'])+v@(np.array(rows[i]['frame'])@rot).T
     vals.append(env.ce.chemistry(world));weights.append(wt);hits.append(abs(ds)<=1.5)
    vals=np.array(vals);weights=np.array(weights);den=weights.sum(0)
    mean=(vals*weights[:,:,None]).sum(0)/np.maximum(den[:,None],1e-12)
    sd=np.sqrt(np.maximum((vals**2*weights[:,:,None]).sum(0)/np.maximum(den[:,None],1e-12)-mean**2,0))
    support=np.array(hits).mean(0);valid=(support>=.2)&(den>.01)
    stats,pay=m.export_mesh(sub,name,v,f,n,mean,sd,support,valid);payload.append(pay)
    for channel,k,target in [('HBA',1,np.array([.15,.55,.85])),('HBD',2,np.array([.7,.25,.7]))]:
     col=.93+(target-.93)*mean[:,k,None];col[~valid]=.55
     payload.append(dict(pay,name=name+'_'+channel,c=col.round(4).tolist()))
   orig=decode(d/(mat+'_'+cid+'_coverage'+str(pct)+'_view.py'))
   payload.append(next(x for x in orig['meshes'] if x['name']=='enzyme_contact'))
   fn=mat+'_'+cid+'_coverage'+str(pct)+'_material_view.py'
   m.viewer(sub/fn,payload,orig['core'],'material_water_consensus')
   ast.parse((sub/fn).read_text())
   t['water_surface_viewer']=str(sub/fn)
  summary.append(rec)
  print('MATERIAL_WATER_SURFACES',mat,cid,flush=True)
with (OUT/'all_core_results.json').open('x') as f:json.dump(summary,f,indent=2)
lines=['# 核心条件下的90/70/50%完整位点覆盖界面','',
'本轮使用同一固定核心、同一组完整位点，不使用逐点q分位数。每档配对给出材料朝水外表面、材料化学分布、真实代表位点和候选酶接触面。','',
'## 分母与验收','',
'分母是该固定核心原来的合格材料位点，不是全PET/PA6/PA66位点。目标位点数取ceil(N*比例)。先从每个位点作为起点贪心选择同一子集，在空间中尽量保留局部蛋白体积；各档验证前3个候选子集，再选择通过覆盖目标且接触面积较多的一个。这是启发式搜索，未证明全局最优。',
'原核心与供体坐标保持不变。沿用原重原子范德华半径和0.4 A重叠容差；对整张局部接触面顶点、边中点、面中心、局部1 A蛋白占据网格及固定核心逐位点检查。PASS是该离散几何口径的通过，不是完整蛋白可折叠或有活性。','',
'## 全部结果','',
'|材料|核心|分母|90%实际通过/要求|70%实际通过/要求|50%实际通过/要求|90/70/50%所选位点平均接触面积%|','|---|---|---:|---|---|---|---|']
for rec in summary:
 t=rec['tiers'];lines.append('|%s|%s|%d|%s|%s|%s|%s|'%(rec['material'],rec['core'],rec['denominator'],*['%d/%d (%s)'%(a['actual_coverage_n'],a['required_n'],'PASS' if a['subset_all_pass'] and a['actual_coverage_n']>=a['required_n'] else 'NOT_MET') for a in t],'/'.join('%.1f'%(a['selected_mean_contact_fraction']*100) for a in t)))
lines+=['','## 材料亲疏水及查看方式','',
'材料展示使用此前基于外部连通水得到的近似SES朝水表面，保持固定核心坐标系；碰撞检查则使用原筛选的重原子范德华体积。两者用途不同，不能拿SES图上的相交代替原子碰撞判定。每档的材料面、化学统计和酶面都来自同一份selected_site_ids。',
'片段亲脂性按完整化学图计算Wildman-Crippen原子贡献（显式H归入相邻重原子），8近邻1.5 A平滑。棕正/蓝负是经验片段指数，不能当成真实水合自由能；材料HBA/HBD有独立图层。灰色为共同近表面支持不足，不是疏水。','',
'PyMOL查看文件默认显示material_water_consensus。可禁用后启用material_water_real查看真实代表，*_HBA和*_HBD查看材料氢键特征；enzyme_contact是酶面。黄色D1/D2和球棍为固定供体/核心，青色点为羰基C，不是暴露面积。','',
'## 当前边界','',
'- 没有新增MD或供体采样，基于既有单快照和原位点名单。',
'- 共识面仍可能因平均而变平；真实代表面保留凹凸。没有人为要求平面，也没有为了显示凹凸强制制造曲率。',
'- 平面只是形状表现，主要目标是完整位点覆盖与接触；90%和50%子集不保证嵌套，三档名单均保存。',
'- 接触面积定义为界面顶点到材料重原子VDW表面的距离在[-0.4,1.5] A内的面积比例，随后按所选位点平均；这不是结合能。',
'- 蛋白骨架支撑、完整外形、折叠、动态稳健性、界面氢键方向与能量仍为NOT_EVALUATED。',
'- 与旧q系列不仅覆盖定义不同，碰撞表示也已统一回原筛选口径，旧新数值不可直接当作纯算法提升比较。',
'- 全部原始网格与54个成对查看文件保留在远程；GitHub发布摘要、图和M329三个代表查看文件，不是全量备份。','']
with (OUT/'README.md').open('x') as f:f.write('\n'.join(lines))
print('WHOLE_SITE_MATERIAL_PUBLICATION_READY',flush=True)
