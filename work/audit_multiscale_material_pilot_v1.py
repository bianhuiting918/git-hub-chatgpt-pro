import pathlib,json,sys,ast,hashlib,numpy as np
from scipy.ndimage import map_coordinates
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
R=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(R))
import multiscale_material_match_v1 as a
OUT=R/"multiscale_material_match_pilot_report_v1";OUT.mkdir(exist_ok=True)
roots={mat:R/("multiscale_material_match_pilot_v1" if mat=="PET" else "multiscale_material_match_"+mat+"_pilot_v1") for mat in ["PET","PA6","PA66"]}
doc=[];checks=[];manifest=[]
for mat,root in roots.items():
 d=json.load(open(root/"summary.json"));assert d["pilot"] and d["evaluated_sites"]==8 and len(d["results"])==18
 rows=json.load(open(root/"sites.json"))
 for frac in [.9,.7,.5]:
  for group in ["all","current_core_accessible"]:
   rr=[x for x in d["results"] if x["cohort"]==group and x["target_fraction"]==frac]
   assert len({tuple(x["chosen_site_ids"]) for x in rr})==1
 for x in d["results"]:
  key=f"{mat}_{x['cohort']}_coverage{int(x['target_fraction']*100)}_R{x['radius_A']}"
  folder=root/key;z=np.load(folder/(key+".npz"));v=z["vertices_A"];f=z["faces"]
  assert np.isfinite(v).all() and np.isfinite(z["colors_rgb"]).all() and len(f)>0
  assert np.linalg.norm(v,axis=1).max()<=x["radius_A"]+1e-5
  ii=np.r_[f[:,0],f[:,1],f[:,2]];jj=np.r_[f[:,1],f[:,2],f[:,0]]
  n,_=connected_components(coo_matrix((np.ones(len(ii)),(ii,jj)),shape=(len(v),len(v))),directed=False)
  assert n==1
  err={e["site_id"]:e["D95_A"] for e in x["per_site_errors"]}
  assert all(err[s]<=x["tolerance_for_same_subset_A"]+1e-8 for s in x["chosen_site_ids"])
  assert x["actual_n_at_tolerance"]>=x["required_n"]
  if mat=="PET":assert z["material_HBD_weight"].max()==0
  p=folder/(key+"_view.py");ast.parse(p.read_text())
  manifest.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
  checks.append(dict(key=key,finite=True,one_component=True,within_radius=True,chosen_subset_matches=True))
 doc.append(d)
# Existing-field PET controls: separate support filtering from frame choice.
root=roots["PET"];z=np.load(root/"fields_R22.npz");a.m.AX=z["axis_A"];rows=json.load(open(root/"sites.json"))
panel={s["site_id"]:s for s in json.load(open(R/"panel_surface_scan_v1/sites.json")) if s["surface"]=="PET_original"}
changed=[]
for s,field in zip(rows,z["fields"]):
 ps=panel[s["site_id"]];co=np.array(ps["frame"])@np.array(ps["amide"])[1];leave=np.array(ps["frame"])@np.array(ps["amide"])[2]
 new=a.m.old.make_frame([0,0,1 if s["side"]=="top" else -1],co,leave)
 q=a.GRID@new.T@np.array(s["surface_frame"])
 changed.append(map_coordinates(field.astype(float),(q-a.AX[0]).T,order=1,mode="constant",cval=-30).reshape(field.shape))
controls=[]
for frame,fields in [("local_density_normal",z["fields"]),("slab_normal_control",changed)]:
 for radius in [10,15,20]:
  rep=[a.mesh_rep(*a.m.mesh_patch(f.astype(float),radius=radius)[:2]) for f in fields]
  errors=[a.distance(rep[i],rep[j])[0] for i in range(len(rep)) for j in range(i)]
  controls.append(dict(frame=frame,support_filter=False,radius_A=radius,pairwise_D95_quantiles_A=np.quantile(errors,[0,.25,.5,.75,1]).tolist()))
audit=dict(status="TECHNICAL_PILOT_CHECKS_PASS_NOT_SCIENTIFIC_ACCEPTANCE",materials=3,sites_evaluated=24,meshes_and_viewers_checked=len(checks),checks=checks,control_results=controls,manifest=manifest,unverified=["PyMOL GUI","full population fitting","global optimality","grid convergence","independent validation","enzyme-interface fit at R20"])
for name,data in [("audit.json",audit),("all_pilot_results.json",doc)]:
 with (OUT/name).open("x") as f:json.dump(data,f,indent=2,allow_nan=False)
lines=["# 三尺度材料表面：实际形状匹配试算","",
"**状态：试算完成，非全量结果；未证明单一共性面能在小误差下覆盖50/70/90%的总体位点。**",
"PET、PA6、PA66各确定性抽取8个位点，其中4个属于当前核心面板可接近集合。各材料2集合×3尺度×3档，共54项。全体集合的8个位点已包含那4个可接近位点，不是每种材料12个。",
"",
"![匹配所需容差](pilot_required_tolerance.svg)",
"",
"## 70%档实际需要的形状容差",
"",
"|材料|试算集合|半径10 A|半径15 A|半径20 A|",
"|---|---|---:|---:|---:|"]
for d in doc:
 for group in ["all","current_core_accessible"]:
  vals=[next(x for x in d["results"] if x["cohort"]==group and x["target_fraction"]==.7 and x["radius_A"]==r)["tolerance_for_same_subset_A"] for r in [10,15,20]]
  lines.append("|"+d["material"]+"|"+("全部抽样位点 N=8" if group=="all" else "可接近抽样位点 N=4")+"|"+"|".join(f"{v:.2f} A" for v in vals)+"|")
lines+=["","这里70%要求ceil(8×0.7)=6或ceil(4×0.7)=3；相同R20选中位点用于R10/R15，较小尺度上可能额外匹配其他位点。这个容差不是先验合格阈值，更不能因为把容差放大后计数达到目标就称为良好匹配。",
"","## 图与化学图层",""]
for mat in roots:
 for group in ["all","current_core_accessible"]:
  fn=mat+"_"+group+"_pilot_surfaces.svg"
  lines.append(f"- [{mat} {group} 三尺度亲脂性/HBA/HBD]({fn})")
lines+=["","所有图均为70%档候选的试算展示。亲脂性、HBA、HBD独立显示；颜色为后注释，没有参与形状优化，也不是水合自由能或表面电势。",
"","## 方法与局限","",
"- 原始完整材料快照、XY周期、外部水连通近似SES，保持0.5 A全局网格；局部网格1 A，扩展至±22 A。未改变MD或供体采样。",
"- 原版只有中心位点的块体筛选；本试算增加网格点内侧密度支撑>=70%的筛选（R10圆盘81点、向内5 A、沿固定局部法向）。这是剔除松散表面部分的几何近似，不是链拓扑判定；过滤前后面积保存在远程mesh_scope_audit.json。",
"- 同一材料局部法向+C=O投影坐标，三尺度使用同一候选场和选中位点，不做逐位点自由旋转来美化匹配。",
"- D95为两个方向各自按表面积加权的95%最近采样面距离之较大者；另存最大采样偏差。不是全部表面最大距离，也不是原子碰撞。",
"- 候选仅有样本真实面和中位数面，有限库内按R20所需容差选择；不是全局最优，也不足以否定其他拟合/分群方法。没有人为削平，也没有强制制造凹凸。",
"- 1/2/3 A下的实际计数均保存在all_pilot_results.json；缺失/无法评估必须保留，不能暗中改变分母。",
"- 全量位点未计算：PET578/当前面板可接近131；PA6 413/47；PA66 445/58。PA6全量2个位点法向不确定需保留NOT_EVALUATED。本试算不代表这些总体分布。",
"- 历史473套PET双残基、3021套尼龙双残基、5151932套单残基耦合候选与当前101套面板不同；本轮没有完成历史全库的供体条件界面拟合。",
"",
"## 诊断与下一步",
"",
"PET不做网格点支撑过滤时，R10/R15/R20真实面两两D95中位数约4.64/7.86/10.22 A；改用统一slab法向的对照约4.50/6.51/7.97 A。框架选择影响部分差异，但大偏差不只由新增过滤导致。该对照没有替换主方法。",
"建议先在全量位点上建立形状类别，分别报告每类共识、类内误差与总体占比，再研究多个界面共同覆盖；不得将多个类别覆盖之和标为同一酶面覆盖率。单面全量拟合未启动，不能宣称科学目标已完成。",
"",
"## 文件与复现",
"",
"服务器项目根目录：/data/bht2/polymer_material_reference_20260827/simulation_slabs/interface_surface_robustness_20260908_v1",
"脚本：multiscale_material_match_v1.py、test_multiscale_material_match_v1.py、audit_multiscale_material_pilot_v1.py、render_multiscale_material_pilot_v1.py。运行日志在RUN_LOG.jsonl，说明追加于RUNBOOK.md。",
"GitHub发布本报告、汇总JSON与图，不宣称全量NPZ/PLY/PyMOL文件均已备份。本地未下载，新PyMOL GUI未验证。",
""]
with (OUT/"README.md").open("x") as f:f.write("\n".join(lines))
print(json.dumps(dict(checked=len(checks),report=str(OUT),status=audit["status"])))
