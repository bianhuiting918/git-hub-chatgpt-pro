import pathlib,json,sys,datetime,hashlib,numpy as np
R=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(R))
import tangent_probe_v1 as p
OUT=R/"tangent_probe_publication_20260909_v1";OUT.mkdir(exist_ok=False)
roots=[R/"tangent_probe_pilot_20260909_v1",R/"tangent_probe_surface_controls_20260909_v1"]
docs=[json.load(open(x/"summary.json")) for x in roots]
oldrows=json.load(open(R/"panel_surface_scan_v1/sites.json"))
checks=[];rows=[];corrections=[];manifest=[]
for root,doc in zip(roots,docs):
 for mat in doc["materials"]:
  oldids={s["site_id"].rsplit("_",1)[0]:s["site_id"] for s in oldrows if s["surface"]==mat["material"]+"_original"}
  for rec in mat["results"]:
   z=np.load(root/mat["material"]/(rec["site_id"]+"_directions.npz"))
   a=z["accessible"];e=z["endpoint"];u=z["unresolved"];n=len(z["directions_chemical"])
   assert a.shape==e.shape==u.shape==(10,1024)
   assert not np.any(a&~e) and not np.any(a&u) and not np.any(u&~e)
   assert not np.any(a[1:]&~a[:-1]) and not np.any(e[1:]&~e[:-1])
   assert np.allclose(np.linalg.norm(z["directions_chemical"],axis=1),1)
   aa=a.mean(1)
   assert np.allclose(aa,[x["A"] for x in rec["radii"]])
   assert abs(p.score(p.RADII,aa)-rec["S_log"])<1e-12
   assert 0<=rec["S_log"]<=1
   for k,x in enumerate(rec["radii"]):
    assert x["endpoint_n"]==int(e[k].sum()) and x["accessible_n"]==int(a[k].sum()) and x["unresolved_n"]==int(u[k].sum())
    assert x["largest_fraction"]<=x["A"]+1e-10
   legacy=oldids.get(rec["site_id"].rsplit("_",1)[0])
   if "in_legacy_ge70_cohort" in rec and rec["in_legacy_ge70_cohort"]!=(legacy is not None):
    corrections.append(dict(material=mat["material"],site_id=rec["site_id"],legacy_site_id=legacy,field="in_legacy_ge70_cohort",reason="old PET labels use leaving O; new labels use carbonyl O; match carbon identity"))
   rec=dict(rec,legacy_site_id=legacy,in_legacy_ge70_cohort=legacy is not None,cohort="parent_pilot" if root==roots[0] else "surface_positive_control")
   rows.append(rec);checks.append(dict(material=mat["material"],site=rec["site_id"],cohort=rec["cohort"],checks="PASS"))
   path=root/mat["material"]/(rec["site_id"]+"_directions.npz");manifest.append(dict(path=str(path),sha256=p.old.sha(path)))
# Exact collision segmentation should not alter answer.
seg=[]
for mat in docs[1]["materials"]:
 ce=p.chem.ChemEnv(mat["material"]);cat={x["site_id"]:x for x in p.catalogue(mat["material"],ce)}
 rec=max(mat["results"],key=lambda x:x["S_log"]);s=cat[rec["site_id"]];u=p.old.fibonacci(1024)@np.array(s["frame"]).T
 for r in [.5,2.,6.]:
  a=p.trace(ce.env,np.array(s["center"]),u,r,step=2)
  b=p.trace(ce.env,np.array(s["center"]),u,r,step=1)
  assert all(np.array_equal(a[k],b[k]) for k in ["endpoint","reached","blocked","unresolved"])
  seg.append(dict(material=mat["material"],site=s["site_id"],radius_A=r,step1_vs_step2="EXACT_MASK_MATCH"))
ref=[]
for label,doc in zip(["parent_pilot","surface_control"],docs):
 for mat in doc["materials"]:
  rs=mat["refinement"] if isinstance(mat["refinement"],list) else [mat["refinement"]]
  ref.extend([dict(material=mat["material"],cohort=label,**r) for r in rs])
audit=dict(status="PILOT_TECHNICAL_AUDIT_PASS",site_records=len(rows),unique_material_carbon_sites=len({(x["material"],x["site_id"].rsplit("_",1)[0]) for x in rows}),radius_records=len(rows)*10,checks=checks,segmentation_checks=seg,refinement=ref,largest_refinement_absolute_A_change=max(x["max_absolute_A_change"] for x in ref),metadata_corrections=corrections,manifest=manifest,initial_control_process_exit=1,logging_error="NameError Path; numerical data persisted; logging fix separately recorded; no numerical recompute",not_evaluated=["all-population1024-direction classification","loose-tail/model-end exclusion","curved approach paths","triad-donor scan"])
with (OUT/"audit.json").open("x") as f:json.dump(audit,f,indent=2)
with (OUT/"results.json").open("x") as f:json.dump(dict(definitions=docs[0]["definitions"],radii_A=p.RADII.tolist(),results=rows),f,indent=2)
lines=["# 羰基碳相切球探针：可及性试算","",
"**试算与技术核验已完成；不是全量表面分级，也没有新增三联体—供体扫描。**",
"",
"## 已完成范围","",
"- 完整模型枚举：PET8400个C(=O)-O-C键、PA6和PA66各6800个C(=O)-N-C键。没有支撑/密度阈值筛选；这是模型键数，包含模型链端，不能当成块体表面键数。",
"- 全模型64方向、r=0.5 A的终点初筛仅用于选取试算样本，不是最终暴露分级。",
"- 每材料12个完整模型分层样本，另加12个原外表面阳性校验点，共72条位点计算记录。表面校验点来自旧支撑筛选集合，只用于验证，不能作为全体表面的无偏抽样。",
"- 每个位点10种球半径、1024个方向；部分加密2048方向；9项路径分段1 A与2 A对照完全一致。",
"",
"## 几何定义","",
"球心=C+(1.7+r)u；目标C及所有邻近原子（包括羰基O、酯O/酰胺N和显式H）均保留。严格相切允许，穿透不允许，数值容差1e-7 A。",
"先检查相切终点，再检查球沿同一直线到外部的扫掠体。每段解析计算球与原子的碰撞，不是只查离散路径点；XY周期性保持。外部判据为超过整个原始模型所有原子VDW范围的上/下自由半空间；路径最长120 A，未到达但未碰撞者记unresolved、不计成功。",
"半径：0.5、1、1.5、2、3、4、5、6、8、10 A。方向分母为完整4pi球面，化学局部坐标x沿C=O、y位于离去原子平面。没有施加亲核进攻角。",
"",
"A(r)=成功到达外部的方向数/1024。S_log为A(r)对log(r)积分后归一化，范围0..1；它没有密度或水暴露项，也不是活性或表面积百分数。保留A(r)全曲线比单一S更重要。",
"最大半径只是此次离散采样中的最大可达值。10 A通过只能写>=10 A；零方向为采样未检出，不代表严格无通路。最大连续开放区域由球面邻接图估计，不能视为已证明的连续开放圆锥。",
"",
"## 表面校验点结果","",
"|材料|旧位点ID|S_log|r1 A可及方向%|r2 A可及方向%|r4 A可及方向%|采样最大半径|",
"|---|---|---:|---:|---:|---:|---|"]
for mat in ["PET","PA6","PA66"]:
 for x in sorted([q for q in rows if q["material"]==mat and q["cohort"]=="surface_positive_control"],key=lambda q:-q["S_log"]):
  vals=[next(a["A"] for a in x["radii"] if a["radius_A"]==r)*100 for r in [1,2,4]]
  mx=(">=10 A" if x["Rmax_right_censored"] else (str(x["sampled_Rmax_A"])+" A" if x["sampled_Rmax_A"] is not None else "未检出"))
  lines.append("|"+mat+"|"+str(x["legacy_site_id"])+"|"+f'{x["S_log"]:.5f}'+"|"+"|".join(f"{v:.3f}" for v in vals)+"|"+mx+"|")
lines+=["","## 核验与解释边界","",
f'- {len(rows)}条位点记录、{len(rows)*10}项半径结果完成掩码/计数/得分/方向归一化核验；同方向的大球可及集合没有超出小球集合。',
f'- 1024→2048方向加密的最大可及比例变化为{audit["largest_refinement_absolute_A_change"]:.6f}（绝对比例）；窄入口仍有采样不确定性。',
"- PET新旧ID使用的氧名称不同，已按同一链的羰基碳及坐标核对别名；原始数值不变。结果JSON中的legacy_site_id用于关联旧扫描。",
"- 阳性校验仍使用旧名单，只证明本方法在已知表面点上能识别开放/受阻，不证明全量分级已经完成。",
"- 球形探针是形状与方向开放度的代理；非球形催化核心可能表现不同。后续必须扫描保持固定相对关系的真实三联体—双供体。",
"- 此轮未排除所有松散链尾/模型端基，未检验弯曲路径，也没有做MD、静电或能量计算。",
"",
"## 文件与复现","",
"[逐位点全部曲线、得分与方向区统计](results.json) · [核验与原始方向文件SHA256](audit.json)",
"",
"服务器根目录：/data/bht2/polymer_material_reference_20260827/simulation_slabs/interface_surface_robustness_20260908_v1",
"脚本：tangent_probe_v1.py --output NEW_UNIQUE_NAME --pilot 12 --directions 1024；表面校验脚本tangent_probe_surface_controls_v1.py，需使用新的OUT目录复跑以保护旧结果。运行与日志恢复记录已保留在RUN_LOG.jsonl。",
"原始方向掩码保存在服务器，GitHub仅发布脚本、结果和核验摘要，不声称全部NPZ已备份。本地未下载。",
""]
with (OUT/"README.md").open("x") as f:f.write("\n".join(lines))
with (R/"RUN_LOG.jsonl").open("a") as f:f.write(json.dumps(dict(time=datetime.datetime.now(datetime.timezone.utc).isoformat(),script=str(pathlib.Path(__file__)),output=str(OUT),exit_status=0,checks=len(checks),segmentation_checks=len(seg)))+"\n")
print(json.dumps({k:v for k,v in audit.items() if k in ["status","site_records","unique_material_carbon_sites","radius_records","largest_refinement_absolute_A_change","metadata_corrections"]}))
