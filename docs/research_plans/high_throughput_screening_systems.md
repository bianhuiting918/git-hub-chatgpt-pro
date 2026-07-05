# 塑料降解酶与尼龙酶高通量筛选体系

_Last updated: 2026-07-05_

本文件补充两个课题都必须建立的实验闭环：

```text
计算筛选 / 设计
        ↓
一级高通量活性初筛
        ↓
二级产物谱确认
        ↓
三级真实材料验证
        ↓
数据回流到模型 / 机制假设
```

核心原则：**一级筛选追求通量，二级验证追求化学特异性，三级验证追求真实材料相关性。** 不建议只依赖单一信号，例如 halo、pNP、总胺、总荧光或单一 LC-MS 产物。

---

## 1. 塑料降解酶高通量筛选体系

这里的“塑料降解酶”主要指 PETase / cutinase / polyester hydrolase / polyurethane esterase 等水解酶。不同塑料底物的筛选体系成熟度不同：PET/polyester 最成熟，PU 次之，PA/nylon 需要单独处理。

### 1.1 一级筛选 A：模型底物显色 / 荧光法

#### 原理

用可溶小分子模型底物代替真实塑料。水解后释放可见或荧光产物。

常见读出：

```text
p-nitrophenyl esters → p-nitrophenol, A405
fluorescein diacetate / fluorescein esters → fluorescein fluorescence
MpNPT / mono-p-nitrophenyl terephthalate → p-nitrophenol, 用于 MHETase-like assays
```

#### 流程

```text
1. 96/384-well 表达候选酶或突变体；
2. 加入模型酯底物；
3. 30–60 min 内读 A405 或 fluorescence；
4. 选高信号变体进入真实聚合物验证。
```

#### 优点

```text
通量高；
成本低；
适合快速排除完全无活性或表达失败的变体；
适合初筛 esterase-like 活性。
```

#### 缺点

```text
模型底物太小，不能代表真实塑料表面；
容易筛到普通 esterase，而不是 polymer-active enzyme；
无法反映聚合物吸附、链段进入、结晶度和产物释放；
对 PET/PA 真实解聚能力相关性有限。
```

#### 在本课题中的定位

```text
只作为表达/基本催化能力预筛，不能作为最终活性标准。
```

---

### 1.2 一级筛选 B：透明圈 / halo assay

#### 原理

把乳化或纳米颗粒形式的聚酯、PU 或其他可浑浊化聚合物加入 agar plate。分泌型酶或菌落周围发生降解后产生透明圈。

#### 流程

```text
1. 制备含 polymer emulsion / nanoparticles 的 agar plate；
2. 点接菌株、文库克隆或纯化酶；
3. 培养或孵育；
4. 观察 clear zone / halo；
5. 挑选 halo 最大的克隆进入二级验证。
```

#### 优点

```text
可直接用于菌落级筛选；
不需要复杂仪器；
适合发现分泌型 polyesterase / cutinase / PU esterase。
```

#### 缺点

```text
主要反映浑浊度变化，不一定等于化学解聚；
halo 受酶扩散、表达、分泌、聚合物颗粒大小影响很大；
定量性差；
对高结晶 PET、PA6、PA66 等真实材料相关性有限。
```

#### 在本课题中的定位

```text
适合作为 metagenomic discovery 或 secretion library 的粗筛，不能作为最终定量活性数据。
```

---

### 1.3 一级 / 二级筛选 C：真实聚合物微孔板 + LC/UPLC/HPLC 定量

#### 原理

直接用 PET film、PET powder、PET nanoparticles 或真实塑料片作为底物。反应后用 HPLC/UPLC/LC-MS 定量 TPA、MHET、BHET、EG 或其他单体/寡聚体产物。

#### 流程

```text
1. 将标准化 PET film / powder / textile particles 放入 96-well plate；
2. 加入纯化酶、粗酶液或细胞裂解液；
3. 在目标温度反应数小时到数天；
4. 离心或过滤去除固体；
5. 用 HPLC/UPLC/LC-MS 定量产物；
6. 根据产物总量和产物分布排名。
```

#### 优点

```text
底物真实；
能直接测化学产物；
能区分 TPA/MHET/BHET 等产物；
适合最终验证和小规模文库筛选。
```

#### 缺点

```text
通量低于显色法；
仪器和样品处理成本高；
对上千到上万突变体不友好；
需要严格标准化 film 厚度、粒径、结晶度和预洗步骤。
```

#### 文献状态

已有 PETase directed evolution 使用真实聚合物底物并通过 UPLC 定量降解产物；RWTH AMIDE 论文在引言中也指出，当时成功使用 authentic polymer 进行 directed evolution 的典型案例是 HotPETase，其 variant fitness 通过 UPLC 定量产物获得。

#### 在本课题中的定位

```text
二级验证和最终实验标签来源。
```

---

### 1.4 高通量 MS / 自动化质谱筛选

#### 原理

用自动进样或液滴系统把 96/384-well 反应液快速送入 MS，直接检测真实产物。

对 PET/polyester 可检测：

```text
TPA
MHET
BHET
oligomeric PET hydrolysis products
```

对 PA/nylon 可检测：

```text
PA6 Ahx oligomers
PA66 L1/L2/C1/C2
```

#### 优点

```text
化学特异性高；
可同时看多个产物；
比 HPLC/UPLC 快；
适合产物选择性筛选。
```

#### 缺点

```text
需要 MS 平台和自动化；
基质效应明显；
需要内标和标准品；
不是所有实验室容易建立。
```

#### 在本课题中的定位

```text
对课题一：用于二级确认新骨架真实产物。
对课题二：用于 PA6/PA66 产物分布筛选，尤其是区分 L1/L2/Ahx oligomers。
```

---

### 1.5 稳定性预筛：DSF / nanoDSF / heat challenge

#### 原理

塑料降解通常需要高温以提高聚合物链段运动性，因此热稳定性是必须筛的性状。

流程：

```text
1. 表达或纯化候选酶；
2. 用 DSF/nanoDSF 测 Tm；
3. 或在目标温度预孵育后测残余活性；
4. 过滤低稳定性候选。
```

#### 优点

```text
通量高；
可快速排除无法在 60–80 °C 工作的酶；
适合和活性筛选并行。
```

#### 缺点

```text
热稳定性不等于塑料降解活性；
有些高 Tm 酶仍然不能形成 productive binding。
```

#### 在本课题中的定位

```text
作为活性筛选的并行过滤条件，而不是主筛选信号。
```

---

## 2. 尼龙酶高通量筛选体系

尼龙酶筛选和 PET 不同。PA6/PA66 水解释放的是胺类或含胺寡聚体，因此可以用总伯胺显色，也可以用 MS 直接看 L1/L2/Ahx 分布。

---

### 2.1 AMIDE / MAFC 总伯胺显色筛选

#### 文献体系

RWTH Aachen 的 AMIDE assay 是目前最明确的尼龙酶高通量 directed evolution 平台。它用 Meldrum’s acid furfural conjugate, MAFC, 和 PA/PUR 水解释放的伯胺反应，形成 494 nm 可读的显色产物。该体系报告可检测 PA6、PA66 和 PU 降解产物，并被用于 NylC directed evolution。文献中用真实 PA6 film 作为底物，筛选 1700 个 NylC random mutagenesis clones，获得 turnover frequency 提高 1.9 倍的 NylC_TS^P27Q/F301L。

#### 实验流程

```text
1. 构建 NylC_TS error-prone PCR 文库；
2. 96-well E. coli BL21(DE3) 表达；
3. 裂解细胞，取 clarified lysate；
4. 每孔加入标准化 PA6 film disc；
5. 60 °C, 4 h 反应；
6. 离心去除固体和变性蛋白；
7. 取 100 μL 上清 + 100 μL MAFC solution；
8. 37 °C, 30 min 显色；
9. 读 A494；
10. 高信号克隆进入 HPLC/LC-MS 验证。
```

#### 关键条件

```text
PA6 film：0.2 mm thick，6.4 mm diameter discs，约 7.6 mg/disc
反应体积：200 μL
底物浓度：约 38 g/L PA6 film
buffer：50 mM Bicine, 100 mM NaCl, pH 9.0
酶源：10 μL clarified cell lysate
反应：60 °C, 800 rpm, 4 h
MAFC：4 mM in ethanol
读数：A494
```

#### 优点

```text
使用真实 PA6 film；
96-well 格式；
适合 directed evolution；
可用粗裂解液；
成本低于 LC-MS；
能快速筛出总水解释放增强的变体。
```

#### 缺点

```text
测的是总伯胺，不区分 Ahx1/Ahx2/Ahx3/Ahx4 或 PA66 L1/L2；
不能检测 caprolactam；
不能直接判断 endo/exo 机制；
可能受细胞裂解液背景和其他胺干扰；
hit 必须用 HPLC/LC-MS 验证产物分布。
```

#### 在本课题中的定位

```text
尼龙酶一级初筛。
尤其适合：
- NylC/Nyl50/new scaffold 文库；
- 口袋/通道工程文库；
- 酸预处理寡聚体水解初筛；
- 耐盐/耐残余酸环境活性筛选。
```

---

### 2.2 I.DOT/OPSI-MS 高通量产物谱筛选

#### 文献体系

ORNL 建立了 PAL/I.DOT 液体处理系统与 open port sampling interface mass spectrometry, OPSI-MS, 结合的方法，用于 PA6/PA66 hydrolysis 的高通量产物检测。该方法可在 96-well plate format 中以约 8–20 s/sample 的速度检测 PA6/PA66 水解产物，并能区分 PA6 linear dimer、PA66 linear monomer、cyclic oligomers 等产物。

#### 实验流程

```text
1. PA6/PA66 pellet 或 powder 与 NylC 或候选酶反应；
2. 反应后简单稀释或甲醇/甲酸淬灭；
3. 反应液装入 96-well plate；
4. PAL 或 I.DOT 自动取样；
5. OPSI 将 nL–μL 级样品引入 MS；
6. MS 直接检测 PA6/PA66 oligomer masses；
7. 用内标和标准品定量 L1/L2/Ahx oligomers。
```

#### 可检测产物

```text
PA6:
caprolactam
cyclic dimer/trimer/tetramer...
linear dimer L2
linear trimer L3

PA66:
adipic acid
hexamethylene diamine
cyclic monomer C1
cyclic dimer C2
linear monomer L1
linear dimer L2
```

#### 优点

```text
化学特异性高；
可同时看多个产物；
能区分 PA6 和 PA66 产物分布；
速度显著高于 HPLC；
适合 100–1000 个样品的产物谱筛选。
```

#### 缺点

```text
需要 MS 和自动液滴平台；
标准品不一定齐全；
对大寡聚体或离子化效率差的产物可能有偏差；
不是所有实验室都容易搭建。
```

#### 在本课题中的定位

```text
尼龙酶二级筛选和产物选择性筛选。
特别适合回答：
- Nyl50-like 酶是否偏 PA66 L2？
- NylC-like 酶是否偏 PA66 L1 或 PA6 Ahx2？
- 口袋/通道突变是否改变产物分布？
- 酸低聚化产物经酶处理后是否向目标产物收敛？
```

---

### 2.3 LC-MS / HPLC 二级确认

#### 原理

用标准 LC-MS/HPLC 方法对 hit 做定量验证。

#### 适合检测

```text
PA6:
Ahx1 / 6-AHA
Ahx2
Ahx3
Ahx4
Ahx5
cyclic dimer/trimer
caprolactam, if relevant

PA66:
L1
L2
C1
C2
HMD
adipic acid
longer oligomers, if standards and method allow
```

#### 优点

```text
定量可靠；
适合最终报告数据；
可做标准曲线；
可作为模型训练的实验标签。
```

#### 缺点

```text
通量低；
样品处理复杂；
标准品和分离方法要求高。
```

#### 在本课题中的定位

```text
所有 AMIDE 或 OPSI-MS hits 的最终确认。
```

---

### 2.4 GPC/SEC、DSC/WAXS、SEM/AFM 材料级验证

高通量产物筛选只能说明有可溶产物，不一定说明真实聚合物主链大量断裂。

因此最终需要材料层面验证：

```text
GPC/SEC：残余 polymer Mw 是否下降；
DSC/WAXS：结晶度是否改变；
SEM/AFM：表面是否发生侵蚀或粗糙化；
mass loss：总失重；
tensile testing：力学性能变化，可选。
```

特别注意：

```text
材料碎裂 ≠ 化学解聚；
总胺释放 ↑ ≠ 产物分布理想；
低聚物预存污染必须 blank wash 扣除。
```

---

## 3. 建议为两个课题建立的筛选平台

### 3.1 课题一：通用塑料降解酶平台

建议三层筛选：

```text
一级：模型底物 / 小分子快速筛
- pNP esters
- fluorescein esters
- soluble oligomers
- 快速排除无表达或无基本水解能力的设计

二级：真实或半真实聚合物底物
- PET/PET-like powder, film, nanoparticles
- PU model polymer / polyester-PU substrate
- HPLC/UPLC/LC-MS 定量产物

三级：材料和产物验证
- LC-MS product spectrum
- GPC/SEC
- SEM/AFM
- heat stability / residual activity
```

对于你课题一的计算模型，实验标签应分层：

```text
fast label：模型底物活性
mechanistic label：真实产物 LC-MS/HPLC
application label：材料转化率 / product distribution
```

### 3.2 课题二：尼龙酶平台

建议三层筛选：

```text
一级：AMIDE / MAFC 总伯胺筛选
- 适合大文库
- 用真实 PA6 film 或酸低聚化 PA6/PA66 寡聚体
- 输出 total amine release

二级：I.DOT/OPSI-MS 或 LC-MS 产物谱
- 区分 PA6 Ahx1–Ahx5
- 区分 PA66 L1/L2/C1/C2
- 判断 product selectivity

三级：材料级验证
- GPC/SEC 看 Mw 是否下降
- DSC/WAXS 看结晶度
- SEM/AFM 看表面侵蚀
- blank wash 扣除预存低聚物
```

---

## 4. 与计算课题的连接

高通量筛选体系不是单独的实验模块，而应直接服务于计算假设。

### 4.1 对课题一

```text
QM/MM energy surrogate 预测：
某些序列更稳定 TS-like ensemble，降低 ΔE‡。

实验验证：
1. 模型底物是否提高；
2. 真实聚合物产物是否增加；
3. 产物谱是否符合模型预测；
4. 是否出现 product trapping 或只提高 substrate binding。
```

### 4.2 对课题二

```text
MD/enhanced sampling 预测：
某个 Nyl50-like 突变降低 PA66 chain-entry barrier，偏向 L2-producing pose。

实验验证：
1. AMIDE 总胺是否增加；
2. OPSI-MS / LC-MS 中 L2/L1 比例是否改变；
3. GPC 是否显示 polymer Mw 变化；
4. 表面验证是否支持更强 productive adsorption。
```

---

## 5. 最终推荐

```text
塑料降解酶课题：
一级用模型底物和可溶低聚物快速筛；二级用真实聚合物 + LC/UPLC/MS；最终以产物谱和材料变化作为实验标签。

尼龙酶课题：
一级用 AMIDE/MAFC 总伯胺筛；二级用 I.DOT/OPSI-MS 或 LC-MS 区分 Ahx/L1/L2 产物；三级用 GPC/DSC/WAXS/SEM 验证真实 PA6/PA66 主链变化。
```

最关键原则：

```text
高通量筛选不能只追求“信号强”。
你的两个课题都需要把筛选信号拆成：

1. total activity；
2. product selectivity；
3. material-level depolymerization；
4. process compatibility；
5. 和计算预测的一致性。
```

---

## 6. 关键文献笔记

1. **AMIDE / MAFC nylon assay**  
   Puetz et al., *ACS Sustainable Chemistry & Engineering* 2023.  
   Validated high-throughput screening system for directed evolution of nylon-depolymerizing enzymes.  
   关键信息：96-well PA6 film + cell lysate + MAFC A494；筛 1700 clones；NylC_TS^P27Q/F301L turnover frequency 提高 1.9 倍。

2. **PAL/I.DOT/OPSI-MS polyamide assay**  
   Cahill et al., *Journal of the American Society for Mass Spectrometry* manuscript / ORNL public manuscript.  
   关键信息：96-well compatible；8–20 s/sample；检测 PA6 linear dimer 和 PA66 linear monomer；可同时分析 linear/cyclic oligomers。

3. **Nyl10/Nyl50 diversity panel screen**  
   Drufva et al., 2025, substrate- and product-selective nylon hydrolases.  
   关键信息：95 个 Ntn hydrolase homologs；crude lysate + washed PA6/PA66 powder；I.DOT/OPSI-MS 产物定量；Nyl10/Nyl50 PA66-selective。

4. **PET/polyester authentic polymer screening**  
   PETase/cutinase directed evolution 文献通常以真实 PET film/powder 反应后 HPLC/UPLC/LC-MS 定量 TPA/MHET/BHET 作为可靠二级或终点筛选；这种方法化学特异性强，但通量低于 plate colorimetric assays。
