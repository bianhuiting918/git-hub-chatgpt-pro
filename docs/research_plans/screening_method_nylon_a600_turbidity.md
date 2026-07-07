# 尼龙酶筛选补充方法：PA6/PA66 颗粒悬液 A600 浊度法

_Last updated: 2026-07-07_

## 1. 方法定位

PET hydrolase 文献中常用 PET nanoparticles 的 600 nm 浊度变化来做动力学分析。这个方法可以借鉴到尼龙体系，但在 PA6/PA66 中应作为**辅助筛选**，不能替代 AMIDE/MAFC 和 LC-MS/OPSI-MS。

推荐定位：

```text
PA6/PA66 颗粒 A600 浊度法
= 材料颗粒悬液形态变化的辅助读出

AMIDE/MAFC
= 总伯胺释放一级活性读出

LC-MS/OPSI-MS
= Ahx/L1/L2/C1/C2 产物谱确认

GPC/SEC、SEM、DSC/WAXS
= 材料级验证
```

核心原则：

```text
不能只凭 A600 下降判断尼龙发生化学解聚；
必须同时看到 MAFC 总胺上升和 LC-MS 中 Ahx/L1/L2 等产物增加。
```

---

## 2. 原理

如果 PA6/PA66 被制备成稳定的纳米颗粒或亚微米颗粒悬液，颗粒会散射光，因此在 600 nm 处有浊度信号。

```text
PA6/PA66 颗粒悬液
→ 颗粒散射光
→ A600 / OD600 较高
```

如果尼龙酶侵蚀颗粒表面或使颗粒变小，理论上会出现：

```text
颗粒尺寸减小
颗粒数量减少
表面被水解
部分材料转为可溶寡聚体
        ↓
悬液浊度下降
        ↓
A600 随时间下降
```

但 A600 下降也可能来自非酶原因：

```text
颗粒沉降；
颗粒团聚；
分散状态改变；
气泡或边缘效应；
温度造成悬液稳定性变化。
```

因此必须设置严格对照。

---

## 3. 推荐粒径

用于 A600 连续读数的 PA6/PA66 底物应能稳定悬浮。

推荐粒径：

```text
最佳：100–500 nm
可用：0.5–2 μm
勉强可用：2–5 μm，需要强摇晃或搅拌
不推荐：>10 μm，容易沉降，A600 不稳定
```

不建议把文献中用于真实材料酶解的几百微米 cryomilled PA6/PA66 powder 用于 A600 实时浊度法，因为颗粒太大，沉降会严重干扰读数。

---

## 4. 材料来源和制备建议

### 4.1 直接购买

如果供应商能提供：

```text
PA6 nanoparticles
PA66 nanoparticles
PA6/PA66 submicron powder
```

优先选择粒径：

```text
100–500 nm 或 0.5–2 μm
```

购买时需要记录：

```text
粒径分布；
结晶度；
是否含添加剂；
是否表面改性；
是否有低聚物残留。
```

### 4.2 自制 PA6/PA66 微粒或纳米颗粒

可采用溶解-反沉淀策略：

```text
1. 用合适溶剂溶解 PA6/PA66；
2. 快速加入反溶剂中形成颗粒；
3. 离心或透析洗去残余溶剂；
4. 超声分散；
5. 用 DLS / 激光粒度 / SEM 测粒径；
6. 调整悬液浓度，使初始 A600 约 0.2–1.0。
```

注意：

```text
残余甲酸、HFIP、TFE 或其他溶剂会影响酶活，必须充分去除。
```

---

## 5. 推荐实验流程

```text
1. 制备 PA6/PA66 颗粒悬液；
2. 超声分散，保证无明显大团聚；
3. 调整初始 A600 到 0.2–1.0；
4. 分装到 96-well plate；
5. 加入纯化酶或 crude lysate；
6. 设定目标温度，例如 50–70 °C；
7. 600 nm 连续读数，记录 A600(t)；
8. 同一反应体系取样做 MAFC 和 LC-MS；
9. 反应后用 SEM/DLS 检查颗粒尺寸变化。
```

建议使用：

```text
连续摇板；
同一孔内多点读数；
边缘孔不用或加 buffer；
每组至少 triplicate；
反应前后均测 DLS 或 SEM。
```

---

## 6. 必须设置的对照

```text
1. no enzyme control
只有 PA6/PA66 颗粒，无酶。
用于判断自然沉降或温度导致的 A600 变化。

2. heat-denatured enzyme control
热失活酶。
用于判断蛋白或裂解液对颗粒分散状态的非催化影响。

3. buffer-only control
无颗粒、无酶。
用于扣背景。

4. active enzyme + soluble oligomer positive control
确认酶本身有活性。

5. MAFC parallel assay
同一时间点检测总伯胺释放。

6. LC-MS / OPSI-MS confirmation
确认 Ahx1–Ahx5 或 PA66 L1/L2/C1/C2 产物。

7. particle-size control
反应前后 DLS/SEM，区分水解和团聚/沉降。
```

---

## 7. 数据解释

### 7.1 有效降解的理想模式

```text
A600 下降
+
MAFC A494 上升
+
LC-MS 检测到 Ahx/L1/L2 产物增加
+
DLS/SEM 显示颗粒变小或表面侵蚀
```

这可以支持：

```text
颗粒材料发生酶促水解和形态变化。
```

### 7.2 不能单独解释为降解的情况

```text
只有 A600 下降，但 MAFC 和 LC-MS 不变
```

可能原因：

```text
颗粒沉降；
颗粒团聚；
蛋白改变悬液稳定性；
气泡或光学伪影。
```

```text
MAFC 上升但 A600 不变
```

可能说明：

```text
有少量表面水解，但材料整体颗粒形态变化不明显。
```

```text
LC-MS 有 L1/L2，但 GPC/SEM 无变化
```

可能说明：

```text
主要处理预存低聚物或表面少量可及链段，不代表大量高分子主链断裂。
```

---

## 8. 与现有筛选体系的组合

推荐组合：

```text
一级主筛：AMIDE/MAFC A494
辅助主筛：PA6/PA66 particle A600
二级确认：LC-MS/OPSI-MS
三级验证：GPC/SEC、SEM、DSC/WAXS
```

实际筛选流程：

```text
PA6/PA66 颗粒悬液
        ↓
加 NylC/Nyl10/Nyl50/new scaffold
        ↓
同时记录 A600(t)
        ↓
定时取样做 MAFC A494
        ↓
Top hits 做 LC-MS/OPSI-MS
        ↓
Top hits 做 GPC/SEM/DSC/WAXS
```

---

## 9. 对尼龙酶课题的具体用途

A600 浊度法适合回答：

```text
1. Nyl50/Nyl10 是否比 NylC 更能改变 PA66 颗粒悬液；
2. NylC-V3 或 NylC-HP 是否改善 PA66 颗粒材料侵蚀；
3. 界面突变是否改变酶-颗粒表面相互作用；
4. 酸预处理或低结晶颗粒是否更容易出现 A600 下降；
5. A600 下降是否与 L2/L1 产物比例相关。
```

但它不能单独回答：

```text
1. 产物是 L1 还是 L2；
2. 是 chain-end cleavage 还是 internal-loop cleavage；
3. 是否有大量主链内切；
4. 是否可以再聚合回 PA66。
```

这些仍需 LC-MS/OPSI-MS 和 GPC/SEC。

---

## 10. 采购和制备建议

### 优先购买

```text
PA6 film，0.2 mm，用于 AMIDE 标准体系；
additive-free PA6 powder，用于真实材料验证；
additive-free PA66 pellets/powder，用于 Nyl10/Nyl50 筛选；
6-AHA，PA6 产物标准；
HMD 和 adipic acid，PA66 单体标准。
```

### A600 浊度法专用

```text
PA6 nanoparticles / submicron particles，100–500 nm 或 0.5–2 μm；
PA66 nanoparticles / submicron particles，100–500 nm 或 0.5–2 μm；
或自制 PA6/PA66 溶解-反沉淀颗粒。
```

### 不推荐直接用于 A600 的材料

```text
5–50 μm PA6 powder：可用于反应，但 A600 稳定性一般；
几百 μm cryomilled PA6/PA66 powder：适合 LC-MS 真实材料反应，不适合浊度动力学；
film/disc：适合 MAFC/HPLC，不适合悬液浊度。
```

---

## 11. 最终定位

```text
PA6/PA66 颗粒 A600 浊度法可以作为 PET-NP turbidity assay 的尼龙版本，
用于辅助判断颗粒材料是否发生形态变化。

但尼龙酶筛选不能只用 A600。
最稳妥的体系是：
A600 + MAFC + LC-MS/OPSI-MS + GPC/SEM。
```
