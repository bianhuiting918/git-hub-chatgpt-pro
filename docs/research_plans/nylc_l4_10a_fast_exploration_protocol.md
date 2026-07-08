# NylC A/D cleft + PA66 L4：10 Å 起点快速自由探索模拟 protocol

_Last updated: 2026-07-08_

## 0. 本文件定位

本文件补充 `nylc_tetramer_l4_l8_slab_simulation_plan.md`，专门细化一个实际执行问题：

```text
L4 初始放在 NylC A/D-associated cleft 外约 8–10 Å，
如何在 50 ns 左右尽快接触表面，
同时尽量保留“自由探索”的意义，
并统计 release 阶段的 productive-like pose 比例。
```

核心原则：

```text
可以加弱引导帮助 L4 接触 A/D 表面；
但 productive-like pose fraction 只能在去掉引导后的 release MD 中统计。
```

---

## 1. 目标问题

本 protocol 不以直接得到反应路径为目标，而是回答：

```text
1. L4 从 A/D cleft 外 8–10 Å 起点出发，是否能在 50 ns 内接触 NylC 表面；
2. 接触后是否会进入 A/D-associated cleft；
3. 是否能自发形成 productive-like pose；
4. Register A 或 Register B 哪个更容易在 release MD 中保持；
5. 观察到的构象是否只是端基盐桥导致的非生产性吸附。
```

该 protocol 用于筛选，不替代后续：

```text
长时间 MD；
umbrella / metadynamics；
QM/MM。
```

---

## 2. 模型选择

### 2.1 酶结构

第一版可使用：

```text
NylC A/D dimer
```

用于快速测试 A/D-associated cleft。

后续必须迁移到：

```text
NylC A/B/C/D mature tetramer
```

因为 NylC 的真实装配环境可能影响 A/D cleft 的可及性。

### 2.2 底物结构

优先：

```text
L4-Nfree
H3N+—[U1]—[U2]—[U3]—[U4]—NHMe
```

第二组：

```text
L4-Cfree
Ac—[U1]—[U2]—[U3]—[U4]—COO-
```

对照：

```text
L4-capped
Ac—[U1]—[U2]—[U3]—[U4]—NHMe
```

### 2.3 L4 初始构象

不要只用自由 relax 后的卷曲最低能构象。每个体系建议准备：

```text
1. semi-extended L4；
2. pre-bent L4；
3. collapsed L4 作为对照。
```

第一版优先使用：

```text
semi-extended
pre-bent
```

因为链端底物在真实聚合物环境下不应完全像孤立小分子那样塌缩。

---

## 3. 起点设置

### 3.1 初始距离

把 L4 放在 A/D-associated cleft 外侧：

```text
L4 COM 到 A/D cleft center：约 8–10 Å
或 L4 最近重原子到蛋白表面：约 8–10 Å
```

不要一开始把 L4 直接塞到 Thr267 附近，否则会失去“接触/进入探索”的意义。

### 3.2 初始方向

至少设置两种链端方向：

```text
Direction 1：N-free end toward A/D cleft
Direction 2：C-capped side / opposite orientation toward A/D cleft
```

如果使用 L4-Cfree，则对应：

```text
C-free end toward A/D cleft
N-capped side / opposite orientation toward A/D cleft
```

### 3.3 Register 不强制作为初始约束

如果目标是看“自由探索是否能产生 productive-like pose”，第一阶段不要强制 scissile bond 到 Thr267。

可以在后处理时识别：

```text
Register A：U1–U2 靠近 Thr267；
Register B：U2–U3 靠近 Thr267；
Register C：U3–U4 靠近 Thr267。
```

---

## 4. 为什么不直接用普通 MD

如果完全无引导，50 ns 内 L4 可能不会有效接触 A/D cleft。常见情况：

```text
1. L4 在水中卷曲；
2. L4 漂离 cleft；
3. L4 贴到非活性表面；
4. L4 与蛋白短暂碰撞后离开；
5. 即使靠近 Thr，也不一定是正确 register。
```

因此第一版推荐：

```text
weak encounter guidance
+
release MD
```

而不是：

```text
无引导直接跑 50–200 ns。
```

---

## 5. 推荐 protocol A：弱引导接触 + release MD

这是第一版最推荐方案。

### 5.1 定义两个 groups

```text
Group 1：L4 COM
可用 L4 全部重原子，或 L4 中部重原子。

Group 2：A/D cleft center
可用 Thr267 附近 pocket residues 的 COM，或 A/D cleft 入口残基 COM。
```

不要在第一阶段约束：

```text
Thr267 Oγ — scissile carbonyl C
```

因为这会强迫反应构象，污染 productive fraction。

### 5.2 三阶段 50 ns protocol

```text
Stage 1：0–5 ns，moving restraint
COM_L4 到 A/D cleft center：10 Å → 6 Å
目的：让 L4 快速接触 A/D 表面。

Stage 2：5–10 ns，flat-bottom surface restraint
COM_L4 到 A/D cleft center 保持在 4–8 Å 区间。
目的：让 L4 在表面附近自由转动和寻找接触方式。

Stage 3：10–50 ns，release MD
去掉 encounter restraint。
目的：统计无引导条件下 productive-like pose fraction。
```

### 5.3 restraint 强度

推荐弱约束：

```text
1–5 kcal/mol/Å²
```

如果用 SI 单位或 GROMACS，需要转换。关键原则：

```text
只让 L4 靠近 surface，
不强迫特定 scissile amide bond 到 Thr267。
```

---

## 6. 推荐 protocol B：soft confinement + 轻微接触引导

如果希望更接近自由探索，可以使用 soft wall / spherical confinement。

### 6.1 soft confinement

定义：

```text
L4 COM 不允许离 A/D cleft center 超过 12–15 Å。
```

作用：

```text
防止 L4 漂到水盒远处；
但不规定 L4 必须以什么方向接触 cleft。
```

### 6.2 组合策略

```text
0–5 ns：轻微 COM 接触引导，让 L4 从 10 Å 靠近到 6–8 Å；
5–50 ns：只保留 soft confinement 或完全释放。
```

这种方案比直接 Thr-distance restraint 更适合估计“附近自由探索”的 productive-like pose 比例。

---

## 7. 温度和扩散加速

可以略升温来增加 L4 构象搜索速度，但不要过高。

推荐比较：

```text
300 K：基准；
323 K：轻度加速；
333 K：高温探索。
```

注意：

```text
如果 productive-like pose 只在 333 K 出现，而 300–323 K 不出现，结论要谨慎。
```

不建议第一版直接使用很高温度，例如：

```text
>350 K
```

因为可能引入非真实构象、蛋白 loop 异常波动或底物过度塌缩/展开。

---

## 8. 多起点与 replicate 设计

一条 50 ns 轨迹不够。推荐多起点、多速度 replicate。

第一版最小组：

```text
底物：L4-Nfree
构象：semi-extended / pre-bent
方向：N-end first / opposite orientation
replicate：每种 5 条
```

总数：

```text
2 conformations × 2 orientations × 5 replicates = 20 trajectories
```

每条：

```text
5 ns moving restraint
5 ns flat-bottom surface restraint
40 ns release MD
```

累计：

```text
20 × 50 ns = 1 μs cumulative sampling
```

这通常比一条 1 μs 轨迹更适合探索初始结合与 productive-like pose 形成概率。

---

## 9. productive-like pose 统计只在 release 阶段做

不要在 Stage 1 / Stage 2 统计 productive-like fraction，因为这些阶段有引导力。

正式统计区间：

```text
Stage 3：10–50 ns release MD
```

如果后续延长：

```text
50–100 ns 或 50–200 ns 也可统计为 release extension。
```

---

## 10. productive-like pose 定义

对 NylC + L4，推荐判定标准：

```text
1. 任一 PA66 scissile amide carbonyl C 到 Thr267 Oγ 距离在 3.0–4.5 Å；
2. Thr267 Oγ 对 carbonyl C 的攻击角合理；
3. carbonyl O 指向 oxyanion-stabilizing region；
4. leaving amide N 附近存在可能质子转移路径；
5. register 在 release MD 中不频繁滑移；
6. L4 没有严重卷曲脱离 active cleft；
7. 不是只靠端基 -NH3+ / -COO- 与蛋白形成盐桥固定。
```

后处理时应分别识别：

```text
Register A productive-like：U1–U2 carbonyl C 满足几何标准；
Register B productive-like：U2–U3 carbonyl C 满足几何标准；
Register C productive-like：U3–U4 carbonyl C 满足几何标准。
```

---

## 11. 需要输出的统计表

每条轨迹至少输出：

```text
trajectory_id
enzyme_model: NylC_AD_dimer / NylC_tetramer_AD_zone
ligand_model: L4-Nfree / L4-Cfree / L4-capped
initial_conformation: semi-extended / pre-bent / collapsed
initial_orientation: N-end-first / C-end-first / opposite / random
protocol: moving_restraint + flat_bottom + release
simulation_temperature
release_time_ns
surface_contact_time_fraction
AD_cleft_contact_time_fraction
nearest_Thr267_distance_min
nearest_Thr267_distance_mean
Register_A_productive_fraction
Register_B_productive_fraction
Register_C_productive_fraction
terminal_salt_bridge_fraction
ligand_collapse_metric
classification: no-contact / surface-binding / pocket-binding / productive-like
notes
```

建议额外保存：

```text
representative productive-like snapshots；
representative surface-binding snapshots；
representative nonproductive adsorption snapshots；
L4 contact map；
Thr267 distance time series；
register assignment time series。
```

---

## 12. 如果 50 ns 内仍不能接触怎么办

如果多数 release 前都无法有效接触 A/D 表面，优先调整：

```text
1. Stage 1 moving restraint 目标从 6 Å 改为 5 Å；
2. Stage 2 flat-bottom 区间从 4–8 Å 改为 3.5–7 Å；
3. 增加 Stage 2 到 10–20 ns；
4. 增加 initial orientations；
5. 轻度升温到 323 K；
6. 使用 soft confinement 防止 L4 漂走。
```

不要直接改成：

```text
强拉 scissile carbonyl C 到 Thr267 Oγ。
```

因为这会让“自由探索 productive pose 比例”失去意义。

---

## 13. 如果只出现端基盐桥怎么办

如果 L4-Nfree 或 L4-Cfree 主要通过自由端基形成强盐桥，而非主链酰胺进入 pocket，应追加：

```text
L4-capped
```

并比较：

```text
one-end-free L4 是否只因端基电荷吸附；
capped L4 是否还能形成 pocket-binding 或 productive-like pose。
```

如果只有 one-end-free 有结合，而 capped 完全没有，需谨慎解释为：

```text
可能是 chain-end recognition；
也可能是端基盐桥假象。
```

需要通过 L8 / slab exposed chain-end 模型进一步验证。

---

## 14. 与 docking 的关系

本 protocol 是“不先做严格 docking”的自由探索版本。

如果目标变成快速生成高质量 register pose，则应改用：

```text
A/D-focused constrained docking
+
几何过滤
+
minimization
+
release MD
```

两者区别：

```text
10 Å 自由探索 protocol：
更适合估计 encounter 和 spontaneous productive-like pose 形成概率。

constrained docking protocol：
更适合快速生成 Register A/B/C 候选反应构象。
```

二者可以互补，不应互相替代。

---

## 15. 最小执行版本

如果只想今天先启动一版：

```text
Model：NylC A/D dimer
Ligand：L4-Nfree
Initial distance：8–10 Å from A/D cleft
Initial conformations：semi-extended, pre-bent
Initial orientations：N-end-first, opposite
Replicates：5 each
Temperature：323 K
Protocol per trajectory：
  0–5 ns moving restraint: COM_L4 → AD_cleft_center, 10 Å → 6 Å
  5–10 ns flat-bottom: COM_L4 within 4–8 Å of AD_cleft_center
  10–50 ns release MD: no encounter restraint
Analysis：release-only productive-like fraction
```

第一批总计：

```text
2 × 2 × 5 × 50 ns = 1 μs cumulative MD
```

---

## 16. 一句话总结

```text
为了在 50 ns 内让 L4 从 8–10 Å 外更快接触 NylC A/D cleft，
可以使用弱 COM/contact encounter restraint 或 soft confinement，
但不能用 Thr–carbonyl 强约束来统计 productive-like fraction。
真正的 productive-like 构象比例必须在去掉引导后的 release MD 阶段统计。
多起点、多 replicate、轻度升温，比一条长时间无偏轨迹更适合当前问题。
```
