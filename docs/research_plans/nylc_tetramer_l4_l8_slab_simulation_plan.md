# NylC 四聚体模拟补充计划：L4 多方向 register 扫描与 L8 / PA66 slab 吸附-链端进入模拟

_Last updated: 2026-07-08_

## 0. 本文件定位

本文件补充 `project_2_nylon_two_track_strategy.md` 中的代表性酶机制模拟部分，专门细化两个方向：

```text
方向 A：NylC 四聚体 + PA66 L4 多方向 / 多 register 扫描
目的：判断 NylC 四聚体中 A/D、A/B 或其他表面区域，哪个能形成真正有反应意义的 productive-like pose。

方向 B：NylC 四聚体 + tethered L8 / PA66 slab 多取向模拟
目的：判断 NylC 四聚体面对 PA66 材料表面时，倾向用哪个区域吸附，并判断表面链端是否有机会进入 active-site cleft。
```

这里的核心区分是：

```text
L4 = register 判别模型
重点比较不同界面和不同 register 是否能形成反应构象。

L8 / slab = 长链段进入和材料表面吸附模型
重点看四聚体在 PA66 表面的吸附取向，以及链端是否能进入口袋并形成预反应构象。
```

---

## 1. 为什么 NylC 主模型应使用四聚体

文献支持 NylC / NylC-GYAQ 在相关浓度下可以形成 A/B/C/D tetramer，其中 A/B interface 面积大，是主要稳定装配界面。但 putative active-site / substrate-access 相关区域更可能跨 A/D interface。

因此，NylC 模拟不应只取 A/B dimer 或 A/D dimer。第一优先级应使用：

```text
NylC mature A/B/C/D tetramer
```

在四聚体中同时区分：

```text
A/B interface：structural assembly interface，主要负责装配稳定。
A/D-associated cleft：functional active-site / substrate-access candidate。
其他 monomer-local pocket 或 surface groove：探索性对照。
```

关键判断不是“哪个界面面积最大”，而是：

```text
哪个界面能让 PA66 可切酰胺键靠近 catalytic Thr，
并保持合理亲核攻击几何和 oxyanion stabilization。
```

---

## 2. 方向 A：NylC 四聚体 + L4 多方向 / 多 register 扫描

### 2.1 目标

L4 是机制判别模型，用于回答：

```text
1. NylC 四聚体中 A/D-associated cleft 是否比 A/B interface 更容易形成 productive-like pose；
2. NylC 是否偏 L1-producing Register A；
3. L2-producing Register B 在 NylC 中是否不稳定；
4. A/B interface 是否只是材料/底物非生产性吸附界面；
5. 端基电荷是否造成人工吸附假象。
```

L4 写作：

```text
[U1]—[U2]—[U3]—[U4]
```

其中 U 是一个 PA66 repeat unit。

---

### 2.2 L4 端基版本

第一版优先：

```text
L4-Nfree：胺端自由，一端封端
H3N+—[U1]—[U2]—[U3]—[U4]—NHMe
代表 PA66 胺端链端。

L4-Cfree：羧酸端自由，一端封端
Ac—[U1]—[U2]—[U3]—[U4]—COO-
代表 PA66 羧酸端链端。
```

第二优先级：

```text
L4-capped：两端封端
Ac—[U1]—[U2]—[U3]—[U4]—NHMe
代表 PA66 内部链段。
```

第三优先级：

```text
L4-free-free：两端自由
H3N+—[U1]—[U2]—[U3]—[U4]—COO-
代表完全游离可溶寡聚体。
```

---

### 2.3 L4 register 定义

每个 L4 端基版本都应构建三类 register：

```text
Register A：切 U1–U2
[U1] | [U2]—[U3]—[U4]
预测释放 L1。

Register B：切 U2–U3
[U1]—[U2] | [U3]—[U4]
预测释放 L2。

Register C：切 U3–U4
[U1]—[U2]—[U3] | [U4]
作为更深切割对照。
```

NylC 的预期主假设：

```text
A/D-associated cleft + Register A 更可能形成 productive-like pose。
```

对照假设：

```text
A/B interface 即使能吸附 L4，也不一定能形成 productive-like pose。
```

---

### 2.4 L4 初始放置位置

不要简单按“前后左右”定义，而应按 NylC 四聚体的结构区域定义：

```text
Zone 1：A/D-associated active-site cleft
主假设，可能是底物进入/反应相关区域。

Zone 2：B/C 或等价 A/D-like active-site cleft
用于检查四聚体中对称等价区域。

Zone 3：A/B structural interface
大面积稳定装配界面，作为对照。

Zone 4：C/D structural interface
A/B 的对称或类似结构装配界面对照。

Zone 5：monomer-local active-site opening
判断是否存在不依赖界面的单体口袋入口。

Zone 6：其他 surface groove / charged patch
探索非预设吸附位点。
```

第一版最小测试只做：

```text
A/D-associated cleft × Register A/B/C
A/B structural interface × Register A/B/C
```

---

### 2.5 是否同时放 4 个 L4

第一版不建议在同一个体系里同时放 4 个 L4。

原因：

```text
1. 多个 L4 之间可能互相接触、聚集或竞争；
2. 会改变局部离子和氢键环境；
3. 分析时难以判断结合事件是否独立；
4. 容易把“多底物竞争”与“单底物 productive-like pose”混在一起。
```

推荐：

```text
一个体系只放一个 L4；
但做多个体系，每个体系把 L4 放在不同 zone 和不同 register。
```

多 L4 竞争体系可以作为后续补充，用于探索不同表面区域的吸附竞争，但不能作为 productive pose 的主证据。

---

### 2.6 如何加引导让结合更快

可以加引导，但引导只能用于采样，不能作为最终证据。最终必须做 release MD，即去掉关键引导后，看 L4 是否仍能保持 productive-like pose。

推荐引导方式：

#### 方式 1：constrained docking / 手工 register pose

对每个 register，构建初始反应构象：

```text
Thr Oγ 到 scissile amide carbonyl C：约 3.0–4.0 Å；
carbonyl O 朝向 oxyanion-stabilizing region；
L4 链段沿 pocket / cleft / groove 放置；
避免严重 clash。
```

适合第一轮 L4 register scan。

#### 方式 2：flat-bottom restraint

在短 MD 平衡中设置弱距离约束：

```text
CV = distance(Thr Oγ, scissile carbonyl C)
flat-bottom 区间：3.0–5.0 Å 不加力
超过 5.0 Å 后轻微拉回
```

作用：

```text
防止 L4 初始阶段立即跑掉；
允许口袋侧链和 L4 局部调整；
不把底物硬锁死。
```

#### 方式 3：steered MD

如果从表面到口袋的进入路径未知，可做轻量 steered MD：

```text
从表面 position 拉向 active-site cleft；
CV 可用 L4 中心到 cleft center 的距离，或 scissile carbonyl C 到 Thr Oγ 的距离。
```

steered MD 只用于生成路径候选。正确用法：

```text
1. 选取多个中间构象；
2. 去掉拉力；
3. 分别做 release MD；
4. 判断哪些构象能自然稳定。
```

#### 方式 4：umbrella sampling / metadynamics

只对少数关键体系做，例如：

```text
A/D cleft + Register A
A/D cleft + Register B
A/B interface + Register A
A/B interface + Register B
```

CV 可选：

```text
Thr Oγ 到 scissile carbonyl C 距离；
L4 插入口袋深度；
L4 与 active-site cleft 接触数；
U1/U2 subpocket occupancy。
```

---

### 2.7 L4 推荐第一版流程

```text
Step 1：准备 NylC mature tetramer。

Step 2：选择 L4-Nfree 作为第一版底物。

Step 3：构建 8 组初始体系：
1. A/D cleft + Register A
2. A/D cleft + Register B
3. A/D cleft + Register C
4. A/B interface + Register A
5. A/B interface + Register B
6. A/B interface + Register C
7. random surface placement 1
8. random surface placement 2

Step 4：constrained placement / local docking。

Step 5：minimization + side-chain repacking。

Step 6：1–5 ns weak restraint equilibration。

Step 7：20–50 ns release MD，必要时延长到 100 ns。

Step 8：统计 productive-like pose 指标。
```

后续再追加：

```text
L4-Cfree；
L4-capped；
更多 zone；
多 replicate；
少量 umbrella / metadynamics。
```

---

### 2.8 L4 判定标准

不要只看 L4 是否停留在蛋白表面。分三级判定：

```text
Level 1：surface-binding
L4 是否稳定接触某个 zone。

Level 2：pocket-binding
L4 可切酰胺键是否进入 active-site cleft。

Level 3：productive-like pose
是否形成接近反应前构象。
```

productive-like pose 指标：

```text
1. scissile amide carbonyl C 到 Thr267 Oγ：约 3.0–4.0 Å；
2. Thr Oγ 对 carbonyl C 的攻击角合理；
3. carbonyl O 指向 oxyanion-stabilizing region；
4. leaving amide N 附近存在可能质子转移路径；
5. register 在 release MD 中保持；
6. L4 没有严重扭曲；
7. 不是只靠端基 -NH3+ / -COO- 与蛋白盐桥固定。
```

核心输出表：

```text
system_id
zone
register
terminus_model
restraint_protocol
productive_like_frame_fraction
ThrOgamma_to_carbonylC_distance_mean
attack_angle_mean
oxyanion_Hbond_occupancy
register_stability
terminal_salt_bridge_fraction
L4_residence_time
classification: non-binding / surface-binding / pocket-binding / productive-like
```

---

## 3. 方向 B：NylC 四聚体 + L8 / PA66 slab 多取向模拟

### 3.1 目标

L8 / slab 不是用来一开始做化学反应，也不是直接判断产物比例。它的主要目标是：

```text
1. 判断 NylC 四聚体面对 PA66 材料表面时倾向用哪个区域吸附；
2. 判断 A/D cleft、A/B interface 或其他 zone 是否更常接触 PA66；
3. 判断 PA66 表面链端是否能进入 active-site cleft；
4. 判断是否形成接近 L4 推断的 productive-like pose；
5. 区分 productive chain-entry 和 nonproductive adsorption。
```

---

### 3.2 L8 / slab 的层级

#### 模型 1：free / tethered L8

```text
free L8：可溶长寡聚体对照。
one-end-free / tethered L8：模拟 PA66 表面链端，主模型。
```

示意：

```text
自由链端
  ↓
[U1]—[U2]—[U3]—[U4]—[U5]—[U6]—[U7]—[U8]—tether
```

用途：

```text
验证 L4 推断出的 preferred register 是否能在长链段中形成。
```

#### 模型 2：PA66 amorphous slab

```text
NylC tetramer + PA66 amorphous slab
```

用途：

```text
看 NylC 四聚体整体面对 PA66 材料时更倾向用哪个表面区域吸附。
```

#### 模型 3：PA66 slab + exposed chain end

```text
PA66 slab — [U8]—[U7]—...—[U2]—[U1]
                                      ↑
                                  可进入酶的 surface chain end
```

用途：

```text
把材料表面吸附和 chain-entry 连接起来，
判断表面链端是否能进入 A/D cleft 并形成 productive-like pose。
```

---

### 3.3 slab 模拟应跑很多次

slab 体系高度依赖初始取向，因此不能只跑一条长轨迹。推荐 multi-start simulations：

```text
同一个 NylC tetramer
+
同一个或多个 PA66 slab
+
不同酶取向
+
多个随机速度 replicate
+
多条短 MD
```

原因：

```text
1. 酶初始朝向会决定最先接触哪个 surface zone；
2. PA66 slab 表面不均一；
3. PA66 链段进入口袋是低概率事件；
4. 普通 MD 容易卡在 nonproductive adsorption；
5. 多条短轨迹比一条超长轨迹更适合筛选吸附偏好。
```

第一版推荐：

```text
8 个初始取向 × 2–3 个随机速度 replicate = 16–24 条轨迹
```

每条：

```text
20–50 ns：快速初筛吸附；
50–200 ns：延长有价值轨迹。
```

---

### 3.4 slab 初始取向设计

至少设置以下取向：

```text
Orientation 1：A/D-associated active-site cleft facing slab
测试功能性 cleft 是否能直接面向材料。

Orientation 2：A/B structural interface facing slab
测试大装配界面是否主导非生产性吸附。

Orientation 3：monomer active-site opening facing slab
测试单体口袋入口是否能接触材料。

Orientation 4：opposite/back side facing slab
负向对照。

Orientation 5–8：random orientations
测试非预设吸附偏好。
```

对 slab + exposed chain end，优先做：

```text
A/D cleft facing exposed chain end；
A/B interface facing exposed chain end；
random orientation control。
```

---

### 3.5 slab 模拟的判定标准

不能只看哪个区域“贴得最多”。应分为两层：

```text
第一层：adsorption orientation
判断 NylC 四聚体哪个 zone 倾向接触 PA66 surface。

第二层：chain-entry / productive-like potential
判断表面链段是否能靠近 catalytic Thr 并形成预反应构象。
```

吸附指标：

```text
enzyme-slab contact count；
enzyme-slab contact area；
每条 chain 的接触比例；
每个 zone 的接触比例；
接触残基 frequency；
蛋白取向角；
active cleft 到 slab 的距离；
active cleft 是否朝向 slab。
```

A/D vs A/B 指标：

```text
A/D cleft residues 与 slab 的接触频率；
A/B interface residues 与 slab 的接触频率；
Thr267 到最近 PA66 链段距离；
A/D cleft 是否被 PA66 链段接近；
A/B interface 是否只是表面贴附。
```

productive-like 指标：

```text
最近 PA66 amide carbonyl C 到 Thr267 Oγ 的最小距离；
是否有 amide carbonyl 进入 3–5 Å 范围；
carbonyl O 是否朝向 oxyanion region；
链段是否沿 pocket / cleft 方向放置；
是否出现可延续到 L4/L8 的 register。
```

nonproductive adsorption 指标：

```text
PA66 大面积贴在蛋白非活性表面；
active cleft 远离 slab；
PA66 链段稳定接触但没有任何可切键靠近 Thr；
蛋白被 slab 表面吸附后取向锁死。
```

---

### 3.6 slab 输出表

每条轨迹至少输出：

```text
trajectory_id
orientation_id
replicate_id
slab_model: amorphous / exposed_chain_end
dominant_contact_zone
contact_fraction_AD
contact_fraction_AB
contact_fraction_monomer_pocket
contact_fraction_other_surface
Thr267_to_nearest_PA66_amide_min_distance
active_cleft_facing_slab: yes/no
productive_like_event: yes/no
nonproductive_adsorption: yes/no
final_orientation_cluster
notes
```

解释示例：

```text
如果 A/B 接触频率高，但 Thr267 始终远离 PA66 链段，
说明 A/B 可能是 nonproductive adsorption surface。

如果 A/D 接触频率不一定最高，但一旦接触更容易形成 Thr-near pose，
说明 A/D 更可能是 functional substrate-entry route。
```

---

### 3.7 是否需要增强采样

第一版不一定需要。先做多起点短 MD。

如果出现以下情况，再做增强采样：

```text
酶能吸附 slab；
A/D cleft 偶尔接近表面链段；
但 chain end 很难自然进入口袋。
```

增强采样应针对：

```text
surface chain end → A/D cleft
```

而不是让整个 slab 随机探索。

CV 可选：

```text
CV1：chain-end amide carbonyl C 到 Thr267 Oγ 距离；
CV2：chain-end 插入 A/D cleft 深度；
CV3：enzyme-slab orientation angle；
CV4：PA66-enzyme contact count。
```

方法：

```text
umbrella sampling；
metadynamics；
steered MD followed by release MD。
```

---

## 4. 推荐总顺序

不要一开始直接做最复杂的 slab。推荐顺序：

```text
Step 1：NylC tetramer apo 检查
标记 A/D、A/B、monomer pocket、surface groove。

Step 2：L4 多方向 / 多 register 扫描
找出哪个 zone 和 register 能形成 productive-like pose。

Step 3：tethered L8 chain-entry
验证长链段是否能复现 L4 推断的 productive-like pose。

Step 4：PA66 amorphous slab 多取向吸附
判断 NylC 四聚体整体在 PA66 surface 上倾向怎么吸附。

Step 5：PA66 slab + exposed chain end
验证材料表面链端是否能进入 A/D cleft 或其他 functional route。

Step 6：少数关键体系做增强采样 / QM/MM
只对 L4 或 L2 的关键 productive pose 做 QM/MM；L8/slab 主要做到预反应构象即可。
```

---

## 5. 与 Nyl50 / Nyl10 的连接

虽然本文件以 NylC 四聚体为主，但同样逻辑可以迁移：

```text
Nyl50：
优先使用 A/D-like dimer，重点看 two-subpocket tunnel 和 L2-producing register。

Nyl10：
没有确定 A/B 或 A/D 标准界面，先用 Nyl50-like 和 NylC-like 两套 grafted model 对比。

NylC-V3 / NylC-HP：
在 NylC tetramer 模型上建突变，复查 L4 register 和 L8/slab chain-entry 是否改变。
```

最关键对比仍然是：

```text
NylC tetramer：是否偏 A/D + Register A，释放 L1；
Nyl50 dimer：是否偏 two-subpocket + Register B，释放 L2；
Nyl10：更像 Nyl50，还是代表另一类 PA66-selective route。
```

---

## 6. 最小执行版本

第一版只做 NylC：

```text
1. NylC mature tetramer apo 检查。

2. NylC tetramer + L4-Nfree：
   A/D Register A/B/C；
   A/B Register A/B/C；
   2 个 random surface placement。

3. 每组：
   constrained placement → minimization → 1–5 ns weak restraint → 20–50 ns release MD。

4. 统计 productive-like frame fraction 和端基盐桥假象。

5. 选 2–3 个最有希望的 zone / register。

6. 对这些 zone 做 tethered L8。

7. 建 PA66 amorphous slab，做 8 个取向 × 2–3 replicate 的短 MD。

8. 若 A/D cleft 有链端接近迹象，再做 slab + exposed chain end。
```

---

## 7. 一句话总结

```text
L4 多方向 / 多 register 扫描负责回答：NylC 四聚体上哪个界面能形成真正有反应意义的构象。

L8 / PA66 slab 多取向模拟负责回答：NylC 四聚体面对真实材料表面时倾向用哪个区域吸附，以及表面链端是否能进入这些 functional cleft。

两者必须配合：L4 定义 productive-like pose 和 register，L8/slab 验证长链段和材料表面是否能到达这个 pose。
```
