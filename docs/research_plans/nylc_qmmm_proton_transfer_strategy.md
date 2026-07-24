# NylC PA66-L2 QM/MM 质子转移机制研究方案

## 核心目标
不是预设某一种催化机制，而是比较多个竞争机制，确定NylC催化PA66-L2时最可能的质子转移路径。

## 背景
NylC属于N-terminal nucleophile (Ntn) hydrolase。目前Ntn家族并不存在完全统一的质子转移模型。

已报道机制主要包括：
1. 协同机制（亲核攻击与质子转移同时发生）
2. N端氨基作为proton shuttle
3. Asp/Glu参与质子接受
4. 水桥(proton relay)

目前尚未发现针对NylC催化PA66底物的QM/MM机制研究，因此具有研究价值。

## 不建议的策略
不要直接计算一条预设反应路径。

## 推荐策略
建立同一Michaelis complex后，分别测试三条竞争路径：

### Path A
Thr267攻击羰基，随后质子转移至leaving amide nitrogen。

### Path B
Thr267质子首先转移至Asp306/Asp308，再完成后续步骤。

### Path C
催化水形成proton relay，完成质子转移。

比较三条路径的：
- TS结构
- 活化自由能
- IRC
- 四面体中间体稳定性

## 推荐计算流程
1. 优化Michaelis complex。
2. 100-500 ps QM/MM MD观察天然氢键网络。
3. 统计Thr-H与Asp、N端、水、leaving N的距离及占据率。
4. 根据天然构象分别建立Path A/B/C。
5. 每条路径进行relaxed scan。
6. 对最低能垒路径进行TS优化与IRC。
7. 最终进行自由能比较。

## QM区域建议
- Thr267完整侧链及N端氨基
- Asp306
- Asp308
- 被切割酰胺键附近两个repeat
- 2-3个催化水

## 可借鉴的QM/MM模板
- Penicillin G acylase
- Human Asparaginase III
- Proteasome
- γ-Glutamyltranspeptidase

重点不是照搬其机制，而是借鉴其比较竞争机制的研究框架。

## 最终目标
回答三个问题：
1. NylC真正的质子受体是谁？
2. PA66长链底物是否改变经典Ntn机制？
3. NylC、Nyl10、Nyl50是否采用不同的质子转移策略。