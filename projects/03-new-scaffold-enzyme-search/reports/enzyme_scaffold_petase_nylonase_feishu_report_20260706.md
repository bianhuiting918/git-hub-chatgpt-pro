# PETase 与尼龙酶骨架搜索筛选方法及初步统计

生成日期：2026-07-06  
更新日期：2026-07-06  
项目服务器路径：`/Dell/Dell14/bianht/enzyme_scaffold_search_v2`  
GPU 项目路径：`/data/bht/enzyme_scaffold_search_v2_gpu`  
适用范围：PETase 与尼龙酶（NylC / Nyl50 / Nyl10）骨架搜索、序列候选合并、缺结构序列预测、口袋分层统计、Pythia-pocket 评分、cluster 筛选、二聚体/四聚体界面口袋筛选前的交接报告。

## 1. 研究目标

本项目的目标不是只找到与已知酶序列相似的同源蛋白，而是尽可能系统地召回可能具有相似催化口袋或可迁移功能环境的新骨架。当前策略以“高召回、分层描述、生成缺失结构、后续再排序”为原则：

1. 用序列搜索扩大同源和远缘同源候选。
2. 用 Foldseek 从整体结构相似性角度召回已有结构或预测结构。
3. 用 Folddisco 从活性残基和周围口袋 motif 的结构几何角度召回候选。
4. 将 sequence / Foldseek / Folddisco 的候选合并去重，避免过早用单一方法排除潜在新骨架。
5. 对“有序列但没有结构”的候选建立 ESMFold 预测队列，生成结构后再回流到结构搜索和打分流程。
6. 对合并候选计算全局序列相似性、整体结构 RMSD、Pythia-pocket embedding cosine similarity，再进行 cluster、oligomer-aware pocket 和人工审查。
7. 对尼龙酶候选增加 oligomer-aware 判断：二聚体/四聚体倾向、界面口袋、跨链 substrate tunnel 与 PA6/PA66 productive docking。

关键原则：Folddisco、Foldseek、BLAST/MMseqs2 都是召回工具。工具返回 hit 不能直接等同于功能确认；下载失败、缺结构、缺序列、没有 predicted pocket、尚未排队预测结构，都不能写成生物学失败，应标记为 `NOT_EVALUATED` 或对应技术状态。

## 2. 证据宇宙与 denominator

当前项目至少有三类不同 denominator，必须分开报告：

| Universe | 含义 | 能做什么 | 不能直接做什么 |
|---|---|---|---|
| sequence recall candidates | BLAST/MMseqs/HMMER 从大库召回的序列 | 全局序列比对、cluster、判断是否需要结构预测 | 没有结构时不能做 pocket RMSD、Foldseek、Folddisco、Pythia-pocket |
| structure recall candidates | Foldseek / direct Folddisco 返回的结构候选 | 结构 RMSD、layer motif、Pythia-pocket、Foldseek/Folddisco 证据 | 不代表覆盖了所有序列搜索结果 |
| predicted-structure candidates | sequence recall 中有序列但缺结构，经过 ESMFold/ColabFold 生成结构的候选 | 回流到结构搜索、Pythia-pocket、RMSD、Folddisco 评分 | 预测前不能当作结构搜索已评估 |

因此，报告中不能写“所有序列都已经结构评估”。准确说法应是：

```text
所有候选先合并为 candidate universe；
已有结构的候选进入结构评估；
只有序列但缺结构的候选进入结构预测队列；
结构预测完成后再回流到结构评估。
```

## 3. 完整数据流

当前应按下面的数据流推进，而不是把序列搜索和结构搜索看成两个互不相干的结果表。

```text
Seed set
  ↓
Sequence recall: BLASTP / MMseqs2 / HMMER
  ↓
标准化 ID: UniProt / UniParc / MGYP / PDB / AFDB / sequence hash
  ↓
与 Foldseek / Folddisco structure candidates 合并去重
  ↓
分三类：
  A. 已有结构：PDB / AFDB / Foldseek / Folddisco 可用结构
  B. 有序列但没有结构：进入 ESMFold / ColabFold 预测队列
  C. 没有序列或 ID 无法解析：NOT_EVALUATED，进入 sequence recovery
  ↓
A 和 B 生成的结构统一进入：
  Foldseek / Folddisco after-recall scoring
  L0 / L0+L1 / L0+L2 / L0+L1+L2 layer 统计
  Pythia-pocket embedding cosine
  global RMSD
  ↓
所有有序列候选统一计算：
  Needleman-Wunsch global PID
  cluster90 / cluster50 / cluster30
  ↓
尼龙酶额外进入：
  oligomer-aware interface pocket / tunnel / docking
  ↓
candidate-level final table
```

## 4. 去重与回流规则

### 4.1 合并去重 key

不同来源的 ID 格式不一致，不能只按 accession 字符串去重。建议按以下优先级建立 candidate key：

1. 精确 accession：UniProt / UniParc / MGYP / PDB chain / AFDB accession。
2. sequence hash：同一氨基酸序列来自不同数据库时合并为同一 candidate。
3. structure ID mapping：Foldseek / Folddisco 返回结构 ID 能映射到 accession 时并入已有 candidate。
4. 近重复 cluster：不能替代 exact dedupe，但可用于后续 90% 去冗余。

### 4.2 去重后的状态标签

| 状态 | 含义 | 后续动作 |
|---|---|---|
| `STRUCTURE_AVAILABLE` | 已有 PDB/AFDB/预测结构可用 | 进入 Pythia-pocket、RMSD、Folddisco scoring |
| `SEQUENCE_ONLY_NEEDS_STRUCTURE` | 有序列，但没有可用结构 | 进入 ESMFold/ColabFold 队列 |
| `SEQUENCE_RECOVERY_NEEDED` | 有候选 ID，但序列未取回 | 继续 UniProt/UniParc/MGYP recovery |
| `PREDICTION_RUNNING` | 结构预测已排队或正在运行 | 等待 PDB 产出，不能写成 fail |
| `PREDICTION_DONE` | ESMFold/ColabFold PDB 已生成 | 回流结构评估 |
| `NOT_EVALUATED` | 缺结构、缺序列、下载/解析失败或尚未排队 | 技术未评估，不是生物学失败 |

### 4.3 结构生成后的回流

ESMFold 生成的 PDB 不是终点。每一批新 PDB 同步回 CPU 后，需要重新执行：

1. 更新 structure manifest。
2. 运行 Pythia-pocket predicted-pocket embedding。
3. 计算与所有对应 seed 的 pocket cosine。
4. 计算整体 backbone RMSD。
5. 对该结构运行 layer motif scoring：L0、L0+L1、L0+L2、L0+L1+L2。
6. 将结果并回 candidate-level 主表。
7. 对新增序列纳入 cluster90 / cluster50 / cluster30。

## 5. 总体筛选流程

| 阶段 | 工具/方法 | 目的 | 当前解释口径 |
|---|---|---|---|
| 序列大库召回 | BLASTP / MMseqs2 / HMMER | 从 UniRef90、UniProtKB/TrEMBL 等大库召回同源和远缘同源序列 | tool-native hit，只作为候选来源 |
| 序列标准化 | accession mapping + sequence hash | 把不同数据库命中的同一蛋白合并 | exact dedupe，不能替代 cluster |
| 结构召回 | Foldseek | 通过整体或局部结构相似性召回结构候选 | raw alignment 与 unique target 分开统计 |
| 口袋/基序召回 | Folddisco | 用催化残基和周围 residue layer 定义 motif，检索结构库 | Folddisco hit 是检索命中，不是最终功能 pass |
| 缺结构预测 | ESMFold / ColabFold | 给 sequence-only 候选生成结构 | 预测完成后再进入结构评估 |
| Layer 统计 | L0、L0+L1、L0+L2、L0+L1+L2 | 描述候选对核心催化残基及周围口袋的支持程度 | 先统计，不作为进入 Pythia 前的硬筛 |
| Pocket embedding | Pythia-pocket | 由模型预测 pocket residues，并计算 pocket embedding cosine | 当前使用 `prob_threshold = 0.6` |
| 全局序列相似性 | Needleman-Wunsch global alignment | 计算全长序列差异 | 使用 `PID_global = identical_residues / alignment_columns_including_gaps` |
| Cluster | MMseqs2 / CD-HIT | 去冗余与 scaffold diversity | 建议 cluster90 / cluster50 / cluster30 |
| 整体骨架差异 | Foldseek/TM-align 类 RMSD 统计 | 评估候选整体骨架是否与 seed 接近或偏离 | 用于散点图和排序解释 |
| 寡聚体界面 | ColabFold-Multimer / AlphaFold-Multimer / template assembly | 判断尼龙酶候选是否形成 NylC/Nyl50-like functional interface | 不问“任意二聚体”，而问“正确界面是否形成” |
| 界面口袋/通道 | P2Rank / fpocket / CAVER / POVME / docking | 判断二聚体/四聚体界面是否构成 PA6/PA66 可进入的 catalytic pocket/tunnel | 放在 monomer 初筛之后、实验候选之前 |

## 5A. 关键搜索与打分参数

本节记录当前已经在服务器脚本或日志中确认的参数。没有日志证据的参数不写成已执行标准，统一标为 `TO_VERIFY`。

### 5A.1 序列大库搜索参数

数据来源：`scripts/07_run_broad_sequence_recall.sh`、`scripts/05_run_hmmer_swissprot.sh`、`scripts/06_prepare_large_sequence_dbs.sh`。

| 项目 | 当前参数 |
|---|---|
| 大库 | UniRef90；UniProtKB/TrEMBL complete FASTA |
| 初始 HMM 构建库 | UniProt Swiss-Prot |
| BLASTP | `-evalue 1e-3 -max_target_seqs 4000 -num_threads $THREADS` |
| BLASTP output | `qseqid sseqid pident length qlen slen evalue bitscore qcovs stitle` |
| MMseqs2 | `mmseqs search -s 9.5 -e 1e-3 --max-seqs 4000 --threads $THREADS` |
| MMseqs2 output | `query,target,pident,alnlen,qlen,tlen,evalue,bits,qcov,tcov,theader` |
| HMM 构建 | seed FASTA 先用 `mafft --auto --quiet` 比对，再 `hmmbuild` |
| Swiss-Prot HMMER 初筛 | `hmmsearch --cpu 16 -E 1e-5 --domE 1e-5` |
| UniRef90/TrEMBL HMMER broad recall | `hmmsearch --cpu $THREADS --noali -E 1e-3 --domE 1e-3` |
| broad recall 合并 | 按 class + normalized target 合并，保留 method、database、best e-value、best bitscore、best qcov/tcov |

解释：`max_target_seqs / max-seqs = 4000` 是每个 seed/class/database/method 的召回上限，不是最终候选上限。最终 unique candidates 需要跨 seed、method、database 去重后统计。

### 5A.2 Foldseek 参数

数据来源：`scripts/13_prepare_foldseek_target_dbs.sh`、`scripts/14_run_foldseek_seed_search.sh`、`scripts/16_run_online_foldseek_seed_search.py`。

| 场景 | 当前参数 |
|---|---|
| 本地 Foldseek DB 下载/构建 | `foldseek databases <DB_NAME> <OUT_DB> <TMP_DIR> --threads ${THREADS:-32}` |
| 本地 Foldseek target DB | PDB；Alphafold/Swiss-Prot；Alphafold/UniProt50-minimal；以及按需其他 Foldseek databases |
| 本地 Foldseek search | `foldseek search <query_db> <target_db> <out_db> <tmp> --threads ${THREADS:-32} -e 1e-2 --max-seqs 10000` |
| 本地 Foldseek output | `query,target,evalue,bits,alntmscore,qtmscore,ttmscore,lddt,prob,qcov,tcov,qaln,taln,taxid,taxname` |
| 在线 Foldseek endpoint | `https://search.foldseek.com/api/ticket` |
| 在线 Foldseek mode | `3diaa` |
| 在线 Foldseek databases | `afdb50,mgnify_esm30,BFVD,gmgcl_id,pdb100,cath50` |
| 在线 Foldseek submit sleep | 8 s |
| 在线 Foldseek poll interval | 30 s |
| 在线 Foldseek timeout | 120 s |
| 在线 Foldseek rate-limit sleep | 300 s |

解释：Foldseek hit 是结构召回证据；后续表格需要区分 raw alignment、unique target、可下载结构、可评估结构四个层级。

### 5A.3 Folddisco direct search 参数

数据来源：`scripts/19_run_online_folddisco_pet_triads.py`、`results/layered_folddisco/runs/direct_search/online_missing_l0_l3/submissions/online_folddisco_tickets.tsv`。

| 项目 | 当前参数 |
|---|---|
| 在线 Folddisco endpoint | `https://search.foldseek.com/api/ticket/folddisco` |
| 提交字段 | query PDB file + `motif` + `database[]` |
| direct Folddisco databases | `afdb50_folddisco,afdb-proteome_folddisco,esm30_folddisco,BFVD_folddisco,pdb_folddisco` |
| submit sleep | 60 s |
| poll interval | 60 s |
| timeout | 120 s |
| result download | `/api/result/folddisco/download/{ticket}` |
| direct query motif | 使用 L0、L0+L1、L0+L2、L0+L1+L2 分层 residue 列表分别提交 |

解释：direct Folddisco 的 hit 是“在线数据库中按 motif 返回的结构命中”。这不是最终 pass/fail；后续还要合并去重、计算 layer support、Pythia-pocket 和 RMSD。

### 5A.4 Folddisco after-recall scoring 参数

数据来源：`results/layered_folddisco/runs/after_recall/query_result_manifest.tsv`、`scripts/21_run_foldseek_only_folddisco_batch.py`。

| 项目 | 当前参数 |
|---|---|
| index command | `folddisco index -p <candidate_structure_dir> -i <index_path> -r -t 16 --id relpath` |
| query command | `folddisco query -p <seed_pdb> -q <motif_residues> -i <index_path> -t 16 --per-structure` |
| output fields | `tid,total_match_count,node_count,idf,nres,plddt,max_node_cov,min_rmsd,matching_residues,query_residues` |
| sort rule | `--sort-by max_node_count:desc,idf:desc,min_rmsd:asc` |
| after-recall candidate_count example | PETase Foldseek-derived set: 4,059 structures per class-level index |
| foldseek-only batch download workers | `--download-workers 16` default |
| foldseek-only batch timeout | 60 s default |
| foldseek-only batch all mode | `--limit-per-class 0` means all rows |

解释：after-recall Folddisco 是“拿我们已有候选结构作为数据库，再用 Folddisco 计算 motif support”。它不受在线数据库 top 返回条数影响，适合对 sequence/Foldseek 得到的候选做统一结构 motif 打分。

### 5A.5 Layer motif 参数

当前项目实际执行时使用分层 motif，而不是单一 strict motif。代表性 PETase 参数示例：

| Seed | Layer | Motif residue list | Motif size |
|---|---|---|---:|
| PET_01_6ILW | L0 | `A160,A206,A237` | 3 |
| PET_01_6ILW | L1 | `A87,A161,A236` | 3 |
| PET_01_6ILW | L2 | `A159,A185,A238` | 3 |
| PET_01_6ILW | L0+L1 | `A160,A206,A237,A87,A161,A236` | 6 |
| PET_01_6ILW | L0+L2 | `A160,A206,A237,A159,A185,A238` | 6 |
| PET_01_6ILW | L0+L1+L2 | `A160,A206,A237,A87,A161,A236,A159,A185,A238` | 9 |

所有 seed 都按同一思想生成 L0、L1、L2、L0+L1、L0+L2、L0+L1+L2；具体 residue list 以 `results/layered_folddisco/runs/after_recall/query_result_manifest.tsv` 和 layered definitions 文件为准。

### 5A.6 ESMFold 结构预测参数

数据来源：`scripts/predict_failed_targets_esmfold.py` 和 GPU 运行命令。

| 项目 | 当前参数 |
|---|---|
| 模型 | `esm.pretrained.esmfold_v1()` |
| 设备 | `.eval().cuda()`，GPU A100 80GB 当前在跑 |
| 输入 | FASTA records |
| 输出 | `<outdir>/pdb/<candidate_id>.pdb` + `esmfold_prediction_manifest.tsv` |
| 当前尼龙酶第一批运行参数 | `--sort-by-length --chunk-size 64` |
| 脚本默认 chunk-size | 128 |
| 当前 batch size | 脚本逐条 `model.infer_pdb(seq)`，不是多序列 batch |
| 可选限制 | `--max-count`、`--max-len`，默认 0 表示不限制 |

解释：ESMFold 只负责把 sequence-only candidate 变成 predicted structure；结构质量、pLDDT、pocket 是否可用，要在后续 Pythia/RMSD/Folddisco 阶段再判断。

### 5A.7 Pythia-pocket 参数

数据来源：`scripts/run_predicted_pocket_pythia_cosine.py`、`scripts/run_structure_merged_layer_pythia_pocket_cosine.py`。

| 场景 | 当前参数 |
|---|---|
| predicted-pocket 模式 | Pythia-pocket 预测 residue probability，按阈值选 pocket residues |
| predicted-pocket threshold | `--threshold 0.6` 默认 |
| predicted-pocket model | `tools/src/Pythia/pythia-pocket/model-fl-b.pt` |
| cosine | 对 seed pocket embedding 与 candidate pocket embedding 做 cosine similarity |
| no pocket label | seed 或 candidate 在阈值下无 residue 时标记 `NO_PREDICTED_POCKET` |
| resume | `--resume` 可跳过已完成 candidate/seed pairs |
| flush | `--flush-every 25` 默认 |
| CPU threads | `PYTHIA_TORCH_THREADS`，脚本默认 64 或 32；实际高核脚本曾设置 96，但当前 CPU 使用应受 64 核限制 |
| layer-specific 模式 | `PYTHIA_POCKET_LAYER=L0_L1_L2` 时使用 mapped seed/target residues pooled embedding |

解释：predicted-pocket Pythia 与 L0-L3 手工 layer 是两种不同口袋定义。前者由模型预测 pocket residue；后者由我们指定 catalytic/environment residues。报告和主表必须分别保留这两类字段。

### 5A.8 全局序列相似性参数

当前用户要求的最终序列相似性指标是全局比对，而不是 BLAST 局部相似性：

```text
PID_global = identical_residues / alignment_columns_including_gaps
global_sequence_difference = 1 - PID_global
query_coverage = aligned_query_residues / query_length
target_coverage = aligned_target_residues / target_length
gap_fraction = gap_columns / alignment_columns
```

其中用于散点图颜色的核心值是 `PID_global`，避免局部高相似片段导致“序列相似度很高但 coverage 很低”的误读。

### 5A.9 Cluster 参数

当前报告建议但尚需统一执行的 cluster 参数：

| Cluster level | 参数 | 目的 |
|---|---|---|
| cluster90 | sequence identity 90% | 工程去冗余，删除近重复 |
| cluster50 | sequence identity 50% | subfamily 代表选择，对应 NylC 文献 diversity panel 的 30%-50% 区间 |
| cluster30 | sequence identity 30% | family/scaffold SSN，识别远缘 scaffold |

建议执行时使用 MMseqs2 cluster/linclust 或 CD-HIT，并在结果表中写入具体命令，例如 `--min-seq-id`、coverage mode、coverage cutoff。当前文档中 cluster90/50/30 是项目决策参数，不应误写成已经全部完成的过滤结果。

## 6. Layer 与口袋定义

当前 layer 是项目定义的分层口袋描述，不是 Folddisco 自带的唯一标准。Folddisco 自身负责按给定 motif 和数据库参数返回结构命中；我们另外统计不同 layer 下候选对 motif 的支持情况。

| Layer | 定义 | 用途 |
|---|---|---|
| L0 | 核心催化残基或最关键的反应相关残基 | 作为 Pythia-pocket 对比和结构 motif 的核心 anchor |
| L1 | L0 周围更近的口袋支撑残基，通常直接参与底物定位、催化几何或局部环境 | 用来描述核心催化环境是否完整 |
| L2 | 更外层的口袋环境或二级支撑残基，可能影响底物通道、形状和稳定性 | 用来描述口袋上下文是否接近 |
| L0+L1 | 核心催化残基加近邻口袋 | 比 L0 更严格，但仍避免过度收缩 |
| L0+L2 | 核心催化残基加外层环境 | 用来捕获可能保留核心和外层形状但近邻发生替换的候选 |
| L0+L1+L2 | 当前最完整的项目口袋层 | 用于后续优先级解释，不作为唯一进入 Pythia 的门槛 |

这里的“通过 layer”应理解为“该候选在对应 layer 下有结构映射/几何支持”，而不是功能已确认。后续 Pythia-pocket cosine、整体 RMSD、全局序列差异、cluster 与 oligomer-aware pocket geometry 会一起用于排序。

## 7. PETase：seed、搜索与缺结构序列状态

### 7.1 PETase seed 来源

当前 PETase 有 5 个结构 seed。早期序列大库搜索使用了 4 个有明确 UniProt accession 的 seed；后续已补跑 PET_05_3VIS 的 PDB-derived sequence。五个 seed 不是互相替代关系，而是用于覆盖不同 PET hydrolase / cutinase scaffold。

| Seed | Accession / 来源 | 代表结构 | 说明 | 状态 |
|---|---|---|---|---|
| PET_01 | A0A0K8P6T7 | 6ILW | IsPETase，典型 mesophilic PETase | accepted |
| PET_02 | G9BY57 | 4EB0 | LCC，高温 PET hydrolase/cutinase | accepted |
| PET_03 | Q6A0I4 | 4CG1 | TfCut2 / Thermobifida cutinase scaffold | accepted |
| PET_04 | W0TJ64 | 4WFI / 4WFJ | Cut190，Ca2+ regulated PET-active scaffold | accepted |
| PET_05 | PDB-derived sequence | 3VIS | Est119/TaCut-like cutinase，用于补充 thermophilic actinobacterial cutinase 多样性 | accession 待最终确认 |

### 7.2 PETase 序列搜索统计

| 搜索集合 | 方法/数据库 | 初步统计 |
|---|---|---:|
| 早期 PETase broad recall，4 个 accession seed | BLASTP / MMseqs2 / HMMER against UniRef90 + TrEMBL | 7,306 unique candidates |
| PET_05_3VIS 补充搜索 | BLASTP + MMseqs2 against UniRef90 + TrEMBL | 6,358 unique targets |
| 已有候选 FASTA 中可用于 seed 相似性计算的 PETase 序列 | candidate FASTA vs seed FASTA | 25,367 candidate sequences；25,258 有 seed hit |

PET_05_3VIS 补充搜索的细分结果：

| 方法 | 数据库 | Raw hits | Unique targets |
|---|---|---:|---:|
| BLASTP | UniRef90 | 3,433 | 1,730 |
| BLASTP | TrEMBL | 4,043 | 2,033 |
| MMseqs2 | UniRef90 | 1,936 | 1,936 |
| MMseqs2 | TrEMBL | 2,218 | 2,218 |
| 合并去重 | UniRef90 + TrEMBL | - | 6,358 |

### 7.3 PETase 结构候选统计

当前 PETase 结构候选来自 Foldseek 与 direct Folddisco 的合并结果。RMSD 报告中每个 PETase seed 对应 10,089 个候选行。

| 来源 | 行数 |
|---|---:|
| Foldseek-only | 3,959 |
| Folddisco-only | 6,030 |
| Foldseek and Folddisco overlap | 100 |
| 每个 PETase seed 可评估结构候选行 | 10,089 |

这里的行数是结构/RMSD 可评估集合，不等同于序列搜索全集。

### 7.4 PETase sequence-only no-structure 当前状态

截至 2026-07-06 22:34 +08:00，PETase 的 sequence-only no-structure 预测队列还没有真正建立完成：

| 检查项 | 状态 |
|---|---|
| CPU `results/petase_sequence_available_no_structure_prediction/` | 目录存在 |
| CPU `esmfold_top5000/pdb` | 0 个 PDB |
| CPU `top5000_sequence_available_no_structure.fasta` | 不存在 |
| CPU `petase_sequence_available_no_structure.fasta` | 不存在 |
| GPU PETase ESMFold 输入 FASTA | 不存在 |
| GPU PETase ESMFold 进程 | 未运行 |

因此 PETase 需要补做一个明确步骤：

```text
PETase sequence recall candidates
  ↓
与 PETase structure merged candidates exact dedupe
  ↓
筛出 sequence_available && !structure_available
  ↓
生成 PETase sequence-only no-structure FASTA
  ↓
上传 GPU，排队 ESMFold
  ↓
PDB 回流 Pythia-pocket / RMSD / Folddisco scoring
```

## 8. 尼龙酶：seed、搜索与缺结构序列状态

### 8.1 尼龙酶 seed 来源

尼龙酶当前以 NylC、Nyl50、Nyl10 三类为核心。NylC 的 WT / GYAQ / HP 多数属于同一大骨架背景；Nyl50 和 Nyl10 是更需要关注的新 scaffold / 功能方向。

| Seed | Scaffold group | 用途 | 说明 |
|---|---|---|---|
| NylC_WT_Q79F77 / NylC_WT_5XYO | NylC-like | sequence + structure | PA6-biased baseline NylC |
| NylC_GYAQ_5Y0M | NylC-like | sequence + structure | thermostable GYAQ background |
| NylC_HP_from_GYAQ | NylC-like | sequence + backbone model | HP mutation background，主要用于 backbone retrieval |
| Nyl50_9CXR_9DYS | Nyl50-like | sequence + structure | PA66-selective new scaffold candidate |
| Nyl10_A0A1M5P6R3 | Nyl10-like | sequence + structure/model | PA66-selective purified enzyme；AFDB v6 model 用于结构搜索 |

### 8.2 尼龙酶序列搜索统计

| 搜索集合 | 方法/数据库 | 初步统计 |
|---|---|---:|
| 早期 nylonase broad recall | BLASTP / MMseqs2 / HMMER against UniRef90 + TrEMBL | 33,675 unique candidates |
| PA6/PA66 扩展 seed：NylC WT/GYAQ/HP、Nyl50 | BLASTP + MMseqs2 against UniRef90 + TrEMBL | 已完成 |
| Nyl10 追加搜索 | BLASTP + MMseqs2 against UniRef90 + TrEMBL | 已完成 |
| 已有候选 FASTA 中可用于 seed 相似性计算的尼龙酶序列 | candidate FASTA vs seed FASTA | 15,083 candidate sequences；14,803 有 seed hit |

Nyl10 不是完全搜不到。它在大库序列搜索中有 self-hit，Nyl50 与 NylC 搜索也能命中 Nyl10，但相似度层级不同：Nyl50 到 Nyl10 的 BLAST pident 约 40.6%，query coverage 约 93%；NylC 到 Nyl10 的 BLAST pident 约 30-31%，query coverage 约 76%。这说明 Nyl10 与 NylC-like seed 的序列距离较远，不能只用序列相似性判断是否同类。

### 8.3 尼龙酶结构合并集合与 layer 统计

当前重新定义 layer 后的 nylonase merged candidate universe：

| 指标 | 数量 |
|---|---:|
| structure_merged_candidates_nylonase | 21,211 |
| with existing available structure | 4,035 |
| redefined Folddisco evidence rows | 43,306 |
| redefined Folddisco unique accessions | 24,010 |
| predicted-structure manifest 当前行数 | 3,192 rows / 3,166 unique candidates |
| predicted all-seed embedding 当前行数 | 6,200 rows / 3,106 unique candidates |

Layer 支持统计：

| Layer | Any support | Full support |
|---|---:|---:|
| L0 | 15,284 | 5,383 |
| L0+L1 | 6,701 | 159 |
| L0+L2 | 7,320 | 265 |
| L0+L1+L2 | 5,301 | 51 |

候选来源 route 统计：

| Route | Candidates |
|---|---:|
| direct_folddisco | 17,134 |
| foldseek | 3,555 |
| foldseek;direct_folddisco | 163 |
| sequence;direct_folddisco | 63 |
| sequence;foldseek | 275 |
| sequence;foldseek;direct_folddisco | 21 |

### 8.4 尼龙酶 sequence-only no-structure 当前状态

截至 2026-07-06 22:34 +08:00，尼龙酶“序列有了但没有结构”的候选已经进入 GPU ESMFold 队列，且第一批正在运行。

| 队列 | 输入序列数 | GPU 已生成 PDB | CPU 已同步 PDB | 状态 |
|---|---:|---:|---:|---|
| `failed_structure_targets_with_sequences.fasta` | 5,606 | 4,296 | 3,876 | GPU 正在跑 |
| `remote_fetched_missing_sequences.fasta` | 3,633 | 0 | 0 | waiter 等第一批结束后自动开始 |
| MGYP partial / incremental ready | 若干批 | 0 | 0 | waiter 等前序 ESMFold 结束 |

当前 GPU 进程：

```text
predict_failed_targets_esmfold.py
--fasta failed_structure_targets_with_sequences.fasta
--outdir structure_prediction/esmfold_failed_structure_targets
--sort-by-length
--chunk-size 64
```

解释：尼龙酶这边不是没排队，而是 GPU 正在跑第一批；后续 remote-fetched / MGYP 补回序列已经有 waiter，但要等前序任务结束。

## 9. Pythia-pocket embedding 与相似性

当前 Pythia-pocket 的使用方式是：不直接用我们手工定义的 L0-L3 residue 作为最终 pocket，而是让 Pythia 按模型预测 pocket residues，然后提取 predicted-pocket embedding。我们使用 seed 的核心口袋作为比较对象，计算候选与 seed pocket embedding 的 cosine similarity。

当前参数与解释：

| 项目 | 当前设置 |
|---|---|
| Pocket residue 来源 | Pythia-pocket predicted pocket |
| Probability threshold | 0.6 |
| 相似性指标 | cosine similarity |
| 没有 predicted pocket | 标记为 `NO_PREDICTED_POCKET`，不等同功能失败 |
| 缺结构或解析失败 | 标记为 `NOT_EVALUATED` |

尼龙酶 all-seed predicted-pocket embedding 当前状态：

| 状态 | 数量 |
|---|---:|
| OK | 5,570 |
| NO_PREDICTED_POCKET | 630 |

按 seed 的 all-seed predicted-pocket rows：

| Seed | Rows |
|---|---:|
| Nyl10_A0A1M5P6R3 | 3,106 |
| Nyl50_9DYS | 2,790 |
| NylC_WT_5XYO | 304 |

已读取到的 NylC 与 Nyl10/Nyl50 的 Pythia-pocket cosine：

| 比较 | Cosine | 口径 |
|---|---:|---|
| NylC_WT_5XYO vs Nyl10_A0A1M5P6R3 | 0.860980 | Pythia predicted pocket, prob_threshold 0.6 |
| NylC_WT_5XYO vs Nyl50_9DYS | 0.925020 | Pythia predicted pocket, prob_threshold 0.6 |
| NylC_WT_5XYO vs NylC_WT_5XYO | 1.000000 | self |

### 9.1 Nylonase-specific pocket cosine 分层

GRASE 在 Aes72-centered urethanase discovery 中使用的是 key catalytic residue embedding cosine 1.00-0.97 四个 tier。但 nylonase 已知有效 seed 之间本身更分散：NylC vs Nyl50 为 0.925020，NylC vs Nyl10 为 0.860980。因此不能把 0.97-1.00 直接作为 nylonase 硬阈值，否则会把 Nyl50-like 和 Nyl10-like 已知有效方向排除。

尼龙酶建议使用 multi-seed max cosine，而不是只用 NylC seed：

```text
max_pocket_cosine = max(
  cosine_to_NylC,
  cosine_to_Nyl50,
  cosine_to_Nyl10
)
```

建议分层：

| Tier | 阈值 | 解释 | 处理建议 |
|---|---:|---|---|
| Nyl-pCOS-high | `max_pocket_cosine >= 0.93` | NylC-like / Nyl50-like 近口袋候选 | 高优先级，但仍需 cluster 去冗余和二聚体界面检查 |
| Nyl-pCOS-mid | `0.86 <= max_pocket_cosine < 0.93` | Nyl10-like / Nyl50-like 远缘有效候选区间 | 必须保留；不能按 GRASE 阈值丢弃 |
| Nyl-pCOS-low-but-salvageable | `0.80 <= max_pocket_cosine < 0.86` | 探索区间 | 需要 strong motif、Folddisco、二聚体界面或 CAVER tunnel 支持 |
| Nyl-pCOS-exploratory | `< 0.80` | 低相似探索区间 | 一般不进首轮实验，除非 cluster novelty 或几何证据很强 |

## 10. 文献筛选策略对本项目的约束

### 10.1 GRASE 查找聚氨酯酶：高 embedding 阈值适合 Aes72-centered urethanase，不宜直接套到 nylonase

GRASE 的 template 是 Aes72。作者先实验测试 14 个文献或专利来源的 PURase/urethane hydrolase，只有 Aes72、GatA、SP1 对 N-aryl carbamates 有明显活性，其中 Aes72 对 2,4-TDA-DEG 的活性最高，因此被选为后续搜索模板。

GRASE 前先做了两个 baseline：传统序列/结构搜索，以及 cluster + complex prediction。最终有效路线是 GRASE：Pythia-Pocket 对 Aes72 和候选结构的 key catalytic residues 生成 residue-level embeddings，计算 cosine similarity；Pythia 同时用 NLL 估计稳定性。作者把候选与 Aes72 的 Pythia-Pocket cosine similarity 分成 1.00-0.97 的四个 tier，最终得到 24 个 GRASE candidates，其中 21 个表达且有 urethanase activity。

对本项目的关键判断：**GRASE 的 0.97-1.00 threshold 不能直接套到 nylonase。** 当前 nylonase 统计中 NylC-WT vs Nyl10 的 predicted-pocket cosine 为 0.860980，NylC-WT vs Nyl50 为 0.925020。如果使用 `cosine >= 0.97` 作为硬阈值，Nyl10 和 Nyl50 这两个已知有效方向都会被排除。nylonase 应使用 multi-seed tier，而不是单 seed 高阈值过滤。

### 10.2 Metal-coordination mining：机制几何 motif 先定义功能类型，cluster 和底物假说决定实验序列

该 Nature 工作研究的是 Fe(II)/alphaKG radical halogenase，核心不是 embedding，而是机制驱动的金属配位几何筛选。其可迁移逻辑是：

```text
motif / pocket / geometry
负责判断“可能有目标反应能力”

30% cluster
负责判断“属于哪个 family / scaffold group”

50% cluster
负责判断“哪个 subfamily 值得选代表”

annotation / closest characterized homolog / genome context
负责判断“实验应该测 PA6、PA66 还是更具体的 oligomer/product panel”
```

### 10.3 NylC / Nyl10 / Nyl50 文献：不是 GRASE 式高 embedding 筛选，而是 diversity panel + 实验读出

Nyl10 / Nyl50 文献的路线和 GRASE 完全不同。作者以 NylC-GYAQ 氨基酸序列作为 query，在 ColabFold 的 MMseqs2-based sequence search 中搜索 UniRef100，得到 2,839 条 homologous sequences；随后用 HHfilter 要求与 NylC-GYAQ 的 coverage >= 50%，留下 2,643 条。

为了构建 96-well plate，作者从 2,643 条中选择 95 条 maximally diverse homologs，核心标准是 lowest maximum pairwise sequence identity，而不是最像 NylC。最终 panel 的两两 identity 在 30%-50%。随后全部合成表达，用 PA6 / PA66 粗筛直接读出活性和产物谱，最终识别 Nyl10、Nyl12、Nyl50 等代表。

对本项目的含义：当前 Nyl10 与 NylC 的 BLAST pident 只有约 30-31%，Nyl50 与 Nyl10 约 40.6%，而 NylC vs Nyl10 pocket cosine 也只有 0.861。因此 Nyl10-like scaffold 很容易被高 identity 或高 pocket cosine 筛法排除。必须保留 30%-50% identity 的远缘 cluster 代表。

## 11. Nylonase 需要引入二聚体/四聚体-aware 筛选

NylC/Nyl50-like nylonase 与 PETase 不同。PETase 通常可按 monomer active-site pocket 处理；NylC/Nyl50-like 系统可能依赖 oligomeric interface 形成完整底物通道或活性口袋。

Nyl50 结构分析显示，Nyl50 的 putative substrate access tunnel 是跨两个 monomers 的 U-shaped tunnel；入口半径约 2.5 Å，active-site 附近 bottleneck radius 约 1.5 Å，tunnel volume 约 631 Å³，并且一个 large subpocket 位于 monomer interface，含有来自两个亚基的残基。文献还指出 NylC-GYAQ 和 Nyl50 的 putative active sites 跨越 A/D dimer interface，提示 oligomerization 对活性很可能重要。

因此，本项目后续不应只做 monomer Pythia-pocket embedding、monomer CAVER、monomer docking。这些可以作为高召回初筛，但不能作为最终实验排序依据。nylonase 最终应新增：

```text
oligomer-aware pocket / tunnel / interface screening
```

核心判断不是“是否形成任意二聚体”，而是：

```text
是否形成 NylC/Nyl50-like 正确二聚体/四聚体界面；
该界面是否参与 catalytic pocket / substrate access tunnel。
```

## 12. 下一步筛选决策

### 12.1 P0：补齐 sequence recall 到 structure prediction 的缺口

| 酶类 | 当前状态 | P0 动作 |
|---|---|---|
| PETase | sequence-only no-structure FASTA 未生成，GPU 未排队 | 生成 PETase sequence-only no-structure FASTA，和 structure merged candidates 去重后上传 GPU 排队 ESMFold |
| 尼龙酶 | 第一批 5,606 正在 GPU ESMFold；第二批 3,633 已 waiter | 等第一批完成并同步 CPU；第一批 PDB 回流 embedding/RMSD/Folddisco；第二批自动开始后继续同步 |

PETase 不能只停留在已有结构候选。需要把序列搜索得到但缺结构的候选也纳入 ESMFold，否则 sequence recall 的新增多样性不会进入 Pythia-pocket 或 RMSD 图。

### 12.2 必须进行 cluster，建议使用 90%、50%、30% 三层 identity

当前 nylonase 候选规模已经较大：broad recall 33,675 unique candidates，candidate FASTA 15,083 sequences，structure merged candidates 21,211，predicted all-seed embedding 3,106 unique candidates。这个规模下如果不 cluster，近重复序列会占据 top rank，实验候选会缺乏 family/scaffold 多样性。

建议三层 cluster：

| Identity | 用途 | 解释 |
|---:|---|---|
| 90% | 工程去冗余 | 删除近重复、strain-level redundancy、几乎相同序列 |
| 50% | subfamily 代表选择 | 对应 NylC 文献 diversity panel 的 30%-50% identity 区间，适合挑实验候选 |
| 30% | family / scaffold 级别 SSN | 对应 metal-coordination mining 的 family-level cluster 逻辑，用于识别远缘 scaffold 和 novel cluster |

可以加 70% identity 作为过渡层，但第一版建议先做：

```text
cluster90 + cluster50 + cluster30
```

### 12.3 不建议直接使用 GRASE 的 0.97-1.00 embedding 阈值

GRASE 的 0.97-1.00 是 Aes72-centered urethanase discovery 的经验阈值。nylonase 已知有效 seed 更分散，因此应使用：

```text
max cosine to NylC / Nyl50 / Nyl10
```

而不是只看 NylC seed。

建议 candidate-level 主表新增：

```text
cosine_to_NylC
cosine_to_Nyl50
cosine_to_Nyl10
max_pocket_cosine
best_pocket_seed
pocket_cosine_tier
```

### 12.4 引入二聚体/四聚体预测，但不要全量跑 AlphaFold-Multimer

建议漏斗：

```text
Step 1: sequence QC
去 fragment、去 X、去长度异常、保留 Ntn-hydrolase / Nyl-like motif

Step 2: MMseqs2 cluster
90%, 50%, 30%

Step 3: 每个 50% cluster 选 1-3 条代表
优先 Pythia-pocket OK
优先 layer support 强
优先 structure pLDDT 高
优先多 route 支持 sequence + Foldseek + Folddisco

Step 4: 对 top 500-1000 representative candidates 做 homodimer prediction
必要时做 homotetramer prediction

Step 5: 计算 oligomer confidence
ipTM、interface PAE、pDockQ / pDockQ2、buried SASA、interface contact number、multi-seed consistency、A/D or A/B contact-map similarity to Nyl50/NylC

Step 6: 只对 oligomer confidence 高，或 confidence 中等但 motif/tunnel 很强的候选，做 interface-pocket analysis
```

### 12.5 口袋几何匹配应从 monomer pocket 升级到 interface-pocket geometry

建议新增一个 `interface_pocket_geometry_score`，并与 Pythia-pocket cosine 并列，而不是替代它。建议指标：

```text
1. catalytic Thr / Ntn-hydrolase autocleavage motif 是否完整
2. catalytic residue RMSD to NylC/Nyl50/Nyl10
3. interface pocket 是否包含来自两个 chain 的 residues
4. candidate 是否有跨链 substrate tunnel
5. CAVER bottleneck radius
6. CAVER tunnel volume
7. tunnel-lining residues composition
8. PA6 / PA66 oligomer fragment docking 是否形成 productive pose
9. PA66-like / PA6-like substrate channel 是否有差异
```

### 12.6 首轮实验候选应按 group 分层抽样，而不是按单一总分 top N

如果首轮实验规模为 24-48 条，建议：

| Group | 选择逻辑 | 建议数量 |
|---|---|---:|
| A. NylC-like 高把握候选 | max cosine high，motif 完整，二聚体/四聚体界面合理 | 6-10 |
| B. Nyl50-like / PA66-selective 候选 | cosine_to_Nyl50 高，interface tunnel 接近 Nyl50，PA66 docking pose 合理 | 6-10 |
| C. Nyl10-like 远缘候选 | cosine 不一定高，但 sequence / structure / motif 指向 Nyl10-like | 6-10 |
| D. low-cosine but strong-geometry candidates | cosine 0.80-0.86，但 Folddisco motif、CAVER tunnel、dimer interface 强 | 4-8 |
| E. cluster novelty representatives | 30% cluster 中远离 NylC/Nyl50/Nyl10，但 motif 和结构可信 | 4-8 |

## 13. Candidate-level 主表建议字段

后续统一主表不应只包括 route、RMSD 和 pocket cosine。建议新增 sequence/structure 队列状态、cluster、multi-seed pocket、oligomer 与 interface pocket 字段：

```text
candidate_id
enzyme_class
route
sequence_available
structure_available
sequence_source
structure_source
structure_prediction_status
structure_prediction_batch
structure_prediction_pdb_path
cluster90_id
cluster50_id
cluster30_id
best_seed_by_sequence
best_seed_by_structure
best_seed_by_pocket_cosine
global_PID
global_RMSD
cosine_to_NylC
cosine_to_Nyl50
cosine_to_Nyl10
max_pocket_cosine
best_pocket_seed
pocket_cosine_tier
L0_support
L0L1_support
L0L2_support
L0L1L2_support
oligomer_model_status
predicted_oligomer_state
dimer_ipTM
dimer_interface_PAE
dimer_pDockQ2
interface_BSA
interface_contact_count
interface_contactmap_similarity_to_Nyl50
CAVER_tunnel_found
CAVER_bottleneck_radius
CAVER_tunnel_volume
CAVER_cross_chain
interface_pocket_residue_count
dock_PA6_productive_pose
dock_PA66_productive_pose
dock_best_substrate
dock_reactive_distance
dock_attack_angle
final_priority_group
status_label
```

## 14. 推荐状态标签

| 标签 | 含义 |
|---|---|
| TOOL_HIT | 工具原生检索命中，例如 Foldseek/Folddisco/BLAST/MMseqs 返回 hit |
| SEQUENCE_AVAILABLE | 已有可用于全局序列比对的序列 |
| STRUCTURE_AVAILABLE | 已有可用于 RMSD/Pythia 的结构 |
| SEQUENCE_ONLY_NEEDS_STRUCTURE | 序列可用但无结构，需要 ESMFold/ColabFold |
| PREDICTION_RUNNING | 结构预测已排队或正在运行 |
| PREDICTION_DONE | 预测结构已生成，等待或已回流结构评估 |
| PREDICTED_POCKET_OK | Pythia-pocket 在阈值下得到 pocket embedding |
| NO_PREDICTED_POCKET | Pythia-pocket 阈值下无 pocket residue，不等于生物学失败 |
| NOT_EVALUATED | 缺结构、缺序列、下载/解析失败或尚未进入队列 |
| CLUSTER_REPRESENTATIVE | 该候选是 90%/50%/30% cluster 的代表序列之一 |
| OLIGOMER_MODEL_OK | 二聚体/四聚体模型置信度可接受 |
| INTERFACE_POCKET_OK | oligomer interface 上存在与 catalytic site/tunnel 相关的 pocket |
| CROSS_CHAIN_TUNNEL_OK | CAVER/MOLE 显示存在跨链 substrate access tunnel |
| PRODUCTIVE_DOCKING_POSE | PA6/PA66 fragment docking 形成可反应构象 |
| PROJECT_PRIORITY | 后续根据多指标综合排序得到的候选优先级，不是工具原生标签 |

## 15. 当前结论

1. **必须把 sequence recall 的新增候选回流到结构流程。** 已有结构候选只能覆盖结构数据库中的部分 candidate；序列搜索得到但缺结构的候选需要先 ESMFold/ColabFold，再做 Pythia-pocket、RMSD、Folddisco layer scoring。
2. **尼龙酶 sequence-only no-structure 已经排队并正在跑。** 第一批 5,606 条中 GPU 已生成约 4,296 个 PDB；第二批 3,633 条已排 waiter。
3. **PETase sequence-only no-structure 还没有排队成功。** 需要立刻生成 PETase sequence-only no-structure FASTA，与已有 structure merged candidates 去重后上传 GPU。
4. **GRASE 的 0.97-1.00 pocket cosine 阈值不能直接用于 nylonase。** NylC vs Nyl50 为 0.925020，NylC vs Nyl10 为 0.860980；直接套用会漏掉已知有效 scaffold。
5. **Nylonase 应使用 multi-seed pocket cosine。** 以 NylC、Nyl50、Nyl10 三类 seed 的 max cosine 分层，保留 0.86-0.93 的远缘有效候选区间。
6. **必须 cluster。** 建议同时做 90%、50%、30% identity：90% 去冗余，50% 选 subfamily 代表，30% 做 family/scaffold SSN。
7. **monomer-only pocket 只能作为初筛。** NylC/Nyl50-like enzyme 很可能依赖 A/D dimer interface 或四聚体中的跨链 tunnel，最终排序必须引入 oligomer-aware pocket / tunnel / docking。
8. **实验候选不应按单一总分 top N。** 应按 NylC-like、Nyl50-like、Nyl10-like、low-cosine strong-geometry、cluster novelty 分组抽样。
