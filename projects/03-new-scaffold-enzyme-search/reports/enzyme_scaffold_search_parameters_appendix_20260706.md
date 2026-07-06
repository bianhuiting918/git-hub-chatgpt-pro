# Enzyme scaffold search 参数附录（PETase / 尼龙酶）

日期：2026-07-06  
主报告：`projects/03-new-scaffold-enzyme-search/reports/enzyme_scaffold_petase_nylonase_feishu_report_20260706.md`  
服务器路径：`/Dell/Dell14/bianht/enzyme_scaffold_search_v2`  
GPU 路径：`/data/bht/enzyme_scaffold_search_v2_gpu`

本附录只记录“搜索和打分参数”。研究思路、结果统计、下一步任务以主报告为准。

## 1. 总体口径

序列搜索、Foldseek、Folddisco 都作为高召回候选来源，不直接作为功能确认。后续统一合并去重后，再按 sequence / structure / predicted-structure 三类状态进入不同流程。

| 状态 | 含义 | 是否能做结构打分 |
|---|---|---|
| `SEQUENCE_AVAILABLE` | 已有氨基酸序列 | 不能，除非已有或预测出结构 |
| `STRUCTURE_AVAILABLE` | 已有 PDB / AFDB / Foldseek / Folddisco 可用结构 | 可以 |
| `SEQUENCE_ONLY_NEEDS_STRUCTURE` | 有序列但缺结构 | 先排队 ESMFold / ColabFold |
| `NOT_EVALUATED` | 缺序列、缺结构、下载失败、解析失败或尚未排队 | 不能写成生物学失败 |

## 2. 序列搜索参数

来源脚本：`scripts/07_run_broad_sequence_recall.sh`、`scripts/05_run_hmmer_swissprot.sh`、`scripts/06_prepare_large_sequence_dbs.sh`。

| 项目 | 参数 |
|---|---|
| 搜索库 | UniRef90；UniProtKB/TrEMBL complete FASTA |
| 初始 HMM 构建库 | UniProt Swiss-Prot |
| BLASTP | `-evalue 1e-3 -max_target_seqs 4000 -num_threads $THREADS` |
| BLASTP 输出字段 | `qseqid sseqid pident length qlen slen evalue bitscore qcovs stitle` |
| MMseqs2 | `mmseqs search -s 9.5 -e 1e-3 --max-seqs 4000 --threads $THREADS` |
| MMseqs2 输出字段 | `query,target,pident,alnlen,qlen,tlen,evalue,bits,qcov,tcov,theader` |
| HMM 构建 | seed FASTA 先 `mafft --auto --quiet`，再 `hmmbuild` |
| Swiss-Prot HMMER 初筛 | `hmmsearch --cpu 16 -E 1e-5 --domE 1e-5` |
| UniRef90 / TrEMBL HMMER broad recall | `hmmsearch --cpu $THREADS --noali -E 1e-3 --domE 1e-3` |
| 合并规则 | class + normalized target 合并，保留 method、database、best e-value、best bitscore、best qcov/tcov |

注意：`max_target_seqs` / `--max-seqs = 4000` 是每个 seed / database / method 的召回上限，不是最终候选数。最终候选要跨 seed、方法和数据库去重。

## 3. Foldseek 结构搜索参数

来源脚本：`scripts/13_prepare_foldseek_target_dbs.sh`、`scripts/14_run_foldseek_seed_search.sh`、`scripts/16_run_online_foldseek_seed_search.py`。

| 场景 | 参数 |
|---|---|
| 本地数据库构建 | `foldseek databases <DB_NAME> <OUT_DB> <TMP_DIR> --threads ${THREADS:-32}` |
| 本地 target DB | PDB；Alphafold/Swiss-Prot；Alphafold/UniProt50-minimal；其他按需 Foldseek databases |
| 本地搜索 | `foldseek search <query_db> <target_db> <out_db> <tmp> --threads ${THREADS:-32} -e 1e-2 --max-seqs 10000` |
| 输出字段 | `query,target,evalue,bits,alntmscore,qtmscore,ttmscore,lddt,prob,qcov,tcov,qaln,taln,taxid,taxname` |
| 在线 endpoint | `https://search.foldseek.com/api/ticket` |
| 在线 mode | `3diaa` |
| 在线数据库 | `afdb50,mgnify_esm30,BFVD,gmgcl_id,pdb100,cath50` |
| 提交间隔 | 8 s |
| 轮询间隔 | 30 s |
| ticket timeout | 120 s |
| rate-limit 后等待 | 300 s |

报告时必须区分：raw alignment、unique target、可下载结构、可评估结构。Foldseek hit 只是结构召回证据。

## 4. Folddisco direct 在线搜索参数

来源脚本：`scripts/19_run_online_folddisco_pet_triads.py`，以及 `results/layered_folddisco/runs/direct_search/online_missing_l0_l3/submissions/online_folddisco_tickets.tsv`。

| 项目 | 参数 |
|---|---|
| endpoint | `https://search.foldseek.com/api/ticket/folddisco` |
| 提交字段 | query PDB file + `motif` + `database[]` |
| 数据库 | `afdb50_folddisco,afdb-proteome_folddisco,esm30_folddisco,BFVD_folddisco,pdb_folddisco` |
| submit sleep | 60 s |
| poll interval | 60 s |
| timeout | 120 s |
| 结果下载 | `/api/result/folddisco/download/{ticket}` |
| motif 策略 | L0、L0+L1、L0+L2、L0+L1+L2 分层 residue list 分别提交 |

Folddisco direct hit 是在线数据库按 motif 返回的结构命中，不是最终功能 pass。由于在线服务有返回限制，宽松 motif 和严格 motif 都需要分别查询并合并去重。

## 5. Folddisco after-recall 打分参数

来源：`results/layered_folddisco/runs/after_recall/query_result_manifest.tsv`、`scripts/21_run_foldseek_only_folddisco_batch.py`。

| 项目 | 参数 |
|---|---|
| index | `folddisco index -p <candidate_structure_dir> -i <index_path> -r -t 16 --id relpath` |
| query | `folddisco query -p <seed_pdb> -q <motif_residues> -i <index_path> -t 16 --per-structure` |
| 输出字段 | `tid,total_match_count,node_count,idf,nres,plddt,max_node_cov,min_rmsd,matching_residues,query_residues` |
| 排序规则 | `max_node_count:desc,idf:desc,min_rmsd:asc` |
| Foldseek-only batch workers | `--download-workers 16` default |
| Foldseek-only timeout | 60 s default |
| 全量模式 | `--limit-per-class 0` 表示不设 class 内限制 |

这个步骤用于把已经由 sequence / Foldseek / direct Folddisco 得到的候选结构重新放进自建结构集合，对每个候选做统一 motif support 评分。

## 6. Layer motif 搜索与评分

Layer 是项目定义，不是 Folddisco 原生唯一标准。

| Layer | 含义 | 作用 |
|---|---|---|
| L0 | 核心催化残基或最关键反应残基 | 核心 anchor |
| L1 | L0 近邻口袋支撑残基 | 描述核心催化环境 |
| L2 | 更外层口袋环境残基 | 描述口袋上下文 |
| L0+L1 | 核心 + 近邻口袋 | 中等严格 |
| L0+L2 | 核心 + 外层环境 | 捕获近邻发生替换但整体环境相似的候选 |
| L0+L1+L2 | 当前最完整项目口袋层 | 用于排序解释，不作为唯一硬筛 |

PETase 示例：`PET_01_6ILW`

| Layer | Residues |
|---|---|
| L0 | `A160,A206,A237` |
| L1 | `A87,A161,A236` |
| L2 | `A159,A185,A238` |
| L0+L1 | `A160,A206,A237,A87,A161,A236` |
| L0+L2 | `A160,A206,A237,A159,A185,A238` |
| L0+L1+L2 | `A160,A206,A237,A87,A161,A236,A159,A185,A238` |

所有 seed 都按同一原则生成具体 residue list。最终 residue list 以服务器上的 layered definition / query manifest 为准。

## 7. ESMFold 缺结构预测参数

来源脚本：`scripts/predict_failed_targets_esmfold.py`。

| 项目 | 参数 |
|---|---|
| 模型 | `esm.pretrained.esmfold_v1()` |
| 设备 | GPU `.eval().cuda()` |
| 输入 | FASTA records |
| 输出 | `<outdir>/pdb/<candidate_id>.pdb` + `esmfold_prediction_manifest.tsv` |
| 当前尼龙酶第一批 | `--sort-by-length --chunk-size 64` |
| 脚本默认 chunk-size | 128 |
| batch 方式 | 逐条 `model.infer_pdb(seq)`，不是多序列 batch |
| 可选限制 | `--max-count`、`--max-len`，默认 0 表示不限制 |

ESMFold 输出结构后必须回流：更新 manifest、跑 Pythia-pocket、整体 RMSD、Folddisco layer scoring、cluster 主表。

## 8. Pythia-pocket 参数

来源脚本：`scripts/run_predicted_pocket_pythia_cosine.py`、`scripts/run_structure_merged_layer_pythia_pocket_cosine.py`。

| 项目 | 参数 |
|---|---|
| 模式 | predicted-pocket：Pythia-pocket 预测 pocket residue probability |
| 阈值 | `--threshold 0.6` |
| 模型 | `tools/src/Pythia/pythia-pocket/model-fl-b.pt` |
| 相似性 | seed pocket embedding vs candidate pocket embedding 的 cosine similarity |
| 无 pocket | 标记 `NO_PREDICTED_POCKET`，不是生物学失败 |
| resume | `--resume` 跳过已完成 pairs |
| flush | `--flush-every 25` |
| CPU threads | `PYTHIA_TORCH_THREADS`，当前应受不超过 64 核限制 |
| layer-specific | `PYTHIA_POCKET_LAYER=L0_L1_L2` 时使用 mapped residues pooled embedding |

Pythia predicted pocket 与手工 L0-L3 layer 是两种口袋定义，主表需要分别保留字段。

## 9. 全局序列相似性参数

最终序列相似性使用 Needleman-Wunsch 全局比对，不使用 BLAST 局部 identity 作为最终散点图颜色。

```text
PID_global = identical_residues / alignment_columns_including_gaps
global_sequence_difference = 1 - PID_global
query_coverage = aligned_query_residues / query_length
target_coverage = aligned_target_residues / target_length
gap_fraction = gap_columns / alignment_columns
```

散点图颜色使用 `PID_global`，避免局部短片段高度相似造成误读。

## 10. Cluster 参数

目前是项目决策参数，执行时还需要把实际命令写入结果日志。

| 层级 | 参数 | 用途 |
|---|---|---|
| cluster90 | 90% sequence identity | 工程去冗余 |
| cluster50 | 50% sequence identity | subfamily 代表选择 |
| cluster30 | 30% sequence identity | family / scaffold SSN |

建议用 MMseqs2 cluster / linclust 或 CD-HIT 执行，并记录 `--min-seq-id`、coverage mode、coverage cutoff。

## 11. 报告中不能混淆的点

1. 序列搜索结果不等于结构搜索结果；有序列但缺结构的候选要先预测结构。
2. Foldseek hit 不等于 Folddisco hit；Foldseek 是整体结构召回，Folddisco 是 motif/geometry 召回或评分。
3. direct Folddisco 是在线库召回；after-recall Folddisco 是对我们已有候选结构重新评分。
4. `NO_HIT`、`NO_PREDICTED_POCKET`、`NOT_EVALUATED` 必须分开。
5. L0-L3 是项目分层定义；Folddisco 自身只负责按给定 motif 返回结构命中。
6. PETase 和尼龙酶后续排序不能只看单一分数，需要结合 global PID、global RMSD、pocket cosine、layer support、cluster 和 nylonase oligomer/interface 证据。
