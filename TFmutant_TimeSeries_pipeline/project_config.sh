#!/bin/bash
# 文件名: project_config.sh

# ================= 1. 项目基础配置 =================
PROJECT_NAME="CrDof_Phosphate_TimeCourse_RNAseq"  # 项目名称
BASE_DIR="/home/cfy/Data/CrDof_project"           # 项目根目录
METADATA_FILE="/home/cfy/Data/CrDof_project/RNAseq/samples.txt"    # 样本信息表

# ================= 2. 实验设计参数 (核心控制) =================
# 控制条件 (基准基因型，用于Logic 2比较)
CONTROL_CONDITION="WT"

# 控制时间点 (基准时间，用于Logic 1比较)
# 注意：必须与samples.txt中的时间格式严格一致 (例如 00h 或 00)
CONTROL_TIME="00h"

# 统计阈值
FDR_CUTOFF=0.05
LFC_THRESHOLD=1.0

# ================= 3. 软件与资源 =================

# --- 核心并行数 ---
# 控制 P1(HISAT2) 和 P3(StringTie) 同时处理的样本数
PARALLEL_JOBS=8

# --- P1: HISAT2 ---
HISAT2_THREADS=8
SAMTOOLS_THREADS=2

# --- P2: featureCounts ---
# 单任务逻辑: I/O密集型 (最大64)
FEATURECOUNTS_THREADS=24

# --- P3: StringTie ---
STRINGTIE_THREADS=10

# 软件路径
HISAT2_BIN="/home/jxq/biosoft/hisat2-2.2.1/hisat2"
SAMTOOLS_BIN="samtools"
FEATURECOUNTS_BIN="featureCounts"
STRINGTIE_BIN="/home/cfy/biosoft/stringtie/stringtie-2.2.3.Linux_x86_64/stringtie"

# ================= 4. 参考基因组 =================
GENOME_FASTA="/home/cfy/Data/CrDof_project/ref/PhytozomeV13/CreinhardtiiCC_4532/v6.1/CreinhardtiiCC_4532_707_v6.0.fa"
INDEX_PREFIX="/home/cfy/Data/CrDof_project/ref/PhytozomeV13/CreinhardtiiCC_4532/v6.1/hisat2/Creinhardtii_v6.1.hisat2"
GTF_FILE="/home/cfy/Data/CrDof_project/ref/PhytozomeV13/CreinhardtiiCC_4532/v6.1/CreinhardtiiCC_4532_707_v6.1.gene_exons.gtf"
GFF_FILE="/home/cfy/Data/CrDof_project/ref/PhytozomeV13/CreinhardtiiCC_4532/v6.1/CreinhardtiiCC_4532_707_v6.1.gene_exons.gff3"

# 基因ID清理正则 (可选, 留空则不做任何剥离)
# 衣藻示例: GENE_ID_STRIP_PATTERN="_4532(\\.v6\\.1)?$"
# 地钱: 通常不需要, 留空即可
GENE_ID_STRIP_PATTERN="_4532.*$"
# ================= 5. 可视化工程化配置 (Visualization) =================

# --- A. Condition 颜色 (用于分面图/Venn图) ---
# 核心逻辑：用于区分基因型。
# 格式: Key:Value (Value可以是HEX或颜色名)
STYLE_COND_COLORS="WT:#2E86AB,dof:#A23B72"

# --- B. Condition 形状 (用于全局图) ---
# 核心逻辑：用于区分基因型。因为要配合 Time 的颜色，建议使用实心形状。
# R语言形状代码推荐: 16=圆, 17=三角, 15=方, 18=菱形
STYLE_COND_SHAPES="WT:16,dof:17"

# --- C. Time 颜色 (用于全局图) ---
# 核心逻辑：用于展示时序。建议使用渐变色系（如浅蓝到深蓝）或彩虹色。
# 必须覆盖所有时间点，否则脚本会自动生成。
STYLE_TIME_COLORS="00h:#E5D65C,03h:#E9B730,06h:#E9992C,12h:#E47B43,24h:#DA595B"

# --- D. 标签美化 ---
# 格式: RawName:DisplayName
STYLE_LABELS="WT:WT,dof:italic(dof),00h:0hr,03h:3hr,06h:6hr,12h:12hr,24h:24hr"

# --- E. 统计阈值 ---
VIS_FDR=0.01
VIS_LFC=1.0

# --- F. 排序逻辑 ---
# 用于控制 UpSet 图侧边栏的显示顺序
CONDITION_ORDER="WT,dof"

# 时间排序通常自动处理，但如果需要强制顺序：
# TIME_ORDER="00h,03h,06h,12h,24h"


# ================= 6. Pipeline 5: Mfuzz 时序聚类配置 (v11.1 Standardized) =================

# --- A. 条件标签与显示配置 ---
# 逻辑：Display 变量支持 "italic(Name)" 语法，Shell脚本会自动解析并传递给R脚本

# 1. 真实标签 (必须与 samples.txt 中的 Condition 列严格一致)
CONDITION_WT="WT"
CONDITION_MUT="dof"          # 统一变量名逻辑，对应 R 脚本的 condition_mut

# 2. 图形显示名称 (支持 italic() 语法)
CONDITION_WT_DISPLAY="WT"              # WT 通常不斜体
CONDITION_MUT_DISPLAY="italic(dof)"    # 基因名通常斜体，脚本会自动识别 wrapper

# --- B. 聚类数控制 (可选) ---
# 留空则自动通过 PP-XB Algorithm 选择最优聚类数
# MFUZZ_C_WT=12
# MFUZZ_C_MUT=10   # 对应 R 脚本的 c_mut

# --- C. 聚类数搜索范围 ---
MFUZZ_C_RANGE_MIN=6
MFUZZ_C_RANGE_MAX=24

# --- D. 模糊因子 m 的约束范围 ---
MFUZZ_M_LOWER=1.25
MFUZZ_M_UPPER=2.00

# --- E. 聚类成员阈值 ---
MFUZZ_MEMBERSHIP_CUTOFF=0.5   # 核心基因的隶属度阈值
MFUZZ_MIN_STD=0.4            # 最小标准差过滤 (建议 0.25-0.5,控制基因数目在3k-5k)

# --- F. DDTW 成本函数权重 ---
# 用于 WT-MUT 聚类匹配 (Pipeline 5a)
MFUZZ_W_SHAPE=0.5             # 形状权重
MFUZZ_W_AMP=0.3               # 振幅权重
MFUZZ_W_PHASE=0.2             # 相位权重

# --- G. 分类阈值 ---
# 用于判断 Conserved/Shape_Altered/Amp_Altered/Phase_Shifted
MFUZZ_CLASSIFY_SHAPE=0.15
MFUZZ_CLASSIFY_AMP=0.15
MFUZZ_CLASSIFY_PHASE=0.25

# --- H. Pipeline5b 子聚类参数 ---
MFUZZ_P5B_C_MIN=2
MFUZZ_P5B_C_MAX=6
MFUZZ_P5B_MIN_GENES=10

# --- I. 运行模式控制 ---
RUN_PIPELINE5A=TRUE           # 运行 WT-MUT 匹配分析
RUN_PIPELINE5B=TRUE           # 运行分化分析

# ================= 7. 高级配置 (v11.0 PP-XB Algorithm) =================

# --- Cluster Selection Parameters (PP-XB Algorithm) ---
# v11.0 uses Kwon Index instead of XB to solve monotonic decreasing issue
# Prominence threshold for valley detection
XB_GLOBAL_PROMINENCE=0.10      # Pipeline5a: genome-wide clustering
XB_SUBCLUSTER_PROMINENCE=0.15  # Pipeline5b: sub-clustering (more conservative)
XB_KNEEDLE_S=1.0               # Kneedle sensitivity (legacy, rarely used)

# --- DDTW Algorithm Parameters ---
DDTW_INTERP_POINTS=100                # Interpolation points for DTW
DDTW_DERIVATIVE_METHOD="central"      # "central", "forward", "backward"
# DDTW_WINDOW=NULL                    # Sakoe-Chiba bandwidth (NULL=unconstrained)

# --- Sample Size Constraint (v11.0) ---
# c_max <= min(sqrt(N), N/10) for Pipeline5a
# Pipeline5b uses more conservative constraints for small samples

# ================= 8. Pipeline 5c: 聚类时序区分  =================
RUN_PIPELINE5C=TRUE

# --- A. 阶段区间 (Stage Ranges) ---
# 核心修改: 使用 "Start-End" 的格式定义区间
# 逻辑: 只要 Cluster 的加权时间落在该区间内，就属于该阶段
# 注意: 区间必须覆盖所有时间点
# 格式: 逗号分隔的 "Min-Max" 字符串，数字格式需严格对应samples.txt文件中的数值部分，mfuzz脚本处理后会去除字母前后缀
P5C_STAGE_RANGES="00-03,03-12,12-24"

# --- B. 阶段命名 (Stage Names) ---
# 必须与上面的区间数量严格一致
P5C_STAGE_NAMES="Early,Middle,Late"

# --- C. 阶段颜色 (Stage Colors) ---
# 必须与上面的区间数量严格一致
P5C_STAGE_COLORS="#D1C4E9,#7E57C2,#311B92"

# 其他参数保持不变
P5C_N_PSEUDOTIME=100
P5C_LOESS_SPAN=0.4

# ================= 9. Pipeline 6: pCRE Discovery & Machine Learning =================

# --- A. 通用控制 ---
RUN_PIPELINE6A=TRUE     # 运行 6a (特征提取与富集)
RUN_PIPELINE6B=TRUE     # 运行 6b (随机森林建模)
RUN_PIPELINE6C=TRUE     # 运行 6c (Motif聚类与合并)
RUN_PIPELINE6D=TRUE     # 运行 6d (Publication-Quality Visualization v2.0)
RUN_PIPELINE6E=TRUE     # 运行 6e (标签置换检验)
P6_GENOTYPE="both"      # 处理基因型: WT, Mutant, 或 both
P6_CORES=100             # 并行核心数

# --- B. Pipeline 6a: 序列特征提取与富集 ---
# 区域模式: TRUE=分区模式(Core/Proximal/Distal), FALSE=全长模式
P6_PARTITIONED=TRUE

# 全长模式参数 (仅当 P6_PARTITIONED=FALSE 时生效)
P6_FULL_UP=1500
P6_FULL_DOWN=200

# 分区模式参数 (仅当 P6_PARTITIONED=TRUE 时生效)
P6_CORE_START=-50
P6_CORE_END=200
P6_PROX_START=-500
P6_PROX_END=-50
P6_DIST_START=-1500
P6_DIST_END=-500

# 特征提取参数
P6_K_VALUES="6"         # k-mer长度, 逗号分隔 (例如 "6" 或 "6,7")
P6_MIN_OCCURRENCE=5     # k-mer 最小出现次数
P6_TOP_N=100            # 每个Cluster/Region提取的特征数 (无FDR过滤，保留Top N)

# 富集背景策略
# P6_ENRICH_LOG2FC=""   # 留空则使用全基因组背景，设置数值则使用非响应基因背景
P6_ENRICH_GC_MATCH=FALSE # 是否在富集分析阶段进行GC匹配 (通常选FALSE以提高速度，6b会有GC匹配)

# ML定义参数
P6_ML_LOG2FC=0.8        # 用于定义ML非响应基因(负样本池)的阈值

# --- C. Pipeline 6b: Random Forest 建模 ---
P6_RF_NTREE=500         # 决策树数量
P6_RF_REPLICATES=100    # 重复次数 (Replicates)
P6_RF_FOLDS=10          # 交叉验证折数 (CV Folds)

# 负样本策略
P6_RF_GC_MATCH=TRUE    # 是否在RF训练时对负样本进行GC匹配 (推荐TRUE，但需消耗更多资源)
P6_RF_GC_BINS=10        # GC分箱数
P6_RF_NEG_RATIO=3       # 负样本:正样本 比例

# ================= Pipeline 6c v3.0: Motif Integration & Annotation =================
# 核心功能: Levenshtein聚类 + Shift-Aware PWM合并 + PCC数据库注释
# 注意: v3.0 取代了旧版 6c 和 6d 的功能

P6_CLUST_DIST=2           # Levenshtein 编辑距离阈值 (定义 Motif Family)
P6_CLUST_MIN_OVERLAP=4    # Shift-Aware 合并时的最小重叠长度 (建议 4)
P6_CLUST_TRIM_IC=0.25     # PWM 边缘修剪的信息含量(IC)阈值 (防止合并后边缘出现噪声)

# 资源控制
P6_CORES=60               # 并行核心数 (universalmotif)
P6_GENOTYPE="both"        # 处理基因型: WT, Mutant, 或 both

MEME_DB_PATH="/home/cfy/biosoft/motif_databases/ARABD/ArabidopsisDAPv1.meme"

# ================= Pipeline 6d v2.0: Publication-Quality Visualization =================
# 用于生成 AUC 柱状图、正交下三角热图以及每个 Cluster 的详细 Motif 报告 (DB-Anchored 对齐)

# --- A. 显示标签 (Label Parsing) ---
# 逻辑: 复用 Pipeline 5 中已定义的 CONDITION_WT_DISPLAY / CONDITION_MUT_DISPLAY
# 支持 italic(text) 语法，R脚本会自动解析为斜体表达式
# 如果需要独立于 Pipeline 5 的标签，可在此覆盖:
# P6_VIS_WT_LABEL="WT"
# P6_VIS_MUT_LABEL="italic(dof)"

# --- B. 颜色一致性 ---
# 逻辑: 从 STYLE_COND_COLORS 中自动解析 WT / Mutant 的颜色
# 也可在此手动覆盖:
# P6_VIS_WT_COLOR="#2E86AB"
# P6_VIS_MUT_COLOR="#A23B72"

# --- C. 图1: AUC-ROC 柱状图 ---
P6_VIS_AUC_THRESHOLD=0.7    # AUC 参考阈值线

# --- D. 图2: 正交相似性热图 (下三角布局) ---
P6_VIS_HEATMAP_METRIC="overlap" # 相似性指标: 'overlap' (Simpson) 或 'jaccard'
P6_VIS_TOP_N_FEATURES=50       # 仅使用每个Cluster中Importance排名前N的特征计算相似性 (提升视觉对比度)

# --- E. 图3: Per-Cluster Scatter + DB-Anchored Motif Logos ---
P6_VIS_TOP_N=3                  # 每个Cluster报告中展示Top N个最重要的Motif
P6_VIS_TARGET_KEYWORD="dof"        # 目标关键词 (逗号分隔,例如 "dof,WRKY,bHLH")，在scatter中高亮匹配motif; 留空则不启用
P6_VIS_TARGET_SCOPE="ALL"          # 关键词高亮作用域: "ALL"(所有Cluster), "WT"(仅WT), "MUT"(仅Mutant), 或逗号分隔的内部Cluster ID (如 "WT_C1,MUT_C2")
                                    # 注意: 使用内部ID (WT_C1/MUT_C1), 不要使用显示标签 (如 dof_C1)

# --- F. 绘图外观 ---
P6_VIS_FIG_WIDTH=12         # PDF 宽度
P6_VIS_FIG_HEIGHT=8         # PDF 高度

# ================= Pipeline 6e: Permutation Baseline Test =================
RUN_PIPELINE6E=TRUE
P6_PERM_N=1000              # 置换次数
P6_PERM_REPLICATES=10       # 每次置换的简化 CV 重复数 (6b 用 100, 这里 10 足够)
P6_PERM_FOLDS=5             # 每次重复的 CV 折数 (6b 用 10, 这里 5 足够)
P6_PERM_NTREE=200           # 置换 RF 树数 (6b 用 500, 置换测试不需要那么多)
P6_PERM_ALPHA=0.05          # 显著性阈值

# 注: 负样本策略 (GC_MATCH, GC_BINS, NEG_RATIO) 复用 6b 的配置以保持一致性


# ========================================================
# Pipeline 7: GRN Inference Configuration
# ========================================================

# --- 核心输入文件 ---
GRN_TF_LIST_FILE="/home/cfy/Data/CrDof_project/ref/Cre_TF_list_v6.1.txt"
GRN_HIGHLIGHT_FILE="/home/cfy/Data/CrDof_project/ref/Cr_PSR_gene_v6.1.txt"

# --- 算法参数 ---
GRN_LAG=1
GRN_N_PERMUTATIONS=1000
GRN_FDR_CUTOFF=0.05
GRN_FDR_CUTOFF_DIFF=0.05
GRN_FDR_MUT_LOOSE=0.10
GRN_COR_METHOD="spearman"
GRN_OUTPUT_MODE="significant"
GRN_N_CORES=100
GRN_SEED=42

# --- v1.3 新增参数 ---
GRN_TARGET_FILTER="mfuzz"
GRN_PREFILTER_R=0

# ========== Pipeline 7b/7c 新增参数 ==========
GRN_FOCUS_TF_FAMILY="Dof"
GRN_FOCUS_GENE="Cre12.g521150_4532.v6.1"
GRN_CLUSTER_PAIRS_FILE=""
GRN_PCRE_TF_MAP=""
GRN_BACKGROUND_MODE="mfuzz"
GRN_W_DIFFCOEXP="1.0"
GRN_W_CLUSTER="1.0"
GRN_W_PCRE="1.5"
GRN_ENRICHMENT_METHOD="fisher"
GRN_MIN_EDGES="5"
GRN_EDGE_CLASSES="Lost,Gained,Rewired,Conserved"

# ========== 7c 可视化参数 ==========
GRN_FIG_FORMAT="both"
GRN_FIG_DPI="300"
GRN_FIG_SCALE="1.0"
GRN_PLOTS="heatmap,volcano,barplot,network,circos,outdegree,alluvial"
GRN_LABEL_TOP_N="15"
GRN_LABEL_FOCUS="TRUE"
GRN_SHOW_HEATMAP_STARS="FALSE"

# --- Fig4 网络图参数 (v4.4) ---
GRN_NETWORK_MAX_EDGES="500"
GRN_NETWORK_KEEP_DEG1_ALL="FALSE"
GRN_NETWORK_KEEP_DEG1_FOCUS="TRUE"
GRN_NETWORK_PRUNE="TRUE"
GRN_NETWORK_LAYOUT="kk"
GRN_NETWORK_SEED="42"
GRN_NETWORK_FIG_W="18"
GRN_NETWORK_FIG_H="14"
GRN_NETWORK_ARROW_MM="1.5"
GRN_NETWORK_EDGE_WIDTH_MIN="0.5"
GRN_NETWORK_EDGE_WIDTH_MAX="2.0"
GRN_NETWORK_NODE_SIZE_MIN="4"
GRN_NETWORK_NODE_SIZE_MAX="12"

# --- v4.5 拓扑约束 ---
# 注意: 以下两个步数参数作用于不同阶段
# max_secondary_steps: Step 1b2, 过滤二级锚点 (基于全基因组图到一级锚点距离)
# ego_steps:           Step 1d2, ego 过滤 (基于已选子图到 focus_gene 距离)
# 当两者相同时 1b2 几乎被 1d2 覆盖; 建议 max_secondary >= ego 以避免双重淘汰
GRN_NETWORK_MAX_SECONDARY_STEPS="2"
GRN_NETWORK_INTRA_COMPLETION="FALSE"   # 默认关闭，恢复 v4.4 hub-spoke 结构
GRN_NETWORK_SHARED_TARGETS="FALSE"     
GRN_NETWORK_EGO_STEPS="2"              # 仅保留 focus_gene N 步内的节点

# --- Fig5 Circos 参数 ---
GRN_CIRCOS_MAX_EDGES="500"
GRN_CIRCOS_TOP_EDGES="500"

# ========== 7d 自动生成 highlight 列表参数 (v4.6) ==========
GRN_AUTO_HIGHLIGHT_CLUSTERS="WT_C3"
GRN_AUTO_HIGHLIGHT_MAX_PER_CLUSTER=50
# MECS 质量过滤门槛 (v4.6 新增)
GRN_AUTO_HIGHLIGHT_MIN_TIER=3             # 基因须参与 Tier <= N 的边
GRN_AUTO_HIGHLIGHT_MIN_EDGE_COUNT=1       # 基因须参与的最小边数
# 已移除: GRN_AUTO_HIGHLIGHT_HUB_N (v4.6 不再提取 hub target)