#!/bin/bash
# script_name: 07_run_grn_all.sh
# 功能: 运行 Pipeline 7 GRN 完整流程 (7a + 7b + 7c)
# 版本: v1.0 (整合版)
# 用法: 
#   默认顺序执行 a→b→c:  ./07_run_grn_all.sh
#   单独执行某步骤:       ./07_run_grn_all.sh -a / -b / -c
#   组合执行:             ./07_run_grn_all.sh -a -b

# ========================================================
# 1. 自动获取脚本所在目录 & 加载配置
# ========================================================
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

if [ -f "${SCRIPT_DIR}/project_config.sh" ]; then
    source "${SCRIPT_DIR}/project_config.sh"
else
    echo "❌ Error: project_config.sh not found in ${SCRIPT_DIR}"
    exit 1
fi

# ========================================================
# 2. 步骤控制标志 (默认全部执行)
# ========================================================
RUN_7A="FALSE"
RUN_7B="FALSE"
RUN_7C="FALSE"
RUN_7C_NET="FALSE"    # v4.5: 只执行 7c 中的 Fig4 网络图
STEP_SPECIFIED="FALSE"

# ========================================================
# 3. 输入输出路径 (全局)
# ========================================================
# Pipeline 4/5 输入
RDATA_P4="${BASE_DIR}/04_DESeq2/RData/Pipeline4_Logic1_TimeCourse_for_Mfuzz.RData"
RDATA_P5="${BASE_DIR}/05_Mfuzz/rdata/Mfuzz_v11_Results.RData"

# Pipeline 7 中间产物
RDATA_7A="${BASE_DIR}/07_GRN/rdata/Pipeline7a_DiffCoexp_Results.RData"
RDATA_7B="${BASE_DIR}/07_GRN/rdata/Pipeline7b_Discovery_Results.RData"

# 输出目录
OUT_DIR="${BASE_DIR}/07_GRN"

# ========================================================
# 4. 从 Config 读取所有 GRN 参数 (带默认值)
# ========================================================

# === 核心输入文件 ===
GRN_TF_LIST="${GRN_TF_LIST_FILE:-${BASE_DIR}/00_metadata/TF_list.txt}"
GRN_HIGHLIGHT="${GRN_HIGHLIGHT_FILE:-}"

# === 条件标签 ===
COND_WT="${CONDITION_WT:-WT}"
COND_MUT="${CONDITION_MUT:-dof}"

# === Pipeline 7a 参数 ===
LAG="${GRN_LAG:-1}"
N_PERM="${GRN_N_PERMUTATIONS:-1000}"
FDR="${GRN_FDR_CUTOFF:-0.05}"
FDR_DIFF="${GRN_FDR_CUTOFF_DIFF:-0.05}"
FDR_MUT_LOOSE="${GRN_FDR_MUT_LOOSE:-0.10}"
COR_METHOD="${GRN_COR_METHOD:-spearman}"
OUTPUT_MODE="${GRN_OUTPUT_MODE:-significant}"
N_CORES="${GRN_N_CORES:-1}"
SEED="${GRN_SEED:-42}"
TARGET_FILTER="${GRN_TARGET_FILTER:-mfuzz}"
PREFILTER_R="${GRN_PREFILTER_R:-0}"

# === Pipeline 7b 参数 ===
FOCUS_TF_FAMILY="${GRN_FOCUS_TF_FAMILY:-none}"
FOCUS_GENE="${GRN_FOCUS_GENE:-none}"
CLUSTER_PAIRS_FILE="${GRN_CLUSTER_PAIRS_FILE:-}"
PCRE_TF_MAP="${GRN_PCRE_TF_MAP:-}"
PCRE_INTEGRATED_DIR="${GRN_PCRE_INTEGRATED_DIR:-${BASE_DIR}/06_pCRE/Pipeline6c_v3.0_MotifIntegration}"
PCRE_FAMILY_ALIAS="${GRN_PCRE_FAMILY_ALIAS:-}"
PCRE_MIN_PCC="${GRN_PCRE_MIN_PCC:-0}"
BACKGROUND_MODE="${GRN_BACKGROUND_MODE:-mfuzz}"
ENRICHMENT_METHOD="${GRN_ENRICHMENT_METHOD:-fisher}"
MIN_EDGES="${GRN_MIN_EDGES:-5}"
W_DIFFCOEXP="${GRN_W_DIFFCOEXP:-1.0}"
W_CLUSTER="${GRN_W_CLUSTER:-1.0}"
W_PCRE="${GRN_W_PCRE:-1.5}"
EDGE_CLASSES="${GRN_EDGE_CLASSES:-Lost,Gained,Rewired,Conserved}"
HIGHLIGHT_FILE="${GRN_HIGHLIGHT_FILE:-}"

# === Pipeline 7c 参数 ===
FIG_FORMAT="${GRN_FIG_FORMAT:-both}"
FIG_DPI="${GRN_FIG_DPI:-300}"
FIG_SCALE="${GRN_FIG_SCALE:-1.0}"
PLOTS="${GRN_PLOTS:-heatmap,volcano,barplot,network,circos,outdegree,evidence_heatmap}"
LABEL_TOP_N="${GRN_LABEL_TOP_N:-15}"
LABEL_FOCUS="${GRN_LABEL_FOCUS:-TRUE}"
SHOW_HEATMAP_STARS="${GRN_SHOW_HEATMAP_STARS:-FALSE}"

NETWORK_MAX_EDGES="${GRN_NETWORK_MAX_EDGES:-600}"
NETWORK_KEEP_DEG1_ALL="${GRN_NETWORK_KEEP_DEG1_ALL:-FALSE}"
NETWORK_KEEP_DEG1_FOCUS="${GRN_NETWORK_KEEP_DEG1_FOCUS:-FALSE}"
NETWORK_LAYOUT="${GRN_NETWORK_LAYOUT:-kk}"
NETWORK_SEED="${GRN_NETWORK_SEED:-42}"
NETWORK_PRUNE="${GRN_NETWORK_PRUNE:-TRUE}"
NETWORK_FIG_W="${GRN_NETWORK_FIG_W:-18}"
NETWORK_FIG_H="${GRN_NETWORK_FIG_H:-14}"
NETWORK_ARROW_MM="${GRN_NETWORK_ARROW_MM:-1.5}"
NETWORK_EDGE_WIDTH_MIN="${GRN_NETWORK_EDGE_WIDTH_MIN:-0.5}"
NETWORK_EDGE_WIDTH_MAX="${GRN_NETWORK_EDGE_WIDTH_MAX:-2.0}"
NETWORK_NODE_SIZE_MIN="${GRN_NETWORK_NODE_SIZE_MIN:-4}"
NETWORK_NODE_SIZE_MAX="${GRN_NETWORK_NODE_SIZE_MAX:-12}"

CIRCOS_MAX_EDGES="${GRN_CIRCOS_MAX_EDGES:-500}"
CIRCOS_TOP_EDGES="${GRN_CIRCOS_TOP_EDGES:-300}"

DESEQ2_RDATA="${GRN_DESEQ2_RDATA:-${BASE_DIR}/04_DESeq2/RData/${PROJECT_NAME}_DESeq2_objects.RData}"
CTRL_TIME="${CONTROL_TIME:-00h}"

CLUSTER_COLORS_FILE="${GRN_CLUSTER_COLORS_FILE:-}"
EDGE_CLASS_COLORS_FILE="${GRN_EDGE_CLASS_COLORS_FILE:-}"

AUTO_HL_CLUSTERS="${GRN_AUTO_HIGHLIGHT_CLUSTERS:-}"
# v4.6: MECS quality filter parameters (replaces AUTO_HL_HUB_N)
AUTO_HL_MIN_TIER="${GRN_AUTO_HIGHLIGHT_MIN_TIER:-3}"
AUTO_HL_MIN_EDGE_COUNT="${GRN_AUTO_HIGHLIGHT_MIN_EDGE_COUNT:-1}"

CLI_EGO_STEPS=""

# v4.5: 分级锚点 + 共享靶标
NETWORK_MAX_SECONDARY_STEPS="${GRN_NETWORK_MAX_SECONDARY_STEPS:-2}"
NETWORK_INTRA_COMPLETION="${GRN_NETWORK_INTRA_COMPLETION:-FALSE}"
NETWORK_SHARED_TARGETS="${GRN_NETWORK_SHARED_TARGETS:-TRUE}"

# ========================================================
# 5. 帮助信息
# ========================================================
usage() {
    echo "Usage: $0 [step flags] [options]"
    echo ""
    echo "Pipeline 7: Complete GRN Inference Workflow"
    echo "  7a: Differential Co-expression Analysis"
    echo "  7b: Cluster-Constrained Discovery + MECS Scoring"
    echo "  7c: Network Visualization"
    echo ""
    echo "Step Flags (默认全部执行):"
    echo "  -a, --step-a          Run Pipeline 7a only"
    echo "  -b, --step-b          Run Pipeline 7b only"
    echo "  -c, --step-c          Run Pipeline 7c only"
    echo "  -n, --step-net        Run Pipeline 7c network plot (Fig4) only"
    echo "  (可组合使用, 如: -a -b 仅执行 7a 和 7b)"
    echo ""
    echo "Global Options:"
    echo "  --out_dir PATH        Override output directory"
    echo "  --highlight PATH      Highlight genes file"
    echo ""
    echo "Pipeline 7a Options:"
    echo "  --rdata_p4 PATH       Override Pipeline4 RData path"
    echo "  --rdata_p5 PATH       Override Pipeline5 RData path (or 'none' to skip)"
    echo "  --tf_list PATH        Override TF list file path"
    echo "  --lag N               Lag steps (default: ${LAG})"
    echo "  --n_perm N            Permutation count (default: ${N_PERM})"
    echo "  --fdr N               FDR cutoff (default: ${FDR})"
    echo "  --fdr_diff N          Differential FDR cutoff (default: ${FDR_DIFF})"
    echo "  --output_mode MODE    significant/classified/full (default: ${OUTPUT_MODE})"
    echo "  --n_cores N           Parallel cores (default: ${N_CORES})"
    echo "  --seed N              Random seed (default: ${SEED})"
    echo "  --target_filter MODE  all/mfuzz/non_conserved (default: ${TARGET_FILTER})"
    echo "  --prefilter_r N       |r| pre-filter threshold, 0=disabled (default: ${PREFILTER_R})"
    echo ""
    echo "Pipeline 7b Options:"
    echo "  --rdata_7a PATH       Override Pipeline7a RData path"
    echo "  --focus_tf_family N   TF family to focus (default: ${FOCUS_TF_FAMILY})"
    echo "  --focus_gene ID       Gene ID to focus (default: ${FOCUS_GENE})"
    echo "  --cluster_pairs PATH  Cluster pair definition TSV"
    echo "  --pcre_tf_map PATH    pCRE-TF family mapping TSV (manual, highest priority)"
    echo "  --pcre_integrated_dir PATH  Pipeline 6c output dir for auto pCRE parsing"
    echo "  --pcre_family_alias STR     Manual alias e.g. 'C2C2dof=DOF,C2C2gata=GATA'"
    echo "  --pcre_min_pcc N      Min PCC for pCRE DB match (default: ${PCRE_MIN_PCC})"
    echo "  --background_mode M   mfuzz/all_deg/genome (default: ${BACKGROUND_MODE})"
    echo "  --enrichment_method M fisher/hypergeometric (default: ${ENRICHMENT_METHOD})"
    echo "  --min_edges N         Minimum edges for enrichment (default: ${MIN_EDGES})"
    echo "  --w_diffcoexp N       MECS S1 weight (default: ${W_DIFFCOEXP})"
    echo "  --w_cluster N         MECS S2 weight (default: ${W_CLUSTER})"
    echo "  --w_pcre N            MECS S3 weight (default: ${W_PCRE})"
    echo "  --edge_classes LIST   Comma-separated edge classes (default: ${EDGE_CLASSES})"
    echo ""
    echo "Pipeline 7c Options:"
    echo "  --rdata_7b PATH       Override Pipeline7b RData path"
    echo "  --fig_format FMT      pdf/png/both (default: ${FIG_FORMAT})"
    echo "  --fig_dpi N           DPI for PNG (default: ${FIG_DPI})"
    echo "  --fig_scale N         Figure scale (default: ${FIG_SCALE})"
    echo "  --plots LIST          Comma-separated plot types (default: ${PLOTS})"
    echo "  --label_top_n N       Top N nodes to label (default: ${LABEL_TOP_N})"
    echo "  --label_focus BOOL    Label focus genes (default: ${LABEL_FOCUS})"
    echo "  --show_heatmap_stars  Show stars on heatmap (default: ${SHOW_HEATMAP_STARS})"
    echo "  --network_max_edges N Max edges in network (default: ${NETWORK_MAX_EDGES})"
    echo "  --network_keep_deg1_all BOOL  Keep all TF deg=1 targets (default: ${NETWORK_KEEP_DEG1_ALL})"
    echo "  --network_keep_deg1_focus BOOL Keep focus gene deg=1 targets (default: ${NETWORK_KEEP_DEG1_FOCUS})"
    echo "  --network_layout L    Layout algorithm: kk/fr/stress (default: ${NETWORK_LAYOUT})"
    echo "  --network_seed N      Network layout seed (default: ${NETWORK_SEED})"
    echo "  --network_prune BOOL  Prune degree-1 targets (default: ${NETWORK_PRUNE})"
    echo "  --network_fig_w N     Network figure width (default: ${NETWORK_FIG_W})"
    echo "  --network_fig_h N     Network figure height (default: ${NETWORK_FIG_H})"
    echo "  --network_arrow_mm N  Arrow size in mm (default: ${NETWORK_ARROW_MM})"
    echo "  --network_edge_width_min N  Min edge width (default: ${NETWORK_EDGE_WIDTH_MIN})"
    echo "  --network_edge_width_max N  Max edge width (default: ${NETWORK_EDGE_WIDTH_MAX})"
    echo "  --network_node_size_min N   Min node size (default: ${NETWORK_NODE_SIZE_MIN})"
    echo "  --network_node_size_max N   Max node size (default: ${NETWORK_NODE_SIZE_MAX})"
    echo "  --network_ego_steps N Ego-network topology steps (default: auto)"
    echo "  --network_max_secondary_steps N  Max topological distance for secondary anchors (default: 2)"
    echo "  --network_intra_completion BOOL  Enable intra-network edge completion (default: FALSE)"
    echo "  --network_shared_targets BOOL   Enable shared target discovery (default: TRUE)"
    echo "  --gene_id_strip REGEX  Gene ID cleanup regex (default: from config GENE_ID_STRIP_PATTERN)"
    echo "  --circos_max_edges N  Max edges in circos (default: ${CIRCOS_MAX_EDGES})"
    echo "  --circos_top_edges N  Top edges in circos (default: ${CIRCOS_TOP_EDGES})"
    echo "  --deseq2_rdata PATH   DESeq2 RData for evidence heatmap"
    echo "  --control_time T      Control time point (default: ${CTRL_TIME})"
    echo "  --cluster_colors PATH Cluster colors file"
    echo "  --edge_class_colors PATH  Edge class colors file"
    echo ""
    echo "  -h, --help            Show this help"
    exit 1
}

# ========================================================
# 6. 命令行参数解析
# ========================================================
while [[ $# -gt 0 ]]; do
    case $1 in
        # === 步骤控制 ===
        -a|--step-a)          RUN_7A="TRUE"; STEP_SPECIFIED="TRUE"; shift ;;
        -b|--step-b)          RUN_7B="TRUE"; STEP_SPECIFIED="TRUE"; shift ;;
        -c|--step-c)          RUN_7C="TRUE"; STEP_SPECIFIED="TRUE"; shift ;;
        -n|--step-net)        RUN_7C_NET="TRUE"; STEP_SPECIFIED="TRUE"; shift ;;
        
        # === 全局参数 ===
        --out_dir)            OUT_DIR="$2";              shift 2 ;;
        --highlight)          GRN_HIGHLIGHT="$2"; HIGHLIGHT_FILE="$2"; shift 2 ;;
        
        # === Pipeline 7a 参数 ===
        --rdata_p4)           RDATA_P4="$2";             shift 2 ;;
        --rdata_p5)           RDATA_P5="$2";             shift 2 ;;
        --tf_list)            GRN_TF_LIST="$2";          shift 2 ;;
        --lag)                LAG="$2";                  shift 2 ;;
        --n_perm)             N_PERM="$2";               shift 2 ;;
        --fdr)                FDR="$2";                  shift 2 ;;
        --fdr_diff)           FDR_DIFF="$2";             shift 2 ;;
        --output_mode)        OUTPUT_MODE="$2";          shift 2 ;;
        --n_cores)            N_CORES="$2";              shift 2 ;;
        --seed)               SEED="$2";                 shift 2 ;;
        --target_filter)      TARGET_FILTER="$2";        shift 2 ;;
        --prefilter_r)        PREFILTER_R="$2";          shift 2 ;;
        
        # === Pipeline 7b 参数 ===
        --rdata_7a)           RDATA_7A="$2";             shift 2 ;;
        --focus_tf_family)    FOCUS_TF_FAMILY="$2";      shift 2 ;;
        --focus_gene)         FOCUS_GENE="$2";           shift 2 ;;
        --cluster_pairs)      CLUSTER_PAIRS_FILE="$2";   shift 2 ;;
        --pcre_tf_map)        PCRE_TF_MAP="$2";          shift 2 ;;
        --pcre_integrated_dir) PCRE_INTEGRATED_DIR="$2"; shift 2 ;;
        --pcre_family_alias)  PCRE_FAMILY_ALIAS="$2";    shift 2 ;;
        --pcre_min_pcc)       PCRE_MIN_PCC="$2";         shift 2 ;;
        --background_mode)    BACKGROUND_MODE="$2";      shift 2 ;;
        --enrichment_method)  ENRICHMENT_METHOD="$2";    shift 2 ;;
        --min_edges)          MIN_EDGES="$2";            shift 2 ;;
        --w_diffcoexp)        W_DIFFCOEXP="$2";          shift 2 ;;
        --w_cluster)          W_CLUSTER="$2";            shift 2 ;;
        --w_pcre)             W_PCRE="$2";               shift 2 ;;
        --edge_classes)       EDGE_CLASSES="$2";         shift 2 ;;
        
        # === Pipeline 7c 参数 ===
        --rdata_7b)           RDATA_7B="$2";             shift 2 ;;
        --fig_format)         FIG_FORMAT="$2";           shift 2 ;;
        --fig_dpi)            FIG_DPI="$2";              shift 2 ;;
        --fig_scale)          FIG_SCALE="$2";            shift 2 ;;
        --plots)              PLOTS="$2";                shift 2 ;;
        --label_top_n)        LABEL_TOP_N="$2";          shift 2 ;;
        --label_focus)        LABEL_FOCUS="$2";          shift 2 ;;
        --show_heatmap_stars) SHOW_HEATMAP_STARS="$2";   shift 2 ;;
        --network_max_edges)  NETWORK_MAX_EDGES="$2";    shift 2 ;;
        --network_keep_deg1_all)   NETWORK_KEEP_DEG1_ALL="$2";   shift 2 ;;
        --network_keep_deg1_focus) NETWORK_KEEP_DEG1_FOCUS="$2"; shift 2 ;;
        --network_layout)     NETWORK_LAYOUT="$2";       shift 2 ;;
        --network_seed)       NETWORK_SEED="$2";         shift 2 ;;
        --network_prune)      NETWORK_PRUNE="$2";        shift 2 ;;
        --network_fig_w)      NETWORK_FIG_W="$2";        shift 2 ;;
        --network_fig_h)      NETWORK_FIG_H="$2";        shift 2 ;;
        --network_arrow_mm)   NETWORK_ARROW_MM="$2";     shift 2 ;;
        --network_edge_width_min) NETWORK_EDGE_WIDTH_MIN="$2"; shift 2 ;;
        --network_edge_width_max) NETWORK_EDGE_WIDTH_MAX="$2"; shift 2 ;;
        --network_node_size_min)  NETWORK_NODE_SIZE_MIN="$2";  shift 2 ;;
        --network_node_size_max)  NETWORK_NODE_SIZE_MAX="$2";  shift 2 ;;
        --network_ego_steps)  CLI_EGO_STEPS="$2";        shift 2 ;;
        --network_max_secondary_steps) NETWORK_MAX_SECONDARY_STEPS="$2"; shift 2 ;;
        --network_intra_completion) NETWORK_INTRA_COMPLETION="$2"; shift 2 ;;
        --network_shared_targets) NETWORK_SHARED_TARGETS="$2"; shift 2 ;;
        --gene_id_strip)          GENE_ID_STRIP_PATTERN="$2"; shift 2 ;;
        --circos_max_edges)   CIRCOS_MAX_EDGES="$2";     shift 2 ;;
        --circos_top_edges)   CIRCOS_TOP_EDGES="$2";     shift 2 ;;
        --deseq2_rdata)       DESEQ2_RDATA="$2";         shift 2 ;;
        --control_time)       CTRL_TIME="$2";            shift 2 ;;
        --cluster_colors)     CLUSTER_COLORS_FILE="$2";  shift 2 ;;
        --edge_class_colors)  EDGE_CLASS_COLORS_FILE="$2"; shift 2 ;;
        
        -h|--help)            usage ;;
        *)                    echo "Unknown option: $1"; usage ;;
    esac
done

# 如果未指定步骤，则默认全部执行
if [ "${STEP_SPECIFIED}" == "FALSE" ]; then
    RUN_7A="TRUE"
    RUN_7B="TRUE"
    RUN_7C="TRUE"
fi

# v4.5: -n 等价于 -c + network_only
if [ "${RUN_7C_NET}" == "TRUE" ]; then
    RUN_7C="TRUE"
    PLOTS="network"
fi

# ========================================================
# 7. 输入验证 (按需验证)
# ========================================================
ERRORS=0

# --- 7a 验证 ---
if [ "${RUN_7A}" == "TRUE" ]; then
    if [ ! -f "${RDATA_P4}" ]; then
        echo "❌ Error: Pipeline4 RData not found: ${RDATA_P4}"
        ERRORS=1
    fi
    if [ ! -f "${GRN_TF_LIST}" ]; then
        echo "❌ Error: TF list file not found: ${GRN_TF_LIST}"
        ERRORS=1
    fi
fi

# --- 7b 验证 ---
if [ "${RUN_7B}" == "TRUE" ] && [ "${RUN_7A}" == "FALSE" ]; then
    if [ ! -f "${RDATA_7A}" ]; then
        echo "❌ Error: Pipeline7a RData not found: ${RDATA_7A}"
        echo "   Please run Pipeline 7a first or provide --rdata_7a"
        ERRORS=1
    fi
fi

# --- 7c 验证 ---
if [ "${RUN_7C}" == "TRUE" ] && [ "${RUN_7B}" == "FALSE" ]; then
    if [ ! -f "${RDATA_7B}" ]; then
        echo "❌ Error: Pipeline7b RData not found: ${RDATA_7B}"
        echo "   Please run Pipeline 7b first or provide --rdata_7b"
        ERRORS=1
    fi
fi

# Pipeline 5 RData 检查 (7a 可选)
USE_P5="TRUE"
if [ "${RDATA_P5}" == "none" ] || [ ! -f "${RDATA_P5}" ]; then
    if [ "${RDATA_P5}" != "none" ] && [ ! -f "${RDATA_P5}" ]; then
        echo "⚠  Warning: Pipeline5 RData not found: ${RDATA_P5}"
        echo "   Mfuzz gene filtering will not be available."
        if [ "${TARGET_FILTER}" != "all" ]; then
            echo "   → target_filter will fall back to 'all'"
        fi
    fi
    RDATA_P5=""
    USE_P5="FALSE"
fi

# Cluster pairs 检查 (7b 可选)
if [ -n "${CLUSTER_PAIRS_FILE}" ] && [ ! -f "${CLUSTER_PAIRS_FILE}" ]; then
    echo "⚠  Warning: Cluster pairs file not found: ${CLUSTER_PAIRS_FILE}"
    CLUSTER_PAIRS_FILE=""
fi

# pCRE-TF map 检查 (7b 可选, 三级优先级)
USE_PCRE="FALSE"
PCRE_SOURCE="none"
if [ -n "${PCRE_TF_MAP}" ] && [ -f "${PCRE_TF_MAP}" ]; then
    USE_PCRE="TRUE"
    PCRE_SOURCE="manual"
elif [ -n "${PCRE_TF_MAP}" ] && [ ! -f "${PCRE_TF_MAP}" ]; then
    echo "⚠  Warning: pCRE-TF map not found: ${PCRE_TF_MAP}"
    PCRE_TF_MAP=""
fi

# 优先级 2: 6c 自动解析 (仅当无手动 map 时生效)
USE_PCRE_AUTO="FALSE"
if [ "${USE_PCRE}" == "FALSE" ] && [ -n "${PCRE_INTEGRATED_DIR}" ]; then
    if [ -d "${PCRE_INTEGRATED_DIR}" ]; then
        # 检查 WT 子目录下是否有整合表
        if [ -f "${PCRE_INTEGRATED_DIR}/WT/Integrated_Features_Full.txt" ] || \
           [ -f "${PCRE_INTEGRATED_DIR}/Integrated_Features_Full.txt" ]; then
            USE_PCRE_AUTO="TRUE"
            PCRE_SOURCE="auto_6c"
            echo "ℹ  pCRE auto-parse: will use Pipeline 6c from ${PCRE_INTEGRATED_DIR}"
        else
            echo "⚠  Warning: Pipeline 6c dir exists but no Integrated_Features_Full.txt found"
        fi
    else
        echo "⚠  Warning: Pipeline 6c dir not found: ${PCRE_INTEGRATED_DIR}"
        PCRE_INTEGRATED_DIR=""
    fi
fi

if [ ${ERRORS} -ne 0 ]; then
    echo ""
    echo "Fix the above errors and re-run."
    exit 1
fi

# ========================================================
# 8. 创建输出目录
# ========================================================
mkdir -p "${OUT_DIR}/Pipeline7a"
mkdir -p "${OUT_DIR}/Pipeline7b"
mkdir -p "${OUT_DIR}/Pipeline7b/agriGO_lists"
mkdir -p "${OUT_DIR}/Pipeline7c"
mkdir -p "${OUT_DIR}/rdata"

# ========================================================
# 9. 显示执行计划
# ========================================================
echo ""
echo "╔═══════════════════════════════════════════════════════╗"
echo "║     Pipeline 7: GRN Inference - Unified Workflow      ║"
echo "╚═══════════════════════════════════════════════════════╝"
echo ""
echo " Execution Plan:"
[ "${RUN_7A}" == "TRUE" ] && echo "   ✓ Step A: Differential Co-expression Analysis"
[ "${RUN_7B}" == "TRUE" ] && echo "   ✓ Step B: Cluster-Constrained Discovery + MECS"
[ "${RUN_7C}" == "TRUE" ] && [ "${RUN_7C_NET}" == "TRUE" ] && echo "   ✓ Step C: Network Visualization (Fig4 only)"
[ "${RUN_7C}" == "TRUE" ] && [ "${RUN_7C_NET}" != "TRUE" ] && echo "   ✓ Step C: Network Visualization"
echo ""
echo " Output Directory: ${OUT_DIR}"
echo ""

# ========================================================
# 10. 执行 Pipeline 7a
# ========================================================
if [ "${RUN_7A}" == "TRUE" ]; then
    echo "================================================="
    echo " [Step A] Differential Co-expression Analysis"
    echo "   R script: v1.3 (analytical + Mfuzz filtering)"
    echo "================================================="
    echo ""
    echo " Input:"
    echo "   P4 RData:   ${RDATA_P4}"
    [ "${USE_P5}" == "TRUE" ] && echo "   P5 RData:   ${RDATA_P5}" || echo "   P5 RData:   (not used)"
    echo "   TF List:    ${GRN_TF_LIST}"
    [ -n "${GRN_HIGHLIGHT}" ] && [ -f "${GRN_HIGHLIGHT}" ] && echo "   Highlight:  ${GRN_HIGHLIGHT}"
    echo ""
    echo " Configuration:"
    echo "   Conditions:     WT=${COND_WT}, MUT=${COND_MUT}"
    echo "   Lag: ${LAG} | Method: ${COR_METHOD}"
    echo "   Permutations:   ${N_PERM} | FDR: ${FDR}"
    echo "   Target filter:  ${TARGET_FILTER}"
    echo "   |r| pre-filter: ${PREFILTER_R}"
    echo "   Output Mode:    ${OUTPUT_MODE}"
    echo ""

    # 构建参数
    ARGS_7A=(
        "--rdata_p4"        "${RDATA_P4}"
        "--tf_list"         "${GRN_TF_LIST}"
        "--out_dir"         "${OUT_DIR}"
        "--condition_wt"    "${COND_WT}"
        "--condition_mut"   "${COND_MUT}"
        "--lag"             "${LAG}"
        "--n_permutations"  "${N_PERM}"
        "--fdr_cutoff"      "${FDR}"
        "--fdr_cutoff_diff" "${FDR_DIFF}"
        "--fdr_mut_loose"   "${FDR_MUT_LOOSE}"
        "--cor_method"      "${COR_METHOD}"
        "--output_mode"     "${OUTPUT_MODE}"
        "--n_cores"         "${N_CORES}"
        "--seed"            "${SEED}"
        "--target_filter"   "${TARGET_FILTER}"
        "--prefilter_r"     "${PREFILTER_R}"
    )

    [ "${USE_P5}" == "TRUE" ] && ARGS_7A+=("--rdata_p5" "${RDATA_P5}")
    [ -n "${GRN_HIGHLIGHT}" ] && [ -f "${GRN_HIGHLIGHT}" ] && ARGS_7A+=("--highlight_file" "${GRN_HIGHLIGHT}")
    [ -f "${METADATA_FILE}" ] && ARGS_7A+=("--samples_file" "${METADATA_FILE}")

    Rscript "${SCRIPT_DIR}/RNAseq_pipeline7a_diffcoexp.R" "${ARGS_7A[@]}"
    RET_7A=$?

    echo ""
    if [ ${RET_7A} -eq 0 ]; then
        echo "✅ Pipeline 7a: Success"
        SUMMARY_FILE="${OUT_DIR}/Pipeline7a/Pipeline7a_Summary.txt"
        if [ -f "${SUMMARY_FILE}" ]; then
            echo "--- Quick Summary ---"
            grep -E "^  (Significant|Lost|Gained|Rewired|Conserved|WT_only|MUT_only|NS)" "${SUMMARY_FILE}" 2>/dev/null
        fi
        # 更新 7b 输入路径
        RDATA_7A="${OUT_DIR}/rdata/Pipeline7a_DiffCoexp_Results.RData"
    else
        echo "❌ Pipeline 7a: Failed (exit code: ${RET_7A})"
        echo "   Check log: ${OUT_DIR}/Pipeline7a/Pipeline7a_diffcoexp.log"
        exit 1
    fi
    echo ""
fi

# ========================================================
# 11. 执行 Pipeline 7b
# ========================================================
if [ "${RUN_7B}" == "TRUE" ]; then
    echo "================================================="
    echo " [Step B] Cluster-Constrained Discovery + MECS"
    echo "================================================="
    echo ""
    echo " Input:"
    echo "   7a RData:      ${RDATA_7A}"
    [ -n "${CLUSTER_PAIRS_FILE}" ] && echo "   Cluster Pairs: ${CLUSTER_PAIRS_FILE}" || echo "   Cluster Pairs: (auto-infer from P5)"
    if [ "${USE_PCRE}" == "TRUE" ]; then
        echo "   pCRE-TF Map:   ${PCRE_TF_MAP} (manual)"
    elif [ "${USE_PCRE_AUTO}" == "TRUE" ]; then
        echo "   pCRE Source:    ${PCRE_INTEGRATED_DIR} (auto 6c)"
        [ -n "${PCRE_FAMILY_ALIAS}" ] && echo "   pCRE Alias:     ${PCRE_FAMILY_ALIAS}"
    else
        echo "   pCRE-TF Map:   (not used, MECS = S1+S2)"
    fi
    echo ""
    echo " Focus Configuration:"
    echo "   TF Family:  ${FOCUS_TF_FAMILY}"
    echo "   Gene:       ${FOCUS_GENE}"
    echo ""
    echo " Parameters:"
    echo "   Edge Classes:      ${EDGE_CLASSES}"
    echo "   Enrichment:        ${ENRICHMENT_METHOD} (min_edges=${MIN_EDGES})"
    echo "   Background:        ${BACKGROUND_MODE}"
    echo "   MECS Weights:      S1=${W_DIFFCOEXP}, S2=${W_CLUSTER}, S3=${W_PCRE}"
    echo ""

    # 构建参数
    ARGS_7B=(
        "--rdata_7a"          "${RDATA_7A}"
        "--out_dir"           "${OUT_DIR}"
        "--focus_tf_family"   "${FOCUS_TF_FAMILY}"
        "--focus_gene"        "${FOCUS_GENE}"
        "--background_mode"   "${BACKGROUND_MODE}"
        "--enrichment_method" "${ENRICHMENT_METHOD}"
        "--min_edges"         "${MIN_EDGES}"
        "--w_diffcoexp"       "${W_DIFFCOEXP}"
        "--w_cluster"         "${W_CLUSTER}"
        "--w_pcre"            "${W_PCRE}"
        "--edge_classes"      "${EDGE_CLASSES}"
    )

    [ -n "${CLUSTER_PAIRS_FILE}" ] && [ -f "${CLUSTER_PAIRS_FILE}" ] && ARGS_7B+=("--cluster_pairs" "${CLUSTER_PAIRS_FILE}")
    [ -n "${PCRE_TF_MAP}" ] && [ -f "${PCRE_TF_MAP}" ] && ARGS_7B+=("--pcre_tf_map" "${PCRE_TF_MAP}")
    [ "${USE_PCRE_AUTO}" == "TRUE" ] && [ -n "${PCRE_INTEGRATED_DIR}" ] && ARGS_7B+=("--pcre_integrated_dir" "${PCRE_INTEGRATED_DIR}")
    [ -n "${PCRE_FAMILY_ALIAS}" ] && ARGS_7B+=("--pcre_family_alias" "${PCRE_FAMILY_ALIAS}")
    [ "${PCRE_MIN_PCC}" != "0" ] && ARGS_7B+=("--pcre_min_pcc" "${PCRE_MIN_PCC}")
    [ -n "${HIGHLIGHT_FILE}" ] && [ -f "${HIGHLIGHT_FILE}" ] && ARGS_7B+=("--highlight_file" "${HIGHLIGHT_FILE}")

    Rscript "${SCRIPT_DIR}/RNAseq_pipeline7b_discovery.R" "${ARGS_7B[@]}"
    RET_7B=$?

    echo ""
    if [ ${RET_7B} -eq 0 ]; then
        echo "✅ Pipeline 7b: Success"
        MECS_SUMMARY="${OUT_DIR}/Pipeline7b/Pipeline7b_MECS_Summary.txt"
        if [ -f "${MECS_SUMMARY}" ]; then
            echo "--- MECS Tier Summary ---"
            grep -E "^  Tier" "${MECS_SUMMARY}" 2>/dev/null
        fi
        # 更新 7c 输入路径
        RDATA_7B="${OUT_DIR}/rdata/Pipeline7b_Discovery_Results.RData"
    else
        echo "❌ Pipeline 7b: Failed (exit code: ${RET_7B})"
        echo "   Check log: ${OUT_DIR}/Pipeline7b/Pipeline7b_discovery.log"
        exit 1
    fi
    echo ""
fi

# ========================================================
# 12. 执行 Pipeline 7c
# ========================================================
if [ "${RUN_7C}" == "TRUE" ]; then
    echo "================================================="
    echo " [Step C] Network Visualization (v4.6)"
    echo "================================================="
    echo ""

    # 检查 7c 依赖
    if [ ! -f "${RDATA_7B}" ]; then
        echo "❌ Error: Pipeline7b RData not found: ${RDATA_7B}"
        exit 1
    fi

    # DESeq2 RData 检查 (evidence_heatmap 依赖)
    if echo "${PLOTS}" | grep -q "evidence_heatmap"; then
        if [ ! -f "${DESEQ2_RDATA}" ]; then
            echo "⚠ Warning: DESeq2 RData not found: ${DESEQ2_RDATA}"
        fi
    fi

    # === v4.5: Highlight 来源解析 ===
    HAS_USER_HIGHLIGHT="FALSE"
    HAS_AUTO_HIGHLIGHT="FALSE"
    AUTO_HL_OUTPUT="${OUT_DIR}/Pipeline7c/auto_highlight.txt"
    EFFECTIVE_HIGHLIGHT=""

    # 检查用户提供的 highlight
    if [ -n "${GRN_HIGHLIGHT}" ] && [ -f "${GRN_HIGHLIGHT}" ]; then
        HAS_USER_HIGHLIGHT="TRUE"
        echo " User highlight file: ${GRN_HIGHLIGHT}"
    fi

    # 检查是否需要 cluster auto-highlight (独立于用户 highlight)
    if [ -n "${AUTO_HL_CLUSTERS}" ] && [ "${AUTO_HL_CLUSTERS}" != "none" ]; then
        echo " Auto-highlight clusters: ${AUTO_HL_CLUSTERS}"
        
        # --- 运行 7d 生成 auto-highlight ---
        echo ""
        echo " --- Auto-generating highlight from clusters ---"
        
        WT_MEM="${BASE_DIR}/05_Mfuzz/morphology_analysis/WT_Membership_Matrix.txt"
        MUT_MEM="${BASE_DIR}/05_Mfuzz/morphology_analysis/MUT_Membership_Matrix.txt"
        
        AUTO_HL_ARGS=(
            "--clusters"              "${AUTO_HL_CLUSTERS}"
            "--membership_threshold"  "${MFUZZ_MEMBERSHIP_CUTOFF:-0.5}"
            "--max_genes_per_cluster" "${GRN_AUTO_HIGHLIGHT_MAX_PER_CLUSTER:-100}"
            "--tier_filter"           "1,2"
            "--min_tier"              "${AUTO_HL_MIN_TIER}"
            "--min_edge_count"        "${AUTO_HL_MIN_EDGE_COUNT}"
            "--output"                "${AUTO_HL_OUTPUT}"
        )
        
        [ -f "${WT_MEM}" ]   && AUTO_HL_ARGS+=("--wt_membership" "${WT_MEM}")
        [ -f "${MUT_MEM}" ]  && AUTO_HL_ARGS+=("--mut_membership" "${MUT_MEM}")
        [ -f "${RDATA_7B}" ] && AUTO_HL_ARGS+=("--rdata_7b" "${RDATA_7B}")
        [ -n "${GFF_FILE}" ] && [ -f "${GFF_FILE}" ] && AUTO_HL_ARGS+=("--gff_file" "${GFF_FILE}")
        [ -n "${GENE_ID_STRIP_PATTERN}" ] && AUTO_HL_ARGS+=("--gene_id_strip" "${GENE_ID_STRIP_PATTERN}")
        
        Rscript "${SCRIPT_DIR}/RNAseq_pipeline7d_generate_highlight_from_clusters.R" \
            "${AUTO_HL_ARGS[@]}"
        
        if [ $? -eq 0 ] && [ -f "${AUTO_HL_OUTPUT}" ]; then
            N_GENES=$(wc -l < "${AUTO_HL_OUTPUT}")
            echo " ✔ Auto-highlight generated: ${N_GENES} genes"
            HAS_AUTO_HIGHLIGHT="TRUE"
        else
            echo " ⚠ Auto-highlight generation failed."
        fi
    fi

    # === v4.5: 合并逻辑 ===
    if [ "${HAS_USER_HIGHLIGHT}" == "TRUE" ] && [ "${HAS_AUTO_HIGHLIGHT}" == "TRUE" ]; then
        # 两者共存: 合并并标记来源
        echo " Merging user highlight (primary) + cluster highlight (secondary)..."
        MERGED="${OUT_DIR}/Pipeline7c/merged_highlight.txt"
        
        # v4.6 fix: 跳过用户文件表头行 (NR>1), 避免 "Symbol\tGene\tprimary" 污染合并文件
        # 先写入标准表头
        echo -e "Symbol\tGene\tSource" > "${MERGED}"
        awk -F'\t' 'NR>1{print $1"\t"$2"\tprimary"}' "${GRN_HIGHLIGHT}" >> "${MERGED}"
        awk -F'\t' 'NR>1{print $2}' "${GRN_HIGHLIGHT}" | sort -u > "${MERGED}.user_ids"
        grep -v -F -f "${MERGED}.user_ids" <(awk -F'\t' '{print $1"\t"$2"\tsecondary"}' "${AUTO_HL_OUTPUT}") >> "${MERGED}" 2>/dev/null || true
        rm -f "${MERGED}.user_ids"
        
        N_PRI=$(grep -c "primary" "${MERGED}" || true)
        N_SEC=$(grep -c "secondary" "${MERGED}" || true)
        echo "   Primary: ${N_PRI} | Secondary: ${N_SEC} | Total: $((N_PRI + N_SEC))"
        
        EFFECTIVE_HIGHLIGHT="${MERGED}"

    elif [ "${HAS_USER_HIGHLIGHT}" == "TRUE" ]; then
        # 仅用户列表: 全部标记为 primary (添加第三列)
        TAGGED="${OUT_DIR}/Pipeline7c/tagged_highlight.txt"
        echo -e "Symbol\tGene\tSource" > "${TAGGED}"
        awk -F'\t' 'NR>1{print $1"\t"$2"\tprimary"}' "${GRN_HIGHLIGHT}" >> "${TAGGED}"
        EFFECTIVE_HIGHLIGHT="${TAGGED}"

    elif [ "${HAS_AUTO_HIGHLIGHT}" == "TRUE" ]; then
        # 仅 auto: 全部标记为 primary (无用户列表时 cluster 基因享有完整锚点地位)
        TAGGED="${OUT_DIR}/Pipeline7c/tagged_highlight.txt"
        echo -e "Symbol\tGene\tSource" > "${TAGGED}"
        awk -F'\t' '{print $1"\t"$2"\tprimary"}' "${AUTO_HL_OUTPUT}" >> "${TAGGED}"
        EFFECTIVE_HIGHLIGHT="${TAGGED}"

    else
        echo " ⚠ No highlight genes from any source."
    fi

    # === v4.5: Ego-steps 智能分配 ===
    if [ -n "${CLI_EGO_STEPS}" ]; then
        NETWORK_EGO_STEPS="${CLI_EGO_STEPS}"
    elif [ "${HAS_USER_HIGHLIGHT}" == "TRUE" ]; then
        # 有用户列表 (无论是否合并了 cluster): ego=0, 不做距离截断
        NETWORK_EGO_STEPS="0"
    else
        # 纯 auto 模式: 使用 config 值
        NETWORK_EGO_STEPS="${GRN_NETWORK_EGO_STEPS:-2}"
    fi
    echo " ★ Ego-steps: ${NETWORK_EGO_STEPS}"

    echo ""
    echo " Input:    ${RDATA_7B}"
    echo " Plots:    ${PLOTS}"
    echo " Format:   ${FIG_FORMAT} (DPI: ${FIG_DPI})"
    echo ""

    # 构建参数
    ARGS_7C=(
        "--rdata_7b"                "${RDATA_7B}"
        "--out_dir"                 "${OUT_DIR}"
        "--fig_format"              "${FIG_FORMAT}"
        "--fig_dpi"                 "${FIG_DPI}"
        "--fig_scale"               "${FIG_SCALE}"
        "--plots"                   "${PLOTS}"
        "--focus_tf_family"         "${FOCUS_TF_FAMILY}"
        "--focus_gene"              "${FOCUS_GENE}"
        "--condition_wt"            "${COND_WT}"
        "--condition_mut"           "${COND_MUT}"
        "--label_top_n"             "${LABEL_TOP_N}"
        "--label_focus"             "${LABEL_FOCUS}"
        "--show_heatmap_stars"      "${SHOW_HEATMAP_STARS}"
        "--network_max_edges"       "${NETWORK_MAX_EDGES}"
        "--network_keep_deg1_all"   "${NETWORK_KEEP_DEG1_ALL}"
        "--network_keep_deg1_focus" "${NETWORK_KEEP_DEG1_FOCUS}"
        "--network_layout"          "${NETWORK_LAYOUT}"
        "--network_seed"            "${NETWORK_SEED}"
        "--network_prune"           "${NETWORK_PRUNE}"
        "--network_fig_w"           "${NETWORK_FIG_W}"
        "--network_fig_h"           "${NETWORK_FIG_H}"
        "--network_arrow_mm"        "${NETWORK_ARROW_MM}"
        "--network_edge_width_min"  "${NETWORK_EDGE_WIDTH_MIN}"
        "--network_edge_width_max"  "${NETWORK_EDGE_WIDTH_MAX}"
        "--network_node_size_min"   "${NETWORK_NODE_SIZE_MIN}"
        "--network_node_size_max"   "${NETWORK_NODE_SIZE_MAX}"
        "--network_ego_steps"       "${NETWORK_EGO_STEPS}"
        "--network_max_secondary_steps" "${NETWORK_MAX_SECONDARY_STEPS}"
        "--network_intra_completion" "${NETWORK_INTRA_COMPLETION}"
        "--network_shared_targets" "${NETWORK_SHARED_TARGETS}"
        "--circos_max_edges"        "${CIRCOS_MAX_EDGES}"
        "--circos_top_edges"        "${CIRCOS_TOP_EDGES}"
        "--deseq2_rdata"            "${DESEQ2_RDATA}"
        "--control_time"            "${CTRL_TIME}"
    )

    # ★ v4.6: 传递基因ID清理正则给 7c (统一与 7d 相同的正则化逻辑)
    [ -n "${GENE_ID_STRIP_PATTERN}" ] && ARGS_7C+=("--gene_id_strip" "${GENE_ID_STRIP_PATTERN}")

    [ -n "${CLUSTER_COLORS_FILE}" ] && [ -f "${CLUSTER_COLORS_FILE}" ] && ARGS_7C+=("--cluster_colors" "${CLUSTER_COLORS_FILE}")
    [ -n "${EDGE_CLASS_COLORS_FILE}" ] && [ -f "${EDGE_CLASS_COLORS_FILE}" ] && ARGS_7C+=("--edge_class_colors" "${EDGE_CLASS_COLORS_FILE}")
    [ -n "${EFFECTIVE_HIGHLIGHT}" ] && [ -f "${EFFECTIVE_HIGHLIGHT}" ] && ARGS_7C+=("--highlight_file" "${EFFECTIVE_HIGHLIGHT}")

    Rscript "${SCRIPT_DIR}/RNAseq_pipeline7c_visualization.R" "${ARGS_7C[@]}"
    RET_7C=$?

    echo ""
    if [ ${RET_7C} -eq 0 ]; then
        echo "✅ Pipeline 7c: Success"
    else
        echo "❌ Pipeline 7c: Failed (exit code: ${RET_7C})"
        exit 1
    fi
fi

# ========================================================
# 13. 最终汇总
# ========================================================
echo ""
echo "╔═══════════════════════════════════════════════════════╗"
echo "║            Pipeline 7 Execution Complete              ║"
echo "╚═══════════════════════════════════════════════════════╝"
echo ""
echo " Output Directory: ${OUT_DIR}"
echo ""
echo " Generated Files:"
[ "${RUN_7A}" == "TRUE" ] && echo "   ├── Pipeline7a/  (Differential co-expression results)"
[ "${RUN_7B}" == "TRUE" ] && echo "   ├── Pipeline7b/  (Discovery + MECS scoring)"
[ "${RUN_7C}" == "TRUE" ] && echo "   ├── Pipeline7c/  (Visualizations)"
echo "   └── rdata/       (RData objects for downstream analysis)"
echo ""
echo "✅ All requested steps completed successfully!"