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

NETWORK_MAX_NODES="${GRN_NETWORK_MAX_NODES:-200}"
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
AUTO_HL_HUB_N="${GRN_AUTO_HIGHLIGHT_HUB_N:-30}"

CLI_EGO_STEPS=""

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
    echo "  --pcre_tf_map PATH    pCRE-TF family mapping TSV"
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
    echo "  --network_max_nodes N Max nodes in network (default: ${NETWORK_MAX_NODES})"
    echo "  --network_layout L    Layout algorithm: kk/fr/circle (default: ${NETWORK_LAYOUT})"
    echo "  --network_seed N      Network layout seed (default: ${NETWORK_SEED})"
    echo "  --network_prune BOOL  Prune network (default: ${NETWORK_PRUNE})"
    echo "  --network_fig_w N     Network figure width (default: ${NETWORK_FIG_W})"
    echo "  --network_fig_h N     Network figure height (default: ${NETWORK_FIG_H})"
    echo "  --network_arrow_mm N  Arrow size in mm (default: ${NETWORK_ARROW_MM})"
    echo "  --network_edge_width_min N  Min edge width (default: ${NETWORK_EDGE_WIDTH_MIN})"
    echo "  --network_edge_width_max N  Max edge width (default: ${NETWORK_EDGE_WIDTH_MAX})"
    echo "  --network_node_size_min N   Min node size (default: ${NETWORK_NODE_SIZE_MIN})"
    echo "  --network_node_size_max N   Max node size (default: ${NETWORK_NODE_SIZE_MAX})"
    echo "  --network_ego_steps N Ego-network topology steps (default: auto)"
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
        --network_max_nodes)  NETWORK_MAX_NODES="$2";    shift 2 ;;
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

# pCRE-TF map 检查 (7b 可选)
USE_PCRE="FALSE"
if [ -n "${PCRE_TF_MAP}" ] && [ -f "${PCRE_TF_MAP}" ]; then
    USE_PCRE="TRUE"
elif [ -n "${PCRE_TF_MAP}" ] && [ ! -f "${PCRE_TF_MAP}" ]; then
    echo "⚠  Warning: pCRE-TF map not found: ${PCRE_TF_MAP}"
    PCRE_TF_MAP=""
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
[ "${RUN_7C}" == "TRUE" ] && echo "   ✓ Step C: Network Visualization"
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
    [ "${USE_PCRE}" == "TRUE" ] && echo "   pCRE-TF Map:   ${PCRE_TF_MAP}" || echo "   pCRE-TF Map:   (not used, MECS = S1+S2)"
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
    echo " [Step C] Network Visualization (v4.2.1)"
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

    # === Auto-generate highlight if needed ===
    AUTO_HL_OUTPUT="${OUT_DIR}/Pipeline7c/auto_highlight.txt"
    NEED_AUTO_HIGHLIGHT="FALSE"
    EFFECTIVE_HIGHLIGHT=""

    if [ -n "${GRN_HIGHLIGHT}" ] && [ -f "${GRN_HIGHLIGHT}" ]; then
        echo " Highlight file: ${GRN_HIGHLIGHT} (user-provided)"
        EFFECTIVE_HIGHLIGHT="${GRN_HIGHLIGHT}"
    elif [ -n "${AUTO_HL_CLUSTERS}" ] && [ "${AUTO_HL_CLUSTERS}" != "none" ]; then
        NEED_AUTO_HIGHLIGHT="TRUE"
        echo " No highlight file → auto-generating from clusters: ${AUTO_HL_CLUSTERS}"
    else
        echo " ⚠ No highlight file and no clusters specified."
    fi

    # === Ego-steps 智能分配 ===
    if [ -n "${CLI_EGO_STEPS}" ]; then
        NETWORK_EGO_STEPS="${CLI_EGO_STEPS}"
        echo " ★ Ego-steps explicitly set via CLI: ${NETWORK_EGO_STEPS}"
    elif [ "${NEED_AUTO_HIGHLIGHT}" == "TRUE" ]; then
        NETWORK_EGO_STEPS="${GRN_NETWORK_EGO_STEPS:-2}"
        echo " ★ Auto-mode: using Config/Default NETWORK_EGO_STEPS=${NETWORK_EGO_STEPS}"
    else
        NETWORK_EGO_STEPS="0"
        echo " ★ User-list mode: forcing NETWORK_EGO_STEPS=0"
    fi

    # === Auto-highlight 生成 ===
    if [ "${NEED_AUTO_HIGHLIGHT}" == "TRUE" ]; then
        echo ""
        echo " --- Auto-generating highlight from clusters ---"
        
        WT_MEM="${BASE_DIR}/05_Mfuzz/morphology_analysis/WT_Membership_Matrix.txt"
        MUT_MEM="${BASE_DIR}/05_Mfuzz/morphology_analysis/MUT_Membership_Matrix.txt"
        
        AUTO_HL_ARGS=(
            "--clusters"              "${AUTO_HL_CLUSTERS}"
            "--membership_threshold"  "${MFUZZ_MEMBERSHIP_CUTOFF:-0.5}"
            "--max_genes_per_cluster" "${GRN_AUTO_HIGHLIGHT_MAX_PER_CLUSTER:-100}"
            "--hub_top_n"             "${AUTO_HL_HUB_N}"
            "--tier_filter"           "1,2"
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
            EFFECTIVE_HIGHLIGHT="${AUTO_HL_OUTPUT}"
        else
            echo " ⚠ Auto-highlight generation failed, proceeding without."
        fi
    fi

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
        "--network_max_nodes"       "${NETWORK_MAX_NODES}"
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
        "--circos_max_edges"        "${CIRCOS_MAX_EDGES}"
        "--circos_top_edges"        "${CIRCOS_TOP_EDGES}"
        "--deseq2_rdata"            "${DESEQ2_RDATA}"
        "--control_time"            "${CTRL_TIME}"
    )

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