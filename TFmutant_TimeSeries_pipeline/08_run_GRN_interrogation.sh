#!/bin/bash
# =============================================================================
# script_name: 08_run_GRN_interrogation.sh
# 功能: GRN 锚点网络深度查询 — 衔接 07_run_GRN.sh 流程
# 版本: v1.1
#
# 用法示例:
#   # 单基因查询
#   ./08_run_GRN_interrogation.sh --mode gene_profile \
#       --query_genes Cre16.g659700
#
#   # 带批次标签 (输出到 08_GRN_Interrogation/CrDof_PSR1/)
#   ./08_run_GRN_interrogation.sh --mode shared_targets \
#       --query_genes Cre12.g521150,Cre12.g495100 \
#       --batch CrDof_PSR1
#
#   # 拓扑表 + MADS重叠度
#   ./08_run_GRN_interrogation.sh --mode topology_table \
#       --degree_cutoff 15 --ref_tf Cre11.g467577
#
#   # 全量报告
#   ./08_run_GRN_interrogation.sh --mode full_report \
#       --query_genes Cre12.g521150,Cre12.g495100,Cre11.g467577 \
#       --batch full_20260320
# =============================================================================

# ========================================================
# 1. 定位脚本 & 加载配置
# ========================================================
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

if [ -f "${SCRIPT_DIR}/project_config.sh" ]; then
    source "${SCRIPT_DIR}/project_config.sh"
else
    echo "⚠  Warning: project_config.sh not found in ${SCRIPT_DIR}"
fi

# ========================================================
# 2. 默认参数
# ========================================================
INPUT_DIR="${GRN_P7C_DIR:-${BASE_DIR:-}/07_GRN/Pipeline7c}"
OUT_DIR="${BASE_DIR:-$(pwd)}/08_GRN_Interrogation"

MODE="gene_profile"
QUERY_GENES=""
QUERY_FILE=""
DEGREE_CUTOFF=15
REF_TF=""
OUT_PREFIX="P8"
BATCH=""
COND_WT="${CONDITION_WT:-WT}"
COND_MUT="${CONDITION_MUT:-dof}"

EDGES_CSV=""
NODES_CSV=""
TARGET_LFC=""
TF_LFC=""

# ========================================================
# 3. 帮助
# ========================================================
usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Pipeline 8: GRN Interrogation (v1.1)"
    echo ""
    echo "Modes (--mode):"
    echo "  gene_profile            基因中心拓扑+表达全景"
    echo "  topology_table          高Degree节点边类型分布表"
    echo "  shared_targets          两个TF共享靶基因"
    echo "  overlap_matrix          多TF靶基因重叠矩阵"
    echo "  expression_divergence   聚类×边类型 WT-dof分歧"
    echo "  subspace                调控子空间鉴定"
    echo "  full_report             全部模块"
    echo ""
    echo "Options:"
    echo "  --mode MODE             分析模式 (default: gene_profile)"
    echo "  --query_genes IDs       查询基因, 逗号分隔"
    echo "  --query_file PATH       查询基因列表文件"
    echo "  --input_dir PATH        Pipeline7c 输出目录"
    echo "  --out_dir PATH          输出根目录"
    echo "  --batch NAME            批次/项目标签 → 子目录"
    echo "  --out_prefix STR        文件前缀 (default: P8)"
    echo "  --degree_cutoff N       Degree阈值 (default: 15)"
    echo "  --ref_tf ID             重叠度参考TF"
    echo "  --condition_wt NAME     WT名 (default: WT)"
    echo "  --condition_mut NAME    突变体名 (default: dof)"
    echo "  --edges_csv PATH        覆盖Edges CSV"
    echo "  --nodes_csv PATH        覆盖Nodes CSV"
    echo "  --target_lfc PATH       覆盖Target LFC CSV"
    echo "  --tf_lfc PATH           覆盖TF LFC CSV"
    echo "  -h, --help"
    exit 0
}

# ========================================================
# 4. 参数解析
# ========================================================
while [[ $# -gt 0 ]]; do
    case $1 in
        --mode)           MODE="$2";           shift 2 ;;
        --query_genes)    QUERY_GENES="$2";    shift 2 ;;
        --query_file)     QUERY_FILE="$2";     shift 2 ;;
        --input_dir)      INPUT_DIR="$2";      shift 2 ;;
        --out_dir)        OUT_DIR="$2";        shift 2 ;;
        --batch)          BATCH="$2";          shift 2 ;;
        --out_prefix)     OUT_PREFIX="$2";     shift 2 ;;
        --degree_cutoff)  DEGREE_CUTOFF="$2";  shift 2 ;;
        --ref_tf)         REF_TF="$2";         shift 2 ;;
        --condition_wt)   COND_WT="$2";        shift 2 ;;
        --condition_mut)  COND_MUT="$2";       shift 2 ;;
        --edges_csv)      EDGES_CSV="$2";      shift 2 ;;
        --nodes_csv)      NODES_CSV="$2";      shift 2 ;;
        --target_lfc)     TARGET_LFC="$2";     shift 2 ;;
        --tf_lfc)         TF_LFC="$2";         shift 2 ;;
        -h|--help)        usage ;;
        *)                echo "Unknown: $1"; usage ;;
    esac
done

# ========================================================
# 5. 显示配置
# ========================================================
echo ""
echo "╔══════════════════════════════════════════════════════════╗"
echo "║  Pipeline 8: GRN Interrogation (v1.1)                   ║"
echo "╚══════════════════════════════════════════════════════════╝"
echo ""
echo "  Input dir: ${INPUT_DIR}"
echo "  Mode:      ${MODE}"
[ -n "${BATCH}" ] && echo "  Batch:     ${BATCH}"
echo "  Output:    ${OUT_DIR}"
echo ""

# ========================================================
# 6. 构建 R 参数 & 执行
# ========================================================
mkdir -p "${OUT_DIR}"

R_ARGS=(
    "--mode"          "${MODE}"
    "--out_dir"       "${OUT_DIR}"
    "--out_prefix"    "${OUT_PREFIX}"
    "--degree_cutoff" "${DEGREE_CUTOFF}"
    "--condition_wt"  "${COND_WT}"
    "--condition_mut" "${COND_MUT}"
)

[ -n "${BATCH}" ]      && R_ARGS+=("--batch" "${BATCH}")
[ -n "${EDGES_CSV}" ]  && R_ARGS+=("--edges_csv" "${EDGES_CSV}")
[ -n "${NODES_CSV}" ]  && R_ARGS+=("--nodes_csv" "${NODES_CSV}")
[ -n "${TARGET_LFC}" ] && R_ARGS+=("--target_lfc" "${TARGET_LFC}")
[ -n "${TF_LFC}" ]     && R_ARGS+=("--tf_lfc" "${TF_LFC}")

if [ -z "${EDGES_CSV}" ] && [ -z "${NODES_CSV}" ] && \
   [ -z "${TARGET_LFC}" ] && [ -z "${TF_LFC}" ]; then
    R_ARGS+=("--input_dir" "${INPUT_DIR}")
fi

[ -n "${QUERY_GENES}" ] && R_ARGS+=("--query_genes" "${QUERY_GENES}")
[ -n "${QUERY_FILE}" ]  && R_ARGS+=("--query_file" "${QUERY_FILE}")
[ -n "${REF_TF}" ]      && R_ARGS+=("--ref_tf" "${REF_TF}")

Rscript "${SCRIPT_DIR}/RNAseq_pipeline8_grn_interrogation.R" "${R_ARGS[@]}"
RET=$?

# ========================================================
# 7. 结果汇总
# ========================================================
echo ""
if [ ${RET} -eq 0 ]; then
    # 确定实际输出目录
    if [ -n "${BATCH}" ]; then
        ACTUAL_OUT="${OUT_DIR}/${BATCH}"
    else
        ACTUAL_OUT="${OUT_DIR}/Pipeline8"
    fi

    echo "╔══════════════════════════════════════════════════════════╗"
    echo "║  ✅ Pipeline 8 Complete                                  ║"
    echo "╚══════════════════════════════════════════════════════════╝"
    echo ""
    echo "  Output: ${ACTUAL_OUT}/"
    echo ""

    if [ -d "${ACTUAL_OUT}" ]; then
        echo "  Files:"
        # 人类可读的文件大小
        ls -lhS "${ACTUAL_OUT}/" 2>/dev/null | grep "^-" | \
            awk '{printf "    %-45s %s\n", $NF, $5}'
    fi
else
    echo "❌ Pipeline 8 Failed (exit code: ${RET})"
    exit 1
fi