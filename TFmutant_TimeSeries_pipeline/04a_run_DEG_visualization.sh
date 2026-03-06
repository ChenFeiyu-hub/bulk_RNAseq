#!/bin/bash
# script_name: 05_run_DEG_visualization.sh
# 功能: 执行可视化分析 (P4b/4c)，支持绝对路径调用

set -e

# ========================================================
# 1. 路径解析 (支持 --script-dir 参数)
# ========================================================

# 默认脚本目录为当前脚本所在目录
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# 解析参数，允许外部覆盖 SCRIPT_DIR
for arg in "$@"; do
  case $arg in
    --script-dir=*)
      SCRIPT_DIR="${arg#*=}"
      shift
      ;;
    --script-dir)
      SCRIPT_DIR="$2"
      shift 2
      ;;
  esac
done

# 加载配置 (使用绝对路径)
if [ -f "${SCRIPT_DIR}/project_config.sh" ]; then
    source "${SCRIPT_DIR}/project_config.sh"
else
    echo "❌ Error: project_config.sh not found in ${SCRIPT_DIR}"
    exit 1
fi

# ========================================================
# 2. 定义输入输出 (核心修复)
# ========================================================

# [修复] 显式定义 DESEQ2_DIR
DESEQ2_DIR="${BASE_DIR}/04_DESeq2"

# 定义输入文件
PCA_INPUT="${DESEQ2_DIR}/RData/Pipeline4_Logic1_TimeCourse_for_Mfuzz.RData"
UPSET_INPUT="${DESEQ2_DIR}/RData/Pipeline4_Logic1_TimeCourse_for_Mfuzz.RData"

OUT_DIR="${DESEQ2_DIR}/Visualization"
LOG_DIR="${BASE_DIR}/logs"

# 检查输入
if [ ! -f "$PCA_INPUT" ]; then
    echo "❌ [Error] Input RData not found: $PCA_INPUT"
    exit 1
fi

mkdir -p "${OUT_DIR}/PCA"
mkdir -p "${OUT_DIR}/UpSet"
mkdir -p "${OUT_DIR}/Venn"
mkdir -p "${LOG_DIR}"

echo "========================================================"
echo "   Pipeline 4b/4c: Visualization"
echo "   Script Directory: ${SCRIPT_DIR}"
echo "   Data Directory:   ${DESEQ2_DIR}"
echo "========================================================"

# ========================================================
# 3. 执行 PCA Analysis (Pipeline 4b)
# ========================================================
echo -e "\n>>> [Module 4b] Running PCA Analysis..."

# [修复] 使用 ${SCRIPT_DIR} 调用 R 脚本
Rscript "${SCRIPT_DIR}/RNAseq_pipeline4b_PCA.R" \
    --input "${PCA_INPUT}" \
    --metadata "${METADATA_FILE}" \
    --outdir "${OUT_DIR}/PCA" \
    --prefix "${PROJECT_NAME}" \
    --control_cond "${CONTROL_CONDITION}" \
    --cond_colors "${STYLE_COND_COLORS}" \
    --cond_shapes "${STYLE_COND_SHAPES}" \
    --time_colors "${STYLE_TIME_COLORS}" \
    --label_config "${STYLE_LABELS}" \
    --ellipse_conf 0.95 \
    > "${LOG_DIR}/04b_PCA.log" 2>&1

if [ $? -eq 0 ]; then
    echo "✅ PCA Analysis Done."
else
    echo "❌ PCA Analysis Failed. Check log: ${LOG_DIR}/04b_PCA.log"
fi

# ========================================================
# 4. 执行 UpSet Analysis (Pipeline 4c)
# ========================================================
# 默认开启，除非显式关闭
if [[ "${RUN_UPSET:-true}" == true ]]; then
    echo -e "\n>>> [4c] Running UpSet Analysis (Standardized)..."
    
    # [修复] 使用 ${SCRIPT_DIR} 调用 R 脚本
    Rscript "${SCRIPT_DIR}/RNAseq_pipeline4c_UpSet.R" \
        --input "${UPSET_INPUT}" \
        --outdir "${OUT_DIR}/UpSet" \
        --prefix "${PROJECT_NAME}" \
        --condition_order "${CONDITION_ORDER}" \
        --time_order "${TIME_ORDER}" \
        --cond_colors "${STYLE_COND_COLORS}" \
        --fdr "${FDR_CUTOFF}" \
        --lfc "${LFC_THRESHOLD}" \
        > "${LOG_DIR}/04c_upset.log" 2>&1
    
    if [ $? -eq 0 ]; then
        echo "✅ UpSet Analysis Done."
    else
        echo "❌ UpSet Analysis Failed. Check log: ${LOG_DIR}/04c_upset.log"
        # 打印最后几行错误信息以便调试
        tail -n 3 "${LOG_DIR}/04c_upset.log"
        exit 1
    fi
fi

# ========================================================
# 5. 执行 Venn Panel Analysis (Pipeline 4d - New Design)
# ========================================================
if [[ "${RUN_VENN:-true}" == true ]]; then
    echo -e "\n>>> [4d] Running Venn Panel Analysis (Structured)..."
    
    # 使用包含所有 Logic 的 Master RData (或者 Logic 1 RData 也可以，只要里面有 res_list_time)
    # 注意: 4d 脚本会自动在 RData 中搜寻需要的比较对象
    INPUT_VENN="${DESEQ2_DIR}/RData/Pipeline4_Logic1_TimeCourse_for_Mfuzz.RData"
    
    Rscript "${SCRIPT_DIR}/RNAseq_pipeline4d_Venn.R" \
        --input "$INPUT_VENN" \
        --outdir "${OUT_DIR}/Venn" \
        --prefix "$PROJECT_NAME" \
        --control_cond "${CONTROL_CONDITION}" \
        --condition_order "${CONDITION_ORDER}" \
        --control_time "${CONTROL_TIME}" \
        --time_order "${TIME_ORDER}" \
        --cond_colors "${STYLE_COND_COLORS}" \
        --time_colors "${STYLE_TIME_COLORS}" \
        --fdr "${FDR_CUTOFF}" \
        --lfc "${LFC_THRESHOLD}" \
        > "${LOG_DIR}/04d_venn.log" 2>&1
    
    if [ $? -eq 0 ]; then
        echo "✅ Venn Panel Done."
    else
        echo "❌ Venn Panel Failed. Check log: ${LOG_DIR}/04d_venn.log"
        tail -n 3 "${LOG_DIR}/04d_venn.log"
    fi
fi

echo -e "\n========================================================"
echo "   Visualization Pipeline Completed."
echo "========================================================"