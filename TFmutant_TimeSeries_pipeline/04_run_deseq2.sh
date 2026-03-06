#!/bin/bash
# script_name: 04_run_deseq2.sh
# 功能: 运行 Pipeline 4 (DESeq2)，支持任意目录执行
# 版本: v4.1 (路径健壮性修复 + data目录输出)

# ========================================================
# 1. 自动获取脚本所在目录 (核心修复)
# ========================================================
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# 加载全局配置 (使用绝对路径)
if [ -f "${SCRIPT_DIR}/project_config.sh" ]; then
    source "${SCRIPT_DIR}/project_config.sh"
else
    echo "❌ Error: project_config.sh not found in ${SCRIPT_DIR}"
    exit 1
fi

# ========================================================
# 2. 输入输出路径
# ========================================================
COUNT_MATRIX="${BASE_DIR}/02_counts/${PROJECT_NAME}.count_matrix.txt"
OUT_DIR="${BASE_DIR}/04_DESeq2"

# 检查输入文件
if [ ! -f "${COUNT_MATRIX}" ]; then
    echo "❌ Error: Count matrix not found at ${COUNT_MATRIX}"
    exit 1
fi

# ========================================================
# 3. 创建目录结构
# ========================================================
mkdir -p "${OUT_DIR}/RData"
mkdir -p "${OUT_DIR}/Logic1_TimeCourse"
mkdir -p "${OUT_DIR}/Logic2_GenotypeComp"
mkdir -p "${OUT_DIR}/data"  # <--- 之前补充的 data 目录

echo "================================================="
echo " Pipeline 4: DESeq2 Differential Analysis"
echo " Logic 1: Time-Course (Base: ${CONTROL_TIME})"
echo " Logic 2: Cross-Genotype (Ref: ${CONTROL_CONDITION})"
echo " Script Dir: ${SCRIPT_DIR}"
echo "================================================="

# ========================================================
# 4. 执行 R 脚本 (使用绝对路径)
# ========================================================
# 核心修复: 指定 R 脚本的完整路径
Rscript "${SCRIPT_DIR}/RNAseq_pipeline4_DESeq2.R" \
    --count_file "${COUNT_MATRIX}" \
    --metadata_file "${METADATA_FILE}" \
    --out_dir "${OUT_DIR}" \
    --project_name "${PROJECT_NAME}" \
    --control_cond "${CONTROL_CONDITION}" \
    --control_time "${CONTROL_TIME}" \
    --fdr "${FDR_CUTOFF}" \
    --lfc "${LFC_THRESHOLD}"

# ========================================================
# 5. 结果检查
# ========================================================
if [ $? -eq 0 ]; then
    echo "✅ Pipeline 4 Success."
    echo "   Files generated:"
    echo "   1. Logic 1 RData: ${OUT_DIR}/RData/Pipeline4_Logic1_TimeCourse_for_Mfuzz.RData"
    echo "   2. Master RData:  ${OUT_DIR}/RData/${PROJECT_NAME}_DESeq2_objects.RData"
    echo "   3. VST Matrix:    ${OUT_DIR}/data/${PROJECT_NAME}_VST_transformed.txt"
else
    echo "❌ Pipeline 4 Failed."
    exit 1
fi