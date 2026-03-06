#!/bin/bash
# script_name: 02_run_featurecounts.sh

# ========================================================
# 自动获取脚本所在目录并加载配置
# ========================================================
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

if [ -f "${SCRIPT_DIR}/project_config.sh" ]; then
    source "${SCRIPT_DIR}/project_config.sh"
else
    echo "❌ Error: project_config.sh not found in ${SCRIPT_DIR}"
    exit 1
fi
IN_DIR="${BASE_DIR}/01_assembly"
OUT_DIR="${BASE_DIR}/02_counts"
mkdir -p "${OUT_DIR}"

COUNT_FILE="${OUT_DIR}/${PROJECT_NAME}.featureCounts.txt"
MATRIX_FILE="${OUT_DIR}/${PROJECT_NAME}.count_matrix.txt"

echo "================================================="
echo " Pipeline 2: featureCounts Quantification"
echo " Threads: ${FEATURECOUNTS_THREADS} (Max 64)"
echo "================================================="

BAM_FILES=$(find "${IN_DIR}" -name "*.sort.bam" | sort)

if [ -z "$BAM_FILES" ]; then
    echo "❌ No BAM files found in ${IN_DIR}"
    exit 1
fi

echo ">>> Running featureCounts on $(echo "$BAM_FILES" | wc -w) samples..."


$FEATURECOUNTS_BIN \
    -T ${FEATURECOUNTS_THREADS} \
    -p -t exon -g gene_id \
    -a "${GTF_FILE}" \
    -o "${COUNT_FILE}" \
    $BAM_FILES \
    > "${OUT_DIR}/${PROJECT_NAME}.featureCounts.log" 2>&1

if [ $? -eq 0 ]; then
    echo ">>> Formatting Count Matrix..."
    # 生成纯净矩阵
    grep -v "^#" "${COUNT_FILE}" | \
    cut -f 1,7- | \
    sed "s|${IN_DIR}/||g" | \
    sed "s|.sort.bam||g" \
    > "${MATRIX_FILE}"
    
    echo "✅ Success! Matrix saved to: ${MATRIX_FILE}"
else
    echo "❌ Failed. Check log."
    exit 1
fi