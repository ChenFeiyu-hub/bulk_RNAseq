#!/bin/bash
# script_name: 03_run_stringtie.sh

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
OUT_DIR="${BASE_DIR}/03_tpm"

echo ">>> 启动 Pipeline 3: StringTie TPM 计算 (R Version)"


# 1. 执行 R 脚本
# 注意: 这里的传参风格已经适配了 R 的 optparse
Rscript "${SCRIPT_DIR}/RNAseq_pipeline3_StringTie_TPM.R" \
    --bam_dir "${IN_DIR}" \
    --out_dir "${OUT_DIR}" \
    --gtf "${GTF_FILE}" \
    --stringtie_bin "${STRINGTIE_BIN}" \
    --project_name "${PROJECT_NAME}" \
    --threads "${STRINGTIE_THREADS}" \
    --jobs "${PARALLEL_JOBS}"

# 2. 捕获并检查退出码
EXIT_CODE=$?

if [ $EXIT_CODE -ne 0 ]; then
    echo "❌ [Error] R script failed with exit code: $EXIT_CODE"
    exit $EXIT_CODE
fi

# 3. 只有成功时才打印完成信息
echo ">>> Pipeline 3 完成。"