#!/bin/bash
# script_name: 00_run_all_pipelines.sh
# 功能: RNA-seq 全流程总控 (P1-P5)，支持绝对路径与任意目录运行
# 版本: v2.1 (增强路径健壮性)

# ========================================================
# 1. 自动获取脚本所在目录
# ========================================================
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# 加载全局配置
if [ -f "${SCRIPT_DIR}/project_config.sh" ]; then
    source "${SCRIPT_DIR}/project_config.sh"
else
    echo "❌ Error: project_config.sh not found in ${SCRIPT_DIR}"
    exit 1
fi

# ========================================================
# 2. 帮助与参数解析
# ========================================================
usage() {
    echo "Usage: $0 [options]"
    echo "Options:"
    echo "  -m    Run Mapping & Counting (Pipeline 1 + 2)"
    echo "  -t    Run TPM Quantification (Pipeline 3)"
    echo "  -d    Run Differential Expression (Pipeline 4)"
    echo "  -v    Run Visualization (Pipeline 5: PCA, UpSet)"
    echo "  -a    Run ALL Pipelines (Step 1 to 5)"
    echo "  -h    Show this help message"
    exit 1
}

RUN_MAP=false
RUN_TPM=false
RUN_DEG=false
RUN_VIS=false

if [ $# -eq 0 ]; then usage; fi

while getopts "mtdvah" opt; do
  case $opt in
    m) RUN_MAP=true ;;
    t) RUN_TPM=true ;;
    d) RUN_DEG=true ;;
    v) RUN_VIS=true ;;
    a) RUN_MAP=true; RUN_TPM=true; RUN_DEG=true; RUN_VIS=true ;;
    h) usage ;;
    *) usage ;;
  esac
done

# ========================================================
# 3. 基础环境准备
# ========================================================

# 确保日志目录存在 (使用绝对路径)
LOG_DIR="${BASE_DIR}/logs"
mkdir -p "${LOG_DIR}"

clean_dir() {
    local target_dir="$1"
    if [ -d "$target_dir" ]; then
        echo "🧹 [Clean] Removing old directory: $target_dir"
        rm -rf "$target_dir"
    fi
}

check_status() {
    local step_name="$1"
    local status_code="$?" # 获取上一个命令的退出码
    
    if [ ${status_code} -ne 0 ]; then
        echo "❌ [Error] ${step_name} Failed! (Exit code: ${status_code})"
        echo "   Check logs in ${LOG_DIR}"
        exit 1
    fi
    echo "✅ [Success] ${step_name} Done."
}

echo "========================================================"
echo "   RNA-seq Pipeline Controller"
echo "   Script Dir: ${SCRIPT_DIR}"
echo "   Project: ${PROJECT_NAME}"
echo "   Logs: ${LOG_DIR}"
echo "========================================================"

# ========================================================
# 4. 模块化执行逻辑
# ========================================================

# --- 模块 1: Mapping & Counting ---
if [ "$RUN_MAP" = true ]; then
    echo -e "\n>>> [Module M] Starting Mapping & Counting..."
    clean_dir "${BASE_DIR}/01_assembly"
    clean_dir "${BASE_DIR}/02_counts"
    
    # 使用 bash 调用，确保子脚本在一个新的 shell 中运行，避免环境污染
    # 使用绝对路径调用子脚本
    
    echo "   -> Running Pipeline 1: HISAT2..."
    bash "${SCRIPT_DIR}/01_run_hisat2.sh" > "${LOG_DIR}/01_hisat2.log" 2>&1
    check_status "HISAT2"
    
    echo "   -> Running Pipeline 2: featureCounts..."
    bash "${SCRIPT_DIR}/02_run_featurecounts.sh" > "${LOG_DIR}/02_featurecounts.log" 2>&1
    check_status "featureCounts"
else
    echo -e "\n>>> [Module M] Skipped."
fi

# --- 模块 2: TPM Quantification ---
if [ "$RUN_TPM" = true ]; then
    echo -e "\n>>> [Module T] Starting TPM Quantification..."
    if [ ! -d "${BASE_DIR}/01_assembly" ]; then echo "❌ Error: 01_assembly missing!"; exit 1; fi
    clean_dir "${BASE_DIR}/03_tpm"
    
    echo "   -> Running Pipeline 3: StringTie..."
    bash "${SCRIPT_DIR}/03_run_stringtie.sh" > "${LOG_DIR}/03_stringtie.log" 2>&1
    check_status "StringTie"
else
    echo -e "\n>>> [Module T] Skipped."
fi

# --- 模块 3: Differential Expression ---
if [ "$RUN_DEG" = true ]; then
    echo -e "\n>>> [Module D] Starting Differential Analysis..."
    if [ ! -f "${BASE_DIR}/02_counts/${PROJECT_NAME}.count_matrix.txt" ]; then echo "❌ Error: Count Matrix missing!"; exit 1; fi
    clean_dir "${BASE_DIR}/04_DESeq2"
    
    echo "   -> Running Pipeline 4: DESeq2..."
    # 这里的关键是 04 脚本内部现在能自我定位了
    bash "${SCRIPT_DIR}/04_run_deseq2.sh" > "${LOG_DIR}/04_deseq2.log" 2>&1
    check_status "DESeq2"
else
    echo -e "\n>>> [Module D] Skipped."
fi

# --- 模块 4: Visualization (P5) ---
if [ "$RUN_VIS" = true ]; then
    echo -e "\n>>> [Module V] Starting Visualization..."
    RDATA="${BASE_DIR}/04_DESeq2/RData/Pipeline4_Logic1_TimeCourse_for_Mfuzz.RData"
    if [ ! -f "$RDATA" ]; then echo "❌ Error: DESeq2 RData missing!"; exit 1; fi
    
    echo "   -> Running Pipeline 5: Visualization..."
    # [核心] 传递 --script-dir 参数，确保子脚本能找到 R 代码
    bash "${SCRIPT_DIR}/04a_run_DEG_visualization.sh" --script-dir "${SCRIPT_DIR}" > "${LOG_DIR}/04a_visualization.log" 2>&1
    check_status "Visualization"
else
    echo -e "\n>>> [Module V] Skipped."
fi

echo -e "\n========================================================"
echo "   All Tasks Completed."
echo "   End Time: $(date '+%Y-%m-%d %H:%M:%S')"
echo "========================================================"