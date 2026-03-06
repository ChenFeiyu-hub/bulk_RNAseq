#!/bin/bash
# script_name: 05_run_mfuzz.sh
# 功能: 运行 Mfuzz 时间序列聚类分析 (Pipeline 5a/5b/5c)
# 版本: v11.1 (Standardized: Auto-parse italics & unified variable names)

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
# 2. 智能参数解析函数 (Core Logic for v11.1)
# ========================================================
# 功能: 解析 "italic(Name)" 格式，分离 Display Name 和 Italic Boolean
parse_display_format() {
    local input_str="$1"
    local clean_name=""
    local use_italic="FALSE"

    # 正则匹配 italic(...)
    if [[ "$input_str" =~ ^italic\((.+)\)$ ]]; then
        clean_name="${BASH_REMATCH[1]}"
        use_italic="TRUE"
    else
        clean_name="$input_str"
        use_italic="FALSE"
    fi

    # 返回两个值，以分号分隔 (Shell trick to return multiple values)
    echo "${clean_name};${use_italic}"
}

# ========================================================
# 3. 输入输出路径
# ========================================================
RDATA_FILE="${BASE_DIR}/04_DESeq2/RData/Pipeline4_Logic1_TimeCourse_for_Mfuzz.RData"
OUT_DIR="${BASE_DIR}/05_Mfuzz"

# ========================================================
# 4. 准备参数 (解析 Config)
# ========================================================

# --- A. 解析 WT 显示格式 ---
RAW_WT_DISPLAY="${CONDITION_WT_DISPLAY:-WT}"
PARSED_WT=$(parse_display_format "$RAW_WT_DISPLAY")
# 读取解析结果
VAL_WT_DISPLAY="${PARSED_WT%;*}"
VAL_WT_ITALIC="${PARSED_WT#*;}"

# --- B. 解析 MUT 显示格式 ---
RAW_MUT_DISPLAY="${CONDITION_MUT_DISPLAY:-Mutant}"
PARSED_MUT=$(parse_display_format "$RAW_MUT_DISPLAY")
# 读取解析结果
VAL_MUT_DISPLAY="${PARSED_MUT%;*}"
VAL_MUT_ITALIC="${PARSED_MUT#*;}"

# --- C. 其他参数 ---
COND_WT="${CONDITION_WT:-WT}"
COND_MUT="${CONDITION_MUT:-dof}"  # 对应 Config 中的真实 Sample 标签

# 聚类数 (兼容旧版变量名，优先使用 Config v11.1)
C_WT="${MFUZZ_C_WT:-}"
C_MUT="${MFUZZ_C_MUT:-}"

# 核心 Mfuzz 参数
MIN_STD="${MFUZZ_MIN_STD:-0.25}"
MEMBERSHIP_CUTOFF="${MFUZZ_MEMBERSHIP_CUTOFF:-0.5}"

# 运行模式
RUN_5A="${RUN_PIPELINE5A:-TRUE}"
RUN_5B="${RUN_PIPELINE5B:-TRUE}"
RUN_5C="${RUN_PIPELINE5C:-TRUE}"

# ========================================================
# 5. 帮助信息
# ========================================================
usage() {
    echo "Usage: $0 [options]"
    echo ""
    echo "Pipeline 5: Mfuzz Time-Series Clustering Analysis (v11.1 Standardized)"
    echo ""
    echo "Options:"
    echo "  --rdata PATH        Specify custom RData path"
    echo "  --out_dir PATH      Specify custom output directory"
    echo "  --c_wt N            Force WT cluster number"
    echo "  --c_mut N           Force mutant cluster number"
    echo "  --min_std N         Override min_std threshold (Default: ${MIN_STD})"
    echo "  --no-5a             Skip Pipeline5a (WT-MUT matching)"
    echo "  --no-5b             Skip Pipeline5b (fate divergence)"
    echo "  --no-5c             Skip Pipeline5c (temporal cascade)"
    echo "  -h, --help          Show this help message"
    exit 1
}

# ========================================================
# 6. 命令行参数解析 (允许覆盖 Config)
# ========================================================
while [[ $# -gt 0 ]]; do
    case $1 in
        --rdata)
            RDATA_FILE="$2"
            shift 2
            ;;
        --out_dir)
            OUT_DIR="$2"
            shift 2
            ;;
        --c_wt)
            C_WT="$2"
            shift 2
            ;;
        --c_mut|--c_dof) # 兼容旧参数名
            C_MUT="$2"
            shift 2
            ;;
        --min_std)
            MIN_STD="$2"
            shift 2
            ;;
        --no-5a)
            RUN_5A="FALSE"
            shift
            ;;
        --no-5b)
            RUN_5B="FALSE"
            shift
            ;;
        --no-5c)
            RUN_5C="FALSE"
            shift
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# ========================================================
# 7. 验证输入
# ========================================================
if [ ! -f "${RDATA_FILE}" ]; then
    echo "❌ Error: RData file not found: ${RDATA_FILE}"
    echo "Please ensure Pipeline 4 (DESeq2) has completed successfully."
    exit 1
fi

# ========================================================
# 8. 运行 Pipeline 5a & 5b (Mfuzz Core)
# ========================================================

# 构建参数数组 (适配 R v11.1 接口)
# 构建参数数组 (适配 R v11.1 接口)
ARGS_5A=(
    "--rdata" "${RDATA_FILE}"
    "--out_dir" "${OUT_DIR}"
    
    # --- 条件标签 (Filtering) ---
    "--condition_wt" "${COND_WT}"
    "--condition_mut" "${COND_MUT}"
    
    # --- 显示控制 (Visualization) ---
    "--condition_wt_display" "${VAL_WT_DISPLAY}"
    "--condition_mut_display" "${VAL_MUT_DISPLAY}"
    "--use_italic_wt" "${VAL_WT_ITALIC}"
    "--use_italic_mut" "${VAL_MUT_ITALIC}"
    
    # --- 运行模式 ---
    "--run_pipeline5a" "${RUN_5A}"
    "--run_pipeline5b" "${RUN_5B}"
    
    # --- 核心算法参数 ---
    "--min_std" "${MIN_STD}"
    "--membership_cutoff" "${MEMBERSHIP_CUTOFF}"
    
    # 补充 Pipeline 5a 的聚类搜索范围
    "--c_range_min" "${MFUZZ_C_RANGE_MIN:-6}"
    "--c_range_max" "${MFUZZ_C_RANGE_MAX:-24}"

    # Pipeline 5b 的子聚类搜索范围
    "--p5b_c_range_min" "${MFUZZ_P5B_C_MIN:-2}"
    "--p5b_c_range_max" "${MFUZZ_P5B_C_MAX:-6}"
    
    # --- DDTW 权重 (Config) ---
    "--w_shape" "${MFUZZ_W_SHAPE:-0.5}"
    "--w_amp" "${MFUZZ_W_AMP:-0.3}"
    "--w_phase" "${MFUZZ_W_PHASE:-0.2}"
    
    # --- 分类阈值 (Config) ---
    "--classify_shape_threshold" "${MFUZZ_CLASSIFY_SHAPE:-0.15}"
    "--classify_amp_threshold" "${MFUZZ_CLASSIFY_AMP:-0.12}"
    "--classify_phase_threshold" "${MFUZZ_CLASSIFY_PHASE:-0.25}"
)

# 添加可选的固定聚类数
if [ -n "${C_WT}" ]; then
    ARGS_5A+=("--c_wt" "${C_WT}")
fi

if [ -n "${C_MUT}" ]; then
    # v11.1 R脚本 参数名已标准化为 c_mut
    ARGS_5A+=("--c_mut" "${C_MUT}")
fi

# 添加 samples.txt (如果存在)
if [ -f "${METADATA_FILE}" ]; then
    ARGS_5A+=("--samples_file" "${METADATA_FILE}")
fi

# 创建输出目录
mkdir -p "${OUT_DIR}"

echo "================================================="
echo " Step 1: Mfuzz Clustering & Fate Analysis (5a/5b)"
echo "================================================="
echo " Input RData: ${RDATA_FILE}"
echo " Output Dir:  ${OUT_DIR}"
echo ""
echo " Configuration:"
echo "   WT  Tag: ${COND_WT}  -> Display: ${VAL_WT_DISPLAY} (Italic: ${VAL_WT_ITALIC})"
echo "   MUT Tag: ${COND_MUT} -> Display: ${VAL_MUT_DISPLAY} (Italic: ${VAL_MUT_ITALIC})"
echo "   Params:  min_std=${MIN_STD}, cutoff=${MEMBERSHIP_CUTOFF}"
echo ""

# 运行 R 脚本 5a
Rscript "${SCRIPT_DIR}/RNAseq_pipeline5a_mfuzz.R" "${ARGS_5A[@]}"
RET_5A=$?

if [ ${RET_5A} -ne 0 ]; then
    echo "❌ Pipeline 5a/5b Failed. Stopping."
    exit 1
fi

# ========================================================
# 9. 运行 Pipeline 5c (Temporal Cascade)
# ========================================================

MFUZZ_RESULT_RDATA="${OUT_DIR}/rdata/Mfuzz_v11_Results.RData"

if [ "${RUN_5C}" == "TRUE" ]; then
    echo ""
    echo "================================================="
    echo " Step 2: Temporal Cascade Analysis (Pipeline 5c)"
    echo "================================================="
    
    if [ ! -f "${MFUZZ_RESULT_RDATA}" ]; then
        echo "❌ Error: 5a output not found: ${MFUZZ_RESULT_RDATA}"
        exit 1
    fi

    # 5c 脚本目前主要只需要基础参数
    ARGS_5C=(
        "--rdata" "${MFUZZ_RESULT_RDATA}"
        "--out_dir" "${OUT_DIR}/Pipeline5c_Temporal"
        "--stage_ranges" "${P5C_STAGE_RANGES}"
        "--stage_names" "${P5C_STAGE_NAMES}"
        "--stage_colors" "${P5C_STAGE_COLORS}"
        "--n_pseudotime" "${P5C_N_PSEUDOTIME:-100}"
        "--loess_span" "${P5C_LOESS_SPAN:-0.4}"
        "--membership_cutoff" "${MEMBERSHIP_CUTOFF}"
    )

    Rscript "${SCRIPT_DIR}/RNAseq_pipeline5c_temporal.R" "${ARGS_5C[@]}"
    RET_5C=$?
    
    if [ ${RET_5C} -eq 0 ]; then
        STATUS_5C="Success"
    else
        STATUS_5C="Failed"
    fi
else
    STATUS_5C="Skipped"
fi

# ========================================================
# 10. 最终结果汇总
# ========================================================
echo ""
echo "================================================="
echo " Pipeline 5 Execution Summary (v11.1)"
echo "================================================="
echo ""

if [ ${RET_5A} -eq 0 ]; then
    echo "✅ Pipeline 5a/5b: Success"
    echo "   - Results: ${MFUZZ_RESULT_RDATA}"
else
    echo "❌ Pipeline 5a/5b: Failed"
    exit 1
fi

if [ "${STATUS_5C}" == "Success" ]; then
    echo "✅ Pipeline 5c:    Success"
elif [ "${STATUS_5C}" == "Skipped" ]; then
    echo "⚪ Pipeline 5c:    Skipped"
else
    echo "❌ Pipeline 5c:    Failed"
fi

echo ""