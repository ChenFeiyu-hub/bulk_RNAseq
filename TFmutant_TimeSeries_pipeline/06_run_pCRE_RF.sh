#!/bin/bash
# script_name: 06_run_pCRE_RF.sh
# 功能: pCRE 挖掘、建模、聚类、注释与可视化 (Pipeline 6a-6e)
# 版本: v3.0 (重构 6d/6e, 支持 FULL/CONSENSUS 双轨运行)

# ========================================================
# 1. 环境初始化
# ========================================================
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"

# 加载配置
if [ -f "${SCRIPT_DIR}/project_config.sh" ]; then
    source "${SCRIPT_DIR}/project_config.sh"
else
    echo "❌ Error: project_config.sh not found."
    exit 1
fi

# ========================================================
# 2. 命令行参数解析 (独立控制模块)
# ========================================================
# 用法示例: 
#   bash 06.sh          # 运行 Config 中设置为 TRUE 的所有步骤
#   bash 06.sh -d -e    # 仅运行 6d 和 6e
#   bash 06.sh -b       # 仅运行 6b

# 如果没有提供参数，保持 Config 中的默认设置
# 如果提供了参数，先将所有开关设为 FALSE，再根据参数开启
if [ $# -gt 0 ]; then
    RUN_PIPELINE6A="FALSE"
    RUN_PIPELINE6B="FALSE"
    RUN_PIPELINE6C="FALSE"
    RUN_PIPELINE6D="FALSE"
    RUN_PIPELINE6E="FALSE"
    
    while getopts "abcde" opt; do
        case $opt in
            a) RUN_PIPELINE6A="TRUE" ;;
            b) RUN_PIPELINE6B="TRUE" ;;
            c) RUN_PIPELINE6C="TRUE" ;;
            d) RUN_PIPELINE6D="TRUE" ;;
            e) RUN_PIPELINE6E="TRUE" ;;
            *) echo "Usage: $0 [-a] [-b] [-c] [-d] [-e]"; exit 1 ;;
        esac
    done
fi

echo "================================================="
echo " Pipeline 6: pCRE Discovery System v3.0"
echo " Project: ${PROJECT_NAME}"
echo " Execution Plan:"
echo "   [6a] Enrichment:  ${RUN_PIPELINE6A}"
echo "   [6b] RF Modeling: ${RUN_PIPELINE6B}"
echo "   [6c] Clustering:  ${RUN_PIPELINE6C}"
echo "   [6d] Integration: ${RUN_PIPELINE6D}"
echo "   [6e] Baseline:   ${RUN_PIPELINE6E}"
echo ""
echo " Run Mode: ${P6_RUN_MODE:-BOTH}"
echo "================================================="

# ========================================================
# 3. 路径推断与准备 (全局可用)
# ========================================================
# Mfuzz 路径 (Pipeline 5 输出)
MFUZZ_DIR="${BASE_DIR}/05_Mfuzz"
WT_MEMBERSHIP="${MFUZZ_DIR}/morphology_analysis/WT_Membership_Matrix.txt"
MUT_MEMBERSHIP="${MFUZZ_DIR}/morphology_analysis/MUT_Membership_Matrix.txt"
VST_FILE="${BASE_DIR}/04_DESeq2/data/${PROJECT_NAME}_VST_transformed.txt"

# 基因组 FASTA 推断
if [ -n "${GENOME_FASTA}" ]; then TARGET_FASTA="${GENOME_FASTA}";
elif [ -n "${REF_FASTA}" ]; then TARGET_FASTA="${REF_FASTA}";
else TARGET_FASTA="${INDEX_PREFIX%.hisat2}.fa"; fi

# GFF 推断
if [ -n "${GFF_FILE}" ]; then TARGET_GFF="${GFF_FILE}"; else TARGET_GFF="${GTF_FILE%.gtf}.gff3"; fi

# --- 关键：预先计算各步骤的目录名 (无论是否运行，下游步骤都需要路径) ---
PCRE_BASE_DIR="${BASE_DIR}/06_pCRE"

# 6a 输出目录名逻辑
BG_LABEL=$([ -z "${P6_ENRICH_LOG2FC}" ] && echo "GenomeWideBG" || printf "log2FC%.1f" "${P6_ENRICH_LOG2FC}")
GC_LABEL=$([ "${P6_ENRICH_GC_MATCH}" == "TRUE" ] && echo "_GCmatched" || echo "")
MODE_LABEL=$([ "${P6_PARTITIONED}" == "TRUE" ] && echo "_Partitioned" || echo "_FullLength")
DIR_NAME_6A="Pipeline6a_v8.0_${BG_LABEL}${GC_LABEL}${MODE_LABEL}"
OUTPUT_DIR_6A="${PCRE_BASE_DIR}/${DIR_NAME_6A}"

# 6b 输出目录名逻辑
SUFFIX_6B=$([ "${P6_RF_GC_MATCH}" == "TRUE" ] && echo "_GCmatched" || echo "")
DIR_NAME_6B="Pipeline6b_v4.0${SUFFIX_6B}"
OUTPUT_DIR_6B="${PCRE_BASE_DIR}/${DIR_NAME_6B}"

# 6c 输出目录名逻辑 (通常是固定的后缀)
DIR_NAME_6C="Pipeline6c_v2.0_MotifClustering"
OUTPUT_DIR_6C="${PCRE_BASE_DIR}/${DIR_NAME_6C}"

# 6d 输出目录名逻辑
DIR_NAME_6D="Pipeline6d_v2.0_Visualization"
OUTPUT_DIR_6D="${PCRE_BASE_DIR}/${DIR_NAME_6D}"

# 6e 输出目录名逻辑
DIR_NAME_6E="Pipeline6e1_v4.0_Figures"
OUTPUT_DIR_6E="${PCRE_BASE_DIR}/${DIR_NAME_6E}"

# ========================================================
# 4. Pipeline 6a: 统计富集
# ========================================================
if [ "${RUN_PIPELINE6A}" == "TRUE" ]; then
    echo -e "\n>>> Running Pipeline 6a..."
    
    # 检查输入
    if [ ! -f "${WT_MEMBERSHIP}" ] || [ ! -f "${VST_FILE}" ]; then
        echo "❌ Error: Missing input files from Pipeline 4/5."
        exit 1
    fi

    ARGS_6A=(
        "--base_dir" "${BASE_DIR}"
        "--genotype" "${P6_GENOTYPE}"
        "--condition_wt" "${CONDITION_WT}"          # ← 新增
        "--condition_mut" "${CONDITION_MUT}"         # ← 新增
        "--control_time" "${CONTROL_TIME}"           # ← 新增
        "--wt_membership" "${WT_MEMBERSHIP}"
        "--mutant_membership" "${MUT_MEMBERSHIP}"
        "--vst_file" "${VST_FILE}"
        "--cores" "${P6_CORES}"
        "--top_n" "${P6_TOP_N}"
        "--k_values" "${P6_K_VALUES}"
        "--min_occurrence" "${P6_MIN_OCCURRENCE}"
        "--ml_log2FC" "${P6_ML_LOG2FC}"
        "--membership_threshold" "${MFUZZ_MEMBERSHIP_CUTOFF}"
        "--genome_fasta" "${TARGET_FASTA}"
    )

    if [ "${P6_PARTITIONED}" == "TRUE" ]; then
        ARGS_6A+=("--partitioned")
        ARGS_6A+=("--core_start" "${P6_CORE_START}" "--core_end" "${P6_CORE_END}")
        ARGS_6A+=("--proximal_start" "${P6_PROX_START}" "--proximal_end" "${P6_PROX_END}")
        ARGS_6A+=("--distal_start" "${P6_DIST_START}" "--distal_end" "${P6_DIST_END}")
    else
        ARGS_6A+=("--full_length")
        ARGS_6A+=("--full_upstream" "${P6_FULL_UP}" "--full_downstream" "${P6_FULL_DOWN}")
    fi

    if [ -n "${P6_ENRICH_LOG2FC}" ]; then ARGS_6A+=("--log2FC" "${P6_ENRICH_LOG2FC}"); fi
    if [ "${P6_ENRICH_GC_MATCH}" == "TRUE" ]; then ARGS_6A+=("--gc_match"); fi
    if [ -f "${TARGET_GFF}" ]; then ARGS_6A+=("--gff_file" "${TARGET_GFF}"); fi
    if [ -n "${GENE_ID_STRIP_PATTERN}" ]; then
        ARGS_6A+=("--gene_id_strip" "${GENE_ID_STRIP_PATTERN}")
    fi

    Rscript "${SCRIPT_DIR}/RNAseq_Pipeline6a_pCRE_Statistical_Enrichment.R" "${ARGS_6A[@]}"
    if [ $? -ne 0 ]; then echo "❌ 6a Failed"; exit 1; fi
fi

# ========================================================
# 5. Pipeline 6b: 随机森林建模
# ========================================================
if [ "${RUN_PIPELINE6B}" == "TRUE" ]; then
    echo -e "\n>>> Running Pipeline 6b..."
    
    if [ ! -d "${OUTPUT_DIR_6A}" ]; then
        echo "❌ Error: 6a output directory not found: ${OUTPUT_DIR_6A}"
        exit 1
    fi

    ARGS_6B=(
        "--input_dir" "${OUTPUT_DIR_6A}"
        "--genotype" "${P6_GENOTYPE}"
        "--ntree" "${P6_RF_NTREE}"
        "--replicates" "${P6_RF_REPLICATES}"
        "--folds" "${P6_RF_FOLDS}"
        "--gc_bins" "${P6_RF_GC_BINS}"
        "--neg_ratio" "${P6_RF_NEG_RATIO}"
        "--cores" "${P6_CORES}"
    )
    if [ "${P6_RF_GC_MATCH}" == "TRUE" ]; then ARGS_6B+=("--gc_match"); fi

    Rscript "${SCRIPT_DIR}/RNAseq_Pipeline6b_pCRE_RF_Model.R" "${ARGS_6B[@]}"
    if [ $? -ne 0 ]; then echo "❌ 6b Failed"; exit 1; fi
fi

# ========================================================
# 6. Pipeline 6c v3.0: Unified Motif Integration & Annotation
# ========================================================
# 功能: Levenshtein聚类 -> Shift-Aware加权合并 -> 内部PCC注释
# 取代: 旧版 6c (聚类) 和 6d (Tomtom整合)

if [ "${RUN_PIPELINE6C}" == "TRUE" ]; then
    echo -e "\n>>> Running Pipeline 6c v3.0 (Integration & Annotation)..."
    
    # 检查输入目录
    # 6a 输出目录 (Enrichment)
    if [ ! -d "${OUTPUT_DIR_6A}" ]; then
        echo "❌ Error: 6a output directory not found: ${OUTPUT_DIR_6A}"
        exit 1
    fi
    # 6b 输出目录 (Importance - Master List)
    if [ ! -d "${OUTPUT_DIR_6B}" ]; then
        echo "❌ Error: 6b output directory not found: ${OUTPUT_DIR_6B}"
        exit 1
    fi
    # 检查 MEME 数据库
    if [ -z "${MEME_DB_PATH}" ] || [ ! -f "${MEME_DB_PATH}" ]; then
        echo "❌ Error: MEME database not found at: ${MEME_DB_PATH}"
        echo "   Please check MEME_DB_PATH in project_config.sh"
        exit 1
    fi

    # 更新输出目录名 (v3.0)
    DIR_NAME_6C="Pipeline6c_v3.0_MotifIntegration"
    OUTPUT_DIR_6C="${PCRE_BASE_DIR}/${DIR_NAME_6C}"

    # 构建参数
    ARGS_6C=(
        "--p6a_dir" "${OUTPUT_DIR_6A}"
        "--p6b_dir" "${OUTPUT_DIR_6B}"
        "--output_dir" "${OUTPUT_DIR_6C}"
        "--meme_db" "${MEME_DB_PATH}"
        "--genotype" "${P6_GENOTYPE}"
        "--n_cores" "${P6_CORES}"
        "--cluster_distance" "${P6_CLUST_DIST}"
        "--min_overlap" "${P6_CLUST_MIN_OVERLAP}"
        "--trim_ic" "${P6_CLUST_TRIM_IC}"
        "--seed" "42"
    )

    echo "  Input 6a: ${OUTPUT_DIR_6A}"
    echo "  Input 6b: ${OUTPUT_DIR_6B}"
    echo "  Output:   ${OUTPUT_DIR_6C}"
    echo "  Database: ${MEME_DB_PATH}"

    # 执行 R 脚本
    # 确保脚本文件名与您保存的一致
    Rscript "${SCRIPT_DIR}/RNAseq_Pipeline6c_Motif_Integration.R" "${ARGS_6C[@]}"
    
    if [ $? -ne 0 ]; then 
        echo "❌ 6c v3.0 Failed"; exit 1; 
    else
        echo "✅ 6c v3.0 Completed Successfully"
    fi
fi

# ========================================================
# 7. Pipeline 6d v2.0: Publication-Quality Visualization
# ========================================================
# 功能: AUC柱状图(WT左/Mut右) + 下三角正交热图(TopN特征) + DB-Anchored Motif报告
# 输入: 6b (AUC), 6c (Motifs+PWMs), MEME_DB
# 新增: --wt_display_label, --mut_display_label, --wt_color, --mut_color,
#        --top_n_features, --target_keyword

if [ "${RUN_PIPELINE6D}" == "TRUE" ]; then
    echo -e "\n>>> Running Pipeline 6d (Visualization v2.0)..."
    
    # --- 检查输入目录 ---
    INPUT_DIR_6C="${PCRE_BASE_DIR}/Pipeline6c_v3.0_MotifIntegration"
    
    if [ ! -d "${INPUT_DIR_6C}" ]; then
        echo "❌ Error: 6c output directory not found: ${INPUT_DIR_6C}"
        exit 1
    fi

    if [ ! -d "${OUTPUT_DIR_6B}" ]; then
        echo "❌ Error: 6b output directory not found: ${OUTPUT_DIR_6B}"
        exit 1
    fi

    if [ -z "${MEME_DB_PATH}" ] || [ ! -f "${MEME_DB_PATH}" ]; then
        echo "❌ Error: MEME database not found at: ${MEME_DB_PATH}"
        echo "   Visualization (Fig 3) requires the DB for motif alignment."
        exit 1
    fi

    # --- 定义输出目录 ---
    DIR_NAME_6D="Pipeline6d_v2.0_Visualization"
    OUTPUT_DIR_6D="${PCRE_BASE_DIR}/${DIR_NAME_6D}"
    mkdir -p "${OUTPUT_DIR_6D}"

    # --- 解析显示标签 (优先使用 6d 独立配置, 否则复用 Pipeline 5 的全局配置) ---
    VIS_WT_LABEL="${P6_VIS_WT_LABEL:-${CONDITION_WT_DISPLAY:-WT}}"
    VIS_MUT_LABEL="${P6_VIS_MUT_LABEL:-${CONDITION_MUT_DISPLAY:-Mutant}}"

    # --- 解析颜色 (优先使用 6d 独立配置, 否则从 STYLE_COND_COLORS 自动提取) ---
    # 默认回退颜色
    DEFAULT_WT_COLOR="#2E86AB"
    DEFAULT_MUT_COLOR="#A23B72"

    if [ -n "${P6_VIS_WT_COLOR}" ]; then
        VIS_WT_COLOR="${P6_VIS_WT_COLOR}"
    elif [ -n "${STYLE_COND_COLORS}" ]; then
        # 从 "WT:#2E86AB,dof:#A23B72" 中提取 WT 的颜色
        VIS_WT_COLOR=$(echo "${STYLE_COND_COLORS}" | tr ',' '\n' | grep "^${CONDITION_WT}:" | cut -d: -f2)
        VIS_WT_COLOR="${VIS_WT_COLOR:-${DEFAULT_WT_COLOR}}"
    else
        VIS_WT_COLOR="${DEFAULT_WT_COLOR}"
    fi

    if [ -n "${P6_VIS_MUT_COLOR}" ]; then
        VIS_MUT_COLOR="${P6_VIS_MUT_COLOR}"
    elif [ -n "${STYLE_COND_COLORS}" ]; then
        VIS_MUT_COLOR=$(echo "${STYLE_COND_COLORS}" | tr ',' '\n' | grep "^${CONDITION_MUT}:" | cut -d: -f2)
        VIS_MUT_COLOR="${VIS_MUT_COLOR:-${DEFAULT_MUT_COLOR}}"
    else
        VIS_MUT_COLOR="${DEFAULT_MUT_COLOR}"
    fi

    # --- 构建参数 ---
    ARGS_6D=(
        "--p6c_dir"            "${INPUT_DIR_6C}"
        "--p6b_dir"            "${OUTPUT_DIR_6B}"
        "--output_dir"         "${OUTPUT_DIR_6D}"
        "--meme_db"            "${MEME_DB_PATH}"
        "--genotype"           "${P6_GENOTYPE}"
        "--wt_display_label"   "${VIS_WT_LABEL}"
        "--mut_display_label"  "${VIS_MUT_LABEL}"
        "--wt_color"           "${VIS_WT_COLOR}"
        "--mut_color"          "${VIS_MUT_COLOR}"
        "--plot_top_n"         "${P6_VIS_TOP_N:-3}"
        "--top_n_features"     "${P6_VIS_TOP_N_FEATURES:-100}"
        "--heatmap_metric"     "${P6_VIS_HEATMAP_METRIC:-overlap}"
        "--fig_width"          "${P6_VIS_FIG_WIDTH:-12}"
        "--fig_height"         "${P6_VIS_FIG_HEIGHT:-8}"
        "--seed"               "42"
    )

    # 可选: target_keyword (仅在非空时传递)
    if [ -n "${P6_VIS_TARGET_KEYWORD}" ]; then
        ARGS_6D+=("--target_keyword" "${P6_VIS_TARGET_KEYWORD}")
    fi

    # 可选: target_scope (仅在非空且非默认ALL时传递，或始终传递以明确意图)
    if [ -n "${P6_VIS_TARGET_SCOPE}" ]; then
        ARGS_6D+=("--target_scope" "${P6_VIS_TARGET_SCOPE}")
    fi

    echo "  Input 6c : ${INPUT_DIR_6C}"
    echo "  Input 6b : ${OUTPUT_DIR_6B}"
    echo "  Output   : ${OUTPUT_DIR_6D}"
    echo "  WT Label : ${VIS_WT_LABEL}  (color: ${VIS_WT_COLOR})"
    echo "  Mut Label: ${VIS_MUT_LABEL} (color: ${VIS_MUT_COLOR})"
    echo "  Top N    : scatter=${P6_VIS_TOP_N:-3}, heatmap=${P6_VIS_TOP_N_FEATURES:-100}"
    if [ -n "${P6_VIS_TARGET_KEYWORD}" ]; then
        echo "  Target   : ${P6_VIS_TARGET_KEYWORD}"
    fi
    if [ -n "${P6_VIS_TARGET_SCOPE}" ]; then
        echo "  Scope    : ${P6_VIS_TARGET_SCOPE}"
    fi

    # 执行 R 脚本
    Rscript "${SCRIPT_DIR}/RNAseq_Pipeline6d_Visualization.R" "${ARGS_6D[@]}"
    
    if [ $? -ne 0 ]; then 
        echo "❌ 6d Visualization Failed"; exit 1; 
    else
        echo "✅ 6d Visualization v2.0 Completed Successfully"
        echo "   Figures: ${OUTPUT_DIR_6D}/Fig1_AUC_ROC_Barplot.pdf"
        echo "            ${OUTPUT_DIR_6D}/Fig2_Similarity_Heatmap.pdf"
        echo "   Reports: ${OUTPUT_DIR_6D}/Cluster_Reports/"
        echo "   CSV:     ${OUTPUT_DIR_6D}/Selected_Motifs_List.csv"
        echo "   TSV:     ${OUTPUT_DIR_6D}/Cluster_Reports/*/Top20_Motif_Details.tsv"
    fi
fi

# ========================================================
# 8. Pipeline 6e v1.0: Permutation Baseline Test
# ========================================================
# 功能: 标签置换检验，验证 6b 的 AUC 是否显著优于随机基准
# 输入: 6a (ML_ready 数据) + 6b (真实 AUC)

if [ "${RUN_PIPELINE6E}" == "TRUE" ]; then
    echo -e "\n>>> Running Pipeline 6e (Permutation Baseline Test v1.0)..."
    
    # 检查输入
    if [ ! -d "${OUTPUT_DIR_6A}" ]; then
        echo "❌ Error: 6a output directory not found: ${OUTPUT_DIR_6A}"
        exit 1
    fi
    if [ ! -d "${OUTPUT_DIR_6B}" ]; then
        echo "❌ Error: 6b output directory not found: ${OUTPUT_DIR_6B}"
        exit 1
    fi

    # 定义输出目录
    DIR_NAME_6E="Pipeline6e_v1.0_PermutationTest"
    OUTPUT_DIR_6E="${PCRE_BASE_DIR}/${DIR_NAME_6E}"

    ARGS_6E=(
        "--p6a_dir"      "${OUTPUT_DIR_6A}"
        "--p6b_dir"      "${OUTPUT_DIR_6B}"
        "--output_dir"   "${OUTPUT_DIR_6E}"
        "--genotype"     "${P6_GENOTYPE}"
        "--n_perm"       "${P6_PERM_N:-100}"
        "--n_replicates" "${P6_PERM_REPLICATES:-10}"
        "--n_folds"      "${P6_PERM_FOLDS:-5}"
        "--ntree"        "${P6_PERM_NTREE:-200}"
        "--neg_ratio"    "${P6_RF_NEG_RATIO}"
        "--cores"        "${P6_CORES}"
        "--alpha"        "${P6_PERM_ALPHA:-0.05}"
    )
    if [ "${P6_RF_GC_MATCH}" == "TRUE" ]; then ARGS_6E+=("--gc_match"); fi
    ARGS_6E+=("--gc_bins" "${P6_RF_GC_BINS}")

    echo "  Input 6a: ${OUTPUT_DIR_6A}"
    echo "  Input 6b: ${OUTPUT_DIR_6B}"
    echo "  Output:   ${OUTPUT_DIR_6E}"
    echo "  Permutations: ${P6_PERM_N:-100}"
    echo "  CV Spec: ${P6_PERM_REPLICATES:-10} reps × ${P6_PERM_FOLDS:-5}-fold"
    echo "  ntree: ${P6_PERM_NTREE:-200}"

    Rscript "${SCRIPT_DIR}/RNAseq_Pipeline6e_Permutation_Baseline.R" "${ARGS_6E[@]}"
    
    if [ $? -ne 0 ]; then 
        echo "❌ 6e Permutation Test Failed"; exit 1; 
    else
        echo "✅ 6e Permutation Baseline Test v1.0 Completed Successfully"
        echo "   Summary: ${OUTPUT_DIR_6E}/*/Permutation_Summary.txt"
        echo "   Figures: ${OUTPUT_DIR_6E}/*/Fig_Permutation_Test_*.pdf"
    fi
fi
# ========================================================
# 9. 完成汇总
# ========================================================
echo ""
echo "================================================="
echo " Pipeline 6 Execution Summary"
echo "================================================="
echo " Output Directories:"
echo "   6a: ${OUTPUT_DIR_6A}"
echo "   6b: ${OUTPUT_DIR_6B}"
echo "   6c: ${OUTPUT_DIR_6C}"
echo "   6d: ${OUTPUT_DIR_6D}"
echo "   6e: ${OUTPUT_DIR_6E}"
echo ""
echo " Run Mode: ${P6_RUN_MODE:-BOTH}"
if [ "${P6_RUN_MODE}" == "BOTH" ] || [ -z "${P6_RUN_MODE}" ]; then
    echo "   - Full_Features: All filtered k-mers"
    echo "   - Degenerate_Features: Degenerate merged motifs"
fi
echo ""
echo " Completed at: $(date '+%Y-%m-%d %H:%M:%S')"
echo "================================================="