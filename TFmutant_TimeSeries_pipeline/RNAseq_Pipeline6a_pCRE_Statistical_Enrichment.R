#!/usr/bin/env Rscript
################################################################################
# RNAseq Pipeline 6a v8.0: Unified Multi-Region pCRE Enrichment Analysis
#
# 功能特性：
#   - 支持WT和Mutant基因型（可单独或同时处理）
#   - 支持分区模式(Partitioned)与全长模式(Full-length)，互斥
#   - 可配置的启动子区域边界参数
#   - 标准化输出目录结构，配套6b流程
#   - 完全移除FDR过滤逻辑
#
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(Biostrings)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(Rsamtools)
  library(rtracklayer)
  library(parallel)
  library(doParallel)
  library(foreach)
  library(optparse)
})

################################################################################
# 命令行参数解析
################################################################################

option_list <- list(
  # === 基因型选择 ===
  make_option(c("--genotype"), type = "character", default = "both",
              help = "分析的基因型: 'WT', 'Mutant', 或 'both' [default: %default]"),
  
  # === 【新增】条件名与时间基线 ===
  make_option(c("--condition_wt"), type = "character", default = NULL,
              help = "WT条件名 (必须与samples.txt/VST列名一致, 如 'MpWT', 'CrWT')"),
  make_option(c("--condition_mut"), type = "character", default = NULL,
              help = "Mutant条件名 (如 'Mpstop1', 'Crdof')"),
  make_option(c("--control_time"), type = "character", default = NULL,
              help = "基线时间标签 (如 '0d', '00h'), 必须与样本命名一致"),
  make_option(c("--gene_id_strip"), type = "character", default = "",
              help = "基因ID清理正则 (如 '_4532(\\\\.v6\\\\.1)?$'), 留空不剥离"),
  
  # === 区域模式选择（互斥） ===
  make_option(c("--partitioned"), action = "store_true", default = TRUE,
              help = "启用分区模式（默认）"),
  make_option(c("--full_length"), action = "store_true", default = FALSE,
              help = "启用全长模式（与--partitioned互斥）"),
  
  # === 全长模式参数 ===
  make_option(c("--full_upstream"), type = "integer", default = 1500,
              help = "全长模式: TSS上游距离(bp) [default: %default]"),
  make_option(c("--full_downstream"), type = "integer", default = 500,
              help = "全长模式: TSS下游距离(bp) [default: %default]"),
  
  # === 分区模式参数 ===
  make_option(c("--core_start"), type = "integer", default = -50,
              help = "核心启动子起始位置(相对TSS) [default: %default]"),
  make_option(c("--core_end"), type = "integer", default = 200,
              help = "核心启动子结束位置(相对TSS) [default: %default]"),
  make_option(c("--proximal_start"), type = "integer", default = -500,
              help = "近端启动子起始位置 [default: %default]"),
  make_option(c("--proximal_end"), type = "integer", default = -50,
              help = "近端启动子结束位置 [default: %default]"),
  make_option(c("--distal_start"), type = "integer", default = -1500,
              help = "远端启动子起始位置 [default: %default]"),
  make_option(c("--distal_end"), type = "integer", default = -500,
              help = "远端启动子结束位置 [default: %default]"),
  
  # === 富集分析参数 ===
  make_option(c("--log2FC"), type = "double", default = NULL,
              help = "富集分析的Log2FC阈值。若不指定则使用全基因组背景"),
  make_option(c("--ml_log2FC"), type = "double", default = 0.8,
              help = "ML模型的Log2FC阈值 [default: %default]"),
  make_option(c("--gc_match"), action = "store_true", default = FALSE,
              help = "是否在富集分析中使用GC匹配背景"),
  
  # === 路径配置 ===
  make_option(c("--base_dir"), type = "character", default = NULL,
              help = "项目基础目录 (必需)"),
  make_option(c("--wt_membership"), type = "character", default = NULL,
              help = "WT Membership Matrix文件路径"),
  make_option(c("--mutant_membership"), type = "character", default = NULL,
              help = "Mutant Membership Matrix文件路径"),
  make_option(c("--vst_file"), type = "character", default = NULL,
              help = "VST transformed表达矩阵文件路径"),
  make_option(c("--genome_fasta"), type = "character", default = NULL,
              help = "基因组FASTA文件路径"),
  make_option(c("--gff_file"), type = "character", default = NULL,
              help = "GFF3注释文件路径"),
  
  # === 计算参数 ===
  make_option(c("--cores"), type = "integer", default = 40,
              help = "并行核心数 [default: %default]"),
  make_option(c("--top_n"), type = "integer", default = 100,
              help = "每个Cluster/Region提取的Top N pCREs [default: %default]"),
  make_option(c("--membership_threshold"), type = "double", default = 0.5,
              help = "Membership阈值 [default: %default]"),
  make_option(c("--min_occurrence"), type = "integer", default = 5,
              help = "k-mer最小出现次数 [default: %default]"),
  make_option(c("--k_values"), type = "character", default = "6",
              help = "k-mer长度，逗号分隔 [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

################################################################################
# 参数验证与处理
################################################################################

# --- 必需参数检查 ---
if (is.null(opt$base_dir)) stop("必须提供 --base_dir")
if (is.null(opt$condition_wt)) stop("必须提供 --condition_wt (如 'MpWT')")
if (is.null(opt$condition_mut)) stop("必须提供 --condition_mut (如 'Mpstop1')")
if (is.null(opt$control_time)) stop("必须提供 --control_time (如 '0d' 或 '00h')")

# 处理互斥参数
if (opt$full_length && opt$partitioned) {
  opt$partitioned <- FALSE  # full_length优先
}

# 解析k_values
k_values <- as.integer(strsplit(opt$k_values, ",")[[1]])

# 验证基因型参数
valid_genotypes <- c("WT", "Mutant", "both")
if (!opt$genotype %in% valid_genotypes) {
  stop(sprintf("无效的基因型参数: %s. 有效值: %s", 
               opt$genotype, paste(valid_genotypes, collapse = ", ")))
}

# --- 从 control_time 中提取数值基线 ---
# 约定: 时间标签 = 数字 + 单位字母(s), 如 "0d", "00h", "12h", "3d"
# 提取纯数字部分作为数值比较的基线
CONTROL_TIME_LABEL <- opt$control_time
CONTROL_TIME_NUMERIC <- as.numeric(gsub("[^0-9.]", "", CONTROL_TIME_LABEL))
cat(sprintf("  基线时间: '%s' → 数值 %g\n", CONTROL_TIME_LABEL, CONTROL_TIME_NUMERIC))

################################################################################
# 参数配置 (纯参数驱动, 无硬编码默认路径)
################################################################################

# 1. 区域配置构建函数 (保留了你原本的逻辑)
build_region_config <- function(opt) {
  if (opt$full_length) {
    # 全长模式：单一区域
    return(list(
      FullLength = list(
        name = "FullLength",
        description = sprintf("TSS上游%dbp至下游%dbp", opt$full_upstream, opt$full_downstream),
        upstream = opt$full_upstream,
        downstream = opt$full_downstream
      )
    ))
  } else {
    # 分区模式：多区域
    regions <- list()
    
    # 核心启动子
    core_up <- abs(min(opt$core_start, opt$core_end))
    core_down <- max(0, max(opt$core_start, opt$core_end))
    regions$Core <- list(
      name = "Core",
      description = sprintf("核心启动子(%d to %d)", opt$core_start, opt$core_end),
      upstream = core_up,
      downstream = core_down,
      region_start = opt$core_start,
      region_end = opt$core_end
    )
    
    # 近端启动子
    proximal_up <- abs(min(opt$proximal_start, opt$proximal_end))
    proximal_down <- max(0, max(opt$proximal_start, opt$proximal_end))
    # 检查是否需要trim
    if (opt$proximal_end > opt$distal_end) {
      regions$Proximal <- list(
        name = "Proximal",
        description = sprintf("近端启动子(%d to %d)", opt$proximal_start, opt$proximal_end),
        upstream = proximal_up,
        downstream = proximal_down,
        region_start = opt$proximal_start,
        region_end = opt$proximal_end
      )
    }
    
    # 远端启动子
    distal_up <- abs(min(opt$distal_start, opt$distal_end))
    regions$Distal <- list(
      name = "Distal",
      description = sprintf("远端启动子(%d to %d)", opt$distal_start, opt$distal_end),
      upstream = distal_up,
      downstream = 0,
      trim_start = abs(opt$distal_end),  # 需要trim掉的部分
      region_start = opt$distal_start,
      region_end = opt$distal_end
    )
    
    return(regions)
  }
}

# 2. 文件路径解析函数: 全部由命令行或Shell传入, 引入泛用回退逻辑
resolve_path <- function(opt_value, fallback_dir, fallback_name, label) {
  if (!is.null(opt_value) && opt_value != "") {
    return(opt_value)
  }
  # 通用回退: 基于标准目录结构推断
  fb <- file.path(fallback_dir, fallback_name)
  if (file.exists(fb)) {
    cat(sprintf("  [自动推断] %s → %s\n", label, fb))
    return(fb)
  }
  stop(sprintf("必须通过命令行提供 --%s 或确保文件存在: %s", label, fb))
}

GENE_ID_STRIP <- opt$gene_id_strip  # 全局基因ID剥离模式

# 3. 参数汇总列表 (结合了 Claude 的变量和你的区域函数)
PARAMS <- list(
  base_dir = opt$base_dir,
  genotypes_to_process = if (opt$genotype == "both") c("WT", "Mutant") else opt$genotype,
  
  # 条件名 (由Config传入, 泛用)
  condition_wt = opt$condition_wt,
  condition_mut = opt$condition_mut,
  control_time_label = CONTROL_TIME_LABEL,
  control_time_numeric = CONTROL_TIME_NUMERIC,
  gene_id_strip = GENE_ID_STRIP,
  
  # 文件路径 (全部来自命令行, 无硬编码)
  wt_membership_file = resolve_path(opt$wt_membership, opt$base_dir, 
                                     "05_Mfuzz/morphology_analysis/WT_Membership_Matrix.txt",
                                     "wt_membership"),
  mutant_membership_file = resolve_path(opt$mutant_membership, opt$base_dir,
                                         "05_Mfuzz/morphology_analysis/MUT_Membership_Matrix.txt",
                                         "mutant_membership"),
  vst_file = opt$vst_file,
  genome_fasta = opt$genome_fasta,
  gff_file = opt$gff_file,
  
  # 区域配置
  region_mode = ifelse(opt$full_length, "full_length", "partitioned"),
  genomic_regions = build_region_config(opt),
  
  # 分析参数
  membership_threshold = opt$membership_threshold,
  enrichment_log2fc = opt$log2FC,
  use_gc_match_enrichment = opt$gc_match,
  ml_log2fc_threshold = opt$ml_log2FC,
  min_expression = 1,
  k_values = k_values,
  min_occurrence = opt$min_occurrence,
  top_n_per_cluster_region = opt$top_n,
  n_cores = opt$cores,
  
  seed = 42
)

################################################################################
# 路径设置
################################################################################

PCRE_DIR <- file.path(PARAMS$base_dir, "06_pCRE")

# 构建输出目录名称
bg_label <- if (is.null(PARAMS$enrichment_log2fc)) "GenomeWideBG" else sprintf("log2FC%.1f", PARAMS$enrichment_log2fc)
gc_label <- if (PARAMS$use_gc_match_enrichment) "_GCmatched" else ""
mode_label <- ifelse(PARAMS$region_mode == "full_length", "_FullLength", "_Partitioned")

OUTPUT_BASE <- file.path(PCRE_DIR, paste0("Pipeline6a_v8.0_", bg_label, gc_label, mode_label))

setup_output_directories <- function(output_base, regions, genotypes) {
  dir.create(output_base, recursive = TRUE, showWarnings = FALSE)
  
  genotype_dirs <- list()
  
  for (geno in genotypes) {
    geno_dir <- file.path(output_base, geno)
    dir.create(geno_dir, recursive = TRUE, showWarnings = FALSE)
    
    region_dirs <- list()
    for (region_id in names(regions)) {
      region_dir <- file.path(geno_dir, paste0("Region_", region_id))
      for (subdir in c("promoter_sequences", "cluster_results")) {
        dir.create(file.path(region_dir, subdir), recursive = TRUE, showWarnings = FALSE)
      }
      region_dirs[[region_id]] <- region_dir
    }
    
    dir.create(file.path(geno_dir, "ML_ready"), recursive = TRUE, showWarnings = FALSE)
    
    genotype_dirs[[geno]] <- list(
      base = geno_dir,
      regions = region_dirs,
      ml_ready = file.path(geno_dir, "ML_ready")
    )
  }
  
  dir.create(file.path(output_base, "logs"), recursive = TRUE, showWarnings = FALSE)
  
  return(genotype_dirs)
}

GENOTYPE_DIRS <- setup_output_directories(OUTPUT_BASE, PARAMS$genomic_regions, 
                                           PARAMS$genotypes_to_process)

################################################################################
# 日志设置
################################################################################

log_file <- file.path(OUTPUT_BASE, "logs",
                      paste0("pipeline6a_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("================================================================================\n")
cat("  Pipeline 6a v8.0: Unified Multi-Region pCRE Enrichment Analysis\n")
cat("================================================================================\n\n")
cat(sprintf("开始时间: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("输出目录: %s\n\n", OUTPUT_BASE))

cat("================================================================================\n")
cat("  核心参数配置\n")
cat("================================================================================\n\n")

cat("【分析模式】\n")
cat(sprintf("  ✔ 基因型: %s\n", paste(PARAMS$genotypes_to_process, collapse = ", ")))
cat(sprintf("  ✔ 区域模式: %s\n", PARAMS$region_mode))
cat("\n")

cat("【区域配置】\n")
for (region_id in names(PARAMS$genomic_regions)) {
  region <- PARAMS$genomic_regions[[region_id]]
  cat(sprintf("  - %s: %s\n", region_id, region$description))
}
cat("\n")

cat("【富集分析背景】\n")
if (is.null(PARAMS$enrichment_log2fc)) {
  cat("  ✔ 全基因组背景\n")
} else {
  cat(sprintf("  ✔ Log2FC阈值: |log2FC| < %.1f\n", PARAMS$enrichment_log2fc))
}
cat(sprintf("  ✔ GC匹配: %s\n", ifelse(PARAMS$use_gc_match_enrichment, "启用", "禁用")))
cat("\n")

cat("【ML模型背景】\n")
cat(sprintf("  ✔ 非响应基因定义: |log2FC| < %.1f\n", PARAMS$ml_log2fc_threshold))
cat("\n")

cat("【特征选择】\n")
cat(sprintf("  ✔ 每个Cluster/Region选取Top %d (按Raw P-value排序)\n", 
            PARAMS$top_n_per_cluster_region))
cat("  ✔ 注意: 不进行FDR过滤，保留所有候选特征供RF筛选\n")
cat("\n")

set.seed(PARAMS$seed)

################################################################################
# 辅助函数
################################################################################

# 基因ID标准化: 根据用户提供的正则剥离后缀, 留空则原样返回
standardize_gene_ids <- function(gene_ids, strip_pattern = "") {
  if (is.null(strip_pattern) || strip_pattern == "") {
    return(gene_ids)
  }
  return(gsub(strip_pattern, "", gene_ids))
}

calculate_gc_content <- function(sequences) {
  gc_content <- sapply(as.character(sequences), function(seq) {
    seq_upper <- toupper(seq)
    gc_count <- sum(strsplit(seq_upper, "")[[1]] %in% c("G", "C"))
    total <- nchar(seq_upper)
    if (total > 0) gc_count / total else NA
  })
  return(gc_content)
}

gc_matched_sampling <- function(pos_genes, neg_pool, gc_info, n_bins = 10, ratio = 1) {
  pos_gc <- gc_info[pos_genes]
  pos_gc <- pos_gc[!is.na(pos_gc)]
  
  if (length(pos_gc) == 0) return(sample(neg_pool, min(length(neg_pool), length(pos_genes) * ratio)))
  
  gc_breaks <- quantile(gc_info[!is.na(gc_info)], probs = seq(0, 1, length.out = n_bins + 1))
  gc_breaks[1] <- 0
  gc_breaks[length(gc_breaks)] <- 1
  
  pos_bins <- cut(pos_gc, breaks = gc_breaks, labels = FALSE, include.lowest = TRUE)
  pos_bin_counts <- table(pos_bins)
  
  neg_gc <- gc_info[neg_pool]
  neg_gc <- neg_gc[!is.na(neg_gc)]
  neg_bins <- cut(neg_gc, breaks = gc_breaks, labels = FALSE, include.lowest = TRUE)
  
  sampled_neg <- c()
  
  for (bin in names(pos_bin_counts)) {
    n_needed <- as.integer(pos_bin_counts[bin] * ratio)
    neg_in_bin <- names(neg_gc)[neg_bins == as.integer(bin)]
    
    if (length(neg_in_bin) >= n_needed) {
      sampled_neg <- c(sampled_neg, sample(neg_in_bin, n_needed))
    } else if (length(neg_in_bin) > 0) {
      sampled_neg <- c(sampled_neg, neg_in_bin)
    }
  }
  
  return(unique(sampled_neg))
}

################################################################################
# 模块 1: 读取 Membership Matrix (通用化)
################################################################################

load_membership_matrix <- function(membership_file, genotype, threshold = 0.5,
                                   condition_wt = NULL, condition_mut = NULL,
                                   gene_id_strip = "") {
  cat(sprintf("\n--- 读取 %s Membership Matrix ---\n", genotype))
  
  if (!file.exists(membership_file)) {
    stop(sprintf("文件不存在: %s", membership_file))
  }
  
  mem_df <- read.table(membership_file, header = TRUE, sep = "\t", 
                       stringsAsFactors = FALSE, check.names = FALSE)
  
  # 识别基因列
  if ("gene" %in% colnames(mem_df)) {
    gene_col <- "gene"
  } else if ("GeneID" %in% colnames(mem_df)) {
    gene_col <- "GeneID"
  } else {
    gene_col <- colnames(mem_df)[1]
  }
  
  genes <- standardize_gene_ids(mem_df[[gene_col]], gene_id_strip)
  
  # 识别Cluster列: 通用策略
  all_candidate_cols <- setdiff(colnames(mem_df), gene_col)
  
  # 策略1: 根据基因型和条件名构建匹配模式
  if (genotype == "WT") {
    # 尝试匹配: WT_C*, <condition_wt>_C*, 或裸 C[0-9]+
    patterns <- c("^WT_C", "^C[0-9]")
    if (!is.null(condition_wt) && condition_wt != "") {
      # 转义条件名中的特殊正则字符
      esc_cond <- gsub("([\\.|()\\[\\]{}^$*+?\\\\])", "\\\\\\1", condition_wt, perl = TRUE)
      patterns <- c(paste0("^", esc_cond, "_C"), patterns)
    }
  } else {
    patterns <- c("^Mutant_C", "^MUT_C", "^C[0-9]")
    if (!is.null(condition_mut) && condition_mut != "") {
      esc_cond <- gsub("([\\.|()\\[\\]{}^$*+?\\\\])", "\\\\\\1", condition_mut, perl = TRUE)
      patterns <- c(paste0("^", esc_cond, "_C"), patterns)
    }
  }
  
  cluster_cols <- character(0)
  for (pat in patterns) {
    cluster_cols <- all_candidate_cols[grepl(pat, all_candidate_cols, ignore.case = TRUE)]
    if (length(cluster_cols) > 0) break
  }
  
  # 回退: 所有数值列
  if (length(cluster_cols) == 0) {
    cluster_cols <- all_candidate_cols[sapply(all_candidate_cols, function(x) is.numeric(mem_df[[x]]))]
    cat("  ⚠ 未能按模式匹配Cluster列, 回退至全部数值列\n")
  }
  
  # 排序
  cluster_nums <- as.numeric(gsub("[^0-9]", "", cluster_cols))
  if (!any(is.na(cluster_nums))) {
    cluster_cols <- cluster_cols[order(cluster_nums)]
  }
  
  mem_matrix <- as.matrix(mem_df[, cluster_cols])
  rownames(mem_matrix) <- genes
  
  cat(sprintf("  检测到 %d 个 Clusters: %s\n", 
              length(cluster_cols), 
              paste(head(cluster_cols, 5), collapse = ", ")))
  cat(sprintf("  总基因数: %d\n", nrow(mem_matrix)))
  
  return(list(
    matrix = mem_matrix,
    genes = genes,
    cluster_cols = cluster_cols,
    genotype = genotype
  ))
}

################################################################################
# 模块 2: 计算表达变化 (通用化)
################################################################################

calculate_expression_changes <- function(vst_file, responsive_genes, genotype,
                                         condition_wt, condition_mut,
                                         control_time_label, control_time_numeric,
                                         gene_id_strip = "") {
  cat(sprintf("\n--- 计算表达变化 (%s 样本) ---\n", genotype))
  
  vst_data <- read.table(vst_file, header = TRUE, row.names = 1, sep = "\t",
                         check.names = FALSE)
  rownames(vst_data) <- standardize_gene_ids(rownames(vst_data), gene_id_strip)
  
  # ========== 通用样本列筛选 ==========
  # 使用Config传入的条件名, 而非硬编码
  if (genotype == "WT") {
    target_condition <- condition_wt
    exclude_condition <- condition_mut
  } else {
    target_condition <- condition_mut
    exclude_condition <- condition_wt
  }
  
  # 转义用于正则
  esc_target <- gsub("([\\.|()\\[\\]{}^$*+?\\\\])", "\\\\\\1", target_condition, perl = TRUE)
  esc_exclude <- gsub("([\\.|()\\[\\]{}^$*+?\\\\])", "\\\\\\1", exclude_condition, perl = TRUE)
  
  sample_cols <- grep(esc_target, colnames(vst_data), value = TRUE, ignore.case = TRUE)
  # 排除另一基因型 (防止 "WT" 匹配到包含 "MpWT" 的同时也匹配到 "Mpstop1_WT_xxx" 之类)
  if (length(exclude_condition) > 0 && exclude_condition != "") {
    sample_cols <- sample_cols[!grepl(esc_exclude, sample_cols, ignore.case = TRUE)]
  }
  
  if (length(sample_cols) == 0) {
    stop(sprintf("错误: 在VST文件中未找到包含 '%s' 的样本列\n  可用列: %s", 
                 target_condition, paste(head(colnames(vst_data), 10), collapse = ", ")))
  }
  
  cat(sprintf("  找到 %d 个样本列\n", length(sample_cols)))
  
  # ========== 通用时间点提取 (基线锚定策略) ==========
  # 策略: 从用户指定的 control_time_label (如 "0d") 中提取时间单位 ("d")
  # 从而在复杂的样本名 (如 "MpWT.DR0d.1") 中精准提取被该单位修饰的数字
  time_unit <- gsub("[0-9.]", "", control_time_label)
  if (nchar(time_unit) == 0) time_unit <- "[a-zA-Z]*" # 兼容无单位的纯数字情况
  
  extract_timepoint <- function(sample_name) {
    # 构造正则捕获组: 匹配 (数字) + 单位 + 非字母分隔符/结尾
    # 例如匹配 "DR0d.1" 中的 "0", "12h_1" 中的 "12"
    pattern <- paste0("([0-9]+(?:\\.[0-9]+)?)", time_unit, "(?:[^a-zA-Z]|$)")
    
    m <- regexec(pattern, sample_name)
    matches <- regmatches(sample_name, m)
    
    # 如果正则捕获成功，matches[[1]][2] 就是纯数字部分
    if (length(matches[[1]]) >= 2) {
      return(as.numeric(matches[[1]][2]))
    }
    
    # 终极回退策略: 暴力拆分，提取被 _ 或 . 隔开的片段中包含数字的部分
    parts <- unlist(strsplit(sample_name, "[_.]"))
    for (part in parts) {
      if (grepl("[0-9]", part)) {
        num <- as.numeric(gsub("[^0-9.]", "", part))
        if (!is.na(num)) return(num)
      }
    }
    return(NA)
  }
  
  sample_timepoints <- sapply(sample_cols, extract_timepoint)
  time_points <- sort(unique(sample_timepoints[!is.na(sample_timepoints)]))
  
  if (length(time_points) == 0) {
    stop(sprintf("错误: 无法从样本名中解析时间点。样本示例: %s\n  期望格式: <Condition>_<数字><单位>_<rep>", 
                 paste(head(sample_cols, 3), collapse = ", ")))
  }
  
  cat(sprintf("  解析到时间点: %s\n", paste(time_points, collapse = ", ")))
  
  # ========== 通用T0识别 ==========
  # 优先使用 control_time_numeric (从 Config 中的 CONTROL_TIME 解析而来)
  t0_cols <- sample_cols[sample_timepoints == control_time_numeric]
  
  if (length(t0_cols) == 0) {
    # 回退: 使用最小时间点
    min_tp <- min(sample_timepoints, na.rm = TRUE)
    t0_cols <- sample_cols[sample_timepoints == min_tp]
    cat(sprintf("  ⚠ 未找到时间=%g的样本, 使用最小时间点 %g 作为基线 (%d 列)\n", 
                control_time_numeric, min_tp, length(t0_cols)))
  } else {
    cat(sprintf("  基线T0: %d 列 (时间=%g)\n", length(t0_cols), control_time_numeric))
  }
  
  t0_mean <- rowMeans(vst_data[, t0_cols, drop = FALSE], na.rm = TRUE)
  
  # 计算最大|log2FC|
  max_abs_log2fc <- rep(0, nrow(vst_data))
  names(max_abs_log2fc) <- rownames(vst_data)
  
  for (tp in time_points) {
    if (tp == control_time_numeric) next
    tp_cols <- sample_cols[sample_timepoints == tp]
    if (length(tp_cols) == 0) next
    tp_mean <- rowMeans(vst_data[, tp_cols, drop = FALSE], na.rm = TRUE)
    log2fc <- tp_mean - t0_mean
    max_abs_log2fc <- pmax(max_abs_log2fc, abs(log2fc), na.rm = TRUE)
  }
  
  mean_expr <- rowMeans(vst_data[, sample_cols], na.rm = TRUE)
  all_vst_genes <- rownames(vst_data)
  
  cat(sprintf("  VST矩阵基因数: %d\n", length(all_vst_genes)))
  cat(sprintf("  响应基因数: %d\n", length(responsive_genes)))
  
  return(list(
    all_genes = all_vst_genes,
    max_log2fc = max_abs_log2fc,
    mean_expr = mean_expr,
    responsive_genes = responsive_genes
  ))
}
get_nonresponsive_genes <- function(expr_data, log2fc_threshold, min_expr = 1) {
  if (is.null(log2fc_threshold)) return(expr_data$all_genes)
  is_stable <- expr_data$max_log2fc < log2fc_threshold & expr_data$mean_expr >= min_expr
  stable_genes <- names(expr_data$max_log2fc)[is_stable]
  stable_genes <- setdiff(stable_genes, expr_data$responsive_genes)
  return(stable_genes[!is.na(stable_genes)])
}

################################################################################
# 模块 3: 序列提取
################################################################################

extract_region_sequences_with_gc <- function(gene_ids, genome_fasta, gff_file, 
                                              region_config, output_file = NULL,
                                              gene_id_strip = "") {
  cat(sprintf("\n--- 提取区域序列: %s ---\n", region_config$name))
  
  fai_file <- paste0(genome_fasta, ".fai")
  if (!file.exists(fai_file)) {
    Rsamtools::indexFa(genome_fasta)
  }
  
  fa_file <- Rsamtools::FaFile(genome_fasta)
  open(fa_file)
  fa_seqinfo <- seqinfo(fa_file)
  
  txdb <- tryCatch({
    suppressMessages(makeTxDbFromGFF(gff_file, format = "gff3"))
  }, error = function(e) {
    gff <- rtracklayer::import(gff_file)
    genes_gr <- gff[gff$type == "gene"]
    names(genes_gr) <- genes_gr$ID
    return(genes_gr)
  })
  
  if (inherits(txdb, "TxDb")) {
    genes_gr <- genes(txdb)
  } else {
    genes_gr <- txdb
  }
  
  # 通用化基因ID匹配
  ref_ids <- names(genes_gr)
  ref_short <- standardize_gene_ids(ref_ids, gene_id_strip)
  short_to_full <- setNames(ref_ids, ref_short)
  
  query_short <- standardize_gene_ids(gene_ids, gene_id_strip)
  matched_full_ids <- short_to_full[query_short]
  valid_idx <- !is.na(matched_full_ids)
  gene_ids_matched <- matched_full_ids[valid_idx]
  
  if (length(gene_ids_matched) == 0) {
    # 尝试直接匹配 (无需标准化的情况)
    direct_match <- intersect(gene_ids, ref_ids)
    if (length(direct_match) > 0) {
      cat(sprintf("  ⚠ 标准化匹配失败, 直接匹配到 %d 基因\n", length(direct_match)))
      gene_ids_matched <- setNames(direct_match, direct_match)
    } else {
      cat("  ❌ 无法匹配任何基因ID\n")
      cat(sprintf("    查询示例: %s\n", paste(head(gene_ids, 3), collapse = ", ")))
      cat(sprintf("    参考示例: %s\n", paste(head(ref_ids, 3), collapse = ", ")))
      close(fa_file)
      return(NULL)
    }
  }
  
  genes_subset <- genes_gr[gene_ids_matched]
  
  # 构建区域 (不变)
  if (!is.null(region_config$trim_start)) {
    full_upstream_gr <- promoters(genes_subset, upstream = region_config$upstream, downstream = 0)
    distal_length <- region_config$upstream - region_config$trim_start
    region_gr <- resize(full_upstream_gr, width = distal_length, fix = "start")
  } else {
    region_gr <- promoters(genes_subset, upstream = region_config$upstream, 
                           downstream = region_config$downstream)
  }
  
  genome_chroms <- seqnames(fa_seqinfo)
  seqlengths_genome <- seqlengths(fa_seqinfo)
  
  common_chroms <- intersect(as.character(seqnames(region_gr)), genome_chroms)
  region_gr <- region_gr[as.character(seqnames(region_gr)) %in% common_chroms]
  seqlevels(region_gr) <- common_chroms
  seqlengths(region_gr) <- seqlengths_genome[common_chroms]
  region_gr <- trim(region_gr)
  
  if (!is.null(region_config$trim_start)) {
    expected_length <- region_config$upstream - region_config$trim_start
  } else {
    expected_length <- region_config$upstream + region_config$downstream
  }
  min_length <- expected_length * 0.5
  region_gr <- region_gr[width(region_gr) >= min_length]
  
  region_seqs <- getSeq(fa_file, region_gr)
  close(fa_file)
  
  names(region_seqs) <- standardize_gene_ids(names(region_gr), gene_id_strip)
  gc_content <- calculate_gc_content(region_seqs)
  
  if (!is.null(output_file)) {
    writeXStringSet(region_seqs, output_file)
  }
  
  cat(sprintf("  ✔ 有效序列: %d, 平均GC: %.1f%%\n", 
              length(region_seqs), mean(gc_content, na.rm = TRUE) * 100))
  
  return(list(
    sequences = region_seqs,
    gc_content = gc_content
  ))
}

################################################################################
# 模块 4: 二值化k-mer特征矩阵
################################################################################

generate_binary_kmer_matrix <- function(sequences, k_values = c(6), n_cores = 40) {
  cat("\n--- 生成二值化 k-mer 特征矩阵 ---\n")
  
  all_kmers <- unlist(lapply(k_values, function(k) {
    as.character(mkAllStrings(c("A", "C", "G", "T"), k))
  }))
  
  seq_list <- as.character(sequences)
  seq_names <- names(sequences)
  n_seqs <- length(seq_list)
  
  cl <- parallel::makeForkCluster(min(n_cores, detectCores() - 1))
  registerDoParallel(cl)
  
  chunk_size <- ceiling(n_seqs / n_cores)
  chunks <- split(1:n_seqs, ceiling(seq_along(1:n_seqs) / chunk_size))
  
  binary_matrix <- foreach(
    chunk_idx = seq_along(chunks),
    .combine = rbind,
    .errorhandling = "pass"
  ) %dopar% {
    gene_indices <- chunks[[chunk_idx]]
    chunk_mat <- matrix(0L, nrow = length(gene_indices), ncol = length(all_kmers))
    
    for (i in seq_along(gene_indices)) {
      idx <- gene_indices[i]
      seq_str <- seq_list[[idx]]
      
      for (j in seq_along(all_kmers)) {
        chunk_mat[i, j] <- as.integer(grepl(all_kmers[j], seq_str, fixed = TRUE))
      }
    }
    
    rownames(chunk_mat) <- seq_names[gene_indices]
    chunk_mat
  }
  
  stopCluster(cl)
  colnames(binary_matrix) <- all_kmers
  
  cat(sprintf("  ✔ 特征矩阵: %d x %d\n", nrow(binary_matrix), ncol(binary_matrix)))
  
  return(binary_matrix)
}

################################################################################
# 模块 5: Fisher's Exact Test 富集分析
################################################################################

run_fisher_enrichment <- function(pos_genes, neg_genes, kmer_matrix, 
                                  min_occurrence = 5, cluster_name = "Cluster") {
  
  pos_genes <- intersect(pos_genes, rownames(kmer_matrix))
  neg_genes <- intersect(neg_genes, rownames(kmer_matrix))
  
  n_pos <- length(pos_genes)
  n_neg <- length(neg_genes)
  
  if (n_pos == 0 || n_neg == 0) return(NULL)
  
  pos_matrix <- kmer_matrix[pos_genes, , drop = FALSE]
  neg_matrix <- kmer_matrix[neg_genes, , drop = FALSE]
  
  pos_counts <- colSums(pos_matrix)
  neg_counts <- colSums(neg_matrix)
  total_counts <- pos_counts + neg_counts
  
  keep_idx <- total_counts >= min_occurrence
  filtered_kmers <- colnames(kmer_matrix)[keep_idx]
  
  if (length(filtered_kmers) == 0) return(NULL)
  
  pos_counts <- pos_counts[keep_idx]
  neg_counts <- neg_counts[keep_idx]
  
  results <- data.frame(
    Cluster = cluster_name,
    kmer = filtered_kmers,
    k = nchar(filtered_kmers),
    pos_count = as.integer(pos_counts),
    neg_count = as.integer(neg_counts),
    n_pos_total = n_pos,
    n_neg_total = n_neg,
    pos_ratio = round(pos_counts / n_pos, 4),
    neg_ratio = round(neg_counts / n_neg, 4),
    stringsAsFactors = FALSE
  )
  
  results$fold_enrichment <- with(results, {
    fe <- ifelse(neg_ratio > 0, pos_ratio / neg_ratio, ifelse(pos_ratio > 0, 10, 1))
    pmin(fe, 100)
  })
  
  results$p_value <- sapply(1:nrow(results), function(i) {
    a <- results$pos_count[i]
    b <- n_pos - a
    c <- results$neg_count[i]
    d <- n_neg - c
    contingency <- matrix(c(a, c, b, d), nrow = 2)
    tryCatch({
      fisher.test(contingency, alternative = "greater")$p.value
    }, error = function(e) 1.0)
  })
  
  results <- results %>% 
    arrange(p_value) %>%
    mutate(rank = row_number())
  
  return(results)
}

################################################################################
# 模块 6: 单区域分析
################################################################################

run_single_region_analysis <- function(region_id, region_config, membership_data,
                                       expr_data, all_genes, output_dir, params) {
  cat("\n")
  cat("################################################################################\n")
  cat(sprintf("  区域分析: %s (%s)\n", region_id, membership_data$genotype))
  cat("################################################################################\n")
  
  # 1. 提取序列和GC
  seq_file <- file.path(output_dir, "promoter_sequences", 
                        sprintf("%s_sequences.fa", region_config$name))
  gc_file <- file.path(output_dir, "promoter_sequences", 
                       sprintf("%s_gc_content.rds", region_config$name))
  
  if (file.exists(seq_file) && file.exists(gc_file)) {
    region_seqs <- readDNAStringSet(seq_file)
    gc_content <- readRDS(gc_file)
    cat("  (使用缓存序列)\n")
  } else {
    result <- extract_region_sequences_with_gc(
      all_genes, params$genome_fasta, params$gff_file,
      region_config, seq_file,
      gene_id_strip = params$gene_id_strip    # ← 新增
    )
    if (is.null(result)) return(NULL)
    region_seqs <- result$sequences
    gc_content <- result$gc_content
    saveRDS(gc_content, gc_file)
  }
  
  # 2. 生成k-mer矩阵
  kmer_matrix_file <- file.path(output_dir, 
                                 sprintf("kmer_matrix_%s.rds", region_config$name))
  
  if (file.exists(kmer_matrix_file)) {
    kmer_matrix <- readRDS(kmer_matrix_file)
    cat("  (使用缓存k-mer矩阵)\n")
  } else {
    kmer_matrix <- generate_binary_kmer_matrix(region_seqs, params$k_values, params$n_cores)
    saveRDS(kmer_matrix, kmer_matrix_file)
  }
  
  # 3. 富集分析
  cat("\n--- Cluster-Specific Fisher's Exact Test ---\n")
  
  all_enrichment_results <- data.frame()
  top_n <- params$top_n_per_cluster_region
  
  for (cluster_name in membership_data$cluster_cols) {
    
    core_genes <- membership_data$genes[
      membership_data$matrix[, cluster_name] >= params$membership_threshold
    ]
    
    if (length(core_genes) < 10) {
      cat(sprintf("  %s: 跳过 (仅%d基因)\n", cluster_name, length(core_genes)))
      next
    }
    
    # 确定富集分析的背景基因
    if (is.null(params$enrichment_log2fc)) {
      bg_pool <- setdiff(names(region_seqs), core_genes)
    } else {
      bg_pool <- get_nonresponsive_genes(expr_data, params$enrichment_log2fc, 
                                          params$min_expression)
      bg_pool <- intersect(bg_pool, names(region_seqs))
    }
    
    # GC匹配 (可选)
    if (params$use_gc_match_enrichment && length(bg_pool) > length(core_genes) * 5) {
      background_genes <- gc_matched_sampling(core_genes, bg_pool, gc_content, 
                                               n_bins = 10, ratio = 5)
    } else {
      background_genes <- bg_pool
    }
    
    enrichment <- run_fisher_enrichment(
      pos_genes = core_genes,
      neg_genes = background_genes,
      kmer_matrix = kmer_matrix,
      min_occurrence = params$min_occurrence,
      cluster_name = cluster_name
    )
    
    if (!is.null(enrichment) && nrow(enrichment) > 0) {
      enrichment$Region <- region_id
      enrichment$Genotype <- membership_data$genotype
      all_enrichment_results <- rbind(all_enrichment_results, enrichment)
      
      write.table(enrichment,
                  file.path(output_dir, "cluster_results", 
                            paste0(cluster_name, "_Enrichment.txt")),
                  sep = "\t", quote = FALSE, row.names = FALSE)
      
      cat(sprintf("  %s: %d genes, Top%d selected\n", 
                  cluster_name, length(core_genes), min(top_n, nrow(enrichment))))
    }
  }
  
  if (nrow(all_enrichment_results) > 0) {
    write.table(all_enrichment_results,
                file.path(output_dir, "All_Clusters_Enrichment.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
    
    return(list(
      enrichment = all_enrichment_results,
      kmer_matrix = kmer_matrix,
      gc_content = gc_content,
      region_id = region_id
    ))
  }
  
  return(NULL)
}

################################################################################
# 模块 7: 生成ML_ready数据
################################################################################

generate_ml_ready_data <- function(all_region_results, membership_data, 
                                    expr_data, params, ml_ready_dir, genotype) {
  cat("\n================================================================================\n")
  cat(sprintf("  生成ML_ready数据 (%s)\n", genotype))
  cat("================================================================================\n\n")
  
  top_n <- params$top_n_per_cluster_region
  
  # 1. 合并富集结果
  combined_enrichment <- do.call(rbind, lapply(all_region_results, function(x) x$enrichment))
  cat(sprintf("  合并富集结果: %d 行\n", nrow(combined_enrichment)))
  
  # 2. 生成白名单 (基于Raw P-value排序，不做FDR过滤)
  whitelist_full <- combined_enrichment %>%
    group_by(Cluster, Region) %>%
    arrange(p_value) %>%
    slice_head(n = top_n) %>%
    ungroup() %>%
    mutate(feature_id = paste(Region, kmer, sep = "_"))
  
  whitelist_dedup <- whitelist_full %>%
    group_by(feature_id) %>%
    summarise(
      kmer = dplyr::first(kmer),
      Region = dplyr::first(Region),
      N_Sources = n(),
      Source_Clusters = paste(unique(Cluster), collapse = ";"),
      Best_FE = max(fold_enrichment, na.rm = TRUE),
      Best_pval = min(p_value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Best_pval)
  
  cat(sprintf("  白名单特征: %d (来自%d个Cluster-Region组合)\n", 
              nrow(whitelist_dedup), nrow(whitelist_full)))
  
  # 3. 合并k-mer矩阵和GC信息
  common_genes <- Reduce(intersect, lapply(all_region_results, function(x) rownames(x$kmer_matrix)))
  
  region_matrices <- list()
  combined_gc <- list()
  
  for (region_result in all_region_results) {
    region_id <- region_result$region_id
    region_kmers <- whitelist_full %>%
      filter(Region == region_id) %>%
      pull(kmer) %>%
      unique()
    region_kmers <- intersect(region_kmers, colnames(region_result$kmer_matrix))
    
    if (length(region_kmers) > 0) {
      sub_matrix <- region_result$kmer_matrix[common_genes, region_kmers, drop = FALSE]
      colnames(sub_matrix) <- paste(region_id, colnames(sub_matrix), sep = "_")
      region_matrices[[region_id]] <- sub_matrix
    }
    
    combined_gc[[region_id]] <- region_result$gc_content[common_genes]
  }
  
  combined_matrix <- do.call(cbind, region_matrices)
  
  gc_df <- as.data.frame(combined_gc)
  avg_gc <- rowMeans(gc_df, na.rm = TRUE)
  names(avg_gc) <- common_genes
  
  cat(sprintf("  合并矩阵: %d genes x %d features\n", 
              nrow(combined_matrix), ncol(combined_matrix)))
  
  # 4. 为ML准备独立的非响应基因背景
  ml_nonresponsive <- get_nonresponsive_genes(expr_data, params$ml_log2fc_threshold, 
                                               params$min_expression)
  ml_nonresponsive <- intersect(ml_nonresponsive, common_genes)
  
  cat(sprintf("  ML非响应基因: %d (|log2FC| < %.1f)\n", 
              length(ml_nonresponsive), params$ml_log2fc_threshold))
  
  # 5. 准备Cluster数据
  cluster_data <- list()
  
  for (cluster_name in membership_data$cluster_cols) {
    core_genes <- membership_data$genes[
      membership_data$matrix[, cluster_name] >= params$membership_threshold
    ]
    core_genes <- intersect(core_genes, common_genes)
    
    if (length(core_genes) >= 10) {
      cluster_data[[cluster_name]] <- list(
        cluster = cluster_name,
        pos_genes = core_genes,
        n_pos = length(core_genes)
      )
    }
  }
  
  cat(sprintf("  有效Cluster: %d\n", length(cluster_data)))
  
  # 6. 保存
  write.table(whitelist_full %>%
                dplyr::select(feature_id, kmer, Region, Cluster, 
                              fold_enrichment, p_value, rank),
              file.path(ml_ready_dir, "Features_Full_Details.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  write.table(whitelist_dedup,
              file.path(ml_ready_dir, "Features_Whitelist.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  saveRDS(combined_matrix, file.path(ml_ready_dir, "Combined_Kmer_Matrix.rds"))
  
  ml_input_data <- list(
    genotype = genotype,
    cluster_data = cluster_data,
    kmer_matrix = combined_matrix,
    whitelist = whitelist_dedup$feature_id,
    whitelist_details = whitelist_dedup,
    
    ml_nonresponsive_genes = ml_nonresponsive,
    ml_log2fc_threshold = params$ml_log2fc_threshold,
    
    gc_content = avg_gc,
    gc_by_region = combined_gc,
    max_log2fc = expr_data$max_log2fc[common_genes],
    
    common_genes = common_genes,
    params = list(
      enrichment_background = ifelse(is.null(params$enrichment_log2fc), 
                                     "genome_wide", 
                                     sprintf("log2FC<%.1f", params$enrichment_log2fc)),
      gc_match_enrichment = params$use_gc_match_enrichment,
      ml_log2fc_threshold = params$ml_log2fc_threshold,
      k_values = params$k_values,
      top_n = params$top_n_per_cluster_region,
      regions = names(params$genomic_regions),
      region_mode = params$region_mode
    )
  )
  
  saveRDS(ml_input_data, file.path(ml_ready_dir, "Pipeline6b_Input_Data.rds"))
  
  cat("\n  ✔ 输出完成:\n")
  cat("    - Features_Full_Details.tsv\n")
  cat("    - Features_Whitelist.txt\n")
  cat("    - Combined_Kmer_Matrix.rds\n")
  cat("    - Pipeline6b_Input_Data.rds\n")
  
  return(ml_input_data)
}

################################################################################
# 模块 8: 处理单个基因型 (通用化)
################################################################################

process_genotype <- function(genotype, genotype_dirs, params) {
  cat("\n")
  cat("================================================================================\n")
  cat(sprintf("  处理基因型: %s\n", genotype))
  cat("================================================================================\n")
  
  # 1. 读取Membership Matrix — 传入条件名
  membership_file <- if (genotype == "WT") params$wt_membership_file else params$mutant_membership_file
  membership_data <- load_membership_matrix(
    membership_file, genotype, params$membership_threshold,
    condition_wt = params$condition_wt,
    condition_mut = params$condition_mut,
    gene_id_strip = params$gene_id_strip
  )
  
  all_responsive_genes <- membership_data$genes[
    apply(membership_data$matrix, 1, max) >= params$membership_threshold
  ]
  
  # 2. 计算表达变化 — 传入条件名和时间基线
  expr_data <- calculate_expression_changes(
    params$vst_file, all_responsive_genes, genotype,
    condition_wt = params$condition_wt,
    condition_mut = params$condition_mut,
    control_time_label = params$control_time_label,
    control_time_numeric = params$control_time_numeric,
    gene_id_strip = params$gene_id_strip
  )
  all_genes <- expr_data$all_genes
  
  # 3. 多区域分析 (不变)
  all_region_results <- list()
  
  for (region_id in names(params$genomic_regions)) {
    region_config <- params$genomic_regions[[region_id]]
    output_dir <- genotype_dirs$regions[[region_id]]
    
    result <- run_single_region_analysis(
      region_id = region_id,
      region_config = region_config,
      membership_data = membership_data,
      expr_data = expr_data,
      all_genes = all_genes,
      output_dir = output_dir,
      params = params
    )
    
    if (!is.null(result)) {
      all_region_results[[region_id]] <- result
    }
  }
  
  # 4. 生成ML_ready数据 (不变)
  ml_data <- NULL
  if (length(all_region_results) > 0) {
    ml_data <- generate_ml_ready_data(
      all_region_results, membership_data, expr_data, params,
      genotype_dirs$ml_ready, genotype
    )
  }
  
  return(list(
    genotype = genotype,
    membership_data = membership_data,
    expr_data = expr_data,
    region_results = all_region_results,
    ml_data = ml_data
  ))
}

################################################################################
# 主程序
################################################################################

cat("\n================================================================================\n")
cat("  开始主分析\n")
cat("================================================================================\n")

all_genotype_results <- list()

for (genotype in PARAMS$genotypes_to_process) {
  result <- process_genotype(genotype, GENOTYPE_DIRS[[genotype]], PARAMS)
  all_genotype_results[[genotype]] <- result
}

# 最终汇总
cat("\n================================================================================\n")
cat("  Pipeline 6a v8.0 完成\n")
cat("================================================================================\n")
cat(sprintf("结束时间: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

cat("【处理汇总】\n")
for (genotype in names(all_genotype_results)) {
  result <- all_genotype_results[[genotype]]
  n_clusters <- length(result$membership_data$cluster_cols)
  n_regions <- length(result$region_results)
  n_features <- if (!is.null(result$ml_data)) ncol(result$ml_data$kmer_matrix) else 0
  
  cat(sprintf("  %s: %d Clusters, %d Regions, %d Features\n", 
              genotype, n_clusters, n_regions, n_features))
}

cat(sprintf("\n输出目录: %s\n", OUTPUT_BASE))

sink(type = "message")
sink()
close(log_con)

cat("\n分析完成。\n")