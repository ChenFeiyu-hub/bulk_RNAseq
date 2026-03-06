#!/usr/bin/env Rscript
################################################################################
# RNAseq Pipeline 6c v3.1: Unified Motif Integration & Annotation
#
# 本脚本统一替代旧版 6c (Levenshtein聚类 + IUPAC简并) 与 6d (TOMTOM + 数据整合),
# 采用全 R 原生 universalmotif 方案完成:
#   1. 动态加载并合并 6a (Enrichment) 和 6b (RF Importance) 数据
#   2. Levenshtein 层次聚类定义 Motif Family
#   3. Shift-Aware Weighted Consensus PWM (universalmotif::merge_motifs)
#   4. 内部 PCC 注释 (vs 用户提供的 MEME DB)
#
# 关键设计原则:
#   - 6b 是 Master List — 保留 6b 的全部行, 不做任何过滤
#   - 6a Enrichment 数据以 left_join 方式补充到 6b 上
#   - 移除旧版拐点/阈值过滤逻辑
#   - 用 universalmotif PCC 替代外部 TOMTOM 依赖
#   - 不对 PCC 进行阈值过滤, 直接取 Top 2 注释
#
# 输入:
#   --p6a_dir : Pipeline 6a 输出根目录 (含 {Genotype}/Region_* 子目录)
#   --p6b_dir : Pipeline 6b 输出根目录 (含 {Genotype}/importance/ 子目录)
#   --meme_db : MEME 格式 motif 数据库 (必需)
#
# 输出:
#   - Integrated_Features_Full.txt   : 全表 (行 = Raw K-mer)
#   - Consensus_Motif_Summary.txt    : 汇总 (行 = Motif Family)
#   - Motif_PWMs.rds                 : universalmotif PWM 列表
#
# 版本: 3.1 (修复注释模块接口 + 移除PCC阈值过滤)
# 日期: 2026-02
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(stringdist)
  library(universalmotif)
  library(parallel)
})

################################################################################
# CLI 参数
################################################################################

option_list <- list(
  make_option(c("--p6a_dir"), type = "character", default = NULL,
              help = "Pipeline 6a 输出目录 (含 WT/Mutant 子目录及 Region_* 文件夹)"),
  make_option(c("--p6b_dir"), type = "character", default = NULL,
              help = "Pipeline 6b 输出目录 (含 WT/Mutant/importance/ 子目录)"),
  make_option(c("--output_dir"), type = "character", default = NULL,
              help = "输出目录"),
  make_option(c("--meme_db"), type = "character", default = NULL,
              help = "MEME 格式的 motif 参考数据库文件 (REQUIRED)"),
  make_option(c("--n_cores"), type = "integer", default = 20,
              help = "并行核心数 [default: %default]"),
  make_option(c("--genotype"), type = "character", default = "both",
              help = "处理的基因型: 'WT', 'Mutant', 或 'both' [default: %default]"),
  make_option(c("--cluster_distance"), type = "integer", default = 2,
              help = "Levenshtein 聚类切割高度 [default: %default]"),
  make_option(c("--min_overlap"), type = "integer", default = 4,
              help = "merge_motifs 最小 overlap [default: %default]"),
  make_option(c("--trim_ic"), type = "double", default = 0.25,
              help = "trim_motifs 最低 IC 阈值 [default: %default]"),
  make_option(c("--seed"), type = "integer", default = 42,
              help = "随机种子 [default: %default]")
)

opt_parser <- OptionParser(
  option_list = option_list,
  description = "Pipeline 6c v3.1: Unified Motif Integration & Annotation (universalmotif)"
)
opt <- parse_args(opt_parser)

# --- 参数验证 ---
if (is.null(opt$p6b_dir))  stop("Error: --p6b_dir is required (Pipeline 6b output)")
if (is.null(opt$meme_db))  stop("Error: --meme_db is required (MEME format motif database)")
if (!file.exists(opt$meme_db)) stop(sprintf("Error: MEME database not found: %s", opt$meme_db))

if (is.null(opt$output_dir)) {
  opt$output_dir <- file.path(dirname(opt$p6b_dir), "Pipeline6c_v3.1_MotifIntegration")
}

genotypes_to_process <- if (tolower(opt$genotype) == "both") {
  c("WT", "Mutant")
} else {
  opt$genotype
}

set.seed(opt$seed)

################################################################################
# 创建输出目录
################################################################################

for (geno in genotypes_to_process) {
  dir.create(file.path(opt$output_dir, geno), recursive = TRUE, showWarnings = FALSE)
}
dir.create(file.path(opt$output_dir, "logs"), recursive = TRUE, showWarnings = FALSE)

################################################################################
# 日志设置
################################################################################

log_file <- file.path(opt$output_dir, "logs",
                      paste0("pipeline6c_v3_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("================================================================================\n")
cat("  Pipeline 6c v3.1: Unified Motif Integration & Annotation\n")
cat("================================================================================\n\n")
cat(sprintf("开始时间  : %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("6a 目录   : %s\n", ifelse(is.null(opt$p6a_dir), "(未提供)", opt$p6a_dir)))
cat(sprintf("6b 目录   : %s\n", opt$p6b_dir))
cat(sprintf("MEME DB   : %s\n", opt$meme_db))
cat(sprintf("输出目录  : %s\n", opt$output_dir))
cat(sprintf("基因型    : %s\n", paste(genotypes_to_process, collapse = ", ")))
cat(sprintf("并行核心  : %d\n", opt$n_cores))
cat(sprintf("聚类距离  : Levenshtein h = %d\n", opt$cluster_distance))
cat(sprintf("合并参数  : min.overlap = %d, tryRC = TRUE, method = PCC\n", opt$min_overlap))
cat(sprintf("修剪 IC   : %.2f\n", opt$trim_ic))
cat(sprintf("PCC 过滤  : 无 (直接取 Top 2)\n"))
cat("\n")

################################################################################
# 辅助函数
################################################################################

#' 安全读取制表符分隔文件
safe_read_tsv <- function(path, ...) {
  tryCatch(
    read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
               check.names = FALSE, quote = "", fill = TRUE, ...),
    error = function(e) {
      warning(sprintf("  ⚠ 读取失败: %s — %s", path, e$message))
      NULL
    }
  )
}

#' 将单条 k-mer 序列转为 universalmotif PPM 对象
#' @param kmer   纯 ACGT 序列字符串
#' @param name   motif 名称
#' @param nsites 加权站点数 (用于 merge 时的权重体现)
kmer_to_motif <- function(kmer, name = kmer, nsites = 100) {
  bases <- strsplit(toupper(kmer), "")[[1]]
  len   <- length(bases)
  mat   <- matrix(0, nrow = 4, ncol = len,
                  dimnames = list(c("A", "C", "G", "T"), NULL))
  for (j in seq_len(len)) {
    b <- bases[j]
    if (b %in% c("A", "C", "G", "T")) {
      mat[b, j] <- 1.0
    } else {
      mat[, j] <- 0.25  # ambiguous -> uniform
    }
  }
  create_motif(mat, name = name, type = "PPM", alphabet = "DNA",
               nsites = as.integer(max(1, round(nsites))))
}

#' 从 universalmotif 对象提取 IUPAC 一致性序列
#' (修复: 确保永远返回长度为 1 的字符向量)
safe_consensus <- function(motif_obj) {
  tryCatch({
    res <- NA_character_
    if (!is.null(motif_obj@consensus)) {
      res <- motif_obj@consensus
    } else {
      # Fallback: build from PPM
      ppm <- convert_type(motif_obj, "PPM")@motif
      bases <- c("A", "C", "G", "T")
      res <- paste0(apply(ppm, 2, function(col) {
        idx <- which.max(col)
        if (col[idx] >= 0.5) bases[idx] else "N"
      }), collapse = "")
    }
    
    # 最终防线：确保由 vector 变为 scalar
    if (length(res) > 1) res <- res[1]
    if (length(res) == 0) res <- NA_character_
    return(res)
    
  }, error = function(e) NA_character_)
}

################################################################################
# 模块 1: 加载 6b Importance 数据 (Master List)
################################################################################

load_6b_data <- function(p6b_dir, genotype) {
  cat(sprintf("\n--- [6b] 加载 %s Importance 数据 ---\n", genotype))

  imp_dir <- file.path(p6b_dir, genotype, "importance")
  if (!dir.exists(imp_dir)) {
    cat(sprintf("  ⚠ 目录不存在: %s\n", imp_dir))
    return(NULL)
  }

  imp_files <- list.files(imp_dir, pattern = "_importance\\.txt$", full.names = TRUE)
  if (length(imp_files) == 0) {
    cat(sprintf("  ⚠ 未找到 importance 文件\n"))
    return(NULL)
  }

  cat(sprintf("  找到 %d 个 importance 文件\n", length(imp_files)))

  all_dfs <- list()
  for (f in imp_files) {
    cluster_name <- gsub("_importance\\.txt$", "", basename(f))
    df <- safe_read_tsv(f)
    if (is.null(df) || nrow(df) == 0) next

    # 确保必需列
    if (!"feature" %in% names(df) || !"mean_importance" %in% names(df)) {
      warning(sprintf("    跳过 %s: 缺少 feature/mean_importance 列", basename(f)))
      next
    }

    # 解析 Region 和 kmer (如果文件中未提供)
    if (!"Region" %in% names(df) || !"kmer" %in% names(df)) {
      parsed <- str_match(df$feature, "^([^_]+)_([ACGT]+)$")
      df$Region <- ifelse(is.na(parsed[, 2]), NA_character_, parsed[, 2])
      df$kmer   <- ifelse(is.na(parsed[, 3]), df$feature, parsed[, 3])
    }

    df$Source_Cluster <- cluster_name
    df$Genotype <- genotype
    all_dfs[[cluster_name]] <- df
    cat(sprintf("    %s: %d features\n", cluster_name, nrow(df)))
  }

  if (length(all_dfs) == 0) return(NULL)

  combined <- bind_rows(all_dfs)

  # 仅保留有效 k-mer (纯 ACGT)
  valid <- grepl("^[ACGT]+$", combined$kmer)
  if (sum(!valid) > 0) {
    cat(sprintf("  ⚠ 移除 %d 个非 ACGT k-mer\n", sum(!valid)))
    combined <- combined[valid, ]
  }

  cat(sprintf("  ✔ 6b 总计: %d clusters, %d features\n",
              length(unique(combined$Source_Cluster)), nrow(combined)))
  return(combined)
}

################################################################################
# 模块 2: 加载 6a Enrichment 数据
################################################################################

load_6a_data <- function(p6a_dir, genotype) {
  cat(sprintf("\n--- [6a] 加载 %s Enrichment 数据 ---\n", genotype))

  if (is.null(p6a_dir) || !dir.exists(p6a_dir)) {
    cat("  ⚠ 6a 目录未提供或不存在, 跳过 Enrichment 整合\n")
    return(NULL)
  }

  # 在 6a 目录下搜索基因型子目录
  geno_dir <- file.path(p6a_dir, genotype)
  if (!dir.exists(geno_dir)) {
    # 模糊匹配
    cands <- list.dirs(p6a_dir, recursive = FALSE)
    hit <- cands[grepl(genotype, basename(cands), ignore.case = TRUE)]
    if (length(hit) > 0) {
      geno_dir <- hit[1]
    } else {
      cat(sprintf("  ⚠ 未找到 %s 基因型目录\n", genotype))
      return(NULL)
    }
  }

  # 递归搜索所有 Region_* 目录下的 All_Clusters_Enrichment.txt
  region_dirs <- list.dirs(geno_dir, recursive = TRUE)
  region_dirs <- region_dirs[grepl("Region_", basename(region_dirs))]

  if (length(region_dirs) == 0) {
    cat("  ⚠ 未找到任何 Region_* 目录\n")
    return(NULL)
  }

  cat(sprintf("  找到 %d 个 Region 目录: %s\n",
              length(region_dirs),
              paste(basename(region_dirs), collapse = ", ")))

  all_enrich <- list()
  for (rd in region_dirs) {
    target <- file.path(rd, "All_Clusters_Enrichment.txt")
    if (!file.exists(target)) {
      cat(sprintf("    ⚠ 缺失: %s/All_Clusters_Enrichment.txt\n", basename(rd)))
      next
    }
    df <- tryCatch(
      safe_read_tsv(target),
      error = function(e) NULL
    )
    if (!is.null(df) && nrow(df) > 0) {
      # 确保 Region 列存在 (从目录名提取作为备选)
      dir_region <- gsub("^Region_", "", basename(rd))
      if (!"Region" %in% names(df)) df$Region <- dir_region
      if (!"Genotype" %in% names(df)) df$Genotype <- genotype
      all_enrich[[basename(rd)]] <- df
      cat(sprintf("    %s: %d 行\n", basename(rd), nrow(df)))
    }
  }

  if (length(all_enrich) == 0) {
    cat("  ⚠ 未成功加载任何 Enrichment 数据\n")
    return(NULL)
  }

  combined <- bind_rows(all_enrich)
  cat(sprintf("  ✔ 6a 总计: %d 行 (覆盖 %d clusters × %d regions)\n",
              nrow(combined),
              length(unique(combined$Cluster)),
              length(unique(combined$Region))))
  return(combined)
}

################################################################################
# 模块 3: 合并 6a + 6b 数据
################################################################################

merge_6a_6b <- function(df_6b, df_6a) {
  cat("\n--- 合并 6a Enrichment → 6b Master List ---\n")

  if (is.null(df_6a) || nrow(df_6a) == 0) {
    cat("  6a 数据不可用, 6b 保持原样 (Enrichment 列设为 NA)\n")
    df_6b$Enrichment_P   <- NA_real_
    df_6b$Fold_Enrichment <- NA_real_
    return(df_6b)
  }

  # 6a 列名标准化
  enrich_slim <- df_6a %>%
    select(
      Cluster, kmer, Region, Genotype,
      any_of(c("p_value", "fold_enrichment", "fdr", "rank"))
    ) %>%
    rename(Enrichment_P = p_value) %>%
    {
      if ("fold_enrichment" %in% names(.)) rename(., Fold_Enrichment = fold_enrichment) else mutate(., Fold_Enrichment = NA_real_)
    }

  # 去重: 同一 (Cluster, Region, kmer, Genotype) 可能有多行, 取最显著
  enrich_slim <- enrich_slim %>%
    group_by(Cluster, Region, kmer, Genotype) %>%
    arrange(Enrichment_P) %>%
    slice_head(n = 1) %>%
    ungroup()

  n_before <- nrow(df_6b)
  merged <- df_6b %>%
    left_join(
      enrich_slim,
      by = c("Source_Cluster" = "Cluster", "kmer" = "kmer",
             "Region" = "Region", "Genotype" = "Genotype")
    )

  n_matched <- sum(!is.na(merged$Enrichment_P))
  cat(sprintf("  ✔ 合并完成: %d / %d 行有 Enrichment 匹配 (%.1f%%)\n",
              n_matched, n_before, 100 * n_matched / max(1, n_before)))
  cat(sprintf("  ℹ 保留全部 %d 行 (6b 为 Master)\n", nrow(merged)))

  # 填充缺失
  if (!"Enrichment_P"    %in% names(merged)) merged$Enrichment_P    <- NA_real_
  if (!"Fold_Enrichment" %in% names(merged)) merged$Fold_Enrichment <- NA_real_

  return(merged)
}

################################################################################
# 模块 4: Levenshtein 层次聚类 → 定义 Motif Family
################################################################################

assign_motif_families <- function(master_df, max_dist = 2) {
  cat("\n--- Levenshtein 聚类 (定义 Motif Family) ---\n")

  # 按 Source_Cluster + Region 分组
  master_df$group_key <- paste(master_df$Source_Cluster, master_df$Region, sep = "||")

  group_keys <- unique(master_df$group_key)
  cat(sprintf("  分组数: %d (Source_Cluster × Region)\n", length(group_keys)))

  # 全局 Motif Family 计数器
  global_fam_counter <- 0L
  result_list <- list()

  for (gk in group_keys) {
    sub <- master_df[master_df$group_key == gk, , drop = FALSE]
    kmers <- sub$kmer
    n <- length(kmers)

    if (n == 0) next

    if (n == 1) {
      local_fam <- 1L
    } else {
      dist_mat <- stringdistmatrix(kmers, kmers, method = "lv")
      hc <- hclust(as.dist(dist_mat), method = "complete")
      local_fam <- cutree(hc, h = max_dist)
    }

    sub$Motif_Family_Local <- local_fam
    sub$Motif_Family <- local_fam + global_fam_counter

    # 为 family 构建可读标签
    parts <- strsplit(gk, "\\|\\|")[[1]]
    src_cl <- parts[1]
    region <- parts[2]
    sub$Motif_Family_Label <- paste0(src_cl, "_", region, "_MF", local_fam)

    global_fam_counter <- global_fam_counter + max(local_fam)
    result_list[[gk]] <- sub
  }

  result <- bind_rows(result_list)
  n_families <- length(unique(result$Motif_Family))
  cat(sprintf("  ✔ 总 Motif Family 数: %d (cut h=%d)\n", n_families, max_dist))

  # 分布统计
  fam_sizes <- result %>% count(Motif_Family) %>% pull(n)
  cat(sprintf("    Family 大小: min=%d, median=%.0f, max=%d, mean=%.1f\n",
              min(fam_sizes), median(fam_sizes), max(fam_sizes), mean(fam_sizes)))
  cat(sprintf("    单例 Family: %d (%.1f%%)\n",
              sum(fam_sizes == 1), 100 * sum(fam_sizes == 1) / length(fam_sizes)))

  return(result)
}

################################################################################
# 模块 5: Shift-Aware Weighted Consensus PWM
################################################################################

#' 对单个 Motif Family 构建加权 consensus PWM
build_family_consensus <- function(kmers, importances, family_label,
                                   min_overlap = 4, trim_ic = 0.25) {
  n <- length(kmers)

  # 边界情况: 单条 k-mer
  if (n == 1) {
    motif <- kmer_to_motif(kmers[1], name = family_label,
                           nsites = max(1, round(importances[1] * 1000)))
    return(motif)
  }

  # 权重归一化: 将 importance 映射到 nsites (最小 10, 最大 10000)
  imp_norm <- importances / max(importances, na.rm = TRUE)
  nsites_vec <- pmax(10, round(imp_norm * 1000))

  # 创建 universalmotif 对象列表
  motif_list <- lapply(seq_len(n), function(i) {
    kmer_to_motif(kmers[i], name = paste0(family_label, "_m", i),
                  nsites = nsites_vec[i])
  })

  # 尝试 merge_motifs (shift-aware, RC-aware)
  consensus <- tryCatch({
    merged <- merge_motifs(motif_list,
                           method = "PCC",
                           min.overlap = min(min_overlap, min(nchar(kmers)) - 1),
                           tryRC = TRUE)
    # 设置名称
    merged@name <- family_label
    merged@nsites <- as.integer(sum(nsites_vec))
    merged
  }, error = function(e) {
    # Fallback: 用最高 importance 的 k-mer
    best_idx <- which.max(importances)
    fallback <- kmer_to_motif(kmers[best_idx], name = family_label,
                              nsites = sum(nsites_vec))
    fallback@extrainfo <- c(merge_error = e$message)
    fallback
  })

  # Trim low-IC edges
  consensus <- tryCatch({
    trim_motifs(consensus, min.ic = trim_ic)
  }, error = function(e) consensus)

  return(consensus)
}

#' 对所有 Motif Family 批量构建 consensus PWM
build_all_consensus <- function(master_df, min_overlap = 4, trim_ic = 0.25, n_cores = 1) {
  cat("\n--- 构建 Shift-Aware Weighted Consensus PWMs ---\n")

  families <- split(master_df, master_df$Motif_Family)
  n_fam <- length(families)
  cat(sprintf("  Motif Family 数: %d\n", n_fam))

  # 并行或串行处理
  if (n_cores > 1 && n_fam > 10) {
    cat(sprintf("  使用 %d 核并行...\n", min(n_cores, n_fam)))
    pwm_list <- mclapply(families, function(fam_df) {
      build_family_consensus(
        kmers = fam_df$kmer,
        importances = fam_df$mean_importance,
        family_label = fam_df$Motif_Family_Label[1],
        min_overlap = min_overlap,
        trim_ic = trim_ic
      )
    }, mc.cores = min(n_cores, n_fam))
  } else {
    pwm_list <- lapply(families, function(fam_df) {
      build_family_consensus(
        kmers = fam_df$kmer,
        importances = fam_df$mean_importance,
        family_label = fam_df$Motif_Family_Label[1],
        min_overlap = min_overlap,
        trim_ic = trim_ic
      )
    })
  }

  # 命名
  names(pwm_list) <- sapply(pwm_list, function(m) m@name)

  cat(sprintf("  ✔ 成功构建 %d 个 consensus PWM\n", length(pwm_list)))

  # PWM 长度统计
  widths <- sapply(pwm_list, function(m) ncol(m@motif))
  cat(sprintf("    PWM 宽度: min=%d, median=%.0f, max=%d\n",
              min(widths), median(widths), max(widths)))

  return(pwm_list)
}

################################################################################
# 模块 6: PCC 注释 — 高性能并行版 (Native C++ Multithreading)
# 策略: 分批构建矩阵 + OpenMP 并行 + 向量化提取 Top 2
################################################################################

#' 加载 MEME 数据库
load_meme_database <- function(meme_db_path) {
  cat(sprintf("\n--- 加载 MEME 数据库: %s ---\n", basename(meme_db_path)))
  db <- read_meme(meme_db_path)
  if (inherits(db, "universalmotif")) db <- list(db)
  cat(sprintf("  ✔ 加载 %d 个参考 motif\n", length(db)))
  return(db)
}

#' 强制将 motif 的 alphabet 设为 "DNA"
force_dna <- function(motif_obj) {
  if (motif_obj@alphabet != "DNA") {
    motif_obj@alphabet <- "DNA"
  }
  return(motif_obj)
}

#' 核心比对函数: 批量并行处理
#' @param queries    Query motifs 列表
#' @param db_motifs  Database motifs 列表
#' @param db_names   Database motifs 名称向量
#' @param n_cores    并行线程数
annotate_batch_parallel <- function(queries, db_motifs, db_names, 
                                    n_cores = 1, min_overlap = 4) {
  
  # 1. 预处理: 强制 DNA 模式
  queries <- lapply(queries, force_dna)
  # 注意: db_motifs 在外部已经处理过 force_dna，这里不再重复以节省时间
  
  # 2. 合并列表: c(Queries, DB)
  # compare_motifs 计算的是 All-vs-All 矩阵
  combined_motifs <- c(queries, db_motifs)
  
  n_q <- length(queries)
  n_db <- length(db_motifs)
  
  # 3. 高性能并行比对 (Native OpenMP)
  # 这里的 nthreads 是 C++ 线程，不会产生 R 的 fork 问题，非常稳定且快
  sim_matrix <- compare_motifs(
    combined_motifs,
    method = "PCC",
    min.overlap = min_overlap,
    tryRC = TRUE,
    nthreads = n_cores 
  )
  
  # 4. 提取子矩阵: Rows = Queries, Cols = DB
  # 矩阵结构: 前 n_q 行是 Queries，后 n_db 行是 DB
  # 我们只需要 [1:n_q, (n_q+1):end] 这部分
  sub_mat <- sim_matrix[1:n_q, (n_q + 1):(n_q + n_db), drop = FALSE]
  
  # 5. 向量化提取 Top 2
  # 使用 apply 遍历每一行(每个Query)
  results_list <- apply(sub_mat, 1, function(scores) {
    # 部分 NA 处理
    scores[is.na(scores)] <- -1 
    
    # 快速排序找最大索引 (partial sort 更快)
    # n-1, n 是最后两个（最大）
    n <- length(scores)
    if (n < 2) {
       return(c(db_names[1], round(scores[1], 4), NA, NA))
    }
    
    # order 默认升序，取最后两个
    top_idx <- order(scores, decreasing = TRUE)[1:2]
    
    val1 <- scores[top_idx[1]]
    val2 <- scores[top_idx[2]]
    
    # 如果 score 是 -1 (原 NA)，则还原为 NA
    c(
      db_names[top_idx[1]], 
      ifelse(val1 == -1, NA, sprintf("%.4f", val1)),
      db_names[top_idx[2]], 
      ifelse(val2 == -1, NA, sprintf("%.4f", val2))
    )
  })
  
  # 结果转置并整理为 DF
  res_df <- as.data.frame(t(results_list), stringsAsFactors = FALSE)
  colnames(res_df) <- c("DB_Match_1", "PCC_1", "DB_Match_2", "PCC_2")
  
  # 恢复数值类型
  res_df$PCC_1 <- as.numeric(res_df$PCC_1)
  res_df$PCC_2 <- as.numeric(res_df$PCC_2)
  
  return(res_df)
}

#' 主注释流程 (分批调度)
run_annotation <- function(master_df, pwm_list, db_motifs, 
                           n_cores = 1, min_overlap = 4) {
  
  cat(sprintf("\n--- PCC 注释 (Native Parallel, Threads=%d) ---\n", n_cores))
  
  # 预处理 DB (只做一次)
  db_motifs <- lapply(db_motifs, force_dna)
  db_names <- sapply(db_motifs, function(m) m@name)
  
  # 定义分批大小 (防止内存爆炸，每批 1000 个 query 比较合适)
  # 1000 queries x 1000 db = 1M elements matrix (8MB RAM), 非常安全
  BATCH_SIZE <- 1000
  
  # ========== Task A: Raw K-mers ==========
  cat("\n[Task A] 注释 Raw K-mers...\n")
  unique_kmers <- unique(master_df$kmer)
  n_kmers <- length(unique_kmers)
  cat(sprintf("  唯一 k-mer 数: %d\n", n_kmers))
  
  # 将 k-mer 转为对象列表
  kmer_motifs <- lapply(seq_along(unique_kmers), function(i) {
    kmer_to_motif(unique_kmers[i], name = unique_kmers[i])
  })
  
  # 分批循环
  n_batches <- ceiling(n_kmers / BATCH_SIZE)
  res_chunks <- list()
  
  for (b in 1:n_batches) {
    start_idx <- (b - 1) * BATCH_SIZE + 1
    end_idx <- min(b * BATCH_SIZE, n_kmers)
    
    cat(sprintf("\r  处理批次 %d/%d (k-mers %d-%d)...", b, n_batches, start_idx, end_idx))
    
    batch_queries <- kmer_motifs[start_idx:end_idx]
    
    res_chunks[[b]] <- annotate_batch_parallel(
      batch_queries, db_motifs, db_names, 
      n_cores = n_cores, min_overlap = min_overlap
    )
  }
  cat("\n")
  
  kmer_annot <- bind_rows(res_chunks)
  kmer_annot$kmer <- unique_kmers # 绑定 ID
  
  # 整理列顺序
  kmer_annot <- kmer_annot %>% select(kmer, everything())
  
  cat(sprintf("  ✔ Task A 完成: %d 个 k-mer 注释\n", nrow(kmer_annot)))
  
  # ========== Task B: Consensus PWMs ==========
  cat("\n[Task B] 注释 Consensus PWMs...\n")
  n_pwm <- length(pwm_list)
  cat(sprintf("  Consensus PWM 数: %d\n", n_pwm))
  
  # 分批循环 (PWM 数量也可能较多)
  n_batches_pwm <- ceiling(n_pwm / BATCH_SIZE)
  res_chunks_pwm <- list()
  
  for (b in 1:n_batches_pwm) {
    start_idx <- (b - 1) * BATCH_SIZE + 1
    end_idx <- min(b * BATCH_SIZE, n_pwm)
    
    cat(sprintf("\r  处理批次 %d/%d (PWMs %d-%d)...", b, n_batches_pwm, start_idx, end_idx))
    
    batch_queries <- pwm_list[start_idx:end_idx]
    
    res_chunks_pwm[[b]] <- annotate_batch_parallel(
      batch_queries, db_motifs, db_names, 
      n_cores = n_cores, min_overlap = min_overlap
    )
  }
  cat("\n")
  
  pwm_annot <- bind_rows(res_chunks_pwm)
  pwm_annot$Motif_Family_Label <- names(pwm_list)
  
  # 整理列顺序
  pwm_annot <- pwm_annot %>% select(Motif_Family_Label, everything())
  
  rename_map <- c(
    "DB_Match_1" = "Consensus_DB_Match_1", "PCC_1" = "Consensus_PCC_1",
    "DB_Match_2" = "Consensus_DB_Match_2", "PCC_2" = "Consensus_PCC_2"
  )
  colnames(pwm_annot)[match(names(rename_map), colnames(pwm_annot))] <- rename_map
  
  cat(sprintf("  ✔ Task B 完成: %d 个 consensus PWM 注释\n", nrow(pwm_annot)))
  
  return(list(
    kmer_annot = kmer_annot %>% 
      rename(Kmer_DB_Match_1 = DB_Match_1, Kmer_PCC_1 = PCC_1,
             Kmer_DB_Match_2 = DB_Match_2, Kmer_PCC_2 = PCC_2),
    pwm_annot = pwm_annot
  ))
}

################################################################################
# 模块 7: 生成输出表格
################################################################################

#' 生成 Output 1: Integrated_Features_Full.txt
build_integrated_table <- function(master_df, kmer_annot_df) {
  cat("\n--- 生成 Integrated_Features_Full.txt ---\n")

  # left_join k-mer 注释 (join key = kmer)
  result <- master_df %>%
    left_join(kmer_annot_df, by = "kmer")

  # 选择/重命名最终列
  final_cols <- c(
    "Genotype", "Source_Cluster", "Region", "Motif_Family", "Motif_Family_Label",
    "kmer", "feature", "mean_importance", "sd_importance", "selection_freq", "rank",
    "Enrichment_P", "Fold_Enrichment",
    "Kmer_DB_Match_1", "Kmer_PCC_1", "Kmer_DB_Match_2", "Kmer_PCC_2"
  )

  # 仅保留存在的列
  available <- intersect(final_cols, names(result))
  result <- result[, available, drop = FALSE]

  # 美化列名
  rename_map <- c(
    "kmer"             = "Kmer",
    "feature"          = "Feature_ID",
    "mean_importance"  = "Importance",
    "sd_importance"    = "Importance_SD",
    "selection_freq"   = "Selection_Freq",
    "rank"             = "RF_Rank"
  )
  for (old_name in names(rename_map)) {
    if (old_name %in% names(result)) {
      names(result)[names(result) == old_name] <- rename_map[[old_name]]
    }
  }

  result <- result %>% arrange(Genotype, Source_Cluster, Region, Motif_Family)

  cat(sprintf("  ✔ %d 行 × %d 列\n", nrow(result), ncol(result)))
  return(result)
}

#' 生成 Output 2: Consensus_Motif_Summary.txt (修复 List Column 报错版)
build_consensus_summary <- function(master_df, pwm_list, pwm_annot_df) {
  cat("\n--- 生成 Consensus_Motif_Summary.txt ---\n")

  # 1. 按 Motif_Family 聚合统计
  fam_stats <- master_df %>%
    group_by(Genotype, Source_Cluster, Region,
             Motif_Family, Motif_Family_Label) %>%
    summarise(
      Num_Members      = n(),
      Total_Importance = sum(mean_importance, na.rm = TRUE),
      Mean_Importance  = mean(mean_importance, na.rm = TRUE),
      Max_Importance   = max(mean_importance, na.rm = TRUE),
      Members          = paste(kmer, collapse = ";"),
      .groups = "drop"
    )

  # 2. 提取 consensus 序列 (使用 strict 模式防止生成 List)
  # 使用 vapply 强制要求返回字符型，如果出错会直接暴露，而不是生成 list
  fam_stats$Consensus_Sequence <- vapply(fam_stats$Motif_Family_Label, function(label) {
    if (label %in% names(pwm_list)) {
      val <- safe_consensus(pwm_list[[label]])
      # 确保只返回长度为1的字符串
      if (length(val) == 0) return(NA_character_)
      if (length(val) > 1) return(val[1])
      return(val)
    } else {
      return(NA_character_)
    }
  }, FUN.VALUE = character(1))

  # 3. PWM 宽度 (使用 strict 模式)
  fam_stats$PWM_Width <- vapply(fam_stats$Motif_Family_Label, function(label) {
    if (label %in% names(pwm_list)) {
      return(as.integer(ncol(pwm_list[[label]]@motif)))
    } else {
      return(NA_integer_)
    }
  }, FUN.VALUE = integer(1))

  # 4. 合并 consensus 注释
  fam_stats <- fam_stats %>%
    left_join(pwm_annot_df, by = "Motif_Family_Label")

  fam_stats <- fam_stats %>%
    arrange(desc(Total_Importance))

  # [关键修复] 5. 最终安全检查：检测并展平所有 List 列
  # 这是解决 'unimplemented type list' 错误的终极保险
  list_cols <- sapply(fam_stats, is.list)
  if (any(list_cols)) {
    bad_cols <- names(fam_stats)[list_cols]
    cat(sprintf("  ⚠ 检测到潜在的 List 列，正在强制修复: %s\n", paste(bad_cols, collapse=", ")))
    
    for (col in bad_cols) {
      # 将 list 强制转换为字符向量，NULL 转为 NA
      fam_stats[[col]] <- sapply(fam_stats[[col]], function(x) {
        if (is.null(x) || length(x) == 0) return(NA_character_)
        # 如果里面有多个元素，用逗号连接
        paste(as.character(x), collapse = ",") 
      })
    }
  }

  cat(sprintf("  ✔ %d 行 (Motif Family)\n", nrow(fam_stats)))
  return(fam_stats)
}

################################################################################
# 模块 8: 处理单个基因型
################################################################################

process_genotype <- function(genotype, opt) {
  cat("\n")
  cat("================================================================================\n")
  cat(sprintf("  处理基因型: %s\n", genotype))
  cat("================================================================================\n")

  # --- Step A: 加载 6b (Master) ---
  df_6b <- load_6b_data(opt$p6b_dir, genotype)
  if (is.null(df_6b) || nrow(df_6b) == 0) {
    cat("  ⚠ 6b 数据为空, 跳过\n")
    return(NULL)
  }

  # --- Step B: 加载 6a (Enrichment) ---
  df_6a <- load_6a_data(opt$p6a_dir, genotype)

  # --- Step C: 合并 ---
  master <- merge_6a_6b(df_6b, df_6a)

  # --- Step D: Levenshtein 聚类 ---
  master <- assign_motif_families(master, max_dist = opt$cluster_distance)

  # --- Step E: 构建 Consensus PWMs ---
  pwm_list <- build_all_consensus(master,
                                   min_overlap = opt$min_overlap,
                                   trim_ic = opt$trim_ic,
                                   n_cores = opt$n_cores)

  # --- Step F: PCC 注释 (修复: 正确调用 run_annotation) ---
  db_motifs <- load_meme_database(opt$meme_db)
  annot_results <- run_annotation(master, pwm_list, db_motifs,
                                  n_cores = opt$n_cores,
                                  min_overlap = opt$min_overlap)

  # --- Step G: 组装输出 ---
  geno_out <- file.path(opt$output_dir, genotype)

  # Output 1: Integrated_Features_Full.txt
  integrated_tbl <- build_integrated_table(master, annot_results$kmer_annot)
  out_file_1 <- file.path(geno_out, "Integrated_Features_Full.txt")
  write.table(integrated_tbl, out_file_1,
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  ✔ 输出: %s\n", basename(out_file_1)))

  # Output 2: Consensus_Motif_Summary.txt
  consensus_tbl <- build_consensus_summary(master, pwm_list, annot_results$pwm_annot)
  out_file_2 <- file.path(geno_out, "Consensus_Motif_Summary.txt")
  write.table(consensus_tbl, out_file_2,
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  ✔ 输出: %s\n", basename(out_file_2)))

  # Output 3: Motif_PWMs.rds
  out_file_3 <- file.path(geno_out, "Motif_PWMs.rds")
  saveRDS(pwm_list, out_file_3)
  cat(sprintf("  ✔ 输出: %s (%d PWMs)\n", basename(out_file_3), length(pwm_list)))

  return(list(
    genotype      = genotype,
    integrated    = integrated_tbl,
    consensus     = consensus_tbl,
    pwm_list      = pwm_list,
    n_features    = nrow(integrated_tbl),
    n_families    = nrow(consensus_tbl)
  ))
}

################################################################################
# 主程序
################################################################################

main <- function() {
  all_results <- list()

  for (genotype in genotypes_to_process) {
    result <- tryCatch(
      process_genotype(genotype, opt),
      error = function(e) {
        cat(sprintf("\n  ✖ 处理 %s 时出错: %s\n", genotype, e$message))
        NULL
      }
    )
    if (!is.null(result)) {
      all_results[[genotype]] <- result
    }
  }

  # ===== 汇总报告 =====
  cat("\n")
  cat("================================================================================\n")
  cat("  Pipeline 6c v3.1 完成\n")
  cat("================================================================================\n\n")
  cat(sprintf("结束时间: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

  cat("【处理汇总】\n\n")
  for (genotype in names(all_results)) {
    r <- all_results[[genotype]]
    cat(sprintf("  %s:\n", genotype))
    cat(sprintf("    总特征数 (Raw K-mers): %d\n", r$n_features))
    cat(sprintf("    Motif Family 数      : %d\n", r$n_families))
    cat(sprintf("    Consensus PWM 数     : %d\n", length(r$pwm_list)))

    # 注释命中率
    if ("Kmer_PCC_1" %in% names(r$integrated)) {
      n_annotated <- sum(!is.na(r$integrated$Kmer_PCC_1))
      cat(sprintf("    K-mer 注释数         : %d (%.1f%%)\n",
                  n_annotated, 100 * n_annotated / max(1, r$n_features)))
    }
    if ("Consensus_PCC_1" %in% names(r$consensus)) {
      n_con_ann <- sum(!is.na(r$consensus$Consensus_PCC_1))
      cat(sprintf("    Consensus 注释数     : %d (%.1f%%)\n",
                  n_con_ann, 100 * n_con_ann / max(1, r$n_families)))
    }
    cat("\n")
  }

  cat("【输出文件说明】\n\n")
  cat("每个基因型目录下包含:\n")
  cat("  ├── Integrated_Features_Full.txt    # 全表 (行=Raw K-mer, 含6a+6b+注释)\n")
  cat("  ├── Consensus_Motif_Summary.txt     # 汇总 (行=Motif Family, 含consensus+注释)\n")
  cat("  └── Motif_PWMs.rds                  # universalmotif PWM 列表 (供下游绘图)\n")
  cat("\n")

  cat(sprintf("输出目录: %s\n", opt$output_dir))
  cat(sprintf("日志文件: %s\n", log_file))
  cat("\n分析完成。\n")

  invisible(all_results)
}

# 执行
results <- main()

# 关闭日志
sink(type = "message")
sink()
close(log_con)