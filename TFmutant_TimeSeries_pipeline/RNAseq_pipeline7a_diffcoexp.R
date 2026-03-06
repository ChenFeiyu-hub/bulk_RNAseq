#!/usr/bin/env Rscript
################################################################################
# RNAseq Pipeline 7a: 时滞差异共表达分析 (Lag-k Differential Co-expression)
#
# 版本: v1.3  (解析条件内检验 + Mfuzz过滤 + 检查点缓存)
#
# ============================================================================
# v1.3 更新日志 (基于 v1.2)
# ============================================================================
# 1. [CRITICAL FIX] 步骤 6: 置换检验 → 解析 t-test p 值
#    - 问题: v1.2 参数化 null 拟合仍受 p_floor=1e-04 限制,
#            且 null 分布 SD≈0.18 导致 |r|=0.95 → z≈3.84 → p≈1.2e-04,
#            经 BH 校正 2.28M 对后, 最好情况 FDR≈280, 永远无法 <0.05
#    - 方案: 使用 Spearman rho 的 t-approximation:
#            t = r × sqrt((n-2)/(1-r²)), df = n-2
#            产生连续 p 值, 无下界限制, 计算成本 O(1)
#    - 效果: 不再需要条件内置换循环, 运行时间从 ~25 min → <1 sec
#
# 2. [NEW] 步骤 4b: Mfuzz 共同基因过滤
#    - 将目标基因集从全部 DEGs (11660) 缩减为 Mfuzz 共同基因 (3268)
#    - 检验对数从 ~2.28M → ~640K, BH 校正负担降低 3.6 倍
#    - BH rank 1 门槛: 2.19e-08 → 7.81e-08
#    - 需要的 |r| 从 ~0.94 降至 ~0.93 (解析 t-test)
#
# 3. [NEW] 步骤 6b: 可选 |r| 预过滤
#    - 仅对 max(|r_WT|, |r_MUT|) >= 0.3 的对做 BH 校正
#    - 进一步降低多重检验负担 (通常 ~200K 对)
#    - 通过 --prefilter_r 参数控制 (默认 0, 即不启用)
#
# 4. [IMPROVED] 步骤 7: 保留参数化置换 (v1.2 已验证有效)
#    - v1.2 步骤 7 已检测到 56,694 条显著差异边
#    - 移除 p_floor 限制, 允许更精确的参数化 p 值
#    - 检查点机制保持不变
#
# 5. [IMPROVED] 步骤 8: 边分类逻辑优化
#    - 解析 p 值使 sig_WT/sig_MUT 能正常产生 TRUE 值
#    - 分类逻辑与 v1.2 完全一致, 但现在能实际触发所有分支
#
# ============================================================================
# 方法学基础
# ============================================================================
# 条件内显著性 (步骤 6):
#   - Spearman rho 的 t-approximation (解析, 无置换)
#   - t = r × sqrt(df / (1 - r²)), df = n_eff - 2
#   - 双侧 p 值: 2 × P(T > |t|), T ~ t(df)
#   - BH FDR 校正
#
# 差异显著性 (步骤 7):
#   - 跨条件标签置换 + 参数化 null 拟合 (v1.2 方法, 已验证有效)
#   - Fisher Z-test: z_diff = (atanh(r_WT) - atanh(r_MUT)) / SE
#   - B 次置换估计 null 分布 → 正态拟合 → 连续 p 值
#   - BH FDR 校正
#
# ============================================================================
# 分析流程
# ============================================================================
# 1. 加载 Pipeline 4 VST 数据 + Pipeline 5 Mfuzz 结果
# 2. 加载 TF 列表, 过滤为数据集中实际表达的 TF
# 3. 检测 metadata 列名
# 4. 构建 Replicate-Aware Lag-k Design-Response 矩阵
# 4b. [NEW] Mfuzz 共同基因过滤 (缩减目标基因集)
# 5. 计算观测 Spearman lag-k 相关 (WT, MUT) + Fisher Z 差异
# 6. [NEW] 条件内解析 t-test 显著性 (替代置换检验)
# 6b. [NEW] 可选 |r| 预过滤 + 重新 BH 校正
# 7. 跨条件差异置换检验 (参数化 null, v1.2 方法)  ← 检查点
# 8. 边分类 (Lost / Gained / Rewired / Conserved / NS)
# 9. 输出边表 + TF 摘要 + RData
#
# ============================================================================
# 使用方式
# ============================================================================
#   Rscript RNAseq_pipeline7a_diffcoexp.R \
#     --rdata_p4 /path/to/Pipeline4_Logic1_TimeCourse_for_Mfuzz.RData \
#     --rdata_p5 /path/to/Mfuzz_v11_Results.RData \
#     --tf_list /path/to/TF_list.txt \
#     --out_dir /path/to/07_GRN \
#     --condition_wt WT --condition_mut dof \
#     --n_cores 40
#
#   新增参数:
#     --prefilter_r 0.3    # 启用 |r| 预过滤 (默认 0 = 不启用)
#     --target_filter mfuzz # 目标基因过滤策略: all/mfuzz/non_conserved
#
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(parallel)
})

################################################################################
# 参数解析
################################################################################

parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  params <- list(
    # ========== 数据输入 ==========
    rdata_p4       = NULL,     # [必需] Pipeline4 RData (含 VST + metadata + 时间响应基因)
    rdata_p5       = NULL,     # [可选] Pipeline5 Mfuzz RData (用于注释聚类归属 + 基因过滤)
    tf_list        = NULL,     # [必需] TF 列表文件 (TSV: TF_ID, Gene_ID, Family)
    highlight_file = NULL,     # [可选] 重点基因标注文件 (TSV: Symbol, Gene)
    samples_file   = NULL,     # [可选] 覆盖 metadata
    
    # ========== 输出控制 ==========
    out_dir      = NULL,       # [可选] 输出目录 (自动推断)
    project_name = NULL,       # [可选] 项目名称
    
    # ========== 条件标签 ==========
    condition_wt  = "WT",
    condition_mut = "mut",
    
    # ========== 算法参数 ==========
    lag            = 1,        # 时滞步数 (默认 1)
    n_permutations = 1000,     # 置换次数 (仅用于步骤 7 差异检验)
    fdr_cutoff     = 0.05,     # FDR 阈值 (WT/MUT 显著性)
    fdr_cutoff_diff = 0.05,    # FDR 阈值 (差异显著性)
    fdr_mut_loose  = 0.10,     # MUT 宽松阈值 (Lost Edge 定义)
    cor_method     = "spearman", # 相关方法
    
    # ========== v1.3 新增参数 ==========
    target_filter  = "mfuzz",  # 目标基因过滤: "all" / "mfuzz" / "non_conserved"
    prefilter_r    = 0,        # |r| 预过滤阈值 (0 = 不启用, 推荐 0.3)
    
    # ========== 输出模式 ==========
    output_mode = "significant",
    
    # ========== 计算资源 ==========
    n_cores   = 1,             # 并行核数 (仅步骤 7 使用)
    seed      = 42,            # 随机种子
    
    # ========== 检查点控制 ==========
    force_rerun = FALSE,       # 强制重算 (忽略检查点)
    
    # ========== 调试 ==========
    verbose = TRUE
  )
  
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (grepl("^--", arg)) {
      key <- gsub("-", "_", sub("^--", "", arg))
      if (grepl("=", key)) {
        parts <- strsplit(key, "=")[[1]]
        key <- parts[1]; value <- parts[2]
      } else if (i + 1 <= length(args) && !grepl("^--", args[i + 1])) {
        i <- i + 1; value <- args[i]
      } else {
        value <- "TRUE"
      }
      if (key %in% names(params)) {
        if (is.logical(params[[key]])) {
          params[[key]] <- as.logical(value)
        } else if (is.numeric(params[[key]]) && !is.null(params[[key]])) {
          params[[key]] <- as.numeric(value)
        } else if (is.null(params[[key]]) && !is.na(suppressWarnings(as.numeric(value)))) {
          params[[key]] <- as.numeric(value)
        } else {
          params[[key]] <- value
        }
      }
    }
    i <- i + 1
  }
  
  return(params)
}

PARAMS <- parse_arguments()

################################################################################
# 辅助函数: 时间点解析 (复用 Pipeline 5 逻辑)
################################################################################

extract_time_numeric <- function(time_label) {
  if (is.na(time_label) || time_label == "") return(NA)
  cleaned <- gsub("[^0-9.]", "", as.character(time_label))
  if (cleaned == "" || cleaned == ".") return(NA)
  return(suppressWarnings(as.numeric(cleaned)))
}

################################################################################
# 核心函数 1: 构建 Replicate-Aware Lag-k 矩阵
################################################################################

build_lag_matrices <- function(vst_matrix, metadata, condition, gene_set,
                                time_col, cond_col, lag = 1) {
  
  cond_mask <- as.character(metadata[[cond_col]]) == condition
  cond_samples <- rownames(metadata)[cond_mask]
  
  if (length(cond_samples) == 0) {
    stop(sprintf("No samples found for condition '%s'", condition))
  }
  
  sample_times <- data.frame(
    sample = cond_samples,
    time_label = as.character(metadata[cond_samples, time_col]),
    stringsAsFactors = FALSE
  )
  sample_times$time_numeric <- sapply(sample_times$time_label, extract_time_numeric)
  sample_times <- sample_times[!is.na(sample_times$time_numeric), ]
  
  sorted_times <- sort(unique(sample_times$time_numeric))
  n_times <- length(sorted_times)
  
  if (n_times <= lag) {
    stop(sprintf("Not enough timepoints (%d) for lag=%d", n_times, lag))
  }
  
  sample_times$rep_idx <- NA
  for (t_val in sorted_times) {
    t_mask <- abs(sample_times$time_numeric - t_val) < 0.1
    t_samples <- sample_times$sample[t_mask]
    sample_times$rep_idx[t_mask] <- seq_along(t_samples)
  }
  
  n_reps <- max(sample_times$rep_idx)
  n_lag_pairs <- n_times - lag
  n_eff <- n_reps * n_lag_pairs
  
  available_genes <- intersect(gene_set, rownames(vst_matrix))
  if (length(available_genes) == 0) {
    stop("No genes from gene_set found in VST matrix")
  }
  
  design_mat  <- matrix(NA, nrow = n_eff, ncol = length(available_genes))
  response_mat <- matrix(NA, nrow = n_eff, ncol = length(available_genes))
  colnames(design_mat) <- available_genes
  colnames(response_mat) <- available_genes
  
  rep_blocks <- integer(n_eff)
  row_info <- data.frame(
    row = integer(n_eff), rep = integer(n_eff),
    t_design = numeric(n_eff), t_response = numeric(n_eff),
    stringsAsFactors = FALSE
  )
  
  row_idx <- 0
  for (r in 1:n_reps) {
    for (p in 1:n_lag_pairs) {
      row_idx <- row_idx + 1
      
      t_design  <- sorted_times[p]
      t_response <- sorted_times[p + lag]
      
      design_sample <- sample_times$sample[
        abs(sample_times$time_numeric - t_design) < 0.1 & sample_times$rep_idx == r
      ]
      response_sample <- sample_times$sample[
        abs(sample_times$time_numeric - t_response) < 0.1 & sample_times$rep_idx == r
      ]
      
      if (length(design_sample) == 1 && length(response_sample) == 1) {
        design_mat[row_idx, ]  <- vst_matrix[available_genes, design_sample]
        response_mat[row_idx, ] <- vst_matrix[available_genes, response_sample]
      }
      
      rep_blocks[row_idx] <- r
      row_info$row[row_idx] <- row_idx
      row_info$rep[row_idx] <- r
      row_info$t_design[row_idx] <- t_design
      row_info$t_response[row_idx] <- t_response
    }
  }
  
  complete_rows <- complete.cases(design_mat) & complete.cases(response_mat)
  if (sum(complete_rows) < n_eff) {
    cat(sprintf("  Warning: %d/%d rows incomplete, removing\n", 
                n_eff - sum(complete_rows), n_eff))
  }
  
  design_mat  <- design_mat[complete_rows, , drop = FALSE]
  response_mat <- response_mat[complete_rows, , drop = FALSE]
  rep_blocks   <- rep_blocks[complete_rows]
  row_info     <- row_info[complete_rows, ]
  n_eff        <- sum(complete_rows)
  
  return(list(
    design   = design_mat,
    response = response_mat,
    n_eff    = n_eff,
    n_reps   = n_reps,
    n_times  = n_times,
    n_lag_pairs = n_lag_pairs,
    rep_blocks  = rep_blocks,
    row_info    = row_info,
    sorted_times = sorted_times,
    genes    = available_genes
  ))
}

################################################################################
# 核心函数 2: 矩阵级相关计算
################################################################################

calc_cor_matrix <- function(tf_mat, target_mat, method = "spearman") {
  r_mat <- cor(tf_mat, target_mat, method = method, use = "pairwise.complete.obs")
  return(r_mat)
}

################################################################################
# 核心函数 3: 解析条件内 p 值 (t-approximation)  ★ v1.3 NEW ★
################################################################################

#' 条件内解析 p 值 (Spearman rho t-approximation)
#'
#' v1.3 替代方案 (替换 v1.2 置换检验):
#'   - 对 Spearman rho 使用 t-approximation:
#'     t = r × sqrt(df / (1 - r²)), df = n_eff - 2
#'   - 双侧 p 值: 2 × pt(-|t|, df)
#'   - 优势: 连续 p 值 (无 1/B 下限), O(1) 计算成本
#'   - 当 n_eff >= 10 时, t-approximation 精度极佳
#'
#' @param r_mat 相关系数矩阵 (TF × Target)
#' @param n_eff 有效样本量
#' @return list(p_values, fdr_values, method)
calc_analytical_pvalues <- function(r_mat, n_eff, fdr_method = "BH") {
  
  df <- n_eff - 2
  
  if (df < 1) {
    warning("Degrees of freedom < 1 (n_eff=", n_eff, "), all p-values set to 1")
    p_mat <- matrix(1, nrow = nrow(r_mat), ncol = ncol(r_mat),
                    dimnames = dimnames(r_mat))
    fdr_mat <- p_mat
    return(list(p_values = p_mat, fdr_values = fdr_mat,
                method = "analytical_t_test", df = df, n_eff = n_eff))
  }
  
  # Clip to avoid Inf at |r| = 1
  r_clipped <- pmin(pmax(r_mat, -0.9999), 0.9999)
  
  # t-statistic
  t_stat <- r_clipped * sqrt(df / (1 - r_clipped^2))
  
  # Two-sided p-value from t-distribution
  p_mat <- 2 * pt(-abs(t_stat), df = df)
  dimnames(p_mat) <- dimnames(r_mat)
  
  # Handle NA (from NA correlations)
  p_mat[is.na(p_mat)] <- 1
  
  # BH FDR correction
  fdr_mat <- matrix(p.adjust(as.vector(p_mat), method = fdr_method),
                    nrow = nrow(p_mat), ncol = ncol(p_mat),
                    dimnames = dimnames(p_mat))
  
  return(list(
    p_values   = p_mat,
    fdr_values = fdr_mat,
    t_stats    = t_stat,
    method     = "analytical_t_test",
    df         = df,
    n_eff      = n_eff
  ))
}

################################################################################
# 核心函数 4: 跨条件差异置换检验  (v1.2 方法, 已验证有效)
################################################################################

#' v1.3 变更: 移除 p_floor 限制, 允许更精确的参数化 p 值
run_differential_permutation <- function(vst_matrix, metadata,
                                          tf_genes, target_genes,
                                          n_eff_wt, n_eff_mut,
                                          z_diff_obs,
                                          condition_wt, condition_mut,
                                          time_col, cond_col,
                                          lag = 1, B = 1000,
                                          method = "spearman",
                                          seed = 42, n_cores = 1,
                                          verbose = TRUE) {
  
  n_tf <- length(tf_genes)
  n_target <- length(target_genes)
  
  abs_z_diff_obs <- abs(z_diff_obs)
  
  # ---- 预计算时间点分组信息 ----
  cond_mask_wt <- as.character(metadata[[cond_col]]) == condition_wt
  cond_mask_mut <- as.character(metadata[[cond_col]]) == condition_mut
  
  samples_wt <- rownames(metadata)[cond_mask_wt]
  samples_mut <- rownames(metadata)[cond_mask_mut]
  
  all_samples <- c(samples_wt, samples_mut)
  sample_info <- data.frame(
    sample = all_samples,
    condition = c(rep(condition_wt, length(samples_wt)),
                  rep(condition_mut, length(samples_mut))),
    time_numeric = sapply(metadata[all_samples, time_col], extract_time_numeric),
    stringsAsFactors = FALSE
  )
  sample_info <- sample_info[!is.na(sample_info$time_numeric), ]
  
  sorted_times <- sort(unique(sample_info$time_numeric))
  
  samples_by_time <- list()
  for (t_val in sorted_times) {
    t_mask <- abs(sample_info$time_numeric - t_val) < 0.1
    samples_by_time[[as.character(t_val)]] <- sample_info[t_mask, ]
  }
  
  n_wt_per_time <- sum(abs(sample_info$time_numeric[sample_info$condition == condition_wt] -
                           sorted_times[1]) < 0.1)
  n_mut_per_time <- sum(abs(sample_info$time_numeric[sample_info$condition == condition_mut] -
                            sorted_times[1]) < 0.1)
  
  all_genes <- union(tf_genes, target_genes)
  
  available_cores <- detectCores(logical = FALSE)
  if (is.na(available_cores)) available_cores <- 1
  n_cores <- min(n_cores, available_cores, B)
  n_cores <- max(n_cores, 1)
  
  if (verbose) {
    cat(sprintf("    Running %d differential permutations on %d cores...\n", B, n_cores))
  }
  
  chunk_sizes <- rep(B %/% n_cores, n_cores)
  remainder <- B %% n_cores
  if (remainder > 0) chunk_sizes[1:remainder] <- chunk_sizes[1:remainder] + 1
  
  set.seed(seed + 9999)
  chunk_seeds <- sample.int(.Machine$integer.max, n_cores)
  
  # ---- Worker: 返回 count, sum, sum_sq ----
  worker_fn <- function(chunk_idx) {
    local_B <- chunk_sizes[chunk_idx]
    local_count  <- matrix(0L,  nrow = n_tf, ncol = n_target,
                           dimnames = list(tf_genes, target_genes))
    local_sum    <- matrix(0.0, nrow = n_tf, ncol = n_target,
                           dimnames = list(tf_genes, target_genes))
    local_sum_sq <- matrix(0.0, nrow = n_tf, ncol = n_target,
                           dimnames = list(tf_genes, target_genes))
    
    set.seed(chunk_seeds[chunk_idx])
    
    for (b in 1:local_B) {
      
      pseudo_wt_samples <- character(0)
      pseudo_mut_samples <- character(0)
      
      for (t_key in names(samples_by_time)) {
        t_info <- samples_by_time[[t_key]]
        t_samples <- t_info$sample
        shuffled <- sample(t_samples)
        pseudo_wt_samples <- c(pseudo_wt_samples, shuffled[1:n_wt_per_time])
        pseudo_mut_samples <- c(pseudo_mut_samples, shuffled[(n_wt_per_time + 1):length(shuffled)])
      }
      
      pseudo_meta <- metadata[c(pseudo_wt_samples, pseudo_mut_samples), ]
      pseudo_meta$pseudo_cond <- c(rep("pseudo_WT", length(pseudo_wt_samples)),
                                    rep("pseudo_MUT", length(pseudo_mut_samples)))
      rownames(pseudo_meta) <- c(pseudo_wt_samples, pseudo_mut_samples)
      
      lag_pseudo_wt <- tryCatch(
        build_lag_matrices(vst_matrix, pseudo_meta, "pseudo_WT", all_genes,
                           time_col, "pseudo_cond", lag),
        error = function(e) NULL
      )
      lag_pseudo_mut <- tryCatch(
        build_lag_matrices(vst_matrix, pseudo_meta, "pseudo_MUT", all_genes,
                           time_col, "pseudo_cond", lag),
        error = function(e) NULL
      )
      
      if (is.null(lag_pseudo_wt) || is.null(lag_pseudo_mut)) next
      
      tf_in_pwt <- intersect(tf_genes, colnames(lag_pseudo_wt$design))
      tg_in_pwt <- intersect(target_genes, colnames(lag_pseudo_wt$response))
      tf_in_pmut <- intersect(tf_genes, colnames(lag_pseudo_mut$design))
      tg_in_pmut <- intersect(target_genes, colnames(lag_pseudo_mut$response))
      
      common_tf <- intersect(tf_in_pwt, tf_in_pmut)
      common_tg <- intersect(tg_in_pwt, tg_in_pmut)
      
      if (length(common_tf) == 0 || length(common_tg) == 0) next
      
      r_pwt <- calc_cor_matrix(lag_pseudo_wt$design[, common_tf, drop = FALSE],
                                lag_pseudo_wt$response[, common_tg, drop = FALSE],
                                method)
      r_pmut <- calc_cor_matrix(lag_pseudo_mut$design[, common_tf, drop = FALSE],
                                 lag_pseudo_mut$response[, common_tg, drop = FALSE],
                                 method)
      
      z_pwt <- atanh(pmin(pmax(r_pwt, -0.9999), 0.9999))
      z_pmut <- atanh(pmin(pmax(r_pmut, -0.9999), 0.9999))
      se <- sqrt(1 / max(1, lag_pseudo_wt$n_eff - 3) +
                 1 / max(1, lag_pseudo_mut$n_eff - 3))
      z_diff_perm <- (z_pwt - z_pmut) / se
      
      shared_tf <- intersect(common_tf, tf_genes)
      shared_tg <- intersect(common_tg, target_genes)
      if (length(shared_tf) > 0 && length(shared_tg) > 0) {
        z_perm_sub <- z_diff_perm[shared_tf, shared_tg, drop = FALSE]
        z_obs_sub  <- abs_z_diff_obs[shared_tf, shared_tg, drop = FALSE]
        
        abs_z_perm <- abs(z_perm_sub)
        abs_z_perm[is.na(abs_z_perm)] <- 0
        
        # 经验计数
        exceed_mask <- abs_z_perm >= z_obs_sub
        exceed_mask[is.na(exceed_mask)] <- FALSE
        local_count[shared_tf, shared_tg] <- local_count[shared_tf, shared_tg] + exceed_mask
        
        # 累计矩 (用于参数化 null)
        local_sum[shared_tf, shared_tg]    <- local_sum[shared_tf, shared_tg]    + abs_z_perm
        local_sum_sq[shared_tf, shared_tg] <- local_sum_sq[shared_tf, shared_tg] + abs_z_perm^2
      }
    }
    
    return(list(count = local_count, sum = local_sum, sum_sq = local_sum_sq,
                n_perm = local_B))
  }
  
  t_start <- Sys.time()
  
  if (n_cores == 1) {
    count_exceed <- matrix(0L,  nrow = n_tf, ncol = n_target,
                           dimnames = list(tf_genes, target_genes))
    null_sum     <- matrix(0.0, nrow = n_tf, ncol = n_target,
                           dimnames = list(tf_genes, target_genes))
    null_sum_sq  <- matrix(0.0, nrow = n_tf, ncol = n_target,
                           dimnames = list(tf_genes, target_genes))
    set.seed(seed + 9999)
    report_interval <- max(1, B %/% 10)
    
    for (b in 1:B) {
      pseudo_wt_samples <- character(0)
      pseudo_mut_samples <- character(0)
      
      for (t_key in names(samples_by_time)) {
        t_info <- samples_by_time[[t_key]]
        t_samples <- t_info$sample
        shuffled <- sample(t_samples)
        pseudo_wt_samples <- c(pseudo_wt_samples, shuffled[1:n_wt_per_time])
        pseudo_mut_samples <- c(pseudo_mut_samples, shuffled[(n_wt_per_time + 1):length(shuffled)])
      }
      
      pseudo_meta <- metadata[c(pseudo_wt_samples, pseudo_mut_samples), ]
      pseudo_meta$pseudo_cond <- c(rep("pseudo_WT", length(pseudo_wt_samples)),
                                    rep("pseudo_MUT", length(pseudo_mut_samples)))
      rownames(pseudo_meta) <- c(pseudo_wt_samples, pseudo_mut_samples)
      
      lag_pseudo_wt <- tryCatch(
        build_lag_matrices(vst_matrix, pseudo_meta, "pseudo_WT", all_genes,
                           time_col, "pseudo_cond", lag),
        error = function(e) NULL
      )
      lag_pseudo_mut <- tryCatch(
        build_lag_matrices(vst_matrix, pseudo_meta, "pseudo_MUT", all_genes,
                           time_col, "pseudo_cond", lag),
        error = function(e) NULL
      )
      
      if (is.null(lag_pseudo_wt) || is.null(lag_pseudo_mut)) next
      
      tf_in_pwt <- intersect(tf_genes, colnames(lag_pseudo_wt$design))
      tg_in_pwt <- intersect(target_genes, colnames(lag_pseudo_wt$response))
      tf_in_pmut <- intersect(tf_genes, colnames(lag_pseudo_mut$design))
      tg_in_pmut <- intersect(target_genes, colnames(lag_pseudo_mut$response))
      
      common_tf <- intersect(tf_in_pwt, tf_in_pmut)
      common_tg <- intersect(tg_in_pwt, tg_in_pmut)
      
      if (length(common_tf) == 0 || length(common_tg) == 0) next
      
      r_pwt <- calc_cor_matrix(lag_pseudo_wt$design[, common_tf, drop = FALSE],
                                lag_pseudo_wt$response[, common_tg, drop = FALSE],
                                method)
      r_pmut <- calc_cor_matrix(lag_pseudo_mut$design[, common_tf, drop = FALSE],
                                 lag_pseudo_mut$response[, common_tg, drop = FALSE],
                                 method)
      
      z_pwt <- atanh(pmin(pmax(r_pwt, -0.9999), 0.9999))
      z_pmut <- atanh(pmin(pmax(r_pmut, -0.9999), 0.9999))
      se <- sqrt(1 / max(1, lag_pseudo_wt$n_eff - 3) +
                 1 / max(1, lag_pseudo_mut$n_eff - 3))
      z_diff_perm <- (z_pwt - z_pmut) / se
      
      shared_tf <- intersect(common_tf, tf_genes)
      shared_tg <- intersect(common_tg, target_genes)
      if (length(shared_tf) > 0 && length(shared_tg) > 0) {
        z_perm_sub <- z_diff_perm[shared_tf, shared_tg, drop = FALSE]
        z_obs_sub  <- abs_z_diff_obs[shared_tf, shared_tg, drop = FALSE]
        
        abs_z_perm <- abs(z_perm_sub)
        abs_z_perm[is.na(abs_z_perm)] <- 0
        
        exceed_mask <- abs_z_perm >= z_obs_sub
        exceed_mask[is.na(exceed_mask)] <- FALSE
        count_exceed[shared_tf, shared_tg] <- count_exceed[shared_tf, shared_tg] + exceed_mask
        null_sum[shared_tf, shared_tg]     <- null_sum[shared_tf, shared_tg]     + abs_z_perm
        null_sum_sq[shared_tf, shared_tg]  <- null_sum_sq[shared_tf, shared_tg]  + abs_z_perm^2
      }
      
      if (verbose && b %% report_interval == 0) {
        cat(sprintf("      Differential permutation %d/%d (%.0f%%)\n", b, B, 100 * b / B))
      }
    }
  } else {
    if (verbose) cat(sprintf("    Dispatching to %d workers...\n", n_cores))
    
    results <- mclapply(1:n_cores, worker_fn, mc.cores = n_cores,
                        mc.preschedule = TRUE)
    
    errors <- sapply(results, function(x) inherits(x, "try-error"))
    if (any(errors)) {
      warning(sprintf("  %d/%d workers failed", sum(errors), n_cores))
      results <- results[!errors]
    }
    
    count_exceed <- Reduce(`+`, lapply(results, `[[`, "count"))
    null_sum     <- Reduce(`+`, lapply(results, `[[`, "sum"))
    null_sum_sq  <- Reduce(`+`, lapply(results, `[[`, "sum_sq"))
  }
  
  t_elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  if (verbose) cat(sprintf("    Completed in %.1f seconds\n", t_elapsed))
  
  # ---- 参数化 null 拟合 ----
  null_mean <- null_sum / B
  null_var  <- (null_sum_sq / B) - null_mean^2
  null_var[null_var < 0] <- 0
  null_sd   <- sqrt(null_var)
  
  EPSILON <- 1e-10
  sd_valid <- null_sd > EPSILON
  
  z_scores <- matrix(NA_real_, nrow = n_tf, ncol = n_target,
                     dimnames = list(tf_genes, target_genes))
  z_scores[sd_valid] <- (abs_z_diff_obs[sd_valid] - null_mean[sd_valid]) / null_sd[sd_valid]
  
  p_parametric <- matrix(1.0, nrow = n_tf, ncol = n_target,
                         dimnames = list(tf_genes, target_genes))
  p_parametric[sd_valid] <- pnorm(z_scores[sd_valid], lower.tail = FALSE)
  
  p_empirical <- (count_exceed + 1L) / (B + 1L)
  p_parametric[!sd_valid] <- p_empirical[!sd_valid]
  
  # v1.3: 移除 p_floor 限制 (v1.2 中 p_floor = 1/(B*10) 会截断有效信号)
  # 参数化 p 值来自正态拟合, 其精度不受 B 限制, 无需人为设下限
  
  dimnames(p_parametric) <- dimnames(z_diff_obs)
  dimnames(p_empirical)  <- dimnames(z_diff_obs)
  
  fdr_parametric <- matrix(p.adjust(as.vector(p_parametric), method = "BH"),
                           nrow = n_tf, ncol = n_target)
  fdr_empirical  <- matrix(p.adjust(as.vector(p_empirical), method = "BH"),
                           nrow = n_tf, ncol = n_target)
  dimnames(fdr_parametric) <- dimnames(z_diff_obs)
  dimnames(fdr_empirical)  <- dimnames(z_diff_obs)
  
  if (verbose) {
    cat(sprintf("    Null distribution |z_diff|: mean=%.4f, sd=%.4f (median across pairs)\n",
                median(null_mean, na.rm = TRUE), median(null_sd, na.rm = TRUE)))
    cat(sprintf("    Parametric p-value range: [%.2e, %.4f]\n",
                min(p_parametric, na.rm = TRUE), max(p_parametric, na.rm = TRUE)))
    cat(sprintf("    Empirical p-value range:  [%.2e, %.4f]\n",
                min(p_empirical, na.rm = TRUE), max(p_empirical, na.rm = TRUE)))
    n_sig_param <- sum(fdr_parametric < 0.05, na.rm = TRUE)
    n_sig_emp   <- sum(fdr_empirical < 0.05, na.rm = TRUE)
    cat(sprintf("    Significant edges (FDR<0.05): parametric=%d, empirical=%d\n",
                n_sig_param, n_sig_emp))
  }
  
  return(list(
    p_values      = p_parametric,
    fdr_values    = fdr_parametric,
    p_empirical   = p_empirical,
    fdr_empirical = fdr_empirical,
    null_mean     = null_mean,
    null_sd       = null_sd,
    count_exceed  = count_exceed
  ))
}

################################################################################
# 核心函数 5: 边分类
################################################################################

classify_edges <- function(sig_wt, sig_mut, sig_diff, sig_mut_loose) {
  
  edge_class <- rep("NS", length(sig_wt))
  
  edge_class[sig_wt & !sig_mut & !sig_diff]  <- "WT_only"
  edge_class[!sig_wt & sig_mut & !sig_diff]  <- "MUT_only"
  edge_class[sig_wt & sig_mut & !sig_diff] <- "Conserved"
  edge_class[sig_wt & sig_diff & !sig_mut_loose] <- "Lost"
  edge_class[!sig_wt & sig_diff & sig_mut] <- "Gained"
  edge_class[sig_wt & sig_diff & sig_mut_loose] <- "Rewired"
  
  return(edge_class)
}

################################################################################
# 验证参数
################################################################################

if (is.null(PARAMS$rdata_p4) || is.null(PARAMS$tf_list)) {
  cat("================================================================================\n")
  cat("  Pipeline 7a: 时滞差异共表达分析 (v1.3)\n")
  cat("================================================================================\n\n")
  cat("错误: 必须指定 --rdata_p4 和 --tf_list\n\n")
  cat("使用方法:\n")
  cat("  Rscript RNAseq_pipeline7a_diffcoexp.R \\\n")
  cat("    --rdata_p4 /path/to/Pipeline4_Logic1_TimeCourse_for_Mfuzz.RData \\\n")
  cat("    --rdata_p5 /path/to/Mfuzz_v11_Results.RData \\\n")
  cat("    --tf_list /path/to/TF_list.txt \\\n")
  cat("    --condition_wt WT --condition_mut dof \\\n")
  cat("    --n_cores 40\n\n")
  cat("核心参数:\n")
  cat("  --rdata_p4        [必需] Pipeline4 RData\n")
  cat("  --tf_list         [必需] TF 列表 (TSV: TF_ID, Gene_ID, Family)\n")
  cat("  --rdata_p5        [推荐] Pipeline5 Mfuzz RData (基因过滤 + 聚类注释)\n")
  cat("  --out_dir         [可选] 输出目录\n")
  cat("  --lag             [可选] 时滞步数 (默认 1)\n")
  cat("  --n_permutations  [可选] 置换次数 (默认 1000, 仅步骤 7)\n")
  cat("  --fdr_cutoff      [可选] FDR 阈值 (默认 0.05)\n")
  cat("  --target_filter   [可选] 目标基因过滤: all/mfuzz/non_conserved (默认 mfuzz)\n")
  cat("  --prefilter_r     [可选] |r| 预过滤阈值 (默认 0 = 不启用, 推荐 0.3)\n")
  cat("  --output_mode     [可选] significant/classified/full\n")
  cat("  --n_cores         [可选] 并行核数 (默认 1)\n")
  cat("  --force_rerun     [可选] 忽略检查点, 强制重算\n")
  quit(status = 1)
}

################################################################################
# 主程序
################################################################################

cat("================================================================================\n")
cat("  Pipeline 7a: 时滞差异共表达分析\n")
cat("  v1.3 - 解析条件内检验 + Mfuzz过滤 + 检查点缓存\n")
cat("================================================================================\n\n")
cat(sprintf("开始时间: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("并行核数: %d (请求) / %d (可用)\n",
            PARAMS$n_cores, max(1, detectCores(logical = FALSE), na.rm = TRUE)))
cat(sprintf("检查点模式: %s\n",
            ifelse(PARAMS$force_rerun, "强制重算 (忽略缓存)", "启用 (自动跳过已完成步骤)")))
cat(sprintf("目标基因过滤: %s\n", PARAMS$target_filter))
cat(sprintf("|r| 预过滤阈值: %s\n\n",
            ifelse(PARAMS$prefilter_r > 0, sprintf("%.2f", PARAMS$prefilter_r), "不启用")))

# ========================================================================
# 步骤 1: 加载数据
# ========================================================================

cat("--- 步骤 1: 加载数据 ---\n")

load(PARAMS$rdata_p4)
cat(sprintf("  Loaded Pipeline4 RData: %s\n", basename(PARAMS$rdata_p4)))

if (exists("mfuzz_input")) {
  vst_matrix <- mfuzz_input$vst_matrix
  metadata <- mfuzz_input$metadata
  time_responsive_genes <- mfuzz_input$responsive_genes
  cat(sprintf("  Source: mfuzz_input (Pipeline4 Logic1)\n"))
} else if (exists("vsd")) {
  vst_matrix <- assay(vsd)
  if (!exists("metadata")) metadata <- meta
  time_responsive_genes <- rownames(vst_matrix)
  cat(sprintf("  Source: vsd object (Master RData)\n"))
} else {
  stop("RData 格式不匹配: 未找到 mfuzz_input 或 vsd")
}

if (exists("vsd")) {
  vst_matrix_full <- SummarizedExperiment::assay(vsd)
  cat(sprintf("  Full VST matrix available: %d genes x %d samples\n",
              nrow(vst_matrix_full), ncol(vst_matrix_full)))
} else {
  vst_matrix_full <- vst_matrix
  cat(sprintf("  Using subset VST: %d genes x %d samples\n",
              nrow(vst_matrix_full), ncol(vst_matrix_full)))
}

cat(sprintf("  Time-responsive DEGs: %d\n", length(time_responsive_genes)))

if (!is.null(PARAMS$samples_file) && file.exists(PARAMS$samples_file)) {
  metadata <- read.table(PARAMS$samples_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (!"SampleID" %in% colnames(metadata)) colnames(metadata)[1] <- "SampleID"
  rownames(metadata) <- metadata$SampleID
  cat(sprintf("  Metadata overridden from: %s (%d samples)\n",
              basename(PARAMS$samples_file), nrow(metadata)))
}

tf_info <- read.table(PARAMS$tf_list, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "")
cat(sprintf("  TF list loaded: %d entries\n", nrow(tf_info)))

if (!"Gene_ID" %in% colnames(tf_info)) {
  for (alt in c("GeneID", "gene_id", "Gene", "gene")) {
    if (alt %in% colnames(tf_info)) {
      colnames(tf_info)[colnames(tf_info) == alt] <- "Gene_ID"
      break
    }
  }
}
if (!"Family" %in% colnames(tf_info)) {
  for (alt in c("family", "TF_Family", "tf_family")) {
    if (alt %in% colnames(tf_info)) {
      colnames(tf_info)[colnames(tf_info) == alt] <- "Family"
      break
    }
  }
}

tf_info_dedup <- tf_info %>%
  group_by(Gene_ID) %>%
  summarise(Family = first(Family), .groups = "drop") %>%
  as.data.frame()

tf_genes_all <- unique(tf_info_dedup$Gene_ID)
tf_genes <- intersect(tf_genes_all, rownames(vst_matrix_full))
tf_info_dedup <- tf_info_dedup[tf_info_dedup$Gene_ID %in% tf_genes, ]
cat(sprintf("  TFs in dataset: %d / %d\n", length(tf_genes), length(tf_genes_all)))

if (length(tf_genes) == 0) {
  stop("No TF genes found in VST matrix. Check Gene_ID format consistency.")
}

tf_family_map <- setNames(tf_info_dedup$Family, tf_info_dedup$Gene_ID)

target_genes <- time_responsive_genes
cat(sprintf("  Target gene set (all DEGs): %d\n", length(target_genes)))

highlight_map <- NULL
if (!is.null(PARAMS$highlight_file) && file.exists(PARAMS$highlight_file)) {
  highlight_df <- read.table(PARAMS$highlight_file, header = TRUE, sep = "\t",
                              stringsAsFactors = FALSE, quote = "")
  if (!"Gene" %in% colnames(highlight_df)) {
    for (alt in c("GeneID", "Gene_ID", "gene_id", "gene")) {
      if (alt %in% colnames(highlight_df)) {
        colnames(highlight_df)[colnames(highlight_df) == alt] <- "Gene"
        break
      }
    }
  }
  if (!"Symbol" %in% colnames(highlight_df)) {
    for (alt in c("symbol", "Name", "name")) {
      if (alt %in% colnames(highlight_df)) {
        colnames(highlight_df)[colnames(highlight_df) == alt] <- "Symbol"
        break
      }
    }
  }
  highlight_map <- setNames(highlight_df$Symbol, highlight_df$Gene)
  cat(sprintf("  Highlight genes loaded: %d\n", length(highlight_map)))
}

# ---- 加载 Pipeline 5 结果 (聚类注释 + Mfuzz 基因集) ----
cluster_annotation <- NULL
fate_annotation <- NULL
mfuzz_common_genes <- NULL  # v1.3: 用于目标基因过滤
non_conserved_genes <- NULL # v1.3: 非保守基因

if (!is.null(PARAMS$rdata_p5) && file.exists(PARAMS$rdata_p5)) {
  p5_env <- new.env()
  load(PARAMS$rdata_p5, envir = p5_env)
  
  if (exists("wt_clusters", envir = p5_env) && exists("mut_clusters", envir = p5_env)) {
    wt_cl <- get("wt_clusters", envir = p5_env)
    mut_cl <- get("mut_clusters", envir = p5_env)
    
    membership_cutoff <- 0.5
    if (exists("MFUZZ_PARAMS", envir = p5_env)) {
      mp <- get("MFUZZ_PARAMS", envir = p5_env)
      if (!is.null(mp$membership_cutoff)) membership_cutoff <- mp$membership_cutoff
    }
    
    wt_best <- apply(wt_cl$membership, 1, which.max)
    wt_max_mem <- apply(wt_cl$membership, 1, max)
    mut_best <- apply(mut_cl$membership, 1, which.max)
    mut_max_mem <- apply(mut_cl$membership, 1, max)
    
    cluster_annotation <- data.frame(
      Gene_ID = rownames(wt_cl$membership),
      WT_Cluster = wt_best,
      WT_Membership = round(wt_max_mem, 4),
      MUT_Cluster = mut_best[rownames(wt_cl$membership)],
      MUT_Membership = round(mut_max_mem[rownames(wt_cl$membership)], 4),
      stringsAsFactors = FALSE
    )
    cat(sprintf("  Mfuzz cluster annotation loaded: %d genes\n", nrow(cluster_annotation)))
    
    # v1.3: 提取 Mfuzz 共同基因集 (用于目标过滤)
    mfuzz_common_genes <- rownames(wt_cl$membership)
    cat(sprintf("  Mfuzz common genes extracted: %d\n", length(mfuzz_common_genes)))
  }
  
  if (exists("fate_results", envir = p5_env)) {
    fate_res <- get("fate_results", envir = p5_env)
    fate_gene_list <- list()
    non_conserved_list <- c()
    
    for (fr in fate_res) {
      if (is.null(fr$status) || fr$status != "analyzed") next
      wt_cl_id <- fr$wt_cluster
      
      for (k in seq_along(fr$fate_labels)) {
        mem <- fr$subcluster_membership[, k]
        core_genes <- names(which(mem >= membership_cutoff))
        if (length(core_genes) > 0) {
          fate_gene_list[[length(fate_gene_list) + 1]] <- data.frame(
            Gene_ID = core_genes,
            WT_Cluster_P5b = wt_cl_id,
            SubCluster = k,
            Fate_Label = fr$fate_labels[k],
            stringsAsFactors = FALSE
          )
          # v1.3: 收集非保守基因
          if (fr$fate_labels[k] != "Conserved") {
            non_conserved_list <- c(non_conserved_list, core_genes)
          }
        }
      }
    }
    
    if (length(fate_gene_list) > 0) {
      fate_annotation <- do.call(rbind, fate_gene_list)
      fate_annotation <- fate_annotation %>%
        group_by(Gene_ID) %>%
        slice_head(n = 1) %>%
        ungroup() %>%
        as.data.frame()
      cat(sprintf("  Fate annotation loaded: %d genes\n", nrow(fate_annotation)))
    }
    
    non_conserved_genes <- unique(non_conserved_list)
    cat(sprintf("  Non-conserved genes (Phase_Shifted + Shape_Diverged): %d\n",
                length(non_conserved_genes)))
  }
  
  rm(p5_env)
} else {
  if (PARAMS$target_filter != "all") {
    cat("  ⚠ Pipeline5 RData 不可用, 目标基因过滤回退为 'all'\n")
    PARAMS$target_filter <- "all"
  }
}

# ========================================================================
# 步骤 2: 路径与输出准备
# ========================================================================

cat("\n--- 步骤 2: 输出准备 ---\n")

if (is.null(PARAMS$out_dir)) {
  rdata_dir <- dirname(dirname(PARAMS$rdata_p4))
  OUTPUT_DIR <- file.path(rdata_dir, "07_GRN")
} else {
  OUTPUT_DIR <- PARAMS$out_dir
}

dir.create(file.path(OUTPUT_DIR, "Pipeline7a"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUTPUT_DIR, "rdata"), recursive = TRUE, showWarnings = FALSE)

# 检查点目录
CHECKPOINT_DIR <- file.path(OUTPUT_DIR, "Pipeline7a", ".checkpoints")
dir.create(CHECKPOINT_DIR, recursive = TRUE, showWarnings = FALSE)

cat(sprintf("  Output directory: %s\n", OUTPUT_DIR))
cat(sprintf("  Checkpoint directory: %s\n", CHECKPOINT_DIR))

log_file <- file.path(OUTPUT_DIR, "Pipeline7a", "Pipeline7a_diffcoexp.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

# ========================================================================
# 步骤 3: 检测 metadata 列名
# ========================================================================

cat("\n--- 步骤 3: 元数据解析 ---\n")

time_col <- NULL
for (cn in c("Time", "time", "Timepoint", "timepoint", "TIME")) {
  if (cn %in% colnames(metadata)) { time_col <- cn; break }
}
if (is.null(time_col)) stop("Cannot find time column in metadata")

cond_col <- NULL
for (cn in c("Condition", "condition", "Genotype", "genotype")) {
  if (cn %in% colnames(metadata)) { cond_col <- cn; break }
}
if (is.null(cond_col)) stop("Cannot find condition column in metadata")

cat(sprintf("  Time column: %s\n", time_col))
cat(sprintf("  Condition column: %s\n", cond_col))
cat(sprintf("  WT condition: %s\n", PARAMS$condition_wt))
cat(sprintf("  MUT condition: %s\n", PARAMS$condition_mut))

# ========================================================================
# 步骤 4: 构建 Lag-k 矩阵
# ========================================================================

cat(sprintf("\n--- 步骤 4: 构建 Lag-%d Replicate-Aware 矩阵 ---\n", PARAMS$lag))

all_analysis_genes <- union(tf_genes, target_genes)
cat(sprintf("  Total analysis genes (TF U Target): %d\n", length(all_analysis_genes)))

cat(sprintf("\n  [WT] Building lag-%d matrices...\n", PARAMS$lag))
lag_wt <- build_lag_matrices(vst_matrix_full, metadata, PARAMS$condition_wt,
                              all_analysis_genes, time_col, cond_col, PARAMS$lag)
cat(sprintf("    n_eff = %d (reps=%d x lag_pairs=%d)\n",
            lag_wt$n_eff, lag_wt$n_reps, lag_wt$n_lag_pairs))

cat(sprintf("\n  [MUT] Building lag-%d matrices...\n", PARAMS$lag))
lag_mut <- build_lag_matrices(vst_matrix_full, metadata, PARAMS$condition_mut,
                              all_analysis_genes, time_col, cond_col, PARAMS$lag)
cat(sprintf("    n_eff = %d (reps=%d x lag_pairs=%d)\n",
            lag_mut$n_eff, lag_mut$n_reps, lag_mut$n_lag_pairs))

tf_in_data <- intersect(tf_genes, intersect(colnames(lag_wt$design), colnames(lag_mut$design)))
tg_in_data_full <- intersect(target_genes, intersect(colnames(lag_wt$response), colnames(lag_mut$response)))

cat(sprintf("  TFs available in both conditions: %d\n", length(tf_in_data)))
cat(sprintf("  Targets available in both conditions (unfiltered): %d\n", length(tg_in_data_full)))

# ========================================================================
# 步骤 4b: 目标基因过滤  ★ v1.3 NEW ★
# ========================================================================

cat(sprintf("\n--- 步骤 4b: 目标基因过滤 (策略: %s) ---\n", PARAMS$target_filter))

tg_in_data <- tg_in_data_full  # 默认使用全部

if (PARAMS$target_filter == "mfuzz" && !is.null(mfuzz_common_genes)) {
  # 策略 A: Mfuzz 共同基因 (推荐, 平衡覆盖与统计功效)
  tg_in_data <- intersect(tg_in_data_full, mfuzz_common_genes)
  cat(sprintf("  Mfuzz common genes: %d\n", length(mfuzz_common_genes)))
  cat(sprintf("  Targets after Mfuzz filter: %d -> %d (%.1f%% retained)\n",
              length(tg_in_data_full), length(tg_in_data),
              100 * length(tg_in_data) / length(tg_in_data_full)))
  
} else if (PARAMS$target_filter == "non_conserved" && !is.null(non_conserved_genes)) {
  # 策略 B: 仅非保守基因 (最大统计功效, 但覆盖有限)
  tg_in_data <- intersect(tg_in_data_full, non_conserved_genes)
  cat(sprintf("  Non-conserved genes: %d\n", length(non_conserved_genes)))
  cat(sprintf("  Targets after non-conserved filter: %d -> %d (%.1f%% retained)\n",
              length(tg_in_data_full), length(tg_in_data),
              100 * length(tg_in_data) / length(tg_in_data_full)))
  if (length(tg_in_data) < 50) {
    cat("  ⚠ 警告: 非保守基因数过少, 可能遗漏重要调控关系. 建议使用 --target_filter mfuzz\n")
  }
  
} else {
  # 策略 C: 使用全部 DEGs (原 v1.2 行为)
  cat(sprintf("  Using all DEGs as targets: %d\n", length(tg_in_data)))
  if (PARAMS$target_filter != "all") {
    cat(sprintf("  ⚠ 请求的过滤策略 '%s' 不可用 (缺少 P5 数据), 回退为 'all'\n",
                PARAMS$target_filter))
  }
}

n_total_pairs <- length(tf_in_data) * length(tg_in_data) - length(intersect(tf_in_data, tg_in_data))
cat(sprintf("\n  最终检验规模:\n"))
cat(sprintf("    TFs: %d\n", length(tf_in_data)))
cat(sprintf("    Targets: %d\n", length(tg_in_data)))
cat(sprintf("    Total pairs (excl. self-loops): %d\n", n_total_pairs))
cat(sprintf("    BH rank-1 门槛 (FDR<0.05): p < %.2e\n", 0.05 / n_total_pairs))

# 预估解析 t-test 所需最小 |r|
df_wt <- lag_wt$n_eff - 2
if (df_wt >= 1) {
  # 使用 qt 反推: 对于 rank-1 的 p 值门槛, 需要多大的 t 值
  p_threshold_rank1 <- 0.05 / n_total_pairs
  t_required <- qt(1 - p_threshold_rank1 / 2, df = df_wt)
  # t = r * sqrt(df / (1 - r^2)) → r = t / sqrt(t^2 + df)
  r_required <- t_required / sqrt(t_required^2 + df_wt)
  cat(sprintf("    解析 t-test (df=%d): rank-1 需要 |r| > %.3f, rank-100 需要 |r| > %.3f\n",
              df_wt, r_required,
              {
                p100 <- 0.05 * 100 / n_total_pairs
                t100 <- qt(1 - p100 / 2, df = df_wt)
                t100 / sqrt(t100^2 + df_wt)
              }))
}

# ========================================================================
# 步骤 5: 计算观测相关矩阵
# ========================================================================

cat(sprintf("\n--- 步骤 5: 计算观测 Spearman Lag-%d 相关 ---\n", PARAMS$lag))

r_wt <- calc_cor_matrix(
  lag_wt$design[, tf_in_data, drop = FALSE],
  lag_wt$response[, tg_in_data, drop = FALSE],
  PARAMS$cor_method
)
cat(sprintf("  WT correlation matrix: %d x %d\n", nrow(r_wt), ncol(r_wt)))
cat(sprintf("  WT |r| distribution: mean=%.3f, sd=%.3f, max=%.3f\n",
            mean(abs(r_wt), na.rm = TRUE),
            sd(abs(r_wt), na.rm = TRUE),
            max(abs(r_wt), na.rm = TRUE)))

r_mut <- calc_cor_matrix(
  lag_mut$design[, tf_in_data, drop = FALSE],
  lag_mut$response[, tg_in_data, drop = FALSE],
  PARAMS$cor_method
)
cat(sprintf("  MUT correlation matrix: %d x %d\n", nrow(r_mut), ncol(r_mut)))
cat(sprintf("  MUT |r| distribution: mean=%.3f, sd=%.3f, max=%.3f\n",
            mean(abs(r_mut), na.rm = TRUE),
            sd(abs(r_mut), na.rm = TRUE),
            max(abs(r_mut), na.rm = TRUE)))

clip_r <- function(r) pmin(pmax(r, -0.9999), 0.9999)
z_wt <- atanh(clip_r(r_wt))
z_mut <- atanh(clip_r(r_mut))

se_diff <- sqrt(1 / max(1, lag_wt$n_eff - 3) + 1 / max(1, lag_mut$n_eff - 3))
z_diff <- (z_wt - z_mut) / se_diff

cat(sprintf("  Fisher Z SE: %.4f (n_wt=%d, n_mut=%d)\n", se_diff, lag_wt$n_eff, lag_mut$n_eff))

# ========================================================================
# 步骤 6: 条件内解析 p 值  ★ v1.3 CRITICAL FIX ★
# ========================================================================
#
# v1.2 问题诊断:
#   置换参数化 null: mean(|r|)=0.2543, sd=0.1809
#   即使 |r_obs|=0.95: z=(0.95-0.2543)/0.1809=3.84 → p=1.22e-04
#   p_floor 截断至 1e-04
#   BH rank-1 需要 p < 2.19e-08 → FDR=280 ≫ 0.05
#   → 数学上不可能产生显著结果
#
# v1.3 解决方案:
#   解析 t-test: t = r × sqrt(df/(1-r²)), df = n_eff - 2
#   |r|=0.95 → t = 0.95 × sqrt(10/0.0975) = 9.62 → p = 2.2e-06
#   连续 p 值, 无下限限制, 计算成本 O(1)
#
# ========================================================================

cat(sprintf("\n--- 步骤 6: 条件内解析显著性检验 (t-approximation, df=%d) ---\n", df_wt))
cat(sprintf("  方法: Spearman rho t-approximation (替代 v1.2 置换检验)\n"))
cat(sprintf("  优势: 连续 p 值 (无 1/B 下限), 计算 O(1), 无需并行\n"))

t_step6_start <- Sys.time()

# WT 解析 p 值
cat("\n  [WT] Analytical t-test...\n")
analytical_wt <- calc_analytical_pvalues(r_wt, lag_wt$n_eff)
n_sig_wt <- sum(analytical_wt$fdr_values < PARAMS$fdr_cutoff, na.rm = TRUE)
cat(sprintf("    n_eff=%d, df=%d\n", lag_wt$n_eff, analytical_wt$df))
cat(sprintf("    p-value range: [%.2e, %.4f]\n",
            min(analytical_wt$p_values, na.rm = TRUE),
            max(analytical_wt$p_values, na.rm = TRUE)))
cat(sprintf("    Significant edges (FDR<%.3f): %d\n", PARAMS$fdr_cutoff, n_sig_wt))

# MUT 解析 p 值
cat("\n  [MUT] Analytical t-test...\n")
analytical_mut <- calc_analytical_pvalues(r_mut, lag_mut$n_eff)
n_sig_mut <- sum(analytical_mut$fdr_values < PARAMS$fdr_cutoff, na.rm = TRUE)
cat(sprintf("    n_eff=%d, df=%d\n", lag_mut$n_eff, analytical_mut$df))
cat(sprintf("    p-value range: [%.2e, %.4f]\n",
            min(analytical_mut$p_values, na.rm = TRUE),
            max(analytical_mut$p_values, na.rm = TRUE)))
cat(sprintf("    Significant edges (FDR<%.3f): %d\n", PARAMS$fdr_cutoff, n_sig_mut))

# 封装为兼容结构 (与 v1.2 perm_wt/perm_mut 接口一致)
perm_wt <- list(
  p_values      = analytical_wt$p_values,
  fdr_values    = analytical_wt$fdr_values,
  p_empirical   = analytical_wt$p_values,   # 解析方法无经验 p 值, 复制
  fdr_empirical = analytical_wt$fdr_values,  # 同上
  null_mean     = matrix(0, nrow = nrow(r_wt), ncol = ncol(r_wt),
                         dimnames = dimnames(r_wt)),
  null_sd       = matrix(1, nrow = nrow(r_wt), ncol = ncol(r_wt),
                         dimnames = dimnames(r_wt)),
  method        = "analytical_t_test",
  df            = analytical_wt$df
)

perm_mut <- list(
  p_values      = analytical_mut$p_values,
  fdr_values    = analytical_mut$fdr_values,
  p_empirical   = analytical_mut$p_values,
  fdr_empirical = analytical_mut$fdr_values,
  null_mean     = matrix(0, nrow = nrow(r_mut), ncol = ncol(r_mut),
                         dimnames = dimnames(r_mut)),
  null_sd       = matrix(1, nrow = nrow(r_mut), ncol = ncol(r_mut),
                         dimnames = dimnames(r_mut)),
  method        = "analytical_t_test",
  df            = analytical_mut$df
)

t_step6_elapsed <- as.numeric(difftime(Sys.time(), t_step6_start, units = "secs"))
cat(sprintf("\n  步骤 6 完成: %.1f 秒 (解析方法, 无需置换)\n", t_step6_elapsed))

# ========================================================================
# 步骤 6b: 可选 |r| 预过滤 + 重新 BH 校正  ★ v1.3 NEW ★
# ========================================================================

if (PARAMS$prefilter_r > 0) {
  cat(sprintf("\n--- 步骤 6b: |r| 预过滤 (阈值=%.2f) + 重新 BH 校正 ---\n",
              PARAMS$prefilter_r))
  
  r_max <- pmax(abs(r_wt), abs(r_mut), na.rm = TRUE)
  keep_mask <- r_max >= PARAMS$prefilter_r
  n_kept <- sum(keep_mask, na.rm = TRUE)
  n_total <- length(r_max)
  
  cat(sprintf("  Pre-filter: %d / %d pairs with max(|r_WT|, |r_MUT|) >= %.2f (%.1f%%)\n",
              n_kept, n_total, PARAMS$prefilter_r, 100 * n_kept / n_total))
  
  if (n_kept > 0 && n_kept < n_total) {
    # WT: 仅对过滤后子集做 BH
    fdr_wt_filtered <- matrix(1, nrow = nrow(analytical_wt$p_values),
                               ncol = ncol(analytical_wt$p_values),
                               dimnames = dimnames(analytical_wt$p_values))
    fdr_wt_filtered[keep_mask] <- p.adjust(analytical_wt$p_values[keep_mask], method = "BH")
    
    # MUT
    fdr_mut_filtered <- matrix(1, nrow = nrow(analytical_mut$p_values),
                                ncol = ncol(analytical_mut$p_values),
                                dimnames = dimnames(analytical_mut$p_values))
    fdr_mut_filtered[keep_mask] <- p.adjust(analytical_mut$p_values[keep_mask], method = "BH")
    
    # 更新 perm 对象
    perm_wt$fdr_values  <- fdr_wt_filtered
    perm_mut$fdr_values <- fdr_mut_filtered
    
    n_sig_wt_filtered <- sum(fdr_wt_filtered < PARAMS$fdr_cutoff, na.rm = TRUE)
    n_sig_mut_filtered <- sum(fdr_mut_filtered < PARAMS$fdr_cutoff, na.rm = TRUE)
    
    cat(sprintf("  After pre-filter BH:\n"))
    cat(sprintf("    WT sig edges: %d -> %d\n", n_sig_wt, n_sig_wt_filtered))
    cat(sprintf("    MUT sig edges: %d -> %d\n", n_sig_mut, n_sig_mut_filtered))
    cat(sprintf("    BH rank-1 门槛: %.2e -> %.2e (%.1fx 更宽松)\n",
                0.05 / n_total, 0.05 / n_kept,
                n_total / n_kept))
    
    # 保存预过滤 mask 供步骤 7 使用
    prefilter_mask <- keep_mask
  } else {
    cat("  ⚠ 预过滤未改变基因对数, 跳过\n")
    prefilter_mask <- NULL
  }
} else {
  prefilter_mask <- NULL
}

# ========================================================================
# 步骤 7: 跨条件差异置换检验  ★ 检查点支持 ★
# ========================================================================

cat(sprintf("\n--- 步骤 7: 跨条件差异置换检验 (B=%d, cores=%d) ---\n",
            PARAMS$n_permutations, PARAMS$n_cores))

# v1.3: 检查点文件名包含目标基因数以区分不同过滤策略
checkpoint_tag <- sprintf("step7_%s_%dtg", PARAMS$target_filter, length(tg_in_data))
checkpoint_step7 <- file.path(CHECKPOINT_DIR, sprintf("Pipeline7a_checkpoint_%s.RData", checkpoint_tag))

if (!PARAMS$force_rerun && file.exists(checkpoint_step7)) {
  cat(sprintf("  ★ 检测到步骤 7 检查点 (%s), 加载缓存结果...\n", checkpoint_tag))
  load(checkpoint_step7)
  cat(sprintf("    Differential: %d sig edges (parametric FDR<0.05)\n",
              sum(perm_diff$fdr_values < 0.05, na.rm = TRUE)))
  cat("  ★ 步骤 7 跳过 (使用缓存). 如需重算请删除检查点文件或使用 --force_rerun\n")
} else {
  perm_diff <- run_differential_permutation(
    vst_matrix    = vst_matrix_full,
    metadata      = metadata,
    tf_genes      = tf_in_data,
    target_genes  = tg_in_data,
    n_eff_wt      = lag_wt$n_eff,
    n_eff_mut     = lag_mut$n_eff,
    z_diff_obs    = z_diff,
    condition_wt  = PARAMS$condition_wt,
    condition_mut = PARAMS$condition_mut,
    time_col      = time_col,
    cond_col      = cond_col,
    lag           = PARAMS$lag,
    B             = PARAMS$n_permutations,
    method        = PARAMS$cor_method,
    seed          = PARAMS$seed,
    n_cores       = PARAMS$n_cores,
    verbose       = PARAMS$verbose
  )
  
  cat(sprintf("\n  ★ 保存步骤 7 检查点: %s\n", basename(checkpoint_step7)))
  save(perm_diff, file = checkpoint_step7)
}

# v1.3: 如果启用了 |r| 预过滤, 对差异 p 值也重新做 BH
if (!is.null(prefilter_mask)) {
  cat("\n  步骤 7b: 对差异 p 值应用预过滤 BH 校正...\n")
  fdr_diff_filtered <- matrix(1, nrow = nrow(perm_diff$p_values),
                               ncol = ncol(perm_diff$p_values),
                               dimnames = dimnames(perm_diff$p_values))
  fdr_diff_filtered[prefilter_mask] <- p.adjust(perm_diff$p_values[prefilter_mask], method = "BH")
  
  n_sig_diff_before <- sum(perm_diff$fdr_values < PARAMS$fdr_cutoff_diff, na.rm = TRUE)
  n_sig_diff_after  <- sum(fdr_diff_filtered < PARAMS$fdr_cutoff_diff, na.rm = TRUE)
  
  perm_diff$fdr_values <- fdr_diff_filtered
  
  cat(sprintf("    Diff sig edges: %d -> %d\n", n_sig_diff_before, n_sig_diff_after))
}

# ========================================================================
# 步骤 8: 构建边表
# ========================================================================

cat("\n--- 步骤 8: 构建边表与边分类 ---\n")

edge_table <- expand.grid(
  TF_Gene = tf_in_data,
  Target_Gene = tg_in_data,
  stringsAsFactors = FALSE
)

edge_table <- edge_table[edge_table$TF_Gene != edge_table$Target_Gene, ]
cat(sprintf("  Total edges (excl. self-loops): %d\n", nrow(edge_table)))

get_mat_val <- function(mat, tf_col, tg_col) {
  mapply(function(tf, tg) {
    if (tf %in% rownames(mat) && tg %in% colnames(mat)) {
      mat[tf, tg]
    } else {
      NA_real_
    }
  }, tf_col, tg_col, USE.NAMES = FALSE)
}

edge_table$r_WT     <- get_mat_val(r_wt, edge_table$TF_Gene, edge_table$Target_Gene)
edge_table$r_MUT    <- get_mat_val(r_mut, edge_table$TF_Gene, edge_table$Target_Gene)
edge_table$delta_r  <- edge_table$r_WT - edge_table$r_MUT
edge_table$abs_delta_r <- abs(edge_table$delta_r)
edge_table$z_WT     <- get_mat_val(z_wt, edge_table$TF_Gene, edge_table$Target_Gene)
edge_table$z_MUT    <- get_mat_val(z_mut, edge_table$TF_Gene, edge_table$Target_Gene)
edge_table$z_diff   <- get_mat_val(z_diff, edge_table$TF_Gene, edge_table$Target_Gene)

# 条件内 p 值 (v1.3: 解析 t-test)
edge_table$p_WT     <- get_mat_val(perm_wt$p_values, edge_table$TF_Gene, edge_table$Target_Gene)
edge_table$p_MUT    <- get_mat_val(perm_mut$p_values, edge_table$TF_Gene, edge_table$Target_Gene)
edge_table$p_diff   <- get_mat_val(perm_diff$p_values, edge_table$TF_Gene, edge_table$Target_Gene)

edge_table$fdr_WT   <- get_mat_val(perm_wt$fdr_values, edge_table$TF_Gene, edge_table$Target_Gene)
edge_table$fdr_MUT  <- get_mat_val(perm_mut$fdr_values, edge_table$TF_Gene, edge_table$Target_Gene)
edge_table$fdr_diff <- get_mat_val(perm_diff$fdr_values, edge_table$TF_Gene, edge_table$Target_Gene)

# 差异经验 p 值 (诊断列, 仅步骤 7 有)
edge_table$p_diff_empirical <- get_mat_val(perm_diff$p_empirical, edge_table$TF_Gene, edge_table$Target_Gene)

edge_table$sig_WT   <- edge_table$fdr_WT < PARAMS$fdr_cutoff
edge_table$sig_MUT  <- edge_table$fdr_MUT < PARAMS$fdr_cutoff
edge_table$sig_diff <- edge_table$fdr_diff < PARAMS$fdr_cutoff_diff

sig_mut_loose <- edge_table$fdr_MUT < PARAMS$fdr_mut_loose

edge_table$edge_class <- classify_edges(
  sig_wt = edge_table$sig_WT,
  sig_mut = edge_table$sig_MUT,
  sig_diff = edge_table$sig_diff,
  sig_mut_loose = sig_mut_loose
)

edge_table$direction_WT <- ifelse(
  edge_table$sig_WT,
  ifelse(edge_table$r_WT > 0, "positive", "negative"),
  "NS"
)
edge_table$direction_MUT <- ifelse(
  edge_table$sig_MUT,
  ifelse(edge_table$r_MUT > 0, "positive", "negative"),
  "NS"
)

edge_table$TF_Family <- tf_family_map[edge_table$TF_Gene]

if (!is.null(highlight_map)) {
  edge_table$TF_Symbol <- highlight_map[edge_table$TF_Gene]
  edge_table$Target_Symbol <- highlight_map[edge_table$Target_Gene]
} else {
  edge_table$TF_Symbol <- NA_character_
  edge_table$Target_Symbol <- NA_character_
}

if (!is.null(cluster_annotation)) {
  edge_table <- edge_table %>%
    left_join(cluster_annotation %>% select(Gene_ID, WT_Cluster, MUT_Cluster) %>%
                rename(Target_WT_Cluster = WT_Cluster, Target_MUT_Cluster = MUT_Cluster),
              by = c("Target_Gene" = "Gene_ID"))
} else {
  edge_table$Target_WT_Cluster <- NA_integer_
  edge_table$Target_MUT_Cluster <- NA_integer_
}

if (!is.null(fate_annotation)) {
  edge_table <- edge_table %>%
    left_join(fate_annotation %>% select(Gene_ID, Fate_Label) %>%
                rename(Target_Fate = Fate_Label),
              by = c("Target_Gene" = "Gene_ID"))
} else {
  edge_table$Target_Fate <- NA_character_
}

col_order <- c(
  "TF_Gene", "Target_Gene", "TF_Family", "TF_Symbol", "Target_Symbol",
  "edge_class",
  "r_WT", "r_MUT", "delta_r", "abs_delta_r",
  "z_WT", "z_MUT", "z_diff",
  "p_WT", "p_MUT", "p_diff",
  "fdr_WT", "fdr_MUT", "fdr_diff",
  "p_diff_empirical",
  "sig_WT", "sig_MUT", "sig_diff",
  "direction_WT", "direction_MUT",
  "Target_WT_Cluster", "Target_MUT_Cluster", "Target_Fate"
)
col_order <- intersect(col_order, colnames(edge_table))
extra_cols <- setdiff(colnames(edge_table), col_order)
edge_table <- edge_table[, c(col_order, extra_cols)]

# ========================================================================
# 步骤 9: 输出
# ========================================================================

cat("\n--- 步骤 9: 输出结果 ---\n")

if (PARAMS$output_mode == "full") {
  edge_output <- edge_table
} else if (PARAMS$output_mode == "classified") {
  edge_output <- edge_table %>% filter(edge_class != "NS")
} else {
  edge_output <- edge_table %>% filter(sig_WT | sig_MUT)
}

edge_file <- file.path(OUTPUT_DIR, "Pipeline7a", "Pipeline7a_DiffCoexp_Edges.tsv")
write.table(edge_output, edge_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  Edge table: %s (%d edges, mode=%s)\n",
            basename(edge_file), nrow(edge_output), PARAMS$output_mode))

tf_summary <- edge_table %>%
  filter(sig_WT | sig_MUT) %>%
  group_by(TF_Gene, TF_Family) %>%
  summarise(
    n_targets_WT        = sum(sig_WT, na.rm = TRUE),
    n_targets_MUT       = sum(sig_MUT, na.rm = TRUE),
    n_lost              = sum(edge_class == "Lost", na.rm = TRUE),
    n_gained            = sum(edge_class == "Gained", na.rm = TRUE),
    n_rewired           = sum(edge_class == "Rewired", na.rm = TRUE),
    n_conserved         = sum(edge_class == "Conserved", na.rm = TRUE),
    delta_out_degree    = sum(sig_WT, na.rm = TRUE) - sum(sig_MUT, na.rm = TRUE),
    mean_r_WT           = mean(r_WT[sig_WT], na.rm = TRUE),
    mean_r_MUT          = mean(r_MUT[sig_MUT], na.rm = TRUE),
    mean_abs_delta_r    = mean(abs_delta_r[sig_WT | sig_MUT], na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(abs(delta_out_degree)))

if (!is.null(highlight_map)) {
  tf_summary$TF_Symbol <- highlight_map[tf_summary$TF_Gene]
}

tf_file <- file.path(OUTPUT_DIR, "Pipeline7a", "Pipeline7a_TF_Summary.tsv")
write.table(tf_summary, tf_file, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("  TF summary: %s (%d TFs)\n", basename(tf_file), nrow(tf_summary)))

class_counts <- table(edge_table$edge_class)
class_df <- data.frame(
  Edge_Class = names(class_counts),
  Count = as.integer(class_counts),
  stringsAsFactors = FALSE
)

summary_lines <- c(
  "================================================================================",
  "  Pipeline 7a: Differential Co-expression Analysis Summary (v1.3)",
  "================================================================================",
  "",
  sprintf("Analysis Time: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("Lag: %d", PARAMS$lag),
  sprintf("Correlation Method: %s", PARAMS$cor_method),
  "",
  "--- v1.3 Method Configuration ---",
  sprintf("  Within-condition significance: analytical t-test (df=%d)", df_wt),
  sprintf("  Differential significance: parametric permutation (B=%d)", PARAMS$n_permutations),
  sprintf("  Target gene filter: %s", PARAMS$target_filter),
  sprintf("  |r| pre-filter: %s", ifelse(PARAMS$prefilter_r > 0,
                                          sprintf("%.2f", PARAMS$prefilter_r), "disabled")),
  sprintf("  Parallel Cores: %d (differential permutation only)", PARAMS$n_cores),
  sprintf("  FDR Cutoff (WT/MUT sig): %.3f", PARAMS$fdr_cutoff),
  sprintf("  FDR Cutoff (differential): %.3f", PARAMS$fdr_cutoff_diff),
  sprintf("  FDR MUT Loose (Lost Edge): %.3f", PARAMS$fdr_mut_loose),
  "",
  "--- Data Dimensions ---",
  sprintf("  TFs tested: %d", length(tf_in_data)),
  sprintf("  Targets (before filter): %d", length(tg_in_data_full)),
  sprintf("  Targets (after filter): %d", length(tg_in_data)),
  sprintf("  Total pairs (excl. self-loops): %d", nrow(edge_table)),
  sprintf("  n_eff WT: %d (reps=%d x lag_pairs=%d)", lag_wt$n_eff, lag_wt$n_reps, lag_wt$n_lag_pairs),
  sprintf("  n_eff MUT: %d (reps=%d x lag_pairs=%d)", lag_mut$n_eff, lag_mut$n_reps, lag_mut$n_lag_pairs),
  "",
  "--- Significance Summary ---",
  sprintf("  Significant WT edges (FDR < %.3f): %d", PARAMS$fdr_cutoff,
          sum(edge_table$sig_WT, na.rm = TRUE)),
  sprintf("  Significant MUT edges (FDR < %.3f): %d", PARAMS$fdr_cutoff,
          sum(edge_table$sig_MUT, na.rm = TRUE)),
  sprintf("  Significant differential (FDR < %.3f): %d", PARAMS$fdr_cutoff_diff,
          sum(edge_table$sig_diff, na.rm = TRUE)),
  "",
  "--- Edge Classification ---"
)
for (i in 1:nrow(class_df)) {
  summary_lines <- c(summary_lines,
                     sprintf("  %s: %d", class_df$Edge_Class[i], class_df$Count[i]))
}

# v1.3 诊断信息
summary_lines <- c(summary_lines,
  "",
  "--- v1.3 Diagnostic: Analytical p-value Statistics ---",
  sprintf("  WT p-value range:  [%.2e, %.4f]",
          min(perm_wt$p_values, na.rm = TRUE),
          max(perm_wt$p_values, na.rm = TRUE)),
  sprintf("  MUT p-value range: [%.2e, %.4f]",
          min(perm_mut$p_values, na.rm = TRUE),
          max(perm_mut$p_values, na.rm = TRUE)),
  sprintf("  Diff p-value range (parametric): [%.2e, %.4f]",
          min(perm_diff$p_values, na.rm = TRUE),
          max(perm_diff$p_values, na.rm = TRUE)),
  sprintf("  BH correction factor: %d tests", nrow(edge_table)),
  sprintf("  Required raw p for FDR<0.05 (approx): < %.2e", 0.05 / nrow(edge_table)),
  "",
  "--- v1.3 vs v1.2 Comparison ---",
  sprintf("  v1.2: Step 6 method = permutation + parametric null (p_floor=1e-04)"),
  sprintf("  v1.3: Step 6 method = analytical t-test (continuous p, no floor)"),
  sprintf("  v1.2: Targets = %d (all DEGs)", length(tg_in_data_full)),
  sprintf("  v1.3: Targets = %d (%s filter)", length(tg_in_data), PARAMS$target_filter)
)

if (!is.null(prefilter_mask)) {
  summary_lines <- c(summary_lines,
    sprintf("  v1.3: Pre-filter pairs: %d (|r| >= %.2f)",
            sum(prefilter_mask, na.rm = TRUE), PARAMS$prefilter_r))
}

summary_file <- file.path(OUTPUT_DIR, "Pipeline7a", "Pipeline7a_Summary.txt")
writeLines(summary_lines, summary_file)
cat(sprintf("  Summary: %s\n", basename(summary_file)))

cat("\n")
for (line in summary_lines) cat(line, "\n")

# RData
save(
  r_wt, r_mut, z_wt, z_mut, z_diff,
  perm_wt, perm_mut, perm_diff,
  edge_table, edge_output,
  tf_summary,
  tf_in_data, tg_in_data, tg_in_data_full, tf_family_map,
  lag_wt, lag_mut,
  cluster_annotation, fate_annotation, highlight_map,
  mfuzz_common_genes, non_conserved_genes,
  PARAMS,
  file = file.path(OUTPUT_DIR, "rdata", "Pipeline7a_DiffCoexp_Results.RData")
)
cat(sprintf("  RData: %s\n", "Pipeline7a_DiffCoexp_Results.RData"))

# Parameters JSON
params_json <- list(
  pipeline = "7a",
  version = "1.3",
  description = "Lag-k differential co-expression with analytical within-condition tests + Mfuzz filtering",
  method = list(
    correlation = PARAMS$cor_method,
    lag = PARAMS$lag,
    within_condition = sprintf("Analytical t-test (t = r*sqrt(df/(1-r^2)), df=%d)", df_wt),
    differential = sprintf("Cross-condition permutation with parametric null fitting (B=%d)",
                          PARAMS$n_permutations),
    fdr_correction = "Benjamini-Hochberg",
    target_filter = PARAMS$target_filter,
    prefilter_r = PARAMS$prefilter_r,
    replicate_handling = "Replicate-Aware matrix concatenation",
    parallelization = sprintf("mclapply with %d cores (differential permutation only)", PARAMS$n_cores),
    checkpoint = "Step 7 cached in .checkpoints/ directory"
  ),
  parameters = list(
    n_permutations = PARAMS$n_permutations,
    fdr_cutoff = PARAMS$fdr_cutoff,
    fdr_cutoff_diff = PARAMS$fdr_cutoff_diff,
    fdr_mut_loose = PARAMS$fdr_mut_loose,
    seed = PARAMS$seed,
    n_cores = PARAMS$n_cores,
    output_mode = PARAMS$output_mode,
    target_filter = PARAMS$target_filter,
    prefilter_r = PARAMS$prefilter_r
  ),
  data_dimensions = list(
    n_tf = length(tf_in_data),
    n_targets_before_filter = length(tg_in_data_full),
    n_targets_after_filter = length(tg_in_data),
    n_eff_wt = lag_wt$n_eff,
    n_eff_mut = lag_mut$n_eff,
    n_timepoints = lag_wt$n_times,
    analytical_df = df_wt
  ),
  edge_classification = list(
    Lost = "sig_WT & sig_diff & !sig_MUT_loose",
    Gained = "!sig_WT & sig_diff & sig_MUT",
    Rewired = "sig_WT & sig_diff & sig_MUT_loose",
    Conserved = "sig_WT & sig_MUT & !sig_diff",
    WT_only = "sig_WT & !sig_MUT & !sig_diff",
    MUT_only = "!sig_WT & sig_MUT & !sig_diff",
    NS = "neither significant"
  ),
  v13_changes = list(
    step6 = "Replaced permutation test with analytical t-test (eliminates p_floor bottleneck)",
    step4b = "Added Mfuzz common gene filtering (reduces multiple testing burden ~3.6x)",
    step6b = "Optional |r| pre-filtering for further BH correction reduction",
    step7 = "Removed p_floor from differential permutation (allows more precise parametric p)"
  ),
  input_files = list(
    rdata_p4 = PARAMS$rdata_p4,
    rdata_p5 = PARAMS$rdata_p5,
    tf_list = PARAMS$tf_list,
    highlight_file = PARAMS$highlight_file
  ),
  analysis_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

writeLines(
  jsonlite::toJSON(params_json, pretty = TRUE, auto_unbox = TRUE),
  file.path(OUTPUT_DIR, "Pipeline7a", "Pipeline7a_Parameters.json")
)
cat(sprintf("  Parameters JSON: %s\n", "Pipeline7a_Parameters.json"))

################################################################################
# 完成
################################################################################

cat("\n================================================================================\n")
cat("  Pipeline 7a 完成 (v1.3 解析条件内检验 + Mfuzz过滤)\n")
cat("================================================================================\n\n")
cat(sprintf("结束时间: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

cat("输出文件:\n")
cat(sprintf("  %s/Pipeline7a/\n", OUTPUT_DIR))
cat("    ├── Pipeline7a_DiffCoexp_Edges.tsv   - 差异共表达边表 (核心输出)\n")
cat("    ├── Pipeline7a_TF_Summary.tsv        - TF 级统计摘要\n")
cat("    ├── Pipeline7a_Summary.txt           - 运行摘要\n")
cat("    ├── Pipeline7a_Parameters.json       - 参数记录\n")
cat("    ├── Pipeline7a_diffcoexp.log         - 完整日志\n")
cat("    └── .checkpoints/                    - 步骤7置换检验缓存\n")
cat(sprintf("  %s/rdata/\n", OUTPUT_DIR))
cat("    └── Pipeline7a_DiffCoexp_Results.RData - 完整中间结果\n")

cat("\nv1.3 关键变更:\n")
cat("  1. 步骤 6: 置换检验 → 解析 t-test (连续 p 值, <1 秒)\n")
cat(sprintf("  2. 步骤 4b: 目标基因过滤 (%s): %d → %d 基因\n",
            PARAMS$target_filter, length(tg_in_data_full), length(tg_in_data)))
if (PARAMS$prefilter_r > 0) {
  cat(sprintf("  3. 步骤 6b: |r| >= %.2f 预过滤 (进一步减少 BH 负担)\n",
              PARAMS$prefilter_r))
}
cat("  4. 步骤 7: 移除 p_floor, 保留参数化置换 (已验证有效)\n")

cat("\n检查点说明:\n")
cat("  - 步骤 6 不再需要检查点 (解析计算 <1 秒)\n")
cat("  - 步骤 7 置换检验结果缓存至 .checkpoints/ 目录\n")
cat("  - 检查点文件名包含过滤策略标签, 不同策略不会互相覆盖\n")
cat("  - 如需强制重算: 删除 .checkpoints/ 下的文件, 或使用 --force_rerun\n")

cat("\n下游使用指南:\n")
cat("  边表列说明:\n")
cat("    edge_class    = Lost/Gained/Rewired/Conserved/WT_only/MUT_only/NS\n")
cat("    p_WT/p_MUT    = 解析 t-test p 值 (v1.3, 连续分布)\n")
cat("    p_diff        = 参数化置换 p 值 (v1.2 方法)\n")
cat("    p_diff_empirical = 经验 p 值 (离散, 用于诊断对比)\n")
cat("    Target_Fate   = Pipeline5b 命运分化标签 (若加载了 P5 RData)\n")
cat("  过滤示例 (在 7b/7c 中):\n")
cat("    Lost Edges:      filter(edge_class == 'Lost')\n")
cat("    Fate-filtered:   filter(edge_class == 'Lost', Target_Fate != 'Conserved')\n")
cat("    High-confidence: filter(edge_class == 'Lost', fdr_diff < 0.01)\n")

sink(type = "message")
sink(type = "output")
close(log_con)

cat("\n日志已保存至: ", log_file, "\n")
