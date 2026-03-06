#!/usr/bin/env Rscript
################################################################################
# RNAseq Pipeline 6e v1.0: Permutation Baseline Test
#
# 功能特性：
#   - 标签置换检验：验证 6b 报告的 AUC 是否显著优于随机预期
#   - 每次置换独立采样负样本 + 打乱正/负标签 → 简化 RF CV → 记录 null AUC
#   - 经验 p-value 计算
#   - 自动生成 Null Distribution 可视化 (直方图 + 真实 AUC 标注)
#   - 参数化控制：置换次数、CV 规格、核心数均通过命令行/Config 配置
#
# 输入:
#   - 6a 的 ML_ready 数据 (Pipeline6b_Input_Data.rds)
#   - 6b 的 RF 结果 (Pipeline6b_RF_Results.rds) → 提取 real AUC
#
# 输出:
#   - Permutation_Summary.txt         (per-cluster real AUC + null stats + p-value)
#   - Permutation_NullDistributions.rds (完整 null 分布，供后续分析)
#   - Fig_Permutation_Test_*.pdf      (per-genotype 可视化)
#
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(randomForest)
  library(pROC)
  library(doParallel)
  library(foreach)
  library(optparse)
})

################################################################################
# 命令行参数
################################################################################

option_list <- list(
  # === 输入输出 ===
  make_option(c("--p6a_dir"), type = "character", default = NULL,
              help = "Pipeline 6a 输出目录 (含 ML_ready 数据)"),
  make_option(c("--p6b_dir"), type = "character", default = NULL,
              help = "Pipeline 6b 输出目录 (含 RF_Results.rds)"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "输出目录 [default: auto]"),
  make_option(c("--genotype"), type = "character", default = "both",
              help = "处理的基因型: 'WT', 'Mutant', 或 'both' [default: %default]"),
  
  # === 置换参数 ===
  make_option(c("--n_perm"), type = "integer", default = 100,
              help = "置换次数 [default: %default]"),
  make_option(c("--n_replicates"), type = "integer", default = 10,
              help = "每次置换的 CV 重复次数 (简化版) [default: %default]"),
  make_option(c("--n_folds"), type = "integer", default = 5,
              help = "每次重复的 CV 折数 (简化版) [default: %default]"),
  make_option(c("--ntree"), type = "integer", default = 200,
              help = "置换 RF 树数量 (可小于 6b) [default: %default]"),
  
  # === 负样本策略 (与 6b 保持一致) ===
  make_option(c("--gc_match"), action = "store_true", default = FALSE,
              help = "是否对负样本进行 GC 匹配采样"),
  make_option(c("--gc_bins"), type = "integer", default = 10,
              help = "GC 匹配的分层数量 [default: %default]"),
  make_option(c("--neg_ratio"), type = "double", default = 3,
              help = "负样本与正样本的比例 [default: %default]"),
  
  # === 资源 ===
  make_option(c("--cores"), type = "integer", default = 40,
              help = "并行核心数 [default: %default]"),
  make_option(c("--seed"), type = "integer", default = 42,
              help = "随机种子 [default: %default]"),
  
  # === 可视化 ===
  make_option(c("--fig_width"), type = "double", default = 10,
              help = "Figure 宽度 (inches) [default: %default]"),
  make_option(c("--fig_height"), type = "double", default = 6,
              help = "Figure 高度 (inches) [default: %default]"),
  make_option(c("--alpha"), type = "double", default = 0.05,
              help = "显著性阈值 [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 默认输出路径
if (is.null(opt$output_dir)) {
  opt$output_dir <- file.path(dirname(opt$p6a_dir), "Pipeline6e_v1.0_PermutationTest")
}

################################################################################
# 参数配置
################################################################################

PARAMS <- list(
  p6a_dir    = opt$p6a_dir,
  p6b_dir    = opt$p6b_dir,
  output_dir = opt$output_dir,
  genotypes  = if (opt$genotype == "both") c("WT", "Mutant") else opt$genotype,
  
  # 置换
  n_perm       = opt$n_perm,
  n_replicates = opt$n_replicates,
  n_cv_folds   = opt$n_folds,
  rf_ntree     = opt$ntree,
  
  # 负样本
  use_gc_match = opt$gc_match,
  gc_bins      = opt$gc_bins,
  neg_ratio    = opt$neg_ratio,
  
  # 资源
  n_cores = opt$cores,
  seed    = opt$seed,
  
  # 可视化
  fig_width  = opt$fig_width,
  fig_height = opt$fig_height,
  alpha      = opt$alpha
)

# 创建输出目录
for (geno in PARAMS$genotypes) {
  dir.create(file.path(PARAMS$output_dir, geno), recursive = TRUE, showWarnings = FALSE)
}

################################################################################
# 日志设置
################################################################################

log_file <- file.path(PARAMS$output_dir,
                      paste0("pipeline6e_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("================================================================================\n")
cat("  Pipeline 6e v1.0: Permutation Baseline Test\n")
cat("================================================================================\n\n")
cat(sprintf("开始时间: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("6a 目录: %s\n", PARAMS$p6a_dir))
cat(sprintf("6b 目录: %s\n", PARAMS$p6b_dir))
cat(sprintf("输出目录: %s\n\n", PARAMS$output_dir))

cat("核心参数:\n")
cat(sprintf("  基因型: %s\n", paste(PARAMS$genotypes, collapse = ", ")))
cat(sprintf("  置换次数: %d\n", PARAMS$n_perm))
cat(sprintf("  简化 CV: %d replicates × %d-fold\n", PARAMS$n_replicates, PARAMS$n_cv_folds))
cat(sprintf("  RF ntree: %d\n", PARAMS$rf_ntree))
cat(sprintf("  GC 匹配: %s\n", ifelse(PARAMS$use_gc_match,
                                       sprintf("启用 (%d bins, ratio=%.1f)", PARAMS$gc_bins, PARAMS$neg_ratio),
                                       sprintf("禁用 (ratio=%.1f)", PARAMS$neg_ratio))))
cat(sprintf("  显著性阈值: %.2f\n", PARAMS$alpha))
cat(sprintf("  并行核心: %d\n", PARAMS$n_cores))
cat("\n")

set.seed(PARAMS$seed)

################################################################################
# 辅助函数: GC 匹配采样 (复用 6b 逻辑)
################################################################################

gc_matched_sampling <- function(pos_genes, neg_pool, gc_info, n_bins = 10, ratio = 1) {
  pos_gc <- gc_info[pos_genes]
  pos_gc <- pos_gc[!is.na(pos_gc)]
  
  if (length(pos_gc) == 0) {
    return(sample(neg_pool, min(length(neg_pool), length(pos_genes) * ratio)))
  }
  
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
# 核心函数: 单次置换的简化 RF CV
################################################################################

#' 在标签置换后运行简化版 RF CV, 只返回 AUC (不追踪 importance)
#' @param X_combined 合并后的特征矩阵 (pos + neg)
#' @param y_permuted 置换后的标签向量 (factor)
#' @param n_reps 简化 CV 重复数
#' @param n_folds CV 折数
#' @param ntree RF 树数
#' @param seed_offset 种子偏移
#' @return 单个 mean_auc 值
run_permuted_rf <- function(X_combined, y_permuted, n_reps, n_folds, 
                            ntree, seed_offset) {
  
  n_total <- nrow(X_combined)
  n_features <- ncol(X_combined)
  best_mtry <- max(1, floor(sqrt(n_features)))
  
  rep_aucs <- numeric(n_reps)
  
  for (rep_id in 1:n_reps) {
    set.seed(seed_offset + rep_id)
    
    fold_assignments <- sample(rep(1:n_folds, length.out = n_total))
    fold_aucs <- numeric(n_folds)
    
    for (fold in 1:n_folds) {
      test_idx  <- which(fold_assignments == fold)
      train_idx <- which(fold_assignments != fold)
      
      # 跳过单类 fold
      y_train <- y_permuted[train_idx]
      y_test  <- as.numeric(as.character(y_permuted[test_idx]))
      
      if (length(unique(y_train)) < 2 || length(unique(y_test)) < 2) {
        fold_aucs[fold] <- 0.5
        next
      }
      
      tryCatch({
        rf_model <- randomForest(
          x = X_combined[train_idx, , drop = FALSE],
          y = y_train,
          ntree = ntree,
          mtry = best_mtry,
          importance = FALSE  # 置换测试不需要 importance
        )
        
        y_prob <- predict(rf_model, X_combined[test_idx, , drop = FALSE], type = "prob")[, "1"]
        
        fold_aucs[fold] <- tryCatch({
          as.numeric(pROC::auc(pROC::roc(y_test, y_prob, quiet = TRUE,
                                          levels = c(0, 1), direction = "<")))
        }, error = function(e) 0.5)
        
      }, error = function(e) {
        fold_aucs[fold] <- 0.5
      })
    }
    
    rep_aucs[rep_id] <- mean(fold_aucs, na.rm = TRUE)
  }
  
  mean(rep_aucs, na.rm = TRUE)
}

################################################################################
# 模块: 单个 Cluster 的置换检验
################################################################################

permutation_test_cluster <- function(cluster_name, cluster_data, kmer_matrix,
                                     ml_nonresponsive, gc_content, real_auc,
                                     params) {
  
  cat(sprintf("\n【%s】 (real AUC = %.3f)\n", cluster_name, real_auc))
  
  pos_genes <- intersect(cluster_data$pos_genes, rownames(kmer_matrix))
  neg_pool  <- intersect(ml_nonresponsive, rownames(kmer_matrix))
  
  n_pos <- length(pos_genes)
  
  if (n_pos < 30) {
    cat("  ⚠ 正样本不足 30，跳过\n")
    return(NULL)
  }
  
  feature_names <- colnames(kmer_matrix)
  n_features <- length(feature_names)
  
  if (n_features < 10) {
    cat("  ⚠ 特征不足 10，跳过\n")
    return(NULL)
  }
  
  X_pos <- as.matrix(kmer_matrix[pos_genes, feature_names, drop = FALSE])
  
  cat(sprintf("  正样本: %d | 负样本池: %d | 特征: %d\n",
              n_pos, length(neg_pool), n_features))
  cat(sprintf("  运行 %d 次置换 (%d reps × %d-fold)...\n",
              params$n_perm, params$n_replicates, params$n_cv_folds))
  
  # 并行执行置换
  registerDoParallel(cores = min(params$n_cores, params$n_perm))
  
  null_aucs <- foreach(
    perm_id = 1:params$n_perm,
    .combine = c,
    .packages = c("randomForest", "pROC")
  ) %dopar% {
    
    # 1. 独立采样负样本 (与 6b v4.1 一致)
    set.seed(params$seed + perm_id * 13)
    
    if (params$use_gc_match && !is.null(gc_content)) {
      neg_genes <- gc_matched_sampling(pos_genes, neg_pool, gc_content,
                                        n_bins = params$gc_bins,
                                        ratio = params$neg_ratio)
    } else {
      max_neg <- min(length(neg_pool), n_pos * params$neg_ratio)
      neg_genes <- sample(neg_pool, max_neg)
    }
    
    X_neg <- as.matrix(kmer_matrix[neg_genes, feature_names, drop = FALSE])
    
    # 2. 合并数据集
    X_combined <- rbind(X_pos, X_neg)
    y_true <- factor(c(rep(1, n_pos), rep(0, nrow(X_neg))))
    
    # 3. 置换标签
    y_permuted <- factor(sample(y_true))
    
    # 4. 简化 RF CV
    run_permuted_rf(
      X_combined  = X_combined,
      y_permuted  = y_permuted,
      n_reps      = params$n_replicates,
      n_folds     = params$n_cv_folds,
      ntree       = params$rf_ntree,
      seed_offset = params$seed + perm_id * 1000
    )
  }
  
  stopImplicitCluster()
  
  # 经验 p-value (加 1 修正)
  p_value <- (sum(null_aucs >= real_auc) + 1) / (params$n_perm + 1)
  
  # Null 分布统计
  null_mean <- mean(null_aucs)
  null_sd   <- sd(null_aucs)
  null_q95  <- quantile(null_aucs, 0.95)
  null_q99  <- quantile(null_aucs, 0.99)
  
  # Effect size (标准化距离)
  effect_size <- (real_auc - null_mean) / null_sd
  
  sig_label <- ifelse(p_value < 0.001, "***",
               ifelse(p_value < 0.01, "**",
               ifelse(p_value < params$alpha, "*", "ns")))
  
  cat(sprintf("  ✔ Null AUC = %.3f ± %.3f | Real AUC = %.3f | p = %.4f %s | d = %.2f\n",
              null_mean, null_sd, real_auc, p_value, sig_label, effect_size))
  
  list(
    cluster      = cluster_name,
    n_pos        = n_pos,
    n_features   = n_features,
    real_auc     = real_auc,
    null_mean    = null_mean,
    null_sd      = null_sd,
    null_q95     = as.numeric(null_q95),
    null_q99     = as.numeric(null_q99),
    p_value      = p_value,
    effect_size  = effect_size,
    sig          = sig_label,
    null_aucs    = null_aucs,
    n_perm       = params$n_perm
  )
}

################################################################################
# 模块: 可视化 — Null Distribution 面板图
################################################################################

plot_permutation_results <- function(results_list, genotype, output_dir, params) {
  
  n_clusters <- length(results_list)
  if (n_clusters == 0) return(invisible(NULL))
  
  # 准备绘图数据
  plot_data <- do.call(rbind, lapply(results_list, function(res) {
    data.frame(
      cluster   = res$cluster,
      null_auc  = res$null_aucs,
      real_auc  = res$real_auc,
      p_value   = res$p_value,
      sig       = res$sig,
      stringsAsFactors = FALSE
    )
  }))
  
  # 真实 AUC 标注数据
  real_df <- data.frame(
    cluster  = sapply(results_list, function(r) r$cluster),
    real_auc = sapply(results_list, function(r) r$real_auc),
    p_value  = sapply(results_list, function(r) r$p_value),
    sig      = sapply(results_list, function(r) r$sig),
    label    = sapply(results_list, function(r) {
      sprintf("AUC=%.3f\np=%.4f %s", r$real_auc, r$p_value, r$sig)
    }),
    stringsAsFactors = FALSE
  )
  
  # 按真实 AUC 降序排列 cluster
  cluster_order <- real_df$cluster[order(real_df$real_auc, decreasing = TRUE)]
  plot_data$cluster <- factor(plot_data$cluster, levels = cluster_order)
  real_df$cluster   <- factor(real_df$cluster, levels = cluster_order)
  
  # 自适应面板布局
  n_cols <- min(4, n_clusters)
  n_rows <- ceiling(n_clusters / n_cols)
  
  fig_w <- params$fig_width
  fig_h <- max(params$fig_height, n_rows * 2.5)
  
  p <- ggplot(plot_data, aes(x = null_auc)) +
    geom_histogram(aes(y = after_stat(density)),
                   bins = 30, fill = "grey70", color = "grey40", linewidth = 0.3) +
    geom_density(color = "steelblue", linewidth = 0.8) +
    geom_vline(data = real_df, aes(xintercept = real_auc),
               color = "firebrick", linewidth = 1, linetype = "dashed") +
    geom_vline(xintercept = 0.5, color = "grey50", linewidth = 0.5, linetype = "dotted") +
    geom_label(data = real_df,
               aes(x = real_auc, y = Inf, label = label),
               vjust = 1.2, hjust = -0.05, size = 2.8,
               fill = "white", alpha = 0.8, label.size = 0.3) +
    facet_wrap(~ cluster, scales = "free_y", ncol = n_cols) +
    scale_x_continuous(limits = c(0.3, 1.0)) +
    labs(
      title = sprintf("Permutation Baseline Test — %s (%d permutations)",
                       genotype, params$n_perm),
      subtitle = sprintf("Simplified CV: %d reps × %d-fold, ntree=%d",
                          params$n_replicates, params$n_cv_folds, params$rf_ntree),
      x = "AUC-ROC (Null Distribution)",
      y = "Density"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(face = "bold", size = 9),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(color = "grey40", size = 10),
      panel.grid.minor = element_blank()
    )
  
  pdf_path <- file.path(output_dir, genotype,
                        sprintf("Fig_Permutation_Test_%s.pdf", genotype))
  ggsave(pdf_path, plot = p, width = fig_w, height = fig_h, limitsize = FALSE)
  cat(sprintf("  📊 Figure saved: %s\n", pdf_path))
  
  # 同时输出 PNG
  png_path <- sub("\\.pdf$", ".png", pdf_path)
  ggsave(png_path, plot = p, width = fig_w, height = fig_h, dpi = 200, limitsize = FALSE)
  
  invisible(p)
}

################################################################################
# 模块: 处理单个基因型
################################################################################

process_genotype <- function(genotype, params) {
  cat("\n================================================================================\n")
  cat(sprintf("  处理基因型: %s\n", genotype))
  cat("================================================================================\n")
  
  # ---- 1. 加载 6a ML_ready 数据 ----
  ml_ready_file <- file.path(params$p6a_dir, genotype, "ML_ready", "Pipeline6b_Input_Data.rds")
  
  if (!file.exists(ml_ready_file)) {
    cat(sprintf("  ⚠ ML_ready 数据不存在: %s\n", ml_ready_file))
    return(NULL)
  }
  
  data_6a <- readRDS(ml_ready_file)
  cat(sprintf("  6a 数据: %d clusters, %d genes × %d features\n",
              length(data_6a$cluster_data),
              nrow(data_6a$kmer_matrix),
              ncol(data_6a$kmer_matrix)))
  
  # ---- 2. 加载 6b 真实 AUC ----
  rf_results_file <- file.path(params$p6b_dir, genotype, "results", "Pipeline6b_RF_Results.rds")
  
  if (!file.exists(rf_results_file)) {
    cat(sprintf("  ⚠ 6b 结果不存在: %s\n", rf_results_file))
    return(NULL)
  }
  
  rf_results <- readRDS(rf_results_file)
  real_aucs <- sapply(rf_results, function(x) x$mean_auc)
  cat(sprintf("  6b 结果: %d clusters, AUC range [%.3f, %.3f]\n",
              length(real_aucs), min(real_aucs), max(real_aucs)))
  
  # ---- 3. 对每个 Cluster 运行置换检验 ----
  perm_results <- list()
  
  for (cluster_name in names(data_6a$cluster_data)) {
    
    # 必须在 6b 中有对应的真实 AUC
    if (!(cluster_name %in% names(rf_results))) {
      cat(sprintf("\n【%s】⚠ 6b 无对应结果，跳过\n", cluster_name))
      next
    }
    
    real_auc <- rf_results[[cluster_name]]$mean_auc
    
    result <- permutation_test_cluster(
      cluster_name     = cluster_name,
      cluster_data     = data_6a$cluster_data[[cluster_name]],
      kmer_matrix      = data_6a$kmer_matrix,
      ml_nonresponsive = data_6a$ml_nonresponsive_genes,
      gc_content       = data_6a$gc_content,
      real_auc         = real_auc,
      params           = params
    )
    
    if (!is.null(result)) {
      perm_results[[cluster_name]] <- result
    }
    
    gc(verbose = FALSE)
  }
  
  if (length(perm_results) == 0) {
    cat("\n  ⚠ 没有成功的置换结果\n")
    return(NULL)
  }
  
  # ---- 4. 输出汇总表 ----
  geno_output_dir <- file.path(params$output_dir, genotype)
  
  summary_df <- data.frame(
    Cluster     = sapply(perm_results, function(r) r$cluster),
    N_Pos       = sapply(perm_results, function(r) r$n_pos),
    N_Features  = sapply(perm_results, function(r) r$n_features),
    Real_AUC    = sapply(perm_results, function(r) round(r$real_auc, 4)),
    Null_Mean   = sapply(perm_results, function(r) round(r$null_mean, 4)),
    Null_SD     = sapply(perm_results, function(r) round(r$null_sd, 4)),
    Null_Q95    = sapply(perm_results, function(r) round(r$null_q95, 4)),
    Null_Q99    = sapply(perm_results, function(r) round(r$null_q99, 4)),
    P_value     = sapply(perm_results, function(r) round(r$p_value, 5)),
    Effect_Size = sapply(perm_results, function(r) round(r$effect_size, 2)),
    Sig         = sapply(perm_results, function(r) r$sig),
    N_Perm      = sapply(perm_results, function(r) r$n_perm),
    stringsAsFactors = FALSE,
    row.names = NULL
  ) %>%
    arrange(P_value)
  
  write.table(summary_df,
              file.path(geno_output_dir, "Permutation_Summary.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # 完整 null 分布 (供后续分析)
  saveRDS(perm_results,
          file.path(geno_output_dir, "Permutation_NullDistributions.rds"))
  
  # ---- 5. 可视化 ----
  plot_permutation_results(perm_results, genotype, params$output_dir, params)
  
  # ---- 6. 打印汇总 ----
  cat(sprintf("\n--- %s 置换检验汇总 ---\n", genotype))
  n_sig <- sum(summary_df$P_value < params$alpha)
  cat(sprintf("  显著 Cluster (p < %.2f): %d / %d\n",
              params$alpha, n_sig, nrow(summary_df)))
  print(summary_df, row.names = FALSE)
  
  return(list(genotype = genotype, summary = summary_df, results = perm_results))
}

################################################################################
# 主程序
################################################################################

main <- function() {
  
  all_results <- list()
  
  for (genotype in PARAMS$genotypes) {
    result <- process_genotype(genotype, PARAMS)
    if (!is.null(result)) {
      all_results[[genotype]] <- result
    }
  }
  
  # 全局汇总
  cat("\n================================================================================\n")
  cat("  Pipeline 6e v1.0 完成: Permutation Baseline Test\n")
  cat("================================================================================\n\n")
  cat(sprintf("结束时间: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  
  cat("【全局汇总】\n")
  cat(sprintf("  置换次数: %d | CV: %d reps × %d-fold | ntree: %d\n",
              PARAMS$n_perm, PARAMS$n_replicates, PARAMS$n_cv_folds, PARAMS$rf_ntree))
  
  for (geno in names(all_results)) {
    s <- all_results[[geno]]$summary
    n_sig <- sum(s$P_value < PARAMS$alpha)
    cat(sprintf("\n  %s: %d/%d clusters 通过置换检验 (p < %.2f)\n",
                geno, n_sig, nrow(s), PARAMS$alpha))
    
    # 标注未通过的 cluster
    failed <- s %>% filter(P_value >= PARAMS$alpha)
    if (nrow(failed) > 0) {
      cat("    ⚠ 未通过:\n")
      for (i in 1:nrow(failed)) {
        cat(sprintf("      - %s: AUC=%.3f, null=%.3f±%.3f, p=%.4f\n",
                    failed$Cluster[i], failed$Real_AUC[i],
                    failed$Null_Mean[i], failed$Null_SD[i], failed$P_value[i]))
      }
    }
  }
  
  # 合并全基因型汇总表
  if (length(all_results) > 1) {
    combined_summary <- do.call(rbind, lapply(names(all_results), function(geno) {
      all_results[[geno]]$summary %>% mutate(Genotype = geno, .before = 1)
    }))
    write.table(combined_summary,
                file.path(PARAMS$output_dir, "Permutation_Summary_All.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  cat(sprintf("\n输出目录: %s\n", PARAMS$output_dir))
  cat("\n分析完成。\n")
}

main()

sink(type = "message")
sink()
close(log_con)
