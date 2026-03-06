#!/usr/bin/env Rscript
################################################################################
# RNAseq Pipeline 6b v4.1: Random Forest Model (Unified Version)
#
# 功能特性：
#   - 配套Pipeline 6a v8.0的输出结构
#   - 支持处理单个或多个基因型
#   - 完全移除FDR过滤逻辑（使用所有Top N特征）
#   - 可选GC匹配负样本
#   - v4.1: 修复方差塌陷 — 每个 Replicate 独立采样负样本
#
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(randomForest)
  library(caret)
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
  make_option(c("-i", "--input_dir"), type = "character", default = NULL,
              help = "Pipeline 6a 输出目录（包含WT/Mutant子目录）"),
  make_option(c("-o", "--output_dir"), type = "character", default = NULL,
              help = "输出目录"),
  make_option(c("--genotype"), type = "character", default = "both",
              help = "处理的基因型: 'WT', 'Mutant', 或 'both' [default: %default]"),
  
  # === 负样本策略 ===
  make_option(c("--gc_match"), action = "store_true", default = FALSE,
              help = "是否对负样本进行GC匹配采样"),
  make_option(c("--gc_bins"), type = "integer", default = 10,
              help = "GC匹配的分层数量 [default: %default]"),
  make_option(c("--neg_ratio"), type = "double", default = 3,
              help = "负样本与正样本的比例 [default: %default]"),
  
  # === RF参数 ===
  make_option(c("--ntree"), type = "integer", default = 500,
              help = "随机森林树数量 [default: %default]"),
  make_option(c("--replicates"), type = "integer", default = 100,
              help = "CV重复次数 [default: %default]"),
  make_option(c("--folds"), type = "integer", default = 10,
              help = "CV折数 [default: %default]"),
  make_option(c("--cores"), type = "integer", default = 40,
              help = "并行核心数 [default: %default]")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# 设置默认路径
if (is.null(opt$input_dir)) {
  opt$input_dir <- "/home/cfy/Data/DOF/06_pCRE/Pipeline6a_v8.0_GenomeWideBG_Partitioned"
}
if (is.null(opt$output_dir)) {
  gc_suffix <- ifelse(opt$gc_match, "_GCmatched", "")
  opt$output_dir <- paste0(dirname(opt$input_dir), "/Pipeline6b_v4.0", gc_suffix)
}

################################################################################
# 参数配置
################################################################################

PARAMS <- list(
  input_dir = opt$input_dir,
  output_dir = opt$output_dir,
  genotypes_to_process = if (opt$genotype == "both") c("WT", "Mutant") else opt$genotype,
  
  # GC匹配
  use_gc_match = opt$gc_match,
  gc_bins = opt$gc_bins,
  neg_ratio = opt$neg_ratio,
  
  # RF参数
  n_replicates = opt$replicates,
  n_cv_folds = opt$folds,
  rf_ntree = opt$ntree,
  n_cores = opt$cores,
  
  seed = 42
)

################################################################################
# 创建输出目录
################################################################################

for (geno in PARAMS$genotypes_to_process) {
  for (subdir in c("logs", "results", "importance")) {
    dir.create(file.path(PARAMS$output_dir, geno, subdir), 
               recursive = TRUE, showWarnings = FALSE)
  }
}

################################################################################
# 日志设置
################################################################################

log_file <- file.path(PARAMS$output_dir, 
                      paste0("pipeline6b_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
dir.create(dirname(log_file), recursive = TRUE, showWarnings = FALSE)
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("================================================================================\n")
cat("  Pipeline 6b v4.1: Random Forest Model (Unified Version)\n")
cat("================================================================================\n\n")
cat(sprintf("开始时间: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("输入目录: %s\n", PARAMS$input_dir))
cat(sprintf("输出目录: %s\n\n", PARAMS$output_dir))

cat("核心参数:\n")
cat(sprintf("  基因型: %s\n", paste(PARAMS$genotypes_to_process, collapse = ", ")))
cat(sprintf("  GC匹配: %s\n", ifelse(PARAMS$use_gc_match, 
                                    sprintf("启用 (%d bins, ratio=%.1f)", 
                                            PARAMS$gc_bins, PARAMS$neg_ratio), 
                                    "禁用")))
cat(sprintf("  CV策略: %d replicates × %d-fold balanced CV\n", 
            PARAMS$n_replicates, PARAMS$n_cv_folds))
cat(sprintf("  随机森林: ntree=%d\n", PARAMS$rf_ntree))
cat(sprintf("  特征选择: 使用6a输出的所有白名单特征 (无FDR过滤)\n"))
cat(sprintf("  负样本策略: 每个Replicate独立采样 (v4.1, 消除方差塌陷)\n"))
cat("\n")

set.seed(PARAMS$seed)

################################################################################
# 辅助函数
################################################################################

#' GC匹配采样
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
# 模块 1: 加载数据
################################################################################

load_genotype_data <- function(input_dir, genotype) {
  cat(sprintf("\n--- 加载 %s 数据 ---\n", genotype))
  
  ml_ready_dir <- file.path(input_dir, genotype, "ML_ready")
  input_file <- file.path(ml_ready_dir, "Pipeline6b_Input_Data.rds")
  
  if (!file.exists(input_file)) {
    cat(sprintf("  ⚠ 文件不存在: %s\n", input_file))
    return(NULL)
  }
  
  data <- readRDS(input_file)
  
  cat(sprintf("  Cluster数量: %d\n", length(data$cluster_data)))
  cat(sprintf("  k-mer矩阵: %d genes x %d features\n", 
              nrow(data$kmer_matrix), ncol(data$kmer_matrix)))
  cat(sprintf("  ML非响应基因: %d\n", length(data$ml_nonresponsive_genes)))
  
  if (!is.null(data$gc_content)) {
    cat(sprintf("  GC含量信息: 可用 (均值 %.1f%%)\n", 
                mean(data$gc_content, na.rm = TRUE) * 100))
  }
  
  return(data)
}

################################################################################
# 模块 2: 单次Replicate的Balanced CV
################################################################################

run_single_replicate <- function(rep_id, X_pos, X_neg, feature_names, 
                                  n_folds, ntree, mtry, seed_base) {
  set.seed(seed_base + rep_id)
  
  n_pos <- nrow(X_pos)
  n_neg <- nrow(X_neg)
  
  pos_folds <- sample(rep(1:n_folds, length.out = n_pos))
  
  fold_aucs <- numeric(n_folds)
  fold_importance <- matrix(0, nrow = n_folds, ncol = length(feature_names))
  colnames(fold_importance) <- feature_names
  
  for (fold in 1:n_folds) {
    pos_test_idx <- which(pos_folds == fold)
    pos_train_idx <- which(pos_folds != fold)
    
    n_pos_train <- length(pos_train_idx)
    n_pos_test <- length(pos_test_idx)
    
    # Balanced: 每次重新采样负样本
    neg_train_idx <- sample(1:n_neg, min(n_pos_train, n_neg), replace = FALSE)
    remaining_neg <- setdiff(1:n_neg, neg_train_idx)
    neg_test_idx <- sample(remaining_neg, min(n_pos_test, length(remaining_neg)), replace = FALSE)
    
    X_train <- rbind(X_pos[pos_train_idx, , drop = FALSE],
                     X_neg[neg_train_idx, , drop = FALSE])
    y_train <- factor(c(rep(1, n_pos_train), rep(0, length(neg_train_idx))))
    
    X_test <- rbind(X_pos[pos_test_idx, , drop = FALSE],
                    X_neg[neg_test_idx, , drop = FALSE])
    y_test <- c(rep(1, n_pos_test), rep(0, length(neg_test_idx)))
    
    tryCatch({
      rf_model <- randomForest(
        x = X_train,
        y = y_train,
        ntree = ntree,
        mtry = mtry,
        importance = TRUE
      )
      
      y_prob <- predict(rf_model, X_test, type = "prob")[, "1"]
      
      fold_aucs[fold] <- tryCatch({
        as.numeric(pROC::auc(pROC::roc(y_test, y_prob, quiet = TRUE, 
                                        levels = c(0, 1), direction = "<")))
      }, error = function(e) 0.5)
      
      fold_importance[fold, ] <- importance(rf_model)[, "MeanDecreaseGini"]
      
    }, error = function(e) {
      fold_aucs[fold] <- NA
      fold_importance[fold, ] <- 0
    })
  }
  
  list(
    mean_auc = mean(fold_aucs, na.rm = TRUE),
    fold_aucs = fold_aucs,
    mean_importance = colMeans(fold_importance, na.rm = TRUE)
  )
}

################################################################################
# 模块 3: Cluster级别RF分析
################################################################################

run_rf_for_cluster <- function(cluster_name, cluster_data, kmer_matrix, 
                                ml_nonresponsive, gc_content, params) {
  
  cat(sprintf("\n【%s】\n", cluster_name))
  
  pos_genes <- intersect(cluster_data$pos_genes, rownames(kmer_matrix))
  neg_pool <- intersect(ml_nonresponsive, rownames(kmer_matrix))
  
  n_pos <- length(pos_genes)
  cat(sprintf("  正样本: %d\n", n_pos))
  cat(sprintf("  负样本池: %d\n", length(neg_pool)))
  
  if (n_pos < 30) {
    cat("  ⚠ 正样本不足30，跳过\n")
    return(NULL)
  }
  
  # 使用所有可用特征（不做FDR过滤）
  feature_names <- colnames(kmer_matrix)
  n_features <- length(feature_names)
  cat(sprintf("  特征数: %d (使用全部白名单特征)\n", n_features))
  
  if (n_features < 10) {
    cat("  ⚠ 特征数不足10，跳过\n")
    return(NULL)
  }
  
  # ---- v4.1 修复: 每个 Replicate 独立采样负样本 (消除方差塌陷) ----
  # 预生成 N 组独立的负样本基因集，使 sd_auc 能正确反映负样本选择的方差
  cat(sprintf("  预生成 %d 组独立负样本...\n", params$n_replicates))
  
  neg_samples_list <- lapply(1:params$n_replicates, function(rep_id) {
    set.seed(params$seed + rep_id * 7)  # 每组不同的种子
    if (params$use_gc_match && !is.null(gc_content)) {
      gc_matched_sampling(pos_genes, neg_pool, gc_content,
                          n_bins = params$gc_bins,
                          ratio = params$neg_ratio)
    } else {
      max_neg <- min(length(neg_pool), n_pos * params$neg_ratio)
      sample(neg_pool, max_neg)
    }
  })
  
  # 报告负样本统计 (基于全部组的汇总)
  neg_sizes <- sapply(neg_samples_list, length)
  cat(sprintf("  负样本/组: median=%d, range=[%d, %d]\n",
              median(neg_sizes), min(neg_sizes), max(neg_sizes)))
  
  if (params$use_gc_match && !is.null(gc_content)) {
    pos_gc_mean <- mean(gc_content[pos_genes], na.rm = TRUE)
    # 汇报所有组的 neg GC 均值分布
    neg_gc_means <- sapply(neg_samples_list, function(ng) mean(gc_content[ng], na.rm = TRUE))
    cat(sprintf("  GC含量: 正=%.1f%%, 负=%.1f%%±%.1f%%\n",
                pos_gc_mean * 100, mean(neg_gc_means) * 100, sd(neg_gc_means) * 100))
  }
  
  n_neg <- median(neg_sizes)  # 代表值用于汇总报告
  
  # 提取正样本特征矩阵 (所有 replicate 共享)
  X_pos <- as.matrix(kmer_matrix[pos_genes, feature_names, drop = FALSE])
  
  # 确定mtry
  best_mtry <- max(1, floor(sqrt(n_features)))
  
  # 并行运行 (每个 replicate 使用独立的负样本集)
  cat(sprintf("  运行 %d replicates × %d-fold CV (mtry=%d)...\n", 
              params$n_replicates, params$n_cv_folds, best_mtry))
  
  registerDoParallel(cores = min(params$n_cores, params$n_replicates))
  
  all_results <- foreach(
    rep_id = 1:params$n_replicates,
    .combine = list,
    .multicombine = TRUE,
    .maxcombine = params$n_replicates
  ) %dopar% {
    # 每个 replicate 从预生成列表中获取自己的负样本
    neg_genes_rep <- neg_samples_list[[rep_id]]
    X_neg_rep <- as.matrix(kmer_matrix[neg_genes_rep, feature_names, drop = FALSE])
    
    run_single_replicate(
      rep_id = rep_id,
      X_pos = X_pos,
      X_neg = X_neg_rep,
      feature_names = feature_names,
      n_folds = params$n_cv_folds,
      ntree = params$rf_ntree,
      mtry = best_mtry,
      seed_base = params$seed
    )
  }
  
  stopImplicitCluster()
  
  # 聚合结果
  all_aucs <- sapply(all_results, function(x) x$mean_auc)
  valid_aucs <- all_aucs[!is.na(all_aucs)]
  
  if (length(valid_aucs) < 10) {
    cat("  ⚠ 有效结果不足\n")
    return(NULL)
  }
  
  mean_auc <- mean(valid_aucs)
  sd_auc <- sd(valid_aucs)
  
  # 聚合特征重要性
  importance_matrix <- do.call(rbind, lapply(all_results, function(x) x$mean_importance))
  mean_importance <- colMeans(importance_matrix, na.rm = TRUE)
  sd_importance <- apply(importance_matrix, 2, sd, na.rm = TRUE)
  selection_freq <- colMeans(importance_matrix > 0, na.rm = TRUE)
  
  importance_df <- data.frame(
    feature = feature_names,
    mean_importance = mean_importance,
    sd_importance = sd_importance,
    selection_freq = selection_freq,
    stringsAsFactors = FALSE
  ) %>%
    arrange(desc(mean_importance)) %>%
    mutate(
      rank = row_number(),
      Region = sub("_[ACGT]+$", "", feature),
      kmer = sub("^[^_]+_", "", feature)
    )
  
  cat(sprintf("  ✔ AUC-ROC = %.3f ± %.3f\n", mean_auc, sd_auc))
  
  # 显示Top特征
  cat("  Top 10 重要pCRE:\n")
  top10 <- head(importance_df, 10)
  for (i in 1:nrow(top10)) {
    cat(sprintf("    %2d. %s (imp=%.3f)\n", 
                i, top10$feature[i], top10$mean_importance[i]))
  }
  
  list(
    cluster = cluster_name,
    n_pos = n_pos,
    n_neg = n_neg,
    n_features = n_features,
    best_mtry = best_mtry,
    mean_auc = mean_auc,
    sd_auc = sd_auc,
    all_aucs = valid_aucs,
    importance = importance_df,
    gc_match_used = params$use_gc_match
  )
}

################################################################################
# 模块 4: 处理单个基因型
################################################################################

process_genotype <- function(genotype, data, output_dir, params) {
  cat("\n================================================================================\n")
  cat(sprintf("  处理基因型: %s\n", genotype))
  cat("================================================================================\n")
  
  if (is.null(data)) {
    cat("  ⚠ 无数据，跳过\n")
    return(NULL)
  }
  
  all_results <- list()
  
  for (cluster_name in names(data$cluster_data)) {
    result <- run_rf_for_cluster(
      cluster_name = cluster_name,
      cluster_data = data$cluster_data[[cluster_name]],
      kmer_matrix = data$kmer_matrix,
      ml_nonresponsive = data$ml_nonresponsive_genes,
      gc_content = data$gc_content,
      params = params
    )
    
    if (!is.null(result)) {
      all_results[[cluster_name]] <- result
    }
    
    gc(verbose = FALSE)
  }
  
  if (length(all_results) == 0) {
    cat("\n⚠ 没有成功的分析结果\n")
    return(NULL)
  }
  
  # 保存结果
  geno_output_dir <- file.path(output_dir, genotype)
  
  # 1. 性能汇总
  perf_summary <- data.frame(
    Cluster = names(all_results),
    N_Pos = sapply(all_results, function(x) x$n_pos),
    N_Neg = sapply(all_results, function(x) x$n_neg),
    N_Features = sapply(all_results, function(x) x$n_features),
    Best_mtry = sapply(all_results, function(x) x$best_mtry),
    Mean_AUC = sapply(all_results, function(x) round(x$mean_auc, 3)),
    SD_AUC = sapply(all_results, function(x) round(x$sd_auc, 3)),
    GC_Match = sapply(all_results, function(x) x$gc_match_used),
    stringsAsFactors = FALSE
  )
  
  write.table(perf_summary,
              file.path(geno_output_dir, "results", "RF_Performance_Summary.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # 2. 特征重要性
  for (cl_name in names(all_results)) {
    write.table(all_results[[cl_name]]$importance,
                file.path(geno_output_dir, "importance", 
                          sprintf("%s_importance.txt", cl_name)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
  
  # 3. 完整结果
  saveRDS(all_results, file.path(geno_output_dir, "results", "Pipeline6b_RF_Results.rds"))
  
  cat("\n  ✔ 结果已保存\n")
  
  return(list(
    genotype = genotype,
    results = all_results,
    summary = perf_summary
  ))
}

################################################################################
# 主程序
################################################################################

main <- function() {
  
  all_genotype_results <- list()
  
  for (genotype in PARAMS$genotypes_to_process) {
    # 加载数据
    data <- load_genotype_data(PARAMS$input_dir, genotype)
    
    # 运行RF分析
    result <- process_genotype(genotype, data, PARAMS$output_dir, PARAMS)
    
    if (!is.null(result)) {
      all_genotype_results[[genotype]] <- result
    }
  }
  
  # 汇总报告
  cat("\n================================================================================\n")
  cat("  Pipeline 6b v4.1 完成\n")
  cat("================================================================================\n\n")
  cat(sprintf("结束时间: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  
  cat("【性能汇总】\n")
  for (genotype in names(all_genotype_results)) {
    cat(sprintf("\n%s:\n", genotype))
    print(all_genotype_results[[genotype]]$summary %>% arrange(desc(Mean_AUC)))
  }
  
  cat(sprintf("\n输出目录: %s\n", PARAMS$output_dir))
  cat("\n分析完成。\n")
}

main()

sink(type = "message")
sink()
close(log_con)