#!/usr/bin/env Rscript
################################################################################
# RNAseq Pipeline 7b: Cluster-Constrained Discovery + MECS Scoring
#
# 版本: v1.0
#
# ============================================================================
# 模块概览
# ============================================================================
# Module A: Cluster-Constrained TF Hub Discovery
#           (超几何检验 per TF per Cluster, TF Family 级 Fisher 检验)
# Module B: agriGO Gene List Export
#           (按 Cluster / EdgeClass / Focus TF 组织的基因列表)
# Module C: TF Family × Cluster × Edge Class Enrichment Matrix
#           (Fisher 精确检验, 长/宽格式输出)
# Module D: Focus TF Rewiring Profiles
#           (聚焦分析指定 TF 家族或特定基因)
# Module E: MECS (Multi-Evidence Confidence Score) + Tier Classification
#           (S1: 差异共表达 + S2: 聚类约束 + S3: pCRE一致性[可选])
#
# ============================================================================
# 工程化泛化设计
# ============================================================================
# - --focus_tf_family / --focus_gene 参数化, 不硬编码任何基因/家族
# - Cluster 数量自动适配 (从数据推断, 不假设固定数量)
# - 每层证据 (P5 聚类, P5b fate, P6 pCRE) 均可选, 缺失时优雅退化
# - 与 7a 一致的 optparse + JSON + 日志 + RData 风格
#
# ============================================================================
# 使用方式
# ============================================================================
#   Rscript RNAseq_pipeline7b_discovery.R \
#     --rdata_7a /path/to/Pipeline7a_DiffCoexp_Results.RData \
#     --out_dir /path/to/07_GRN \
#     --focus_tf_family DOF \
#     --focus_gene Cre06.g250050
#
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(jsonlite)
})

################################################################################
# 参数解析
################################################################################

parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  params <- list(
    # ========== 数据输入 ==========
    rdata_7a          = NULL,        # [必需] Pipeline7a RData
    cluster_pairs     = NULL,        # [可选] Cluster pair TSV
    pcre_tf_map       = NULL,        # [可选] pCRE-TF mapping TSV
    highlight_file    = NULL,        # [可选] Highlight genes (覆盖 7a 中的)

    # ========== 输出 ==========
    out_dir           = NULL,        # [可选] 输出目录

    # ========== 核心泛化参数 ==========
    focus_tf_family   = "none",      # TF 家族聚焦 (逗号分隔多家族)
    focus_gene        = "none",      # 突变基因 ID 聚焦

    # ========== 富集检验 ==========
    enrichment_method = "fisher",    # fisher / hypergeometric
    min_edges         = 5,           # 最小边数阈值
    edge_classes      = "Lost,Gained,Rewired,Conserved",

    # ========== agriGO ==========
    background_mode   = "mfuzz",     # mfuzz / all_deg / genome
    gene_id_format    = "auto",

    # ========== MECS ==========
    w_diffcoexp       = 1.0,         # S1 权重
    w_cluster         = 1.0,         # S2 权重
    w_pcre            = 1.5,         # S3 权重
    output_tiers      = "1,2,3,4",   # 输出的 Tier 级别

    # ========== 调试 ==========
    verbose           = TRUE
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
        } else {
          params[[key]] <- value
        }
      }
    }
    i <- i + 1
  }

  # 解析逗号分隔参数
  params$edge_classes_vec <- trimws(unlist(strsplit(params$edge_classes, ",")))
  params$output_tiers_vec <- as.integer(trimws(unlist(strsplit(params$output_tiers, ","))))

  if (params$focus_tf_family != "none") {
    params$focus_tf_family_vec <- trimws(unlist(strsplit(params$focus_tf_family, ",")))
  } else {
    params$focus_tf_family_vec <- character(0)
  }

  return(params)
}

PARAMS <- parse_arguments()

################################################################################
# 验证参数
################################################################################

if (is.null(PARAMS$rdata_7a)) {
  cat("================================================================================\n")
  cat("  Pipeline 7b: Cluster-Constrained Discovery + MECS (v1.0)\n")
  cat("================================================================================\n\n")
  cat("错误: 必须指定 --rdata_7a\n\n")
  cat("使用方法:\n")
  cat("  Rscript RNAseq_pipeline7b_discovery.R \\\n")
  cat("    --rdata_7a /path/to/Pipeline7a_DiffCoexp_Results.RData \\\n")
  cat("    --focus_tf_family DOF --focus_gene Cre06.g250050\n")
  quit(status = 1)
}

################################################################################
# 主程序
################################################################################

cat("================================================================================\n")
cat("  Pipeline 7b: Cluster-Constrained Discovery + MECS Scoring\n")
cat("  v1.0 - Generalized Engineering Framework\n")
cat("================================================================================\n\n")
cat(sprintf("开始时间: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# ========================================================================
# Part 0: 加载数据 + 初始化
# ========================================================================

cat("--- Part 0: 加载 Pipeline 7a 数据 ---\n")

# 注意：load() 会覆盖名为 PARAMS 的对象 (来自 7a), 所以必须在 load 之前备份
PARAMS_7B_BACKUP <- PARAMS
rdata_7a_path <- PARAMS$rdata_7a

load(rdata_7a_path)
cat(sprintf("  Loaded 7a RData: %s\n", basename(rdata_7a_path)))

# 此时全局环境中的 PARAMS 是从 7a 继承来的
if (exists("PARAMS") && !identical(PARAMS, PARAMS_7B_BACKUP)) {
  PARAMS_7a_inherited <- PARAMS    # 存下 7a 的原始参数
  PARAMS <- PARAMS_7B_BACKUP       # 恢复 7b 自己的参数
} else {
  PARAMS_7a_inherited <- list()
}

# 输出目录
if (is.null(PARAMS$out_dir)) {
  PARAMS$out_dir <- dirname(dirname(PARAMS$rdata_7a))
}
OUT_DIR <- PARAMS$out_dir
OUT_7B <- file.path(OUT_DIR, "Pipeline7b")
dir.create(OUT_7B, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "rdata"), recursive = TRUE, showWarnings = FALSE)

# 日志
log_file <- file.path(OUT_7B, "Pipeline7b_discovery.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

# 更新 highlight_map (如果 7b 提供了新的)
if (!is.null(PARAMS$highlight_file) && file.exists(PARAMS$highlight_file)) {
  hl_df <- read.table(PARAMS$highlight_file, header = TRUE, sep = "\t",
                       stringsAsFactors = FALSE, quote = "")
  if (!"Gene" %in% colnames(hl_df)) {
    for (alt in c("GeneID", "Gene_ID", "gene_id", "gene")) {
      if (alt %in% colnames(hl_df)) { colnames(hl_df)[colnames(hl_df) == alt] <- "Gene"; break }
    }
  }
  if (!"Symbol" %in% colnames(hl_df)) {
    for (alt in c("symbol", "Name", "name")) {
      if (alt %in% colnames(hl_df)) { colnames(hl_df)[colnames(hl_df) == alt] <- "Symbol"; break }
    }
  }
  highlight_map <- setNames(hl_df$Symbol, hl_df$Gene)
  cat(sprintf("  Highlight genes updated: %d\n", length(highlight_map)))
}

# 获取 cluster 信息
has_cluster_info <- !is.null(cluster_annotation) && nrow(cluster_annotation) > 0
has_fate_info <- !is.null(fate_annotation) && nrow(fate_annotation) > 0

if (has_cluster_info) {
  n_wt_clusters <- max(cluster_annotation$WT_Cluster, na.rm = TRUE)
  n_mut_clusters <- max(cluster_annotation$MUT_Cluster, na.rm = TRUE)
  cat(sprintf("  Cluster info: %d WT clusters, %d MUT clusters\n",
              n_wt_clusters, n_mut_clusters))
} else {
  n_wt_clusters <- 0
  n_mut_clusters <- 0
  cat("  ⚠ No cluster annotation available\n")
}

if (has_fate_info) {
  cat(sprintf("  Fate annotation: %d genes\n", nrow(fate_annotation)))
}

# ========================================================================
# Part 0b: 构建 Cluster Pair Info
# ========================================================================

cat("\n--- Part 0b: 构建 Cluster Pair 信息 ---\n")

cluster_pair_info <- NULL

if (!is.null(PARAMS$cluster_pairs) && file.exists(PARAMS$cluster_pairs)) {
  # 模式 1: 从外部文件读取
  cluster_pair_info <- read.table(PARAMS$cluster_pairs, header = TRUE, sep = "\t",
                                   stringsAsFactors = FALSE)
  cat(sprintf("  Loaded cluster pairs from file: %d pairs\n", nrow(cluster_pair_info)))

} else if (has_cluster_info) {
  # 模式 2: 自动推断
  cat("  Auto-inferring cluster pairs from P5 annotation...\n")

  # 构建 confusion matrix
  conf_mat <- table(
    WT = cluster_annotation$WT_Cluster,
    MUT = cluster_annotation$MUT_Cluster
  )

  # 对每个 WT cluster 取 MUT cluster 中 argmax
  pairs_list <- list()
  for (wt_c in 1:n_wt_clusters) {
    if (as.character(wt_c) %in% rownames(conf_mat)) {
      row_vals <- conf_mat[as.character(wt_c), ]
      best_mut <- as.integer(names(which.max(row_vals)))
      overlap_n <- max(row_vals)
      total_n <- sum(row_vals)

      pairs_list[[length(pairs_list) + 1]] <- data.frame(
        WT_Cluster = wt_c,
        MUT_Cluster = best_mut,
        Overlap_Genes = overlap_n,
        Total_WT_Genes = total_n,
        Overlap_Ratio = round(overlap_n / total_n, 3),
        stringsAsFactors = FALSE
      )
    }
  }

  cluster_pair_info <- do.call(rbind, pairs_list)

  # 计算 divergence ratio (从 fate_annotation)
  if (has_fate_info) {
    div_ratios <- fate_annotation %>%
      group_by(WT_Cluster_P5b) %>%
      summarise(
        n_subclusters_total = n_distinct(SubCluster),
        n_subclusters_diverged = n_distinct(SubCluster[Fate_Label != "Conserved"]),
        .groups = "drop"
      ) %>%
      mutate(Divergence_Ratio = round(n_subclusters_diverged / n_subclusters_total, 3))

    cluster_pair_info <- cluster_pair_info %>%
      left_join(div_ratios, by = c("WT_Cluster" = "WT_Cluster_P5b"))

    cat(sprintf("  Divergence ratios calculated from fate annotation\n"))
  } else {
    cluster_pair_info$n_subclusters_total <- NA
    cluster_pair_info$n_subclusters_diverged <- NA
    cluster_pair_info$Divergence_Ratio <- NA
    cat("  ⚠ No fate annotation → Divergence_Ratio = NA for all clusters\n")
  }

  cat(sprintf("  Auto-inferred %d cluster pairs:\n", nrow(cluster_pair_info)))
  for (i in 1:nrow(cluster_pair_info)) {
    cp <- cluster_pair_info[i, ]
    div_str <- if (!is.na(cp$Divergence_Ratio)) sprintf("div=%.2f", cp$Divergence_Ratio) else "div=NA"
    cat(sprintf("    WT_C%d → MUT_C%d (overlap=%.1f%%, %s)\n",
                cp$WT_Cluster, cp$MUT_Cluster,
                cp$Overlap_Ratio * 100, div_str))
  }
} else {
  cat("  ⚠ No cluster pair info available (no P5 data or external file)\n")
}

# 保存 cluster pair info
if (!is.null(cluster_pair_info)) {
  write.table(cluster_pair_info,
              file.path(OUT_7B, "Pipeline7b_ModA_ClusterPair_Info.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}

# ========================================================================
# Part 0c: 解析 Focus 参数
# ========================================================================

cat("\n--- Part 0c: 解析 Focus 参数 ---\n")

focus_tfs <- character(0)
focus_gene_is_tf <- FALSE
focus_gene_is_target <- FALSE

# Focus TF Family
if (length(PARAMS$focus_tf_family_vec) > 0) {
  for (fam in PARAMS$focus_tf_family_vec) {
    fam_tfs <- names(tf_family_map)[tf_family_map == fam]
    cat(sprintf("  Focus family '%s': %d TFs in dataset\n", fam, length(fam_tfs)))
    focus_tfs <- union(focus_tfs, fam_tfs)
  }
}

# Focus Gene
if (PARAMS$focus_gene != "none") {
  fg <- PARAMS$focus_gene
  if (fg %in% tf_in_data) {
    focus_gene_is_tf <- TRUE
    focus_tfs <- union(focus_tfs, fg)
    # 自动将其家族加入 focus 列表
    fg_family <- tf_family_map[fg]
    if (!is.na(fg_family) && !(fg_family %in% PARAMS$focus_tf_family_vec)) {
      PARAMS$focus_tf_family_vec <- c(PARAMS$focus_tf_family_vec, fg_family)
      cat(sprintf("  Focus gene %s is a TF (family: %s) → added to focus\n", fg, fg_family))
    }
  } else if (fg %in% tg_in_data) {
    focus_gene_is_target <- TRUE
    cat(sprintf("  Focus gene %s is a target → will analyze upstream regulators\n", fg))
  } else {
    cat(sprintf("  ⚠ Focus gene %s not found in 7a analysis set\n", fg))
  }
}

has_focus <- length(focus_tfs) > 0 || focus_gene_is_target
cat(sprintf("  Total focus TFs: %d\n", length(focus_tfs)))
cat(sprintf("  Module D will %s\n", ifelse(has_focus, "run", "be skipped")))

# ========================================================================
# Part 0d: 加载 pCRE 信息 (可选)
# ========================================================================

pcre_map <- NULL
use_pcre <- FALSE

if (!is.null(PARAMS$pcre_tf_map) && file.exists(PARAMS$pcre_tf_map)) {
  pcre_map <- read.table(PARAMS$pcre_tf_map, header = TRUE, sep = "\t",
                          stringsAsFactors = FALSE)
  use_pcre <- TRUE
  cat(sprintf("\n  pCRE-TF map loaded: %d entries\n", nrow(pcre_map)))
} else {
  cat(sprintf("\n  pCRE-TF map not available → MECS degrades to S1+S2\n"))
  PARAMS$w_pcre <- 0
}

# ========================================================================
# Part 0e: 为边表添加 cluster pair 信息
# ========================================================================

cat("\n--- Part 0e: 增强边表 (cluster pair + membership) ---\n")

if (has_cluster_info && !is.null(cluster_pair_info)) {
  # 为每个 target 添加其所在 cluster pair 的 divergence ratio
  div_lookup <- setNames(cluster_pair_info$Divergence_Ratio, cluster_pair_info$WT_Cluster)

  if (!"Target_WT_Cluster" %in% colnames(edge_table)) {
    edge_table <- edge_table %>%
      left_join(cluster_annotation %>% select(Gene_ID, WT_Cluster, MUT_Cluster) %>%
                  rename(Target_WT_Cluster = WT_Cluster, Target_MUT_Cluster = MUT_Cluster),
                by = c("Target_Gene" = "Gene_ID"))
  }

  edge_table$Target_Divergence_Ratio <- div_lookup[as.character(edge_table$Target_WT_Cluster)]

  # 添加 membership 值
  mem_lookup <- setNames(cluster_annotation$WT_Membership, cluster_annotation$Gene_ID)
  edge_table$Target_WT_Membership <- mem_lookup[edge_table$Target_Gene]

  cat(sprintf("  Added Divergence_Ratio and Membership to edge table\n"))
} else {
  edge_table$Target_Divergence_Ratio <- NA_real_
  edge_table$Target_WT_Membership <- NA_real_
  cat("  ⚠ Cluster pair info not available → S2 will be uniform\n")
}


# ########################################################################
# MODULE A: Cluster-Constrained TF Hub Discovery
# ########################################################################

cat("\n========================================================================\n")
cat("  Module A: Cluster-Constrained TF Hub Discovery\n")
cat("========================================================================\n")

run_module_A <- function(edge_tbl, cluster_info, tf_fam_map,
                          focus_fams, edge_classes, min_edges, method) {

  # ---- Step 1: Per-TF per-Cluster edge count ----
  cat("\n  [A.1] Computing per-TF per-Cluster edge counts...\n")

  # 过滤有 cluster 信息的边
  tbl_with_cl <- edge_tbl %>%
    filter(!is.na(Target_WT_Cluster))

  if (nrow(tbl_with_cl) == 0) {
    cat("  ⚠ No edges with cluster info → Module A skipped\n")
    return(list(hub_table = NULL, family_enrich = NULL, focus_report = NULL))
  }

  # 对每个 (TF, Cluster, EdgeClass) 组合计数
  hub_counts <- tbl_with_cl %>%
    filter(edge_class %in% edge_classes) %>%
    group_by(TF_Gene, TF_Family, Target_WT_Cluster, edge_class) %>%
    summarise(n_edges = n(), .groups = "drop")

  cat(sprintf("    %d (TF, Cluster, EdgeClass) combinations\n", nrow(hub_counts)))

  # ---- Step 2: 超几何检验 ----
  cat("  [A.2] Running hypergeometric tests...\n")

  # 全局统计
  global_counts <- tbl_with_cl %>%
    filter(edge_class %in% edge_classes) %>%
    group_by(TF_Gene, edge_class) %>%
    summarise(tf_total = n(), .groups = "drop")

  cluster_totals <- tbl_with_cl %>%
    filter(edge_class %in% edge_classes) %>%
    group_by(Target_WT_Cluster, edge_class) %>%
    summarise(cluster_total = n(), .groups = "drop")

  overall_totals <- tbl_with_cl %>%
    filter(edge_class %in% edge_classes) %>%
    group_by(edge_class) %>%
    summarise(grand_total = n(), .groups = "drop")

  # 合并所有信息
  hub_test <- hub_counts %>%
    left_join(global_counts, by = c("TF_Gene", "edge_class")) %>%
    left_join(cluster_totals, by = c("Target_WT_Cluster", "edge_class")) %>%
    left_join(overall_totals, by = "edge_class")

  # 超几何检验
  hub_test$p_value <- with(hub_test, {
    mapply(function(q, m, n_black, k) {
      if (any(is.na(c(q, m, n_black, k)))) return(NA_real_)
      if (q < 1 || m < 1) return(1.0)
      phyper(q - 1, m, n_black, k, lower.tail = FALSE)
    }, n_edges, tf_total, grand_total - tf_total, cluster_total)
  })

  # BH 联合校正
  hub_test$fdr <- p.adjust(hub_test$p_value, method = "BH")

  n_sig_hub <- sum(hub_test$fdr < 0.05, na.rm = TRUE)
  cat(sprintf("    Significant (FDR<0.05): %d / %d tests\n", n_sig_hub, nrow(hub_test)))

  # ---- Step 3: TF Family 级 Fisher 检验 ----
  cat("  [A.3] TF Family-level Fisher tests...\n")

  family_counts <- tbl_with_cl %>%
    filter(edge_class %in% edge_classes) %>%
    mutate(TF_Family = tf_fam_map[TF_Gene]) %>%
    filter(!is.na(TF_Family)) %>%
    group_by(TF_Family, Target_WT_Cluster, edge_class) %>%
    summarise(family_edge_count = n(), .groups = "drop")

  # 对每个 (Family, Cluster): 该家族 vs 其他家族 在此 edge_class 上的 Fisher
  family_total_by_cluster <- tbl_with_cl %>%
    filter(edge_class %in% edge_classes) %>%
    group_by(Target_WT_Cluster, edge_class) %>%
    summarise(total_edges = n(), .groups = "drop")

  family_enrich <- family_counts %>%
    left_join(family_total_by_cluster,
              by = c("Target_WT_Cluster", "edge_class")) %>%
    mutate(other_edges = total_edges - family_edge_count)

  # 需要全局该家族非此 edge_class 的边数
  family_global <- tbl_with_cl %>%
    mutate(TF_Family = tf_fam_map[TF_Gene]) %>%
    filter(!is.na(TF_Family)) %>%
    group_by(TF_Family, Target_WT_Cluster) %>%
    summarise(family_all_edges = n(), .groups = "drop")

  family_enrich <- family_enrich %>%
    left_join(family_global, by = c("TF_Family", "Target_WT_Cluster")) %>%
    mutate(family_non_class = family_all_edges - family_edge_count)

  cluster_all_totals <- tbl_with_cl %>%
    group_by(Target_WT_Cluster) %>%
    summarise(cluster_all_total = n(), .groups = "drop")

  family_enrich <- family_enrich %>%
    left_join(cluster_all_totals, by = "Target_WT_Cluster") %>%
    mutate(other_non_class = cluster_all_total - family_all_edges - other_edges + family_edge_count)

  # Fisher 检验 (只对够边的)
  family_enrich$fisher_p <- NA_real_
  family_enrich$fisher_OR <- NA_real_

  for (i in 1:nrow(family_enrich)) {
    fe <- family_enrich[i, ]
    if (fe$family_edge_count < min_edges) next

    a <- fe$family_edge_count
    b <- max(0, fe$family_non_class)
    c <- max(0, fe$other_edges)
    d <- max(0, fe$other_non_class)

    if (any(is.na(c(a, b, c, d)))) next

    mat <- matrix(c(a, c, b, d), nrow = 2)
    ft <- tryCatch(fisher.test(mat), error = function(e) NULL)
    if (!is.null(ft)) {
      family_enrich$fisher_p[i] <- ft$p.value
      family_enrich$fisher_OR[i] <- ft$estimate
    }
  }

  family_enrich$fisher_fdr <- p.adjust(family_enrich$fisher_p, method = "BH")
  family_enrich$log2OR <- log2(pmax(family_enrich$fisher_OR, 1e-10))

  n_sig_fam <- sum(family_enrich$fisher_fdr < 0.05, na.rm = TRUE)
  cat(sprintf("    Family enrichment sig (FDR<0.05): %d / %d tests\n",
              n_sig_fam, sum(!is.na(family_enrich$fisher_p))))

  # ---- Step 4: Focus Family Report ----
  focus_report <- NULL
  if (length(focus_fams) > 0 && nrow(hub_test) > 0) {
    cat("  [A.4] Generating Focus Family Report...\n")

    focus_report <- hub_test %>%
      filter(TF_Family %in% focus_fams) %>%
      arrange(edge_class, fdr) %>%
      mutate(Is_Significant = fdr < 0.05)

    cat(sprintf("    Focus families: %s → %d test results\n",
                paste(focus_fams, collapse = ", "), nrow(focus_report)))
  }

  return(list(hub_table = hub_test, family_enrich = family_enrich,
              focus_report = focus_report))
}

mod_a <- run_module_A(
  edge_tbl = edge_table,
  cluster_info = cluster_pair_info,
  tf_fam_map = tf_family_map,
  focus_fams = PARAMS$focus_tf_family_vec,
  edge_classes = PARAMS$edge_classes_vec,
  min_edges = PARAMS$min_edges,
  method = PARAMS$enrichment_method
)

# 保存 Module A 输出
if (!is.null(mod_a$hub_table)) {
  write.table(mod_a$hub_table,
              file.path(OUT_7B, "Pipeline7b_ModA_TF_Hub_ByCluster.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
if (!is.null(mod_a$family_enrich)) {
  write.table(mod_a$family_enrich,
              file.path(OUT_7B, "Pipeline7b_ModA_TFFamily_Enrichment.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
if (!is.null(mod_a$focus_report)) {
  for (fam in PARAMS$focus_tf_family_vec) {
    fam_report <- mod_a$focus_report %>% filter(TF_Family == fam)
    if (nrow(fam_report) > 0) {
      write.table(fam_report,
                  file.path(OUT_7B, sprintf("Pipeline7b_ModA_Focus_%s_Report.tsv", fam)),
                  sep = "\t", quote = FALSE, row.names = FALSE)
    }
  }
}


# ########################################################################
# MODULE B: agriGO Gene List Export
# ########################################################################

cat("\n========================================================================\n")
cat("  Module B: agriGO Gene List Export\n")
cat("========================================================================\n")

run_module_B <- function(edge_tbl, cluster_info, focus_tfs_list, focus_fams,
                          tf_fam_map, background_mode,
                          mfuzz_genes, tg_full, out_dir) {

  agriGO_dir <- file.path(out_dir, "agriGO_lists")
  dir.create(file.path(agriGO_dir, "background"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(agriGO_dir, "by_cluster"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(agriGO_dir, "by_edge_class"), recursive = TRUE, showWarnings = FALSE)

  # ---- 背景基因集 ----
  cat("  [B.1] Exporting background gene set...\n")
  bg_genes <- switch(background_mode,
    "mfuzz" = if (!is.null(mfuzz_genes)) mfuzz_genes else tg_full,
    "all_deg" = tg_full,
    "genome" = NULL
  )

  if (!is.null(bg_genes) && length(bg_genes) > 0) {
    writeLines(bg_genes, file.path(agriGO_dir, "background", "background_genes.txt"))
    cat(sprintf("    Background: %s mode, %d genes\n", background_mode, length(bg_genes)))
  } else {
    cat("    Background: genome mode (no file, use species default in agriGO)\n")
  }

  # ---- 按 Edge Class 导出 (全基因组) ----
  cat("  [B.2] Exporting by edge class (genome-wide)...\n")

  for (ec in unique(edge_tbl$edge_class)) {
    genes <- unique(edge_tbl$Target_Gene[edge_tbl$edge_class == ec])
    if (length(genes) >= 1) {
      writeLines(genes, file.path(agriGO_dir, "by_edge_class",
                                   sprintf("all_%s_targets.txt", ec)))
      cat(sprintf("    %s: %d targets\n", ec, length(genes)))
    }
  }

  # ---- 按 Cluster 导出 ----
  if (!is.null(cluster_info) && "Target_WT_Cluster" %in% colnames(edge_tbl)) {
    cat("  [B.3] Exporting by cluster...\n")

    wt_clusters <- sort(unique(edge_tbl$Target_WT_Cluster[!is.na(edge_tbl$Target_WT_Cluster)]))

    for (cl in wt_clusters) {
      cl_dir <- file.path(agriGO_dir, "by_cluster", sprintf("WT_C%d", cl))
      dir.create(cl_dir, showWarnings = FALSE)

      cl_edges <- edge_tbl %>% filter(Target_WT_Cluster == cl)

      for (ec in unique(cl_edges$edge_class)) {
        genes <- unique(cl_edges$Target_Gene[cl_edges$edge_class == ec])
        if (length(genes) >= 1) {
          writeLines(genes, file.path(cl_dir, sprintf("%s_targets.txt", ec)))
        }
      }

      # All cluster genes
      all_cl_genes <- unique(cl_edges$Target_Gene)
      if (length(all_cl_genes) >= 1) {
        writeLines(all_cl_genes, file.path(cl_dir, "all_cluster_genes.txt"))
      }
    }
    cat(sprintf("    Exported lists for %d WT clusters\n", length(wt_clusters)))
  }

  # ---- 按 Focus TF 导出 ----
  if (length(focus_tfs_list) > 0 && length(focus_fams) > 0) {
    cat("  [B.4] Exporting by focus TF...\n")
    focus_dir <- file.path(agriGO_dir, "by_focus_tf")
    dir.create(focus_dir, showWarnings = FALSE)

    for (fam in focus_fams) {
      fam_tfs <- names(tf_fam_map)[tf_fam_map == fam]
      fam_tfs <- intersect(fam_tfs, focus_tfs_list)

      # 家族全体 Lost targets
      fam_lost <- unique(edge_tbl$Target_Gene[
        edge_tbl$TF_Gene %in% fam_tfs & edge_tbl$edge_class == "Lost"
      ])
      if (length(fam_lost) >= 1) {
        writeLines(fam_lost, file.path(focus_dir, sprintf("%s_lost_targets.txt", fam)))
        cat(sprintf("    %s family Lost targets: %d\n", fam, length(fam_lost)))
      }

      # 单个 TF
      for (tf in fam_tfs) {
        tf_lost <- unique(edge_tbl$Target_Gene[
          edge_tbl$TF_Gene == tf & edge_tbl$edge_class == "Lost"
        ])
        if (length(tf_lost) >= 1) {
          tf_label <- gsub("[^A-Za-z0-9._]", "_", tf)
          writeLines(tf_lost, file.path(focus_dir,
                                         sprintf("%s_%s_lost_targets.txt", fam, tf_label)))
        }
      }
    }
  }

  # ---- README ----
  readme_lines <- c(
    sprintf("# agriGO Gene Lists - Pipeline 7b"),
    sprintf("# Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
    "",
    sprintf("Background gene set: %s (%s genes)",
            background_mode,
            ifelse(!is.null(bg_genes), as.character(length(bg_genes)), "genome default")),
    "",
    "Gene ID format: auto-detected from upstream data",
    "",
    "Directory structure:",
    "  background/     - Background gene set for enrichment analysis",
    "  by_edge_class/  - Gene lists by edge classification (genome-wide)",
    "  by_cluster/     - Gene lists organized by WT cluster",
    "  by_focus_tf/    - Gene lists for focus TF family (if specified)",
    "",
    "Usage in agriGO:",
    "  1. Select 'Singular Enrichment Analysis (SEA)'",
    "  2. Upload target gene list from appropriate subdirectory",
    "  3. Upload background file from background/ directory",
    "  4. Select your species and gene ID type",
    "",
    "Note: Files with < 10 genes may not produce meaningful enrichment results."
  )
  writeLines(readme_lines, file.path(agriGO_dir, "README.txt"))
  cat("  README.txt generated\n")

  invisible(NULL)
}

run_module_B(
  edge_tbl = edge_table,
  cluster_info = cluster_pair_info,
  focus_tfs_list = focus_tfs,
  focus_fams = PARAMS$focus_tf_family_vec,
  tf_fam_map = tf_family_map,
  background_mode = PARAMS$background_mode,
  mfuzz_genes = mfuzz_common_genes,
  tg_full = tg_in_data_full,
  out_dir = OUT_7B
)


# ########################################################################
# MODULE C: TF Family × Cluster × Edge Class Enrichment Matrix
# ########################################################################

cat("\n========================================================================\n")
cat("  Module C: TF Family × Cluster Enrichment Matrix\n")
cat("========================================================================\n")

run_module_C <- function(edge_tbl, tf_fam_map, edge_classes, min_edges) {

  tbl <- edge_tbl %>%
    filter(!is.na(Target_WT_Cluster)) %>%
    mutate(TF_Family = tf_fam_map[TF_Gene]) %>%
    filter(!is.na(TF_Family))

  if (nrow(tbl) == 0) {
    cat("  ⚠ No edges with cluster + family info → Module C skipped\n")
    return(list(enrich_long = NULL, enrich_wide_list = NULL))
  }

  results_list <- list()

  for (ec in edge_classes) {
    cat(sprintf("  [C] Enrichment for edge class: %s\n", ec))

    families <- sort(unique(tbl$TF_Family))
    clusters <- sort(unique(tbl$Target_WT_Cluster))

    for (fam in families) {
      for (cl in clusters) {
        # 2×2 table: (this family vs others) × (this edge_class vs others)
        a <- sum(tbl$TF_Family == fam & tbl$Target_WT_Cluster == cl & tbl$edge_class == ec)
        b <- sum(tbl$TF_Family == fam & tbl$Target_WT_Cluster == cl & tbl$edge_class != ec)
        c_val <- sum(tbl$TF_Family != fam & tbl$Target_WT_Cluster == cl & tbl$edge_class == ec)
        d <- sum(tbl$TF_Family != fam & tbl$Target_WT_Cluster == cl & tbl$edge_class != ec)

        if (a < min_edges) {
          results_list[[length(results_list) + 1]] <- data.frame(
            TF_Family = fam, WT_Cluster = cl, Edge_Class = ec,
            Count = a, fisher_OR = NA, fisher_p = NA,
            stringsAsFactors = FALSE
          )
          next
        }

        mat <- matrix(c(a, c_val, b, d), nrow = 2)
        ft <- tryCatch(fisher.test(mat), error = function(e) NULL)

        results_list[[length(results_list) + 1]] <- data.frame(
          TF_Family = fam, WT_Cluster = cl, Edge_Class = ec,
          Count = a,
          fisher_OR = if (!is.null(ft)) ft$estimate else NA,
          fisher_p = if (!is.null(ft)) ft$p.value else NA,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  enrich_long <- do.call(rbind, results_list)

  # BH across all tests
  enrich_long$fisher_fdr <- p.adjust(enrich_long$fisher_p, method = "BH")
  enrich_long$log2OR <- log2(pmax(enrich_long$fisher_OR, 1e-10))
  enrich_long$neg_log10_fdr <- -log10(pmax(enrich_long$fisher_fdr, 1e-300))
  enrich_long$neg_log10_fdr <- pmin(enrich_long$neg_log10_fdr, 10)  # cap at 10

  # 显著性标记
  enrich_long$sig_mark <- ifelse(
    enrich_long$fisher_fdr < 0.01, "**",
    ifelse(enrich_long$fisher_fdr < 0.05, "*", "")
  )

  n_sig <- sum(enrich_long$fisher_fdr < 0.05, na.rm = TRUE)
  cat(sprintf("  Total enrichment tests: %d | Significant: %d\n", nrow(enrich_long), n_sig))

  # 宽格式 (per edge class)
  enrich_wide_list <- list()
  for (ec in edge_classes) {
    ec_data <- enrich_long %>% filter(Edge_Class == ec)
    if (nrow(ec_data) == 0) next

    # 全基因组该 edge class 总数 per family (用于排序)
    fam_total <- edge_tbl %>%
      filter(edge_class == ec) %>%
      mutate(TF_Family = tf_fam_map[TF_Gene]) %>%
      filter(!is.na(TF_Family)) %>%
      count(TF_Family, name = "global_count")

    wide_or <- ec_data %>%
      select(TF_Family, WT_Cluster, log2OR) %>%
      pivot_wider(names_from = WT_Cluster, values_from = log2OR,
                  names_prefix = sprintf("C%s_log2OR_", ""))

    wide_fdr <- ec_data %>%
      select(TF_Family, WT_Cluster, fisher_fdr) %>%
      pivot_wider(names_from = WT_Cluster, values_from = fisher_fdr,
                  names_prefix = sprintf("C%s_FDR_", ""))

    wide <- wide_or %>%
      left_join(fam_total, by = "TF_Family") %>%
      arrange(desc(global_count))

    enrich_wide_list[[ec]] <- wide
  }

  return(list(enrich_long = enrich_long, enrich_wide_list = enrich_wide_list))
}

mod_c <- run_module_C(
  edge_tbl = edge_table,
  tf_fam_map = tf_family_map,
  edge_classes = PARAMS$edge_classes_vec,
  min_edges = PARAMS$min_edges
)

# 保存 Module C 输出
if (!is.null(mod_c$enrich_long)) {
  write.table(mod_c$enrich_long,
              file.path(OUT_7B, "Pipeline7b_ModC_Enrichment_Long.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
if (!is.null(mod_c$enrich_wide_list)) {
  for (ec in names(mod_c$enrich_wide_list)) {
    write.table(mod_c$enrich_wide_list[[ec]],
                file.path(OUT_7B, sprintf("Pipeline7b_ModC_Enrichment_Wide_%s.tsv", ec)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}


# ########################################################################
# MODULE D: Focus TF Rewiring Profiles
# ########################################################################

cat("\n========================================================================\n")
cat("  Module D: Focus TF Rewiring Profiles\n")
cat("========================================================================\n")

run_module_D <- function(edge_tbl, focus_tfs_list, focus_gene, focus_gene_is_tf,
                          focus_gene_is_target, cluster_pair_info,
                          tf_fam_map, highlight_map_local) {

  if (length(focus_tfs_list) == 0 && !focus_gene_is_target) {
    cat("  Module D skipped (no focus specified)\n")
    return(list(profiles = NULL, details = NULL, upstream = NULL))
  }

  profiles_list <- list()
  detail_tables <- list()

  # ---- Focus TFs: 下游调控变化分析 ----
  if (length(focus_tfs_list) > 0) {
    cat(sprintf("  [D.1] Analyzing %d focus TFs as regulators...\n", length(focus_tfs_list)))

    for (tf in focus_tfs_list) {
      tf_edges <- edge_tbl %>% filter(TF_Gene == tf)

      if (nrow(tf_edges) == 0) next

      # 各类边统计
      n_by_class <- tf_edges %>%
        count(edge_class, name = "n") %>%
        pivot_wider(names_from = edge_class, values_from = n, values_fill = 0)

      # Cluster 分布
      lost_by_cluster <- tf_edges %>%
        filter(edge_class == "Lost", !is.na(Target_WT_Cluster)) %>%
        count(Target_WT_Cluster, name = "lost_count")

      # Cluster Specificity Index (CSI)
      total_lost <- sum(tf_edges$edge_class == "Lost")
      total_targets <- nrow(tf_edges %>% filter(sig_WT | sig_MUT))

      if (total_lost > 0 && total_targets > 0 && !is.null(cluster_pair_info)) {
        diverged_clusters <- cluster_pair_info$WT_Cluster[
          !is.na(cluster_pair_info$Divergence_Ratio) & cluster_pair_info$Divergence_Ratio > 0.5
        ]
        lost_in_diverged <- sum(tf_edges$edge_class == "Lost" &
                                 tf_edges$Target_WT_Cluster %in% diverged_clusters, na.rm = TRUE)
        targets_in_diverged <- sum((tf_edges$sig_WT | tf_edges$sig_MUT) &
                                    tf_edges$Target_WT_Cluster %in% diverged_clusters, na.rm = TRUE)

        csi <- if (targets_in_diverged > 0 && total_targets > 0) {
          (lost_in_diverged / total_lost) / (targets_in_diverged / total_targets)
        } else { NA_real_ }
      } else {
        csi <- NA_real_
      }

      # Direction consistency
      lost_directions <- tf_edges %>%
        filter(edge_class == "Lost") %>%
        pull(direction_WT)
      n_positive <- sum(lost_directions == "positive")
      n_negative <- sum(lost_directions == "negative")
      direction_consistency <- max(n_positive, n_negative) / max(1, n_positive + n_negative)

      # Rewiring magnitude
      rewiring_mag <- mean(abs(tf_edges$delta_r[tf_edges$edge_class %in%
                               c("Lost", "Gained", "Rewired")]), na.rm = TRUE)

      # Symbol
      tf_sym <- if (!is.null(highlight_map_local) && tf %in% names(highlight_map_local)) {
        highlight_map_local[tf]
      } else { NA_character_ }

      profiles_list[[tf]] <- data.frame(
        TF_Gene = tf,
        TF_Symbol = tf_sym,
        TF_Family = tf_fam_map[tf],
        n_targets_WT = sum(tf_edges$sig_WT),
        n_targets_MUT = sum(tf_edges$sig_MUT),
        n_lost = sum(tf_edges$edge_class == "Lost"),
        n_gained = sum(tf_edges$edge_class == "Gained"),
        n_rewired = sum(tf_edges$edge_class == "Rewired"),
        n_conserved = sum(tf_edges$edge_class == "Conserved"),
        delta_out_degree = sum(tf_edges$sig_WT) - sum(tf_edges$sig_MUT),
        CSI = round(csi, 3),
        Direction_Consistency = round(direction_consistency, 3),
        Rewiring_Magnitude = round(rewiring_mag, 4),
        stringsAsFactors = FALSE
      )

      # Detail table (for this TF)
      detail_tables[[tf]] <- tf_edges %>%
        filter(edge_class %in% c("Lost", "Gained", "Rewired", "Conserved")) %>%
        arrange(edge_class, fdr_diff)
    }
  }

  # ---- Focus Gene as Target: 上游调控分析 ----
  upstream_table <- NULL

  if (focus_gene_is_target) {
    fg <- focus_gene
    cat(sprintf("  [D.2] Analyzing focus gene %s as a target (upstream regulators)...\n", fg))

    upstream_edges <- edge_tbl %>%
      filter(Target_Gene == fg)

    if (nrow(upstream_edges) > 0) {
      upstream_table <- upstream_edges %>%
        select(TF_Gene, TF_Family, TF_Symbol, edge_class,
               r_WT, r_MUT, delta_r, fdr_WT, fdr_MUT, fdr_diff,
               direction_WT, direction_MUT) %>%
        arrange(edge_class, fdr_diff)

      cat(sprintf("    Found %d upstream TF edges for %s\n", nrow(upstream_table), fg))

      # 上游 TF 家族分布
      upstream_fam_counts <- upstream_table %>%
        filter(edge_class %in% c("Lost", "Gained")) %>%
        count(TF_Family, edge_class) %>%
        arrange(desc(n))

      if (nrow(upstream_fam_counts) > 0) {
        cat("    Upstream TF family distribution (Lost/Gained):\n")
        for (i in 1:min(10, nrow(upstream_fam_counts))) {
          row <- upstream_fam_counts[i, ]
          cat(sprintf("      %s (%s): %d edges\n", row$TF_Family, row$edge_class, row$n))
        }
      }
    } else {
      cat(sprintf("    ⚠ No edges found for target %s\n", fg))
    }
  }

  # 合并 profiles
  profiles <- if (length(profiles_list) > 0) do.call(rbind, profiles_list) else NULL

  return(list(profiles = profiles, details = detail_tables, upstream = upstream_table))
}

mod_d <- run_module_D(
  edge_tbl = edge_table,
  focus_tfs_list = focus_tfs,
  focus_gene = PARAMS$focus_gene,
  focus_gene_is_tf = focus_gene_is_tf,
  focus_gene_is_target = focus_gene_is_target,
  cluster_pair_info = cluster_pair_info,
  tf_fam_map = tf_family_map,
  highlight_map_local = highlight_map
)

# 保存 Module D 输出
if (!is.null(mod_d$profiles)) {
  write.table(mod_d$profiles,
              file.path(OUT_7B, "Pipeline7b_ModD_Focus_TF_Profiles.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
if (length(mod_d$details) > 0) {
  for (tf_name in names(mod_d$details)) {
    tf_label <- gsub("[^A-Za-z0-9._]", "_", tf_name)
    write.table(mod_d$details[[tf_name]],
                file.path(OUT_7B, sprintf("Pipeline7b_ModD_%s_Detail.tsv", tf_label)),
                sep = "\t", quote = FALSE, row.names = FALSE)
  }
}
if (!is.null(mod_d$upstream)) {
  write.table(mod_d$upstream,
              file.path(OUT_7B, "Pipeline7b_ModD_FocusGene_Upstream.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
}


# ########################################################################
# MODULE E: MECS Scoring & Tier Classification
# ########################################################################

cat("\n========================================================================\n")
cat("  Module E: MECS Multi-Evidence Confidence Scoring\n")
cat("========================================================================\n")

run_module_E <- function(edge_tbl, cluster_pair_info, pcre_map, use_pcre,
                          w1, w2, w3, output_tiers) {

  cat("  [E.1] Computing S1 (differential co-expression strength)...\n")

  # S1: -log10(fdr_diff), normalized by 99th percentile
  s1_raw <- -log10(pmax(edge_tbl$fdr_diff, 1e-300))
  s1_raw[is.na(s1_raw)] <- 0
  s1_raw[!is.finite(s1_raw)] <- 0

  # 归一化: 99th percentile capping
  s1_cap <- quantile(s1_raw[s1_raw > 0], 0.99, na.rm = TRUE)
  if (is.na(s1_cap) || s1_cap == 0) s1_cap <- 1
  edge_tbl$S1 <- pmin(s1_raw / s1_cap, 1.0)

  cat(sprintf("    S1 range: [%.4f, %.4f], 99th percentile cap: %.2f\n",
              min(edge_tbl$S1), max(edge_tbl$S1), s1_cap))

  # S2: Cluster constraint (membership × divergence_ratio)
  cat("  [E.2] Computing S2 (cluster constraint)...\n")

  if (!all(is.na(edge_tbl$Target_Divergence_Ratio)) &&
      !all(is.na(edge_tbl$Target_WT_Membership))) {
    edge_tbl$S2 <- ifelse(
      !is.na(edge_tbl$Target_WT_Membership) & !is.na(edge_tbl$Target_Divergence_Ratio),
      edge_tbl$Target_WT_Membership * edge_tbl$Target_Divergence_Ratio,
      0
    )
    cat(sprintf("    S2 range: [%.4f, %.4f]\n",
                min(edge_tbl$S2, na.rm = TRUE), max(edge_tbl$S2, na.rm = TRUE)))
  } else {
    edge_tbl$S2 <- 0
    cat("    S2 = 0 for all edges (no cluster/divergence info)\n")
    if (w2 > 0) {
      cat("    ⚠ w_cluster set to 0 (no S2 data)\n")
      w2 <- 0
    }
  }

  # S3: pCRE family consistency (optional)
  cat("  [E.3] Computing S3 (pCRE consistency)...\n")

  if (use_pcre && !is.null(pcre_map)) {
    # pcre_map 格式: WT_Cluster | TF_Family | Importance_Rank
    edge_tbl$S3 <- 0

    # 对每条边检查: TF_Family 是否在 target 的 cluster 的 pCRE map 中
    for (i in 1:nrow(edge_tbl)) {
      tf_fam <- edge_tbl$TF_Family[i]
      tgt_cl <- edge_tbl$Target_WT_Cluster[i]
      if (is.na(tf_fam) || is.na(tgt_cl)) next

      match_rows <- pcre_map$WT_Cluster == tgt_cl & pcre_map$TF_Family == tf_fam
      if (any(match_rows)) {
        if ("Importance_Rank" %in% colnames(pcre_map)) {
          rank_val <- min(pcre_map$Importance_Rank[match_rows])
          max_rank <- max(pcre_map$Importance_Rank, na.rm = TRUE)
          edge_tbl$S3[i] <- 1 - rank_val / max_rank
        } else {
          edge_tbl$S3[i] <- 1.0
        }
      }
    }
    cat(sprintf("    S3 > 0: %d edges\n", sum(edge_tbl$S3 > 0)))
  } else {
    edge_tbl$S3 <- NA_real_
    w3 <- 0
    cat("    S3 = NA (no pCRE data)\n")
  }

  # MECS composite score
  cat("  [E.4] Computing MECS composite score...\n")

  active_weights <- w1 + w2 + w3
  if (active_weights == 0) active_weights <- 1  # 安全防护

  edge_tbl$MECS_score <- (
    w1 * edge_tbl$S1 +
    w2 * edge_tbl$S2 +
    w3 * ifelse(is.na(edge_tbl$S3), 0, edge_tbl$S3)
  ) / active_weights

  cat(sprintf("    MECS weights: S1=%.1f, S2=%.1f, S3=%.1f (active sum=%.1f)\n",
              w1, w2, w3, active_weights))
  cat(sprintf("    MECS range: [%.4f, %.4f]\n",
              min(edge_tbl$MECS_score, na.rm = TRUE),
              max(edge_tbl$MECS_score, na.rm = TRUE)))

  # Tier classification
  cat("  [E.5] Tier classification...\n")

  edge_tbl$Tier <- 5L  # default

  # Tier 1: 三层 (或两层, 若无 pCRE) 证据齐全
  if (use_pcre) {
    tier1_mask <- edge_tbl$edge_class %in% c("Lost", "Gained", "Rewired") &
                  edge_tbl$S2 > 0.5 &
                  !is.na(edge_tbl$S3) & edge_tbl$S3 > 0
    tier2_mask <- edge_tbl$edge_class %in% c("Lost", "Gained", "Rewired") &
                  edge_tbl$S2 > 0.5 &
                  (is.na(edge_tbl$S3) | edge_tbl$S3 == 0)
  } else {
    # 无 pCRE 时: Tier 1 = diffcoexp + cluster, Tier 2 无法区分
    tier1_mask <- edge_tbl$edge_class %in% c("Lost", "Gained", "Rewired") &
                  edge_tbl$S2 > 0.5
    tier2_mask <- rep(FALSE, nrow(edge_tbl))  # 合并到 Tier 1
  }

  tier3_mask <- edge_tbl$edge_class %in% c("Lost", "Gained", "Rewired") &
                edge_tbl$S2 <= 0.5

  tier4_mask <- edge_tbl$edge_class == "Conserved" &
                edge_tbl$S2 < 0.2

  edge_tbl$Tier[tier1_mask] <- 1L
  edge_tbl$Tier[tier2_mask] <- 2L
  edge_tbl$Tier[tier3_mask] <- 3L
  edge_tbl$Tier[tier4_mask] <- 4L

  tier_counts <- table(edge_tbl$Tier)
  cat("    Tier distribution:\n")
  for (t in sort(as.integer(names(tier_counts)))) {
    cat(sprintf("      Tier %d: %d edges\n", t, tier_counts[as.character(t)]))
  }

  # 过滤输出
  scored_edges <- edge_tbl %>%
    filter(Tier %in% output_tiers) %>%
    arrange(Tier, desc(MECS_score))

  cat(sprintf("    Output edges (Tiers %s): %d\n",
              paste(output_tiers, collapse = ","), nrow(scored_edges)))

  return(list(scored_edges = scored_edges, edge_table_full = edge_tbl,
              tier_counts = tier_counts, weights_used = c(S1 = w1, S2 = w2, S3 = w3)))
}

mod_e <- run_module_E(
  edge_tbl = edge_table,
  cluster_pair_info = cluster_pair_info,
  pcre_map = pcre_map,
  use_pcre = use_pcre,
  w1 = PARAMS$w_diffcoexp,
  w2 = PARAMS$w_cluster,
  w3 = PARAMS$w_pcre,
  output_tiers = PARAMS$output_tiers_vec
)

# 保存 Module E 输出
write.table(mod_e$scored_edges,
            file.path(OUT_7B, "Pipeline7b_MECS_Scored_Edges.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("\n  MECS scored edges saved: %d edges\n", nrow(mod_e$scored_edges)))

# MECS Summary
mecs_summary_lines <- c(
  "================================================================================",
  "  Pipeline 7b: MECS Scoring Summary",
  "================================================================================",
  "",
  sprintf("Analysis Time: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  sprintf("MECS Weights: S1=%.1f, S2=%.1f, S3=%.1f",
          mod_e$weights_used["S1"], mod_e$weights_used["S2"], mod_e$weights_used["S3"]),
  sprintf("pCRE evidence: %s", ifelse(use_pcre, "available", "not used (S3=NA)")),
  "",
  "--- Tier Distribution ---"
)

for (t in sort(as.integer(names(mod_e$tier_counts)))) {
  mecs_summary_lines <- c(mecs_summary_lines,
    sprintf("  Tier %d: %d edges", t, mod_e$tier_counts[as.character(t)]))
}

# Focus family in tiers
if (length(PARAMS$focus_tf_family_vec) > 0) {
  mecs_summary_lines <- c(mecs_summary_lines, "",
    "--- Focus Family Tier Distribution ---")
  for (fam in PARAMS$focus_tf_family_vec) {
    fam_tfs <- names(tf_family_map)[tf_family_map == fam]
    fam_edges <- mod_e$edge_table_full %>%
      filter(TF_Gene %in% fam_tfs)
    fam_tier <- table(fam_edges$Tier)
    mecs_summary_lines <- c(mecs_summary_lines,
      sprintf("  %s family:", fam))
    for (t in sort(as.integer(names(fam_tier)))) {
      mecs_summary_lines <- c(mecs_summary_lines,
        sprintf("    Tier %d: %d edges", t, fam_tier[as.character(t)]))
    }
  }
}

# Top MECS edges
mecs_summary_lines <- c(mecs_summary_lines, "",
  "--- Top 20 MECS Edges ---")
top20 <- mod_e$scored_edges %>%
  head(20)
for (i in 1:min(20, nrow(top20))) {
  row <- top20[i, ]
  tf_label <- ifelse(!is.na(row$TF_Symbol), row$TF_Symbol, row$TF_Gene)
  tg_label <- ifelse(!is.na(row$Target_Symbol), row$Target_Symbol, row$Target_Gene)
  mecs_summary_lines <- c(mecs_summary_lines,
    sprintf("  %2d. %s → %s | %s | MECS=%.3f | Tier %d",
            i, tf_label, tg_label, row$edge_class, row$MECS_score, row$Tier))
}

writeLines(mecs_summary_lines, file.path(OUT_7B, "Pipeline7b_MECS_Summary.txt"))


# ########################################################################
# SAVE: RData + Summary + Parameters JSON
# ########################################################################

cat("\n========================================================================\n")
cat("  Saving outputs\n")
cat("========================================================================\n")

# 全局摘要
summary_lines <- c(
  "================================================================================",
  "  Pipeline 7b: Cluster-Constrained Discovery Summary (v1.0)",
  "================================================================================",
  "",
  sprintf("Analysis Time: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("Focus TF Family: %s", PARAMS$focus_tf_family),
  sprintf("Focus Gene: %s", PARAMS$focus_gene),
  "",
  "--- Module A: TF Hub Discovery ---",
  sprintf("  Hub tests: %d",
          ifelse(!is.null(mod_a$hub_table), nrow(mod_a$hub_table), 0)),
  sprintf("  Significant (FDR<0.05): %d",
          ifelse(!is.null(mod_a$hub_table),
                 sum(mod_a$hub_table$fdr < 0.05, na.rm = TRUE), 0)),
  sprintf("  Family enrichment tests: %d",
          ifelse(!is.null(mod_a$family_enrich), nrow(mod_a$family_enrich), 0)),
  "",
  "--- Module B: agriGO Export ---",
  sprintf("  Background mode: %s", PARAMS$background_mode),
  "",
  "--- Module C: Enrichment Matrix ---",
  sprintf("  Total tests: %d",
          ifelse(!is.null(mod_c$enrich_long), nrow(mod_c$enrich_long), 0)),
  sprintf("  Significant (FDR<0.05): %d",
          ifelse(!is.null(mod_c$enrich_long),
                 sum(mod_c$enrich_long$fisher_fdr < 0.05, na.rm = TRUE), 0)),
  "",
  "--- Module D: Focus TF Profiles ---",
  sprintf("  Focus TFs analyzed: %d", length(mod_d$details)),
  sprintf("  Focus gene upstream: %s",
          ifelse(!is.null(mod_d$upstream), "yes", "no")),
  "",
  "--- Module E: MECS Scoring ---",
  sprintf("  Scored edges output: %d", nrow(mod_e$scored_edges))
)

writeLines(summary_lines, file.path(OUT_7B, "Pipeline7b_Summary.txt"))
cat("  Summary saved\n")

# Parameters JSON
params_json <- list(
  pipeline = "7b",
  version = "1.0",
  description = "Cluster-Constrained Discovery + MECS Multi-Evidence Scoring",
  input = list(
    rdata_7a = PARAMS$rdata_7a,
    cluster_pairs = PARAMS$cluster_pairs,
    pcre_tf_map = PARAMS$pcre_tf_map,
    highlight_file = PARAMS$highlight_file
  ),
  focus = list(
    tf_family = PARAMS$focus_tf_family,
    gene = PARAMS$focus_gene,
    gene_is_tf = focus_gene_is_tf,
    gene_is_target = focus_gene_is_target,
    n_focus_tfs = length(focus_tfs)
  ),
  parameters = list(
    enrichment_method = PARAMS$enrichment_method,
    min_edges = PARAMS$min_edges,
    edge_classes = PARAMS$edge_classes_vec,
    background_mode = PARAMS$background_mode,
    mecs_weights = list(
      S1_diffcoexp = PARAMS$w_diffcoexp,
      S2_cluster = PARAMS$w_cluster,
      S3_pcre = PARAMS$w_pcre
    ),
    output_tiers = PARAMS$output_tiers_vec
  ),
  data_dimensions = list(
    n_edges_input = nrow(edge_table),
    n_scored_edges = nrow(mod_e$scored_edges),
    n_wt_clusters = n_wt_clusters,
    n_mut_clusters = n_mut_clusters,
    n_tf_families = length(unique(tf_family_map)),
    cluster_pairs_source = ifelse(!is.null(PARAMS$cluster_pairs),
                                   "external_file", "auto_inferred")
  ),
  analysis_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

writeLines(
  jsonlite::toJSON(params_json, pretty = TRUE, auto_unbox = TRUE),
  file.path(OUT_7B, "Pipeline7b_Parameters.json")
)
cat("  Parameters JSON saved\n")

# RData
scored_edges <- mod_e$scored_edges
edge_table_scored <- mod_e$edge_table_full
hub_table <- mod_a$hub_table
family_enrichment <- mod_a$family_enrich
enrichment_long <- mod_c$enrich_long
enrichment_wide <- mod_c$enrich_wide_list
focus_profiles <- mod_d$profiles
focus_details <- mod_d$details
focus_upstream <- mod_d$upstream
mecs_weights <- mod_e$weights_used
tier_counts <- mod_e$tier_counts

# 保留从 7a 继承的关键对象
save(
  # 7b 核心输出
  scored_edges, edge_table_scored,
  hub_table, family_enrichment,
  enrichment_long, enrichment_wide,
  focus_profiles, focus_details, focus_upstream,
  mecs_weights, tier_counts,
  cluster_pair_info,
  # 从 7a 继承 (7c 需要)
  edge_table, tf_summary,
  tf_in_data, tg_in_data, tg_in_data_full, tf_family_map,
  cluster_annotation, fate_annotation, highlight_map,
  mfuzz_common_genes, non_conserved_genes,
  focus_tfs, focus_gene_is_tf, focus_gene_is_target,
  n_wt_clusters, n_mut_clusters,
  use_pcre,
  # 参数
  PARAMS, PARAMS_7a_inherited,
  file = file.path(OUT_DIR, "rdata", "Pipeline7b_Discovery_Results.RData")
)
cat("  RData saved: Pipeline7b_Discovery_Results.RData\n")

################################################################################
# 完成
################################################################################

cat("\n================================================================================\n")
cat("  Pipeline 7b 完成 (v1.0)\n")
cat("================================================================================\n\n")
cat(sprintf("结束时间: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

cat("输出文件:\n")
cat(sprintf("  %s/Pipeline7b/\n", OUT_DIR))
cat("    ├── Pipeline7b_ModA_TF_Hub_ByCluster.tsv\n")
cat("    ├── Pipeline7b_ModA_TFFamily_Enrichment.tsv\n")
cat("    ├── Pipeline7b_ModA_ClusterPair_Info.tsv\n")
cat("    ├── Pipeline7b_ModC_Enrichment_Long.tsv\n")
cat("    ├── Pipeline7b_ModC_Enrichment_Wide_{EdgeClass}.tsv\n")
if (has_focus) {
  cat("    ├── Pipeline7b_ModD_Focus_TF_Profiles.tsv\n")
  cat("    ├── Pipeline7b_ModD_{TF}_Detail.tsv\n")
}
cat("    ├── Pipeline7b_MECS_Scored_Edges.tsv\n")
cat("    ├── Pipeline7b_MECS_Summary.txt\n")
cat("    ├── Pipeline7b_Summary.txt\n")
cat("    ├── Pipeline7b_Parameters.json\n")
cat("    ├── Pipeline7b_discovery.log\n")
cat("    └── agriGO_lists/\n")
cat(sprintf("  %s/rdata/\n", OUT_DIR))
cat("    └── Pipeline7b_Discovery_Results.RData\n")

cat("\n下游:\n")
cat("  → Pipeline 7c: Visualization (07_run_grn_7c.sh)\n")

sink(type = "message")
sink(type = "output")
close(log_con)

cat("\n日志已保存至: ", log_file, "\n")