#!/usr/bin/env Rscript
################################################################################
# generate_highlight_from_clusters.R
#
# 版本: v4.6
#
# 功能: 当用户没有预定义 highlight 列表时，自动从 Mfuzz cluster 生成 highlight:
#   来源 A: 用户指定的 Mfuzz cluster 中的基因 (pCRE 驱动)
#   来源 B: MECS 质量过滤 — 保留在 scored_edges 中有足够证据支撑的 Source A 基因
#
# 设计原则 (v4.6):
#   - 输出中不包含任何 TF。TF 在 7c 的 anchor-edge 机制中以 "Other TF"
#     身份自然进入网络，但不会辐射出 cluster 无关的靶标。
#   - 如果将 TF 放入 anchor 集，Step 1c 会选取该 TF 的全部显著边，
#     引入大量与指定 cluster 无关的靶标，导致网络膨胀。
#   - v4.5 的 Source B (hub target 补充) 和 Source C (全局 TF 补充) 已移除，
#     因为它们在合并模式下导致 secondary anchor 过度膨胀 (~295 genes)，
#     偏离了 "cluster 补充视角" 的设计初衷。
#
# 输出: 标准 highlight 格式 (TSV: Symbol \t Gene)
#
# v4.6 变更 (相对 v4.5):
#   - 移除 --mode 参数 (不再区分 full/cluster_filter)
#   - 移除 --hub_top_n 参数 (不再提取 hub target)
#   - 新增 --min_tier / --min_edge_count 参数 (MECS 质量门槛)
#   - Source B 统一为 MECS 质量过滤
#   - Source C 永久跳过 (TF 不进入 auto-highlight)
################################################################################

suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
})

option_list <- list(
  # === 来源 A: Cluster-based ===
  make_option("--wt_membership", type = "character", default = NULL,
              help = "WT Membership Matrix 文件路径"),
  make_option("--mut_membership", type = "character", default = NULL,
              help = "MUT Membership Matrix 文件路径"),
  make_option("--clusters", type = "character", default = NULL,
              help = "逗号分隔的 cluster ID 列表 (如 'WT_C8,WT_C9,MUT_C3')"),
  make_option("--membership_threshold", type = "double", default = 0.5,
              help = "Membership 阈值 [default: %default]"),
  make_option("--max_genes_per_cluster", type = "integer", default = 100,
              help = "每个 cluster 最多提取的基因数 [default: %default]"),

  # === 来源 B: MECS 质量过滤 ===
  make_option("--rdata_7b", type = "character", default = NULL,
              help = "Pipeline7b RData 路径 (可选，用于 MECS 质量过滤)"),
  make_option("--tier_filter", type = "character", default = "1,2",
              help = "纳入的 Tier 等级 [default: %default]"),
  make_option("--min_tier", type = "integer", default = 3,
              help = "Source A 基因须参与至少1条 Tier <= N 的边 [default: %default]"),
  make_option("--min_edge_count", type = "integer", default = 1,
              help = "Source A 基因须参与的最小边数 [default: %default]"),

  # === 注释 ===
  make_option("--gff_file", type = "character", default = NULL,
              help = "GFF3 文件路径 (用于提取 Gene → Symbol 映射, 可选)"),
  make_option("--gene_id_strip", type = "character", default = "",
              help = "基因ID清理正则"),

  # === 输出 ===
  make_option("--output", type = "character", default = "auto_highlight.txt",
              help = "输出文件路径 [default: %default]"),
  make_option("--quiet", action = "store_true", default = FALSE)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (!opt$quiet) {
  cat("=============================================\n")
  cat(" Auto-Generate Highlight Gene List (v4.6)\n")
  cat("=============================================\n\n")
}

# ================================================================
# 辅助: 基因ID标准化
# ================================================================
standardize_ids <- function(ids, pattern = "") {
  if (is.null(pattern) || pattern == "") return(ids)
  gsub(pattern, "", ids)
}

# ================================================================
# 辅助: 从 GFF 提取 Symbol 映射
# ================================================================
build_symbol_map <- function(gff_path, strip_pattern = "") {
  if (is.null(gff_path) || !file.exists(gff_path)) return(NULL)

  cat("  Loading GFF for symbol mapping...\n")
  tryCatch({
    gff_lines <- readLines(gff_path)
    gene_lines <- gff_lines[grepl("\tgene\t", gff_lines)]

    # 提取 ID 和 Name 属性
    extract_attr <- function(line, key) {
      m <- regmatches(line, regexec(paste0(key, "=([^;]+)"), line))
      if (length(m[[1]]) >= 2) return(m[[1]][2])
      return(NA_character_)
    }

    ids <- sapply(gene_lines, function(l) extract_attr(l, "ID"), USE.NAMES = FALSE)
    names_attr <- sapply(gene_lines, function(l) extract_attr(l, "Name"), USE.NAMES = FALSE)

    ids <- standardize_ids(ids, strip_pattern)
    names_attr <- standardize_ids(names_attr, strip_pattern)  # ★ fix: Symbol 也需要 strip
    symbol_map <- setNames(names_attr, ids)
    symbol_map <- symbol_map[!is.na(symbol_map) & !is.na(names(symbol_map))]

    cat(sprintf("  Loaded %d gene symbols from GFF\n", length(symbol_map)))
    return(symbol_map)
  }, error = function(e) {
    cat(sprintf("  Warning: GFF parsing failed: %s\n", e$message))
    return(NULL)
  })
}

# ================================================================
# 来源 A: 从 Mfuzz Cluster 提取基因
# ================================================================
genes_from_clusters <- character(0)

if (!is.null(opt$clusters) && opt$clusters != "" && opt$clusters != "none") {

  cluster_ids <- trimws(strsplit(opt$clusters, ",")[[1]])
  cluster_ids <- cluster_ids[cluster_ids != ""]
  cat(sprintf("Source A: Extracting genes from %d clusters: %s\n",
              length(cluster_ids), paste(cluster_ids, collapse = ", ")))

  # 分离 WT 和 MUT clusters
  wt_clusters  <- cluster_ids[grepl("^WT_", cluster_ids)]
  mut_clusters <- cluster_ids[grepl("^MUT_", cluster_ids)]

  # 也支持不带前缀的 cluster (如 "C8") — 此时同时在 WT 和 MUT 中查找
  bare_clusters <- cluster_ids[!grepl("^(WT_|MUT_)", cluster_ids)]

  process_membership <- function(mem_file, prefix, target_clusters, bare_cls) {
    if (is.null(mem_file) || !file.exists(mem_file)) return(character(0))

    mem_df <- read.table(mem_file, header = TRUE, sep = "\t",
                         stringsAsFactors = FALSE, check.names = FALSE)

    gene_col <- intersect(c("gene", "GeneID"), names(mem_df))[1]
    if (is.na(gene_col)) gene_col <- names(mem_df)[1]

    genes_raw <- standardize_ids(mem_df[[gene_col]], opt$gene_id_strip)
    cluster_cols <- setdiff(names(mem_df), gene_col)

    collected <- character(0)

    # 安全地构建待查找的 cluster 列表，过滤掉可能潜伏的空字符串
    clusters_to_check <- target_clusters
    if (length(bare_cls) > 0) {
      valid_bare_cls <- bare_cls[bare_cls != ""]
      if (length(valid_bare_cls) > 0) {
        clusters_to_check <- c(clusters_to_check, paste0(prefix, "_", valid_bare_cls))
      }
    }

    for (cl in clusters_to_check) {
      # 拦截异常的裸前缀
      if (cl == "WT_" || cl == "MUT_") {
        next
      }

      # 尝试精确匹配
      matched_col <- NULL
      if (cl %in% cluster_cols) {
        matched_col <- cl
      } else {
        # 模糊匹配: 去掉前缀后比对
        cl_num <- gsub("^(WT_|MUT_)", "", cl)

        # 防御 cl_num 为空字符串
        if (cl_num == "") {
          next
        }

        candidates <- cluster_cols[grepl(cl_num, cluster_cols)]
        if (length(candidates) > 0) matched_col <- candidates[1]
      }

      if (!is.null(matched_col)) {
        membership_vals <- as.numeric(mem_df[[matched_col]])
        high_mem_idx <- which(membership_vals >= opt$membership_threshold)

        # 按 membership 排序取 top N
        if (length(high_mem_idx) > opt$max_genes_per_cluster) {
          high_mem_idx <- high_mem_idx[order(membership_vals[high_mem_idx],
                                              decreasing = TRUE)]
          high_mem_idx <- high_mem_idx[1:opt$max_genes_per_cluster]
        }

        cl_genes <- genes_raw[high_mem_idx]
        collected <- c(collected, cl_genes)
        cat(sprintf("  %s [%s]: %d genes (threshold=%.2f)\n",
                    cl, matched_col, length(cl_genes), opt$membership_threshold))
      } else {
        cat(sprintf("  %s: column not found, skipping\n", cl))
      }
    }
    return(collected)
  }

  # 处理 WT
  if (length(wt_clusters) > 0 || length(bare_clusters) > 0) {
    wt_genes <- process_membership(opt$wt_membership, "WT",
                                    wt_clusters, bare_clusters)
    genes_from_clusters <- c(genes_from_clusters, wt_genes)
  }

  # 处理 MUT
  if (length(mut_clusters) > 0 || length(bare_clusters) > 0) {
    mut_genes <- process_membership(opt$mut_membership, "MUT",
                                     mut_clusters, bare_clusters)
    genes_from_clusters <- c(genes_from_clusters, mut_genes)
  }

  genes_from_clusters <- unique(genes_from_clusters)
  cat(sprintf("  Source A total: %d unique genes\n\n", length(genes_from_clusters)))
}

# ================================================================
# 来源 B: MECS 质量过滤 (v4.6)
# ================================================================
# 检查每个 Source A 基因在 7b scored_edges 中的参与情况。
# 基因必须参与至少 min_edge_count 条 Tier <= min_tier 的边才保留。
# 未提供 7b 数据时跳过过滤，保留全部 Source A 基因。
#
# 注意: 不再提取 hub target 或全局 TF (v4.5 的 Source B/C 已废弃)。
# TF 在 7c 中通过 anchor-edge 机制自然进入网络。
# ================================================================

genes_filtered <- genes_from_clusters  # 默认: 无过滤

if (!is.null(opt$rdata_7b) && file.exists(opt$rdata_7b)) {

  cat("Source B: MECS quality filtering...\n")
  cat(sprintf("  Criteria: gene must appear in >= %d edges of Tier <= %d\n",
              opt$min_edge_count, opt$min_tier))

  tryCatch({
    env_7b <- new.env()
    load(opt$rdata_7b, envir = env_7b)

    # ★ Bug fix v4.6: 优先按名称查找 scored_edges (含 Tier 列)
    #   旧逻辑用 ls() 遍历取第一个含 TF_Gene 的 data.frame，
    #   但 ls() 按字母排序，edge_table (e) 排在 scored_edges (s) 前面，
    #   导致加载了无 Tier 列的原始全边表 (640k rows)。
    scored_edges <- NULL
    
    # 优先: 直接按名称获取 scored_edges
    for (preferred in c("scored_edges", "edge_table_scored")) {
      if (exists(preferred, envir = env_7b)) {
        obj <- get(preferred, envir = env_7b)
        if (is.data.frame(obj) && "TF_Gene" %in% names(obj) &&
            "Target_Gene" %in% names(obj) && "Tier" %in% names(obj)) {
          scored_edges <- obj
          cat(sprintf("  Found scored edges: '%s' (%d rows, has Tier column)\n",
                      preferred, nrow(obj)))
          break
        }
      }
    }
    
    # 回退: 通用搜索 (仅当按名称找不到时)
    if (is.null(scored_edges)) {
      cat("  Preferred objects not found, searching generically...\n")
      for (obj_name in ls(env_7b)) {
        obj <- get(obj_name, envir = env_7b)
        if (is.data.frame(obj) && "TF_Gene" %in% names(obj) &&
            "Target_Gene" %in% names(obj)) {
          scored_edges <- obj
          cat(sprintf("  Found scored edges (fallback): '%s' (%d rows)\n",
                      obj_name, nrow(obj)))
          break
        }
      }
    }

    if (!is.null(scored_edges)) {
      n_before <- length(genes_from_clusters)

      if (n_before > 0) {
        max_tier <- opt$min_tier
        min_n    <- opt$min_edge_count

        # 筛选合格边 (Tier <= 阈值)
        if ("Tier" %in% names(scored_edges)) {
          qualified_edges <- scored_edges %>%
            filter(Tier <= max_tier)
        } else {
          # 无 Tier 列时保留所有边
          qualified_edges <- scored_edges
          cat("  Warning: No 'Tier' column in scored_edges, using all edges\n")
        }

        cat(sprintf("  Qualified edges (Tier <= %d): %d\n",
                    max_tier, nrow(qualified_edges)))

        # ★ Bug fix v4.6: 对 edge table 中的基因 ID 应用相同的标准化
        #   Source A 基因经过 gene_id_strip 处理 (如去掉 _4532.v6.1 后缀)，
        #   但 scored_edges 中的 TF_Gene/Target_Gene 保留原始格式，
        #   导致比较时两边 ID 永远不一致。
        edge_tf_ids     <- standardize_ids(qualified_edges$TF_Gene, opt$gene_id_strip)
        edge_target_ids <- standardize_ids(qualified_edges$Target_Gene, opt$gene_id_strip)

        cat(sprintf("  ID format check: Source A[1]='%s' | Edge TF[1]='%s' | Edge Target[1]='%s'\n",
                    genes_from_clusters[1], edge_tf_ids[1], edge_target_ids[1]))

        # 统计每个 Source A 基因参与的合格边数
        # (作为 TF 调控端 或 作为 Target 被调控端)
        gene_edge_counts <- tibble(Gene = genes_from_clusters) %>%
          rowwise() %>%
          mutate(
            n_as_target = sum(edge_target_ids == Gene),
            n_as_tf     = sum(edge_tf_ids == Gene),
            n_total     = n_as_target + n_as_tf
          ) %>%
          ungroup()

        # 过滤
        genes_filtered <- gene_edge_counts %>%
          filter(n_total >= min_n) %>%
          pull(Gene)

        n_after <- length(genes_filtered)
        n_removed <- n_before - n_after

        cat(sprintf("  Source A genes: %d -> %d (removed %d with insufficient edges)\n",
                    n_before, n_after, n_removed))

        # 被移除基因的统计
        if (n_removed > 0) {
          cat(sprintf("  Removed genes (no Tier<=%d edges): %d\n",
                      max_tier, n_removed))
        }

        # 通过基因的 edge 分布
        if (n_after > 0) {
          passed <- gene_edge_counts %>% filter(n_total >= min_n)
          cat(sprintf("  Passed genes edge count: median=%g, range=[%d, %d]\n",
                      median(passed$n_total),
                      min(passed$n_total),
                      max(passed$n_total)))
        }
      } else {
        cat("  No Source A genes to filter.\n")
      }
    } else {
      cat("  Warning: No scored_edges found in 7b RData. Skipping filter.\n")
    }

    cat("Source C: Skipped (TFs enter 7c network via anchor-edge mechanism)\n\n")

  }, error = function(e) {
    cat(sprintf("  Warning: Failed to load 7b data: %s\n", e$message))
    cat("  All Source A genes retained without MECS filtering.\n\n")
  })
} else {
  cat("Source B: 7b RData not provided, skipping MECS quality filter.\n")
  cat("  All Source A genes retained without filtering.\n\n")
}

# ================================================================
# 合并 (v4.6: 仅输出过滤后的 Source A 基因)
# ================================================================
all_highlight_genes <- unique(genes_filtered)

if (length(all_highlight_genes) == 0) {
  cat("WARNING: No genes collected from any source!\n")
  cat("  Check --clusters and membership files.\n")
  cat("  If using MECS filter, consider relaxing --min_tier or --min_edge_count.\n")
  # 写入空文件避免下游报错
  writeLines(c("Symbol\tGene"), opt$output)
  quit(status = 1)
}

cat(sprintf("Combined: %d genes (Source A filtered)\n",
            length(all_highlight_genes)))
if (length(genes_from_clusters) > 0) {
  cat(sprintf("  Source A (raw):      %d\n", length(genes_from_clusters)))
  cat(sprintf("  After MECS filter:   %d\n", length(all_highlight_genes)))
  cat(sprintf("  Removed by filter:   %d\n",
              length(genes_from_clusters) - length(all_highlight_genes)))
}

# ================================================================
# Symbol 映射
# ================================================================
symbol_map <- build_symbol_map(opt$gff_file, opt$gene_id_strip)

if (!is.null(symbol_map)) {
  symbols <- ifelse(all_highlight_genes %in% names(symbol_map),
                    symbol_map[all_highlight_genes],
                    all_highlight_genes)
} else {
  # 没有 GFF: 用基因 ID 本身作为 Symbol (去掉常见后缀)
  symbols <- gsub("\\.[0-9]+$", "", all_highlight_genes)
  symbols <- gsub("_.*$", "", symbols)
}

# ================================================================
# 输出
# ================================================================
output_df <- data.frame(
  Symbol = symbols,
  Gene   = all_highlight_genes,
  stringsAsFactors = FALSE
) %>%
  distinct(Gene, .keep_all = TRUE) %>%
  arrange(Symbol)

# v4.6: Source 列统一为 cluster_filtered
output_df$Source <- "cluster_filtered"

write.table(output_df[, c("Symbol", "Gene")], opt$output,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 同时输出一份带来源标注的详细版
detail_path <- sub("\\.txt$", "_detailed.txt", opt$output)
write.table(output_df, detail_path,
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\nOutput: %s (%d genes)\n", opt$output, nrow(output_df)))
cat(sprintf("  Detailed: %s\n", detail_path))
cat("Done.\n")