#!/usr/bin/env Rscript
################################################################################
# generate_highlight_from_clusters.R
# 
# 功能: 当用户没有预定义 highlight 列表时，自动从以下来源生成:
#   来源 A: 用户指定的 Mfuzz cluster 中的基因 (pCRE 驱动)
#   来源 B: Pipeline 7b 中的 hub target 基因 (高 in-degree)
#   来源 C: Pipeline 7b 中所有 Tier 1-2 的 TF
#
# 输出: 标准 highlight 格式 (TSV: Symbol \t Gene)
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
  
  # === 来源 B/C: 7b hub genes ===
  make_option("--rdata_7b", type = "character", default = NULL,
              help = "Pipeline7b RData 路径 (可选，用于提取 hub genes)"),
  make_option("--hub_top_n", type = "integer", default = 30,
              help = "提取 top N hub target 基因 [default: %default]"),
  make_option("--tier_filter", type = "character", default = "1,2",
              help = "纳入的 Tier 等级 [default: %default]"),
  
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
  cat("=====================================\n")
  cat(" Auto-Generate Highlight Gene List\n")
  cat("=====================================\n\n")
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
    
    # 修复点 1：安全地构建待查找的 cluster 列表，过滤掉可能潜伏的空字符串
    clusters_to_check <- target_clusters
    if (length(bare_cls) > 0) {
      valid_bare_cls <- bare_cls[bare_cls != ""]
      if (length(valid_bare_cls) > 0) {
        clusters_to_check <- c(clusters_to_check, paste0(prefix, "_", valid_bare_cls))
      }
    }
    
    for (cl in clusters_to_check) {
      # 修复点 2：拦截异常的裸前缀，防止后续逻辑暴走
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
        
        # 修复点 3：绝对防御，防止 cl_num 为空字符串时 grepl 匹配所有列
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
# 来源 B/C: 从 Pipeline 7b 提取 hub genes 和 TFs
# ================================================================
genes_from_7b <- character(0)

if (!is.null(opt$rdata_7b) && file.exists(opt$rdata_7b)) {
  cat("Source B/C: Extracting hub genes from Pipeline 7b...\n")
  
  tryCatch({
    env_7b <- new.env()
    load(opt$rdata_7b, envir = env_7b)
    
    # 查找 scored_edges
    scored_edges <- NULL
    for (obj_name in ls(env_7b)) {
      obj <- get(obj_name, envir = env_7b)
      if (is.data.frame(obj) && "TF_Gene" %in% names(obj) && 
          "Target_Gene" %in% names(obj)) {
        scored_edges <- obj
        cat(sprintf("  Found scored edges: '%s' (%d rows)\n", obj_name, nrow(obj)))
        break
      }
    }
    
    if (!is.null(scored_edges)) {
      tier_vals <- as.integer(trimws(strsplit(opt$tier_filter, ",")[[1]]))
      
      if ("Tier" %in% names(scored_edges)) {
        tiered <- scored_edges %>% filter(Tier %in% tier_vals)
      } else {
        tiered <- scored_edges
      }
      
      # B: Hub targets (高 in-degree)
      hub_targets <- tiered %>%
        count(Target_Gene, name = "in_degree") %>%
        arrange(desc(in_degree)) %>%
        head(opt$hub_top_n) %>%
        pull(Target_Gene)
      
      # C: 所有涉及的 TF
      hub_tfs <- tiered %>%
        distinct(TF_Gene) %>%
        pull(TF_Gene)
      
      genes_from_7b <- unique(c(hub_targets, hub_tfs))
      cat(sprintf("  Source B: %d hub targets (top %d by in-degree)\n",
                  length(hub_targets), opt$hub_top_n))
      cat(sprintf("  Source C: %d TFs from Tier %s edges\n",
                  length(hub_tfs), paste(tier_vals, collapse = ",")))
      cat(sprintf("  Source B+C total: %d unique genes\n\n", length(genes_from_7b)))
    }
  }, error = function(e) {
    cat(sprintf("  Warning: Failed to load 7b data: %s\n", e$message))
  })
} else {
  cat("Source B/C: 7b RData not provided or not found, skipping.\n\n")
}

# ================================================================
# 合并所有来源
# ================================================================
all_highlight_genes <- unique(c(genes_from_clusters, genes_from_7b))

if (length(all_highlight_genes) == 0) {
  cat("⚠ WARNING: No genes collected from any source!\n")
  cat("  Check --clusters and --rdata_7b parameters.\n")
  # 写入空文件避免下游报错
  writeLines(c("Symbol\tGene"), opt$output)
  quit(status = 1)
}

cat(sprintf("Combined: %d unique highlight genes\n", length(all_highlight_genes)))
cat(sprintf("  From clusters: %d\n", length(genes_from_clusters)))
cat(sprintf("  From 7b hubs:  %d\n", length(genes_from_7b)))
cat(sprintf("  Overlap:       %d\n",
            length(intersect(genes_from_clusters, genes_from_7b))))

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

# 添加来源标注列 (可选，方便调试)
output_df$Source <- ifelse(
  output_df$Gene %in% genes_from_clusters & output_df$Gene %in% genes_from_7b,
  "cluster+7b",
  ifelse(output_df$Gene %in% genes_from_clusters, "cluster", "7b_hub")
)

write.table(output_df[, c("Symbol", "Gene")], opt$output,
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 同时输出一份带来源标注的详细版
detail_path <- sub("\\.txt$", "_detailed.txt", opt$output)
write.table(output_df, detail_path,
            sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("\n✔ Output: %s (%d genes)\n", opt$output, nrow(output_df)))
cat(sprintf("  Detailed: %s\n", detail_path))
cat("Done.\n")