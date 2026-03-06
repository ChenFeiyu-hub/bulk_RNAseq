#!/usr/bin/env Rscript
################################################################################
# RNAseq_pipeline4b_PCA.R (v5.0 - DEG Driven)
# 
# 功能: 标准化 PCA 分析 (Global Overview + Recalculated Panel)
# 核心变更:
#   1. [基因选择] 移除 Top N Variance，改为 "Union of DEGs" (所有显著差异基因的并集)
#   2. [视觉] 移除透明度 (Alpha=1.0)，确保颜色完全匹配 Config
#   3. [输入] 解析 Logic 1 或 Master RData 中的比较结果
################################################################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(ggplot2)
  library(ggrepel)
  library(matrixStats)
  library(optparse)
  library(RColorBrewer)
  library(patchwork)
})

# ==============================================================================
# 1. 参数定义
# ==============================================================================
option_list <- list(
  # 输入输出
  make_option(c("-i", "--input"), type = "character", help = "Input RData (must contain vsd AND res_list_time/all_comparisons)"),
  make_option(c("-m", "--metadata"), type = "character", help = "Metadata file"),
  make_option(c("-o", "--outdir"), type = "character", default = "PCA_Output"),
  make_option(c("-p", "--prefix"), type = "character", default = "Project"),
  
  # 配置参数
  make_option(c("--control_cond"), type = "character", default = "WT"),
  
  # 视觉参数 (Key:Value)
  make_option(c("--cond_colors"), type = "character"),
  make_option(c("--cond_shapes"), type = "character"),
  make_option(c("--time_colors"), type = "character"),
  make_option(c("--label_config"), type = "character"),
  
  # 核心统计参数 (用于筛选 DEG)
  make_option(c("--fdr"), type = "double", default = 0.05, help = "FDR cutoff for gene selection"),
  make_option(c("--lfc"), type = "double", default = 1.0, help = "Log2FC cutoff for gene selection"),
  
  # 绘图微调
  make_option(c("--ellipse_conf"), type = "double", default = 0.95)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) stop("Error: --input is required.")

# ==============================================================================
# 2. 辅助函数
# ==============================================================================

parse_config <- function(config_str, as_numeric=FALSE) {
  if (is.null(config_str) || config_str == "") return(NULL)
  clean_str <- gsub('"', '', config_str)
  pairs <- strsplit(clean_str, ",")[[1]]
  keys <- trimws(sapply(strsplit(pairs, ":"), `[`, 1))
  vals <- trimws(sapply(strsplit(pairs, ":"), `[`, 2))
  if (as_numeric) vals <- as.numeric(vals)
  valid <- which(!is.na(keys) & !is.na(vals) & keys != "")
  if(length(valid) == 0) return(NULL)
  setNames(vals[valid], keys[valid])
}

sort_time_points <- function(times) {
  nums <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", times)))
  if (all(is.na(nums))) return(sort(times))
  times[order(nums)]
}

get_colors <- function(items, manual_map, palette_name="Set1") {
  final_cols <- setNames(rep(NA, length(items)), items)
  if (!is.null(manual_map)) {
    matched <- intersect(names(manual_map), items)
    final_cols[matched] <- manual_map[matched]
  }
  missing <- names(final_cols)[is.na(final_cols)]
  if (length(missing) > 0) {
    pool_size <- length(missing)
    if(pool_size <= 8) pool <- brewer.pal(max(3, pool_size), palette_name)
    else pool <- colorRampPalette(brewer.pal(8, palette_name))(pool_size)
    final_cols[missing] <- pool[1:pool_size]
  }
  return(final_cols)
}

get_shapes <- function(items, manual_map) {
  default_shapes <- c(16, 17, 15, 18, 19, 8) 
  final_shapes <- setNames(rep(NA, length(items)), items)
  if (!is.null(manual_map)) {
    matched <- intersect(names(manual_map), items)
    final_shapes[matched] <- manual_map[matched]
  }
  missing <- names(final_shapes)[is.na(final_shapes)]
  if (length(missing) > 0) {
    final_shapes[missing] <- rep(default_shapes, length.out=length(missing))
  }
  return(final_shapes)
}

format_labels <- function(labels, manual_map = NULL) {
  if (is.null(manual_map)) return(labels)
  lapply(labels, function(x) {
    disp <- if(x %in% names(manual_map)) manual_map[[x]] else x
    if (grepl("italic|bold|~|\\^|\\[", disp)) {
      tryCatch(parse(text=disp), error=function(e) disp)
    } else {
      disp
    }
  })
}

# --- 核心函数: 提取 DEG ---
extract_degs <- function(res_obj, fdr_th, lfc_th) {
  df <- as.data.frame(res_obj)
  # 兼容 list 结构
  if("results" %in% names(res_obj)) df <- as.data.frame(res_obj$results)
  
  # 确保有列
  if(!all(c("padj", "log2FoldChange") %in% colnames(df))) return(character(0))
  
  sig <- df %>% 
    filter(!is.na(padj) & padj < fdr_th & abs(log2FoldChange) > lfc_th)
  
  return(rownames(sig))
}

# ==============================================================================
# 3. 数据加载
# ==============================================================================
cat(">>> [Step 1] Loading Data...\n")
load(opt$input)

if (!exists("vsd")) {
  if (exists("dds")) {
    cat("  Applying VST transformation...\n")
    vsd <- vst(dds, blind=FALSE)
  } else {
    stop("Error: RData must contain 'vsd' or 'dds'.")
  }
}

# 读取元数据
meta_df <- read.table(opt$metadata, header=TRUE, sep="\t", stringsAsFactors=FALSE)
if(! "sample" %in% colnames(meta_df)) colnames(meta_df)[1] <- "sample"

# 对齐样本
vst_mat <- assay(vsd)
common <- intersect(colnames(vst_mat), meta_df$sample)
if(length(common) == 0) stop("Error: No common samples!")
vst_mat <- vst_mat[, common]
meta_df <- meta_df[match(common, meta_df$sample), ]

# 自动探测列名
cond_col <- grep("condition|genotype|group", colnames(meta_df), ignore.case=TRUE, value=TRUE)[1]
time_col <- grep("time|hour|stage", colnames(meta_df), ignore.case=TRUE, value=TRUE)[1]
cat(sprintf("  Detected Columns -> Condition: [%s], Time: [%s]\n", cond_col, time_col))

# ==============================================================================
# 4. 基因选择 (基于 DEG)
# ==============================================================================
cat(">>> [Step 2] Selecting Genes based on DEG thresholds...\n")
cat(sprintf("    Thresholds: FDR < %s, |LFC| > %s\n", opt$fdr, opt$lfc))

# 寻找结果列表
res_source <- NULL
if (exists("res_list_time")) res_source <- res_list_time
if (is.null(res_source) && exists("all_comparisons")) res_source <- all_comparisons
if (is.null(res_source) && exists("res_list_genotype")) res_source <- res_list_genotype

if (is.null(res_source)) {
  stop("Error: No differential expression results (res_list_time / all_comparisons) found in RData. Cannot filter by DEGs.")
}

all_degs <- c()
for (n in names(res_source)) {
  degs <- extract_degs(res_source[[n]], opt$fdr, opt$lfc)
  all_degs <- c(all_degs, degs)
}
selected_genes <- unique(all_degs)
# 确保这些基因在矩阵里
selected_genes <- intersect(selected_genes, rownames(vst_mat))

if (length(selected_genes) < 3) {
  stop("Error: Less than 3 DEGs found! Please check your FDR/LFC thresholds or upstream analysis.")
}

cat(sprintf("    Found %d unique DEGs from %d comparisons.\n", length(selected_genes), length(res_source)))

# ==============================================================================
# 5. 视觉映射
# ==============================================================================
# Time Sorting
sorted_times <- sort_time_points(unique(meta_df[[time_col]]))
meta_df[[time_col]] <- factor(meta_df[[time_col]], levels = sorted_times)

# Condition Sorting
raw_conds <- unique(meta_df[[cond_col]])
others <- sort(setdiff(raw_conds, opt$control_cond))
final_conds <- c(opt$control_cond, others)
final_conds <- final_conds[final_conds %in% raw_conds]
meta_df[[cond_col]] <- factor(meta_df[[cond_col]], levels = final_conds)

# Config 解析
manual_cond_cols <- parse_config(opt$cond_colors)
manual_cond_shapes <- parse_config(opt$cond_shapes, as_numeric=TRUE)
manual_time_cols <- parse_config(opt$time_colors)
manual_labels <- parse_config(opt$label_config)

# 生成映射
MAP_COND_COLS <- get_colors(final_conds, manual_cond_cols, "Set1")
MAP_COND_SHAPES <- get_shapes(final_conds, manual_cond_shapes)
MAP_TIME_COLS <- get_colors(sorted_times, manual_time_cols, "Blues")

cat(">>> Visual Mappings Configured.\n")
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$outdir, "Snapshots"), showWarnings = FALSE)

# ==============================================================================
# 6. Global PCA (DEG Set)
# ==============================================================================
cat("\n>>> [Step 3] Plotting Global PCA (Using DEGs)...\n")

pca <- prcomp(t(vst_mat[selected_genes, ]), center=TRUE, scale.=FALSE)
percentVar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)

d <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], sample=rownames(pca$x)) %>%
  merge(meta_df, by="sample")

p_global <- ggplot(d, aes(x=PC1, y=PC2)) +
  stat_ellipse(aes(color=!!sym(time_col)), level=opt$ellipse_conf, 
               type="norm", linetype="dashed", alpha=0.5, show.legend=FALSE) +
  
  # Global Point: Color=Time, Shape=Condition (无透明度)
  geom_point(aes(color=!!sym(time_col), shape=!!sym(cond_col)), 
             size=5, stroke=0.8) +
  
  geom_text_repel(aes(label=sample), size=2.5, max.overlaps=10, color="grey40") +
  
  scale_color_manual(values=MAP_TIME_COLS, 
                     labels=format_labels(names(MAP_TIME_COLS), manual_labels), name="Time") +
  scale_shape_manual(values=MAP_COND_SHAPES, 
                     labels=format_labels(names(MAP_COND_SHAPES), manual_labels), name="Condition") +
  
  labs(title = paste0("Global PCA: ", opt$prefix),
       subtitle = sprintf("Based on %d DEGs (Union) | Color=Time, Shape=Condition", length(selected_genes)),
       x = paste0("PC1: ", percentVar[1], "% variance"),
       y = paste0("PC2: ", percentVar[2], "% variance")) +
  theme_bw(base_size=14)

ggsave(file.path(opt$outdir, paste0(opt$prefix, "_PCA_Global_Overview.pdf")), p_global, width=10, height=8)

# ==============================================================================
# 7. Recalculated Panel (Local PCA on DEG Subset)
# ==============================================================================
cat("\n>>> [Step 4] Plotting Recalculated Snapshots...\n")

plot_list <- list()

for (tp in sorted_times) {
  sub_meta <- meta_df[meta_df[[time_col]] == tp, ]
  if(nrow(sub_meta) < 3) {
    cat(sprintf("  [Skip] %s: too few samples\n", tp)); next
  }
  sub_mat <- vst_mat[selected_genes, sub_meta$sample] # 依然只使用 DEG 集合
  

  # 移除在局部没有变异的基因 (防止 prcomp 报错)
  sub_mat <- sub_mat[rowVars(sub_mat) > 0, ]
  
  sub_pca <- prcomp(t(sub_mat), center=TRUE, scale.=FALSE)
  sub_pct <- round(100 * sub_pca$sdev^2 / sum(sub_pca$sdev^2), 1)
  
  sub_d <- data.frame(sample=rownames(sub_pca$x), PC1=sub_pca$x[,1], PC2=sub_pca$x[,2]) %>%
    merge(sub_meta, by="sample")
  
  p <- ggplot(sub_d, aes(x=PC1, y=PC2)) +
    # 移除 Alpha，使用纯色
    geom_point(aes(color=!!sym(cond_col)), shape=16, size=5) +
    
    geom_text_repel(aes(label=sample), size=3) +
    scale_color_manual(values=MAP_COND_COLS, 
                       labels=format_labels(names(MAP_COND_COLS), manual_labels)) +
    
    labs(title = format_labels(tp, manual_labels)[[1]], 
         subtitle = sprintf("PC1: %.1f%% | PC2: %.1f%%", sub_pct[1], sub_pct[2]),
         x=NULL, y=NULL) +
    theme_bw(base_size=12) +
    theme(legend.position="none", 
          plot.title=element_text(hjust=0.5, face="bold"),
          plot.subtitle=element_text(hjust=0.5, size=9, color="grey30"))
  
  ggsave(file.path(opt$outdir, "Snapshots", paste0(opt$prefix, "_PCA_", tp, ".pdf")), p, width=5, height=5)
  plot_list[[tp]] <- p
}

if (length(plot_list) > 0) {
  # 共享图例
  legend_plot <- ggplot(meta_df, aes(x=1, y=1, color=!!sym(cond_col))) +
    geom_point(size=5, shape=16) + # 纯色图例
    scale_color_manual(values=MAP_COND_COLS, 
                       labels=format_labels(names(MAP_COND_COLS), manual_labels), 
                       name="Condition") +
    theme_void() + theme(legend.position="bottom", legend.text=element_text(size=12))
  
  shared_legend <- cowplot::get_legend(legend_plot)
  
  n_col <- min(3, length(plot_list))
  n_row <- ceiling(length(plot_list) / n_col)
  
  final_panel <- wrap_plots(plot_list, ncol=n_col) + 
    plot_layout(guides = "collect") & theme(legend.position = "bottom")
  
  final_panel <- final_panel + plot_annotation(
    title = "PCA Snapshots (Recalculated on DEGs)",
    subtitle = sprintf("PCA calculated per timepoint using %d global DEGs", length(selected_genes)),
    theme = theme(plot.title = element_text(size=16, face="bold", hjust=0.5))
  )
  
  out_panel <- file.path(opt$outdir, paste0(opt$prefix, "_PCA_Panel_Recalculated.pdf"))
  ggsave(out_panel, final_panel, width=4*n_col, height=4*n_row+1)
  cat(sprintf("  ✓ Panel Saved: %s\n", basename(out_panel)))
}

cat("\nDone.\n")