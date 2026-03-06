#!/usr/bin/env Rscript
################################################################################
# RNAseq_pipeline4d_Venn.R (v3.1 - Color-bindind bug fix)
# 
# 功能: 
#   1. 可视化: Logic 1 两两比较韦恩图面板 (剔除 Control Time)
#   2. 数据资产: 结构化的基因列表输出 (Time_Comparison_Direction 目录隔离)
#
# v3.1 修复:
#   - [BUG FIX] VennDiagram 在 scaled/euler.d 模式下会根据集合大小
#     自动交换左右圆圈位置，但不会同步交换 fill/col 颜色向量，
#     导致颜色与基因型错位。
#   - 修复方案: 在调用 venn.diagram 前，对输入 list 和颜色向量做
#     动态预排序，始终将较大集合置于第一位，使 VennDiagram 内部
#     的交换逻辑不再触发，从而保证颜色-数据严格绑定。
#
# 视觉设计:
#   - 极简标题，无遮挡
#   - 透明填充 + 实心边框
#   - 彩色时间轴 + 虚线网格
################################################################################

suppressPackageStartupMessages({
  library(VennDiagram)
  library(grid)
  library(gridExtra)
  library(tidyverse)
  library(optparse)
  library(RColorBrewer)
})

futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

# ==============================================================================
# 1. 参数定义
# ==============================================================================
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input Logic 1 RData"),
  make_option(c("-o", "--outdir"), type = "character", default = "Venn_Output"),
  make_option(c("-p", "--prefix"), type = "character", default = "Project"),
  
  # 关键配置
  make_option(c("--control_cond"), type = "character", default = "WT"),
  make_option(c("--control_time"), type = "character", default = "00h"),
  make_option(c("--condition_order"), type = "character", default = NULL),
  make_option(c("--time_order"), type = "character", default = NULL),
  
  # 颜色配置
  make_option(c("--cond_colors"), type = "character"),
  make_option(c("--time_colors"), type = "character"),
  
  # 阈值
  make_option(c("--fdr"), type = "double", default = 0.05),
  make_option(c("--lfc"), type = "double", default = 1.0)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) stop("Error: --input is required.")

# ==============================================================================
# 2. 辅助函数
# ==============================================================================

parse_config <- function(config_str) {
  if (is.null(config_str) || config_str == "") return(NULL)
  clean_str <- gsub('"', '', config_str)
  pairs <- strsplit(clean_str, ",")[[1]]
  keys <- trimws(sapply(strsplit(pairs, ":"), `[`, 1))
  vals <- trimws(sapply(strsplit(pairs, ":"), `[`, 2))
  valid <- which(!is.na(keys) & !is.na(vals) & keys != "")
  if(length(valid) == 0) return(NULL)
  setNames(vals[valid], keys[valid])
}

get_deg_ids <- function(res, fdr_th, lfc_th, direction) {
  df <- as.data.frame(res)
  if("results" %in% names(res)) df <- as.data.frame(res$results)
  
  if(!all(c("padj", "log2FoldChange") %in% colnames(df))) return(character(0))
  
  if (direction == "Up") {
    idx <- which(df$padj < fdr_th & df$log2FoldChange > lfc_th)
  } else {
    idx <- which(df$padj < fdr_th & df$log2FoldChange < -lfc_th)
  }
  
  if("gene" %in% colnames(df)) return(unique(as.character(df$gene[idx])))
  return(unique(rownames(df)[idx]))
}

get_color <- function(key, map, default="grey50") {
  if(!is.null(map) && key %in% names(map)) return(map[[key]])
  return(default)
}

add_alpha <- function(col, alpha=0.5) {
  rgb_vals <- col2rgb(col)
  rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], max = 255, alpha = alpha * 255)
}

# ==============================================================================
# 3. 核心功能: 导出基因列表 (修复目录结构)
# ==============================================================================
export_venn_lists <- function(set_A, set_B, name_A, name_B, timepoint, direction, root_dir) {
  intersect_genes <- intersect(set_A, set_B)
  only_A <- setdiff(set_A, set_B)
  only_B <- setdiff(set_B, set_A)
  
  folder_name <- paste0(timepoint, "_", name_B, "_vs_", name_A, "_", direction)
  sub_dir <- file.path(root_dir, "Venn_GeneLists", folder_name)
  
  dir.create(sub_dir, recursive = TRUE, showWarnings = FALSE)
  
  file_intersect <- file.path(sub_dir, paste0("Intersection_", name_A, "_vs_", name_B, ".txt"))
  file_only_A    <- file.path(sub_dir, paste0("Exclusive_", name_A, ".txt"))
  file_only_B    <- file.path(sub_dir, paste0("Exclusive_", name_B, ".txt"))
  
  if(length(intersect_genes) > 0) writeLines(intersect_genes, file_intersect)
  if(length(only_A) > 0)          writeLines(only_A, file_only_A)
  if(length(only_B) > 0)          writeLines(only_B, file_only_B)
  
  return(list(
    n_intersect = length(intersect_genes),
    n_A = length(only_A),
    n_B = length(only_B)
  ))
}

# ==============================================================================
# 4. 绘图对象生成 [v3.1 核心修复]
# ==============================================================================
create_venn_grob <- function(set_list, set_colors) {
  if (length(set_list) != 2) return(textGrob("Error", gp=gpar(col="red")))
  
  # ==========================================================================
  # [v3.1 BUG FIX] 动态预排序: 始终将较大集合置于 list 第一位
  #
  # 原理: VennDiagram::draw.pairwise.venn 在 scaled/euler.d 模式下,
  #        内部会检查 area1 vs area2, 若 area1 < area2 则交换两个圆的
  #        几何位置, 但 fill/col 等视觉参数按原始位置索引, 不跟随交换.
  #        通过预排序保证 length(set_list[[1]]) >= length(set_list[[2]]),
  #        VennDiagram 的交换分支永远不会触发, 颜色-数据保持严格绑定.
  # ==========================================================================
  sizes <- sapply(set_list, length)
  if (sizes[1] < sizes[2]) {
    original_names <- names(set_list)
    set_list   <- set_list[c(2, 1)]
    set_colors <- set_colors[c(2, 1)]
    cat(sprintf("    [Pre-sort] Swapped order: '%s'(%d) now before '%s'(%d)\n",
                original_names[2], sizes[2], original_names[1], sizes[1]))
  }
  
  fill_cols   <- sapply(set_colors, add_alpha, alpha=0.5)
  border_cols <- set_colors
  
  tryCatch({
    venn_obj <- venn.diagram(
      x = set_list,
      filename = NULL,
      fill = fill_cols,
      col = border_cols,
      lwd = 2,                 
      scaled = TRUE,           
      euler.d = TRUE,
      category.names = c("", ""), 
      cex = 1.0,
      fontface = "bold",
      fontfamily = "sans",
      margin = 0.05,
      alpha = 0.5
    )
    
    # -----------------------------------------------------------------------
    # [v3.1 增强] 添加图例, 使颜色-基因型对应关系一目了然
    # 无论内部顺序如何, 图例始终以 name:color 配对展示, 不依赖位置
    # -----------------------------------------------------------------------
    legend_labels <- names(set_list)
    legend_colors <- unname(set_colors)
    
    legend_grob <- gTree(children = gList(
      # 第一个图例项 (左侧)
      circleGrob(x = unit(0.15, "npc"), y = unit(0.5, "npc"), r = unit(0.03, "npc"),
                 gp = gpar(fill = add_alpha(legend_colors[1], 0.5),
                           col = legend_colors[1], lwd = 2)),
      textGrob(sprintf("%s (%d)", legend_labels[1], sizes[max(which(names(set_list) == legend_labels[1]))]),
               x = unit(0.22, "npc"), y = unit(0.5, "npc"),
               just = "left", gp = gpar(fontsize = 7, col = legend_colors[1])),
      # 第二个图例项 (右侧)
      circleGrob(x = unit(0.55, "npc"), y = unit(0.5, "npc"), r = unit(0.03, "npc"),
                 gp = gpar(fill = add_alpha(legend_colors[2], 0.5),
                           col = legend_colors[2], lwd = 2)),
      textGrob(sprintf("%s (%d)", legend_labels[2], sizes[min(which(names(set_list) == legend_labels[2]))]),
               x = unit(0.62, "npc"), y = unit(0.5, "npc"),
               just = "left", gp = gpar(fontsize = 7, col = legend_colors[2]))
    ))
    
    # 组合: 虚线网格 + 韦恩图 + 底部图例
    gTree(children = gList(
      rectGrob(gp = gpar(col = "grey85", lty = "dashed", fill = NA)), 
      venn_obj,
      # 图例放在底部 5% 的区域
      editGrob(legend_grob, vp = viewport(y = 0.02, height = 0.06, just = "bottom"))
    ))
    
  }, error = function(e) {
    textGrob("No Overlap", gp=gpar(fontsize=8, col="grey60"))
  })
}

create_time_header <- function(label, bg_color) {
  gTree(children = gList(
    rectGrob(gp = gpar(fill = bg_color, col = NA)), 
    textGrob(label, gp = gpar(col = "white", fontface = "bold", fontsize=10))
  ))
}

# ==============================================================================
# 5. 主流程
# ==============================================================================
cat(">>> Loading Data...\n")
load(opt$input)

if (!exists("res_list_time")) {
  if (exists("all_comparisons")) {
    res_list_time <- all_comparisons 
  } else {
    stop("Error: 'res_list_time' not found.")
  }
}

clean_split <- function(x) if (is.null(x) || x == "") NULL else trimws(strsplit(x, ",")[[1]])
COND_ORDER <- clean_split(opt$condition_order)
TIME_ORDER <- clean_split(opt$time_order)
COLOR_COND <- parse_config(opt$cond_colors)
COLOR_TIME <- parse_config(opt$time_colors)

# 推断 Mutants
if (is.null(COND_ORDER)) {
  raw_names <- names(res_list_time)
  all_conds <- unique(gsub("Time_([^_]+)_.*", "\\1", raw_names))
  COND_ORDER <- c(opt$control_cond, setdiff(all_conds, opt$control_cond))
}
mutants <- setdiff(COND_ORDER, opt$control_cond)

if (is.null(TIME_ORDER)) {
  raw_names <- names(res_list_time)
  extracted_times <- unique(gsub("^Time_[^_]+_(.+?)_vs_.*", "\\1", raw_names))
  extracted_times <- extracted_times[!is.na(extracted_times)]
  num_vals <- suppressWarnings(as.numeric(gsub("[^0-9]", "", extracted_times)))
  if(any(is.na(num_vals))) {
      TIME_ORDER <- sort(extracted_times) 
  } else {
      TIME_ORDER <- extracted_times[order(num_vals)]
  }
}

# 剔除 Control Time
PLOT_TIME_ORDER <- setdiff(TIME_ORDER, opt$control_time)

cat(sprintf("  Analysis Plan:\n"))
cat(sprintf("    Comparisons: %s vs %s\n", paste(mutants, collapse=","), opt$control_cond))
cat(sprintf("    Timepoints : %s (Control %s removed)\n", paste(PLOT_TIME_ORDER, collapse=","), opt$control_time))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 6. 生成面板
# ==============================================================================

generate_analysis <- function(direction) {
  cat(sprintf("\n>>> Analyzing: %s-regulated Genes\n", direction))
  
  grob_list <- list()
  
  # --- 1. 表头 (Mutants) ---
  grob_list[[1]] <- nullGrob() 
  for (mut in mutants) {
    col_title <- sprintf("%s vs %s", mut, opt$control_cond)
    grob_list[[length(grob_list)+1]] <- textGrob(col_title, gp=gpar(fontface="bold", fontsize=12))
  }
  
  # --- 2. 数据循环 (仅针对 PLOT_TIME_ORDER) ---
  for (tp in PLOT_TIME_ORDER) {
    # 左侧时间轴
    tp_color <- get_color(tp, COLOR_TIME, "#777777")
    grob_list[[length(grob_list)+1]] <- create_time_header(tp, tp_color)
    
    for (mut in mutants) {
      pattern_A <- paste0("Time_", opt$control_cond, ".*", tp, "_vs_")
      pattern_B <- paste0("Time_", mut, ".*", tp, "_vs_")
      
      name_A <- grep(pattern_A, names(res_list_time), value=TRUE)[1]
      name_B <- grep(pattern_B, names(res_list_time), value=TRUE)[1]
      
      gene_set_list <- list()
      has_data <- FALSE
      
      if (!is.na(name_A) && !is.na(name_B)) {
        genes_A <- get_deg_ids(res_list_time[[name_A]], opt$fdr, opt$lfc, direction)
        genes_B <- get_deg_ids(res_list_time[[name_B]], opt$fdr, opt$lfc, direction)
        
        if (length(genes_A) > 0 || length(genes_B) > 0) {
          has_data <- TRUE
          
          # ================================================================
          # [v3.1] 构建 named list 和 named color vector
          #        名称严格绑定, 后续由 create_venn_grob 做预排序
          # ================================================================
          gene_set_list[[opt$control_cond]] <- genes_A
          gene_set_list[[mut]]              <- genes_B
          
          cols <- c(
            get_color(opt$control_cond, COLOR_COND, "blue"),
            get_color(mut, COLOR_COND, "red")
          )
          names(cols) <- c(opt$control_cond, mut)
          
          # 导出列表 (导出逻辑不受排序影响, 使用原始变量)
          export_venn_lists(genes_A, genes_B, opt$control_cond, mut, tp, direction, opt$outdir)
          
          # 绘图 (create_venn_grob 内部会做预排序)
          grob_list[[length(grob_list)+1]] <- create_venn_grob(gene_set_list, cols)
        }
      }
      
      if (!has_data) {
        grob_list[[length(grob_list)+1]] <- textGrob("No DEGs", gp=gpar(col="grey70", fontsize=9))
      }
    }
  }
  
  # --- 3. 布局与保存 ---
  n_cols <- length(mutants) + 1
  widths <- c(0.15, rep(1, length(mutants)))
  heights <- c(0.1, rep(1, length(PLOT_TIME_ORDER)))
  
  pdf_file <- file.path(opt$outdir, paste0(opt$prefix, "_Venn_Panel_", direction, ".pdf"))
  
  pdf(pdf_file, width = 4 * length(mutants) + 1.5, height = 3.5 * length(PLOT_TIME_ORDER) + 1)
  
  grid.arrange(
    grobs = grob_list,
    ncol = n_cols,
    widths = unit(widths, "null"),
    heights = unit(heights, "null"),
    top = textGrob(
      sprintf("%s-regulated", direction),
      gp = gpar(fontsize=20, fontface="bold"),
      vjust = 1
    ),
    bottom = textGrob(
      sprintf("Gene lists saved in: %s/Venn_GeneLists/", basename(opt$outdir)),
      gp = gpar(fontsize=8, col="grey50", fontface="italic")
    )
  )
  
  dev.off()
  cat(sprintf("  ✓ Panel Saved: %s\n", basename(pdf_file)))
}

generate_analysis("Up")
generate_analysis("Down")

cat("\n>>> Venn Analysis & Export Complete.\n")