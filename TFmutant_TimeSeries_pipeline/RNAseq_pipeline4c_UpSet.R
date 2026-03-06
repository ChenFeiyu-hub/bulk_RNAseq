#!/usr/bin/env Rscript
################################################################################
# RNAseq_pipeline4c_UpSet.R
# 
# 功能: 标准化 UpSet 图绘制 (Logic 1: Time Course)
################################################################################

suppressPackageStartupMessages({
  library(UpSetR)
  library(tidyverse)
  library(optparse)
  library(RColorBrewer)
  library(grid)
})

# ==============================================================================
# 1. 参数定义
# ==============================================================================
option_list <- list(
  make_option(c("-i", "--input"), type = "character", help = "Input Logic1 RData"),
  make_option(c("-o", "--outdir"), type = "character", default = "UpSet_Output"),
  make_option(c("-p", "--prefix"), type = "character", default = "Project"),
  make_option(c("--condition_order"), type = "character", default=NULL),
  make_option(c("--time_order"), type = "character", default=NULL),
  make_option(c("--cond_colors"), type = "character", help = "Key:Value colors"),
  make_option(c("--fdr"), type = "double", default = 0.05),
  make_option(c("--lfc"), type = "double", default = 1.0)
)

opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$input)) stop("Error: --input is required.")

# ==============================================================================
# 2. 辅助函数
# ==============================================================================

parse_color_config <- function(config_str) {
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

draw_std_upset <- function(gene_list, set_metadata, title, output_file) {
  if (length(gene_list) < 2) {
    cat(sprintf("  [Skip] %s: Less than 2 sets have DEGs.\n", title))
    return()
  }
  
  upset_data <- fromList(gene_list)
  ordered_sets <- set_metadata$Display_Name
  ordered_colors <- set_metadata$Color
  
  valid_indices <- ordered_sets %in% names(gene_list)
  ordered_sets <- ordered_sets[valid_indices]
  ordered_colors <- ordered_colors[valid_indices]
  
  if(length(ordered_sets) < 2) return()

  max_set_size <- max(sapply(gene_list, length))
  
  # [修复核心] onefile = FALSE 防止出现空白第一页
  pdf(output_file, width = 14, height = 9, onefile = FALSE)
  
  print(
    upset(upset_data,
          sets = ordered_sets,
          keep.order = TRUE,
          sets.bar.color = ordered_colors,
          order.by = "freq",
          decreasing = TRUE,
          mb.ratio = c(0.6, 0.4),
          point.size = 3.5,
          line.size = 1,
          main.bar.color = "black",
          matrix.color = "black",
          text.scale = c(1.5, 1.2, 1.2, 1, 1.2, 1.2),
          set_size.show = TRUE,
          set_size.scale_max = max_set_size * 1.25
    )
  )
  
  grid.text(title, x = 0.5, y = 0.98, gp = gpar(fontsize = 16, fontface = "bold"))
  dev.off()
  cat(sprintf("  ✓ Saved: %s\n", basename(output_file)))
}

# ==============================================================================
# 3. 数据加载与元数据构建
# ==============================================================================
cat(">>> Loading Logic 1 Data...\n")
load(opt$input)

if (!exists("res_list_time")) stop("Error: 'res_list_time' not found in RData.")

# 参数预处理
clean_split <- function(x) {
  if (is.null(x) || x == "") return(NULL)
  return(trimws(strsplit(x, ",")[[1]]))
}

COND_ORDER <- clean_split(opt$condition_order)
TIME_ORDER <- clean_split(opt$time_order)
if(length(TIME_ORDER) == 1 && TIME_ORDER == "auto") TIME_ORDER <- NULL

COLOR_MAP <- parse_color_config(opt$cond_colors)

# 动态解析
raw_names <- names(res_list_time)
meta_list <- list()

cat(">>> Parsing comparison names...\n")
regex_pattern <- "^(?:Time|TimeCourse)_([^_]+)_(.+?)_vs_"

for (rname in raw_names) {
  matches <- str_match(rname, regex_pattern)
  if (all(!is.na(matches[1,]))) {
    cond <- matches[1,2]
    tp   <- matches[1,3]
    display_name <- paste0(cond, "_", tp)
    meta_list[[rname]] <- data.frame(
      Original_Name = rname, Display_Name = display_name,
      Condition = cond, Time = tp, stringsAsFactors = FALSE
    )
  }
}

if (length(meta_list) == 0) stop("Failed to parse comparison names.")
set_info <- do.call(rbind, meta_list)

# 智能排序
if (!is.null(COND_ORDER)) {
  set_info <- set_info[set_info$Condition %in% COND_ORDER, ]
  set_info$Condition <- factor(set_info$Condition, levels = COND_ORDER)
} else {
  cat("  ! Note: --condition_order is empty. Auto-detecting conditions...\n")
  set_info$Condition <- factor(set_info$Condition)
}

if (is.null(TIME_ORDER)) {
  unique_times <- unique(set_info$Time)
  num_vals <- suppressWarnings(as.numeric(gsub("[^0-9]", "", unique_times)))
  if(any(is.na(num_vals))) {
      TIME_ORDER <- sort(unique_times) 
  } else {
      TIME_ORDER <- unique_times[order(num_vals)]
  }
}
set_info$Time <- factor(set_info$Time, levels = TIME_ORDER)

set_info <- set_info[!is.na(set_info$Condition) & !is.na(set_info$Time), ]
set_info <- set_info %>% arrange(Condition, Time)

if (nrow(set_info) == 0) stop("Error: All sets filtered out! Check config vs data names.")

# 颜色分配
set_info$Color <- sapply(as.character(set_info$Condition), function(cond) {
  if (!is.null(COLOR_MAP) && cond %in% names(COLOR_MAP)) return(COLOR_MAP[[cond]])
  return("#999999") 
})

cat(">>> Set Metadata Constructed (Sorted):\n")
print(set_info[, c("Display_Name", "Condition", "Time", "Color")])

# ==============================================================================
# 4. 数据提取与绘图
# ==============================================================================
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

cat("\n>>> Processing Up-regulated Genes...\n")
list_up <- list()
for (i in 1:nrow(set_info)) {
  orig_name <- set_info$Original_Name[i]
  disp_name <- set_info$Display_Name[i]
  res_obj <- res_list_time[[orig_name]]
  if (!is.null(res_obj)) {
    degs <- get_deg_ids(res_obj, opt$fdr, opt$lfc, "Up")
    if (length(degs) > 0) list_up[[disp_name]] <- degs
  }
}
valid_sets_up <- names(list_up)
meta_up <- set_info[set_info$Display_Name %in% valid_sets_up, ]

if (nrow(meta_up) > 0) {
  draw_std_upset(list_up, meta_up, sprintf("Up-regulated (FDR<%.2f, LFC>%.1f)", opt$fdr, opt$lfc),
                 file.path(opt$outdir, paste0(opt$prefix, "_UpSet_Upregulated.pdf")))
} else {
  cat("  [Warning] No Up-regulated genes found.\n")
}

cat("\n>>> Processing Down-regulated Genes...\n")
list_down <- list()
for (i in 1:nrow(set_info)) {
  orig_name <- set_info$Original_Name[i]
  disp_name <- set_info$Display_Name[i]
  res_obj <- res_list_time[[orig_name]]
  if (!is.null(res_obj)) {
    degs <- get_deg_ids(res_obj, opt$fdr, opt$lfc, "Down")
    if (length(degs) > 0) list_down[[disp_name]] <- degs
  }
}
valid_sets_down <- names(list_down)
meta_down <- set_info[set_info$Display_Name %in% valid_sets_down, ]

if (nrow(meta_down) > 0) {
  draw_std_upset(list_down, meta_down, sprintf("Down-regulated (FDR<%.2f, LFC<-%.1f)", opt$fdr, opt$lfc),
                 file.path(opt$outdir, paste0(opt$prefix, "_UpSet_Downregulated.pdf")))
} else {
  cat("  [Warning] No Down-regulated genes found.\n")
}

cat("\n>>> Analysis Complete.\n")