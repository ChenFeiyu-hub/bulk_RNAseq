#!/usr/bin/env Rscript
################################################################################
# RNAseq Pipeline 5c: Temporal Cascade & Fate Composition (v2.1.1 Peak-Sorting)
# 
# 功能:
#   1. [核心升级] 计算 Cluster LOESS 插值后的绝对峰值时间 (Peak Time) 进行排序
#   2. 划分 Stage (Early/Middle/Late) 并与峰值时间对齐
#   3. 整合 Pipeline 5b 的 Fate 结果
#   4. 绘制 ComplexHeatmap:
#      - 主体: LOESS 平滑后的表达趋势 (Cluster Center)
#      - 右侧: 堆叠柱状图展示 Fate 构成 (Conserved vs Changed 基因数)
#
# v2.1.1 修复:
#   - 修复 tryCatch error handler 中变量赋值的作用域 (scoping) Bug
#   - 改用 tryCatch 返回值模式，确保 fallback 结果正确传播
#   - 增加 LOESS degree 自适应逻辑 (4点时降级为 degree=1)
#
# 依赖: Mfuzz_v10_Results.RData (由 Pipeline 5a/5b 生成)
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(ComplexHeatmap)
  library(circlize)
  library(RColorBrewer)
  library(grid)
  library(Biobase)
  library(Mfuzz) 
})

################################################################################
# 1. 参数解析与配置
################################################################################

parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  params <- list()
  
  # 默认值
  params$rdata <- "05_Mfuzz/rdata/Mfuzz_v11_Results.RData"
  params$out_dir <- "05_Mfuzz/Pipeline5c_Temporal"
  params$stage_ranges <- "00-03,03-12,12-24"
  params$stage_names <- "Early,Middle,Late"
  params$stage_colors <- "#D1C4E9,#7E57C2,#311B92"
  params$n_pseudotime <- 100
  params$loess_span <- 0.4
  params$membership_cutoff <- 0.5
  
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (grepl("^--", arg)) {
      key <- sub("^--", "", arg)
      val <- args[i+1]
      params[[key]] <- val
      i <- i + 2
    } else {
      i <- i + 1
    }
  }
  return(params)
}

PARAMS <- parse_arguments()

# 创建输出目录
dir.create(PARAMS$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PARAMS$out_dir, "figures"), recursive = TRUE, showWarnings = FALSE)

# 日志
log_file <- file.path(PARAMS$out_dir, "Pipeline5c.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat("================================================================================\n")
cat("  Pipeline 5c: Temporal Cascade & Fate Composition (v2.1.1 Peak-Sorting)\n")
cat("================================================================================\n")
cat(sprintf("RData: %s\n", PARAMS$rdata))

################################################################################
# 2. 数据加载与处理
################################################################################

if (!file.exists(PARAMS$rdata)) stop("RData file not found!")
load(PARAMS$rdata) 

# 解析 Stage 配置
parse_ranges <- function(range_str) {
  parts <- strsplit(range_str, ",")[[1]]
  ranges <- lapply(parts, function(x) {
    as.numeric(strsplit(x, "-")[[1]])
  })
  return(ranges)
}

stage_ranges_list <- parse_ranges(PARAMS$stage_ranges)
stage_names_vec <- strsplit(PARAMS$stage_names, ",")[[1]]
stage_colors_vec <- strsplit(PARAMS$stage_colors, ",")[[1]]
names(stage_colors_vec) <- stage_names_vec

# Fate 颜色配置
FATE_COLORS <- c(
  "Conserved" = "#98DF8A",
  "Shape_Diverged" = "#FF9896",
  "Amp_Altered" = "#FFBB78",
  "Phase_Shifted" = "#F7B6D2"
)

################################################################################
# 3. 核心计算函数 (v2.1.1 修复)
################################################################################

# ★ [修复] 自适应 LOESS 拟合辅助函数
# 根据唯一时间点数选择合适的拟合策略:
#   >= 5 点: LOESS degree=2 (完整二次局部拟合)
#   == 4 点: LOESS degree=1 (线性局部拟合, 避免 span*n < degree+1 导致失败)
#   <  4 点: 自然样条插值 (splinefun)
fit_smooth_curve <- function(y, time_points, t_interp, span = 0.4) {
  n_unique <- length(unique(time_points))
  
  if (n_unique >= 5) {
    # 足够多的点: 使用 LOESS degree=2
    fit <- loess(y ~ time_points, span = span, degree = 2)
    predict(fit, newdata = data.frame(time_points = t_interp))
  } else if (n_unique >= 4) {
    # 边界情况: 降级为 degree=1, 并适当增大 span 确保稳定性
    adjusted_span <- max(span, 0.5)
    fit <- loess(y ~ time_points, span = adjusted_span, degree = 1)
    predict(fit, newdata = data.frame(time_points = t_interp))
  } else {
    # 点数不足: 使用自然样条
    fit <- splinefun(time_points, y, method = "natural")
    fit(t_interp)
  }
}

# ★ [修复] 计算插值后的绝对峰值时间 (Visual Peak Sorting)
# 核心修复: 用 tryCatch 的返回值捕获 y_pred, 而非在 error handler 内部赋值
calc_peak_time_loess <- function(centers, time_points, n_interp = 100, span = 0.4) {
  n_clusters <- nrow(centers)
  peak_times <- numeric(n_clusters)
  t_interp <- seq(min(time_points), max(time_points), length.out = n_interp)
  
  for (i in 1:n_clusters) {
    y <- as.numeric(centers[i, ])
    
    # ★ 修复: tryCatch 返回值模式 —— 无论成功还是 fallback, 结果都赋给 y_pred
    y_pred <- tryCatch({
      fit_smooth_curve(y, time_points, t_interp, span)
    }, error = function(e) {
      # 终极 fallback: 线性插值 (永远不会失败)
      approx(time_points, y, xout = t_interp, rule = 2)$y
    })
    
    # 处理可能的 NA (LOESS 边界问题)
    if (any(is.na(y_pred))) {
      y_pred[is.na(y_pred)] <- min(y, na.rm = TRUE)
    }
    
    # 找到最大值所在的索引, 映射回时间轴
    max_idx <- which.max(y_pred)
    peak_times[i] <- t_interp[max_idx]
  }
  return(peak_times)
}

# 判定 Stage (基于 Peak Time)
classify_stage <- function(peak_time, ranges, names) {
  for (i in seq_along(ranges)) {
    r <- ranges[[i]]
    if (peak_time >= r[1] && peak_time <= r[2]) {
      return(names[i])
    }
  }
  return("Unclassified") 
}

# ★ [修复] LOESS 插值 (用于绘图数据)
# 核心修复: 用 tryCatch 返回值赋给 interp_matrix[i, ], 而非在 error handler 内部赋值
interpolate_expression <- function(expr_matrix, time_points, n_interp = 100, span = 0.4) {
  t_interp <- seq(min(time_points), max(time_points), length.out = n_interp)
  interp_matrix <- matrix(NA, nrow = nrow(expr_matrix), ncol = n_interp)
  rownames(interp_matrix) <- rownames(expr_matrix)
  
  for (i in 1:nrow(expr_matrix)) {
    y <- as.numeric(expr_matrix[i, ])
    
    # ★ 修复: 赋值在 tryCatch 外部, 接收其返回值
    interp_matrix[i, ] <- tryCatch({
      fit_smooth_curve(y, time_points, t_interp, span)
    }, error = function(e) {
      # 终极 fallback: 线性插值
      approx(time_points, y, xout = t_interp, rule = 2)$y
    })
  }
  return(list(matrix = interp_matrix, time = t_interp))
}

################################################################################
# 4. 构建绘图数据
################################################################################

cat("--- Processing Temporal Order and Fate Composition ---\n")

# 4.1 计算时序 (使用 Peak Time)
time_points <- MFUZZ_PARAMS$time_points
wt_centers <- wt_clusters$centers

cat(sprintf("  Time points detected: %d (%s)\n", 
            length(unique(time_points)), 
            paste(sort(unique(time_points)), collapse = ", ")))

cat("  Using Peak-Time Sorting (finding max of LOESS curve)...\n")
peak_times <- calc_peak_time_loess(
    wt_centers, 
    time_points, 
    n_interp = as.numeric(PARAMS$n_pseudotime), 
    span = as.numeric(PARAMS$loess_span)
)

# 4.2 构建 Cluster Info 表
n_clusters <- nrow(wt_centers)
cluster_info <- data.frame(
  Cluster = 1:n_clusters,
  Peak_Time = peak_times,
  stringsAsFactors = FALSE
)

# 划分 Stage
cluster_info$Stage <- sapply(cluster_info$Peak_Time, function(t) {
  classify_stage(t, stage_ranges_list, stage_names_vec)
})
cluster_info$Stage <- factor(cluster_info$Stage, levels = stage_names_vec)

# 4.3 提取 Fate 计数 (用于堆叠柱状图)
fate_types <- c("Conserved", "Shape_Diverged", "Amp_Altered", "Phase_Shifted")
fate_counts_mat <- matrix(0, nrow = n_clusters, ncol = length(fate_types))
colnames(fate_counts_mat) <- fate_types
rownames(fate_counts_mat) <- paste0("C", 1:n_clusters)

for (i in 1:n_clusters) {
  if (i > length(fate_results)) next
  res <- fate_results[[i]]
  
  if (!is.null(res) && res$status == "analyzed") {
    labels <- res$fate_labels
    sizes <- res$subcluster_sizes
    
    for (k in seq_along(labels)) {
      lbl <- labels[k]
      sz <- sizes[k]
      
      if (lbl == "Conserved") {
        fate_counts_mat[i, "Conserved"] <- fate_counts_mat[i, "Conserved"] + sz
      } else if (grepl("Shape", lbl)) {
        fate_counts_mat[i, "Shape_Diverged"] <- fate_counts_mat[i, "Shape_Diverged"] + sz
      } else if (grepl("Phase", lbl)) {
        fate_counts_mat[i, "Phase_Shifted"] <- fate_counts_mat[i, "Phase_Shifted"] + sz
      } else if (grepl("Amp", lbl) || grepl("Repressed", lbl) || grepl("Enhanced", lbl)) {
        fate_counts_mat[i, "Amp_Altered"] <- fate_counts_mat[i, "Amp_Altered"] + sz
      }
    }
  }
}

# 4.4 排序 Cluster
# 先按 Stage 排序，再按 Peak_Time 排序
cluster_info$Order_Idx <- order(cluster_info$Stage, cluster_info$Peak_Time)
ordered_info <- cluster_info[cluster_info$Order_Idx, ]

cat(sprintf("  Top 3 earliest clusters: C%s (%.2f), C%s (%.2f), C%s (%.2f)\n",
            ordered_info$Cluster[1], ordered_info$Peak_Time[1],
            ordered_info$Cluster[2], ordered_info$Peak_Time[2],
            ordered_info$Cluster[3], ordered_info$Peak_Time[3]))

################################################################################
# 5. 准备 Heatmap 数据
################################################################################

cat("--- Generating Heatmap Data (Cluster Centers) ---\n")

# 计算每个 Cluster 的插值中心 (用于绘制热图主体的平滑曲线)
cluster_centers_interp <- matrix(NA, nrow = n_clusters, ncol = as.numeric(PARAMS$n_pseudotime))
rownames(cluster_centers_interp) <- paste0("C", 1:n_clusters)

# 执行插值
interp_res <- interpolate_expression(
  wt_centers, 
  time_points, 
  n_interp = as.numeric(PARAMS$n_pseudotime), 
  span = as.numeric(PARAMS$loess_span)
)
mat_smooth <- interp_res$matrix

# Z-score standardization per cluster
mat_centers_z <- t(scale(t(mat_smooth)))
mat_centers_z[mat_centers_z > 2.5] <- 2.5
mat_centers_z[mat_centers_z < -2.5] <- -2.5

# 按排序重新排列矩阵
mat_plot <- mat_centers_z[ordered_info$Cluster, ]
fate_plot <- fate_counts_mat[ordered_info$Cluster, ]
stage_plot <- ordered_info$Stage
row_labels <- paste0("C", ordered_info$Cluster)

################################################################################
# 6. 绘制 ComplexHeatmap
################################################################################

cat("--- Plotting ComplexHeatmap ---\n")

# 6.1 左侧: Stage 颜色条
anno_left <- rowAnnotation(
  Stage = stage_plot,
  col = list(Stage = stage_colors_vec),
  show_legend = TRUE,
  show_annotation_name = FALSE,
  width = unit(5, "mm")
)

# 6.2 右侧: Stacked Barplot (Fate Composition)
anno_right <- rowAnnotation(
  "Fate Composition\n(Gene Counts)" = anno_barplot(
    fate_plot,
    gp = gpar(fill = FATE_COLORS[colnames(fate_plot)], col = NA),
    bar_width = 0.8,
    width = unit(4, "cm"), 
    axis_param = list(side = "top", labels_rot = 0, gp = gpar(fontsize = 8))
  ),
  show_annotation_name = TRUE,
  annotation_name_rot = 0,
  annotation_name_side = "top",
  annotation_name_gp = gpar(fontsize = 9, fontface = "bold")
)

# 6.3 绘制
ht <- Heatmap(
  mat_plot,
  name = "Expression\n(Z-score)",
  col = colorRamp2(c(-2, 0, 2), c("#2166AC", "#F7F7F7", "#B2182B")),
  
  # 行配置
  cluster_rows = FALSE,
  show_row_dend = FALSE,
  row_labels = row_labels,
  row_names_side = "left",
  row_names_gp = gpar(fontsize = 9),
  
  # 列配置
  cluster_columns = FALSE,
  show_column_names = FALSE, 
  column_title = "Temporal Cascade (Peak-Time Ordered)",
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  
  # 分割 (按 Stage)
  row_split = stage_plot,
  row_gap = unit(2, "mm"),
  row_title_rot = 0,
  row_title_gp = gpar(fontsize = 10, fontface = "bold"),
  
  # 注释
  left_annotation = anno_left,
  right_annotation = anno_right,
  
  # 边框
  border = TRUE
)

# 6.4 保存
pdf_file <- file.path(PARAMS$out_dir, "figures", "Cascade_Heatmap_Interpolated.pdf")
pdf(pdf_file, width = 11, height = max(6, n_clusters * 0.45))
draw(ht, merge_legend = TRUE)
dev.off()

cat(sprintf("  ✔ Saved: %s\n", pdf_file))

################################################################################
# 7. 导出统计表
################################################################################

cat("--- Exporting Summary Tables ---\n")

summary_df <- ordered_info
summary_df <- cbind(summary_df, fate_counts_mat[ordered_info$Cluster, ])

write.table(summary_df, 
            file.path(PARAMS$out_dir, "Cluster_Temporal_Fate_Summary.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("================================================================================\n")
cat("  Pipeline 5c Completed Successfully\n")
cat("================================================================================\n")

sink(type = "message")
sink(type = "output")
close(log_con)