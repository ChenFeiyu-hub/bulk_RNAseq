#!/usr/bin/env Rscript
################################################################################
# RNAseq Pipeline 5: Mfuzz时间序列聚类分析 (工程化版本)
# 
# 版本: 11.1 (标准化重构版本)
# 
# ============================================================================
# 版本更新记录
# ============================================================================
#
# 【v11.1 标准化重构】
#   1. 变量名标准化：所有 dof 相关变量名替换为 mut
#      - c_dof -> c_mut, m_dof -> m_mut
#      - dof_center -> mut_center, dof_genes -> mut_genes, etc.
#   
#   2. 可视化显示参数化：
#      - 新增 condition_wt_display: WT 显示名称
#      - 新增 use_italic_wt / use_italic_mut: 控制是否使用斜体格式
#      - 所有图形标题、图例由 DISPLAY 参数控制
#      - 斜体格式不再硬编码，而是由参数控制
#   
#   3. 输出文件名保持标准化命名 (WT/MUT)
#      - 例如: WT_C01_MUT_C08_Profile.pdf
#      - 便于下游脚本解析
#   
#   4. 向后兼容：
#      - 仍支持 --c_dof / --m_dof 参数名 (内部转换为 c_mut/m_mut)
#
# 【v11.0 PP-XB算法更新】
#   基于两份独立研究报告的交叉验证，实现了统一的最优聚类数选择算法：
#   
#   1. Kwon指数替代XB：从根本上解决XB单调递减问题
#      - 惩罚项 = (1/c) * Σ||v_j - v̄||² 随c增加而增加
#      - 抵消紧致度的自然下降，使曲线具备U形潜力
#   
#   2. 曲线形态分类重构：三类自适应策略
#      - U_SHAPED_OR_VALLEY: 直接选择全局最小值
#      - MONOTONIC_DECREASING: 使用曲率最大点（替代PAVA+Kneedle）
#      - FLUCTUATING: LOWESS平滑后选择最小值
#   
#   3. √N样本量约束：基于Schwämmle & Jensen (2010)
#      - Pipeline5a: c_max ≤ min(√N, N/10)
#      - Pipeline5b: 更保守，N<50时 c_max ≤ 4
#   
#   4. 双面板诊断图：同时显示原始XB和修正后Kwon曲线
#   
#   5. 始终报告原始指标值：避免PAVA平滑值的误导
#
# 【v10.1 工程化更新】
#   1. 自动化数据衔接：仅需输入 RData 路径即可运行
#   2. 自动时间点解析：从 metadata 自动获取并排序时间点
#   3. 自动输出路径：基于项目结构自动创建输出目录
#   4. 支持任意数量时间点：不再硬编码时间点数目
#   5. 完全保留 v10.0 的 DDTW 正交三维指标体系
#
# 【使用方式】
#   简化模式 (推荐):
#     Rscript RNAseq_pipeline5a_mfuzz.R \
#       --rdata /path/to/04_DESeq2/RData/Pipeline4_Logic1_TimeCourse_for_Mfuzz.RData
#
#   完整参数模式 (包含显示控制):
#     Rscript RNAseq_pipeline5a_mfuzz.R \
#       --rdata /path/to/RData.RData \
#       --samples /path/to/samples.txt \
#       --out_dir /path/to/output \
#       --condition_wt WT --condition_mut dof \
#       --condition_wt_display WT --condition_mut_display dof \
#       --use_italic_wt FALSE --use_italic_mut TRUE \
#       --c_wt 10 --c_mut 12 ...
#
# ============================================================================
# v10.0 DDTW正交三维指标体系 (保留)
# ============================================================================
#
# 【理论基础】
#   基于 Keogh & Pazzani (2001) 的 Derivative DTW 思想，实现形状-振幅-相位
#   的数学正交分解。
#
# 【1. Shape (形状距离) - DDTW】
#   - 对导数序列进行DTW对齐
#   - D_shape = DDTW_distance / 理论最大距离
#   - 正交性：通过比较"变化率"而非原始值，对幅度缩放不敏感
#   - 通过允许时间弹性对齐，剔除相位差异的影响
#
# 【2. Amplitude (振幅距离) - SDR】
#   - 定义：峰谷差(Range)的归一化差异
#   - D_amp = |Range_A - Range_B| / (Range_A + Range_B)
#   - 正交性：时间顺序无关的标量统计量
#   - 天然与Shape和Phase正交
#
# 【3. Phase (相位距离) - Warping Path Deviation】
#   - 定义：DTW路径偏离对角线的程度
#   - D_phase = mean(|path_i - path_j × n/m|) / 最大可能偏离
#   - 正交性：从DTW对齐中自然得出，与Shape解耦
#   - 测量"需要多少时间扭曲才能对齐形状"
#
################################################################################

suppressPackageStartupMessages({
  library(Mfuzz)
  library(Biobase)
  library(SummarizedExperiment)
  library(tidyverse)
  library(RColorBrewer)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(grid)
})

if (!require("clue", quietly = TRUE)) {
  install.packages("clue", repos = "https://cloud.r-project.org")
  library(clue)
}

################################################################################
# 工程化参数解析 (v10.1 新增)
################################################################################

parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  
  params <- list(
    # ========== 工程化核心参数 (v10.1 新增) ==========
    rdata = NULL,                # [必需] 上游 RData 路径
    samples_file = NULL,         # [可选] samples.txt 路径 (覆盖 RData 中的 metadata)
    out_dir = NULL,              # [可选] 输出目录 (自动推断)
    project_name = NULL,         # [可选] 项目名称 (自动推断)
    
    # ========== 条件标签配置 (v11.1 更新 - 添加完整显示控制) ==========
    condition_wt = "WT",           # WT条件的标签 (用于数据过滤)
    condition_mut = "mut",         # 突变体条件的标签 (用于数据过滤)
    condition_wt_display = "WT",   # WT显示名称 (用于图形)
    condition_mut_display = "mut", # 突变体显示名称 (用于图形)
    use_italic_wt = FALSE,         # WT名称是否使用斜体
    use_italic_mut = TRUE,         # 突变体名称是否使用斜体 (基因名通常斜体)
    
    # ========== 运行模式控制 ==========
    run_pipeline5a = TRUE,
    run_pipeline5b = TRUE,
    
    # ========== Pipeline5a Mfuzz参数 ==========
    c_wt = NULL,
    c_mut = NULL,
    c_range_min = 6,
    c_range_max = 24,
    m_wt = NULL,
    m_mut = NULL,
    m_lower = 1.25,
    m_upper = 2.00,
    membership_cutoff = 0.5,
    min_std = 0.5,
    
    # ========== 成本函数权重 (DDTW) ==========
    w_shape = 0.5,
    w_amp = 0.3,
    w_phase = 0.2,
    
    # ========== DDTW参数 ==========
    interp_points = 100,
    dtw_window = NULL,
    derivative_method = "central",
    
    # ========== Pipeline5b参数 ==========
    p5b_c_range_min = 2,
    p5b_c_range_max = 8,
    p5b_min_genes = 20,
    
    # ========== XB-Based参数 ==========
    global_prominence_threshold = 0.10,
    subcluster_prominence_threshold = 0.10,
    kneedle_S = 1.0,
    
    # ========== 统一分类阈值 ==========
    classify_shape_threshold = 0.15,
    classify_amp_threshold = 0.12,
    classify_phase_threshold = 0.25
  )
  
  # 解析命令行参数
  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    
    if (grepl("^--", arg)) {
      key <- sub("^--", "", arg)
      
      # 处理 --key=value 格式
      if (grepl("=", key)) {
        parts <- strsplit(key, "=")[[1]]
        key <- parts[1]
        value <- parts[2]
      } else {
        # 处理 --key value 格式
        if (i + 1 <= length(args) && !grepl("^--", args[i + 1])) {
          i <- i + 1
          value <- args[i]
        } else {
          value <- "TRUE"  # 布尔标志
        }
      }
      
      # 转换键名 (支持下划线和连字符)
      key <- gsub("-", "_", key)
      
      # v11.1: 向后兼容旧参数名 (c_dof -> c_mut, m_dof -> m_mut)
      if (key == "c_dof") key <- "c_mut"
      if (key == "m_dof") key <- "m_mut"
      
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

#' 从时间标签中提取数值（增强版）
#' 
#' 支持格式: "00h", "03h", "T0", "T03", "0", "3", "0h", "3hr", "006h" 等
#' @param time_label 时间标签字符串
#' @return 数值
extract_time_numeric <- function(time_label) {
  if (is.na(time_label) || time_label == "") {
    return(NA)
  }
  
  # 转为字符串
  time_label <- as.character(time_label)
  
  # 方法：移除所有非数字字符，只保留数字和小数点
  cleaned <- gsub("[^0-9.]", "", time_label)
  
  # 处理空字符串（如果原标签只有字母）
  if (cleaned == "" || cleaned == ".") {
    warning(sprintf("无法解析时间标签: '%s'", time_label))
    return(NA)
  }
  
  # 转换为数值（自动处理前导零，如 "00" -> 0, "03" -> 3）
  num <- suppressWarnings(as.numeric(cleaned))
  
  if (is.na(num)) {
    warning(sprintf("无法解析时间标签: '%s' -> '%s'", time_label, cleaned))
  }
  
  return(num)
}

#' 从 metadata 自动获取并排序时间点
#' 
#' @param metadata 元数据 data.frame
#' @return list 包含 time_points (数值向量) 和 time_labels (原始标签)
auto_detect_timepoints <- function(metadata) {
  # Find time column
  time_col <- NULL
  for (col_name in c("Time", "time", "Timepoint", "timepoint", "TIME")) {
    if (col_name %in% colnames(metadata)) {
      time_col <- col_name
      break
    }
  }
  
  if (is.null(time_col)) {
    stop("Cannot find time column in metadata (tried: Time, time, Timepoint)")
  }
  
  # Get unique time labels
  unique_times <- unique(as.character(metadata[[time_col]]))
  
  # Extract numeric values and sort
  time_df <- data.frame(
    label = unique_times,
    numeric = sapply(unique_times, extract_time_numeric),
    stringsAsFactors = FALSE
  )
  
  # Remove unparseable time points
  time_df <- time_df[!is.na(time_df$numeric), ]
  
  # Sort by numeric value
  time_df <- time_df[order(time_df$numeric), ]
  
  cat(sprintf("  Auto-detected %d time points:\n", nrow(time_df)))
  cat(sprintf("    Labels: %s\n", paste(time_df$label, collapse = ", ")))
  cat(sprintf("    Numeric: %s\n", paste(time_df$numeric, collapse = ", ")))
  
  return(list(
    time_points = time_df$numeric,
    time_labels_raw = time_df$label,
    time_labels = paste0("T", time_df$numeric),  # Standardized labels
    time_col = time_col,
    D = nrow(time_df)
  ))
}

#' 生成时间标签格式映射
#' 
#' 用于匹配 metadata 中的时间格式
#' @param time_info 从 auto_detect_timepoints 返回的 list
#' @return 命名向量，名称为标准格式，值为原始格式
create_time_format_map <- function(time_info) {
  map <- setNames(time_info$time_labels_raw, time_info$time_labels)
  return(map)
}

################################################################################
# 路径自动推断 (v10.1 新增)
################################################################################

#' 从 RData 路径自动推断项目结构
#' 
#' @param rdata_path RData 文件路径
#' @return list 包含 base_dir, deseq2_dir, output_dir, project_name
infer_project_paths <- function(rdata_path) {
  # 规范化路径
  rdata_path <- normalizePath(rdata_path, mustWork = TRUE)
  
  # 获取目录结构
  rdata_dir <- dirname(rdata_path)   # .../04_DESeq2/RData
  deseq2_dir <- dirname(rdata_dir)   # .../04_DESeq2
  base_dir <- dirname(deseq2_dir)    # .../项目根目录
  
  # 检查目录结构
  if (!grepl("RData$", rdata_dir)) {
    warning("RData 文件不在标准的 'RData' 子目录中，路径推断可能不准确")
  }
  
  # 推断项目名称 (从 RData 文件名或目录名)
  rdata_basename <- basename(rdata_path)
  if (grepl("_DESeq2_objects\\.RData$", rdata_basename)) {
    project_name <- sub("_DESeq2_objects\\.RData$", "", rdata_basename)
  } else if (grepl("Pipeline4_Logic1", rdata_basename)) {
    # 从父目录名推断
    project_name <- basename(base_dir)
  } else {
    project_name <- tools::file_path_sans_ext(rdata_basename)
  }
  
  # 输出目录
  output_dir <- file.path(base_dir, "05_Mfuzz")
  
  return(list(
    rdata_path = rdata_path,
    rdata_dir = rdata_dir,
    deseq2_dir = deseq2_dir,
    base_dir = base_dir,
    output_dir = output_dir,
    project_name = project_name
  ))
}

################################################################################
# 验证和初始化
################################################################################

# 验证必需参数
if (is.null(PARAMS$rdata)) {
  cat("================================================================================\n")
  cat("  Pipeline 5: Mfuzz 时间序列聚类分析 (v11.1 标准化版本)\n")
  cat("================================================================================\n\n")
  cat("错误: 必须指定 --rdata 参数\n\n")
  cat("使用方法:\n")
  cat("  简化模式 (推荐):\n")
  cat("    Rscript RNAseq_pipeline5a_mfuzz.R \\\n")
  cat("      --rdata /path/to/04_DESeq2/RData/Pipeline4_Logic1_TimeCourse_for_Mfuzz.RData\n\n")
  cat("  完整参数 (包含显示控制):\n")
  cat("    Rscript RNAseq_pipeline5a_mfuzz.R \\\n")
  cat("      --rdata /path/to/RData.RData \\\n")
  cat("      --samples_file /path/to/samples.txt \\\n")
  cat("      --out_dir /path/to/output \\\n")
  cat("      --condition_wt WT --condition_mut dof \\\n")
  cat("      --condition_wt_display WT --condition_mut_display dof \\\n")
  cat("      --use_italic_wt FALSE --use_italic_mut TRUE \\\n")
  cat("      --c_wt 10 --c_mut 12 \\\n")
  cat("      --run_pipeline5a TRUE --run_pipeline5b TRUE\n\n")
  cat("核心参数说明:\n")
  cat("  --rdata                 [必需] 上游 DESeq2 输出的 RData 文件\n")
  cat("  --samples_file          [可选] samples.txt 元数据文件\n")
  cat("  --out_dir               [可选] 输出目录 (默认自动推断为 05_Mfuzz)\n")
  cat("  --condition_wt          [可选] WT条件标签 (默认: WT)\n")
  cat("  --condition_mut         [可选] 突变体条件标签 (默认: mut)\n")
  cat("  --condition_wt_display  [可选] WT显示名称 (默认: WT)\n")
  cat("  --condition_mut_display [可选] 突变体显示名称 (默认: mut)\n")
  cat("  --use_italic_wt         [可选] WT名称是否斜体 (默认: FALSE)\n")
  cat("  --use_italic_mut        [可选] 突变体名称是否斜体 (默认: TRUE)\n")
  cat("  --c_wt/--c_mut          [可选] 固定聚类数 (默认自动选择)\n")
  cat("  --run_pipeline5a        [可选] 运行 WT-MUT 匹配分析 (默认: TRUE)\n")
  cat("  --run_pipeline5b        [可选] 运行命运分化分析 (默认: TRUE)\n")
  quit(status = 1)
}

# 检查 RData 文件存在
if (!file.exists(PARAMS$rdata)) {
  stop(sprintf("RData 文件不存在: %s", PARAMS$rdata))
}

# 推断路径
paths <- infer_project_paths(PARAMS$rdata)

# 应用用户覆盖或使用推断值
BASE_DIR <- paths$base_dir
DESEQ2_DIR <- paths$deseq2_dir
OUTPUT_DIR <- if (!is.null(PARAMS$out_dir)) PARAMS$out_dir else paths$output_dir
PROJECT_NAME <- if (!is.null(PARAMS$project_name)) PARAMS$project_name else paths$project_name

# 条件标签 (v11.1: 完整的显示控制参数)
CONDITION_WT <- PARAMS$condition_wt
CONDITION_MUT <- PARAMS$condition_mut
CONDITION_WT_DISPLAY <- PARAMS$condition_wt_display
CONDITION_MUT_DISPLAY <- PARAMS$condition_mut_display
USE_ITALIC_WT <- PARAMS$use_italic_wt
USE_ITALIC_MUT <- PARAMS$use_italic_mut

################################################################################
# 创建输出目录
################################################################################

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
for (subdir in c("plots", "morphology_analysis", "rdata", 
                 "cluster_selection", "Pipeline5a_Figures", "Pipeline5b_Figures",
                 "AgriGO_Ready/Background", "AgriGO_Ready/Pipeline5a", 
                 "AgriGO_Ready/Pipeline5b")) {
  dir.create(file.path(OUTPUT_DIR, subdir), recursive = TRUE, showWarnings = FALSE)
}

################################################################################
# 日志设置
################################################################################

log_file <- file.path(OUTPUT_DIR, paste0(PROJECT_NAME, "_mfuzz_v11.log"))
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

################################################################################
# 颜色配置
################################################################################

UNIFIED_COLORS <- list(
  Conserved = "#98DF8A",
  Shape     = "#FF9896",
  Amp       = "#FFBB78",
  Phase     = "#F7B6D2"
)

PRIORITY_COLORS <- UNIFIED_COLORS

FATE_COLORS <- list(
  Conserved      = UNIFIED_COLORS$Conserved,
  Shape_Diverged = UNIFIED_COLORS$Shape,
  Amp_Altered    = UNIFIED_COLORS$Amp,
  Phase_Shifted  = UNIFIED_COLORS$Phase
)

COLOR_WT  <- "#377EB8"
COLOR_MUT <- "#984EA3"
COLOR_WT_RIBBON  <- adjustcolor(COLOR_WT, alpha.f = 0.2)
COLOR_MUT_RIBBON <- adjustcolor(COLOR_MUT, alpha.f = 0.2)

SUBCLUSTER_COLORS <- c(
  "#E74C3C", "#3498DB", "#2ECC71", "#F39C12", 
  "#9B59B6", "#1ABC9C", "#E67E22", "#34495E"
)

################################################################################
# 核心函数: DDTW正交三维指标 (v10.0 完整保留)
################################################################################

#' 计算数值导数
calc_derivative <- function(y, dt = 1, method = "central") {
  n <- length(y)
  dy <- numeric(n)
  
  if (method == "central") {
    for (i in 2:(n-1)) {
      dy[i] <- (y[i+1] - y[i-1]) / (2 * dt)
    }
    dy[1] <- (y[2] - y[1]) / dt
    dy[n] <- (y[n] - y[n-1]) / dt
    
  } else if (method == "forward") {
    for (i in 1:(n-1)) {
      dy[i] <- (y[i+1] - y[i]) / dt
    }
    dy[n] <- dy[n-1]
    
  } else if (method == "backward") {
    dy[1] <- (y[2] - y[1]) / dt
    for (i in 2:n) {
      dy[i] <- (y[i] - y[i-1]) / dt
    }
  }
  
  return(dy)
}

#' DTW算法实现
dtw_core <- function(x, y, window = NULL) {
  n <- length(x)
  m <- length(y)
  
  D <- matrix(Inf, n + 1, m + 1)
  D[1, 1] <- 0
  
  d <- matrix(0, n, m)
  for (i in 1:n) {
    for (j in 1:m) {
      d[i, j] <- (x[i] - y[j])^2
    }
  }
  
  for (i in 1:n) {
    if (!is.null(window)) {
      j_start <- max(1, round(i * m / n) - window)
      j_end <- min(m, round(i * m / n) + window)
    } else {
      j_start <- 1
      j_end <- m
    }
    
    for (j in j_start:j_end) {
      D[i + 1, j + 1] <- d[i, j] + min(
        D[i, j + 1],
        D[i + 1, j],
        D[i, j]
      )
    }
  }
  
  path_i <- c(n)
  path_j <- c(m)
  i <- n
  j <- m
  
  while (i > 1 || j > 1) {
    if (i == 1) {
      j <- j - 1
    } else if (j == 1) {
      i <- i - 1
    } else {
      candidates <- c(D[i, j + 1], D[i + 1, j], D[i, j])
      which_min <- which.min(candidates)
      
      if (which_min == 1) {
        i <- i - 1
      } else if (which_min == 2) {
        j <- j - 1
      } else {
        i <- i - 1
        j <- j - 1
      }
    }
    
    path_i <- c(i, path_i)
    path_j <- c(j, path_j)
  }
  
  dtw_distance <- sqrt(D[n + 1, m + 1])
  path_length <- length(path_i)
  normalized_distance <- dtw_distance / sqrt(path_length)
  
  return(list(
    distance = dtw_distance,
    normalized_distance = normalized_distance,
    path_i = path_i,
    path_j = path_j,
    path_length = path_length,
    cost_matrix = D[-1, -1],
    local_cost = d
  ))
}

#' 计算DTW路径偏离度（Phase指标）
calc_path_deviation <- function(path_i, path_j, n, m) {
  j_expected <- path_i * m / n
  deviations <- abs(path_j - j_expected)
  mean_deviation <- mean(deviations)
  max_possible_deviation <- max(m, n) / 2
  normalized_deviation <- min(mean_deviation / max_possible_deviation, 1)
  
  return(list(
    deviation = normalized_deviation,
    mean_raw_deviation = mean_deviation,
    max_deviation = max(deviations),
    path_length = length(path_i)
  ))
}

#' 计算基于DDTW的正交三维指标 (v10.0核心函数)
calc_ddtw_metrics <- function(center_A, center_B, 
                               t,
                               n_interp = COST_PARAMS$interp_points,
                               dtw_window = COST_PARAMS$dtw_window,
                               deriv_method = COST_PARAMS$derivative_method) {
  
  # Step 1: 插值到均匀细密网格
  t_range <- range(t)
  t_fine <- seq(t_range[1], t_range[2], length.out = n_interp)
  dt <- t_fine[2] - t_fine[1]
  
  spline_A <- splinefun(t, center_A, method = "natural")
  spline_B <- splinefun(t, center_B, method = "natural")
  
  z_A <- spline_A(t_fine)
  z_B <- spline_B(t_fine)
  
  # Step 2: 计算导数序列
  dz_A <- calc_derivative(z_A, dt, method = deriv_method)
  dz_B <- calc_derivative(z_B, dt, method = deriv_method)
  
  # Step 3: DDTW - 对导数序列进行DTW
  ddtw_result <- dtw_core(dz_A, dz_B, window = dtw_window)
  
  # Step 4: 计算Shape距离
  dz_range_A <- max(dz_A) - min(dz_A)
  dz_range_B <- max(dz_B) - min(dz_B)
  max_deriv_range <- max(dz_range_A, dz_range_B, 1e-6)
  
  shape_dist <- min(ddtw_result$normalized_distance / (max_deriv_range + 1e-6), 1)
  
  # Step 5: 计算Phase距离 (路径偏离度)
  path_dev <- calc_path_deviation(
    ddtw_result$path_i, 
    ddtw_result$path_j,
    n_interp, n_interp
  )
  phase_dist <- path_dev$deviation
  
  # Step 6: 计算Amplitude距离 (SDR)
  range_A <- max(center_A) - min(center_A)
  range_B <- max(center_B) - min(center_B)
  amp_dist <- abs(range_A - range_B) / (range_A + range_B + 1e-10)
  
  # Step 7: 计算Pearson相关
  pearson_r <- cor(center_A, center_B)
  deriv_r <- cor(dz_A, dz_B)
  
  # Step 8: 计算总成本
  cost <- COST_PARAMS$w_shape * shape_dist + 
          COST_PARAMS$w_amp * amp_dist + 
          COST_PARAMS$w_phase * phase_dist
  
  return(list(
    shape_dist = shape_dist,
    amp_dist = amp_dist,
    phase_dist = phase_dist,
    cost = cost,
    ddtw_distance = ddtw_result$distance,
    ddtw_normalized = ddtw_result$normalized_distance,
    path_length = ddtw_result$path_length,
    mean_path_deviation = path_dev$mean_raw_deviation,
    max_path_deviation = path_dev$max_deviation,
    range_A = range_A,
    range_B = range_B,
    pearson_r = pearson_r,
    deriv_r = deriv_r,
    dtw_path = list(i = ddtw_result$path_i, j = ddtw_result$path_j),
    derivatives = list(dz_A = dz_A, dz_B = dz_B),
    interpolated = list(z_A = z_A, z_B = z_B, t = t_fine)
  ))
}

#' 计算3-正交维度复合成本矩阵
calc_composite_cost_matrix <- function(wt_centers, mut_centers, time_points) {
  n_wt <- nrow(wt_centers)
  n_mut <- nrow(mut_centers)
  
  cost_matrix <- matrix(0, n_wt, n_mut)
  shape_dist_matrix <- matrix(0, n_wt, n_mut)
  amp_dist_matrix <- matrix(0, n_wt, n_mut)
  phase_dist_matrix <- matrix(0, n_wt, n_mut)
  pearson_r_matrix <- matrix(0, n_wt, n_mut)
  deriv_r_matrix <- matrix(0, n_wt, n_mut)
  
  for (i in 1:n_wt) {
    for (j in 1:n_mut) {
      metrics <- calc_ddtw_metrics(
        center_A = wt_centers[i, ],
        center_B = mut_centers[j, ],
        t = time_points
      )
      
      cost_matrix[i, j] <- metrics$cost
      shape_dist_matrix[i, j] <- metrics$shape_dist
      amp_dist_matrix[i, j] <- metrics$amp_dist
      phase_dist_matrix[i, j] <- metrics$phase_dist
      pearson_r_matrix[i, j] <- metrics$pearson_r
      deriv_r_matrix[i, j] <- metrics$deriv_r
    }
  }
  
  rownames(cost_matrix) <- paste0("WT_C", 1:n_wt)
  colnames(cost_matrix) <- paste0("MUT_C", 1:n_mut)
  dimnames(shape_dist_matrix) <- dimnames(cost_matrix)
  dimnames(amp_dist_matrix) <- dimnames(cost_matrix)
  dimnames(phase_dist_matrix) <- dimnames(cost_matrix)
  dimnames(pearson_r_matrix) <- dimnames(cost_matrix)
  dimnames(deriv_r_matrix) <- dimnames(cost_matrix)
  
  return(list(
    cost_matrix = cost_matrix,
    shape_dist_matrix = shape_dist_matrix,
    amp_dist_matrix = amp_dist_matrix,
    phase_dist_matrix = phase_dist_matrix,
    pearson_r_matrix = pearson_r_matrix,
    deriv_r_matrix = deriv_r_matrix,
    weights = list(
      w_shape = COST_PARAMS$w_shape,
      w_amp = COST_PARAMS$w_amp,
      w_phase = COST_PARAMS$w_phase
    )
  ))
}

################################################################################
# Schwämmle-Jensen公式估算模糊因子m
################################################################################

estimate_m_sj <- function(N, D, m_lower = 1.25, m_upper = 2.00, verbose = TRUE) {
  term1 <- (1418/N + 22.05) * D^(-2)
  term2 <- (12.33/N + 0.243) * D^(-0.0406 * log(N))
  term3 <- -0.1134
  m_raw <- 1 + term1 + term2 + term3
  
  m_constrained <- max(min(m_raw, m_upper), m_lower)
  truncated <- m_raw < m_lower || m_raw > m_upper
  
  if (verbose) {
    cat(sprintf("  Schwämmle-Jensen原始计算: m = %.4f\n", m_raw))
    if (truncated) {
      cat(sprintf("  ! 原始m=%.4f超出约束区间[%.2f, %.2f]，截断至m=%.4f\n", 
                  m_raw, m_lower, m_upper, m_constrained))
    } else {
      cat(sprintf("  ✔ m=%.4f 在约束区间[%.2f, %.2f]内\n", m_constrained, m_lower, m_upper))
    }
  }
  
  return(list(
    m = m_constrained,
    m_raw = m_raw,
    N = N,
    D = D,
    truncated = truncated,
    bounds = c(lower = m_lower, upper = m_upper)
  ))
}

################################################################################
# Xie-Beni Index计算
################################################################################

calc_xie_beni <- function(data, centers, membership, m) {
  N <- nrow(data)
  c <- nrow(centers)
  
  compactness <- 0
  for (i in 1:c) {
    for (j in 1:N) {
      dist_sq <- sum((data[j, ] - centers[i, ])^2)
      compactness <- compactness + (membership[j, i]^m) * dist_sq
    }
  }
  
  if (c < 2) return(Inf)
  
  min_center_dist_sq <- Inf
  for (i in 1:(c-1)) {
    for (k in (i+1):c) {
      dist_sq <- sum((centers[i, ] - centers[k, ])^2)
      if (dist_sq < min_center_dist_sq) {
        min_center_dist_sq <- dist_sq
      }
    }
  }
  
  if (min_center_dist_sq == 0) return(Inf)
  
  separation <- N * min_center_dist_sq
  xb <- compactness / separation
  
  return(xb)
}

################################################################################
# Kwon指数计算 (v11.0 新增 - 修正XB单调递减问题)
# 
# 理论基础：Kwon (1998) 提出的修正指数，通过添加惩罚项解决XB单调递减问题
# 惩罚项 = (1/c) * Σ||v_j - v̄||² 随c增加而增加，抵消紧致度的自然下降
################################################################################

calc_kwon_index <- function(data, centers, membership, m) {
  N <- nrow(data)
  c <- nrow(centers)
  
  # Term 1: 紧致度 (原XB分子)
  compactness <- 0
  for (i in 1:c) {
    for (j in 1:N) {
      dist_sq <- sum((data[j, ] - centers[i, ])^2)
      compactness <- compactness + (membership[j, i]^m) * dist_sq
    }
  }
  
  # Term 2: Kwon惩罚项 (质心到全局质心的距离)
  # 当c增加时，质心会向数据空间边缘扩散，使此项增大
  global_center <- colMeans(centers)
  penalty <- sum(apply(centers, 1, function(v) sum((v - global_center)^2))) / c
  
  # Term 3: 分离度 (最小质心间距，不乘N)
  if (c < 2) return(Inf)
  
  min_center_dist_sq <- Inf
  for (i in 1:(c-1)) {
    for (k in (i+1):c) {
      dist_sq <- sum((centers[i, ] - centers[k, ])^2)
      if (dist_sq < min_center_dist_sq) {
        min_center_dist_sq <- dist_sq
      }
    }
  }
  
  if (min_center_dist_sq == 0) return(Inf)
  
  # Kwon指数（注意：分母不乘N，与XB不同）
  kwon <- (compactness + penalty) / min_center_dist_sq
  
  return(kwon)
}

################################################################################
# 样本量约束函数 (v11.0 新增 - 基于√N规则)
#
# 理论基础：
# - √N规则：Schwämmle & Jensen (2010) 确认 c ≤ √N 作为理论上限
# - N/10规则：每聚类至少需要10个样本才能可靠估计质心
# - Pipeline5b更保守：小样本情况使用 √N/2
################################################################################

apply_sample_size_constraint <- function(c_range, N, pipeline_type = "5a") {
  # 理论上限
  c_max_sqrt <- floor(sqrt(N))
  c_max_ratio <- floor(N / 10)  # 每聚类至少10个样本
  
  if (pipeline_type == "5b") {
    # Pipeline5b小样本更保守
    if (N < 50) {
      c_max_theoretical <- min(4, c_max_sqrt, c_max_ratio)
    } else if (N < 100) {
      c_max_theoretical <- min(floor(sqrt(N) / 1.5), c_max_ratio, 6)
    } else {
      c_max_theoretical <- min(c_max_sqrt, c_max_ratio)
    }
  } else {
    # Pipeline5a
    c_max_theoretical <- min(c_max_sqrt, c_max_ratio)
  }
  
  # 应用约束，确保至少保留用户设定的最小值
  c_min_user <- min(c_range)
  c_max_constrained <- max(c_max_theoretical, c_min_user + 2)  # 至少3个点
  
  # 过滤c_range
  c_range_constrained <- c_range[c_range <= c_max_constrained]
  
  # 确保至少有3个c值用于曲线分析
  if (length(c_range_constrained) < 3) {
    c_range_constrained <- c_min_user:min(c_min_user + 2, max(c_range))
    warning(sprintf("样本量N=%d过小，c_range被约束为[%d, %d]", 
                    N, min(c_range_constrained), max(c_range_constrained)))
  }
  
  return(list(
    c_range = c_range_constrained,
    c_max_theoretical = c_max_theoretical,
    c_max_sqrt = c_max_sqrt,
    c_max_ratio = c_max_ratio,
    was_constrained = max(c_range) > c_max_constrained
  ))
}

################################################################################
# 离散曲率计算 (v11.0 新增 - 用于单调递减曲线的停止点检测)
#
# 基于L-method思想：曲率最大点是最佳停止点
################################################################################

calc_discrete_curvature <- function(y) {
  n <- length(y)
  if (n < 3) return(rep(0, n))
  
  curvature <- rep(0, n)
  
  for (i in 2:(n-1)) {
    # 一阶差分（前后差分的平均）
    d1 <- y[i] - y[i-1]
    d2 <- y[i+1] - y[i]
    
    # 二阶差分近似曲率
    # 使用标准曲率公式: κ = |y''| / (1 + y'^2)^(3/2)
    y_prime <- (d1 + d2) / 2
    y_double_prime <- d2 - d1
    
    curvature[i] <- abs(y_double_prime) / (1 + y_prime^2)^1.5
  }
  
  return(curvature)
}

################################################################################
# XB-Based聚类数选择算法 (v11.0 重构 - PP-XB算法)
################################################################################

find_valleys_with_prominence <- function(x) {
  n <- length(x)
  data_range <- max(x) - min(x)
  if (data_range == 0) data_range <- 1
  
  valley_indices <- c()
  for (i in 2:(n-1)) {
    if (x[i] < x[i-1] && x[i] <= x[i+1]) {
      valley_indices <- c(valley_indices, i)
    }
  }
  
  prominences <- numeric(length(valley_indices))
  for (j in seq_along(valley_indices)) {
    valley_idx <- valley_indices[j]
    valley_val <- x[valley_idx]
    
    left_max <- valley_val
    if (valley_idx > 1) {
      for (i in (valley_idx - 1):1) {
        if (x[i] < valley_val) break
        left_max <- max(left_max, x[i])
      }
    }
    
    right_max <- valley_val
    if (valley_idx < n) {
      for (i in (valley_idx + 1):n) {
        if (x[i] < valley_val) break
        right_max <- max(right_max, x[i])
      }
    }
    
    base <- min(left_max, right_max)
    prominences[j] <- base - valley_val
  }
  
  rel_prominences <- prominences / data_range
  
  return(list(
    valley_indices = valley_indices,
    prominences = prominences,
    rel_prominences = rel_prominences,
    data_range = data_range
  ))
}

#' 曲线形态分类 (v11.0 重构版)
#' 
#' 改进要点：
#' 1. 增加U形专门检测（首尾比中间高）
#' 2. 分离显著谷底检测与单调性检测
#' 3. 避免PAVA破坏U形信息
#' 
#' @param values 指标值向量（Kwon或XB）
#' @param c_values 对应的c值向量
#' @param prominence_threshold 显著性阈值（默认0.10）
#' @return list 包含曲线类型、原因和推荐策略
classify_curve_morphology <- function(values, c_values, prominence_threshold = 0.10) {
  n <- length(values)
  valley_info <- find_valleys_with_prominence(values)
  
  # ========== 检测1: U形曲线（首尾比中间高）==========
  # 将曲线分为三段，检查首尾是否显著高于中间
  first_third_idx <- 1:max(2, n %/% 3)
  last_third_idx <- max(1, n - n %/% 3 + 1):n
  middle_idx <- max(2, n %/% 4):min(n-1, 3*n %/% 4)
  
  first_third_mean <- mean(values[first_third_idx])
  last_third_mean <- mean(values[last_third_idx])
  middle_min <- min(values[middle_idx])
  middle_mean <- mean(values[middle_idx])
  
  # U形判据：首尾均值都比中间最小值高至少5%
  is_u_shaped <- (first_third_mean > middle_min * 1.05) && 
                 (last_third_mean > middle_min * 1.05)
  
  # ========== 检测2: 显著谷底 ==========
  has_significant_valley <- length(valley_info$valley_indices) > 0 && 
                            max(valley_info$rel_prominences, na.rm = TRUE) >= prominence_threshold
  
  # ========== 检测3: 单调性 ==========
  n_decreasing <- sum(diff(values) < 0)
  monotonicity_ratio <- n_decreasing / (n - 1)
  is_strongly_monotonic <- monotonicity_ratio > 0.85
  
  # ========== 分类决策 ==========
  if (is_u_shaped || has_significant_valley) {
    # 类型1: U形或有显著谷底 → 直接选择全局最小值
    reason_parts <- c()
    if (is_u_shaped) reason_parts <- c(reason_parts, "U_shape_detected")
    if (has_significant_valley) {
      max_prom <- max(valley_info$rel_prominences, na.rm = TRUE)
      reason_parts <- c(reason_parts, sprintf("significant_valley(%.1f%%)", max_prom * 100))
    }
    
    return(list(
      curve_type = "U_SHAPED_OR_VALLEY",
      reason = paste(reason_parts, collapse = " + "),
      strategy = "select_global_minimum",
      is_u_shaped = is_u_shaped,
      has_significant_valley = has_significant_valley,
      max_rel_prominence = ifelse(length(valley_info$rel_prominences) > 0, 
                                   max(valley_info$rel_prominences), 0),
      valley_info = valley_info,
      monotonicity_ratio = monotonicity_ratio
    ))
  } else if (is_strongly_monotonic) {
    # 类型2: 严格单调递减 → 使用曲率最大点
    return(list(
      curve_type = "MONOTONIC_DECREASING",
      reason = sprintf("monotonicity=%.1f%%", monotonicity_ratio * 100),
      strategy = "max_curvature_or_elbow",
      is_u_shaped = FALSE,
      has_significant_valley = FALSE,
      max_rel_prominence = ifelse(length(valley_info$rel_prominences) > 0, 
                                   max(valley_info$rel_prominences), 0),
      valley_info = valley_info,
      monotonicity_ratio = monotonicity_ratio
    ))
  } else {
    # 类型3: 波动型 → 平滑后重新分析
    return(list(
      curve_type = "FLUCTUATING",
      reason = sprintf("no_clear_pattern(monotonicity=%.1f%%, max_prom=%.1f%%)", 
                       monotonicity_ratio * 100,
                       ifelse(length(valley_info$rel_prominences) > 0, 
                              max(valley_info$rel_prominences) * 100, 0)),
      strategy = "smooth_then_select",
      is_u_shaped = FALSE,
      has_significant_valley = FALSE,
      max_rel_prominence = ifelse(length(valley_info$rel_prominences) > 0, 
                                   max(valley_info$rel_prominences), 0),
      valley_info = valley_info,
      monotonicity_ratio = monotonicity_ratio
    ))
  }
}

# 保留旧函数名以兼容，但内部调用新函数
classify_xb_curve_shape <- function(xb_values, threshold = 0.03) {
  # 为了向后兼容，保留原函数签名
  result <- classify_curve_morphology(xb_values, seq_along(xb_values), threshold)
  
  # 转换为旧格式
  if (result$curve_type == "U_SHAPED_OR_VALLEY") {
    return(list(
      curve_type = "NON_CONVEX",
      reason = result$reason,
      max_rel_prominence = result$max_rel_prominence,
      valley_info = result$valley_info
    ))
  } else {
    return(list(
      curve_type = "QUASI_MONOTONIC",
      reason = result$reason,
      max_rel_prominence = result$max_rel_prominence,
      valley_info = result$valley_info
    ))
  }
}

pava_decreasing <- function(x) {
  n <- length(x)
  if (n <= 1) return(x)
  
  y <- -x
  result <- y
  
  i <- 1
  while (i < n) {
    if (result[i] > result[i + 1]) {
      j <- i + 1
      pool_sum <- result[i] + result[j]
      pool_count <- 2
      
      while (j < n && (pool_sum / pool_count) > result[j + 1]) {
        j <- j + 1
        pool_sum <- pool_sum + result[j]
        pool_count <- pool_count + 1
      }
      
      pool_mean <- pool_sum / pool_count
      for (k in i:j) {
        result[k] <- pool_mean
      }
      
      while (i > 1 && result[i - 1] > result[i]) {
        i_start <- i - 1
        while (i_start > 1 && result[i_start - 1] == result[i_start]) {
          i_start <- i_start - 1
        }
        
        pool_sum <- sum(result[i_start:j])
        pool_count <- j - i_start + 1
        pool_mean <- pool_sum / pool_count
        
        for (k in i_start:j) {
          result[k] <- pool_mean
        }
        
        i <- i_start
      }
      
      i <- j + 1
    } else {
      i <- i + 1
    }
  }
  
  return(-result)
}

kneedle_detect <- function(x, y, S = 1.0) {
  n <- length(x)
  
  if (n < 3) {
    return(list(
      knee_idx = 1,
      knee_x = x[1],
      knee_y = y[1],
      differences = rep(0, n),
      at_boundary = TRUE
    ))
  }
  
  x_range <- max(x) - min(x)
  y_range <- max(y) - min(y)
  if (x_range == 0) x_range <- 1
  if (y_range == 0) y_range <- 1
  
  x_norm <- (x - min(x)) / x_range
  y_norm <- (y - min(y)) / y_range
  
  differences <- y_norm + x_norm - 1
  abs_diff <- abs(differences)
  
  interior_indices <- 2:(n-1)
  if (length(interior_indices) == 0) interior_indices <- 1:n
  
  knee_idx <- interior_indices[which.max(abs_diff[interior_indices])]
  at_boundary <- (knee_idx == 1 || knee_idx == n)
  
  return(list(
    knee_idx = knee_idx,
    knee_x = x[knee_idx],
    knee_y = y[knee_idx],
    differences = differences,
    at_boundary = at_boundary
  ))
}

# 保留旧函数以兼容
find_optimal_c_xb_minimum <- function(c_values, xb_values) {
  n <- length(c_values)
  
  valley_indices <- c()
  for (i in 2:(n-1)) {
    if (xb_values[i] < xb_values[i-1] && xb_values[i] <= xb_values[i+1]) {
      valley_indices <- c(valley_indices, i)
    }
  }
  
  if (length(valley_indices) > 0) {
    candidate_xbs <- xb_values[valley_indices]
    best_valley_idx <- valley_indices[which.min(candidate_xbs)]
    optimal_idx <- best_valley_idx
    selection_method <- "valley_xb_minimum"
  } else {
    interior_indices <- 2:(n-1)
    if (length(interior_indices) > 0) {
      optimal_idx <- interior_indices[which.min(xb_values[interior_indices])]
      selection_method <- "interior_xb_minimum"
    } else {
      optimal_idx <- which.min(xb_values)
      selection_method <- "global_xb_minimum"
    }
  }
  
  optimal_c <- c_values[optimal_idx]
  at_boundary <- (optimal_idx == 1 || optimal_idx == n)
  
  return(list(
    optimal_c = optimal_c,
    optimal_idx = optimal_idx,
    optimal_xb = xb_values[optimal_idx],
    at_boundary = isTRUE(at_boundary),
    valley_indices = valley_indices,
    selection_method = selection_method
  ))
}

find_optimal_c_kneedle_xb <- function(c_values, xb_values, S = 1.0) {
  xb_smoothed <- pava_decreasing(xb_values)
  kneedle_result <- kneedle_detect(c_values, xb_smoothed, S = S)
  
  optimal_c <- kneedle_result$knee_x
  optimal_idx <- kneedle_result$knee_idx
  at_boundary <- isTRUE(kneedle_result$at_boundary)
  
  return(list(
    optimal_c = optimal_c,
    optimal_idx = optimal_idx,
    optimal_xb = xb_smoothed[optimal_idx],
    xb_smoothed = xb_smoothed,
    kneedle_result = kneedle_result,
    at_boundary = at_boundary,
    selection_method = "pava_kneedle"
  ))
}

# 保留旧接口以兼容
find_optimal_c_xb_based <- function(c_values, xb_values,
                                     prominence_threshold = 0.03,
                                     kneedle_S = 1.0) {
  curve_class <- classify_xb_curve_shape(xb_values, threshold = prominence_threshold)
  
  if (curve_class$curve_type == "NON_CONVEX") {
    selection_result <- find_optimal_c_xb_minimum(c_values, xb_values)
    selection_result$curve_class <- curve_class
    selection_result$strategy <- "NON_CONVEX -> valley_xb_minimum"
  } else {
    selection_result <- find_optimal_c_kneedle_xb(c_values, xb_values, S = kneedle_S)
    selection_result$curve_class <- curve_class
    selection_result$strategy <- "QUASI_MONOTONIC -> pava_kneedle"
  }
  
  return(selection_result)
}

################################################################################
# PP-XB统一算法框架 (v11.0 新增)
#
# 核心改进：
# 1. 使用Kwon指数替代XB（解决单调递减）
# 2. 基于曲线形态分类的自适应选择策略
# 3. 曲率最大点替代PAVA+Kneedle（处理单调递减）
# 4. √N样本量约束
################################################################################

#' 统一的最优聚类数选择函数 (PP-XB算法)
#' 
#' @param eset ExpressionSet对象
#' @param c_range 候选c值范围
#' @param m 模糊因子
#' @param pipeline_type "5a" (全基因组) 或 "5b" (子聚类)
#' @param prominence_threshold 显著性阈值 (Pipeline5a: 0.10, Pipeline5b: 0.15)
#' @param smoothing_span LOWESS平滑带宽 (Pipeline5a: 0.3, Pipeline5b: 0.5)
#' @return list(optimal_c, kwon_values, xb_values, diagnostics)
find_optimal_c_unified <- function(eset, c_range, m, 
                                    pipeline_type = "5a",
                                    prominence_threshold = NULL,
                                    smoothing_span = NULL) {
  
  N <- nrow(exprs(eset))
  D <- ncol(exprs(eset))
  data_matrix <- exprs(eset)
  
  # Set Pipeline-specific parameters
  if (is.null(prominence_threshold)) {
    prominence_threshold <- ifelse(pipeline_type == "5b", 0.15, 0.10)
  }
  if (is.null(smoothing_span)) {
    smoothing_span <- ifelse(pipeline_type == "5b", 0.5, 0.3)
  }
  
  # ========== Phase 0: Sample Size Constraint ==========
  constraint_result <- apply_sample_size_constraint(c_range, N, pipeline_type)
  c_range_constrained <- constraint_result$c_range
  
  cat(sprintf("  [Sample Constraint] N=%d, sqrt(N)=%.1f, c_max_theoretical=%d\n", 
              N, sqrt(N), constraint_result$c_max_theoretical))
  if (constraint_result$was_constrained) {
    cat(sprintf("    Original c_range=[%d,%d] -> Constrained=[%d,%d]\n",
                min(c_range), max(c_range),
                min(c_range_constrained), max(c_range_constrained)))
  }
  
  # ========== Phase 1: Kwon and XB Index Scan ==========
  kwon_values <- numeric(length(c_range_constrained))
  xb_values <- numeric(length(c_range_constrained))
  empty_cluster_detected <- FALSE
  
  cat(sprintf("\n  [Index Scan] c in [%d, %d]\n", 
              min(c_range_constrained), max(c_range_constrained)))
  
  for (i in seq_along(c_range_constrained)) {
    c_val <- c_range_constrained[i]
    set.seed(123)
    cl <- mfuzz(eset, c = c_val, m = m)
    
    # Empty cluster detection
    max_memberships <- apply(cl$membership, 2, max)
    if (any(max_memberships < 0.5)) {
      cat(sprintf("    c=%2d: Warning - Empty cluster detected, stopping scan\n", c_val))
      empty_cluster_detected <- TRUE
      # Truncate data
      if (i > 1) {
        c_range_constrained <- c_range_constrained[1:(i-1)]
        kwon_values <- kwon_values[1:(i-1)]
        xb_values <- xb_values[1:(i-1)]
      }
      break
    }
    
    kwon_values[i] <- calc_kwon_index(data_matrix, cl$centers, cl$membership, m)
    xb_values[i] <- calc_xie_beni(data_matrix, cl$centers, cl$membership, m)
    
    cat(sprintf("    c=%2d: Kwon=%.4f, XB=%.6f\n", c_val, kwon_values[i], xb_values[i]))
  }
  
  # Check if enough points for analysis
  if (length(c_range_constrained) < 3) {
    warning("c_range too small for curve analysis, returning midpoint")
    mid_idx <- ceiling(length(c_range_constrained) / 2)
    return(list(
      optimal_c = c_range_constrained[mid_idx],
      optimal_idx = mid_idx,
      optimal_kwon = kwon_values[mid_idx],
      kwon_values = kwon_values,
      xb_values = xb_values,
      c_range = c_range_constrained,
      curve_class = list(curve_type = "INSUFFICIENT_DATA", reason = "too_few_points"),
      selection_method = "fallback_midpoint",
      at_boundary = FALSE,
      N = N,
      D = D,
      empty_cluster_detected = empty_cluster_detected
    ))
  }
  
  # ========== Phase 2: Curve Morphology Classification (based on Kwon) ==========
  curve_class <- classify_curve_morphology(kwon_values, c_range_constrained, prominence_threshold)
  cat(sprintf("\n  [Curve Classification] %s\n", curve_class$curve_type))
  cat(sprintf("    Reason: %s\n", curve_class$reason))
  cat(sprintf("    Strategy: %s\n", curve_class$strategy))
  
  # ========== Phase 3: Morphology-Specific Selection ==========
  if (curve_class$curve_type == "U_SHAPED_OR_VALLEY") {
    # Strategy 1: U-shaped or significant valley -> select global minimum
    optimal_idx <- which.min(kwon_values)
    selection_method <- "global_kwon_minimum"
    
  } else if (curve_class$curve_type == "MONOTONIC_DECREASING") {
    # Strategy 2: Monotonic decreasing -> use max curvature point
    curvature <- calc_discrete_curvature(kwon_values)
    
    # Find max curvature in interior points
    interior_indices <- 2:(length(kwon_values)-1)
    if (length(interior_indices) > 0) {
      interior_curvatures <- curvature[interior_indices]
      max_curv_interior_idx <- which.max(interior_curvatures)
      optimal_idx <- interior_indices[max_curv_interior_idx]
      selection_method <- "max_curvature"
      
      # If curvature peak not prominent, use conservative strategy
      if (max(curvature, na.rm = TRUE) < 0.01) {
        optimal_idx <- max(2, ceiling(length(c_range_constrained) / 3))
        selection_method <- "conservative_fallback_one_third"
        cat(sprintf("    Warning: Curvature peak not prominent, using conservative strategy\n"))
      }
    } else {
      optimal_idx <- 1
      selection_method <- "fallback_first"
    }
    
  } else {  # FLUCTUATING
    # Strategy 3: Fluctuating -> LOWESS smoothing then select minimum
    if (length(kwon_values) >= 4) {
      kwon_smoothed <- lowess(c_range_constrained, kwon_values, f = smoothing_span)$y
      optimal_idx <- which.min(kwon_smoothed)
      selection_method <- "smoothed_minimum"
    } else {
      optimal_idx <- which.min(kwon_values)
      selection_method <- "global_kwon_minimum_fallback"
    }
  }
  
  optimal_c <- c_range_constrained[optimal_idx]
  
  # ========== Phase 4: Post-hoc Validation ==========
  at_boundary <- (optimal_idx == 1 || optimal_idx == length(c_range_constrained))
  
  # Prominence validation (for valley selection)
  prominence_valid <- NA
  if (curve_class$curve_type == "U_SHAPED_OR_VALLEY") {
    valley_info <- curve_class$valley_info
    if (!is.null(valley_info) && optimal_idx %in% valley_info$valley_indices) {
      prom_idx <- which(valley_info$valley_indices == optimal_idx)
      prominence_valid <- valley_info$rel_prominences[prom_idx] >= prominence_threshold
    }
  }
  
  # ========== Output Diagnostic Info ==========
  cat(sprintf("\n  [Selection Result] c*=%d (Kwon=%.4f, XB=%.6f) [%s]\n", 
              optimal_c, kwon_values[optimal_idx], xb_values[optimal_idx], selection_method))
  
  if (at_boundary) {
    cat("    Warning: Optimal c at boundary, consider adjusting c_range\n")
  }
  
  return(list(
    optimal_c = optimal_c,
    optimal_idx = optimal_idx,
    optimal_kwon = kwon_values[optimal_idx],
    optimal_xb = xb_values[optimal_idx],
    kwon_values = kwon_values,
    xb_values = xb_values,
    c_range = c_range_constrained,
    curve_class = curve_class,
    selection_method = selection_method,
    at_boundary = at_boundary,
    prominence_valid = prominence_valid,
    N = N,
    D = D,
    empty_cluster_detected = empty_cluster_detected,
    constraint_result = constraint_result
  ))
}

################################################################################
# XB-Based聚类数选择算法 (v11.0 重构版 - 双面板诊断图)
################################################################################

#' 扫描最优聚类数 (v11.0 PP-XB算法)
#' 
#' 改进要点：
#' 1. 使用find_optimal_c_unified统一框架
#' 2. 生成双面板诊断图（原始XB + Kwon指数）
#' 3. 始终报告原始指标值
#' 
#' @param eset ExpressionSet对象
#' @param c_range 候选c值范围
#' @param m 模糊因子
#' @param condition_name 条件名称
#' @param output_dir 输出目录
#' @param pipeline_type "5a" 或 "5b"
#' @param prominence_threshold 显著性阈值
#' @param save_plot 是否保存图形
scan_optimal_c <- function(eset, c_range, m, condition_name, output_dir,
                            pipeline_type = "5a",
                            prominence_threshold = NULL,
                            kneedle_S = XB_PARAMS$kneedle_S,
                            save_plot = TRUE) {
  
  # 使用默认阈值（如果未指定）
  if (is.null(prominence_threshold)) {
    if (pipeline_type == "5b") {
      prominence_threshold <- XB_PARAMS$subcluster_prominence_threshold
    } else {
      prominence_threshold <- XB_PARAMS$global_prominence_threshold
    }
  }
  
  cat(sprintf("\n【%s PP-XB最优聚类数选择】(v11.0)\n", condition_name))
  cat(sprintf("  Pipeline类型: %s\n", pipeline_type))
  cat(sprintf("  显著性阈值: %.1f%%\n", prominence_threshold * 100))
  
  # 调用统一框架
  result <- find_optimal_c_unified(
    eset = eset,
    c_range = c_range,
    m = m,
    pipeline_type = pipeline_type,
    prominence_threshold = prominence_threshold
  )
  
  # 生成诊断图
  if (save_plot && length(result$c_range) >= 3) {
    plot_filename <- paste0(condition_name, "_PP-XB_Diagnostic.pdf")
    plot_dual_diagnostic(
      c_range = result$c_range,
      xb_values = result$xb_values,
      kwon_values = result$kwon_values,
      optimal_c = result$optimal_c,
      optimal_idx = result$optimal_idx,
      curve_class = result$curve_class,
      selection_method = result$selection_method,
      condition_name = condition_name,
      output_path = file.path(output_dir, "cluster_selection", plot_filename)
    )
    cat(sprintf("  ✔ 诊断图已保存: %s\n", plot_filename))
  }
  
  # 转换为兼容旧接口的格式
  xb_result <- list(
    optimal_c = result$optimal_c,
    optimal_idx = result$optimal_idx,
    optimal_xb = result$optimal_xb,
    optimal_kwon = result$optimal_kwon,
    curve_class = result$curve_class,
    selection_method = result$selection_method,
    at_boundary = result$at_boundary,
    valley_indices = if (!is.null(result$curve_class$valley_info)) 
                       result$curve_class$valley_info$valley_indices else NULL,
    strategy = sprintf("%s -> %s", result$curve_class$curve_type, result$selection_method)
  )
  
  cat(sprintf("\n  ✔ 最优c=%d (Kwon=%.4f, XB=%.6f) [%s]\n",
              result$optimal_c,
              result$optimal_kwon,
              result$optimal_xb,
              result$selection_method))
  
  if (isTRUE(result$at_boundary)) {
    cat(sprintf("  ⚠ 警告: 最优c在边界，建议调整c_range范围\n"))
  }
  
  return(list(
    optimal_c = result$optimal_c,
    xb_values = result$xb_values,
    kwon_values = result$kwon_values,
    xb_result = xb_result,
    unified_result = result  # 保留完整结果用于诊断
  ))
}

#' 双面板诊断图 (v11.0 新增)
#' 
#' 同时显示原始XB曲线和修正后的Kwon曲线
#' 双面板诊断图 (v11.0 - English Labels)
#' 
plot_dual_diagnostic <- function(c_range, xb_values, kwon_values, 
                                  optimal_c, optimal_idx,
                                  curve_class, selection_method,
                                  condition_name, output_path) {
  
  pdf(output_path, width = 8, height = 8)
  par(mfrow = c(2, 1), mar = c(4, 5, 3, 2))
  
  # ========== Panel 1: Original XB Index ==========
  plot(c_range, xb_values, type = "b", pch = 19, col = "#3498DB",
       xlab = "", ylab = "XB Index (Raw)",
       main = sprintf("%s - Raw Xie-Beni Index", condition_name),
       cex.main = 1.1, cex.lab = 1.0, cex.axis = 0.9)
  
  # Mark optimal point on XB curve

points(optimal_c, xb_values[optimal_idx], pch = 1, col = "#E74C3C", cex = 2, lwd = 2)
  abline(v = optimal_c, lty = 2, col = "#E74C3C", lwd = 1)
  
  # Add annotation
  mtext("(May show monotonic decreasing trend)", side = 3, line = 0.2, cex = 0.8, col = "#7F8C8D")
  
  # ========== Panel 2: Kwon Index ==========
  plot(c_range, kwon_values, type = "b", pch = 17, col = "#27AE60",
       xlab = "Number of Clusters (c)", ylab = "Kwon Index (Corrected)",
       main = sprintf("%s - Kwon Index (%s)", condition_name, curve_class$curve_type),
       cex.main = 1.1, cex.lab = 1.0, cex.axis = 0.9)
  
  # Mark all detected valleys
  if (!is.null(curve_class$valley_info) && length(curve_class$valley_info$valley_indices) > 0) {
    valley_idx <- curve_class$valley_info$valley_indices
    points(c_range[valley_idx], kwon_values[valley_idx],
           pch = 1, col = "#F39C12", cex = 2, lwd = 2)
  }
  
  # Mark optimal point
  points(optimal_c, kwon_values[optimal_idx], pch = 19, col = "#E74C3C", cex = 2.5)
  abline(v = optimal_c, lty = 2, col = "#E74C3C", lwd = 1.5)
  
  # Add legend
  legend_items <- c(sprintf("Optimal c=%d", optimal_c), "Kwon Curve")
  legend_cols <- c("#E74C3C", "#27AE60")
  legend_pch <- c(19, 17)
  
  if (!is.null(curve_class$valley_info) && length(curve_class$valley_info$valley_indices) > 0) {
    legend_items <- c(legend_items, "Detected Valleys")
    legend_cols <- c(legend_cols, "#F39C12")
    legend_pch <- c(legend_pch, 1)
  }
  
  legend("topright", legend = legend_items, col = legend_cols, pch = legend_pch,
         cex = 0.85, bg = "white")
  
  # Bottom annotation
  mtext(sprintf("Strategy: %s | Curve Type: %s", selection_method, curve_class$curve_type), 
        side = 1, line = 2.8, cex = 0.75, col = "#7F8C8D")
  
  dev.off()
}

#' 准备单条件时间序列数据 (v10.1 - English Output)
#' 
#' @param vst_matrix VST转换后的表达矩阵
#' @param metadata 元数据
#' @param condition 条件标签
#' @param gene_list 基因列表
#' @param time_info 时间点信息 (从 auto_detect_timepoints 返回)
#' @param min_std 最小标准差阈值
prepare_single_timeseries <- function(vst_matrix, metadata, condition, gene_list,
                                       time_info, min_std) {
  
  cat(sprintf("\n[Processing Condition: %s]\n", condition))
  
  # ========== 1. Find condition column ==========
  cond_col <- NULL
  for (col_name in c("Condition", "condition", "Genotype", "genotype")) {
    if (col_name %in% colnames(metadata)) {
      cond_col <- col_name
      break
    }
  }
  if (is.null(cond_col)) {
    stop("Metadata missing 'Condition' or 'Genotype' column")
  }
  
  # ========== 2. Filter condition samples ==========
  cond_samples <- rownames(metadata)[as.character(metadata[[cond_col]]) == condition]
  
  if (length(cond_samples) == 0) {
    warning(sprintf("No samples found for condition '%s'", condition))
    return(list(matrix = matrix(nrow=0, ncol=time_info$D), filtered_genes = character(0)))
  }
  cat(sprintf("  Found %d samples for condition '%s'\n", length(cond_samples), condition))

  # ========== 3. Extract expression data subset ==========
  expr_sub <- vst_matrix[, cond_samples, drop = FALSE]
  
  # ========== 4. Initialize mean matrix ==========
  mean_matrix <- matrix(NA, nrow = length(gene_list), ncol = time_info$D)
  rownames(mean_matrix) <- gene_list
  colnames(mean_matrix) <- time_info$time_labels
  
  # ========== 5. Get sample time info ==========
  time_col <- time_info$time_col
  sample_time_labels <- as.character(metadata[cond_samples, time_col])
  
  # Convert sample time labels to numeric
  sample_time_numeric <- sapply(sample_time_labels, extract_time_numeric)
  
  # ========== 6. Debug output ==========
  cat(sprintf("  Time column: %s\n", time_col))
  cat(sprintf("  Sample time labels (unique): %s\n", 
              paste(unique(sample_time_labels), collapse = ", ")))
  cat(sprintf("  Sample time numeric (unique): %s\n", 
              paste(sort(unique(sample_time_numeric)), collapse = ", ")))
  cat(sprintf("  Target time points (numeric): %s\n", 
              paste(time_info$time_points, collapse = ", ")))
  
  # ========== 7. Calculate mean for each time point ==========
  matched_samples_total <- 0
  
  for (i in seq_along(time_info$time_points)) {
    target_time <- time_info$time_points[i]
    
    # Numeric matching (tolerance 0.1)
    tp_mask <- abs(sample_time_numeric - target_time) < 0.1
    tp_samples <- cond_samples[tp_mask]
    
    if (length(tp_samples) > 0) {
      available_genes <- intersect(gene_list, rownames(vst_matrix))
      if (length(available_genes) > 0) {
        mean_matrix[available_genes, i] <- rowMeans(
          expr_sub[available_genes, tp_samples, drop = FALSE], 
          na.rm = TRUE
        )
      }
      matched_samples_total <- matched_samples_total + length(tp_samples)
      cat(sprintf("    T%g: %d samples matched\n", target_time, length(tp_samples)))
    } else {
      cat(sprintf("    T%g: Warning - No samples matched!\n", target_time))
    }
  }
  
  # ========== 8. Check matching completeness ==========
  if (matched_samples_total == 0) {
    cat("  CRITICAL: No samples matched to any time point!\n")
    cat("    Please check time format in metadata.\n")
    return(list(matrix = matrix(nrow=0, ncol=time_info$D), filtered_genes = character(0)))
  }
  
  if (matched_samples_total < length(cond_samples)) {
    cat(sprintf("  Warning: Only %d/%d samples matched to time points\n", 
                matched_samples_total, length(cond_samples)))
  }
  
  # ========== 9. Filter complete data ==========
  complete_genes_idx <- complete.cases(mean_matrix)
  n_complete <- sum(complete_genes_idx)
  
  if (n_complete == 0) {
    cat("  CRITICAL: No genes have data at all time points!\n")
    return(list(matrix = matrix(nrow=0, ncol=time_info$D), filtered_genes = character(0)))
  }
  
  mean_matrix <- mean_matrix[complete_genes_idx, , drop = FALSE]
  cat(sprintf("  Complete data genes: %d / %d (%.1f%%)\n", 
              n_complete, length(gene_list), 100 * n_complete / length(gene_list)))
  
  # ========== 10. Standard deviation filter ==========
  if (nrow(mean_matrix) > 0) {
    gene_sds <- apply(mean_matrix, 1, sd, na.rm = TRUE)
    filtered_genes <- names(which(gene_sds >= min_std & !is.na(gene_sds)))
    result_matrix <- mean_matrix[filtered_genes, , drop = FALSE]
    cat(sprintf("  Std filter (std >= %.2f): %d -> %d genes\n", 
                min_std, n_complete, length(filtered_genes)))
  } else {
    filtered_genes <- character(0)
    result_matrix <- mean_matrix
  }
  
  return(list(
    matrix = result_matrix,
    filtered_genes = filtered_genes
  ))
}

################################################################################
# 匈牙利算法求解
################################################################################

solve_asymmetric_LSAP <- function(cost_matrix) {
  n_wt <- nrow(cost_matrix)
  n_mut <- ncol(cost_matrix)
  
  if (n_wt <= n_mut) {
    cat(sprintf("  执行正向匹配 (N_WT=%d <= N_MUT=%d)\n", n_wt, n_mut))
    
    solution <- solve_LSAP(cost_matrix, maximum = FALSE)
    optimal_matching <- as.integer(solution)
    
    matching_type <- "forward"
    unmatched_mut <- setdiff(1:n_mut, optimal_matching)
    
    result <- list(
      matching = optimal_matching,
      matching_type = matching_type,
      n_matched = n_wt,
      unmatched_wt = integer(0),
      unmatched_mut = unmatched_mut,
      base_condition = "WT"
    )
    
  } else {
    cat(sprintf("  执行转置匹配 (N_WT=%d > N_MUT=%d)\n", n_wt, n_mut))
    
    transposed_matrix <- t(cost_matrix)
    solution <- solve_LSAP(transposed_matrix, maximum = FALSE)
    mut_to_wt_matching <- as.integer(solution)
    
    optimal_matching <- rep(NA_integer_, n_wt)
    for (j in 1:n_mut) {
      wt_idx <- mut_to_wt_matching[j]
      optimal_matching[wt_idx] <- j
    }
    
    matching_type <- "transposed"
    unmatched_wt <- which(is.na(optimal_matching))
    
    result <- list(
      matching = optimal_matching,
      matching_type = matching_type,
      n_matched = n_mut,
      unmatched_wt = unmatched_wt,
      unmatched_mut = integer(0),
      base_condition = "MUT"
    )
  }
  
  return(result)
}

################################################################################
# 形态学评估和分类函数
################################################################################

classify_regulation_status_v10 <- function(shape_dist, amp_dist, phase_dist,
                                           range_wt, range_mut) {
  
  shape_changed <- shape_dist > CLASSIFY_THRESHOLDS$shape
  amp_changed <- amp_dist > CLASSIFY_THRESHOLDS$amp
  phase_changed <- phase_dist > CLASSIFY_THRESHOLDS$phase
  
  if (!shape_changed && !amp_changed && !phase_changed) {
    return("Conserved")
  }
  
  if (shape_changed) {
    return("Shape_Altered")
  }
  
  if (amp_changed) {
    if (range_mut < range_wt) {
      return("Amplitude_Repressed")
    } else {
      return("Amplitude_Enhanced")
    }
  }
  
  if (phase_changed) {
    return("Phase_Shifted")
  }
  
  return("Conserved")
}

assess_cluster_morphology <- function(center_wt, center_mut, 
                                       t) {
  
  metrics <- calc_ddtw_metrics(center_wt, center_mut, t)
  
  regulation_status <- classify_regulation_status_v10(
    metrics$shape_dist, metrics$amp_dist, metrics$phase_dist,
    metrics$range_A, metrics$range_B
  )
  
  return(list(
    cost = metrics$cost,
    shape_dist = metrics$shape_dist,
    amp_dist = metrics$amp_dist,
    phase_dist = metrics$phase_dist,
    ddtw_distance = metrics$ddtw_distance,
    pearson_r = metrics$pearson_r,
    deriv_r = metrics$deriv_r,
    range_wt = metrics$range_A,
    range_mut = metrics$range_B,
    regulation_status = regulation_status,
    dtw_path = metrics$dtw_path
  ))
}

classify_fate_hierarchical <- function(shape_dist, amp_dist, phase_dist,
                                        shape_threshold = CLASSIFY_THRESHOLDS$shape,
                                        amp_threshold = CLASSIFY_THRESHOLDS$amp,
                                        phase_threshold = CLASSIFY_THRESHOLDS$phase) {
  
  if (shape_dist > shape_threshold) {
    return("Shape_Diverged")
  }
  
  if (amp_dist > amp_threshold) {
    return("Amp_Altered")
  }
  
  if (phase_dist > phase_threshold) {
    return("Phase_Shifted")
  }
  
  return("Conserved")
}

calc_raw_fate_metrics <- function(wt_center, sub_center, t) {
  
  metrics <- calc_ddtw_metrics(wt_center, sub_center, t)
  
  return(list(
    shape_dist = metrics$shape_dist,
    amp_dist = metrics$amp_dist,
    phase_dist = metrics$phase_dist,
    ddtw_distance = metrics$ddtw_distance,
    pearson_r = metrics$pearson_r,
    deriv_r = metrics$deriv_r,
    range_wt = metrics$range_A,
    range_sub = metrics$range_B
  ))
}

################################################################################
# 可视化辅助函数
################################################################################

get_priority_category <- function(shape_dist, amp_dist, phase_dist) {
  if (shape_dist > CLASSIFY_THRESHOLDS$shape) return("Shape")
  if (amp_dist > CLASSIFY_THRESHOLDS$amp) return("Amp")
  if (phase_dist > CLASSIFY_THRESHOLDS$phase) return("Phase")
  return("Conserved")
}

get_priority_order <- function(category) {
  priority_map <- c(Conserved = 1, Shape = 2, Amp = 3, Phase = 4)
  order_val <- priority_map[category]
  if (is.na(order_val)) order_val <- 5
  return(order_val)
}

get_priority_color <- function(category) {
  color <- PRIORITY_COLORS[[category]]
  if (is.null(color)) color <- PRIORITY_COLORS$Conserved
  return(color)
}

build_dimension_labels <- function(shape_dist, amp_dist, phase_dist,
                                   range_wt, range_mut) {
  
  shape_changed <- shape_dist > CLASSIFY_THRESHOLDS$shape
  shape_label <- ifelse(shape_changed, "Shape:altered", "Shape:conserved")
  
  amp_changed <- amp_dist > CLASSIFY_THRESHOLDS$amp
  if (amp_changed) {
    amp_direction <- ifelse(range_mut > range_wt, "up", "down")
    amp_label <- paste0("Amp:", amp_direction)
  } else {
    amp_label <- "Amp:conserved"
  }
  
  phase_changed <- phase_dist > CLASSIFY_THRESHOLDS$phase
  if (phase_changed) {
    phase_label <- "Phase:shifted"
  } else {
    phase_label <- "Phase:conserved"
  }
  
  return(list(
    shape    = shape_label,
    amp      = amp_label,
    phase    = phase_label,
    combined = paste(shape_label, "|", amp_label, "|", phase_label)
  ))
}

################################################################################
# 可视化函数 (v11.1 更新 - 完全参数化显示名称和斜体控制)
################################################################################

#' 辅助函数: 生成带格式的图例标签
#' 
#' @param display_name 显示名称
#' @param use_italic 是否使用斜体
#' @return ggplot2 可用的标签表达式
make_formatted_label <- function(display_name, use_italic) {
  if (use_italic) {
    return(bquote(italic(.(display_name))))
  } else {
    return(display_name)
  }
}

#' 辅助函数: 生成带格式的标题
#' 
#' @param wt_cluster_id WT 聚类 ID
#' @param mut_cluster_id MUT 聚类 ID
#' @param wt_display WT 显示名称
#' @param mut_display MUT 显示名称
#' @param use_italic_wt WT 是否斜体
#' @param use_italic_mut MUT 是否斜体
#' @param title_type 标题类型: "profile", "venn", "fate"
make_plot_title <- function(wt_cluster_id, mut_cluster_id = NULL, 
                            wt_display, mut_display,
                            use_italic_wt, use_italic_mut,
                            title_type = "profile") {
  
  if (title_type == "profile") {
    if (use_italic_wt && use_italic_mut) {
      return(bquote(italic(.(wt_display)) ~ "C" * .(wt_cluster_id) ~ "-" ~ italic(.(mut_display)) ~ "C" * .(mut_cluster_id)))
    } else if (!use_italic_wt && use_italic_mut) {
      return(bquote(.(wt_display) ~ "C" * .(wt_cluster_id) ~ "-" ~ italic(.(mut_display)) ~ "C" * .(mut_cluster_id)))
    } else if (use_italic_wt && !use_italic_mut) {
      return(bquote(italic(.(wt_display)) ~ "C" * .(wt_cluster_id) ~ "-" ~ .(mut_display) ~ "C" * .(mut_cluster_id)))
    } else {
      return(bquote(.(wt_display) ~ "C" * .(wt_cluster_id) ~ "-" ~ .(mut_display) ~ "C" * .(mut_cluster_id)))
    }
  } else if (title_type == "venn") {
    if (use_italic_wt && use_italic_mut) {
      return(bquote(italic(.(wt_display)) ~ "C" * .(wt_cluster_id) ~ "\u2229" ~ italic(.(mut_display)) ~ "C" * .(mut_cluster_id)))
    } else if (!use_italic_wt && use_italic_mut) {
      return(bquote(.(wt_display) ~ "C" * .(wt_cluster_id) ~ "\u2229" ~ italic(.(mut_display)) ~ "C" * .(mut_cluster_id)))
    } else if (use_italic_wt && !use_italic_mut) {
      return(bquote(italic(.(wt_display)) ~ "C" * .(wt_cluster_id) ~ "\u2229" ~ .(mut_display) ~ "C" * .(mut_cluster_id)))
    } else {
      return(bquote(.(wt_display) ~ "C" * .(wt_cluster_id) ~ "\u2229" ~ .(mut_display) ~ "C" * .(mut_cluster_id)))
    }
  } else if (title_type == "fate") {
    if (use_italic_wt && use_italic_mut) {
      return(bquote(italic(.(wt_display)) ~ "C" * .(wt_cluster_id) ~ "Fate Divergence in" ~ italic(.(mut_display))))
    } else if (!use_italic_wt && use_italic_mut) {
      return(bquote(.(wt_display) ~ "C" * .(wt_cluster_id) ~ "Fate Divergence in" ~ italic(.(mut_display))))
    } else if (use_italic_wt && !use_italic_mut) {
      return(bquote(italic(.(wt_display)) ~ "C" * .(wt_cluster_id) ~ "Fate Divergence in" ~ .(mut_display)))
    } else {
      return(bquote(.(wt_display) ~ "C" * .(wt_cluster_id) ~ "Fate Divergence in" ~ .(mut_display)))
    }
  }
}

plot_profile_overlay <- function(wt_center, mut_center, 
                                 wt_sd = NULL, mut_sd = NULL,
                                 wt_cluster_id, mut_cluster_id,
                                 shape_dist, amp_dist, phase_dist,
                                 range_wt, range_mut,
                                 priority_category,
                                 time_labels,
                                 wt_display_name = CONDITION_WT_DISPLAY,
                                 mut_display_name = CONDITION_MUT_DISPLAY,
                                 use_italic_wt = USE_ITALIC_WT,
                                 use_italic_mut = USE_ITALIC_MUT) {
  
  n_points <- length(wt_center)
  time_numeric <- 1:n_points
  
  plot_data <- data.frame(
    time       = rep(time_numeric, 2),
    time_label = rep(time_labels, 2),
    expression = c(wt_center, mut_center),
    condition  = rep(c("WT", "MUT"), each = n_points)
  )
  
  has_ribbon <- FALSE
  if (!is.null(wt_sd) && !is.null(mut_sd)) {
    ribbon_wt <- data.frame(
      time = time_numeric, ymin = wt_center - wt_sd, ymax = wt_center + wt_sd, condition = "WT"
    )
    ribbon_mut <- data.frame(
      time = time_numeric, ymin = mut_center - mut_sd, ymax = mut_center + mut_sd, condition = "MUT"
    )
    ribbon_data <- rbind(ribbon_wt, ribbon_mut)
    has_ribbon <- TRUE
  }
  
  border_color <- get_priority_color(priority_category)
  dim_labels <- build_dimension_labels(shape_dist, amp_dist, phase_dist, range_wt, range_mut)
  
  # 生成格式化的图例标签
  wt_label <- make_formatted_label(wt_display_name, use_italic_wt)
  mut_label <- make_formatted_label(mut_display_name, use_italic_mut)
  
  # 生成格式化的标题
  plot_title <- make_plot_title(wt_cluster_id, mut_cluster_id,
                                 wt_display_name, mut_display_name,
                                 use_italic_wt, use_italic_mut,
                                 title_type = "profile")
  
  p <- ggplot(plot_data, aes(x = time, y = expression, color = condition)) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "gray90", linewidth = 0.3),
      panel.border      = element_rect(color = border_color, fill = NA, linewidth = 2.5),
      panel.background  = element_rect(fill = "transparent", color = NA),
      plot.background   = element_rect(fill = "transparent", color = NA),
      plot.title        = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle     = element_text(size = 8, hjust = 0.5, color = "#2C3E50", face = "bold"),
      legend.position      = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.background    = element_rect(fill = alpha("white", 0.85), color = NA),
      legend.key.size      = unit(0.4, "cm"),
      legend.text          = element_text(size = 8),
      legend.title         = element_blank(),
      axis.title           = element_text(size = 9),
      axis.text            = element_text(size = 8),
      plot.margin          = margin(10, 8, 8, 8)
    )
  
  if (has_ribbon) {
    p <- p + 
      geom_ribbon(data = ribbon_data %>% filter(condition == "WT"),
                  aes(x = time, ymin = ymin, ymax = ymax),
                  fill = COLOR_WT_RIBBON, color = NA, inherit.aes = FALSE) +
      geom_ribbon(data = ribbon_data %>% filter(condition == "MUT"),
                  aes(x = time, ymin = ymin, ymax = ymax),
                  fill = COLOR_MUT_RIBBON, color = NA, inherit.aes = FALSE)
  }
  
  p <- p +
    geom_line(linewidth = 1.1) +
    geom_point(size = 2.8, aes(shape = condition)) +
    scale_color_manual(
      values = c("WT" = COLOR_WT, "MUT" = COLOR_MUT),
      labels = c("WT" = wt_label, "MUT" = mut_label)
    ) +
    scale_shape_manual(
      values = c("WT" = 16, "MUT" = 17),
      labels = c("WT" = wt_label, "MUT" = mut_label)
    ) +
    scale_x_continuous(breaks = time_numeric, labels = time_labels) +
    labs(
      title = plot_title,
      subtitle = dim_labels$combined,
      x = "Time",
      y = "Standardized Expression (Z-score)"
    )
  
  return(p)
}

################################################################################
# Venn图可视化函数 (v10.1 更新 - 透明填充样式)
################################################################################

#' 为颜色添加透明度
#' @param col 颜色值
#' @param alpha 透明度 (0-1)
#' @return 带透明度的颜色
add_alpha <- function(col, alpha = 0.5) {
  rgb_vals <- col2rgb(col)
  rgb(rgb_vals[1], rgb_vals[2], rgb_vals[3], max = 255, alpha = alpha * 255)
}

plot_composition_venn <- function(wt_genes, mut_genes, wt_cluster_id, mut_cluster_id,
                                  wt_display_name = CONDITION_WT_DISPLAY,
                                  mut_display_name = CONDITION_MUT_DISPLAY,
                                  use_italic_wt = USE_ITALIC_WT,
                                  use_italic_mut = USE_ITALIC_MUT) {
  
  wt_set <- unique(wt_genes)
  mut_set <- unique(mut_genes)
  
  n_wt <- length(wt_set)
  n_mut <- length(mut_set)
  n_shared <- length(intersect(wt_set, mut_set))
  n_wt_only <- n_wt - n_shared
  n_mut_only <- n_mut - n_shared
  
  # 生成格式化的标题
  plot_title <- make_plot_title(wt_cluster_id, mut_cluster_id,
                                 wt_display_name, mut_display_name,
                                 use_italic_wt, use_italic_mut,
                                 title_type = "venn")
  
  # 确定标签的fontface
  wt_fontface <- if(use_italic_wt) "bold.italic" else "bold"
  mut_fontface <- if(use_italic_mut) "bold.italic" else "bold"
  
  p <- tryCatch({
    # 生成圆形路径
    theta <- seq(0, 2*pi, length.out = 200)
    radius <- 1
    
    # 两个圆的中心位置
    wt_center_x <- -0.5
    mut_center_x <- 0.5
    
    wt_x <- radius * cos(theta) + wt_center_x
    wt_y <- radius * sin(theta)
    mut_x <- radius * cos(theta) + mut_center_x
    mut_y <- radius * sin(theta)
    
    wt_circle_df <- data.frame(x = wt_x, y = wt_y)
    mut_circle_df <- data.frame(x = mut_x, y = mut_y)
    
    # 定义填充颜色（带透明度）
    fill_wt <- add_alpha(COLOR_WT, alpha = 0.5)
    fill_mut <- add_alpha(COLOR_MUT, alpha = 0.5)
    
    # 创建基础图形
    venn_plot <- ggplot() +
      # WT圆形填充（透明）
      geom_polygon(data = wt_circle_df, aes(x = x, y = y), 
                   fill = fill_wt, color = NA) +
      # MUT圆形填充（透明）
      geom_polygon(data = mut_circle_df, aes(x = x, y = y), 
                   fill = fill_mut, color = NA) +
      # WT圆形边框（实心）
      geom_path(data = wt_circle_df, aes(x = x, y = y), 
                color = COLOR_WT, linewidth = 2, lineend = "round") +
      # MUT圆形边框（实心）
      geom_path(data = mut_circle_df, aes(x = x, y = y), 
                color = COLOR_MUT, linewidth = 2, lineend = "round") +
      # 数字标注
      annotate("text", x = -0.85, y = 0, label = n_wt_only, 
               size = 4, fontface = "bold", color = "black") +
      annotate("text", x = 0, y = 0, label = n_shared, 
               size = 4, fontface = "bold", color = "black") +
      annotate("text", x = 0.85, y = 0, label = n_mut_only, 
               size = 4, fontface = "bold", color = "black") +
      # 条件标签 (使用参数化的fontface)
      annotate("text", x = wt_center_x, y = 1.35, label = wt_display_name, 
               size = 5, fontface = wt_fontface, color = COLOR_WT) +
      annotate("text", x = mut_center_x, y = 1.35, label = mut_display_name, 
               size = 5, fontface = mut_fontface, color = COLOR_MUT) +
      # 坐标和主题设置
      coord_fixed(ratio = 1) +
      theme_void() +
      theme(
        legend.position   = "none",
        panel.background  = element_rect(fill = "transparent", color = NA),
        plot.background   = element_rect(fill = "transparent", color = NA),
        plot.title        = element_text(size = 9, face = "bold", hjust = 0.5),
        plot.subtitle     = element_text(size = 7.5, hjust = 0.5, color = "gray40"),
        plot.margin       = margin(8, 8, 8, 8)
      ) +
      labs(
        title = plot_title,
        subtitle = sprintf("%s only: %d | Shared: %d | %s only: %d", 
                          wt_display_name, n_wt_only, n_shared, mut_display_name, n_mut_only)
      ) +
      xlim(-1.8, 1.8) + ylim(-1.3, 1.7)
    
    venn_plot
    
  }, error = function(e) {
    # 错误回退：简单文本显示
    ggplot() +
      annotate("text", x = 0.5, y = 0.5, 
               label = sprintf("WT C%d ∩ %s C%d\nWT only: %d | Shared: %d | %s only: %d",
                              wt_cluster_id, mut_display_name, mut_cluster_id, 
                              n_wt_only, n_shared, mut_display_name, n_mut_only),
               size = 3.5, hjust = 0.5, vjust = 0.5) +
      theme_void() +
      xlim(0, 1) + ylim(0, 1)
  })
  
  return(p)
}

plot_fate_divergence <- function(wt_center, subcluster_centers, 
                                  wt_cluster_id, n_subclusters,
                                  subcluster_sizes, fate_labels,
                                  time_labels,
                                  wt_display_name = CONDITION_WT_DISPLAY,
                                  mut_display_name = CONDITION_MUT_DISPLAY,
                                  use_italic_wt = USE_ITALIC_WT,
                                  use_italic_mut = USE_ITALIC_MUT) {
  
  n_points <- length(wt_center)
  time_numeric <- 1:n_points
  
  # 生成格式化的标题
  plot_title <- make_plot_title(wt_cluster_id, NULL,
                                 wt_display_name, mut_display_name,
                                 use_italic_wt, use_italic_mut,
                                 title_type = "fate")
  
  plot_data <- data.frame(
    time = time_numeric,
    expression = wt_center,
    cluster = sprintf("%s_original", wt_display_name),
    cluster_type = "WT",
    stringsAsFactors = FALSE
  )
  
  for (k in 1:n_subclusters) {
    fate_label <- fate_labels[k]
    subcluster_data <- data.frame(
      time = time_numeric,
      expression = subcluster_centers[k, ],
      cluster = sprintf("Sub_%d: %s (n=%d)", k, fate_label, subcluster_sizes[k]),
      cluster_type = "Fate",
      stringsAsFactors = FALSE
    )
    plot_data <- rbind(plot_data, subcluster_data)
  }
  
  cluster_levels <- c(sprintf("%s_original", wt_display_name), 
                      sprintf("Sub_%d: %s (n=%d)", 1:n_subclusters, fate_labels, subcluster_sizes))
  plot_data$cluster <- factor(plot_data$cluster, levels = cluster_levels)
  plot_data$cluster_type <- factor(plot_data$cluster_type, levels = c("WT", "Fate"))
  
  wt_cluster_name <- sprintf("%s_original", wt_display_name)
  color_vals <- setNames(COLOR_WT, wt_cluster_name)
  for (k in 1:n_subclusters) {
    cl_name <- sprintf("Sub_%d: %s (n=%d)", k, fate_labels[k], subcluster_sizes[k])
    fate_color <- FATE_COLORS[[fate_labels[k]]]
    if (is.null(fate_color)) fate_color <- SUBCLUSTER_COLORS[k]
    color_vals[cl_name] <- fate_color
  }
  
  linewidth_vals <- setNames(1.8, wt_cluster_name)
  size_vals <- setNames(3.5, wt_cluster_name)
  shape_vals <- c("WT" = 16, "Fate" = 17)
  for (k in 1:n_subclusters) {
    cl_name <- sprintf("Sub_%d: %s (n=%d)", k, fate_labels[k], subcluster_sizes[k])
    linewidth_vals[cl_name] <- 1.2
    size_vals[cl_name] <- 2.5
  }
  
  p <- ggplot(plot_data, aes(x = time, y = expression, color = cluster, group = cluster)) +
    geom_line(aes(linewidth = cluster)) +
    geom_point(aes(size = cluster, shape = cluster_type)) +
    scale_linewidth_manual(values = linewidth_vals, guide = "none") +
    scale_size_manual(values = size_vals, guide = "none") +
    scale_shape_manual(values = shape_vals, guide = "none") +
    scale_color_manual(values = color_vals) +
    scale_x_continuous(breaks = time_numeric, labels = time_labels) +
    theme_minimal(base_size = 10) +
    theme(
      panel.grid.minor  = element_blank(),
      panel.grid.major  = element_line(color = "gray90", linewidth = 0.3),
      panel.border      = element_rect(color = "#2C3E50", fill = NA, linewidth = 1.5),
      panel.background  = element_rect(fill = "transparent", color = NA),
      plot.background   = element_rect(fill = "transparent", color = NA),
      plot.title        = element_text(size = 11, face = "bold", hjust = 0.5),
      plot.subtitle     = element_text(size = 9, hjust = 0.5, color = "#7F8C8D"),
      legend.position   = "right",
      legend.background = element_rect(fill = alpha("white", 0.85), color = NA),
      legend.title      = element_blank(),
      axis.title        = element_text(size = 9),
      axis.text         = element_text(size = 8),
      plot.margin       = margin(10, 8, 8, 8)
    ) +
    labs(
      title = plot_title,
      subtitle = sprintf("%d subclusters (DDTW Classification)", n_subclusters),
      x = "Time",
      y = "Standardized Expression (Z-score)"
    )
  
  return(p)
}

################################################################################
# AgriGO输出函数
################################################################################

export_background_genes <- function(common_genes, output_dir) {
  bg_file <- file.path(output_dir, "AgriGO_Ready", "Background", "Universe_GeneID.txt")
  writeLines(common_genes, bg_file)
  cat(sprintf("✔ 背景基因集已保存: %s (%d genes)\n", bg_file, length(common_genes)))
}

export_pipeline5a_genelists <- function(wt_cluster_id, mut_cluster_id,
                                         wt_genes, mut_genes, output_dir) {
  
  pair_name <- sprintf("WT_C%02d_MUT_C%02d", wt_cluster_id, mut_cluster_id)
  pair_dir <- file.path(output_dir, "AgriGO_Ready", "Pipeline5a", pair_name)
  dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
  
  wt_set <- unique(wt_genes)
  mut_set <- unique(mut_genes)
  intersection <- intersect(wt_set, mut_set)
  wt_only <- setdiff(wt_set, mut_set)
  mut_only <- setdiff(mut_set, wt_set)
  
  writeLines(wt_set, file.path(pair_dir, "WT_ClusterGenes.txt"))
  writeLines(mut_set, file.path(pair_dir, "MUT_ClusterGenes.txt"))
  writeLines(intersection, file.path(pair_dir, "Intersection_Genes.txt"))
  writeLines(wt_only, file.path(pair_dir, "WT_Only_Genes.txt"))
  writeLines(mut_only, file.path(pair_dir, "MUT_Only_Genes.txt"))
  
  cat(sprintf("    ✔ %s: WT=%d, MUT=%d, ∩=%d, WT_only=%d, MUT_only=%d\n",
              pair_name, length(wt_set), length(mut_set), length(intersection),
              length(wt_only), length(mut_only)))
  
  return(list(
    pair_name = pair_name,
    n_wt = length(wt_set),
    n_mut = length(mut_set),
    n_intersection = length(intersection),
    n_wt_only = length(wt_only),
    n_mut_only = length(mut_only)
  ))
}

export_pipeline5b_genelists <- function(wt_cluster_id, subcluster_membership,
                                         output_dir, fate_labels = NULL,
                                         membership_cutoff = MFUZZ_PARAMS$membership_cutoff) {
  
  wt_dir <- file.path(output_dir, "AgriGO_Ready", "Pipeline5b", sprintf("WT_C%02d", wt_cluster_id))
  dir.create(wt_dir, recursive = TRUE, showWarnings = FALSE)
  
  n_subclusters <- ncol(subcluster_membership)
  export_summary <- list()
  
  for (k in 1:n_subclusters) {
    core_genes <- names(which(subcluster_membership[, k] >= membership_cutoff))
    
    if (length(core_genes) > 0) {
      file_name <- if (!is.null(fate_labels)) {
        sprintf("SubCluster_%d_%s_Genes.txt", k, fate_labels[k])
      } else {
        sprintf("SubCluster_%d_Genes.txt", k)
      }
      writeLines(core_genes, file.path(wt_dir, file_name))
    }
    
    export_summary[[k]] <- list(
      subcluster = k,
      n_genes = length(core_genes),
      fate_label = if (!is.null(fate_labels)) fate_labels[k] else NA
    )
    
    cat(sprintf("      SubCluster_%d (%s): %d genes\n", 
                k, if (!is.null(fate_labels)) fate_labels[k] else "NA", length(core_genes)))
  }
  
  return(export_summary)
}

################################################################################
# 主程序开始
################################################################################

cat("================================================================================\n")
cat("  Pipeline 5: Mfuzz时间序列聚类分析 (v10.1 工程化版本)\n")
cat("================================================================================\n\n")

cat(sprintf("开始时间: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("项目名称: %s\n", PROJECT_NAME))
cat(sprintf("输入 RData: %s\n", PARAMS$rdata))
cat(sprintf("输出目录: %s\n\n", OUTPUT_DIR))

cat("运行模式:\n")
cat(sprintf("  Pipeline5a (WT-MUT匹配): %s\n", ifelse(PARAMS$run_pipeline5a, "启用", "禁用")))
cat(sprintf("  Pipeline5b (命运分化): %s\n\n", ifelse(PARAMS$run_pipeline5b, "启用", "禁用")))

cat("条件标签配置:\n")
cat(sprintf("  WT条件: %s\n", CONDITION_WT))
cat(sprintf("  突变体条件: %s (显示名: %s)\n\n", CONDITION_MUT, CONDITION_MUT_DISPLAY))

################################################################################
# 步骤1: 加载数据 (v10.1 核心更新)
################################################################################

cat("--- 步骤1: 加载数据 ---\n")

load(PARAMS$rdata)
cat(sprintf("  ✔ 已加载 RData: %s\n", basename(PARAMS$rdata)))

# 检查并获取必要对象
if (exists("mfuzz_input")) {
  # 从 Pipeline4_Logic1_TimeCourse_for_Mfuzz.RData 加载
  vst_matrix <- mfuzz_input$vst_matrix
  metadata <- mfuzz_input$metadata
  all_time_genes <- mfuzz_input$responsive_genes
  cat(sprintf("  数据来源: Pipeline4 Logic1 RData\n"))
} else if (exists("vsd")) {
  # 从 DESeq2_objects.RData 加载
  vst_matrix <- assay(vsd)
  if (!exists("metadata")) metadata <- meta
  
  # 需要从外部文件获取时间响应基因
  pattern_dir <- file.path(DESEQ2_DIR, "Pattern_Feature")
  pattern_file <- file.path(pattern_dir, "Union_Time_Responsive_Genes_with_Cor.txt")
  if (file.exists(pattern_file)) {
    pattern_data <- read.table(pattern_file, header = TRUE, stringsAsFactors = FALSE)
    all_time_genes <- pattern_data$GeneID
  } else {
    # 使用所有基因
    all_time_genes <- rownames(vst_matrix)
    warning("未找到时间响应基因文件，使用全部基因")
  }
  cat(sprintf("  数据来源: DESeq2 Master RData\n"))
} else {
  stop("RData 文件格式不正确，未找到 mfuzz_input 或 vsd 对象")
}

cat(sprintf("  VST数据: %d 基因 x %d 样本\n", nrow(vst_matrix), ncol(vst_matrix)))
cat(sprintf("  时间响应基因: %d 个\n", length(all_time_genes)))

# 如果指定了 samples_file，覆盖 metadata
if (!is.null(PARAMS$samples_file) && file.exists(PARAMS$samples_file)) {
  metadata <- read.table(PARAMS$samples_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  if (!"SampleID" %in% colnames(metadata)) {
    colnames(metadata)[1] <- "SampleID"
  }
  rownames(metadata) <- metadata$SampleID
  cat(sprintf("  ✔ 已从 samples.txt 加载 metadata: %d 样本\n", nrow(metadata)))
}

################################################################################
# 步骤2: 自动检测时间点 (v10.1 核心新增)
################################################################################

cat("\n--- 步骤2: 自动检测时间点 ---\n")

TIME_INFO <- auto_detect_timepoints(metadata)

# 更新全局 MFUZZ_PARAMS
MFUZZ_PARAMS <- list(
  c_wt = PARAMS$c_wt,
  c_mut = PARAMS$c_mut,
  c_range = PARAMS$c_range_min:PARAMS$c_range_max,
  m_wt = PARAMS$m_wt,
  m_mut = PARAMS$m_mut,
  m_bounds = c(lower = PARAMS$m_lower, upper = PARAMS$m_upper),
  membership_cutoff = PARAMS$membership_cutoff,
  min_std = PARAMS$min_std,
  time_points = TIME_INFO$time_points,
  time_labels = TIME_INFO$time_labels,
  time_labels_raw = TIME_INFO$time_labels_raw,
  D = TIME_INFO$D
)

# 更新成本参数
COST_PARAMS <- list(
  w_shape = PARAMS$w_shape,
  w_amp = PARAMS$w_amp,
  w_phase = PARAMS$w_phase,
  interp_points = PARAMS$interp_points,
  dtw_window = PARAMS$dtw_window,
  derivative_method = PARAMS$derivative_method
)

# Pipeline5b参数
P5B_PARAMS <- list(
  c_range = PARAMS$p5b_c_range_min:PARAMS$p5b_c_range_max,
  min_genes = PARAMS$p5b_min_genes
)

# XB-Based参数
XB_PARAMS <- list(
  global_prominence_threshold = PARAMS$global_prominence_threshold,
  subcluster_prominence_threshold = PARAMS$subcluster_prominence_threshold,
  kneedle_S = PARAMS$kneedle_S
)

# 统一分类阈值
CLASSIFY_THRESHOLDS <- list(
  shape = PARAMS$classify_shape_threshold,
  amp = PARAMS$classify_amp_threshold,
  phase = PARAMS$classify_phase_threshold
)

cat(sprintf("\n已配置 %d 个时间点: %s\n", MFUZZ_PARAMS$D, 
            paste(MFUZZ_PARAMS$time_labels, collapse = ", ")))

################################################################################
# 打印参数摘要
################################################################################

cat("\nv10.1 DDTW正交三维指标体系:\n")
cat("  ============================================================\n")
cat("  | 指标  | 定义                          | 正交性保证        |\n")
cat("  |-------|-------------------------------|-------------------|\n")
cat("  | Shape | DDTW(导数序列)归一化距离       | 比较变化率，非原值 |\n")
cat("  | Amp   | |Range_A-Range_B|/(R_A+R_B)   | 时间无关标量      |\n")
cat("  | Phase | DTW路径偏离度                  | 弹性对齐产物      |\n")
cat("  ============================================================\n\n")

cat("DDTW参数:\n")
cat(sprintf("  插值点数: %d\n", COST_PARAMS$interp_points))
cat(sprintf("  导数方法: %s\n", COST_PARAMS$derivative_method))
cat(sprintf("  DTW窗口: %s\n\n", ifelse(is.null(COST_PARAMS$dtw_window), "无约束", COST_PARAMS$dtw_window)))

cat("成本函数权重:\n")
cat(sprintf("  w_shape = %.2f, w_amp = %.2f, w_phase = %.2f\n\n", 
            COST_PARAMS$w_shape, COST_PARAMS$w_amp, COST_PARAMS$w_phase))

cat("统一分类阈值:\n")
cat(sprintf("  Shape阈值: %.2f\n", CLASSIFY_THRESHOLDS$shape))
cat(sprintf("  Amp阈值: %.2f\n", CLASSIFY_THRESHOLDS$amp))
cat(sprintf("  Phase阈值: %.2f\n\n", CLASSIFY_THRESHOLDS$phase))

################################################################################
# 步骤3: 数据准备 (v10.1 更新)
################################################################################

cat("--- 步骤3: 数据准备 ---\n")

wt_ts <- prepare_single_timeseries(vst_matrix, metadata, CONDITION_WT, all_time_genes,
                                   TIME_INFO, MFUZZ_PARAMS$min_std)
mut_ts <- prepare_single_timeseries(vst_matrix, metadata, CONDITION_MUT, all_time_genes,
                                    TIME_INFO, MFUZZ_PARAMS$min_std)

common_genes <- intersect(wt_ts$filtered_genes, mut_ts$filtered_genes)
cat(sprintf("\n共同基因集: %d 基因\n", length(common_genes)))

if (length(common_genes) == 0) {
  stop("共同基因集为0！请检查:\n",
       "  1. metadata 的 Condition/Time 列是否正确\n",
       "  2. 条件标签 (--condition_wt/--condition_mut) 是否与 metadata 匹配")
}

wt_expr <- wt_ts$matrix[common_genes, ]
mut_expr <- mut_ts$matrix[common_genes, ]
N_genes <- length(common_genes)

wt_eset <- ExpressionSet(assayData = wt_expr)
wt_eset_std <- standardise(wt_eset)
mut_eset <- ExpressionSet(assayData = mut_expr)
mut_eset_std <- standardise(mut_eset)

export_background_genes(common_genes, OUTPUT_DIR)

################################################################################
# 步骤4: 估算模糊因子m
################################################################################

cat("\n--- 步骤4: 估算模糊因子m ---\n")

cat("\n【WT模糊因子估算】\n")
if (is.null(MFUZZ_PARAMS$m_wt)) {
  m_wt_result <- estimate_m_sj(N = N_genes, D = MFUZZ_PARAMS$D,
                                m_lower = MFUZZ_PARAMS$m_bounds["lower"],
                                m_upper = MFUZZ_PARAMS$m_bounds["upper"],
                                verbose = TRUE)
  m_wt <- m_wt_result$m
} else {
  m_wt <- MFUZZ_PARAMS$m_wt
  cat(sprintf("  使用用户指定值: m=%.4f\n", m_wt))
}

cat("\n【MUT模糊因子估算】\n")
if (is.null(MFUZZ_PARAMS$m_mut)) {
  m_mut_result <- estimate_m_sj(N = N_genes, D = MFUZZ_PARAMS$D,
                                 m_lower = MFUZZ_PARAMS$m_bounds["lower"],
                                 m_upper = MFUZZ_PARAMS$m_bounds["upper"],
                                 verbose = TRUE)
  m_mut <- m_mut_result$m
} else {
  m_mut <- MFUZZ_PARAMS$m_mut
  cat(sprintf("  使用用户指定值: m=%.4f\n", m_mut))
}

################################################################################
# Pipeline5a: WT-MUT独立聚类与匹配
################################################################################

if (PARAMS$run_pipeline5a) {
  
  cat("\n")
  cat("================================================================================\n")
  cat("  Pipeline5a: WT-MUT独立聚类与匹配分析 (v10.1 DDTW)\n")
  cat("================================================================================\n")
  
  # 步骤5a: 选择最佳聚类数c (PP-XB算法 v11.0)
  cat("\n--- 步骤5a: 选择最佳聚类数c (PP-XB算法 v11.0) ---\n")
  
  if (is.null(MFUZZ_PARAMS$c_wt)) {
    wt_xb_result <- scan_optimal_c(eset = wt_eset_std, c_range = MFUZZ_PARAMS$c_range,
                                    m = m_wt, condition_name = "WT", output_dir = OUTPUT_DIR,
                                    pipeline_type = "5a",
                                    prominence_threshold = XB_PARAMS$global_prominence_threshold)
    c_wt <- wt_xb_result$optimal_c
  } else {
    c_wt <- MFUZZ_PARAMS$c_wt
    cat(sprintf("\n【WT聚类数】使用用户指定值: c=%d\n", c_wt))
  }
  
  if (is.null(MFUZZ_PARAMS$c_mut)) {
    mut_xb_result <- scan_optimal_c(eset = mut_eset_std, c_range = MFUZZ_PARAMS$c_range,
                                     m = m_mut, condition_name = "MUT", output_dir = OUTPUT_DIR,
                                     pipeline_type = "5a",
                                     prominence_threshold = XB_PARAMS$global_prominence_threshold)
    c_mut <- mut_xb_result$optimal_c
  } else {
    c_mut <- MFUZZ_PARAMS$c_mut
    cat(sprintf("\n【MUT聚类数】使用用户指定值: c=%d\n", c_mut))
  }
  
  cat(sprintf("\n【最终聚类数】c_wt=%d, c_mut=%d\n\n", c_wt, c_mut))
  
  # 步骤6a: 执行Mfuzz聚类
  cat("--- 步骤6a: 执行Mfuzz聚类 ---\n")
  
  cat(sprintf("\nWT聚类: c=%d, m=%.4f\n", c_wt, m_wt))
  set.seed(123)
  wt_clusters <- mfuzz(wt_eset_std, c = c_wt, m = m_wt)
  cat(sprintf("  ✔ WT聚类完成\n"))
  
  cat(sprintf("\nMUT聚类: c=%d, m=%.4f\n", c_mut, m_mut))
  set.seed(123)
  mut_clusters <- mfuzz(mut_eset_std, c = c_mut, m = m_mut)
  cat(sprintf("  ✔ MUT聚类完成\n"))
  
  # 生成Mfuzz原生聚类图
  n_cols_plot <- min(4, max(c_wt, c_mut))
  n_rows_wt <- ceiling(c_wt / n_cols_plot)
  n_rows_mut <- ceiling(c_mut / n_cols_plot)
  
  pdf(file.path(OUTPUT_DIR, "plots", "WT_clusters_mfuzz.pdf"), 
      width = 3*n_cols_plot, height = 3*n_rows_wt)
  mfuzz.plot2(wt_eset_std, cl=wt_clusters, mfrow=c(n_rows_wt, n_cols_plot), 
              time.labels=MFUZZ_PARAMS$time_labels, centre=TRUE, x11=FALSE)
  dev.off()
  
  pdf(file.path(OUTPUT_DIR, "plots", "MUT_clusters_mfuzz.pdf"), 
      width = 3*n_cols_plot, height = 3*n_rows_mut)
  mfuzz.plot2(mut_eset_std, cl=mut_clusters, mfrow=c(n_rows_mut, n_cols_plot), 
              time.labels=MFUZZ_PARAMS$time_labels, centre=TRUE, x11=FALSE)
  dev.off()
  
  cat("  ✔ Mfuzz原生聚类图已保存\n")
  
  # 步骤7a: 计算DDTW成本矩阵
  cat("\n--- 步骤7a: 计算DDTW成本矩阵 ---\n")
  
  cost_results <- calc_composite_cost_matrix(
    wt_centers = wt_clusters$centers,
    mut_centers = mut_clusters$centers,
    time_points = MFUZZ_PARAMS$time_points
  )
  
  # 保存特征矩阵
  write.table(cost_results$cost_matrix,
              file.path(OUTPUT_DIR, "morphology_analysis", "Cost_Matrix_v10.txt"),
              sep = "\t", quote = FALSE)
  write.table(cost_results$shape_dist_matrix,
              file.path(OUTPUT_DIR, "morphology_analysis", "Shape_Dist_Matrix_v10.txt"),
              sep = "\t", quote = FALSE)
  write.table(cost_results$amp_dist_matrix,
              file.path(OUTPUT_DIR, "morphology_analysis", "Amp_Dist_Matrix_v10.txt"),
              sep = "\t", quote = FALSE)
  write.table(cost_results$phase_dist_matrix,
              file.path(OUTPUT_DIR, "morphology_analysis", "Phase_Dist_Matrix_v10.txt"),
              sep = "\t", quote = FALSE)
  write.table(cost_results$deriv_r_matrix,
              file.path(OUTPUT_DIR, "morphology_analysis", "Deriv_r_Matrix_v10.txt"),
              sep = "\t", quote = FALSE)
  
  cat("  ✔ DDTW特征矩阵已保存\n")
  
  # 步骤8a: 匈牙利算法匹配
  cat("\n--- 步骤8a: 匈牙利算法最优匹配 ---\n")
  
  matching_result <- solve_asymmetric_LSAP(cost_results$cost_matrix)
  optimal_matching <- matching_result$matching
  
  cat(sprintf("  匹配类型: %s\n", matching_result$matching_type))
  cat(sprintf("  匹配数量: %d\n", matching_result$n_matched))
  
  # 步骤9a: 形态学评估
  cat("\n--- 步骤9a: 形态学评估与分类 ---\n")
  
  morphology_results <- list()
  cluster_features <- data.frame()
  
  for (wt_idx in 1:c_wt) {
    mut_idx <- optimal_matching[wt_idx]
    
    if (is.na(mut_idx)) {
      morphology_results[[wt_idx]] <- list(
        wt_cluster = wt_idx,
        mut_cluster = NA,
        status = "unmatched"
      )
      next
    }
    
    morph <- assess_cluster_morphology(
      center_wt = wt_clusters$centers[wt_idx, ],
      center_mut = mut_clusters$centers[mut_idx, ],
      t = MFUZZ_PARAMS$time_points
    )
    
    morphology_results[[wt_idx]] <- c(
      list(wt_cluster = wt_idx, mut_cluster = mut_idx, status = "matched"),
      morph
    )
    
    priority_cat <- get_priority_category(morph$shape_dist, morph$amp_dist, morph$phase_dist)
    
    cluster_features <- rbind(cluster_features, data.frame(
      WT_Cluster = wt_idx,
      MUT_Cluster = mut_idx,
      Cost = morph$cost,
      Shape_Dist = morph$shape_dist,
      Amp_Dist = morph$amp_dist,
      Phase_Dist = morph$phase_dist,
      DDTW_Distance = morph$ddtw_distance,
      Pearson_r = morph$pearson_r,
      Deriv_r = morph$deriv_r,
      Range_WT = morph$range_wt,
      Range_MUT = morph$range_mut,
      Regulation_Status = morph$regulation_status,
      Priority_Category = priority_cat
    ))
    
    cat(sprintf("  WT_C%02d - MUT_C%02d: Cost=%.4f, Shape=%.3f, Amp=%.3f, Phase=%.3f -> %s\n",
                wt_idx, mut_idx, morph$cost, morph$shape_dist, morph$amp_dist, 
                morph$phase_dist, priority_cat))
  }
  
  write.table(cluster_features,
              file.path(OUTPUT_DIR, "morphology_analysis", "Cluster_Features_v10.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  # 步骤10a: 可视化
  cat("\n--- 步骤10a: 生成Pipeline5a可视化 ---\n")
  
  plots_profile <- list()
  plots_venn <- list()
  agrigo_exports <- list()
  
  for (wt_idx in 1:c_wt) {
    mut_idx <- optimal_matching[wt_idx]
    if (is.na(mut_idx)) next
    
    morph <- morphology_results[[wt_idx]]
    priority_cat <- get_priority_category(morph$shape_dist, morph$amp_dist, morph$phase_dist)
    
    # 获取聚类基因
    wt_genes <- names(which(wt_clusters$membership[, wt_idx] >= MFUZZ_PARAMS$membership_cutoff))
    mut_genes <- names(which(mut_clusters$membership[, mut_idx] >= MFUZZ_PARAMS$membership_cutoff))
    
    # Profile图
    p_profile <- plot_profile_overlay(
      wt_center = wt_clusters$centers[wt_idx, ],
      mut_center = mut_clusters$centers[mut_idx, ],
      wt_cluster_id = wt_idx,
      mut_cluster_id = mut_idx,
      shape_dist = morph$shape_dist,
      amp_dist = morph$amp_dist,
      phase_dist = morph$phase_dist,
      range_wt = morph$range_wt,
      range_mut = morph$range_mut,
      priority_category = priority_cat,
      time_labels = MFUZZ_PARAMS$time_labels,
      mut_display_name = CONDITION_MUT_DISPLAY
    )
    
    plots_profile[[wt_idx]] <- p_profile
    
    # Venn图
    p_venn <- plot_composition_venn(
      wt_genes = wt_genes,
      mut_genes = mut_genes,
      wt_cluster_id = wt_idx,
      mut_cluster_id = mut_idx,
      mut_display_name = CONDITION_MUT_DISPLAY
    )
    
    plots_venn[[wt_idx]] <- p_venn
    
    # AgriGO导出
    agrigo_exports[[wt_idx]] <- export_pipeline5a_genelists(
      wt_cluster_id = wt_idx,
      mut_cluster_id = mut_idx,
      wt_genes = wt_genes,
      mut_genes = mut_genes,
      output_dir = OUTPUT_DIR
    )
    
    # 保存单独图形
    ggsave(file.path(OUTPUT_DIR, "Pipeline5a_Figures", 
                     sprintf("WT_C%02d_MUT_C%02d_Profile.pdf", wt_idx, mut_idx)),
           p_profile, width = 6, height = 4.5, device = cairo_pdf, bg = "transparent")
    
    ggsave(file.path(OUTPUT_DIR, "Pipeline5a_Figures", 
                     sprintf("WT_C%02d_MUT_C%02d_Venn.pdf", wt_idx, mut_idx)),
           p_venn, width = 4, height = 4, device = cairo_pdf, bg = "transparent")
  }
  
  # 组合图
  valid_profile_plots <- plots_profile[!sapply(plots_profile, is.null)]
  valid_venn_plots <- plots_venn[!sapply(plots_venn, is.null)]
  
  if (length(valid_profile_plots) > 0) {
    n_pairs <- length(valid_profile_plots)
    n_cols <- min(3, n_pairs)
    n_rows <- ceiling(n_pairs / n_cols)
    
    composite_profile <- wrap_plots(valid_profile_plots, ncol = n_cols)
    ggsave(file.path(OUTPUT_DIR, "Pipeline5a_Figures", "Pipeline5a_All_Profiles.pdf"),
           composite_profile, width = 6 * n_cols, height = 4.5 * n_rows, 
           device = cairo_pdf, bg = "transparent")
    
    composite_venn <- wrap_plots(valid_venn_plots, ncol = n_cols)
    ggsave(file.path(OUTPUT_DIR, "Pipeline5a_Figures", "Pipeline5a_All_Venns.pdf"),
           composite_venn, width = 4 * n_cols, height = 4 * n_rows, 
           device = cairo_pdf, bg = "transparent")
    
    cat(sprintf("\n✔ Pipeline5a组合图已保存 (%d pairs)\n", n_pairs))
  }
  
  cat("\n✔ Pipeline5a 完成\n")
}

################################################################################
# Pipeline5b: 命运分化分析
################################################################################

if (PARAMS$run_pipeline5b) {
  
  cat("\n")
  cat("================================================================================\n")
  cat("  Pipeline5b: 命运分化分析 (v10.1 DDTW)\n")
  cat("================================================================================\n")
  
  fate_results <- list()
  plots_fate <- list()
  sankey_data_list <- list()
  
  for (wt_idx in 1:c_wt) {
    cat(sprintf("\n--- WT Cluster %d ---\n", wt_idx))
    
    wt_genes <- names(which(wt_clusters$membership[, wt_idx] >= MFUZZ_PARAMS$membership_cutoff))
    n_genes <- length(wt_genes)
    
    if (n_genes < P5B_PARAMS$min_genes) {
      cat(sprintf("  跳过: 基因数 (%d) < 最小阈值 (%d)\n", n_genes, P5B_PARAMS$min_genes))
      fate_results[[wt_idx]] <- list(
        wt_cluster = wt_idx,
        n_genes = n_genes,
        n_subclusters = NA,
        status = "skipped_low_genes"
      )
      next
    }
    
    # 提取MUT中的表达数据
    mut_sub_matrix <- exprs(mut_eset_std)[wt_genes, , drop = FALSE]
    
    if (nrow(mut_sub_matrix) == 0) {
      fate_results[[wt_idx]] <- list(
        wt_cluster = wt_idx,
        n_genes = n_genes,
        n_subclusters = NA,
        status = "no_expression_data"
      )
      next
    }
    
    # 创建子ExpressionSet
    sub_eset <- ExpressionSet(assayData = mut_sub_matrix)
    
    # PP-XB扫描选择最佳c (v11.0 使用Pipeline5b特殊约束)
    sub_xb_result <- scan_optimal_c(
      eset = sub_eset,
      c_range = P5B_PARAMS$c_range,
      m = m_mut,
      condition_name = sprintf("WT_C%02d_Sub", wt_idx),
      output_dir = OUTPUT_DIR,
      pipeline_type = "5b",  # 使用Pipeline5b特殊约束
      prominence_threshold = XB_PARAMS$subcluster_prominence_threshold,
      save_plot = TRUE
    )
    
    c_sub <- sub_xb_result$optimal_c
    
    # Mfuzz聚类
    set.seed(123)
    sub_clusters <- mfuzz(sub_eset, c = c_sub, m = m_mut)
    
    # 计算Fate分类
    fate_labels <- character(c_sub)
    fate_metrics_list <- list()
    
    for (k in 1:c_sub) {
      metrics <- calc_raw_fate_metrics(
        wt_center = wt_clusters$centers[wt_idx, ],
        sub_center = sub_clusters$centers[k, ],
        t = MFUZZ_PARAMS$time_points
      )
      
      fate_labels[k] <- classify_fate_hierarchical(
        metrics$shape_dist, metrics$amp_dist, metrics$phase_dist
      )
      
      fate_metrics_list[[k]] <- metrics
    }
    
    subcluster_sizes <- sapply(1:c_sub, function(k) {
      sum(sub_clusters$membership[, k] >= MFUZZ_PARAMS$membership_cutoff)
    })
    
    # 存储结果
    fate_results[[wt_idx]] <- list(
      wt_cluster = wt_idx,
      n_genes = n_genes,
      n_subclusters = c_sub,
      status = "analyzed",
      xb_result = sub_xb_result$xb_result,
      subcluster_centers = sub_clusters$centers,
      subcluster_membership = sub_clusters$membership,
      fate_labels = fate_labels,
      fate_metrics = fate_metrics_list,
      subcluster_sizes = subcluster_sizes
    )
    
    # Sankey数据
    for (k in 1:c_sub) {
      sankey_data_list[[length(sankey_data_list) + 1]] <- data.frame(
        Source = sprintf("WT_C%02d", wt_idx),
        Target = fate_labels[k],
        Value = subcluster_sizes[k],
        SubCluster = k,
        stringsAsFactors = FALSE
      )
    }
    
    # AgriGO导出
    cat(sprintf("    导出基因列表:\n"))
    export_pipeline5b_genelists(
      wt_cluster_id = wt_idx,
      subcluster_membership = sub_clusters$membership,
      output_dir = OUTPUT_DIR,
      fate_labels = fate_labels,
      membership_cutoff = MFUZZ_PARAMS$membership_cutoff
    )
    
    # Fate Divergence图
    plots_fate[[wt_idx]] <- plot_fate_divergence(
      wt_center = wt_clusters$centers[wt_idx, ],
      subcluster_centers = sub_clusters$centers,
      wt_cluster_id = wt_idx,
      n_subclusters = c_sub,
      subcluster_sizes = subcluster_sizes,
      fate_labels = fate_labels,
      time_labels = MFUZZ_PARAMS$time_labels,
      mut_display_name = CONDITION_MUT_DISPLAY
    )
    
    fate_file <- sprintf("WT_C%02d_Fate_Divergence.pdf", wt_idx)
    ggsave(file.path(OUTPUT_DIR, "Pipeline5b_Figures", fate_file),
           plots_fate[[wt_idx]], width = 8, height = 4, device = cairo_pdf, bg = "transparent")
    cat(sprintf("  ✔ %s\n", fate_file))
  }
  
  # 组合Fate图
  valid_fate_plots <- plots_fate[!sapply(plots_fate, is.null)]
  if (length(valid_fate_plots) > 0) {
    n_rows_fate <- length(valid_fate_plots)
    fate_composite <- wrap_plots(valid_fate_plots, ncol = 1)
    composite_fate_height <- max(6, n_rows_fate * 4)
    
    ggsave(file.path(OUTPUT_DIR, "Pipeline5b_Figures", "Pipeline5b_Composite.pdf"),
           fate_composite, width = 9, height = composite_fate_height, 
           device = cairo_pdf, bg = "transparent")
    cat(sprintf("\n✔ Pipeline5b组合图已保存 (%d clusters)\n", n_rows_fate))
  }
  
  # Sankey源数据
  cat("\n--- 输出Sankey源数据 ---\n")
  
  if (length(sankey_data_list) > 0) {
    sankey_df <- do.call(rbind, sankey_data_list)
    
    sankey_file <- file.path(OUTPUT_DIR, "morphology_analysis", "Pipeline5b_Sankey_Source_Data_v10.csv")
    write.csv(sankey_df, sankey_file, row.names = FALSE)
    cat(sprintf("✔ Sankey源数据已保存: %s (%d rows)\n", sankey_file, nrow(sankey_df)))
    
    fate_summary_table <- sankey_df %>%
      group_by(Target) %>%
      summarise(
        Total_Genes = sum(Value),
        N_SubClusters = n(),
        .groups = "drop"
      ) %>%
      arrange(desc(Total_Genes))
    
    cat("\n  Fate分类汇总:\n")
    for (i in 1:nrow(fate_summary_table)) {
      cat(sprintf("    %s: %d genes (%d subclusters)\n",
                  fate_summary_table$Target[i],
                  fate_summary_table$Total_Genes[i],
                  fate_summary_table$N_SubClusters[i]))
    }
  }
  
  # Fate Summary
  fate_summary <- do.call(rbind, lapply(fate_results, function(x) {
    curve_type <- if (!is.null(x$xb_result)) x$xb_result$curve_class$curve_type else NA
    selection_method <- if (!is.null(x$xb_result)) x$xb_result$selection_method else NA
    fate_dist <- if (!is.null(x$fate_labels)) paste(table(x$fate_labels), collapse="|") else NA
    
    data.frame(
      WT_Cluster = x$wt_cluster,
      N_Genes = x$n_genes,
      N_SubClusters = x$n_subclusters,
      Status = x$status,
      Curve_Type = curve_type,
      Selection_Method = selection_method,
      Fate_Distribution = fate_dist,
      stringsAsFactors = FALSE
    )
  }))
  
  write.table(fate_summary, 
              file.path(OUTPUT_DIR, "morphology_analysis", "Pipeline5b_Fate_Summary_v10.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  
  cat("\n✔ Pipeline5b 完成\n")
}

################################################################################
# 保存所有结果
################################################################################

cat("\n--- 保存所有结果 ---\n")

# 保存Membership矩阵
if (exists("wt_clusters")) {
  wt_membership_df <- as.data.frame(wt_clusters$membership)
  wt_membership_df$gene <- rownames(wt_membership_df)
  colnames(wt_membership_df) <- c(paste0("WT_C", 1:ncol(wt_clusters$membership)), "gene")
  wt_membership_df <- wt_membership_df[, c("gene", paste0("WT_C", 1:ncol(wt_clusters$membership)))]
  write.table(wt_membership_df,
              file.path(OUTPUT_DIR, "morphology_analysis", "WT_Membership_Matrix.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("✔ WT_Membership_Matrix.txt\n")
}

if (exists("mut_clusters")) {
  mut_membership_df <- as.data.frame(mut_clusters$membership)
  mut_membership_df$gene <- rownames(mut_membership_df)
  colnames(mut_membership_df) <- c(paste0("MUT_C", 1:ncol(mut_clusters$membership)), "gene")
  mut_membership_df <- mut_membership_df[, c("gene", paste0("MUT_C", 1:ncol(mut_clusters$membership)))]
  write.table(mut_membership_df,
              file.path(OUTPUT_DIR, "morphology_analysis", "MUT_Membership_Matrix.txt"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("✔ MUT_Membership_Matrix.txt\n")
}

# 保存RData
save_objects <- c("common_genes", "m_wt", "MFUZZ_PARAMS", "COST_PARAMS", 
                  "XB_PARAMS", "CLASSIFY_THRESHOLDS", "TIME_INFO")
if (exists("m_mut")) save_objects <- c(save_objects, "m_mut")
if (exists("wt_clusters")) save_objects <- c(save_objects, "wt_clusters", "wt_eset_std")
if (exists("mut_clusters")) save_objects <- c(save_objects, "mut_clusters", "mut_eset_std")
if (exists("cost_results")) save_objects <- c(save_objects, "cost_results", "optimal_matching", 
                                               "matching_result", "morphology_results", "cluster_features")
if (exists("fate_results")) save_objects <- c(save_objects, "fate_results")

save(list = save_objects,
     file = file.path(OUTPUT_DIR, "rdata", "Mfuzz_v11_Results.RData"))
cat("✔ Mfuzz_v11_Results.RData\n")

# 保存参数JSON
params_json <- list(
  version = "11.1",
  version_update = "标准化重构 - 变量名标准化 + 可视化参数化",
  algorithm = list(
    pipeline5a = "WT-MUT独立聚类 + 匈牙利算法匹配",
    pipeline5b = "WT聚类后分析MUT中的命运分化",
    metrics = "DDTW正交三维指标体系 (Derivative DTW)",
    cluster_selection = "PP-XB统一框架 (Kwon指数 + 曲线形态分类)"
  ),
  ppxb_algorithm = list(
    description = "基于两份独立研究报告交叉验证的统一最优聚类数选择算法",
    kwon_index = "添加惩罚项修正XB单调递减：(compactness + penalty) / min_sep",
    curve_classification = c("U_SHAPED_OR_VALLEY", "MONOTONIC_DECREASING", "FLUCTUATING"),
    sample_constraint = "c_max ≤ min(√N, N/10), Pipeline5b更保守"
  ),
  time_info = list(
    n_timepoints = TIME_INFO$D,
    time_points = TIME_INFO$time_points,
    time_labels = TIME_INFO$time_labels,
    time_labels_raw = TIME_INFO$time_labels_raw
  ),
  conditions = list(
    wt = CONDITION_WT,
    mut = CONDITION_MUT,
    wt_display = CONDITION_WT_DISPLAY,
    mut_display = CONDITION_MUT_DISPLAY,
    use_italic_wt = USE_ITALIC_WT,
    use_italic_mut = USE_ITALIC_MUT
  ),
  ddtw_metrics = list(
    shape = "DDTW(导数序列)归一化距离 - 比较变化率实现与振幅解耦",
    amp = "SDR = |Range_A - Range_B| / (Range_A + Range_B) - 时间无关",
    phase = "DTW路径偏离度 - 从弹性对齐中自然得出"
  ),
  cost_params = list(
    w_shape = COST_PARAMS$w_shape,
    w_amp = COST_PARAMS$w_amp,
    w_phase = COST_PARAMS$w_phase,
    interp_points = COST_PARAMS$interp_points,
    derivative_method = COST_PARAMS$derivative_method
  ),
  classify_thresholds = list(
    shape = CLASSIFY_THRESHOLDS$shape,
    amp = CLASSIFY_THRESHOLDS$amp,
    phase = CLASSIFY_THRESHOLDS$phase
  ),
  unified_colors = UNIFIED_COLORS,
  mfuzz = list(
    c_wt = if(exists("c_wt")) c_wt else NA,
    c_mut = if(exists("c_mut")) c_mut else NA,
    m_wt = m_wt,
    m_mut = if(exists("m_mut")) m_mut else NA
  ),
  xb_params = XB_PARAMS,
  p5b_params = P5B_PARAMS,
  paths = list(
    rdata_input = PARAMS$rdata,
    output_dir = OUTPUT_DIR,
    project_name = PROJECT_NAME
  ),
  analysis_time = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
)

writeLines(
  jsonlite::toJSON(params_json, pretty = TRUE, auto_unbox = TRUE),
  file.path(OUTPUT_DIR, "morphology_analysis", "Analysis_Parameters_v11.json")
)
cat("✔ Analysis_Parameters_v11.json\n")

################################################################################
# 完成
################################################################################

cat("\n================================================================================\n")
cat("  Pipeline 5 (v11.1 标准化版本) 完成\n")
cat("================================================================================\n\n")

cat(sprintf("结束时间: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

cat("结果摘要:\n")
cat(sprintf("  共同基因数: %d\n", length(common_genes)))
cat(sprintf("  时间点数: %d (%s)\n", MFUZZ_PARAMS$D, paste(MFUZZ_PARAMS$time_labels, collapse=", ")))
if (PARAMS$run_pipeline5a && exists("c_wt")) {
  cat(sprintf("  Pipeline5a - WT聚类数: %d, MUT聚类数: %d\n", c_wt, c_mut))
}
if (PARAMS$run_pipeline5b && exists("fate_results")) {
  analyzed_count <- sum(sapply(fate_results, function(x) x$status == "analyzed"))
  cat(sprintf("  Pipeline5b - 分析的WT clusters: %d\n", analyzed_count))
}

cat("\nv11.1 标准化重构更新:\n")
cat("  • 变量名标准化: dof -> mut\n")
cat("  • 可视化完全参数化 (DISPLAY参数控制)\n")
cat("  • 斜体格式由参数控制 (use_italic_wt/mut)\n")
cat("  • 输出文件名保持标准化 (WT/MUT)\n")

cat("\n输出目录结构:\n")
cat(sprintf("  %s/\n", OUTPUT_DIR))
cat("    ├── plots/               - Mfuzz原生聚类图\n")
cat("    ├── cluster_selection/   - PP-XB双面板诊断图\n")
cat("    ├── Pipeline5a_Figures/  - WT-MUT匹配可视化\n")
cat("    ├── Pipeline5b_Figures/  - 命运分化可视化\n")
cat("    ├── morphology_analysis/ - 特征矩阵和参数\n")
cat("    ├── AgriGO_Ready/        - GO分析用基因列表\n")
cat("    └── rdata/               - R数据对象\n")

sink(type = "message")
sink(type = "output")
close(log_con)

cat("\n日志已保存至: ", log_file, "\n")