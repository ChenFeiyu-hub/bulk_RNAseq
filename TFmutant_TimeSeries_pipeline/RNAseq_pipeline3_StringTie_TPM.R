#!/usr/bin/env Rscript

# ==============================================================================
# 脚本名称: RNAseq_pipeline3_StringTie_TPM.R
# 功能: 替代原 Python 脚本，使用 StringTie 计算 TPM 并合并矩阵
# 依赖: optparse, tidyverse, parallel (R base)
# ==============================================================================

# 1. 库加载
suppressPackageStartupMessages({
  library(optparse)
  library(parallel)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
  library(tools)
})

# 2. 参数解析
option_list <- list(
  make_option(c("-b", "--bam_dir"), type="character", help="BAM 文件所在目录"),
  make_option(c("-o", "--out_dir"), type="character", help="输出目录"),
  make_option(c("-g", "--gtf"), type="character", help="GTF 参考文件路径"),
  make_option(c("-s", "--stringtie_bin"), type="character", help="StringTie 可执行文件路径"),
  make_option(c("-n", "--project_name"), type="character", help="项目名称前缀"),
  make_option(c("-t", "--threads"), type="integer", default=4, help="每个 StringTie 任务的线程数"),
  make_option(c("-j", "--jobs"), type="integer", default=8, help="并行任务数 (核心数)")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# 检查必要参数
if (is.null(opt$bam_dir) || is.null(opt$out_dir) || is.null(opt$gtf) || is.null(opt$stringtie_bin)) {
  print_help(opt_parser)
  stop("缺少必要参数！", call.=FALSE)
}

# 创建输出目录
if (!dir.exists(opt$out_dir)) {
  dir.create(opt$out_dir, recursive = TRUE)
}

# 设置日志 (简单的控制台输出)
log_info <- function(msg) {
  cat(sprintf("[%s] [INFO] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

log_error <- function(msg) {
  cat(sprintf("[%s] [ERROR] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg), file = stderr())
}

# ==============================================================================
# 核心函数: 单个样本处理
# ==============================================================================
process_sample <- function(bam_file) {
  tryCatch({
    # 提取样本名 (去除 .sort.bam)
    sample_name <- basename(bam_file) %>% str_remove("\\.sort\\.bam$")
    
    # 定义样本特定目录
    sample_dir <- file.path(opt$out_dir, sample_name)
    if (!dir.exists(sample_dir)) dir.create(sample_dir)
    
    # 定义输出文件路径
    gtf_out <- file.path(sample_dir, paste0(sample_name, ".gtf"))
    gene_abund_tab <- file.path(sample_dir, paste0(sample_name, ".gene_abund.tab"))
    
    # 构建 StringTie 命令
    # 保持与原 Python 脚本一致的参数: -e, -B, -A
    cmd <- paste(
      opt$stringtie_bin,
      "-e -B",
      "-p", opt$threads,
      "-G", opt$gtf,
      "-A", gene_abund_tab,
      "-o", gtf_out,
      bam_file
    )
    
    # 执行命令 (system2 相比 system 更安全)
    # 忽略 stdout, stderr 除非出错，防止刷屏
    exit_code <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    
    if (exit_code != 0) {
      return(list(status = "failed", sample = sample_name, error = "StringTie execution returned non-zero exit code"))
    }
    
    # 检查输出文件是否存在
    if (!file.exists(gene_abund_tab)) {
      return(list(status = "failed", sample = sample_name, error = "Output abundance file not found"))
    }
    
    # 读取结果提取 TPM
    # StringTie -A 输出格式: Gene ID, Gene Name, Reference, Strand, Start, End, Coverage, FPKM, TPM
    df <- read_tsv(gene_abund_tab, show_col_types = FALSE) %>%
      select(`Gene ID`, TPM) %>%
      rename(!!sample_name := TPM) # 动态重命名列为样本名
    
    return(list(status = "success", sample = sample_name, data = df))
    
  }, error = function(e) {
    return(list(status = "failed", sample = basename(bam_file), error = e$message))
  })
}

# ==============================================================================
# 主流程
# ==============================================================================

log_info(paste("Start StringTie Pipeline for project:", opt$project_name))

# 1. 查找 BAM 文件
bam_files <- list.files(opt$bam_dir, pattern = "\\.sort\\.bam$", full.names = TRUE)

if (length(bam_files) == 0) {
  log_error(paste("No BAM files found in", opt$bam_dir))
  quit(status = 1)
}

log_info(paste("Found", length(bam_files), "BAM files. Processing with", opt$jobs, "parallel jobs..."))

# 2. 并行处理 (使用 mclapply，Linux 下最稳定)
results <- mclapply(bam_files, process_sample, mc.cores = opt$jobs)

# 3. 结果检查与合并
success_list <- list()
failed_list <- list()

for (res in results) {
  if (res$status == "success") {
    success_list[[res$sample]] <- res$data
    log_info(paste("Finished:", res$sample))
  } else {
    failed_list[[res$sample]] <- res$error
    log_error(paste("Failed:", res$sample, "-", res$error))
  }
}

# 如果全部失败，退出
if (length(success_list) == 0) {
  log_error("All samples failed processing!")
  quit(status = 1)
}

# 如果部分失败，警告但不退出 (或者你可以选择严格模式退出)
if (length(failed_list) > 0) {
  log_error(paste(length(failed_list), "samples failed. Check logs."))
  # quit(status = 1) # 如果你希望任何一个样本失败就终止流程，取消注释这行
}

# 4. 合并矩阵
log_info("Merging TPM matrices...")

# 使用 purrr::reduce 进行全连接 (Full Join)
merged_df <- success_list %>%
  purrr::reduce(full_join, by = "Gene ID") %>%
  mutate(across(where(is.numeric), ~replace_na(., 0))) # 将 NA 填充为 0

# 5. 输出结果
out_file <- file.path(opt$out_dir, paste0(opt$project_name, ".tpm_matrix.txt"))
write_tsv(merged_df, out_file)

log_info(paste("TPM Matrix saved to:", out_file))
log_info("Pipeline 3 Completed.")