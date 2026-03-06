#!/usr/bin/env Rscript
# script_name: RNAseq_pipeline4_DESeq2.R
# 版本: v4.0 (Engineering Update)
# 功能: DESeq2 标准化流程，输出包含 vsd/dds 的全功能 RData

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(argparser)
})

# ================= 参数解析 =================
p <- arg_parser("DESeq2 Analysis for Time-Course & Multi-Mutants")
p <- add_argument(p, "--count_file", help="FeatureCounts matrix")
p <- add_argument(p, "--metadata_file", help="Metadata file")
p <- add_argument(p, "--out_dir", help="Output directory")
p <- add_argument(p, "--project_name", help="Project Name (for filename prefix)", default="Project")
p <- add_argument(p, "--control_cond", help="Control Condition (e.g., WT)")
p <- add_argument(p, "--control_time", help="Control Time point (e.g., 00h)")
p <- add_argument(p, "--fdr", help="FDR cutoff", default=0.05)
p <- add_argument(p, "--lfc", help="Log2FC threshold", default=1.0)
argv <- parse_args(p)

# ================= 1. 数据加载与预处理 =================
cat(">>> Loading Data...\n")
counts <- read.table(argv$count_file, header=TRUE, row.names=1, sep="\t", check.names=FALSE)
counts <- round(as.matrix(counts))

meta <- read.table(argv$metadata_file, header=TRUE, sep="\t", stringsAsFactors=FALSE)
# 确保样本ID匹配
if(! "SampleID" %in% colnames(meta)) {
    # 尝试第一列作为SampleID
    colnames(meta)[1] <- "SampleID"
}
rownames(meta) <- meta$SampleID

common_samples <- intersect(colnames(counts), rownames(meta))
counts <- counts[, common_samples]
meta <- meta[common_samples, ]

# 因子化处理
# 自动检测列名 (兼容大小写)
col_time <- grep("time", colnames(meta), ignore.case=TRUE, value=TRUE)[1]
col_cond <- grep("condition|genotype", colnames(meta), ignore.case=TRUE, value=TRUE)[1]

meta$Time <- factor(meta[[col_time]], levels = sort(unique(meta[[col_time]]))) 
meta$Condition <- factor(meta[[col_cond]])
meta$Group <- factor(paste(meta$Condition, meta$Time, sep="_"))

cat(sprintf(">>> Samples: %d | Conditions: %d | TimePoints: %d\n", 
            ncol(counts), length(unique(meta$Condition)), length(unique(meta$Time))))

# 构建 DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~ Group)

# 过滤低表达基因 (至少3个样本count >= 10)
keep <- rowSums(counts >= 10) >= 3
dds <- dds[keep, ]
cat(sprintf(">>> Genes retained after filtering: %d\n", nrow(dds)))

# 运行 DESeq
cat(">>> Running DESeq (Wald Test)...\n")
dds <- DESeq(dds)

# VST 转换 (用于PCA和Mfuzz)
cat(">>> Performing VST transformation...\n")
vsd <- vst(dds, blind=FALSE)
vst_mat <- assay(vsd)

# ========================================================
# [新增] 输出 VST 矩阵文本文件 (供 Pipeline 6 使用)
# ========================================================
vst_out_file <- file.path(argv$out_dir, "data", paste0(argv$project_name, "_VST_transformed.txt"))
write.table(vst_mat, file=vst_out_file, sep="\t", quote=FALSE, col.names=NA)
cat(sprintf("   [Saved] VST Matrix text file: %s\n", basename(vst_out_file)))

# ================= 逻辑 1: Time-Course Analysis =================
cat("\n>>> [Logic 1] Running Time-Course Comparisons (vs T0)...\n")

conditions <- unique(as.character(meta$Condition))
times <- unique(as.character(meta$Time))
ctrl_time <- argv$control_time

time_responsive_genes <- c()
res_list_time <- list()

for (cond in conditions) {
  base_group <- paste(cond, ctrl_time, sep="_")
  if (!(base_group %in% levels(meta$Group))) next
  
  for (t in times) {
    if (t == ctrl_time) next 
    target_group <- paste(cond, t, sep="_")
    if (!(target_group %in% levels(meta$Group))) next
    
    contrast_name <- paste0("Time_", cond, "_", t, "_vs_", ctrl_time) # 命名标准化
    
    # Target vs Base
    res <- results(dds, contrast=c("Group", target_group, base_group), alpha=argv$fdr)
    
    sig_genes <- rownames(res)[which(res$padj < argv$fdr & abs(res$log2FoldChange) > argv$lfc)]
    time_responsive_genes <- unique(c(time_responsive_genes, sig_genes))
    
    out_file <- file.path(argv$out_dir, "Logic1_TimeCourse", paste0(contrast_name, ".txt"))
    write.table(as.data.frame(res), out_file, sep="\t", quote=FALSE)
    res_list_time[[contrast_name]] <- res
  }
}

cat(sprintf(">>> Identified %d unique time-responsive genes.\n", length(time_responsive_genes)))

# [关键修复] 保存 Logic 1 RData
# 不仅包含 Mfuzz 数据，还包含 vsd 和 dds，供 PCA 使用
mfuzz_input <- list(
  vst_matrix = vst_mat[time_responsive_genes, ], 
  metadata = meta,
  responsive_genes = time_responsive_genes,
  parameters = list(control_time = ctrl_time)
)

save(mfuzz_input, res_list_time, vsd, dds, meta,
     file = file.path(argv$out_dir, "RData", "Pipeline4_Logic1_TimeCourse_for_Mfuzz.RData"))
cat("   [Saved] Logic 1 RData (Contains vsd/dds for PCA)\n")


# ================= 逻辑 2: Cross-Genotype Analysis =================
cat("\n>>> [Logic 2] Running Cross-Genotype Comparisons...\n")

res_list_genotype <- list()
cond_combinations <- combn(conditions, 2, simplify=FALSE)

for (t in times) {
  for (combo in cond_combinations) {
    c1 <- combo[1]; c2 <- combo[2]
    
    # 确定参考组 (WT优先)
    if (c1 == argv$control_cond) { ref <- c1; target <- c2 }
    else if (c2 == argv$control_cond) { ref <- c2; target <- c1 }
    else { sorted <- sort(c(c1, c2)); ref <- sorted[1]; target <- sorted[2] }
    
    group_ref <- paste(ref, t, sep="_")
    group_target <- paste(target, t, sep="_")
    
    if (!all(c(group_ref, group_target) %in% levels(meta$Group))) next
    
    contrast_name <- paste0("Genotype_", target, "_vs_", ref, "_at_", t)
    
    res <- results(dds, contrast=c("Group", group_target, group_ref), alpha=argv$fdr)
    
    out_file <- file.path(argv$out_dir, "Logic2_GenotypeComp", paste0(contrast_name, ".txt"))
    write.table(as.data.frame(res), out_file, sep="\t", quote=FALSE)
    res_list_genotype[[contrast_name]] <- res
  }
}

# ================= 最终汇总: Master RData =================
# 这是为了满足 Mfuzz 和 Venn 脚本对 "DESeq2_objects.RData" 的需求
# 包含所有比较结果 (all_comparisons)

cat("\n>>> Saving Consolidated Master RData...\n")

all_comparisons <- c(res_list_time, res_list_genotype)
project_rdata <- file.path(argv$out_dir, "RData", paste0(argv$project_name, "_DESeq2_objects.RData"))

save(dds, vsd, meta, all_comparisons, 
     file = project_rdata)

cat(sprintf("   [Saved] Master RData: %s\n", basename(project_rdata)))
cat(">>> Pipeline 4 All Done.\n")