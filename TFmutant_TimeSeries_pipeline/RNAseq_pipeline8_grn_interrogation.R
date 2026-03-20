#!/usr/bin/env Rscript
# =============================================================================
# RNAseq_pipeline8_grn_interrogation.R
# 功能: GRN 锚点网络的交互式深度查询与辅助分析
# 版本: v1.1
# 输入: Pipeline 7c 输出的 4 个 CSV 文件 + 查询基因列表
# 设计: 模块化分析，衔接 07_run_GRN.sh 流程
#
# 分析模块 (--mode):
#   gene_profile       - 基因中心拓扑与表达全景报告
#   topology_table     - 高 Degree 节点边类型分布表 (类似论文表3-7)
#   shared_targets     - 两个 TF 的共享靶基因分析 (类似论文表3-8)
#   overlap_matrix     - 多 TF 靶基因交叉重叠矩阵
#   expression_divergence - 按聚类×边类型的 WT-dof 时序分歧统计
#   subspace           - 调控子空间鉴定 (公共池 vs 私有集)
#   full_report        - 运行所有模块的综合报告
#
# Changelog:
#   v1.1 - 修复: WT_Cluster 误匹配为 LFC 列 (强制 \d+h$ 模式);
#          修复: Symbol == Gene_ID 时视为 NA;
#          移除: delta_outdegree (与 7c Fig6 重复, 7c 使用完整 TF 集);
#          新增: --batch 参数控制输出子目录;
#          优化: 输出日志格式
# =============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# =============================================================================
# 1. 命令行参数
# =============================================================================
option_list <- list(
  make_option("--edges_csv",  type = "character", help = "Pipeline7c_Fig4_Cytoscape_Edges.csv"),
  make_option("--nodes_csv",  type = "character", help = "Pipeline7c_Fig4_Cytoscape_Nodes.csv"),
  make_option("--target_lfc", type = "character", help = "Pipeline7c_Fig7a_Target_LFC_Matrix.csv"),
  make_option("--tf_lfc",     type = "character", help = "Pipeline7c_Fig7b_TF_LFC_Matrix.csv"),
  make_option("--input_dir",  type = "character", default = NULL,
              help = "Pipeline7c 输出目录 (自动定位4个CSV, 替代逐一指定)"),
  make_option("--query_genes", type = "character", default = NULL,
              help = "查询基因ID, 逗号分隔 (如 Cre12.g521150,Cre12.g495100)"),
  make_option("--query_file",  type = "character", default = NULL,
              help = "查询基因列表文件 (每行一个Gene_ID)"),
  make_option("--mode", type = "character", default = "gene_profile",
              help = paste0("分析模式: gene_profile, topology_table, shared_targets,\n",
                           "                           overlap_matrix,\n",
                           "                           expression_divergence, subspace, full_report")),
  make_option("--degree_cutoff", type = "integer", default = 15,
              help = "topology_table 模式的 Degree 阈值 [default: 15]"),
  make_option("--ref_tf", type = "character", default = NULL,
              help = "overlap_matrix/subspace 模式的参考 TF (用于计算重叠度)"),
  make_option("--out_dir",    type = "character", default = ".",
              help = "输出根目录 [default: .]"),
  make_option("--batch",      type = "character", default = NULL,
              help = "批次/项目标签, 作为输出子目录名 (如 CrDof_PSR1_query)"),
  make_option("--out_prefix", type = "character", default = "P8",
              help = "输出文件前缀 [default: P8]"),
  make_option("--condition_wt",  type = "character", default = "WT",
              help = "WT条件名 [default: WT]"),
  make_option("--condition_mut", type = "character", default = "dof",
              help = "突变体条件名 [default: dof]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# =============================================================================
# 2. 加载数据
# =============================================================================
cat("\n╔══════════════════════════════════════════════════════════╗\n")
cat(  "║  Pipeline 8: GRN Interrogation & Deep Query (v1.1)     ║\n")
cat(  "╚══════════════════════════════════════════════════════════╝\n\n")

# 自动定位CSV
if (!is.null(opt$input_dir)) {
  d <- opt$input_dir
  if (is.null(opt$edges_csv))  opt$edges_csv  <- file.path(d, "Pipeline7c_Fig4_Cytoscape_Edges.csv")
  if (is.null(opt$nodes_csv))  opt$nodes_csv  <- file.path(d, "Pipeline7c_Fig4_Cytoscape_Nodes.csv")
  if (is.null(opt$target_lfc)) opt$target_lfc <- file.path(d, "Pipeline7c_Fig7a_Target_LFC_Matrix.csv")
  if (is.null(opt$tf_lfc))     opt$tf_lfc     <- file.path(d, "Pipeline7c_Fig7b_TF_LFC_Matrix.csv")
}

for (f in c("edges_csv", "nodes_csv", "target_lfc", "tf_lfc")) {
  fp <- opt[[f]]
  if (is.null(fp) || !file.exists(fp)) {
    stop(sprintf("❌ 文件缺失: --%s = %s", f, ifelse(is.null(fp), "(未指定)", fp)))
  }
}

cat("  Loading data...\n")
edges  <- read.csv(opt$edges_csv,  stringsAsFactors = FALSE)
nodes  <- read.csv(opt$nodes_csv,  stringsAsFactors = FALSE)
tg_lfc <- read.csv(opt$target_lfc, stringsAsFactors = FALSE)
tf_lfc <- read.csv(opt$tf_lfc,     stringsAsFactors = FALSE)

cat(sprintf("    Edges: %d | Nodes: %d\n", nrow(edges), nrow(nodes)))
cat(sprintf("    Target LFC genes: %d | TF LFC genes: %d\n", nrow(tg_lfc), nrow(tf_lfc)))

# ── v1.1修复: Symbol == Gene_ID → NA ──
if ("Symbol" %in% names(tg_lfc)) {
  tg_lfc$Symbol[!is.na(tg_lfc$Symbol) & tg_lfc$Symbol == tg_lfc$Gene_ID] <- NA_character_
}
if ("Symbol" %in% names(tf_lfc)) {
  tf_lfc$Symbol[!is.na(tf_lfc$Symbol) & tf_lfc$Symbol == tf_lfc$Gene_ID] <- NA_character_
}
if ("Display_Label" %in% names(nodes)) {
  nodes$Display_Label_Clean <- ifelse(
    nodes$Display_Label == nodes$Node_ID, NA_character_, nodes$Display_Label
  )
} else {
  nodes$Display_Label_Clean <- NA_character_
}

# 解析查询基因
query_genes <- c()
if (!is.null(opt$query_genes)) {
  query_genes <- trimws(unlist(strsplit(opt$query_genes, ",")))
}
if (!is.null(opt$query_file) && file.exists(opt$query_file)) {
  file_genes <- trimws(readLines(opt$query_file))
  file_genes <- file_genes[file_genes != "" & !grepl("^#", file_genes)]
  query_genes <- unique(c(query_genes, file_genes))
}

WT  <- opt$condition_wt
MUT <- opt$condition_mut

# ── v1.1修复: 精确识别 LFC 列, 排除 WT_Cluster 等 ──
lfc_time_re <- "\\d+h$"
all_cols_tg <- names(tg_lfc)
all_cols_tf <- names(tf_lfc)

lfc_cols_tg <- all_cols_tg[
  (grepl(paste0("^", WT, "_"), all_cols_tg) | grepl(paste0("^", MUT, "_"), all_cols_tg)) &
  grepl(lfc_time_re, all_cols_tg)
]
lfc_cols_tf <- all_cols_tf[
  (grepl(paste0("^", WT, "_"), all_cols_tf) | grepl(paste0("^", MUT, "_"), all_cols_tf)) &
  grepl(lfc_time_re, all_cols_tf)
]

time_points <- unique(gsub(paste0("^", WT, "_|^", MUT, "_"), "", lfc_cols_tg))
time_num <- as.numeric(gsub("h$", "", time_points))
time_points <- time_points[order(time_num)]

cat(sprintf("    Time points: %s\n", paste(time_points, collapse = ", ")))
cat(sprintf("    LFC columns: %d (target) + %d (TF)\n",
            length(lfc_cols_tg), length(lfc_cols_tf)))
cat(sprintf("    Query genes: %s\n",
            ifelse(length(query_genes) > 0,
                   paste(query_genes, collapse = ", "), "(none)")))
cat(sprintf("    Mode: %s\n", opt$mode))

# ── 输出目录 ──
if (!is.null(opt$batch) && nchar(opt$batch) > 0) {
  OUT <- file.path(opt$out_dir, opt$batch)
  cat(sprintf("    Batch: %s\n", opt$batch))
} else {
  OUT <- file.path(opt$out_dir, "Pipeline8")
}
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
PREFIX <- opt$out_prefix
cat(sprintf("    Output dir: %s\n\n", OUT))


# =============================================================================
# 3. 辅助函数
# =============================================================================

get_gene_edges <- function(gene_id) {
  out_edges <- edges %>% filter(Source == gene_id) %>% mutate(Direction = "Out")
  in_edges  <- edges %>% filter(Target == gene_id) %>% mutate(Direction = "In")
  bind_rows(out_edges, in_edges)
}

get_edge_type_distribution <- function(gene_edges) {
  if (nrow(gene_edges) == 0) return(data.frame())
  total <- nrow(gene_edges)
  gene_edges %>%
    group_by(Edge_Class) %>%
    summarise(Count = n(), Pct = sprintf("%.0f%%", n() / total * 100), .groups = "drop")
}

get_gene_lfc <- function(gene_id) {
  row_tg <- tg_lfc %>% filter(Gene_ID == gene_id)
  if (nrow(row_tg) > 0) {
    lfc_data <- row_tg[1, lfc_cols_tg, drop = FALSE]
    meta <- row_tg[1, intersect(c("Gene_ID","Symbol","WT_Cluster","Edge_Dominant","Is_Highlight"), names(row_tg)), drop = FALSE]
    meta$Source <- "Target"
    return(list(meta = meta, lfc = lfc_data))
  }
  row_tf <- tf_lfc %>% filter(Gene_ID == gene_id)
  if (nrow(row_tf) > 0) {
    lfc_data <- row_tf[1, lfc_cols_tf, drop = FALSE]
    meta <- row_tf[1, intersect(c("Gene_ID","Symbol","Family","Is_Focus_TF","Out_Degree"), names(row_tf)), drop = FALSE]
    meta$Source <- "TF"
    return(list(meta = meta, lfc = lfc_data))
  }
  return(NULL)
}

compute_wt_dof_diff <- function(lfc_data) {
  diffs <- data.frame(Time = time_points, WT_LFC = NA_real_,
                      MUT_LFC = NA_real_, Diff = NA_real_,
                      stringsAsFactors = FALSE)
  for (tp in time_points) {
    wt_col  <- paste0(WT, "_", tp)
    mut_col <- paste0(MUT, "_", tp)
    idx <- which(diffs$Time == tp)
    if (wt_col %in% names(lfc_data) && mut_col %in% names(lfc_data)) {
      wt_val  <- as.numeric(lfc_data[[wt_col]])
      mut_val <- as.numeric(lfc_data[[mut_col]])
      diffs$WT_LFC[idx]  <- wt_val
      diffs$MUT_LFC[idx] <- mut_val
      diffs$Diff[idx]    <- wt_val - mut_val
    }
  }
  diffs
}

get_tf_targets <- function(tf_id) {
  edges %>% filter(Source == tf_id) %>% pull(Target) %>% unique()
}

fmt_gene <- function(gene_id) {
  node_row <- nodes %>% filter(Node_ID == gene_id)
  if (nrow(node_row) > 0) {
    ni <- node_row[1, ]
    sym <- ni$Display_Label_Clean
    if (is.na(sym) || sym == "") {
      tg_row <- tg_lfc %>% filter(Gene_ID == gene_id)
      if (nrow(tg_row) > 0 && !is.na(tg_row$Symbol[1])) sym <- tg_row$Symbol[1]
      if (is.na(sym) || sym == "") {
        tf_row <- tf_lfc %>% filter(Gene_ID == gene_id)
        if (nrow(tf_row) > 0 && !is.na(tf_row$Symbol[1])) sym <- tf_row$Symbol[1]
      }
    }
    fam <- ni$Family
    has_sym <- !is.na(sym) && sym != ""
    has_fam <- !is.na(fam) && fam != "" && fam != "NA"
    if (has_sym && has_fam) return(sprintf("%s (%s, %s)", gene_id, sym, fam))
    if (has_sym) return(sprintf("%s (%s)", gene_id, sym))
    if (has_fam) return(sprintf("%s (%s)", gene_id, fam))
  }
  gene_id
}

write_table_to_report <- function(df, file_conn, title = NULL, indent = "  ") {
  if (!is.null(title)) writeLines(sprintf("\n%s%s", indent, title), file_conn)
  col_widths <- sapply(names(df), function(cn) {
    vals <- nchar(as.character(df[[cn]]))
    max(nchar(cn), max(vals[!is.na(vals)], 2L))  # 2L = "NA" 的宽度兜底
  })
  header <- paste(sprintf(paste0("%-", col_widths, "s"), names(df)), collapse = "  ")
  sep_line <- paste(sapply(col_widths, function(w) paste(rep("-", w), collapse = "")), collapse = "  ")
  writeLines(paste0(indent, header), file_conn)
  writeLines(paste0(indent, sep_line), file_conn)
  for (i in seq_len(nrow(df))) {
    row_str <- paste(sprintf(paste0("%-", col_widths, "s"), as.character(df[i, ])), collapse = "  ")
    writeLines(paste0(indent, row_str), file_conn)
  }
}


# =============================================================================
# 模块 1: gene_profile
# =============================================================================
run_gene_profile <- function(genes) {
  cat("--- Module: Gene Profile ---\n")
  if (length(genes) == 0) { cat("  ⚠ 未指定查询基因\n"); return(invisible(NULL)) }

  report_file <- file.path(OUT, sprintf("%s_Gene_Profile_Report.txt", PREFIX))
  csv_file    <- file.path(OUT, sprintf("%s_Gene_Profile_Edges.csv", PREFIX))
  lfc_file    <- file.path(OUT, sprintf("%s_Gene_Profile_LFC_Diff.csv", PREFIX))

  fc <- file(report_file, "w")
  all_edges_list <- list()
  all_diffs_list <- list()

  writeLines("================================================================", fc)
  writeLines("  GRN Gene Profile Report", fc)
  writeLines(sprintf("  Generated: %s", Sys.time()), fc)
  writeLines(sprintf("  Query genes: %d", length(genes)), fc)
  writeLines("================================================================\n", fc)

  for (gid in genes) {
    writeLines("────────────────────────────────────────────────────────────", fc)
    writeLines(sprintf("  Gene: %s", fmt_gene(gid)), fc)
    writeLines("────────────────────────────────────────────────────────────", fc)

    node_info <- nodes %>% filter(Node_ID == gid)
    if (nrow(node_info) > 0) {
      ni <- node_info[1, ]
      writeLines(sprintf("  NodeType: %s | Family: %s | Degree: %d",
                         ni$NodeType, ni$Family, ni$Degree), fc)
      writeLines(sprintf("  Is_Focus_TF: %s | Is_TF: %s | Is_Highlight: %s | Cluster: %s",
                         ni$Is_Focus_TF, ni$Is_TF, ni$Is_Highlight, ni$Cluster_Label), fc)
    } else {
      writeLines("  ⚠ 不在锚点网络节点表中", fc)
    }

    gene_edges <- get_gene_edges(gid)
    out_n <- sum(gene_edges$Direction == "Out")
    in_n  <- sum(gene_edges$Direction == "In")
    writeLines(sprintf("\n  边统计: Total=%d (Out=%d, In=%d)", nrow(gene_edges), out_n, in_n), fc)

    if (nrow(gene_edges) > 0) {
      dist <- get_edge_type_distribution(gene_edges)
      write_table_to_report(dist, fc, title = "边类型分布:")
      dir_dist <- gene_edges %>%
        group_by(Direction, Edge_Class) %>%
        summarise(Count = n(), .groups = "drop") %>% arrange(Direction, desc(Count))
      write_table_to_report(dir_dist, fc, title = "方向×边类型:")

      if (out_n > 0) {
        writeLines("\n  出边靶基因:", fc)
        out_tgts <- gene_edges %>% filter(Direction == "Out") %>%
          select(Target, Edge_Class, Delta_r, Weight, MECS_Score, Tier)
        for (j in seq_len(nrow(out_tgts))) {
          r <- out_tgts[j, ]
          writeLines(sprintf("    → %s | %s | Δr=%.3f | W=%.3f | MECS=%.3f | T%s",
                             fmt_gene(r$Target), r$Edge_Class,
                             r$Delta_r, r$Weight, r$MECS_Score, r$Tier), fc)
        }
      }
      if (in_n > 0) {
        writeLines("\n  入边来源TF:", fc)
        in_srcs <- gene_edges %>% filter(Direction == "In") %>%
          select(Source, Edge_Class, Delta_r, Weight, MECS_Score, Tier)
        for (j in seq_len(nrow(in_srcs))) {
          r <- in_srcs[j, ]
          writeLines(sprintf("    ← %s | %s | Δr=%.3f | W=%.3f | MECS=%.3f | T%s",
                             fmt_gene(r$Source), r$Edge_Class,
                             r$Delta_r, r$Weight, r$MECS_Score, r$Tier), fc)
        }
      }
      gene_edges$Query_Gene <- gid
      all_edges_list[[gid]] <- gene_edges
    }

    lfc_info <- get_gene_lfc(gid)
    if (!is.null(lfc_info)) {
      writeLines("\n  差异表达时间序列 (log2FC vs 0h):", fc)
      diffs <- compute_wt_dof_diff(lfc_info$lfc)
      diffs_fmt <- diffs %>%
        mutate(across(c(WT_LFC, MUT_LFC, Diff), ~sprintf("%.3f", .)))
      names(diffs_fmt) <- c("Time", paste0(WT,"_LFC"), paste0(MUT,"_LFC"), paste0(WT,"-",MUT))
      write_table_to_report(diffs_fmt, fc)
      diffs$Gene_ID <- gid
      all_diffs_list[[gid]] <- diffs
    } else {
      writeLines("\n  ⚠ 不在LFC矩阵中", fc)
    }
    writeLines("", fc)
  }
  close(fc)

  if (length(all_edges_list) > 0) {
    write.csv(bind_rows(all_edges_list), csv_file, row.names = FALSE)
    cat(sprintf("  ✔ %s\n", basename(csv_file)))
  }
  if (length(all_diffs_list) > 0) {
    write.csv(bind_rows(all_diffs_list), lfc_file, row.names = FALSE)
    cat(sprintf("  ✔ %s\n", basename(lfc_file)))
  }
  cat(sprintf("  ✔ %s\n", basename(report_file)))
  return(report_file)
}


# =============================================================================
# 模块 2: topology_table
# =============================================================================
run_topology_table <- function(degree_cutoff) {
  cat(sprintf("--- Module: Topology Table (Degree >= %d) ---\n", degree_cutoff))
  high_deg <- nodes %>% filter(Degree >= degree_cutoff) %>% arrange(desc(Degree))
  if (nrow(high_deg) == 0) { cat("  ⚠ 无满足条件的节点\n"); return(invisible(NULL)) }

  get_n <- function(ec_df, cls) { v <- ec_df$n[ec_df$Edge_Class == cls]; if(length(v)==0) 0L else v }
  pct_f <- function(x, total) sprintf("%.0f%%", x / total * 100)

  results <- lapply(seq_len(nrow(high_deg)), function(i) {
    gid <- high_deg$Node_ID[i]
    ge <- get_gene_edges(gid)
    ec <- ge %>% group_by(Edge_Class) %>% summarise(n = n(), .groups = "drop")
    total <- nrow(ge)
    data.frame(
      Gene_ID = gid, Family = high_deg$Family[i], Degree = high_deg$Degree[i],
      Out = sum(ge$Direction == "Out"), In = sum(ge$Direction == "In"),
      Conserved = get_n(ec,"Conserved"), Conserved_Pct = pct_f(get_n(ec,"Conserved"), total),
      Lost = get_n(ec,"Lost"), Lost_Pct = pct_f(get_n(ec,"Lost"), total),
      Gained = get_n(ec,"Gained"), Gained_Pct = pct_f(get_n(ec,"Gained"), total),
      Rewired = get_n(ec,"Rewired"),
      stringsAsFactors = FALSE
    )
  })
  topo_df <- bind_rows(results)

  if (!is.null(opt$ref_tf)) {
    ref_targets <- get_tf_targets(opt$ref_tf)
    ref_n <- length(ref_targets)
    cat(sprintf("  Ref TF: %s (%d targets)\n", fmt_gene(opt$ref_tf), ref_n))
    topo_df$Overlap_with_Ref <- sapply(topo_df$Gene_ID, function(gid) {
      if (gid == opt$ref_tf) return("---")
      tf_tgts <- get_tf_targets(gid)
      if (length(tf_tgts) == 0) return("--- (pure target)")
      ov <- length(intersect(tf_tgts, ref_targets))
      sprintf("%d/%d(%.0f%%)", ov, ref_n, ov / ref_n * 100)
    })
  }

  out_csv <- file.path(OUT, sprintf("%s_Topology_Table.csv", PREFIX))
  write.csv(topo_df, out_csv, row.names = FALSE)
  report_file <- file.path(OUT, sprintf("%s_Topology_Table.txt", PREFIX))
  fc <- file(report_file, "w")
  writeLines(sprintf("高Degree节点拓扑特征 (Degree >= %d, N = %d)\n", degree_cutoff, nrow(topo_df)), fc)
  write_table_to_report(topo_df, fc, indent = "")
  close(fc)
  cat(sprintf("  ✔ %d nodes → %s + %s\n", nrow(topo_df), basename(out_csv), basename(report_file)))
  return(topo_df)
}


# =============================================================================
# 模块 3: shared_targets
# =============================================================================
run_shared_targets <- function(genes) {
  cat("--- Module: Shared Targets ---\n")
  if (length(genes) < 2) { cat("  ⚠ 需要至少2个TF\n"); return(invisible(NULL)) }

  tf1 <- genes[1]; tf2 <- genes[2]
  cat(sprintf("  TF1: %s\n  TF2: %s\n", fmt_gene(tf1), fmt_gene(tf2)))

  edges_tf1 <- edges %>% filter(Source == tf1)
  edges_tf2 <- edges %>% filter(Source == tf2)
  targets_tf1 <- edges_tf1$Target
  targets_tf2 <- edges_tf2$Target
  shared <- intersect(targets_tf1, targets_tf2)
  cat(sprintf("  TF1 targets: %d | TF2 targets: %d | Shared: %d\n",
              length(targets_tf1), length(targets_tf2), length(shared)))
  if (length(shared) == 0) { cat("  无共享靶基因\n"); return(invisible(NULL)) }

  shared_df <- data.frame(Shared_Target = shared, stringsAsFactors = FALSE)
  shared_df$TF1_Edge <- edges_tf1$Edge_Class[match(shared, edges_tf1$Target)]
  shared_df$TF2_Edge <- edges_tf2$Edge_Class[match(shared, edges_tf2$Target)]
  shared_df$Edge_Consistent <- shared_df$TF1_Edge == shared_df$TF2_Edge

  shared_df <- shared_df %>%
    left_join(nodes %>% select(Node_ID, Cluster_Label), by = c("Shared_Target" = "Node_ID")) %>%
    left_join(tg_lfc %>% select(Gene_ID, Symbol, WT_Cluster), by = c("Shared_Target" = "Gene_ID"))
  shared_df$Cluster <- ifelse(!is.na(shared_df$WT_Cluster), shared_df$WT_Cluster, shared_df$Cluster_Label)

  last_tp <- time_points[length(time_points)]
  shared_df$Diff_Last_TP <- sapply(shared_df$Shared_Target, function(g) {
    lfc_info <- get_gene_lfc(g)
    if (is.null(lfc_info)) return(NA_real_)
    diffs <- compute_wt_dof_diff(lfc_info$lfc)
    last_row <- diffs %>% filter(Time == last_tp)
    if (nrow(last_row) > 0) last_row$Diff[1] else NA_real_
  })

  tf1_short <- gsub("\\.", "", tf1); tf2_short <- gsub("\\.", "", tf2)
  out_csv <- file.path(OUT, sprintf("%s_Shared_Targets_%s_vs_%s.csv", PREFIX, tf1_short, tf2_short))
  write.csv(shared_df, out_csv, row.names = FALSE)

  report_file <- file.path(OUT, sprintf("%s_Shared_Targets_Report.txt", PREFIX))
  fc <- file(report_file, "w")
  writeLines(sprintf("共享靶基因: %s vs %s", fmt_gene(tf1), fmt_gene(tf2)), fc)
  writeLines(sprintf("共享: %d / TF1总: %d / TF2总: %d",
                     length(shared), length(targets_tf1), length(targets_tf2)), fc)
  writeLines(sprintf("TF2中与TF1重叠: %.0f%%", length(shared)/length(targets_tf2)*100), fc)
  writeLines(sprintf("边类型一致率: %.0f%% (%d/%d)",
                     sum(shared_df$Edge_Consistent)/nrow(shared_df)*100,
                     sum(shared_df$Edge_Consistent), nrow(shared_df)), fc)
  writeLines(sprintf("差值时间点: %s\n", last_tp), fc)

  export_df <- shared_df %>%
    select(Shared_Target, Symbol, Cluster, TF1_Edge, TF2_Edge, Edge_Consistent, Diff_Last_TP) %>%
    mutate(Diff_Last_TP = ifelse(is.na(Diff_Last_TP), "N/A", sprintf("%.3f", Diff_Last_TP)))
  names(export_df)[4:5] <- c(paste0(tf1,"_Edge"), paste0(tf2,"_Edge"))
  names(export_df)[7] <- paste0(WT,"-",MUT,"_",last_tp)
  write_table_to_report(export_df, fc, indent = "")
  close(fc)

  cat(sprintf("  ✔ %s\n  ✔ %s\n", basename(out_csv), basename(report_file)))
  return(shared_df)
}


# =============================================================================
# 模块 4: overlap_matrix
# =============================================================================
run_overlap_matrix <- function(genes) {
  cat("--- Module: Target Overlap Matrix ---\n")
  if (length(genes) == 0) {
    tf_list <- edges %>% group_by(Source) %>% summarise(n=n(),.groups="drop") %>%
      filter(n >= 5) %>% pull(Source)
    cat(sprintf("  Auto: 所有出边>=5的TF (%d)\n", length(tf_list)))
  } else { tf_list <- genes }
  if (length(tf_list) < 2) { cat("  ⚠ 需要>=2个TF\n"); return(invisible(NULL)) }

  target_sets <- lapply(tf_list, get_tf_targets); names(target_sets) <- tf_list
  n <- length(tf_list)
  overlap_mat <- matrix(0, n, n, dimnames = list(tf_list, tf_list))
  overlap_pct <- matrix("", n, n, dimnames = list(tf_list, tf_list))
  for (i in seq_len(n)) for (j in seq_len(n)) {
    ni <- length(target_sets[[i]]); nj <- length(target_sets[[j]])
    common <- length(intersect(target_sets[[i]], target_sets[[j]]))
    overlap_mat[i,j] <- common
    overlap_pct[i,j] <- if(i==j) sprintf("%d",ni) else
      sprintf("%d/%d(%.0f%%)", common, nj, ifelse(nj>0, common/nj*100, 0))
  }

  out1 <- file.path(OUT, sprintf("%s_Overlap_Matrix_Count.csv", PREFIX))
  out2 <- file.path(OUT, sprintf("%s_Overlap_Matrix_Pct.csv", PREFIX))
  write.csv(as.data.frame(overlap_mat), out1)
  write.csv(as.data.frame(overlap_pct), out2)
  cat(sprintf("  ✔ %s\n  ✔ %s\n", basename(out1), basename(out2)))
  return(overlap_mat)
}


# =============================================================================
# 模块 5: expression_divergence
# =============================================================================
run_expression_divergence <- function() {
  cat("--- Module: Expression Divergence ---\n")
  if (nrow(tg_lfc) == 0) { cat("  ⚠ Target LFC为空\n"); return(invisible(NULL)) }

  diff_long <- list()
  for (i in seq_len(nrow(tg_lfc))) {
    gid <- tg_lfc$Gene_ID[i]; cluster <- tg_lfc$WT_Cluster[i]; edge_dom <- tg_lfc$Edge_Dominant[i]
    for (tp in time_points) {
      wt_col <- paste0(WT,"_",tp); mut_col <- paste0(MUT,"_",tp)
      if (wt_col %in% names(tg_lfc) && mut_col %in% names(tg_lfc)) {
        wv <- as.numeric(tg_lfc[[wt_col]][i]); mv <- as.numeric(tg_lfc[[mut_col]][i])
        diff_long[[length(diff_long)+1]] <- data.frame(
          Gene_ID=gid, WT_Cluster=cluster, Edge_Dominant=edge_dom,
          Time=tp, WT_LFC=wv, MUT_LFC=mv, Diff=wv-mv, stringsAsFactors=FALSE)
      }
    }
  }
  diff_all <- bind_rows(diff_long)

  summary_df <- diff_all %>%
    group_by(WT_Cluster, Edge_Dominant, Time) %>%
    summarise(N=n(), Mean_Diff=mean(Diff,na.rm=TRUE), SD_Diff=sd(Diff,na.rm=TRUE),
              Mean_WT=mean(WT_LFC,na.rm=TRUE), Mean_MUT=mean(MUT_LFC,na.rm=TRUE), .groups="drop") %>%
    arrange(WT_Cluster, Edge_Dominant, Time)

  out1 <- file.path(OUT, sprintf("%s_Expression_Divergence.csv", PREFIX))
  out2 <- file.path(OUT, sprintf("%s_Expression_Divergence_ByGene.csv", PREFIX))
  write.csv(summary_df, out1, row.names=FALSE)
  write.csv(diff_all, out2, row.names=FALSE)
  cat(sprintf("  ✔ %s\n  ✔ %s\n", basename(out1), basename(out2)))

  last_tp <- time_points[length(time_points)]
  cat(sprintf("\n  [%s] Cluster × Edge mean(%s-%s):\n", last_tp, WT, MUT))
  last_s <- summary_df %>% filter(Time==last_tp)
  for (i in seq_len(nrow(last_s))) {
    r <- last_s[i,]
    cat(sprintf("    %s × %-10s (N=%2d): %.3f ± %.3f\n",
                r$WT_Cluster, r$Edge_Dominant, r$N,
                r$Mean_Diff, ifelse(is.na(r$SD_Diff),0,r$SD_Diff)))
  }
  return(summary_df)
}


# =============================================================================
# 模块 6: subspace
# =============================================================================
run_subspace <- function(genes) {
  cat("--- Module: Regulatory Subspace ---\n")
  tf_ids <- edges %>% group_by(Source) %>% summarise(Out=n(),.groups="drop") %>%
    filter(Out>=3) %>% pull(Source)
  cat(sprintf("  TFs (Out>=3): %d\n", length(tf_ids)))

  target_sets <- lapply(tf_ids, get_tf_targets); names(target_sets) <- tf_ids
  tf_dom <- sapply(tf_ids, function(tf) {
    ec <- edges %>% filter(Source==tf) %>% group_by(Edge_Class) %>% summarise(n=n(),.groups="drop")
    if(nrow(ec)==0) NA else ec$Edge_Class[which.max(ec$n)]
  })

  cons_tfs <- names(tf_dom[tf_dom=="Conserved" & !is.na(tf_dom)])
  lost_tfs <- names(tf_dom[tf_dom=="Lost" & !is.na(tf_dom)])
  gain_tfs <- names(tf_dom[tf_dom=="Gained" & !is.na(tf_dom)])
  cat(sprintf("  Conserved-dom: %d | Lost-dom: %d | Gained-dom: %d\n",
              length(cons_tfs), length(lost_tfs), length(gain_tfs)))

  pub  <- unique(unlist(target_sets[cons_tfs]))
  priv <- unique(unlist(target_sets[lost_tfs]))
  gain <- unique(unlist(target_sets[gain_tfs]))

  ov_pp <- length(intersect(pub, priv))
  ov_pg <- length(intersect(pub, gain))
  ov_rg <- length(intersect(priv, gain))

  cat(sprintf("  Public pool: %d | Private set: %d | Gained: %d\n",
              length(pub), length(priv), length(gain)))
  cat(sprintf("  Pub∩Priv: %d | Pub∩Gain: %d(%.0f%%) | Priv∩Gain: %d(%.0f%%)\n",
              ov_pp, ov_pg, ifelse(length(gain)>0, ov_pg/length(gain)*100, 0),
              ov_rg, ifelse(length(gain)>0, ov_rg/length(gain)*100, 0)))

  junct <- intersect(pub, priv)
  if (length(junct) > 0) {
    cat(sprintf("\n  交汇靶点 (%d):\n", length(junct)))
    for (jg in junct) {
      ge <- get_gene_edges(jg); deg <- nrow(ge); nl <- sum(ge$Edge_Class=="Lost")
      cat(sprintf("    %s: Deg=%d, Lost=%d(%.0f%%)\n",
                  fmt_gene(jg), deg, nl, ifelse(deg>0, nl/deg*100, 0)))
    }
  }

  rf <- file.path(OUT, sprintf("%s_Subspace_Report.txt", PREFIX))
  fc <- file(rf, "w")
  writeLines("调控子空间分析报告", fc)
  writeLines(sprintf("Generated: %s\n", Sys.time()), fc)
  writeLines(sprintf("保守型TFs (%d): %s", length(cons_tfs), paste(sapply(cons_tfs,fmt_gene),collapse=", ")), fc)
  writeLines(sprintf("Lost型TFs (%d): %s", length(lost_tfs), paste(sapply(lost_tfs,fmt_gene),collapse=", ")), fc)
  writeLines(sprintf("Gained型TFs (%d): %s", length(gain_tfs), paste(sapply(gain_tfs,fmt_gene),collapse=", ")), fc)
  writeLines(sprintf("\n公共池: %d | 私有集: %d | Gained覆盖: %d", length(pub), length(priv), length(gain)), fc)
  writeLines(sprintf("Pub∩Priv: %d | Pub∩Gain: %d | Priv∩Gain: %d", ov_pp, ov_pg, ov_rg), fc)
  if (length(junct)>0) { writeLines(sprintf("\n交汇靶点 (%d):", length(junct)), fc); writeLines(paste(junct,collapse="\n"), fc) }
  close(fc)

  write(pub, file.path(OUT, sprintf("%s_Subspace_Public_Pool.txt", PREFIX)))
  write(priv, file.path(OUT, sprintf("%s_Subspace_Private_Set.txt", PREFIX)))
  write(gain, file.path(OUT, sprintf("%s_Subspace_Gained_Coverage.txt", PREFIX)))
  cat(sprintf("  ✔ %s + gene lists\n", basename(rf)))
  return(invisible(NULL))
}


# =============================================================================
# 4. 主执行逻辑
# =============================================================================
mode <- tolower(opt$mode)
valid_modes <- c("gene_profile","topology_table","shared_targets","overlap_matrix",
                 "expression_divergence","subspace","full_report")
if (!mode %in% valid_modes) stop(sprintf("❌ 未知: %s\n可选: %s", opt$mode, paste(valid_modes,collapse=", ")))

if (mode == "gene_profile")            run_gene_profile(query_genes)
if (mode == "topology_table")          run_topology_table(opt$degree_cutoff)
if (mode == "shared_targets")          run_shared_targets(query_genes)
if (mode == "overlap_matrix")          run_overlap_matrix(query_genes)
if (mode == "expression_divergence")   run_expression_divergence()
if (mode == "subspace")                run_subspace(query_genes)
if (mode == "full_report") {
  cat("========== FULL REPORT ==========\n\n")
  run_topology_table(opt$degree_cutoff); cat("\n")
  run_expression_divergence(); cat("\n")
  run_subspace(query_genes); cat("\n")
  if (length(query_genes) >= 2) { run_shared_targets(query_genes[1:2]); cat("\n") }
  if (length(query_genes) > 0) { run_overlap_matrix(query_genes); cat("\n"); run_gene_profile(query_genes) }
  cat("\n========== COMPLETE ==========\n")
}

cat(sprintf("\n✅ Pipeline 8 完成 (%s)\n   输出: %s\n\n", opt$mode, OUT))