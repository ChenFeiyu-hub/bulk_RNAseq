#!/usr/bin/env Rscript
################################################################################
# RNAseq Pipeline 7c: GRN Network Visualization
#
# 版本: v4.4 (参数化 Degree 控制 + 截断后重剪枝)
#
# ============================================================================
# v4.4 更新日志 (基于 v4.3)
# ============================================================================
# 1. Degree-1 target 保留改为双参数控制:
#    - --network_keep_deg1_all   (默认 FALSE): 保留所有 TF 的 degree=1 target
#    - --network_keep_deg1_focus (默认 FALSE): 单独豁免 focus gene 的 degree=1 target
#    - 默认: 两者均 FALSE → 不显示任何 degree=1 target
#    - highlight 基因不再自动免检 degree 约束
#
# 2. Degree 剪枝逻辑提取为可复用函数 prune_degree1_targets():
#    - Step 1d:  初始剪枝 (anchor 边选取后)
#    - Step 1e2: Safety cap 截断后重新执行 (截断可能产生新的 degree=1 节点)
#
# 3. 处理链路: 1c→1d(degree)→1d2(ego)→1e(cap)→1e2(degree重)→1e3(连通性)
#
# ============================================================================
# v4.3 更新日志 (基于 v4.2)
# ============================================================================
# 1. 参数重命名: --network_max_nodes → --network_max_edges (默认 600)
#    - 直接控制边数上限, 去掉原先 ×3 的经验换算
#    - bash 层面对应: GRN_NETWORK_MAX_EDGES
#
# 2. Safety Cap 后拓扑连通性重校验 (Step 1e2):
#    - 截断边后重新从 focus_gene 计算可达性
#    - 移除与 focus_gene 断联的孤岛节点及其关联边
#    - 确保 Fig4 网络中所有节点与 focus_gene 保持拓扑连接
#
# ============================================================================
# v4.2 更新日志 (基于 v4.1)
# ============================================================================
# 1. 新增 --network_ego_steps 参数:
#    - 0   = 仅去除与 focus_gene 不在同一连通分量的孤岛节点
#    - N>0 = 仅保留 focus_gene N 步拓扑距离以内的节点 (Ego-Network)
#    - 该过滤在 Part 1 的 Smart Pruning 之后执行 (Step 1d2)
#
# 2. 适用于所有模式:
#    - 用户列表模式 (highlight_file): 默认 ego_steps=0 → 去除孤岛
#    - 自动计算模式 (7d clusters):    默认 ego_steps=2 → 2-hop Ego
#    - bash 层面可通过 GRN_NETWORK_EGO_STEPS 配置覆盖
#
# 3. 前向兼容: 未设 focus_gene 时自动跳过 ego 过滤
#
# ============================================================================
# 历史版本
# ============================================================================
# v4.1: Directional boxes, Fig7a/7b 布局优化
# v3.0: Evidence Heatmap + Conserved Edge Fix
#
# ============================================================================
# 图表概览
# ============================================================================
# 1. heatmap          - TF Family × Cluster 富集热图
# 2. volcano          - Volcano Plot: delta_r vs -log10(FDR)
# 3. barplot          - Cluster-wise Edge Class 分布图
# 4. network          - Focus TF Ego Network (tidygraph + ggraph)
# 5. circos           - Circos: TF Family → Cluster
# 6. outdegree        - Delta Out-Degree Bar Plot
# 7. evidence_heatmap - 双热图+Bezier连接: Target×TF 多维证据图
#
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(jsonlite)
})

################################################################################
# 参数解析
################################################################################

parse_arguments <- function() {
  args <- commandArgs(trailingOnly = TRUE)

  params <- list(
    rdata_7b              = NULL,
    out_dir               = NULL,
    fig_format            = "both",
    fig_dpi               = 300,
    fig_scale             = 1.0,
    plots                 = "heatmap,volcano,barplot,network,circos,outdegree,evidence_heatmap",
    focus_tf_family       = "inherit",
    focus_gene            = "inherit",
    condition_wt          = "WT",
    condition_mut         = "dof",
    # Fig4 网络图参数
    network_max_edges     = 600,
    network_layout        = "kk",
    network_seed          = 42,
    network_prune         = TRUE,
    network_fig_w         = 18,
    network_fig_h         = 14,
    network_arrow_mm      = 1.5,
    network_edge_width_min = 0.5,
    network_edge_width_max = 2.0,
    network_node_size_min = 4,
    network_node_size_max = 12,
    # ★ v4.2: Ego-Network 拓扑步数
    network_ego_steps     = 0,
    # ★ v4.4: Degree-1 target 保留控制
    network_keep_deg1_all   = FALSE,   # 保留所有TF的 degree=1 target
    network_keep_deg1_focus = FALSE,   # 单独豁免 focus gene 的 degree=1 target
    # Fig5 Circos
    circos_max_edges      = 500,
    circos_top_edges      = 300,
    # 标签
    label_top_n           = 15,
    label_focus           = TRUE,
    # Fig1
    show_heatmap_stars    = FALSE,
    # 颜色
    cluster_colors        = NULL,
    edge_class_colors     = NULL,
    # DESeq2 差异热图参数
    deseq2_rdata          = NULL,
    control_time          = "00h",
    cond_colors           = "WT:#2E86AB,dof:#A23B72",
    time_colors           = "00h:#E5D65C,03h:#E9B730,06h:#E9992C,12h:#E47B43,24h:#DA595B",
    # Highlight 支持
    highlight_file        = NULL,
    # =========================
    verbose               = TRUE
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (grepl("^--", arg)) {
      key <- gsub("-", "_", sub("^--", "", arg))
      if (grepl("=", key)) {
        parts <- strsplit(key, "=")[[1]]
        key <- parts[1]; value <- paste(parts[-1], collapse = "=")
      } else if (i + 1 <= length(args) && !grepl("^--", args[i + 1])) {
        i <- i + 1; value <- args[i]
      } else {
        value <- "TRUE"
      }
      if (key %in% names(params)) {
        if (is.logical(params[[key]])) {
          params[[key]] <- as.logical(value)
        } else if (is.numeric(params[[key]]) && !is.null(params[[key]])) {
          params[[key]] <- as.numeric(value)
        } else {
          params[[key]] <- value
        }
      }
    }
    i <- i + 1
  }

  params$plots_vec <- trimws(unlist(strsplit(params$plots, ",")))

  # v3.0: 向后兼容 - "alluvial" → "evidence_heatmap"
  params$plots_vec <- gsub("^alluvial$", "evidence_heatmap", params$plots_vec)

  return(params)
}

PARAMS <- parse_arguments()

if (is.null(PARAMS$rdata_7b)) {
  cat("错误: 必须指定 --rdata_7b\n")
  quit(status = 1)
}

################################################################################
# 辅助函数
################################################################################

# --- 保存图表 (双格式) ---
save_plot <- function(plot_obj, filename, out_dir, width, height,
                       fig_format = "both", dpi = 300, scale = 1.0,
                       device_fn = NULL) {
  w <- width * scale
  h <- height * scale

  if (fig_format %in% c("pdf", "both")) {
    if (!is.null(device_fn)) {
      pdf(file.path(out_dir, paste0(filename, ".pdf")), width = w, height = h)
      device_fn()
      dev.off()
    } else {
      ggsave(file.path(out_dir, paste0(filename, ".pdf")),
             plot = plot_obj, width = w, height = h, device = "pdf")
    }
  }
  if (fig_format %in% c("png", "both")) {
    if (!is.null(device_fn)) {
      png(file.path(out_dir, paste0(filename, ".png")),
          width = w, height = h, units = "in", res = dpi)
      device_fn()
      dev.off()
    } else {
      ggsave(file.path(out_dir, paste0(filename, ".png")),
             plot = plot_obj, width = w, height = h, dpi = dpi, device = "png")
    }
  }
  cat(sprintf("    Saved: %s (.%s)\n", filename,
              ifelse(fig_format == "both", "pdf + .png", fig_format)))
}

# --- 颜色字符串解析 ---
parse_color_string <- function(s) {
  if (is.null(s) || s == "") return(NULL)
  pairs <- trimws(unlist(strsplit(s, ",")))
  result <- list()
  for (p in pairs) {
    kv <- strsplit(p, ":")[[1]]
    if (length(kv) == 2) {
      result[[trimws(kv[1])]] <- trimws(kv[2])
    }
  }
  return(unlist(result))
}

# --- 缩短 Gene ID ---
shorten_id <- function(x) sub("_.*$", "", x)

# --- ★ v4.4: Degree-based target pruning (可复用) ---
# 该函数在初始剪枝和 safety cap 后均被调用.
# Essential 节点 = 所有 TF + focus_gene, highlight 基因不再自动免检.
# effective_degree 仅基于 differential edges (排除 Conserved).
prune_degree1_targets <- function(net_edges, tf_family_map,
                                    focus_gene_id, focus_tf_ids,
                                    keep_deg1_all, keep_deg1_focus,
                                    label = "") {

  all_nodes <- unique(c(net_edges$TF_Gene, net_edges$Target_Gene))
  all_tf_ids <- intersect(all_nodes, names(tf_family_map))

  # Essential: TFs + focus_gene (highlight 不再免检)
  essential <- unique(c(all_tf_ids, focus_tf_ids))
  if (!is.null(focus_gene_id) && !is.na(focus_gene_id) && focus_gene_id != "none") {
    essential <- unique(c(essential, focus_gene_id))
  }

  # Effective degree: 仅 differential edges
  diff_df <- net_edges %>% filter(edge_class != "Conserved")
  if (nrow(diff_df) > 0) {
    g_diff <- igraph::graph_from_data_frame(
      diff_df %>% select(TF_Gene, Target_Gene), directed = TRUE
    )
    eff_deg <- igraph::degree(g_diff, mode = "all")
  } else {
    eff_deg <- integer(0)
  }

  # Focus gene 的直接邻居 (全边集, 含 Conserved)
  focus_nbrs <- character(0)
  if (!is.null(focus_gene_id) && !is.na(focus_gene_id) && focus_gene_id != "none") {
    focus_nbrs <- unique(c(
      net_edges$Target_Gene[net_edges$TF_Gene == focus_gene_id],
      net_edges$TF_Gene[net_edges$Target_Gene == focus_gene_id]
    ))
  }

  keep <- character(0)
  for (nid in all_nodes) {
    if (nid %in% essential) {
      keep <- c(keep, nid)
    } else {
      d <- eff_deg[nid]
      if (!is.na(d) && d > 1) {
        keep <- c(keep, nid)
      } else if (keep_deg1_all) {
        keep <- c(keep, nid)
      } else if (keep_deg1_focus && nid %in% focus_nbrs) {
        keep <- c(keep, nid)
      }
      # else: degree ≤ 1 且无豁免 → 剪除
    }
  }
  keep <- unique(keep)

  n_before <- length(all_nodes)
  n_pruned <- n_before - length(keep)

  out_edges <- net_edges %>%
    filter(TF_Gene %in% keep & Target_Gene %in% keep)

  cat(sprintf("       [%s] Degree prune: %d → %d nodes (-%d) | %d edges\n",
              label, n_before, length(keep), n_pruned, nrow(out_edges)))

  out_edges
}

################################################################################
# 主程序
################################################################################

cat("================================================================================\n")
cat("  Pipeline 7c: GRN Network Visualization (v4.4)\n")
cat("================================================================================\n\n")
cat(sprintf("开始时间: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

# ========================================================================
# Part 0: 加载数据
# ========================================================================

cat("--- Part 0: 加载 Pipeline 7b 数据 ---\n")

PARAMS_7c <- PARAMS
load(PARAMS$rdata_7b)
cat(sprintf("  Loaded 7b RData: %s\n", basename(PARAMS_7c$rdata_7b)))
PARAMS <- PARAMS_7c

# ========================================================================
# ★ 读取 7c bash 传来的动态 highlight_file 并覆盖 RData 中的旧值
# ========================================================================
if (!is.null(PARAMS$highlight_file) && file.exists(PARAMS$highlight_file)) {
  cat("\n  [Update] Loading external highlight file generated by 7d...\n")
  hl_df <- read.table(PARAMS$highlight_file, header = FALSE, sep = "\t",
                      stringsAsFactors = FALSE, quote = "")
  
  # 自动识别是否有表头
  if (grepl("^(Symbol|Gene)", hl_df[1,1], ignore.case=TRUE)) {
    colnames(hl_df) <- hl_df[1, ]
    hl_df <- hl_df[-1, ]
  } else {
    colnames(hl_df) <- c("Symbol", "Gene")
  }
  
  # 覆盖从 7b 继承的 highlight_map
  highlight_map <- setNames(hl_df$Symbol, hl_df$Gene)
  cat(sprintf("  Highlight genes updated to: %d genes\n", length(highlight_map)))
}

# 继承 focus 参数
if (PARAMS$focus_tf_family == "inherit" && exists("PARAMS_7a_inherited")) {
  inherited_params <- get("PARAMS", envir = .GlobalEnv)
  if (!is.null(inherited_params$focus_tf_family)) {
    PARAMS$focus_tf_family <- inherited_params$focus_tf_family
  }
}
if (PARAMS$focus_gene == "inherit" && exists("PARAMS_7a_inherited")) {
  inherited_params <- get("PARAMS", envir = .GlobalEnv)
  if (!is.null(inherited_params$focus_gene)) {
    PARAMS$focus_gene <- inherited_params$focus_gene
  }
}

if (PARAMS$focus_tf_family == "inherit") PARAMS$focus_tf_family <- "none"
if (PARAMS$focus_gene == "inherit") PARAMS$focus_gene <- "none"

focus_family_vec <- if (PARAMS$focus_tf_family != "none") {
  trimws(unlist(strsplit(PARAMS$focus_tf_family, ",")))
} else { character(0) }

# 输出目录
if (is.null(PARAMS$out_dir)) {
  PARAMS$out_dir <- dirname(dirname(PARAMS$rdata_7b))
}
OUT_7C <- file.path(PARAMS$out_dir, "Pipeline7c")
dir.create(OUT_7C, recursive = TRUE, showWarnings = FALSE)

# 日志
log_file <- file.path(OUT_7C, "Pipeline7c_visualization.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output", split = TRUE)
sink(log_con, type = "message")

cat(sprintf("  Plots to generate: %s\n", paste(PARAMS$plots_vec, collapse = ", ")))
cat(sprintf("  Focus family: %s | Focus gene: %s\n",
            PARAMS$focus_tf_family, PARAMS$focus_gene))
cat(sprintf("  v4.4 params: ego_steps=%d, max_edges=%d, keep_deg1_all=%s, keep_deg1_focus=%s\n",
            PARAMS$network_ego_steps, PARAMS$network_max_edges,
            PARAMS$network_keep_deg1_all, PARAMS$network_keep_deg1_focus))
cat(sprintf("               show_heatmap_stars=%s, circos_top_edges=%d\n",
            PARAMS$show_heatmap_stars, PARAMS$circos_top_edges))
cat(sprintf("  v4.0 new: deseq2_rdata=%s\n",
            ifelse(is.null(PARAMS$deseq2_rdata), "NULL", basename(PARAMS$deseq2_rdata))))

# ========================================================================
# Part 0b: 颜色方案初始化
# ========================================================================

cat("\n--- Part 0b: 初始化颜色方案 ---\n")

init_colors <- function(n_clusters, n_families, focus_fams,
                         cluster_colors_file, edge_class_colors_file) {

  colors <- list()

  if (!is.null(cluster_colors_file) && file.exists(cluster_colors_file)) {
    cc <- read.table(cluster_colors_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    colors$cluster <- setNames(cc$Color, cc$Cluster)
  } else {
    if (n_clusters <= 8) {
      pal <- RColorBrewer::brewer.pal(max(3, n_clusters), "Set2")
    } else if (n_clusters <= 12) {
      pal <- RColorBrewer::brewer.pal(n_clusters, "Set3")
    } else {
      pal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(n_clusters)
    }
    colors$cluster <- setNames(pal[1:n_clusters], paste0("C", 1:n_clusters))
  }

  if (!is.null(edge_class_colors_file) && file.exists(edge_class_colors_file)) {
    ec <- read.table(edge_class_colors_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    colors$edge_class <- setNames(ec$Color, ec$EdgeClass)
  } else {
    colors$edge_class <- c(
      Lost = "#E74C3C",
      Gained = "#27AE60",
      Rewired = "#8E44AD",
      Conserved = "#3498DB",
      WT_only = "#95A5A6",
      MUT_only = "#BDC3C7",
      NS = "#ECF0F1"
    )
  }

  if (n_families <= 20) {
    fam_pal <- scales::hue_pal()(n_families)
  } else {
    fam_pal <- colorRampPalette(RColorBrewer::brewer.pal(8, "Dark2"))(n_families)
  }
  colors$family <- fam_pal

  colors$focus <- "#FFD700"
  colors$direction <- c(positive = "#C0392B", negative = "#2980B9")

  colors$node_type <- c(
    "Focus TF"    = "#E63946",
    "Other TF"    = "#457B9D",
    "Target"      = "#F4A261",
    "Highlight"   = "#2A9D8F"
  )

  colors$pcre <- c(
    "pCRE_Pos" = "#2ECC71",
    "pCRE_Neg" = "#E67E22",
    "absent"   = "#F0F0F0"
  )

  colors$condition <- parse_color_string(PARAMS$cond_colors)
  colors$time      <- parse_color_string(PARAMS$time_colors)

  return(colors)
}

n_fam <- length(unique(tf_family_map))
colors <- init_colors(
  n_clusters = max(1, n_wt_clusters),
  n_families = n_fam,
  focus_fams = focus_family_vec,
  cluster_colors_file = PARAMS$cluster_colors,
  edge_class_colors_file = PARAMS$edge_class_colors
)

cat(sprintf("  Cluster colors: %d | Family colors: %d | Edge class colors: %d\n",
            length(colors$cluster), n_fam, length(colors$edge_class)))
cat(sprintf("  Condition colors: %s\n", paste(names(colors$condition), colors$condition, sep="=", collapse=", ")))
cat(sprintf("  Time colors: %s\n", paste(names(colors$time), colors$time, sep="=", collapse=", ")))

# 图表计数器
plot_results <- list()

# ========================================================================
# Part 1: 共享网络数据构建 (v4.4: 参数化 Degree + 截断重剪枝)
# ========================================================================

need_network_data <- any(c("network", "evidence_heatmap") %in% PARAMS$plots_vec)

network_data <- NULL

if (need_network_data) {
  cat("\n--- Part 1: 构建共享网络数据 (v4.4: 参数化 Degree + 截断重剪枝) ---\n")

  network_data <- tryCatch({

    # --- Step 1a: 确定边数据来源 ---
    if (exists("edge_table_scored") && is.data.frame(edge_table_scored) &&
        nrow(edge_table_scored) > 0) {
      edge_source <- edge_table_scored
      cat("    Edge source: edge_table_scored\n")
    } else if (exists("scored_edges") && is.data.frame(scored_edges) &&
               nrow(scored_edges) > 0) {
      edge_source <- scored_edges
      cat("    Edge source: scored_edges (fallback)\n")
    } else {
      stop("No edge table found in RData")
    }

    # --- Step 1b: 定义锚定基因集 ---
    anchor_focus_tfs <- if (exists("focus_tfs") && length(focus_tfs) > 0) {
      focus_tfs
    } else { character(0) }

    anchor_focus_gene <- if (PARAMS$focus_gene != "none" && !is.na(PARAMS$focus_gene)) {
      PARAMS$focus_gene
    } else { character(0) }

    anchor_highlight <- if (!is.null(highlight_map) && length(highlight_map) > 0) {
      names(highlight_map)
    } else { character(0) }

    anchor_genes <- unique(c(anchor_focus_tfs, anchor_focus_gene, anchor_highlight))

    if (length(anchor_genes) == 0) {
      stop("No anchor genes defined.")
    }

    cat(sprintf("    Anchor genes: %d (Focus TFs=%d, Focus gene=%d, Highlight=%d)\n",
                length(anchor_genes), length(anchor_focus_tfs),
                length(anchor_focus_gene), length(anchor_highlight)))

    # --- Step 1c: 锚定边选取 ---
    net_edges <- edge_source %>%
      filter(
        (TF_Gene %in% anchor_genes | Target_Gene %in% anchor_genes),
        edge_class %in% c("Lost", "Gained", "Rewired", "Conserved"),
        (sig_WT | sig_MUT)
      )

    cat(sprintf("    Anchor-based edges (raw): %d\n", nrow(net_edges)))

    if (nrow(net_edges) == 0) {
      stop("No significant edges found for anchor genes.")
    }

    raw_class_counts <- table(net_edges$edge_class)
    cat(sprintf("       Lost: %d | Gained: %d | Rewired: %d | Conserved: %d\n",
                raw_class_counts["Lost"] %||% 0,
                raw_class_counts["Gained"] %||% 0,
                raw_class_counts["Rewired"] %||% 0,
                raw_class_counts["Conserved"] %||% 0))

    # --- Step 1d: 智能剪枝 (v4.4: 参数化 degree 控制) ---
    if (PARAMS$network_prune) {
      cat("    Smart pruning (v4.4: parameterized degree control)...\n")
      cat(sprintf("       keep_deg1_all=%s | keep_deg1_focus=%s\n",
                  PARAMS$network_keep_deg1_all, PARAMS$network_keep_deg1_focus))

      net_edges <- prune_degree1_targets(
        net_edges       = net_edges,
        tf_family_map   = tf_family_map,
        focus_gene_id   = PARAMS$focus_gene,
        focus_tf_ids    = anchor_focus_tfs,
        keep_deg1_all   = PARAMS$network_keep_deg1_all,
        keep_deg1_focus = PARAMS$network_keep_deg1_focus,
        label           = "Step1d"
      )

      if (nrow(net_edges) == 0) {
        stop("All edges removed after pruning.")
      }
    }

    # =================================================================
    # ★ Step 1d2: Ego-Network 拓扑约束 (v4.2 新增)
    # =================================================================
    #
    # 逻辑:
    #   ego_steps = 0 : 仅去除与 focus_gene 不在同一连通分量的节点 (去孤岛)
    #   ego_steps > 0 : 仅保留 focus_gene N 步拓扑距离以内的节点
    #   focus_gene = "none" 或不在网络中: 跳过
    #
    # 使用无向图计算最短路径距离, 确保双向可达性
    # =================================================================

    ego_steps <- as.integer(PARAMS$network_ego_steps)
    focus_gene_for_ego <- PARAMS$focus_gene

    ego_filter_applied <- FALSE
    ego_n_removed <- 0L

    if (!is.null(focus_gene_for_ego) &&
        focus_gene_for_ego != "none" &&
        !is.na(focus_gene_for_ego)) {

      all_node_ids_pre_ego <- unique(c(net_edges$TF_Gene, net_edges$Target_Gene))

      if (focus_gene_for_ego %in% all_node_ids_pre_ego) {

        cat(sprintf("    ★ Ego-network filter: center=%s, max_steps=%s\n",
                    focus_gene_for_ego,
                    ifelse(ego_steps > 0,
                           as.character(ego_steps),
                           "\u221E (connected component only)")))

        # 构建无向图用于距离计算
        g_ego_temp <- igraph::graph_from_data_frame(
          net_edges %>% select(TF_Gene, Target_Gene), directed = TRUE
        )
        g_ego_undir <- igraph::as.undirected(g_ego_temp, mode = "collapse")

        # 计算 focus_gene 到所有节点的最短路径距离
        dists_from_focus <- igraph::distances(
          g_ego_undir, v = focus_gene_for_ego, mode = "all"
        )[1, ]

        if (ego_steps > 0) {
          # N-step ego: 仅保留距离 ≤ N 的节点
          ego_keep <- names(dists_from_focus)[
            !is.na(dists_from_focus) &
            is.finite(dists_from_focus) &
            dists_from_focus <= ego_steps
          ]
          cat(sprintf("       %d-step ego: %d / %d nodes reachable within %d hops\n",
                      ego_steps, length(ego_keep),
                      length(all_node_ids_pre_ego), ego_steps))
        } else {
          # 连通分量过滤: 保留所有可达节点 (距离有限)
          ego_keep <- names(dists_from_focus)[
            !is.na(dists_from_focus) & is.finite(dists_from_focus)
          ]
          cat(sprintf("       Connected component: %d / %d nodes reachable\n",
                      length(ego_keep), length(all_node_ids_pre_ego)))
        }

        n_before_ego <- length(all_node_ids_pre_ego)

        # 过滤边: 两端都必须在 ego_keep 中
        net_edges <- net_edges %>%
          filter(TF_Gene %in% ego_keep & Target_Gene %in% ego_keep)

        n_after_ego <- length(unique(c(net_edges$TF_Gene, net_edges$Target_Gene)))
        ego_n_removed <- n_before_ego - n_after_ego
        ego_filter_applied <- TRUE

        cat(sprintf("       Nodes: %d \u2192 %d (removed %d disconnected/distant)\n",
                    n_before_ego, n_after_ego, ego_n_removed))
        cat(sprintf("       Edges after ego filter: %d\n", nrow(net_edges)))

        if (nrow(net_edges) == 0) {
          stop(sprintf(
            "All edges removed after ego filter (ego_steps=%d). %s",
            ego_steps,
            ifelse(ego_steps > 0,
                   "Try increasing --network_ego_steps.",
                   "focus_gene may be isolated in the network.")
          ))
        }

        # 打印距离分布概要
        reachable_dists <- dists_from_focus[is.finite(dists_from_focus)]
        if (length(reachable_dists) > 0) {
          dist_tab <- table(reachable_dists)
          cat("       Distance distribution from focus_gene:\n")
          for (d in sort(as.integer(names(dist_tab)))) {
            cat(sprintf("         %d hop(s): %d nodes\n", d, dist_tab[as.character(d)]))
          }
        }

      } else {
        cat(sprintf("    \u26A0 focus_gene '%s' not found in network nodes, ego filter skipped\n",
                    focus_gene_for_ego))
      }
    } else {
      cat("    \u26A0 No valid focus_gene specified (='none'), ego filter skipped\n")
    }

    # --- Step 1e: 安全截断 ---
    max_edges <- PARAMS$network_max_edges

    if (nrow(net_edges) > max_edges) {
      cat(sprintf("    Safety cap: %d edges > limit %d\n", nrow(net_edges), max_edges))

      anchor_edge_mask <- (net_edges$TF_Gene %in% anchor_genes |
                           net_edges$Target_Gene %in% anchor_genes)
      anchor_edges_df <- net_edges[anchor_edge_mask, ]
      other_edges_df  <- net_edges[!anchor_edge_mask, ]

      if (nrow(anchor_edges_df) < max_edges) {
        n_other_keep <- max_edges - nrow(anchor_edges_df)
        other_edges_df <- other_edges_df %>%
          arrange(desc(abs(MECS_score))) %>%
          head(n_other_keep)
        net_edges <- bind_rows(anchor_edges_df, other_edges_df)
      } else {
        net_edges <- anchor_edges_df %>%
          arrange(desc(abs(MECS_score))) %>%
          head(max_edges)
      }

      cat(sprintf("       Truncated to %d edges\n", nrow(net_edges)))
    }

    # =================================================================
    # ★ Step 1e2: 截断后 Degree 重剪枝 (v4.4 新增)
    # =================================================================
    # Safety cap 截断可能使原来 degree>1 的节点降为 degree=1,
    # 产生新的单线伞状结构. 此步骤用相同参数重新执行 degree 剪枝.
    # =================================================================

    if (PARAMS$network_prune) {
      n_edges_before_reprune <- nrow(net_edges)

      net_edges <- prune_degree1_targets(
        net_edges       = net_edges,
        tf_family_map   = tf_family_map,
        focus_gene_id   = PARAMS$focus_gene,
        focus_tf_ids    = anchor_focus_tfs,
        keep_deg1_all   = PARAMS$network_keep_deg1_all,
        keep_deg1_focus = PARAMS$network_keep_deg1_focus,
        label           = "Step1e2-postCap"
      )

      if (nrow(net_edges) < n_edges_before_reprune) {
        cat(sprintf("       Post-cap degree reprune: %d → %d edges\n",
                    n_edges_before_reprune, nrow(net_edges)))
      }
    }

    # =================================================================
    # ★ Step 1e3: 截断后拓扑连通性重校验 (v4.3)
    # =================================================================
    # Safety cap + degree 重剪枝可能切断 focus_gene 的间接路径,
    # 导致出现与 focus_gene 无拓扑连接的孤岛子图.
    # =================================================================

    if (!is.null(focus_gene_for_ego) &&
        focus_gene_for_ego != "none" &&
        !is.na(focus_gene_for_ego)) {

      post_cap_nodes <- unique(c(net_edges$TF_Gene, net_edges$Target_Gene))

      if (focus_gene_for_ego %in% post_cap_nodes) {
        g_post_cap <- igraph::graph_from_data_frame(
          net_edges %>% select(TF_Gene, Target_Gene), directed = TRUE
        )
        g_post_cap_undir <- igraph::as.undirected(g_post_cap, mode = "collapse")

        dists_post <- igraph::distances(
          g_post_cap_undir, v = focus_gene_for_ego, mode = "all"
        )[1, ]

        reachable_post <- names(dists_post)[
          !is.na(dists_post) & is.finite(dists_post)
        ]

        n_disconnected <- length(post_cap_nodes) - length(reachable_post)

        if (n_disconnected > 0) {
          cat(sprintf("    ★ Step1e3 connectivity check: %d nodes disconnected from %s, removing...\n",
                      n_disconnected, focus_gene_for_ego))

          net_edges <- net_edges %>%
            filter(TF_Gene %in% reachable_post & Target_Gene %in% reachable_post)

          cat(sprintf("       Edges after reconnect filter: %d\n", nrow(net_edges)))
        } else {
          cat("    ★ Step1e3 connectivity check: all nodes reachable, OK\n")
        }
      }
    }

    # --- Step 1f: 构建节点元数据 ---
    all_node_ids <- unique(c(net_edges$TF_Gene, net_edges$Target_Gene))
    highlight_ids <- if (!is.null(highlight_map)) names(highlight_map) else character(0)
    highlight_symbols <- if (!is.null(highlight_map)) highlight_map else setNames(character(0), character(0))

    node_meta <- data.frame(id = all_node_ids, stringsAsFactors = FALSE) %>%
      mutate(
        Is_Focus_TF   = id %in% anchor_focus_tfs,
        Is_TF         = id %in% names(tf_family_map),
        Is_Highlight  = id %in% highlight_ids,
        Is_Focus_Gene = (id == PARAMS$focus_gene),
        Is_Anchor     = id %in% anchor_genes,
        Family        = ifelse(Is_TF, tf_family_map[id], NA_character_),
        Symbol        = ifelse(id %in% names(highlight_symbols),
                               highlight_symbols[id], NA_character_),
        WT_Cluster    = NA_integer_
      )

    # 填充 Cluster 信息
    if (exists("cluster_annotation") && !is.null(cluster_annotation)) {
      ca_match <- match(node_meta$id, cluster_annotation$Gene_ID)
      valid <- !is.na(ca_match)
      node_meta$WT_Cluster[valid] <- cluster_annotation$WT_Cluster[ca_match[valid]]
    }

    # 从边表提取 Symbol
    if ("TF_Symbol" %in% colnames(edge_source)) {
      tf_sym_map <- edge_source %>%
        filter(!is.na(TF_Symbol) & TF_Symbol != "") %>%
        distinct(TF_Gene, TF_Symbol) %>%
        { setNames(.$TF_Symbol, .$TF_Gene) }
    } else {
      tf_sym_map <- setNames(character(0), character(0))
    }

    if ("Target_Symbol" %in% colnames(edge_source)) {
      tg_sym_map <- edge_source %>%
        filter(!is.na(Target_Symbol) & Target_Symbol != "") %>%
        distinct(Target_Gene, Target_Symbol) %>%
        { setNames(.$Target_Symbol, .$Target_Gene) }
    } else {
      tg_sym_map <- setNames(character(0), character(0))
    }

    for (i in seq_len(nrow(node_meta))) {
      if (is.na(node_meta$Symbol[i]) || node_meta$Symbol[i] == "") {
        nid <- node_meta$id[i]
        if (nid %in% names(tf_sym_map)) {
          node_meta$Symbol[i] <- tf_sym_map[nid]
        } else if (nid %in% names(tg_sym_map)) {
          node_meta$Symbol[i] <- tg_sym_map[nid]
        }
      }
    }

    # 节点类型
    node_meta <- node_meta %>%
      mutate(
        NodeType = case_when(
          Is_Focus_TF  ~ "Focus TF",
          Is_TF        ~ "Other TF",
          Is_Highlight ~ "Highlight",
          TRUE         ~ "Target"
        ),
        Label = case_when(
          Is_Focus_TF & !is.na(Symbol) & Symbol != "" ~
            paste0(Symbol, " (", Family, ")"),
          Is_Focus_TF & !is.na(Family) ~
            paste0(shorten_id(id), " (", Family, ")"),
          Is_Focus_TF ~ shorten_id(id),
          Is_Highlight & !is.na(Symbol) & Symbol != "" ~ Symbol,
          Is_Highlight ~ shorten_id(id),
          Is_Focus_Gene & !is.na(Symbol) & Symbol != "" ~ Symbol,
          Is_Focus_Gene ~ shorten_id(id),
          Is_TF & !is.na(Symbol) & Symbol != "" ~ Symbol,
          Is_TF & !is.na(Family) ~ paste0(shorten_id(id), "\n(", Family, ")"),
          Is_TF ~ shorten_id(id),
          TRUE ~ ""
        )
      )

    node_meta$NodeType <- factor(
      node_meta$NodeType,
      levels = c("Focus TF", "Other TF", "Highlight", "Target")
    )

    cat(sprintf("    Network built: %d nodes, %d edges\n",
                nrow(node_meta), nrow(net_edges)))
    cat(sprintf("       Focus TF: %d | Other TF: %d | Highlight: %d | Target: %d\n",
                sum(node_meta$NodeType == "Focus TF"),
                sum(node_meta$NodeType == "Other TF"),
                sum(node_meta$NodeType == "Highlight"),
                sum(node_meta$NodeType == "Target")))

    # ★ v4.2: 将 ego 过滤信息附加到返回结构中
    list(
      net_edges           = net_edges,
      node_meta           = node_meta,
      anchor_genes        = anchor_genes,
      anchor_focus_tfs    = anchor_focus_tfs,
      anchor_focus_gene   = anchor_focus_gene,
      anchor_highlight    = anchor_highlight,
      edge_source         = edge_source,
      tf_sym_map          = tf_sym_map,
      tg_sym_map          = tg_sym_map,
      highlight_symbols   = highlight_symbols,
      ego_filter_applied  = ego_filter_applied,
      ego_steps_used      = ego_steps,
      ego_n_removed       = ego_n_removed
    )
  }, error = function(e) {
    cat(sprintf("    \u26A0 Network data construction failed: %s\n", e$message))
    NULL
  })
}

# ========================================================================
# 图 1: TF Family × Cluster 富集热图
# ========================================================================

if ("heatmap" %in% PARAMS$plots_vec) {
  cat("\n--- 图 1: Enrichment Heatmap ---\n")
  cat(sprintf("    show_heatmap_stars = %s\n", PARAMS$show_heatmap_stars))

  plot_results[["heatmap"]] <- tryCatch({
    suppressPackageStartupMessages(library(ComplexHeatmap))
    suppressPackageStartupMessages(library(circlize))

    if (is.null(enrichment_long) || nrow(enrichment_long) == 0) {
      stop("No enrichment data available")
    }

    hm_data <- enrichment_long %>%
      filter(Edge_Class == "Lost") %>%
      select(TF_Family, WT_Cluster, log2OR, fisher_fdr, sig_mark) %>%
      filter(!is.na(log2OR))

    if (nrow(hm_data) == 0) stop("No Lost enrichment data")

    mat <- hm_data %>%
      select(TF_Family, WT_Cluster, log2OR) %>%
      pivot_wider(names_from = WT_Cluster, values_from = log2OR) %>%
      column_to_rownames("TF_Family") %>%
      as.matrix()

    max_abs <- apply(abs(mat), 1, max, na.rm = TRUE)
    keep_rows <- max_abs >= 0.5
    if (sum(keep_rows) < 2) keep_rows <- rep(TRUE, nrow(mat))
    mat <- mat[keep_rows, , drop = FALSE]

    sig_mat <- NULL
    if (PARAMS$show_heatmap_stars) {
      sig_mat <- hm_data %>%
        select(TF_Family, WT_Cluster, sig_mark) %>%
        pivot_wider(names_from = WT_Cluster, values_from = sig_mark) %>%
        column_to_rownames("TF_Family") %>%
        as.matrix()
      sig_mat <- sig_mat[rownames(mat), colnames(mat), drop = FALSE]
      sig_mat[is.na(sig_mat)] <- ""
    }

    max_abs_val <- max(abs(mat), na.rm = TRUE)
    color_cap <- max(2, ceiling(max_abs_val))
    col_fun <- colorRamp2(
      c(-color_cap, 0, color_cap),
      c("#2166AC", "white", "#B2182B")
    )

    row_labels <- rownames(mat)
    row_highlight <- row_labels %in% focus_family_vec

    if (!is.null(cluster_pair_info)) {
      div_vals <- setNames(cluster_pair_info$Divergence_Ratio, cluster_pair_info$WT_Cluster)
      col_div <- sapply(colnames(mat), function(x) {
        dv <- div_vals[as.character(x)]
        if (is.na(dv)) 0 else dv
      })
      div_col_fun <- colorRamp2(c(0, 0.5, 1), c("#2ECC71", "#F1C40F", "#E74C3C"))
      col_ha <- HeatmapAnnotation(
        Divergence = col_div,
        col = list(Divergence = div_col_fun),
        annotation_name_side = "left"
      )
    } else {
      col_ha <- NULL
    }

    fig_w <- max(8, ncol(mat) * 0.8 + 4)
    fig_h <- max(6, nrow(mat) * 0.4 + 2)

    draw_heatmap <- function() {
      cell_fun_final <- NULL
      if (PARAMS$show_heatmap_stars && !is.null(sig_mat)) {
        cell_fun_final <- function(j, i, x, y, width, height, fill) {
          if (!is.na(sig_mat[i, j]) && sig_mat[i, j] != "") {
            grid::grid.text(sig_mat[i, j], x, y,
                            gp = grid::gpar(fontsize = 10, col = "black"))
          }
        }
      }

      ht <- Heatmap(
        mat,
        name = "log2(OR)",
        col = col_fun,
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        cell_fun = cell_fun_final,
        row_names_gp = grid::gpar(
          fontface = ifelse(row_highlight, "bold", "plain"),
          col = ifelse(row_highlight, "#B8860B", "black")
        ),
        column_title = sprintf(
          "TF Family × Cluster Lost Edge Enrichment%s",
          ifelse(PARAMS$show_heatmap_stars, "", " (effect size only)")
        ),
        top_annotation = col_ha,
        na_col = "grey90",
        border = TRUE
      )
      draw(ht)
    }

    save_plot(NULL, "Fig1_Enrichment_Heatmap", OUT_7C,
              width = fig_w, height = fig_h,
              fig_format = PARAMS$fig_format, dpi = PARAMS$fig_dpi,
              scale = PARAMS$fig_scale, device_fn = draw_heatmap)
    "success"
  }, error = function(e) {
    cat(sprintf("    \u26A0 Heatmap failed: %s\n", e$message))
    paste("failed:", e$message)
  })
}

# ========================================================================
# 图 2: Volcano Plot
# ========================================================================

if ("volcano" %in% PARAMS$plots_vec) {
  cat("\n--- 图 2: Volcano Plot (whitelist labeling) ---\n")

  plot_results[["volcano"]] <- tryCatch({

    vol_data <- scored_edges %>%
      filter(!is.na(delta_r) & !is.na(fdr_diff)) %>%
      mutate(
        neg_log10_fdr = -log10(pmax(fdr_diff, 1e-300)),
        neg_log10_fdr = pmin(neg_log10_fdr, 50),
        point_size = case_when(
          Tier == 1 ~ 3, Tier == 2 ~ 2, Tier == 3 ~ 1.5, TRUE ~ 0.8
        ),
        in_diverged = !is.na(Target_Divergence_Ratio) & Target_Divergence_Ratio > 0.5
      )

    cat("    Building whitelist labels...\n")

    highlight_gene_ids <- if (!is.null(highlight_map)) names(highlight_map) else character(0)
    highlight_symbol_lookup <- if (!is.null(highlight_map)) highlight_map else setNames(character(0), character(0))
    all_tf_ids <- if (exists("tf_in_data")) tf_in_data else character(0)

    cat(sprintf("    Highlight genes: %d | TF genes: %d\n",
                length(highlight_gene_ids), length(all_tf_ids)))

    vol_data$label <- NA_character_

    for (idx in seq_len(nrow(vol_data))) {
      tf_id  <- vol_data$TF_Gene[idx]
      tg_id  <- vol_data$Target_Gene[idx]
      tf_sym <- vol_data$TF_Symbol[idx]
      tg_sym <- vol_data$Target_Symbol[idx]

      tf_label <- NA_character_
      if (tf_id %in% all_tf_ids || tf_id %in% highlight_gene_ids) {
        tf_label <- ifelse(!is.na(tf_sym) && tf_sym != "", tf_sym, tf_id)
      }

      tg_label <- NA_character_
      if (tg_id %in% highlight_gene_ids) {
        sym_from_hl <- highlight_symbol_lookup[tg_id]
        tg_label <- ifelse(!is.na(sym_from_hl) && sym_from_hl != "", sym_from_hl,
                           ifelse(!is.na(tg_sym) && tg_sym != "", tg_sym, tg_id))
      } else if (tg_id %in% all_tf_ids) {
        tg_label <- ifelse(!is.na(tg_sym) && tg_sym != "", tg_sym, tg_id)
      }

      if (!is.na(tf_label) && !is.na(tg_label)) {
        vol_data$label[idx] <- paste0(tf_label, "\u2192", tg_label)
      } else if (!is.na(tf_label) && is.na(tg_label)) {
        if (PARAMS$label_focus && PARAMS$focus_gene != "none" &&
            (tf_id == PARAMS$focus_gene || tg_id == PARAMS$focus_gene)) {
          tg_fallback <- ifelse(!is.na(tg_sym) && tg_sym != "", tg_sym, tg_id)
          vol_data$label[idx] <- paste0(tf_label, "\u2192", tg_fallback)
        }
      } else if (is.na(tf_label) && !is.na(tg_label)) {
        if (PARAMS$label_focus && PARAMS$focus_gene != "none" &&
            (tf_id == PARAMS$focus_gene || tg_id == PARAMS$focus_gene)) {
          tf_fallback <- ifelse(!is.na(tf_sym) && tf_sym != "", tf_sym, tf_id)
          vol_data$label[idx] <- paste0(tf_fallback, "\u2192", tg_label)
        }
      }
    }

    labeled_idx <- which(!is.na(vol_data$label))
    if (length(labeled_idx) > PARAMS$label_top_n * 3) {
      top_keep <- vol_data[labeled_idx, ] %>%
        arrange(desc(MECS_score)) %>%
        head(PARAMS$label_top_n * 3) %>%
        rownames()
      drop_idx <- setdiff(labeled_idx, as.integer(top_keep))
      vol_data$label[drop_idx] <- NA_character_
    }

    n_labeled <- sum(!is.na(vol_data$label))
    cat(sprintf("    Final labeled edges: %d / %d\n", n_labeled, nrow(vol_data)))

    p <- ggplot(vol_data, aes(x = delta_r, y = neg_log10_fdr)) +
      geom_point(aes(color = edge_class, shape = in_diverged, size = point_size),
                 alpha = 0.6) +
      scale_color_manual(values = colors$edge_class, name = "Edge Class") +
      scale_shape_manual(values = c(`TRUE` = 17, `FALSE` = 16),
                         labels = c("Conserved Cluster", "Diverged Cluster"),
                         name = "Target Cluster") +
      scale_size_identity() +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
      labs(
        x = expression(Delta * r ~ "(r_WT - r_MUT)"),
        y = expression(-log[10] * "(FDR"["diff"] * ")"),
        title = "Differential Co-expression Volcano",
        subtitle = "Labels: TFs & highlight genes only (whitelist)"
      ) +
      theme_bw(base_size = 12) +
      theme(legend.position = "right")

    if (requireNamespace("ggrepel", quietly = TRUE) && any(!is.na(vol_data$label))) {
      p <- p + ggrepel::geom_text_repel(
        aes(label = label),
        size = 2.5, max.overlaps = 40, segment.alpha = 0.4,
        fontface = "bold", box.padding = 0.4, point.padding = 0.3,
        min.segment.length = 0.2
      )
    }

    save_plot(p, "Fig2_Volcano", OUT_7C,
              width = 10, height = 7,
              fig_format = PARAMS$fig_format, dpi = PARAMS$fig_dpi,
              scale = PARAMS$fig_scale)
    "success"
  }, error = function(e) {
    cat(sprintf("    \u26A0 Volcano failed: %s\n", e$message))
    paste("failed:", e$message)
  })
}

# ========================================================================
# 图 3: Cluster-wise Edge Class Distribution
# ========================================================================

if ("barplot" %in% PARAMS$plots_vec) {
  cat("\n--- 图 3: Edge Class Barplot ---\n")

  plot_results[["barplot"]] <- tryCatch({

    bar_data <- scored_edges %>%
      filter(!is.na(Target_WT_Cluster)) %>%
      count(Target_WT_Cluster, edge_class, name = "n_edges") %>%
      mutate(Cluster = factor(paste0("C", Target_WT_Cluster),
                               levels = paste0("C", sort(unique(Target_WT_Cluster)))))

    fig_w <- max(8, n_wt_clusters * 1.2 + 2)

    p <- ggplot(bar_data, aes(x = Cluster, y = n_edges, fill = edge_class)) +
      geom_bar(stat = "identity", position = "stack") +
      scale_fill_manual(values = colors$edge_class, name = "Edge Class") +
      labs(x = "WT Cluster", y = "Number of Edges",
           title = "Edge Class Distribution by Cluster") +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    save_plot(p, "Fig3_EdgeClass_Barplot", OUT_7C,
              width = fig_w, height = 5,
              fig_format = PARAMS$fig_format, dpi = PARAMS$fig_dpi,
              scale = PARAMS$fig_scale)
    "success"
  }, error = function(e) {
    cat(sprintf("    \u26A0 Barplot failed: %s\n", e$message))
    paste("failed:", e$message)
  })
}

# ========================================================================
# 图 4: Anchor-Based Directed Network (v4.4: 参数化 Degree)
# ========================================================================

if ("network" %in% PARAMS$plots_vec) {
  cat("\n--- 图 4: Anchor-Based Directed Network (v4.4) ---\n")

  plot_results[["network"]] <- tryCatch({

    if (is.null(network_data)) {
      stop("Network data not available (Part 1 failed)")
    }

    suppressPackageStartupMessages({
      library(igraph)
      library(ggraph)
      library(tidygraph)
    })

    net_edges  <- network_data$net_edges
    node_meta  <- network_data$node_meta
    anchor_genes <- network_data$anchor_genes

    # ===================================================================
    # Step 4g: 构建边数据
    # ===================================================================

    edge_data <- net_edges %>%
      select(
        from = TF_Gene, to = Target_Gene,
        edge_class, delta_r,
        MECS_score = any_of("MECS_score"),
        Tier       = any_of("Tier")
      ) %>%
      mutate(
        edge_width = pmax(abs(delta_r), 0.01),
        edge_alpha = case_when(
          edge_class %in% c("Lost", "Gained", "Rewired") ~ 0.8,
          edge_class == "Conserved" ~ 0.4,
          TRUE ~ 0.5
        )
      )

    pruned_class_counts <- table(edge_data$edge_class)
    cat(sprintf("    Edge summary: %d total\n", nrow(edge_data)))
    cat(sprintf("       Lost: %d | Gained: %d | Rewired: %d | Conserved: %d\n",
                pruned_class_counts["Lost"] %||% 0,
                pruned_class_counts["Gained"] %||% 0,
                pruned_class_counts["Rewired"] %||% 0,
                pruned_class_counts["Conserved"] %||% 0))

    # ===================================================================
    # Step 4h: 构建 tidygraph 对象
    # ===================================================================

    g <- tbl_graph(nodes = node_meta, edges = edge_data, directed = TRUE) %>%
      mutate(degree = centrality_degree(mode = "all"))

    g <- g %>% filter(degree > 0)

    final_n <- igraph::vcount(as.igraph(g))
    final_e <- igraph::ecount(as.igraph(g))
    cat(sprintf("    Final graph: %d nodes, %d edges\n", final_n, final_e))

    # ===================================================================
    # Step 4i: 布局
    # ===================================================================

    layout_algo <- switch(PARAMS$network_layout,
      "kk" = "kk", "fr" = "fr", "stress" = "stress", "kk"
    )
    set.seed(PARAMS$network_seed)
    cat(sprintf("    Layout: %s (seed=%d)\n", layout_algo, PARAMS$network_seed))

    # ===================================================================
    # Step 4j: ggraph 绘图
    # ===================================================================

    # ★ v4.2: 在 subtitle/caption 中注明 ego 过滤信息
    ego_note <- if (isTRUE(network_data$ego_filter_applied)) {
      sprintf("Ego filter: %s (removed %d nodes)",
              ifelse(network_data$ego_steps_used > 0,
                     paste0(network_data$ego_steps_used, "-step"),
                     "connected component"),
              network_data$ego_n_removed)
    } else {
      "Ego filter: not applied"
    }

    p <- ggraph(g, layout = layout_algo) +
      geom_edge_link(
        aes(color = edge_class, width = edge_width,
            linetype = edge_class, alpha = edge_alpha),
        arrow = arrow(length = unit(PARAMS$network_arrow_mm, "mm"), type = "closed"),
        start_cap = circle(3, "mm"), end_cap = circle(5, "mm")
      ) +
      scale_edge_color_manual(values = colors$edge_class, name = "Edge Class") +
      scale_edge_linetype_manual(
        values = c("Lost" = "dashed", "Gained" = "solid",
                   "Rewired" = "longdash", "Conserved" = "dotted"),
        name = "Edge Class"
      ) +
      scale_edge_width(
        range = c(PARAMS$network_edge_width_min, PARAMS$network_edge_width_max),
        guide = "none"
      ) +
      scale_edge_alpha_identity() +
      geom_node_point(
        aes(fill = NodeType, shape = NodeType, size = degree),
        color = "white", stroke = 1.2
      ) +
      scale_fill_manual(values = colors$node_type, name = "Node Type") +
      scale_shape_manual(
        values = c("Focus TF" = 23, "Other TF" = 23,
                   "Highlight" = 21, "Target" = 21),
        name = "Node Type"
      ) +
      scale_size_continuous(
        range = c(PARAMS$network_node_size_min, PARAMS$network_node_size_max),
        name = "Degree",
        guide = guide_legend(
          override.aes = list(shape = 21, fill = "grey50", color = "white", stroke = 0.5)
        )
      ) +
      geom_node_text(
        aes(label = Label), repel = TRUE, size = 3.5, fontface = "bold",
        bg.color = "white", bg.r = 0.15, max.overlaps = 40, show.legend = FALSE
      ) +
      theme_graph(base_family = "sans", background = "white") +
      theme(
        plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
        plot.caption  = element_text(size = 8, hjust = 1, color = "grey60"),
        legend.position = "right", legend.box = "vertical",
        legend.spacing.y = unit(0.3, "cm"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 10, face = "bold")
      ) +
      labs(
        title = sprintf("Anchor-Based GRN: %s + Highlight Genes",
          ifelse(length(focus_family_vec) > 0,
                 paste(focus_family_vec, collapse = " + "), "All TFs")),
        subtitle = sprintf("Nodes: %d | Edges: %d (L=%d, G=%d, R=%d, C=%d)",
          final_n, final_e,
          pruned_class_counts["Lost"] %||% 0,
          pruned_class_counts["Gained"] %||% 0,
          pruned_class_counts["Rewired"] %||% 0,
          pruned_class_counts["Conserved"] %||% 0),
        caption = sprintf("Pruning: Conserved excl. from degree | deg1_all=%s, deg1_focus=%s | %s | Layout: %s | Seed: %d",
          PARAMS$network_keep_deg1_all, PARAMS$network_keep_deg1_focus,
          ego_note, layout_algo, PARAMS$network_seed)
      )

    fig_w <- PARAMS$network_fig_w
    fig_h <- PARAMS$network_fig_h

    save_plot(p, "Fig4_AnchorNetwork", OUT_7C,
              width = fig_w, height = fig_h,
              fig_format = PARAMS$fig_format, dpi = PARAMS$fig_dpi,
              scale = PARAMS$fig_scale)

    # --- Cytoscape CSV 导出 ---
    cat("    Exporting Cytoscape-compatible CSVs...\n")

    cyto_nodes <- as_tibble(g, what = "nodes") %>%
      mutate(
        Display_Label = case_when(!is.na(Symbol) & Symbol != "" ~ Symbol,
                                  TRUE ~ shorten_id(id)),
        Cluster_Label = ifelse(!is.na(WT_Cluster), paste0("C", WT_Cluster), "NA")
      ) %>%
      select(Node_ID = id, Display_Label, NodeType, Family,
             Is_Focus_TF, Is_TF, Is_Highlight, Is_Focus_Gene, Is_Anchor,
             Cluster_Label, Degree = degree)

    write.csv(cyto_nodes,
              file.path(OUT_7C, "Pipeline7c_Fig4_Cytoscape_Nodes.csv"),
              row.names = FALSE)

    cyto_edges <- edge_data %>%
      select(Source = from, Target = to, Edge_Class = edge_class,
             Delta_r = delta_r, Weight = edge_width)
    if ("MECS_score" %in% colnames(edge_data)) cyto_edges$MECS_Score <- edge_data$MECS_score
    if ("Tier" %in% colnames(edge_data)) cyto_edges$Tier <- edge_data$Tier

    write.csv(cyto_edges,
              file.path(OUT_7C, "Pipeline7c_Fig4_Cytoscape_Edges.csv"),
              row.names = FALSE)

    cat(sprintf("    Cytoscape export: %d nodes, %d edges\n",
                nrow(cyto_nodes), nrow(cyto_edges)))

    # --- Network Summary JSON (★ v4.2: 含 ego 信息) ---
    network_summary <- list(
      n_anchor_genes = length(anchor_genes),
      n_final_nodes = final_n, n_final_edges = final_e,
      edge_class_counts = as.list(pruned_class_counts),
      node_type_counts = as.list(table(node_meta$NodeType)),
      pruning_note = sprintf("v4.4: keep_deg1_all=%s, keep_deg1_focus=%s",
                            PARAMS$network_keep_deg1_all, PARAMS$network_keep_deg1_focus),
      ego_filter = list(
        applied       = network_data$ego_filter_applied,
        ego_steps     = network_data$ego_steps_used,
        nodes_removed = network_data$ego_n_removed,
        focus_gene    = PARAMS$focus_gene
      ),
      layout = layout_algo, fig_size = c(fig_w, fig_h)
    )

    writeLines(
      jsonlite::toJSON(network_summary, pretty = TRUE, auto_unbox = TRUE),
      file.path(OUT_7C, "Pipeline7c_Fig4_Network_Summary.json")
    )

    "success"
  }, error = function(e) {
    cat(sprintf("    \u26A0 Network failed: %s\n", e$message))
    paste("failed:", e$message)
  })
}

# ========================================================================
# 图 5: Circos
# ========================================================================

if ("circos" %in% PARAMS$plots_vec) {
  cat("\n--- 图 5: Circos (Edge-Centric Global Ranking) ---\n")
  cat(sprintf("    Strategy: Top %d MECS edges\n", PARAMS$circos_top_edges))

  plot_results[["circos"]] <- tryCatch({
    suppressPackageStartupMessages(library(circlize))

    circ_data <- scored_edges %>%
      filter(Tier %in% c(1, 2, 3),
             edge_class %in% c("Lost", "Gained", "Rewired"),
             !is.na(Target_WT_Cluster)) %>%
      mutate(TF_Family = tf_family_map[TF_Gene]) %>%
      filter(!is.na(TF_Family))

    if (nrow(circ_data) == 0) stop("No data for circos")

    top_n <- min(PARAMS$circos_top_edges, nrow(circ_data))
    circ_data <- circ_data %>% arrange(desc(MECS_score)) %>% head(top_n)

    cat(sprintf("    Selected top %d edges\n", nrow(circ_data)))

    families_in <- sort(unique(circ_data$TF_Family))

    if (nrow(circ_data) > PARAMS$circos_max_edges) {
      circ_data <- head(circ_data, PARAMS$circos_max_edges)
    }

    circ_summary <- circ_data %>%
      count(TF_Family, Target_WT_Cluster, name = "n_edges")

    adj_mat <- circ_summary %>%
      mutate(Cluster = paste0("C", Target_WT_Cluster)) %>%
      select(TF_Family, Cluster, n_edges) %>%
      pivot_wider(names_from = Cluster, values_from = n_edges, values_fill = 0) %>%
      column_to_rownames("TF_Family") %>%
      as.matrix()

    grid_colors <- c()
    for (fam in rownames(adj_mat)) {
      if (fam %in% focus_family_vec) {
        grid_colors[fam] <- colors$focus
      } else {
        fam_idx <- match(fam, families_in)
        n_fams <- length(families_in)
        fam_pal <- if (n_fams <= 8) {
          RColorBrewer::brewer.pal(max(3, n_fams), "Set2")
        } else {
          colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_fams)
        }
        grid_colors[fam] <- fam_pal[fam_idx]
      }
    }
    for (cl in colnames(adj_mat)) {
      grid_colors[cl] <- if (cl %in% names(colors$cluster)) colors$cluster[cl] else "#CCCCCC"
    }

    draw_circos <- function() {
      circos.clear()
      circos.par(start.degree = 90, gap.degree = 2)
      chordDiagram(adj_mat, grid.col = grid_colors,
                   annotationTrack = "grid", preAllocateTracks = 1,
                   transparency = 0.3)
      circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
        sector_name <- CELL_META$sector.index
        xlim <- CELL_META$xlim; ylim <- CELL_META$ylim
        circos.text(mean(xlim), ylim[1], sector_name,
                    facing = "clockwise", niceFacing = TRUE,
                    adj = c(0, 0.5), cex = 0.6)
      }, bg.border = NA)
      title(sprintf("TF Family \u2192 Cluster (Top %d MECS Edges)", nrow(circ_data)))
      circos.clear()
    }

    save_plot(NULL, "Fig5_Circos", OUT_7C,
              width = 10, height = 10,
              fig_format = PARAMS$fig_format, dpi = PARAMS$fig_dpi,
              scale = PARAMS$fig_scale, device_fn = draw_circos)
    "success"
  }, error = function(e) {
    cat(sprintf("    \u26A0 Circos failed: %s\n", e$message))
    paste("failed:", e$message)
  })
}

# ========================================================================
# 图 6: Delta Out-Degree Bar Plot
# ========================================================================

if ("outdegree" %in% PARAMS$plots_vec) {
  cat("\n--- 图 6: Delta Out-Degree ---\n")

  plot_results[["outdegree"]] <- tryCatch({

    if (!exists("tf_summary") || nrow(tf_summary) == 0) stop("No TF summary data")

    if (!"TF_Symbol" %in% colnames(tf_summary)) tf_summary$TF_Symbol <- NA_character_
    if (!"TF_Family" %in% colnames(tf_summary)) tf_summary$TF_Family <- NA_character_

    od_data <- tf_summary %>%
      arrange(desc(abs(delta_out_degree))) %>%
      head(30) %>%
      mutate(
        Name_Part = ifelse(!is.na(TF_Symbol) & TF_Symbol != "",
                           TF_Symbol, TF_Gene),
        Family_Part = ifelse(!is.na(TF_Family) & TF_Family != "",
                             TF_Family, "unknown"),
        TF_label = paste0(Name_Part, " (", Family_Part, ")"),
        is_focus = TF_Family %in% focus_family_vec
      )

    od_data$TF_label <- factor(od_data$TF_label, levels = rev(od_data$TF_label))

    p <- ggplot(od_data, aes(x = TF_label, y = delta_out_degree, fill = is_focus)) +
      geom_bar(stat = "identity") +
      scale_fill_manual(values = c(`TRUE` = colors$focus, `FALSE` = "#3498DB"),
                        labels = c("Other TF", "Focus Family"),
                        name = "TF Category") +
      coord_flip() +
      labs(x = NULL, y = expression(Delta ~ "Out-Degree (WT - MUT)"),
           title = "Top TFs by Regulatory Rewiring Magnitude") +
      theme_bw(base_size = 11) +
      theme(axis.text.y = element_text(size = 9)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50")

    fig_h <- max(6, nrow(od_data) * 0.28 + 1.5)

    save_plot(p, "Fig6_DeltaOutDegree", OUT_7C,
              width = 12, height = fig_h,
              fig_format = PARAMS$fig_format, dpi = PARAMS$fig_dpi,
              scale = PARAMS$fig_scale)
    "success"
  }, error = function(e) {
    cat(sprintf("    \u26A0 OutDegree failed: %s\n", e$message))
    paste("failed:", e$message)
  })
}

# ========================================================================
# 图 7: Differential Expression Evidence Heatmap  ★ v4.2 ★
# ========================================================================

if ("evidence_heatmap" %in% PARAMS$plots_vec) {
  cat("\n--- 图 7: DE Evidence Heatmap (v4.2, ego-filtered) ---\n")

  plot_results[["evidence_heatmap"]] <- tryCatch({

    suppressPackageStartupMessages({
      library(ComplexHeatmap)
      library(circlize)
      library(grid)
      library(RColorBrewer)
    })

    # ===================================================================
    # Step 7a: 前置检查
    # ===================================================================

    if (is.null(network_data)) {
      stop("Network data not available (Part 1 failed or not executed)")
    }

    if (is.null(PARAMS$deseq2_rdata) || !file.exists(PARAMS$deseq2_rdata)) {
      stop(sprintf("DESeq2 RData not found: %s (use --deseq2_rdata)",
                   ifelse(is.null(PARAMS$deseq2_rdata), "NULL", PARAMS$deseq2_rdata)))
    }

    # ===================================================================
    # Step 7b: 加载 DESeq2 结果
    # ===================================================================

    cat("    Loading DESeq2 results...\n")
    deseq2_env <- new.env()
    load(PARAMS$deseq2_rdata, envir = deseq2_env)
    all_comps <- deseq2_env$all_comparisons
    deseq2_meta <- deseq2_env$meta

    cat(sprintf("    DESeq2 comparisons loaded: %d total\n", length(all_comps)))
    cat(sprintf("    Comparison names (first 6): %s\n",
                paste(head(names(all_comps)), collapse = ", ")))

    # ===================================================================
    # Step 7c: 配置参数
    # ===================================================================

    cond_wt  <- PARAMS$condition_wt
    cond_mut <- PARAMS$condition_mut
    ctrl_time <- PARAMS$control_time

    all_times <- sort(unique(as.character(deseq2_meta$Time)))
    non_ctrl_times <- setdiff(all_times, ctrl_time)

    cat(sprintf("    Conditions: %s, %s | Control time: %s\n",
                cond_wt, cond_mut, ctrl_time))
    cat(sprintf("    Non-control timepoints: %s\n",
                paste(non_ctrl_times, collapse = ", ")))

    time_display_map <- setNames(gsub("^0", "", non_ctrl_times), non_ctrl_times)

    # ===================================================================
    # Step 7d: 基因池构建
    # ===================================================================

    nd <- network_data

    target_genes_pool <- nd$node_meta %>%
      filter(NodeType %in% c("Target", "Highlight")) %>%
      pull(id)

    tf_genes_pool <- nd$node_meta %>%
      filter(NodeType %in% c("Focus TF", "Other TF")) %>%
      pull(id)

    first_comp_genes <- rownames(as.data.frame(all_comps[[1]]))
    target_genes_pool <- intersect(target_genes_pool, first_comp_genes)
    tf_genes_pool     <- intersect(tf_genes_pool, first_comp_genes)

    cat(sprintf("    Gene pools (DESeq2-matched): %d targets, %d TFs\n",
                length(target_genes_pool), length(tf_genes_pool)))

    if (length(target_genes_pool) < 2) stop("Too few target genes with DESeq2 data")
    if (length(tf_genes_pool) < 1)     stop("No TF genes with DESeq2 data")

    # ===================================================================
    # Step 7e: 构建 log2FC / 星号 / 方向性方框 矩阵
    # ===================================================================

    build_de_matrices <- function(gene_ids) {

      col_names <- character(0)
      for (tp in non_ctrl_times) {
        col_names <- c(col_names,
                       paste0(cond_wt, "_", tp),
                       paste0(cond_mut, "_", tp))
      }

      n_genes <- length(gene_ids)
      n_cols  <- length(col_names)

      lfc_mat  <- matrix(NA_real_,     nrow = n_genes, ncol = n_cols,
                         dimnames = list(gene_ids, col_names))
      star_mat <- matrix("",           nrow = n_genes, ncol = n_cols,
                         dimnames = list(gene_ids, col_names))
      box_mat  <- matrix(FALSE,        nrow = n_genes, ncol = n_cols,
                         dimnames = list(gene_ids, col_names))

      for (tp in non_ctrl_times) {

        col_wt  <- paste0(cond_wt, "_", tp)
        col_mut <- paste0(cond_mut, "_", tp)

        # --- Logic 1: WT timepoint vs 0h ---
        key_wt <- paste0("Time_", cond_wt, "_", tp, "_vs_", ctrl_time)
        if (key_wt %in% names(all_comps)) {
          res_wt <- as.data.frame(all_comps[[key_wt]])
          matched <- intersect(gene_ids, rownames(res_wt))
          lfc_mat[matched, col_wt] <- res_wt[matched, "log2FoldChange"]
          padj_vec <- res_wt[matched, "padj"]
          star_mat[matched, col_wt] <- ifelse(
            is.na(padj_vec), "",
            ifelse(padj_vec < 0.01, "**",
                   ifelse(padj_vec < 0.05, "*", ""))
          )
        } else {
          cat(sprintf("    \u26A0 Logic 1 comparison not found: %s\n", key_wt))
        }

        # --- Logic 1: MUT timepoint vs 0h ---
        key_mut <- paste0("Time_", cond_mut, "_", tp, "_vs_", ctrl_time)
        if (key_mut %in% names(all_comps)) {
          res_mut <- as.data.frame(all_comps[[key_mut]])
          matched <- intersect(gene_ids, rownames(res_mut))
          lfc_mat[matched, col_mut] <- res_mut[matched, "log2FoldChange"]
          padj_vec <- res_mut[matched, "padj"]
          star_mat[matched, col_mut] <- ifelse(
            is.na(padj_vec), "",
            ifelse(padj_vec < 0.01, "**",
                   ifelse(padj_vec < 0.05, "*", ""))
          )
        } else {
          cat(sprintf("    \u26A0 Logic 1 comparison not found: %s\n", key_mut))
        }

        # --- Logic 2: Genotype comparison at timepoint (方向性方框) ---
        key_geno <- paste0("Genotype_", cond_mut, "_vs_", cond_wt, "_at_", tp)
        if (key_geno %in% names(all_comps)) {
          res_geno <- as.data.frame(all_comps[[key_geno]])
          matched <- intersect(gene_ids, rownames(res_geno))
          padj_geno <- res_geno[matched, "padj"]
          lfc_geno  <- res_geno[matched, "log2FoldChange"]

          for (k in seq_along(matched)) {
            g <- matched[k]
            p_val   <- padj_geno[k]
            lfc_val <- lfc_geno[k]

            if (!is.na(p_val) && p_val < 0.05 && !is.na(lfc_val)) {
              if (lfc_val >= 1) {
                box_mat[g, col_mut] <- TRUE
              } else if (lfc_val <= -1) {
                box_mat[g, col_wt]  <- TRUE
              }
            }
          }
        } else {
          cat(sprintf("    \u26A0 Logic 2 comparison not found: %s\n", key_geno))
        }
      }

      cat(sprintf("      Matrix built: %d genes \u00D7 %d columns\n", n_genes, n_cols))
      cat(sprintf("      Stars: %d cells | Directional boxes: %d cells\n",
                  sum(star_mat != ""), sum(box_mat)))

      list(lfc = lfc_mat, stars = star_mat, boxes = box_mat, col_names = col_names)
    }

    cat("    Building Target matrices...\n")
    target_mats <- build_de_matrices(target_genes_pool)

    cat("    Building TF matrices...\n")
    tf_mats <- build_de_matrices(tf_genes_pool)

    # ===================================================================
    # Step 7f: 配色方案
    # ===================================================================

    brg_pal <- rev(brewer.pal(11, "BrBG"))
    body_col_fun <- colorRamp2(seq(-5, 5, length.out = 11), brg_pal)

    # ===================================================================
    # Step 7g: 列元数据辅助函数
    # ===================================================================

    make_col_split <- function(col_names) {
      tps <- sub(paste0("^(", cond_wt, "|", cond_mut, ")_"), "", col_names)
      factor(time_display_map[tps], levels = time_display_map[non_ctrl_times])
    }

    make_col_labels <- function(col_names) {
      sub("_.*$", "", col_names)
    }

    make_col_fontface <- function(col_names) {
      ifelse(grepl(paste0("^", cond_mut), col_names), "italic", "plain")
    }

    # ===================================================================
    # Step 7h: 通用 cell_fun 工厂
    # ===================================================================

    make_cell_fun <- function(star_mat, box_mat) {
      force(star_mat)
      force(box_mat)
      function(j, i, x, y, width, height, fill) {
        if (!is.na(box_mat[i, j]) && box_mat[i, j]) {
          grid.rect(x, y, width * 0.85, height * 0.85,
                    gp = gpar(col = "#1A1A1A", lwd = 2.2, fill = NA))
        }
        s <- star_mat[i, j]
        if (!is.na(s) && nchar(s) > 0) {
          grid.text(s, x, y,
                    gp = gpar(fontsize = 6, col = "black", fontface = "bold"))
        }
      }
    }

    # ===================================================================
    # Step 7i: 绘制 Target 基因热图 (Fig7a)
    # ===================================================================

    cat("    Drawing Target DE Heatmap (Fig7a)...\n")

    target_meta <- nd$node_meta %>%
      filter(id %in% target_genes_pool) %>%
      select(id, Symbol, WT_Cluster, Is_Highlight)

    target_edge_dom <- nd$net_edges %>%
      filter(Target_Gene %in% target_genes_pool, edge_class != "Conserved") %>%
      count(Target_Gene, edge_class) %>%
      group_by(Target_Gene) %>%
      slice_max(n, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(Target_Gene, Edge_Dominant = edge_class)

    targets_only_conserved <- setdiff(target_genes_pool, target_edge_dom$Target_Gene)
    if (length(targets_only_conserved) > 0) {
      target_edge_dom <- bind_rows(
        target_edge_dom,
        data.frame(Target_Gene = targets_only_conserved,
                   Edge_Dominant = "Conserved", stringsAsFactors = FALSE)
      )
    }

    target_meta <- target_meta %>%
      left_join(target_edge_dom, by = c("id" = "Target_Gene")) %>%
      mutate(
        Edge_Dominant = ifelse(is.na(Edge_Dominant), "Conserved", Edge_Dominant),
        WT_Cluster_label = ifelse(!is.na(WT_Cluster), paste0("C", WT_Cluster), "NA"),
        Gene_Label = case_when(
          !is.na(Symbol) & Symbol != "" ~ Symbol,
          TRUE ~ shorten_id(id)
        )
      )

    edge_order <- c("Lost", "Gained", "Rewired", "Conserved")
    target_meta <- target_meta %>%
      mutate(
        cs = ifelse(is.na(WT_Cluster), 999, WT_Cluster),
        es = match(Edge_Dominant, edge_order, nomatch = 99)
      ) %>%
      arrange(cs, es, Gene_Label)

    target_lfc   <- target_mats$lfc[target_meta$id, , drop = FALSE]
    target_stars <- target_mats$stars[target_meta$id, , drop = FALSE]
    target_boxes <- target_mats$boxes[target_meta$id, , drop = FALSE]

    edge_cls_in_data <- unique(target_meta$Edge_Dominant)
    edge_col_map <- colors$edge_class[intersect(names(colors$edge_class), edge_cls_in_data)]

    target_left_anno <- rowAnnotation(
      Edge = target_meta$Edge_Dominant,
      col  = list(Edge = edge_col_map),
      annotation_name_gp   = gpar(fontsize = 7),
      simple_anno_size     = unit(3, "mm"),
      gap                  = unit(0.5, "mm"),
      annotation_legend_param = list(
        Edge = list(title_gp  = gpar(fontsize = 8, fontface = "bold"),
                    labels_gp = gpar(fontsize = 7))
      )
    )

    cluster_labels_in_data <- unique(target_meta$WT_Cluster_label)
    target_row_split <- factor(
      target_meta$WT_Cluster_label,
      levels = c(paste0("C", sort(unique(
        target_meta$WT_Cluster[!is.na(target_meta$WT_Cluster)]
      ))), if ("NA" %in% cluster_labels_in_data) "NA")
    )

    col_names_all    <- colnames(target_lfc)
    col_split_all    <- make_col_split(col_names_all)
    col_labels_all   <- make_col_labels(col_names_all)
    col_fontface_all <- make_col_fontface(col_names_all)

    n_tgt   <- nrow(target_lfc)
    fig_w_t <- max(8, length(col_names_all) * 0.9 + 6)
    fig_h_t <- max(6, n_tgt * 0.25 + 3)
    fig_h_t <- min(fig_h_t, 35)

    cat(sprintf("    Target heatmap: %d genes \u00D7 %d columns | %.1f \u00D7 %.1f in\n",
                n_tgt, length(col_names_all), fig_w_t, fig_h_t))

    draw_target_heatmap <- function() {

      ht <- Heatmap(
        target_lfc,
        name = "log2(FC)",
        col  = body_col_fun,
        na_col = "grey90",

        cluster_rows    = FALSE,
        cluster_columns = FALSE,

        column_split     = col_split_all,
        column_gap       = unit(3, "mm"),
        column_labels    = col_labels_all,
        column_names_gp  = gpar(fontface = col_fontface_all, fontsize = 9),
        column_names_rot = 45,
        column_title_gp  = gpar(fontsize = 11, fontface = "bold"),

        row_split      = target_row_split,
        row_gap        = unit(1.5, "mm"),
        row_title_gp   = gpar(fontsize = 9, fontface = "bold"),
        row_title_rot  = 0,

        row_labels     = target_meta$Gene_Label,
        row_names_side = "right",
        row_names_gp   = gpar(
          fontsize = ifelse(n_tgt > 80, 4, ifelse(n_tgt > 40, 5.5, 7)),
          fontface = ifelse(target_meta$Is_Highlight, "bold", "plain")
        ),

        left_annotation = target_left_anno,

        cell_fun = make_cell_fun(target_stars, target_boxes),

        border  = TRUE,
        rect_gp = gpar(col = "white", lwd = 0.5),
        use_raster     = (n_tgt * length(col_names_all) > 2000),
        raster_quality = 4,

        heatmap_legend_param = list(
          title     = "log2(FC)",
          at        = c(-5, -2.5, 0, 2.5, 5),
          labels    = c("\u2264-5.0", "-2.5", "0.0", "2.5", "\u22655.0"),
          title_gp  = gpar(fontsize = 9, fontface = "bold"),
          labels_gp = gpar(fontsize = 8),
          legend_height = unit(35, "mm"),
          direction     = "vertical",
          border        = "grey60"
        )
      )

      draw(ht,
           heatmap_legend_side     = "right",
           annotation_legend_side  = "right",
           padding = unit(c(10, 5, 18, 5), "mm"))

      grid.text(
        expression(
          paste("* ", italic(p[adj]), " < 0.05   ** ", italic(p[adj]),
                " < 0.01  (vs 0h)       ",
                "\u25A1 ", italic(p[adj]),
                " < 0.05  (WT vs MUT, marked on higher side)")
        ),
        x = 0.5, y = unit(4, "mm"),
        just = c("center", "bottom"),
        gp = gpar(fontsize = 7.5, col = "grey30")
      )
    }

    save_plot(NULL, "Fig7a_Target_DE_Heatmap", OUT_7C,
              width = fig_w_t, height = fig_h_t,
              fig_format = PARAMS$fig_format, dpi = PARAMS$fig_dpi,
              scale = PARAMS$fig_scale, device_fn = draw_target_heatmap)

    cat("    \u2713 Fig7a saved.\n")

    # ===================================================================
    # Step 7j: 绘制 TF 基因热图 (Fig7b)
    # ===================================================================

    cat("    Drawing TF DE Heatmap (Fig7b)...\n")

    tf_meta <- nd$node_meta %>%
      filter(id %in% tf_genes_pool) %>%
      select(id, Symbol, Family, Is_Focus_TF)

    tf_out_deg <- nd$net_edges %>%
      filter(TF_Gene %in% tf_genes_pool) %>%
      count(TF_Gene, name = "out_degree")

    tf_meta <- tf_meta %>%
      left_join(tf_out_deg, by = c("id" = "TF_Gene")) %>%
      mutate(
        out_degree = ifelse(is.na(out_degree), 0L, out_degree),
        TF_Label = case_when(
          !is.na(Symbol) & Symbol != "" & !is.na(Family) ~
            paste0(Symbol, " (", Family, ")"),
          !is.na(Family) ~ paste0(shorten_id(id), " (", Family, ")"),
          !is.na(Symbol) & Symbol != "" ~ Symbol,
          TRUE ~ shorten_id(id)
        )
      )

    tf_meta <- tf_meta %>%
      arrange(desc(out_degree))

    tf_lfc   <- tf_mats$lfc[tf_meta$id, , drop = FALSE]
    tf_stars <- tf_mats$stars[tf_meta$id, , drop = FALSE]
    tf_boxes <- tf_mats$boxes[tf_meta$id, , drop = FALSE]

    tf_right_anno <- rowAnnotation(
      Out_Degree = anno_barplot(
        tf_meta$out_degree,
        gp    = gpar(fill = ifelse(tf_meta$Is_Focus_TF, colors$focus, "#457B9D")),
        width = unit(15, "mm"),
        border = FALSE
      ),
      annotation_name_gp = gpar(fontsize = 7),
      gap = unit(0.5, "mm")
    )

    n_tf     <- nrow(tf_lfc)
    fig_w_tf <- max(8, length(col_names_all) * 0.9 + 6)
    fig_h_tf <- max(6, n_tf * 0.28 + 3)
    fig_h_tf <- min(fig_h_tf, 35)

    cat(sprintf("    TF heatmap: %d TFs \u00D7 %d columns | %.1f \u00D7 %.1f in\n",
                n_tf, length(col_names_all), fig_w_tf, fig_h_tf))

    draw_tf_heatmap <- function() {

      ht <- Heatmap(
        tf_lfc,
        name = "log2(FC)",
        col  = body_col_fun,
        na_col = "grey90",

        cluster_rows    = FALSE,
        cluster_columns = FALSE,

        column_split     = col_split_all,
        column_gap       = unit(3, "mm"),
        column_labels    = col_labels_all,
        column_names_gp  = gpar(fontface = col_fontface_all, fontsize = 9),
        column_names_rot = 45,
        column_title_gp  = gpar(fontsize = 11, fontface = "bold"),

        row_gap        = unit(1.5, "mm"),
        row_title_gp   = gpar(fontsize = 9, fontface = "bold"),
        row_title_rot  = 0,

        row_labels     = tf_meta$TF_Label,
        row_names_side = "left",
        row_names_gp   = gpar(
          fontsize = ifelse(n_tf > 40, 5.5, ifelse(n_tf > 20, 6.5, 8)),
          fontface = ifelse(tf_meta$Is_Focus_TF, "bold", "plain"),
          col      = ifelse(tf_meta$Is_Focus_TF, "#B8860B", "black")
        ),

        right_annotation = tf_right_anno,

        cell_fun = make_cell_fun(tf_stars, tf_boxes),

        border  = TRUE,
        rect_gp = gpar(col = "white", lwd = 0.5),
        use_raster     = (n_tf * length(col_names_all) > 2000),
        raster_quality = 4,

        heatmap_legend_param = list(
          title     = "log2(FC)",
          at        = c(-5, -2.5, 0, 2.5, 5),
          labels    = c("\u2264-5.0", "-2.5", "0.0", "2.5", "\u22655.0"),
          title_gp  = gpar(fontsize = 9, fontface = "bold"),
          labels_gp = gpar(fontsize = 8),
          legend_height = unit(35, "mm"),
          direction     = "vertical",
          border        = "grey60"
        )
      )

      draw(ht,
           heatmap_legend_side     = "right",
           annotation_legend_side  = "right",
           padding = unit(c(10, 5, 18, 5), "mm"))

      grid.text(
        expression(
          paste("* ", italic(p[adj]), " < 0.05   ** ", italic(p[adj]),
                " < 0.01  (vs 0h)       ",
                "\u25A1 ", italic(p[adj]),
                " < 0.05  (WT vs MUT, marked on higher side)")
        ),
        x = 0.5, y = unit(4, "mm"),
        just = c("center", "bottom"),
        gp = gpar(fontsize = 7.5, col = "grey30")
      )
    }

    save_plot(NULL, "Fig7b_TF_DE_Heatmap", OUT_7C,
              width = fig_w_tf, height = fig_h_tf,
              fig_format = PARAMS$fig_format, dpi = PARAMS$fig_dpi,
              scale = PARAMS$fig_scale, device_fn = draw_tf_heatmap)

    cat("    \u2713 Fig7b saved.\n")

    # ===================================================================
    # Step 7k: 导出辅助数据
    # ===================================================================

    cat("    Exporting supplementary data...\n")

    target_export <- cbind(
      data.frame(
        Gene_ID         = target_meta$id,
        Symbol          = target_meta$Symbol,
        WT_Cluster      = target_meta$WT_Cluster_label,
        Edge_Dominant   = target_meta$Edge_Dominant,
        Is_Highlight    = target_meta$Is_Highlight,
        stringsAsFactors = FALSE
      ),
      as.data.frame(target_lfc)
    )
    write.csv(target_export,
              file.path(OUT_7C, "Pipeline7c_Fig7a_Target_LFC_Matrix.csv"),
              row.names = FALSE)

    tf_export <- cbind(
      data.frame(
        Gene_ID      = tf_meta$id,
        Symbol       = tf_meta$Symbol,
        Family       = tf_meta$Family,
        Is_Focus_TF  = tf_meta$Is_Focus_TF,
        Out_Degree   = tf_meta$out_degree,
        stringsAsFactors = FALSE
      ),
      as.data.frame(tf_lfc)
    )
    write.csv(tf_export,
              file.path(OUT_7C, "Pipeline7c_Fig7b_TF_LFC_Matrix.csv"),
              row.names = FALSE)

    # ★ v4.2: Summary JSON 含 ego 过滤信息
    fig7_summary <- list(
      version           = "v4.4",
      changes           = list(
        "v4.1: directional boxes, Fig7a/7b layout",
        "v4.2: ego-network topological constraint applied before visualization",
        "v4.3: post-cap connectivity re-check; network_max_nodes → network_max_edges",
        "v4.4: parameterized degree-1 control; post-cap degree re-pruning"
      ),
      ego_filter = list(
        applied       = network_data$ego_filter_applied,
        ego_steps     = network_data$ego_steps_used,
        nodes_removed = network_data$ego_n_removed,
        focus_gene    = PARAMS$focus_gene
      ),
      box_logic         = "genotype (Logic 2 LFC direction)",
      n_targets         = n_tgt,
      n_tfs             = n_tf,
      n_timepoints      = length(non_ctrl_times),
      timepoints        = non_ctrl_times,
      control_time      = ctrl_time,
      conditions        = c(cond_wt, cond_mut),
      target_star_count = sum(target_mats$stars != ""),
      target_box_count  = sum(target_mats$boxes),
      tf_star_count     = sum(tf_mats$stars != ""),
      tf_box_count      = sum(tf_mats$boxes),
      color_palette     = "reversed BrBG",
      fig7a_size        = c(fig_w_t, fig_h_t),
      fig7b_size        = c(fig_w_tf, fig_h_tf),
      target_clusters   = as.list(table(target_meta$WT_Cluster_label)),
      tf_out_degree_range = c(min(tf_meta$out_degree), max(tf_meta$out_degree))
    )

    writeLines(
      jsonlite::toJSON(fig7_summary, pretty = TRUE, auto_unbox = TRUE),
      file.path(OUT_7C, "Pipeline7c_Fig7_DE_Evidence_Summary.json")
    )

    cat(sprintf("    Exported: 2 LFC CSVs + 1 JSON\n"))

    "success"
  }, error = function(e) {
    cat(sprintf("    \u26A0 Evidence Heatmap failed: %s\n", e$message))
    cat(sprintf("    Traceback: %s\n",
                paste(capture.output(traceback()), collapse = "\n")))
    paste("failed:", e$message)
  })
}

# ========================================================================
# 最终报告
# ========================================================================

cat("\n================================================================================\n")
cat("  Pipeline 7c: Visualization Complete (v4.4)\n")
cat("================================================================================\n\n")

cat(sprintf("结束时间: %s\n\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))

cat("--- 图表生成结果 ---\n")
for (pname in names(plot_results)) {
  status <- plot_results[[pname]]
  icon <- ifelse(grepl("^success", status), "\u2713", "\u2717")
  cat(sprintf("  %s %s: %s\n", icon, pname, status))
}

n_success <- sum(grepl("^success", unlist(plot_results)))
n_total   <- length(plot_results)
cat(sprintf("\n总计: %d/%d 图表成功生成\n", n_success, n_total))
cat(sprintf("输出目录: %s\n", OUT_7C))

# ★ v4.2: 打印 ego 过滤摘要
if (!is.null(network_data) && isTRUE(network_data$ego_filter_applied)) {
  cat(sprintf("\n★ Ego-network filter: %s, removed %d nodes (focus_gene=%s)\n",
              ifelse(network_data$ego_steps_used > 0,
                     paste0(network_data$ego_steps_used, "-step ego"),
                     "connected component"),
              network_data$ego_n_removed,
              PARAMS$focus_gene))
}

# 保存 session info
writeLines(capture.output(sessionInfo()),
           file.path(OUT_7C, "Pipeline7c_sessionInfo.txt"))

# 关闭日志
sink(type = "message")
sink(type = "output")
close(log_con)

cat(sprintf("\n日志已保存: %s\n", log_file))
cat("Pipeline 7c 完成.\n")