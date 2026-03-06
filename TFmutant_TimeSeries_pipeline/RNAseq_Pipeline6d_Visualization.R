#!/usr/bin/env Rscript
################################################################################
# RNAseq Pipeline 6d v2.1: Publication-Quality Visualization
#
# Final visualization step that consumes outputs from:
#   - Pipeline 6b (RF Performance)
#   - Pipeline 6c v3.1 (Integrated Motif Annotation)
#
# Generates three classes of figures (modeled on Fig 4a-c):
#   Fig 1 — AUC-ROC bar plot per cluster (faceted WT vs Mutant)
#   Fig 2 — Pairwise sequence similarity heatmap (Overlap Coefficient)
#            Upper-triangle layout, WT block first, Mutant block second.
#   Fig 3 — Per-cluster scatter (Importance vs Enrichment) + DB-anchored logos
#
# v2.1 Fixes:
#   - Cluster display names: WT_Cx / MUT_Cx -> {display}_Cx via pre-built dict
#   - Heatmap changed to upper-triangle; axis label colors fixed
#   - Fig 3: fixed logical coercion error with target_keyword grepl
#   - Fig 1: replaced deprecated geom_errorbarh
#
# v2.2 Changes:
#   - Fig 2: Mutant block label repositioned to lower-triangle white space
#   - Fig 3: Column order changed to [pCRE] - [DB Match 1] - [DB Match 2]
#   - Rank in Selected_Motifs_List.csv now reflects True Rank by Importance
#   - New per-cluster output: Top20_Motif_Details.tsv
#   - target_keyword: supports comma-separated multi-keywords (e.g. "dof,WRKY")
#   - New --target_scope: controls which clusters get keyword highlighting
#
# Version : 2.2
# Date    : 2026-02
################################################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(optparse)
  library(universalmotif)
  library(ggseqlogo)
  library(ggrepel)
  library(grid)
  library(gridExtra)
  library(grDevices)
  library(RColorBrewer)
  library(scales)
})

################################################################################
# CLI Arguments
################################################################################

option_list <- list(
  make_option(c("--p6c_dir"), type = "character", default = NULL,
              help = "Pipeline 6c output root (contains WT/ and Mutant/ subfolders)"),
  make_option(c("--p6b_dir"), type = "character", default = NULL,
              help = "Pipeline 6b output root (contains WT/ and Mutant/ subfolders)"),
  make_option(c("--meme_db"), type = "character", default = NULL,
              help = "Path to MEME format motif database (.meme file)"),
  make_option(c("--output_dir"), type = "character", default = NULL,
              help = "Output directory [default: p6c_dir/../Pipeline6d_Visualization]"),
  make_option(c("--genotype"), type = "character", default = "both",
              help = "Genotype to process: 'WT', 'Mutant', or 'both' [default: %default]"),
  make_option(c("--wt_display_label"), type = "character", default = "WT",
              help = "Display label for WT genotype. Use italic(text) for italics [default: %default]"),
  make_option(c("--mut_display_label"), type = "character", default = "Mutant",
              help = "Display label for Mutant genotype. Use italic(text) for italics [default: %default]"),
  make_option(c("--wt_color"), type = "character", default = "#E69F00",
              help = "Hex color for WT genotype [default: %default]"),
  make_option(c("--mut_color"), type = "character", default = "#56B4E9",
              help = "Hex color for Mutant genotype [default: %default]"),
  make_option(c("--plot_top_n"), type = "integer", default = 3,
              help = "Number of top motifs to highlight in scatter plots [default: %default]"),
  make_option(c("--top_n_features"), type = "integer", default = 100,
              help = "Top N features by Importance for heatmap similarity [default: %default]"),
  make_option(c("--heatmap_metric"), type = "character", default = "overlap",
              help = "Similarity metric: 'overlap' (Simpson) or 'jaccard' [default: %default]"),
  make_option(c("--min_overlap"), type = "integer", default = 4,
              help = "min.overlap for compare_motifs alignment [default: %default]"),
  make_option(c("--target_keyword"), type = "character", default = NULL,
              help = "Comma-separated keywords to highlight (e.g. 'dof,WRKY,bHLH') [default: NULL]"),
  make_option(c("--target_scope"), type = "character", default = "ALL",
              help = "Scope for keyword highlighting: 'ALL', 'WT', 'MUT'/'Mutant', or comma-separated internal cluster IDs (e.g. 'WT_C1,MUT_C2') [default: %default]"),
  make_option(c("--fig_width"), type = "double", default = 10,
              help = "Default figure width in inches [default: %default]"),
  make_option(c("--fig_height"), type = "double", default = 8,
              help = "Default figure height in inches [default: %default]"),
  make_option(c("--seed"), type = "integer", default = 42,
              help = "Random seed [default: %default]")
)

opt_parser <- OptionParser(
  option_list  = option_list,
  description  = "Pipeline 6d v2.2: Publication-Quality Visualization"
)
opt <- parse_args(opt_parser)

# --- Validate ---
if (is.null(opt$p6c_dir)) stop("Error: --p6c_dir is required")
if (is.null(opt$p6b_dir)) stop("Error: --p6b_dir is required")

if (is.null(opt$output_dir)) {
  opt$output_dir <- file.path(dirname(opt$p6c_dir), "Pipeline6d_v2.0_Visualization")
}

genotypes_to_process <- if (tolower(opt$genotype) == "both") {
  c("WT", "Mutant")
} else {
  opt$genotype
}

set.seed(opt$seed)

# --- Create output directories ---
dir.create(opt$output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$output_dir, "Cluster_Reports"), recursive = TRUE, showWarnings = FALSE)

cat("================================================================================\n")
cat("  Pipeline 6d v2.2: Publication-Quality Visualization\n")
cat("================================================================================\n\n")
cat(sprintf("Start time       : %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("6b dir           : %s\n", opt$p6b_dir))
cat(sprintf("6c dir           : %s\n", opt$p6c_dir))
cat(sprintf("MEME DB          : %s\n", ifelse(is.null(opt$meme_db), "(not provided)", opt$meme_db)))
cat(sprintf("Output dir       : %s\n", opt$output_dir))
cat(sprintf("Genotypes        : %s\n", paste(genotypes_to_process, collapse = ", ")))
cat(sprintf("WT label         : %s\n", opt$wt_display_label))
cat(sprintf("Mutant label     : %s\n", opt$mut_display_label))
cat(sprintf("Top N (scatter)  : %d\n", opt$plot_top_n))
cat(sprintf("Top N (heatmap)  : %d\n", opt$top_n_features))
cat(sprintf("Target keyword   : %s\n", ifelse(is.null(opt$target_keyword), "(none)", opt$target_keyword)))
cat(sprintf("Target scope     : %s\n", opt$target_scope))
cat("\n")

################################################################################
# Utility: Label Parsing — italic(text) detection
################################################################################

parse_label <- function(label_str) {
  m <- regmatches(label_str, regexec("^italic\\((.+)\\)$", label_str))[[1]]
  if (length(m) == 2) {
    return(list(is_italic = TRUE, text = m[2], raw = label_str))
  }
  return(list(is_italic = FALSE, text = label_str, raw = label_str))
}

make_expr <- function(parsed) {
  if (parsed$is_italic) {
    bquote(italic(.(parsed$text)))
  } else {
    parsed$text
  }
}

make_label_grob <- function(parsed, fontsize = 12, col = "black", fontface = "plain") {
  if (parsed$is_italic) {
    grid::textGrob(parsed$text,
                   gp = grid::gpar(fontsize = fontsize, col = col, fontface = "italic"))
  } else {
    grid::textGrob(parsed$text,
                   gp = grid::gpar(fontsize = fontsize, col = col, fontface = fontface))
  }
}

# --- Parse the user-provided labels ---
wt_label_parsed  <- parse_label(opt$wt_display_label)
mut_label_parsed <- parse_label(opt$mut_display_label)

# --- Build consistent color map ---
geno_colors <- c("WT" = opt$wt_color, "Mutant" = opt$mut_color)

################################################################################
# Utility: Cluster Name Display Dictionary (ONE-TIME, GLOBAL)
#
# Upstream data uses WT_Cx / MUT_Cx format in Source_Cluster.
# We build a replacement function to display user-chosen labels.
# e.g., WT_C1 -> WT_C1 (unchanged), MUT_C5 -> dof_C5
################################################################################

#' Rename a cluster name vector for display purposes.
#' Replaces WT_/MUT_ prefixes with the user-specified display text.
rename_clusters_for_display <- function(x, wt_text, mut_text) {
  x <- gsub("^WT_",  paste0(wt_text, "_"),  x)
  x <- gsub("^MUT_", paste0(mut_text, "_"), x)
  x
}

#' Detect genotype from a raw cluster name (before renaming).
#' Returns "WT" or "Mutant".
detect_geno_from_cluster <- function(x) {
  ifelse(grepl("^WT_", x), "WT", "Mutant")
}

# The plain-text versions (without italic() wrapper)
WT_DISPLAY  <- wt_label_parsed$text
MUT_DISPLAY <- mut_label_parsed$text

cat(sprintf("Display mapping  : WT_Cx -> %s_Cx, MUT_Cx -> %s_Cx\n\n", WT_DISPLAY, MUT_DISPLAY))

################################################################################
# Utility: Safe File Reader
################################################################################

safe_read_tsv <- function(path, ...) {
  if (!file.exists(path)) {
    warning(sprintf("File not found: %s", path))
    return(NULL)
  }
  tryCatch(
    read.table(path, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
               check.names = FALSE, quote = "", fill = TRUE, ...),
    error = function(e) {
      warning(sprintf("Failed to read: %s — %s", path, e$message))
      NULL
    }
  )
}

################################################################################
# MODULE 1: Load All Data
################################################################################

load_all_data <- function(p6b_dir, p6c_dir, genotypes) {
  cat("\n--- Loading Data ---\n")

  all_perf  <- list()
  all_integ <- list()
  all_cons  <- list()
  all_pwms  <- list()

  for (geno in genotypes) {
    cat(sprintf("\n  [%s]\n", geno))

    # 6b: RF Performance Summary
    perf_path <- file.path(p6b_dir, geno, "results", "RF_Performance_Summary.txt")
    perf <- safe_read_tsv(perf_path)
    if (!is.null(perf)) {
      perf$Genotype <- geno
      all_perf[[geno]] <- perf
      cat(sprintf("    6b Performance : %d clusters\n", nrow(perf)))
    }

    # 6c: Integrated Features
    integ_path <- file.path(p6c_dir, geno, "Integrated_Features_Full.txt")
    integ <- safe_read_tsv(integ_path)
    if (!is.null(integ)) {
      if (!"Genotype" %in% names(integ)) integ$Genotype <- geno
      all_integ[[geno]] <- integ
      cat(sprintf("    6c Integrated  : %d rows\n", nrow(integ)))
    }

    # 6c: Consensus Motif Summary
    cons_path <- file.path(p6c_dir, geno, "Consensus_Motif_Summary.txt")
    cons <- safe_read_tsv(cons_path)
    if (!is.null(cons)) {
      if (!"Genotype" %in% names(cons)) cons$Genotype <- geno
      all_cons[[geno]] <- cons
      cat(sprintf("    6c Consensus   : %d families\n", nrow(cons)))
    }

    # 6c: PWMs
    pwm_path <- file.path(p6c_dir, geno, "Motif_PWMs.rds")
    if (file.exists(pwm_path)) {
      pwm <- readRDS(pwm_path)
      all_pwms[[geno]] <- pwm
      cat(sprintf("    6c PWMs        : %d objects\n", length(pwm)))
    }
  }

  list(
    perf  = bind_rows(all_perf),
    integ = bind_rows(all_integ),
    cons  = bind_rows(all_cons),
    pwms  = all_pwms
  )
}

################################################################################
# FIGURE 1: AUC-ROC Bar Plot  (matches Fig 4a)
#
# Layout: WT on Left, Mutant on Right.
# Bar labels: AUC values (e.g., "0.853"), NOT feature counts.
# Cluster names: use display labels.
################################################################################

plot_auc_barplot <- function(perf_df, integ_df, output_dir,
                             wt_parsed, mut_parsed, geno_colors,
                             fig_width = 10, fig_height = 8) {
  cat("\n--- Figure 1: AUC-ROC Bar Plot ---\n")

  if (is.null(perf_df) || nrow(perf_df) == 0) {
    cat("  WARNING: No performance data, skipping AUC plot.\n")
    return(invisible(NULL))
  }

  plot_df <- perf_df

  # --- Rename cluster names for display ---
  plot_df$Display_Cluster <- rename_clusters_for_display(
    plot_df$Cluster, WT_DISPLAY, MUT_DISPLAY
  )

  # --- AUC value label ---
  plot_df$auc_label <- sprintf("%.3f", plot_df$Mean_AUC)

  # --- Natural cluster sort ---
  extract_num <- function(x) {
    nums <- as.numeric(gsub("[^0-9]", "", x))
    ifelse(is.na(nums), 0, nums)
  }
  plot_df <- plot_df %>%
    mutate(cluster_order = extract_num(Cluster)) %>%
    arrange(Genotype, cluster_order)

  plot_df$Display_Cluster <- factor(
    plot_df$Display_Cluster,
    levels = unique(plot_df$Display_Cluster[order(plot_df$cluster_order)])
  )

  # --- Force facet order: WT first (left), Mutant second (right) ---
  plot_df$Genotype <- factor(plot_df$Genotype, levels = c("WT", "Mutant"))

  # --- Build facet labeller with italic support ---
  geno_label_map <- c(
    "WT"     = if (wt_parsed$is_italic)
                 sprintf("italic(\"%s\")", wt_parsed$text) else wt_parsed$text,
    "Mutant" = if (mut_parsed$is_italic)
                 sprintf("italic(\"%s\")", mut_parsed$text) else mut_parsed$text
  )
  facet_labeller <- as_labeller(geno_label_map, default = label_parsed)

  p <- ggplot(plot_df, aes(x = Mean_AUC, y = Display_Cluster, fill = Genotype)) +
    geom_col(width = 0.7, show.legend = FALSE) +
    geom_errorbar(aes(xmin = pmax(0, Mean_AUC - SD_AUC),
                      xmax = pmin(1, Mean_AUC + SD_AUC)),
                  width = 0.25, linewidth = 0.3, color = "grey30",
                  orientation = "y") +
    geom_text(aes(x = pmin(1, Mean_AUC + SD_AUC) + 0.02, label = auc_label),
              hjust = 0, size = 2.8, color = "grey20", fontface = "bold") +
    geom_vline(xintercept = 0.5, linetype = "dashed", color = "red",
               linewidth = 0.5, alpha = 0.7) +
    geom_vline(xintercept = 0.7, linetype = "dashed", color = "grey50",
               linewidth = 0.4, alpha = 0.6) +
    facet_wrap(~ Genotype, scales = "free_y", ncol = 2,
               labeller = facet_labeller) +
    scale_fill_manual(values = geno_colors) +
    scale_x_continuous(limits = c(0, 1.15), breaks = seq(0, 1, 0.25),
                       expand = c(0, 0)) +
    labs(x = "AUC-ROC", y = "Temporally dynamic groups") +
    theme_classic(base_size = 12) +
    theme(
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text       = element_text(size = 11),
      panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
      axis.text.y      = element_text(size = 9),
      plot.margin       = margin(10, 30, 10, 10)
    )

  out_path <- file.path(output_dir, "Fig1_AUC_ROC_Barplot.pdf")
  ggsave(out_path, p, width = fig_width, height = fig_height, dpi = 300)
  cat(sprintf("  Saved: %s\n", out_path))
  ggsave(sub("\\.pdf$", ".png", out_path), p,
         width = fig_width, height = fig_height, dpi = 300)
  return(p)
}

################################################################################
# FIGURE 2: Sequence Similarity Heatmap — Upper-Triangle "Orthogonal" Layout
#
# - Uses Top N features by Importance per cluster (not full feature set).
# - Row/Col order: All WT clusters first, then All Mutant clusters.
#   Top-Left Block:  WT intra-group comparison.
#   Top-Right Block: WT vs Mutant inter-group (orthogonal) comparison.
#   Bottom-Right:    Mutant intra-group comparison.
# - Visual style: UPPER-triangle only (matching reference Fig 4b).
# - Cluster names use display labels; axis colors match genotype.
################################################################################

overlap_coefficient <- function(set_a, set_b) {
  inter <- length(intersect(set_a, set_b))
  min_size <- min(length(set_a), length(set_b))
  if (min_size == 0) return(0)
  inter / min_size
}

plot_similarity_heatmap <- function(integ_df, output_dir,
                                    metric = "overlap",
                                    top_n_features = 100,
                                    wt_parsed, mut_parsed, geno_colors,
                                    fig_width = 10, fig_height = 9) {
  cat("\n--- Figure 2: Sequence Similarity Heatmap ---\n")

  if (is.null(integ_df) || nrow(integ_df) == 0) {
    cat("  WARNING: No integrated data, skipping heatmap.\n")
    return(invisible(NULL))
  }

  # --- Identify columns ---
  kmer_col <- if ("Kmer" %in% names(integ_df)) "Kmer" else
              if ("kmer" %in% names(integ_df)) "kmer" else NULL
  imp_col  <- if ("Importance" %in% names(integ_df)) "Importance" else
              if ("mean_importance" %in% names(integ_df)) "mean_importance" else NULL

  if (is.null(kmer_col)) {
    cat("  WARNING: No Kmer column found, skipping heatmap.\n")
    return(invisible(NULL))
  }

  # --- Build per-cluster kmer sets using ONLY Top N by Importance ---
  if (!is.null(imp_col)) {
    cluster_sets <- integ_df %>%
      group_by(Genotype, Source_Cluster) %>%
      arrange(desc(.data[[imp_col]])) %>%
      slice_head(n = top_n_features) %>%
      summarise(kmers = list(unique(.data[[kmer_col]])), .groups = "drop")
    cat(sprintf("  Using Top %d features per cluster for similarity.\n", top_n_features))
  } else {
    cluster_sets <- integ_df %>%
      group_by(Genotype, Source_Cluster) %>%
      summarise(kmers = list(unique(.data[[kmer_col]])), .groups = "drop")
  }

  n_clusters <- nrow(cluster_sets)
  cat(sprintf("  Clusters: %d\n", n_clusters))
  if (n_clusters < 2) {
    cat("  WARNING: Need at least 2 clusters for heatmap.\n")
    return(invisible(NULL))
  }

  # --- Sort: WT clusters first, then Mutant, natural numeric within ---
  cluster_sets <- cluster_sets %>%
    mutate(
      geno_rank = ifelse(Genotype == "WT", 0, 1),
      clust_num = as.numeric(gsub("[^0-9]", "", Source_Cluster))
    ) %>%
    arrange(geno_rank, clust_num)

  raw_labels <- cluster_sets$Source_Cluster   # e.g. WT_C1, MUT_C5
  geno_vec   <- cluster_sets$Genotype         # WT or Mutant
  kmer_lists <- cluster_sets$kmers

  # Display labels: WT_C1 -> WT_C1, MUT_C1 -> dof_C1
  disp_labels <- rename_clusters_for_display(raw_labels, WT_DISPLAY, MUT_DISPLAY)

  # --- Genotype color for each cluster (aligned to sorted order) ---
  geno_color_vec <- ifelse(geno_vec == "WT",
                           unname(geno_colors["WT"]),
                           unname(geno_colors["Mutant"]))

  # --- Compute pairwise similarity matrix ---
  sim_mat <- matrix(0, nrow = n_clusters, ncol = n_clusters,
                    dimnames = list(disp_labels, disp_labels))

  for (i in seq_len(n_clusters)) {
    sim_mat[i, i] <- 1.0
    if (i < n_clusters) {
      for (j in (i + 1):n_clusters) {
        if (metric == "jaccard") {
          inter <- length(intersect(kmer_lists[[i]], kmer_lists[[j]]))
          union_size <- length(union(kmer_lists[[i]], kmer_lists[[j]]))
          val <- if (union_size == 0) 0 else inter / union_size
        } else {
          val <- overlap_coefficient(kmer_lists[[i]], kmer_lists[[j]])
        }
        sim_mat[i, j] <- val
        sim_mat[j, i] <- val
      }
    }
  }

  # --- UPPER-triangle: set lower triangle to NA ---
  sim_upper <- sim_mat
  sim_upper[lower.tri(sim_upper)] <- NA

  # --- Build long-format for ggplot ---
  melted <- expand.grid(Row = seq_len(n_clusters), Col = seq_len(n_clusters))
  melted$value <- sapply(seq_len(nrow(melted)),
                         function(k) sim_upper[melted$Row[k], melted$Col[k]])
  melted <- melted %>% filter(!is.na(value))

  melted$Row_label <- factor(disp_labels[melted$Row], levels = disp_labels)
  melted$Col_label <- factor(disp_labels[melted$Col], levels = disp_labels)

  # --- Block counts ---
  n_wt  <- sum(geno_vec == "WT")
  n_mut <- sum(geno_vec == "Mutant")

  # --- Axis label colors ---
  # X-axis (top, left to right) follows scale_x_discrete(limits = disp_labels):
  #   ggplot applies the colour vector in the same left->right order.
  x_label_colors <- geno_color_vec

  # Y-axis (right side) with scale_y_discrete(limits = rev(disp_labels)):
  #   Display from top to bottom = rev(disp_labels)[1], ..., rev(disp_labels)[n]
  #   ggplot applies element_text colour from BOTTOM to TOP (i.e., in limits order).
  #   limits order = rev(disp_labels) = [last_cluster, ..., first_cluster]
  #   So colour vector should follow rev(disp_labels) order = rev(geno_color_vec).
  y_label_colors <- rev(geno_color_vec)

  p <- ggplot(melted, aes(x = Col_label, y = Row_label, fill = value)) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradientn(
      colours = c("white", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B", "#67001F"),
      limits  = c(0, 1),
      breaks  = seq(0, 1, 0.2),
      name    = "Correlation"
    ) +
    scale_y_discrete(limits = rev(disp_labels), position = "right") +
    scale_x_discrete(limits = disp_labels, position = "top") +
    coord_fixed() +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x.top   = element_text(angle = 45, hjust = 0, vjust = 0.5,
                                        size = 8, color = x_label_colors),
      axis.text.y.right  = element_text(size = 8, color = y_label_colors),
      axis.title         = element_blank(),
      panel.grid         = element_blank(),
      legend.position    = c(0.15, 0.25),
      legend.direction   = "vertical",
      legend.key.height  = unit(0.8, "cm"),
      legend.key.width   = unit(0.4, "cm"),
      legend.title       = element_text(size = 9),
      legend.text        = element_text(size = 8),
      plot.margin        = margin(10, 10, 30, 10)
    )

  # --- Genotype block rectangles ---
  # Upper-triangle: data visible above diagonal.
  # Matrix row 1 = top of plot (WT_C1), row n = bottom (last MUT).
  # In ggplot tile coords with reversed y:
  #   y-value for row i (1-based) = n_clusters - i + 1 (due to rev)
  #   x-value for col j = j
  # WT block:  rows 1..n_wt, cols 1..n_wt   => top-left
  # MUT block: rows (n_wt+1)..n, cols (n_wt+1)..n => bottom-right
  if (n_wt > 0 && n_mut > 0) {
    p <- p +
      # WT block (top-left)
      annotate("rect",
               xmin = 0.5, xmax = n_wt + 0.5,
               ymin = n_clusters - n_wt + 0.5, ymax = n_clusters + 0.5,
               fill = NA, color = geno_colors["WT"], linewidth = 1.2) +
      # Mutant block (bottom-right)
      annotate("rect",
               xmin = n_wt + 0.5, xmax = n_clusters + 0.5,
               ymin = 0.5, ymax = n_mut + 0.5,
               fill = NA, color = geno_colors["Mutant"], linewidth = 1.2)

    # Block labels
    p <- p +
      annotate("text",
               x = (n_wt + 1) / 2,
               y = n_clusters - n_wt + 0.5 - 0.5,
               label = WT_DISPLAY, color = geno_colors["WT"],
               fontface = if (wt_parsed$is_italic) "italic" else "plain",
               size = 3.5, hjust = 0.5, vjust = 1) +
      annotate("text",
               x = n_wt,              # x坐标: 放在 WT 最后一个 cluster 的位置 (方框左侧空白区)
               y = (n_mut + 1) / 2,   # y坐标: 计算 Mutant 方框的垂直中心高度
               label = MUT_DISPLAY, 
               color = geno_colors["Mutant"],
               fontface = if (mut_parsed$is_italic) "italic" else "plain",
               size = 3.5, 
               hjust = 1,             # 右对齐: 让文字尾部靠近方框左边界
               vjust = 0.5)           # 垂直居中
  }

  out_path <- file.path(output_dir, "Fig2_Similarity_Heatmap.pdf")
  ggsave(out_path, p, width = fig_width, height = fig_height, dpi = 300)
  cat(sprintf("  Saved: %s\n", out_path))
  ggsave(sub("\\.pdf$", ".png", out_path), p,
         width = fig_width, height = fig_height, dpi = 300)

  return(sim_mat)
}

################################################################################
# FIGURE 3: Per-Cluster Scatter + DB-Anchored Motif Logos
#
# Logo Alignment: "Database Anchor" mode.
#   - DB motif is the REFERENCE (drawn fixed, no RC, no shift).
#   - pCRE motif is the QUERY: RC and shifted/padded to align to DB.
#
# Visual: ggseqlogo(method="bits"), Y 0-2, no color legend.
# Target Mode: if --target_keyword provided, highlight matching motif.
# Outputs: PDF report + Selected_Motifs_List.csv
#
# v2.1 FIX: replaced `&&` with vectorized `&` in grepl logic to avoid
#           'length = N in coercion to logical(1)' error.
################################################################################

load_meme_db_safe <- function(meme_db_path) {
  if (is.null(meme_db_path) || !file.exists(meme_db_path)) {
    cat("  MEME DB not available; DB logos will be skipped.\n")
    return(NULL)
  }
  tryCatch({
    db <- read_meme(meme_db_path)
    if (inherits(db, "universalmotif")) db <- list(db)
    names(db) <- sapply(db, function(m) m@name)
    cat(sprintf("  MEME DB loaded: %d motifs\n", length(db)))
    db
  }, error = function(e) {
    warning(sprintf("Failed to load MEME DB: %s", e$message))
    NULL
  })
}

get_db_motif <- function(db, name) {
  if (is.null(db) || is.null(name) || is.na(name)) return(NULL)
  if (name %in% names(db)) return(db[[name]])
  idx <- grep(name, names(db), fixed = TRUE)
  if (length(idx) > 0) return(db[[idx[1]]])
  NULL
}

align_pcre_to_db <- function(db_motif, pcre_motif, min_overlap = 4) {
  if (is.null(db_motif) || is.null(pcre_motif)) return(pcre_motif)

  tryCatch({
    pcc_fwd <- tryCatch(
      compare_motifs(list(db_motif, pcre_motif),
                     method = "PCC", min.overlap = min_overlap,
                     tryRC = FALSE, nthreads = 1)[1, 2],
      error = function(e) -1
    )

    pcre_rc <- motif_rc(pcre_motif)
    pcc_rc <- tryCatch(
      compare_motifs(list(db_motif, pcre_rc),
                     method = "PCC", min.overlap = min_overlap,
                     tryRC = FALSE, nthreads = 1)[1, 2],
      error = function(e) -1
    )

    best_pcre <- if (pcc_rc > pcc_fwd) pcre_rc else pcre_motif
    db_len    <- ncol(db_motif@motif)
    pcre_len  <- ncol(best_pcre@motif)

    best_pcc   <- -1
    best_shift <- 0
    max_neg <- -(pcre_len - min_overlap)
    max_pos <- db_len - min_overlap

    for (shift in max_neg:max_pos) {
      pcre_start <- max(1, 1 - shift)
      pcre_end   <- min(pcre_len, db_len - shift)
      db_start   <- max(1, 1 + shift)
      db_end     <- min(db_len, pcre_len + shift)
      if (pcre_end < pcre_start || db_end < db_start) next
      overlap_len <- pcre_end - pcre_start + 1
      if (overlap_len < min_overlap) next
      pcc_vals <- sapply(seq_len(overlap_len), function(k) {
        p_col <- best_pcre@motif[, pcre_start + k - 1]
        d_col <- db_motif@motif[, db_start + k - 1]
        cor(p_col, d_col)
      })
      mean_pcc <- mean(pcc_vals, na.rm = TRUE)
      if (mean_pcc > best_pcc) {
        best_pcc   <- mean_pcc
        best_shift <- shift
      }
    }

    total_width <- db_len
    if (best_shift > 0) {
      pad_left <- best_shift
      new_mat  <- best_pcre@motif
      if (pad_left > 0) new_mat <- cbind(matrix(0.25, 4, pad_left), new_mat)
      pad_right <- max(0, total_width - ncol(new_mat))
      if (pad_right > 0) new_mat <- cbind(new_mat, matrix(0.25, 4, pad_right))
      new_mat <- new_mat[, seq_len(min(ncol(new_mat), total_width)), drop = FALSE]
    } else {
      new_mat   <- matrix(0.25, 4, total_width)
      dst_start <- max(1, 1 + best_shift)
      src_start <- max(1, 1 - best_shift)
      n_copy    <- min(pcre_len - src_start + 1, total_width - dst_start + 1)
      if (n_copy > 0) {
        new_mat[, dst_start:(dst_start + n_copy - 1)] <-
          best_pcre@motif[, src_start:(src_start + n_copy - 1)]
      }
    }

    rownames(new_mat) <- c("A", "C", "G", "T")
    aligned_pcre <- create_motif(new_mat, name = best_pcre@name,
                                 type = "PPM", alphabet = "DNA")
    return(aligned_pcre)
  }, error = function(e) {
    return(pcre_motif)
  })
}

motif_to_matrix <- function(motif_obj) {
  if (is.null(motif_obj)) return(NULL)
  mat <- motif_obj@motif
  if (is.null(rownames(mat))) rownames(mat) <- c("A", "C", "G", "T")
  mat
}

plot_logo_bits <- function(motif_obj, title = NULL) {
  if (is.null(motif_obj)) {
    return(grid::textGrob("N/A", gp = grid::gpar(col = "grey60", fontsize = 10)))
  }
  tryCatch({
    mat <- motif_to_matrix(motif_obj)
    p <- ggseqlogo(mat, method = "bits") +
      ylim(0, 2) +
      theme_minimal(base_size = 8) +
      theme(
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        axis.title.x    = element_blank(),
        legend.position = "none",
        plot.title      = element_text(size = 7, face = "bold", hjust = 0.5),
        plot.margin     = margin(2, 4, 2, 4)
      )
    if (!is.null(title)) p <- p + ggtitle(title)
    return(ggplotGrob(p))
  }, error = function(e) {
    return(grid::textGrob(paste0("Error: ", e$message),
                          gp = grid::gpar(col = "red", fontsize = 7)))
  })
}

kmer_to_motif <- function(kmer, name = kmer, nsites = 100) {
  bases <- strsplit(toupper(kmer), "")[[1]]
  len   <- length(bases)
  mat   <- matrix(0, nrow = 4, ncol = len,
                  dimnames = list(c("A", "C", "G", "T"), NULL))
  for (j in seq_len(len)) {
    b <- bases[j]
    if (b %in% c("A", "C", "G", "T")) {
      mat[b, j] <- 1.0
    } else {
      mat[, j] <- 0.25
    }
  }
  create_motif(mat, name = name, type = "PPM", alphabet = "DNA",
               nsites = as.integer(max(1, round(nsites))))
}

generate_cluster_reports <- function(integ_df, pwm_lists, meme_db,
                                     output_dir, top_n = 3,
                                     min_overlap = 4,
                                     target_keyword = NULL,
                                     target_scope = "ALL",
                                     wt_parsed, mut_parsed, geno_colors) {
  cat("\n--- Figure 3: Per-Cluster Reports ---\n")

  if (is.null(integ_df) || nrow(integ_df) == 0) {
    cat("  WARNING: No integrated data, skipping cluster reports.\n")
    return(invisible(NULL))
  }

  imp_col <- if ("Importance" %in% names(integ_df)) "Importance" else
             if ("mean_importance" %in% names(integ_df)) "mean_importance" else NULL
  if (is.null(imp_col)) {
    cat("  WARNING: No Importance column found, skipping.\n")
    return(invisible(NULL))
  }

  if ("Enrichment_LogP" %in% names(integ_df)) {
    enrich_col <- "Enrichment_LogP"
  } else if ("Enrichment_P" %in% names(integ_df)) {
    integ_df$Enrichment_LogP <- -log10(pmax(integ_df$Enrichment_P, 1e-300))
    integ_df$Enrichment_LogP[is.na(integ_df$Enrichment_P)] <- NA
    enrich_col <- "Enrichment_LogP"
  } else {
    enrich_col <- NULL
  }

  kmer_col <- if ("Kmer" %in% names(integ_df)) "Kmer" else
              if ("kmer" %in% names(integ_df)) "kmer" else NULL
  cons_col <- if ("Consensus_Sequence" %in% names(integ_df)) "Consensus_Sequence" else NULL

  # Safe column detection — intersect()[1] can return character(0) not NA
  match1_col <- intersect(c("Kmer_DB_Match_1", "DB_Match_1"), names(integ_df))
  match2_col <- intersect(c("Kmer_DB_Match_2", "DB_Match_2"), names(integ_df))
  pcc1_col   <- intersect(c("Kmer_PCC_1", "PCC_1"), names(integ_df))
  pcc2_col   <- intersect(c("Kmer_PCC_2", "PCC_2"), names(integ_df))
  match1_col <- if (length(match1_col) > 0) match1_col[1] else NA_character_
  match2_col <- if (length(match2_col) > 0) match2_col[1] else NA_character_
  pcc1_col   <- if (length(pcc1_col) > 0)   pcc1_col[1]   else NA_character_
  pcc2_col   <- if (length(pcc2_col) > 0)   pcc2_col[1]   else NA_character_

  selected_motifs_list <- list()
  clusters <- unique(paste0(integ_df$Genotype, "||", integ_df$Source_Cluster))
  cat(sprintf("  Processing %d clusters...\n", length(clusters)))

  # --- Parse multi-keyword (v2.2) ---
  keyword_vec <- NULL
  if (!is.null(target_keyword) && nchar(target_keyword) > 0) {
    keyword_vec <- tolower(trimws(strsplit(target_keyword, ",")[[1]]))
    keyword_vec <- keyword_vec[nchar(keyword_vec) > 0]
    if (length(keyword_vec) == 0) keyword_vec <- NULL
    cat(sprintf("  Target keywords: %s\n", paste(keyword_vec, collapse = ", ")))
  }

  # --- Scope checker (v2.2) ---
  # Uses INTERNAL cluster IDs (e.g. WT_C1, MUT_C1), NOT display labels
  check_scope <- function(geno, cl_name, scope_str) {
    scope_upper <- toupper(trimws(scope_str))
    if (scope_upper == "ALL") return(TRUE)
    if (scope_upper == "WT")  return(geno == "WT")
    if (scope_upper %in% c("MUT", "MUTANT")) return(geno == "Mutant")
    # Specific list: comma-separated internal IDs like "WT_C1,MUT_C2"
    scope_ids <- trimws(strsplit(scope_str, ",")[[1]])
    # Build internal ID: Source_Cluster already has prefix (WT_C1, MUT_C1)
    # But also allow Genotype_SourceCluster pattern (e.g. WT_WT_C1 is unlikely,
    # so just check cl_name directly against the list)
    return(cl_name %in% scope_ids)
  }

  cat(sprintf("  Target scope: %s\n", target_scope))

  report_dir <- file.path(output_dir, "Cluster_Reports")

  for (cl_key in clusters) {
    parts   <- strsplit(cl_key, "\\|\\|")[[1]]
    geno    <- parts[1]
    cl_name <- parts[2]

    # Display names
    cl_display   <- rename_clusters_for_display(cl_name, WT_DISPLAY, MUT_DISPLAY)
    geno_parsed  <- if (geno == "WT") wt_parsed else mut_parsed
    geno_display <- geno_parsed$text

    safe_label <- paste0(geno, "_", cl_name)
    cat(sprintf("    [%s -> %s] ", safe_label, cl_display))

    cl_dir <- file.path(report_dir, safe_label)
    dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)

    cl_data <- integ_df %>% filter(Genotype == geno, Source_Cluster == cl_name)
    if (nrow(cl_data) == 0) { cat("empty, skip.\n"); next }

    # --- Compute True Rank by Importance (v2.2) ---
    cl_data <- cl_data %>%
      arrange(desc(.data[[imp_col]])) %>%
      mutate(True_Rank = row_number())

    # --- Export Top20_Motif_Details.tsv (v2.2) ---
    top20_data <- cl_data %>% head(20)
    top20_out <- data.frame(
      Genotype        = geno,
      Cluster         = cl_name,
      True_Rank       = top20_data$True_Rank,
      pCRE_Kmer       = if (!is.null(kmer_col) && kmer_col %in% names(top20_data))
                           top20_data[[kmer_col]] else NA,
      Importance      = top20_data[[imp_col]],
      Enrichment_LogP = if (!is.null(enrich_col) && enrich_col %in% names(top20_data))
                           top20_data[[enrich_col]] else NA,
      DB_Match_1      = if (!is.na(match1_col) && match1_col %in% names(top20_data))
                           top20_data[[match1_col]] else NA,
      PCC_1           = if (!is.na(pcc1_col) && pcc1_col %in% names(top20_data))
                           top20_data[[pcc1_col]] else NA,
      DB_Match_2      = if (!is.na(match2_col) && match2_col %in% names(top20_data))
                           top20_data[[match2_col]] else NA,
      PCC_2           = if (!is.na(pcc2_col) && pcc2_col %in% names(top20_data))
                           top20_data[[pcc2_col]] else NA,
      stringsAsFactors = FALSE
    )
    top20_path <- file.path(cl_dir, "Top20_Motif_Details.tsv")
    write.table(top20_out, top20_path, sep = "\t", row.names = FALSE, quote = FALSE)

    # ====================================================================
    # 3A: Scatter Plot
    # ====================================================================
    scatter_path <- file.path(cl_dir, "Scatter_Importance_Enrichment.pdf")
    has_enrich <- !is.null(enrich_col) && any(!is.na(cl_data[[enrich_col]]))

    # --- Target keyword search (v2.2: multi-keyword + scope) ---
    target_row_idx <- NULL
    if (!is.null(keyword_vec) && check_scope(geno, cl_name, target_scope)) {
      n_rows_cl <- nrow(cl_data)
      hits1 <- rep(FALSE, n_rows_cl)
      hits2 <- rep(FALSE, n_rows_cl)

      if (!is.na(match1_col) && match1_col %in% names(cl_data)) {
        vals1 <- tolower(as.character(cl_data[[match1_col]]))
        for (kw in keyword_vec) {
          hits1 <- hits1 | (!is.na(vals1) & grepl(kw, vals1, fixed = TRUE))
        }
      }
      if (!is.na(match2_col) && match2_col %in% names(cl_data)) {
        vals2 <- tolower(as.character(cl_data[[match2_col]]))
        for (kw in keyword_vec) {
          hits2 <- hits2 | (!is.na(vals2) & grepl(kw, vals2, fixed = TRUE))
        }
      }
      match_hits <- which(hits1 | hits2)
      if (length(match_hits) > 0) {
        target_row_idx <- match_hits[which.max(cl_data[[imp_col]][match_hits])]
      }
    }

    if (has_enrich) {
      top_idx <- order(cl_data[[imp_col]], decreasing = TRUE)
      top_idx <- head(top_idx, top_n)
      cl_data$point_type <- "other"
      cl_data$point_type[top_idx] <- "top"
      if (!is.null(target_row_idx)) {
        cl_data$point_type[target_row_idx] <- "target"
        top_idx <- unique(c(top_idx, target_row_idx))
      }

      cl_data$point_label <- ""
      if (!is.null(cons_col) && any(!is.na(cl_data[[cons_col]]))) {
        cl_data$point_label[top_idx] <- ifelse(
          !is.na(cl_data[[cons_col]][top_idx]),
          cl_data[[cons_col]][top_idx],
          cl_data[[kmer_col]][top_idx]
        )
      } else if (!is.null(kmer_col)) {
        cl_data$point_label[top_idx] <- cl_data[[kmer_col]][top_idx]
      }

      plot_data <- cl_data %>%
        mutate(x_val = .data[[imp_col]], y_val = .data[[enrich_col]]) %>%
        filter(!is.na(x_val), !is.na(y_val))

      point_colors <- c("other" = "grey60", "top" = "red", "target" = "#7B2D8E")
      point_sizes  <- c("other" = 1.5, "top" = 3, "target" = 3.5)

      if (nrow(plot_data) >= 3) {
        p_scatter <- ggplot(plot_data, aes(x = x_val, y = y_val)) +
          geom_point(aes(color = point_type, size = point_type), alpha = 0.7) +
          geom_smooth(method = "lm", formula = y ~ x,
                      color = "grey40", fill = "grey85",
                      linewidth = 0.6, alpha = 0.3, se = TRUE) +
          scale_color_manual(values = point_colors, guide = "none") +
          scale_size_manual(values = point_sizes, guide = "none") +
          geom_text_repel(
            data = filter(plot_data, point_label != ""),
            aes(label = point_label, color = point_type),
            size = 3, box.padding = 0.5, point.padding = 0.3,
            max.overlaps = 20, seed = 42, segment.size = 0.3,
            show.legend = FALSE
          ) +
          labs(
            x     = "Importance (RF Mean Decrease Gini)",
            y     = expression(-log[10](P[enrichment])),
            title = cl_display
          ) +
          theme_classic(base_size = 11) +
          theme(
            plot.title  = element_text(size = 12, face = "bold"),
            plot.margin = margin(10, 10, 10, 10)
          )
        ggsave(scatter_path, p_scatter, width = 6, height = 5, dpi = 300)
      } else {
        pdf(scatter_path, width = 6, height = 5)
        plot.new()
        text(0.5, 0.5,
             sprintf("%s: insufficient data (n=%d)", cl_display, nrow(plot_data)),
             cex = 1.2)
        dev.off()
      }
    } else {
      top_data <- cl_data %>%
        arrange(desc(.data[[imp_col]])) %>%
        head(min(20, nrow(cl_data)))
      bar_label <- if (!is.null(kmer_col)) kmer_col else "Feature_ID"
      if (!bar_label %in% names(top_data)) bar_label <- names(top_data)[1]
      top_data$rank_label <- factor(top_data[[bar_label]],
                                    levels = rev(top_data[[bar_label]]))

      p_bar <- ggplot(top_data, aes(x = .data[[imp_col]], y = rank_label)) +
        geom_col(fill = geno_colors[geno], width = 0.7) +
        labs(x = "Importance", y = NULL,
             title = paste0(cl_display, " (no enrichment data)")) +
        theme_classic(base_size = 10)
      ggsave(scatter_path, p_bar, width = 6, height = 5, dpi = 300)
    }

    # ====================================================================
    # 3B: Motif Logo Alignment — DB-Anchored Mode
    # ====================================================================
    logo_path <- file.path(cl_dir, "Motif_Alignments_TopN.pdf")

    top_rows <- cl_data %>% arrange(desc(.data[[imp_col]])) %>% head(top_n)
    if (!is.null(target_row_idx)) {
      target_row <- cl_data[target_row_idx, , drop = FALSE]
      top_rows <- bind_rows(top_rows, target_row) %>% distinct()
    }

    geno_pwms <- pwm_lists[[geno]]
    if (is.null(geno_pwms) || length(geno_pwms) == 0) cat("no PWMs, ")

    grob_list <- list()

    for (r in seq_len(nrow(top_rows))) {
      row_data <- top_rows[r, ]

      # Get pCRE motif
      family_label <- if ("Motif_Family_Label" %in% names(row_data))
                        row_data$Motif_Family_Label else NA
      pcre_motif <- NULL
      if (!is.na(family_label) && !is.null(geno_pwms) &&
          family_label %in% names(geno_pwms)) {
        pcre_motif <- geno_pwms[[family_label]]
      } else if (!is.null(kmer_col) && !is.na(row_data[[kmer_col]])) {
        pcre_motif <- tryCatch(
          kmer_to_motif(row_data[[kmer_col]], name = row_data[[kmer_col]]),
          error = function(e) NULL
        )
      }

      # Get DB Match 1 (REFERENCE)
      m1_name <- NA_character_
      m1_pcc  <- NA_real_
      if (!is.na(match1_col) && match1_col %in% names(row_data)) {
        m1_name <- as.character(row_data[[match1_col]])
        if (is.na(m1_name) || nchar(m1_name) == 0) m1_name <- NA_character_
      }
      if (!is.na(pcc1_col) && pcc1_col %in% names(row_data)) {
        m1_pcc <- suppressWarnings(as.numeric(row_data[[pcc1_col]]))
      }

      db_m1 <- NULL
      if (!is.na(m1_name) && !is.null(meme_db)) {
        db_m1 <- get_db_motif(meme_db, m1_name)
      }

      # Align pCRE to DB
      aligned_pcre <- pcre_motif
      if (!is.null(db_m1) && !is.null(pcre_motif)) {
        aligned_pcre <- align_pcre_to_db(db_m1, pcre_motif, min_overlap)
      }

      # Column 1: DB Match 1 (fixed reference)
      m1_title <- "DB: N/A"
      if (!is.na(m1_name)) {
        m1_title <- paste0("DB: ", m1_name)
        if (!is.na(m1_pcc)) m1_title <- paste0(m1_title, " (PCC=", round(m1_pcc, 3), ")")
      }
      grob_db1 <- plot_logo_bits(db_m1, title = m1_title)

      # Column 2: pCRE (aligned query)
      pcre_title <- "pCRE"
      if (!is.null(kmer_col) && !is.na(row_data[[kmer_col]]))
        pcre_title <- paste0("pCRE: ", row_data[[kmer_col]])
      grob_pcre <- plot_logo_bits(aligned_pcre, title = pcre_title)

      # Column 3: DB Match 2
      m2_name <- NA_character_
      m2_pcc  <- NA_real_
      if (!is.na(match2_col) && match2_col %in% names(row_data)) {
        m2_name <- as.character(row_data[[match2_col]])
        if (is.na(m2_name) || nchar(m2_name) == 0) m2_name <- NA_character_
      }
      if (!is.na(pcc2_col) && pcc2_col %in% names(row_data)) {
        m2_pcc <- suppressWarnings(as.numeric(row_data[[pcc2_col]]))
      }

      db_m2 <- NULL
      if (!is.na(m2_name) && !is.null(meme_db)) {
        db_m2 <- get_db_motif(meme_db, m2_name)
      }
      m2_title <- "DB: N/A"
      if (!is.na(m2_name)) {
        m2_title <- paste0("DB2: ", m2_name)
        if (!is.na(m2_pcc)) m2_title <- paste0(m2_title, " (PCC=", round(m2_pcc, 3), ")")
      }
      grob_db2 <- plot_logo_bits(db_m2, title = m2_title)

      grob_list <- c(grob_list, list(grob_pcre, grob_db1, grob_db2))

      # Collect for CSV (v2.2: use True_Rank from pre-computed column)
      selected_motifs_list[[length(selected_motifs_list) + 1]] <- data.frame(
        Genotype       = geno,
        Cluster        = cl_name,
        Display_Name   = cl_display,
        Rank           = row_data$True_Rank,
        pCRE_Kmer      = if (!is.null(kmer_col) && !is.na(row_data[[kmer_col]]))
                            row_data[[kmer_col]] else NA,
        Family_Label   = ifelse(!is.na(family_label), family_label, NA),
        DB_Match_1     = ifelse(!is.na(m1_name), m1_name, NA),
        PCC_1          = ifelse(!is.na(m1_pcc), m1_pcc, NA),
        DB_Match_2     = ifelse(!is.na(m2_name), m2_name, NA),
        PCC_2          = ifelse(!is.na(m2_pcc), m2_pcc, NA),
        Importance     = row_data[[imp_col]],
        Is_Target_Hit  = (!is.null(target_row_idx) && r == nrow(top_rows)),
        stringsAsFactors = FALSE
      )
    }

    # Arrange into grid
    if (length(grob_list) > 0) {
      n_logo_rows <- nrow(top_rows)
      header_grobs <- list(
        grid::textGrob("pCRE (Aligned Query)",
                       gp = grid::gpar(fontsize = 10, fontface = "bold")),
        grid::textGrob("DB Match 1 (Reference)",
                       gp = grid::gpar(fontsize = 10, fontface = "bold")),
        grid::textGrob("DB Match 2 (Reference)",
                       gp = grid::gpar(fontsize = 10, fontface = "bold"))
      )
      all_grobs   <- c(header_grobs, grob_list)
      row_heights <- c(0.6, rep(2.5, n_logo_rows))

      pdf(logo_path, width = 14, height = 2.5 * n_logo_rows + 1.5)
      grid.newpage()
      tryCatch({
        grid.arrange(
          grobs   = all_grobs, ncol = 3, nrow = n_logo_rows + 1,
          heights = unit(row_heights, "inches"),
          widths  = unit(c(1, 1, 1), "null"),
          top     = textGrob(
            sprintf("Motif Alignment (DB-Anchored) — %s (Top %d)", cl_display, top_n),
            gp = gpar(fontsize = 13, fontface = "bold")
          )
        )
      }, error = function(e) {
        grid.newpage()
        grid.text(sprintf("Layout error: %s", e$message), gp = gpar(col = "red"))
      })
      dev.off()
    }

    cat(sprintf("%d motifs plotted.\n", nrow(top_rows)))
  }

  # --- Save Selected_Motifs_List.csv ---
  if (length(selected_motifs_list) > 0) {
    csv_df <- bind_rows(selected_motifs_list)
    csv_path <- file.path(output_dir, "Selected_Motifs_List.csv")
    write.csv(csv_df, csv_path, row.names = FALSE)
    cat(sprintf("\n  Selected motifs CSV: %s\n", csv_path))
  }

  cat(sprintf("  Cluster reports saved to: %s\n", report_dir))
  return(invisible(NULL))
}

################################################################################
# MAIN
################################################################################

main <- function() {

  # ===== Step 1: Load all data =====
  data <- load_all_data(opt$p6b_dir, opt$p6c_dir, genotypes_to_process)

  # ===== Step 2: Figure 1 — AUC bar plot =====
  tryCatch(
    plot_auc_barplot(data$perf, data$integ, opt$output_dir,
                     wt_parsed = wt_label_parsed,
                     mut_parsed = mut_label_parsed,
                     geno_colors = geno_colors,
                     fig_width = opt$fig_width, fig_height = opt$fig_height),
    error = function(e) cat(sprintf("  ERROR in Fig 1: %s\n", e$message))
  )

  # ===== Step 3: Figure 2 — Similarity heatmap =====
  tryCatch(
    plot_similarity_heatmap(data$integ, opt$output_dir,
                            metric = opt$heatmap_metric,
                            top_n_features = opt$top_n_features,
                            wt_parsed = wt_label_parsed,
                            mut_parsed = mut_label_parsed,
                            geno_colors = geno_colors,
                            fig_width = opt$fig_width,
                            fig_height = opt$fig_height + 1),
    error = function(e) cat(sprintf("  ERROR in Fig 2: %s\n", e$message))
  )

  # ===== Step 4: Figure 3 — Per-cluster reports =====
  meme_db <- load_meme_db_safe(opt$meme_db)

  tryCatch(
    generate_cluster_reports(
      integ_df       = data$integ,
      pwm_lists      = data$pwms,
      meme_db        = meme_db,
      output_dir     = opt$output_dir,
      top_n          = opt$plot_top_n,
      min_overlap    = opt$min_overlap,
      target_keyword = opt$target_keyword,
      target_scope   = opt$target_scope,
      wt_parsed      = wt_label_parsed,
      mut_parsed     = mut_label_parsed,
      geno_colors    = geno_colors
    ),
    error = function(e) cat(sprintf("  ERROR in Fig 3: %s\n", e$message))
  )

  # ===== Summary =====
  cat("\n")
  cat("================================================================================\n")
  cat("  Pipeline 6d v2.2 Complete\n")
  cat("================================================================================\n\n")
  cat(sprintf("End time   : %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
  cat(sprintf("Output dir : %s\n\n", opt$output_dir))
  cat("Output structure:\n")
  cat("  +-- Fig1_AUC_ROC_Barplot.pdf/png\n")
  cat("  +-- Fig2_Similarity_Heatmap.pdf/png\n")
  cat("  +-- Selected_Motifs_List.csv\n")
  cat("  +-- Cluster_Reports/\n")
  cat("  |   +-- {Genotype}_{Cluster}/\n")
  cat("  |       +-- Scatter_Importance_Enrichment.pdf\n")
  cat("  |       +-- Motif_Alignments_TopN.pdf\n")
  cat("  |       +-- Top20_Motif_Details.tsv\n")
  cat("\nDone.\n")
}

main()