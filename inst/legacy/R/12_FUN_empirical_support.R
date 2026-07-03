###############################################################################
# 12_FUN_empirical_support.R
#
# scONCO v1.1 升級項目 A.4 — Empirical AI / Literature Support
#
# 動機:
#   - 原始 `ai_support` 與 `lit_support` 是主觀 0–1 分數，靠人工標註，難以重現
#   - 改為計算指標:
#       ai_support_emp(g, type)  = 該 marker 出現於幾個 AI 模型版本 / 總版本數
#       lit_support_emp(g, type) = 該 marker 在 reference atlas 上的 Wilcoxon AUC
#                                  rescaled to [0,1]
#
# 對應文獻:
#   - Pullin & McCarthy (2024) Genome Biology 25, 56. — AUC 為 marker selection 最佳指標
#   - mLLMCelltype, AnnDictionary — multi-LLM consensus 方法論
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

###############################################################################
## 1) Compute empirical AI support from multiple AI-generated marker DBs
###############################################################################

#' Compute empirical AI support by tallying across AI version DBs
#'
#' @param markers_df consensus DB with columns gene, cell_type
#' @param ai_dbs named list of data.frames, each with columns gene, cell_type
#'              names() are version IDs (e.g., "ChatGPT5_DR_v1", "Gemini_v2")
#' @return markers_df with new column ai_support_emp in [0,1]
compute_empirical_ai_support <- function(markers_df, ai_dbs) {
  stopifnot(all(c("gene", "cell_type") %in% colnames(markers_df)))
  stopifnot(length(ai_dbs) > 0)
  for (db in ai_dbs) {
    stopifnot(all(c("gene", "cell_type") %in% colnames(db)))
  }

  n_dbs <- length(ai_dbs)
  message(sprintf("[ai_support] Tallying support across %d AI DB versions...", n_dbs))

  # For each (gene, cell_type) row in markers_df, count in how many AI DBs it appears
  keys_consensus <- paste(markers_df$gene, markers_df$cell_type, sep = "@@")

  hit_counts <- rep(0L, nrow(markers_df))
  for (ver in names(ai_dbs)) {
    db <- ai_dbs[[ver]]
    keys_ver <- paste(db$gene, db$cell_type, sep = "@@")
    hit_counts <- hit_counts + as.integer(keys_consensus %in% keys_ver)
  }

  markers_df$ai_support_emp <- hit_counts / n_dbs
  markers_df$ai_support_n_hits <- hit_counts

  message(sprintf("[ai_support] Distribution of AI support:"))
  print(summary(markers_df$ai_support_emp))

  markers_df
}

###############################################################################
## 2) Compute empirical literature support from reference atlas AUC
##    Uses the AUC computation from 09_FUN_marker_specificity_AUC.R
###############################################################################

#' Compute empirical literature support via Wilcoxon AUC on a curated atlas
#'
#' This treats the published, expert-annotated atlas (e.g., Reynolds 2021)
#' as the "literature ground truth" and assigns higher lit_support to markers
#' that empirically discriminate in that atlas.
#'
#' @param markers_df with gene, cell_type
#' @param ref_expr genes × cells reference expression
#' @param ref_labels cell labels matching scONCO ontology
#' @param auc_floor minimum AUC to count as supported (default 0.7)
#' @return markers_df with column lit_support_emp in [0,1]
compute_empirical_lit_support <- function(markers_df, ref_expr, ref_labels,
                                          auc_floor = 0.7) {
  # Reuse compute_marker_AUC_batch if 09 has been sourced
  if (!exists("compute_marker_AUC_batch")) {
    stop("Please source 09_FUN_marker_specificity_AUC.R first")
  }

  message("[lit_support] Computing AUC on reference atlas...")
  tmp <- compute_marker_AUC_batch(ref_expr, ref_labels, markers_df)
  auc <- tmp$auc_empirical

  # Rescale: AUC < floor → 0; floor → 0.5; 1.0 → 1.0
  rescale <- function(x) {
    out <- rep(0, length(x))
    valid <- !is.na(x)
    above <- valid & x >= auc_floor
    # Linear: auc_floor → 0, 1.0 → 1.0
    out[above] <- (x[above] - auc_floor) / (1 - auc_floor)
    out
  }

  markers_df$lit_support_emp <- pmin(1, pmax(0, rescale(auc)))
  markers_df$auc_for_lit <- auc

  message(sprintf("[lit_support] %d/%d markers exceeded AUC floor of %.2f",
                  sum(markers_df$lit_support_emp > 0, na.rm = TRUE),
                  nrow(markers_df), auc_floor))

  markers_df
}

###############################################################################
## 3) Refresh marker weights using empirical AI/lit support
##    Drop-in replacement for sconco_compute_marker_weights() in 01_FUN_*
###############################################################################

#' Recompute scONCO marker weights with empirical (data-driven) AI/lit support
#'
#' @param markers_df with columns gene, cell_type, negative, optionally:
#'   - specificity_empirical (from 09; preferred over default freq-based)
#'   - ai_support_emp (from compute_empirical_ai_support)
#'   - lit_support_emp (from compute_empirical_lit_support)
#' @param alpha,beta,gamma weights for specificity, AI, lit
#' @return markers_df with `weight` column
refresh_marker_weights_empirical <- function(markers_df,
                                             alpha = 0.6,
                                             beta  = 0.2,
                                             gamma = 0.2) {
  stopifnot(all(c("gene", "cell_type") %in% colnames(markers_df)))

  # Specificity (prefer empirical; fall back to default freq-based)
  if ("specificity_empirical" %in% colnames(markers_df)) {
    spec <- markers_df$specificity_empirical
  } else if ("specificity" %in% colnames(markers_df)) {
    spec <- markers_df$specificity
    message("[refresh] Using default `specificity` (no `specificity_empirical` found)")
  } else {
    # Compute simple freq-based fallback
    n_types_total <- length(unique(markers_df$cell_type))
    freq <- markers_df %>%
      group_by(gene) %>%
      summarise(n_types = n_distinct(cell_type), .groups = "drop")
    spec_lookup <- 1 - freq$n_types / n_types_total
    names(spec_lookup) <- freq$gene
    spec <- spec_lookup[markers_df$gene]
    message("[refresh] Computed simple freq-based specificity")
  }
  spec[is.na(spec)] <- 0

  # AI support
  if ("ai_support_emp" %in% colnames(markers_df)) {
    ai <- markers_df$ai_support_emp
  } else {
    ai <- if ("ai_support" %in% colnames(markers_df)) markers_df$ai_support else 0
  }
  ai[is.na(ai)] <- 0

  # Lit support
  if ("lit_support_emp" %in% colnames(markers_df)) {
    lit <- markers_df$lit_support_emp
  } else {
    lit <- if ("lit_support" %in% colnames(markers_df)) markers_df$lit_support else 0
  }
  lit[is.na(lit)] <- 0

  # Composite weight
  markers_df$weight <- alpha * spec + beta * ai + gamma * lit

  # Optional negative-marker handling
  if ("negative" %in% colnames(markers_df)) {
    # Negative markers keep their weight but will be used as penalty downstream
    # (handled in sconco_activity_scores via the `δ` parameter)
  }

  message(sprintf("[refresh] Weights computed: mean=%.3f sd=%.3f min=%.3f max=%.3f",
                  mean(markers_df$weight, na.rm = TRUE),
                  sd(markers_df$weight, na.rm = TRUE),
                  min(markers_df$weight, na.rm = TRUE),
                  max(markers_df$weight, na.rm = TRUE)))

  markers_df
}

###############################################################################
## 4) Full v1.1 marker DB upgrade pipeline (one-shot)
###############################################################################

#' One-shot upgrade: empirical specificity + AI support + lit support + weights
#'
#' @param markers_df starting DB (e.g., DB_pancancer_human_v1.0)
#' @param ai_dbs named list of AI version DBs (for empirical AI support)
#' @param ref_expr reference atlas expression matrix
#' @param ref_labels reference atlas labels (in scONCO ontology)
#' @param weights_for_specificity weights for AUC/tau/mutex composition
#' @param alpha,beta,gamma weights for final marker score
#' @return enriched markers_df ready for scONCO inference
upgrade_to_scONCO_v1_1 <- function(markers_df,
                                    ai_dbs,
                                    ref_expr,
                                    ref_labels,
                                    weights_for_specificity = c(auc = 0.6, tau = 0.3, mutex = 0.1),
                                    alpha = 0.6, beta = 0.2, gamma = 0.2) {
  # Source companion modules if not already in env
  if (!exists("update_specificity_empirical")) {
    stop("Please source 09_FUN_marker_specificity_AUC.R first")
  }

  message("\n=== scONCO v1.1 marker DB upgrade ===")
  message("[Step 1/4] Empirical specificity (AUC + tau + mutex)")
  markers_df <- update_specificity_empirical(
    markers_df, ref_expr, ref_labels,
    weights = weights_for_specificity
  )

  message("\n[Step 2/4] Empirical AI support")
  markers_df <- compute_empirical_ai_support(markers_df, ai_dbs)

  message("\n[Step 3/4] Empirical literature support")
  markers_df <- compute_empirical_lit_support(markers_df, ref_expr, ref_labels)

  message("\n[Step 4/4] Final weight refresh")
  markers_df <- refresh_marker_weights_empirical(markers_df, alpha, beta, gamma)

  message("\n=== Upgrade complete ===")
  message(sprintf("  Markers: %d", nrow(markers_df)))
  message(sprintf("  Cell types covered: %d", length(unique(markers_df$cell_type))))
  message(sprintf("  Columns: %s", paste(colnames(markers_df), collapse = ", ")))

  markers_df
}

###############################################################################
## Example
###############################################################################
# # Load all AI version DBs as list of data.frames
# ai_dbs <- list(
#   DR_v1     = read.csv("database/versions/v_AI_ChatGPT5_DeepResearch/DRV1.csv"),
#   DR_v2     = read.csv("database/versions/v_AI_ChatGPT5_DeepResearch/DRV2.csv"),
#   Pro_v1    = read.csv("database/versions/v_AI_ChatGPT5_Pro/PV1.csv"),
#   Pro_v2    = read.csv("database/versions/v_AI_ChatGPT5_Pro/PV2.csv"),
#   Think_v1  = read.csv("database/versions/v_AI_ChatGPT5_Thinking/TV1.csv"),
#   Inst_v1   = read.csv("database/versions/v_AI_ChatGPT5_Instant/IV1.csv"),
#   Gemini_v1 = read.csv("database/versions/v_AI_Gemini/v1.csv"),
#   Gemini_v2 = read.csv("database/versions/v_AI_Gemini/v2.csv"),
#   Gemini_v3 = read.csv("database/versions/v_AI_Gemini/v3.csv")
# )
#
# # Reference atlas (Reynolds 2021)
# seu_ref <- readRDS("path/to/reynolds_2021.rds")
# ref_expr   <- GetAssayData(seu_ref, layer = "data")
# ref_labels <- seu_ref$cell_type_level3
#
# # One-shot upgrade
# markers_df_v1_1 <- upgrade_to_scONCO_v1_1(
#   markers_df = markers_df,
#   ai_dbs = ai_dbs,
#   ref_expr = ref_expr,
#   ref_labels = ref_labels
# )
#
# saveRDS(markers_df_v1_1, "database/current/markers_df_v1_1.rds")
