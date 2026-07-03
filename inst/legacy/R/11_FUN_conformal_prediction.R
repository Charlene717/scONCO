###############################################################################
# 11_FUN_conformal_prediction.R
#
# scONCO v1.1 升級項目 A.3 — Split Conformal Prediction
#
# 對應文獻:
#   - Boudaoud, R. et al. (2025) Bioinformatics 41(10), btaf521.
#     "Conformal inference for reliable single cell RNA-seq annotation"
#   - Banerjee, P. et al. (2024) Springer CIBB.
#     "Conformal Inference for Cell Type Prediction with Graph-Structured Constraints"
#
# 核心思想:
#   給定校準集 (calibration set) 與目標 coverage 1-α (e.g., 90%):
#   1. 為每個校準細胞計算 "non-conformity score" s_i = 1 - p(true_class_i | x_i)
#   2. 取 (1-α)(n+1)/n 分位數作為閾值 q̂
#   3. 對測試細胞 x_test 輸出 prediction set:
#        {y : p(y | x_test) ≥ 1 - q̂}
#   4. 保證 P(true_label ∈ prediction_set) ≥ 1 - α (邊際 coverage)
#
# Prediction set 大小 = uncertainty 量化:
#   - size = 1 → 高信心
#   - size > 1 → ambiguous (列出所有可能類型)
#   - size = 0 或包含 "Unknown" → out-of-distribution
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

###############################################################################
## 1) Compute softmax probability matrix from raw scores
##    Mirrors compute_confidence() but keeps full matrix
###############################################################################

#' Convert raw activity scores to softmax probabilities
#'
#' @param scores_mat cells × cell_types numeric matrix
#' @param temperature softmax temperature (>1 flatter, <1 sharper). Default 1
#' @return same shape, rows sum to 1
softmax_matrix <- function(scores_mat, temperature = 1) {
  z <- scores_mat / temperature
  z <- z - apply(z, 1, max)
  expz <- exp(z)
  expz / rowSums(expz)
}

###############################################################################
## 2) Split conformal calibration
##
## Input:  calibration matrix + true labels
## Output: threshold q̂ for given confidence level
###############################################################################

#' Calibrate conformal threshold from a held-out calibration set
#'
#' @param cal_probs n_cal × n_types — softmax probabilities for calibration cells
#' @param cal_labels character vector of true cell type labels (length n_cal)
#' @param alpha mis-coverage rate (default 0.1 → 90% coverage)
#' @return list with:
#'   q_hat: scalar threshold
#'   scores: non-conformity score per calibration cell
split_conformal_calibrate <- function(cal_probs, cal_labels, alpha = 0.1) {
  stopifnot(nrow(cal_probs) == length(cal_labels))
  stopifnot(all(cal_labels %in% colnames(cal_probs)))

  # Non-conformity score: 1 - p(true_class)
  scores <- sapply(seq_along(cal_labels), function(i) {
    1 - cal_probs[i, cal_labels[i]]
  })

  # (1-α)(n+1)/n quantile, capped at 1
  n <- length(scores)
  q_level <- min(1, (1 - alpha) * (n + 1) / n)
  q_hat <- quantile(scores, probs = q_level, type = 6, na.rm = TRUE)

  message(sprintf("[Conformal] Calibration n=%d  α=%.2f  q̂=%.4f",
                  n, alpha, q_hat))

  list(q_hat = unname(q_hat), scores = scores, alpha = alpha, n_cal = n)
}

###############################################################################
## 3) Conformal prediction sets for test cells
##
## A type y is included in the prediction set for cell i iff
##   1 - p(y | x_i) ≤ q̂
##   ⇔ p(y | x_i) ≥ 1 - q̂
###############################################################################

#' Generate prediction sets for test cells
#'
#' @param test_probs n_test × n_types softmax probabilities
#' @param q_hat threshold from split_conformal_calibrate()
#' @return list of length n_test; each element = character vector of types in set
conformal_predict <- function(test_probs, q_hat) {
  if (is.null(colnames(test_probs))) {
    stop("test_probs must have column names (cell type names)")
  }
  threshold <- 1 - q_hat
  apply(test_probs, 1, function(p) {
    sel <- which(p >= threshold)
    if (length(sel) == 0) {
      # No type passes — return top-1 (degenerate case)
      top <- which.max(p)
      return(c(colnames(test_probs)[top], "(forced)"))
    }
    colnames(test_probs)[sel]
  }, simplify = FALSE)
}

###############################################################################
## 4) Summarize prediction sets
###############################################################################

#' Summarize prediction set statistics
#'
#' @param pred_sets list from conformal_predict()
#' @return data.frame with cell_id, set_size, top_label, ambiguity_flag
summarize_prediction_sets <- function(pred_sets, cell_ids = NULL) {
  n <- length(pred_sets)
  if (is.null(cell_ids)) cell_ids <- seq_len(n)

  sizes <- sapply(pred_sets, length)
  top_labels <- sapply(pred_sets, function(s) {
    s <- setdiff(s, "(forced)")
    if (length(s) == 0) "Unknown" else s[1]
  })
  ambiguous <- sizes > 1
  forced <- sapply(pred_sets, function(s) "(forced)" %in% s)

  df <- data.frame(
    cell_id = cell_ids,
    set_size = sizes,
    top_label = top_labels,
    is_ambiguous = ambiguous,
    is_oodlike = forced,
    pred_set_str = sapply(pred_sets, paste, collapse = "|"),
    stringsAsFactors = FALSE
  )

  message(sprintf("[Conformal] Summary on n=%d test cells:", n))
  message(sprintf("  Mean set size:        %.2f", mean(sizes)))
  message(sprintf("  Median set size:      %d", median(sizes)))
  message(sprintf("  Singleton (size=1):   %.1f%%", 100 * mean(sizes == 1)))
  message(sprintf("  Ambiguous (size>1):   %.1f%%", 100 * mean(ambiguous)))
  message(sprintf("  OOD-like (forced):    %.1f%%", 100 * mean(forced)))

  df
}

###############################################################################
## 5) Empirical coverage validation
##    Sanity check: on a held-out test set with known labels,
##    actual coverage should match 1 - α
###############################################################################

#' Validate empirical coverage
#'
#' @param pred_sets list from conformal_predict
#' @param true_labels character vector of true labels
#' @return list(coverage, expected, gap)
evaluate_coverage <- function(pred_sets, true_labels, alpha = 0.1) {
  stopifnot(length(pred_sets) == length(true_labels))
  hit <- sapply(seq_along(pred_sets), function(i) {
    true_labels[i] %in% pred_sets[[i]]
  })
  emp <- mean(hit)
  expected <- 1 - alpha
  message(sprintf("[Conformal] Empirical coverage: %.3f (target %.3f, gap %+.3f)",
                  emp, expected, emp - expected))
  list(coverage = emp, expected = expected, gap = emp - expected)
}

###############################################################################
## 6) Apply conformal prediction to Seurat object
###############################################################################

#' Apply conformal prediction to a Seurat object given a calibration sample
#'
#' Workflow:
#'   1) Use a subset of cells with expert/ground-truth labels as calibration
#'   2) Generate prediction sets for remaining (test) cells
#'   3) Store prediction set + size + ambiguity flag in @meta.data
#'
#' @param seurat_obj Seurat object
#' @param scores_matrix cells × cell_types raw activity scores
#'                      (rownames should match colnames of Seurat object)
#' @param calibration_cells character vector of cell barcodes in cal set
#' @param ground_truth_col meta column with ground truth labels for cal cells
#' @param alpha mis-coverage rate
#' @param prefix output column prefix (default "scONCO")
apply_conformal_to_seurat <- function(seurat_obj,
                                      scores_matrix,
                                      calibration_cells,
                                      ground_truth_col,
                                      alpha = 0.1,
                                      prefix = "scONCO") {
  stopifnot(all(calibration_cells %in% rownames(scores_matrix)))
  stopifnot(ground_truth_col %in% colnames(seurat_obj@meta.data))

  probs <- softmax_matrix(scores_matrix)

  cal_idx <- which(rownames(probs) %in% calibration_cells)
  test_idx <- setdiff(seq_len(nrow(probs)), cal_idx)

  cal_probs <- probs[cal_idx, , drop = FALSE]
  cal_labels <- seurat_obj@meta.data[calibration_cells, ground_truth_col]
  # Filter cal cells with labels we cover
  keep <- !is.na(cal_labels) & cal_labels %in% colnames(probs)
  cal_probs <- cal_probs[keep, , drop = FALSE]
  cal_labels <- cal_labels[keep]

  calib <- split_conformal_calibrate(cal_probs, cal_labels, alpha = alpha)

  test_probs <- probs[test_idx, , drop = FALSE]
  pred_sets <- conformal_predict(test_probs, calib$q_hat)
  smry <- summarize_prediction_sets(pred_sets, cell_ids = rownames(test_probs))

  # Write back to Seurat meta (cal cells get NA)
  n <- ncol(seurat_obj)
  meta_set_size <- rep(NA_integer_, n)
  meta_top      <- rep(NA_character_, n)
  meta_amb      <- rep(NA, n)
  meta_set_str  <- rep(NA_character_, n)

  test_barcodes <- rownames(test_probs)
  pos <- match(test_barcodes, colnames(seurat_obj))
  meta_set_size[pos] <- smry$set_size
  meta_top[pos]      <- smry$top_label
  meta_amb[pos]      <- smry$is_ambiguous
  meta_set_str[pos]  <- smry$pred_set_str

  seurat_obj@meta.data[[paste0(prefix, "_pset_size")]]    <- meta_set_size
  seurat_obj@meta.data[[paste0(prefix, "_pset_top")]]     <- meta_top
  seurat_obj@meta.data[[paste0(prefix, "_pset_ambig")]]   <- meta_amb
  seurat_obj@meta.data[[paste0(prefix, "_pset_str")]]     <- meta_set_str
  seurat_obj@misc[[paste0(prefix, "_conformal_calib")]] <- calib

  seurat_obj
}

###############################################################################
## Example
###############################################################################
# # After running scONCO to obtain Level4_Abb scores matrix:
# # `scores_l4` should be cells × cell_types
# probs_l4 <- softmax_matrix(scores_l4)
#
# # Split your labeled data ~50/50 into calibration vs test
# set.seed(42)
# all_cells <- rownames(scores_l4)
# cal <- sample(all_cells, size = floor(length(all_cells) / 2))
#
# seuratObject_Sample <- apply_conformal_to_seurat(
#   seurat_obj = seuratObject_Sample,
#   scores_matrix = scores_l4,
#   calibration_cells = cal,
#   ground_truth_col = "expert_annotation",
#   alpha = 0.1
# )
#
# # Sanity check coverage
# # ... evaluate_coverage(pred_sets, true_labels, alpha = 0.1)
