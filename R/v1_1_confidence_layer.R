# =============================================================
# v1.1 statistical confidence layer — roxygen documentation
# (Implementations live in code/R/09-12_FUN_*.R and are bridged via
#  zzz_source_legacy.R until devtools::document() is run.)
# =============================================================

#' Compute marker AUC specificity (A.1)
#'
#' Wilcoxon AUC of each marker against all other cell types in the reference
#' atlas. Higher AUC means more specific to the target type.
#'
#' @param expr Numeric matrix (genes x cells) of normalised expression.
#' @param labels Character vector of cell-type labels matching the columns of
#'   `expr`.
#' @param genes Optional subset of genes to test (default: all rows of `expr`).
#'
#' @return A long-format data frame with columns `gene`, `celltype`, `auc`.
#'
#' @references Pullin, J.M. & McCarthy, D.J. (2024). A comparison of marker
#'   gene selection methods for scRNA-seq. *Genome Biology* 25, 56.
#'
#' @export
compute_marker_AUC_batch <- function(expr, labels, genes = NULL) {
  stop("Implementation provided by code/R/09_FUN_marker_specificity_AUC.R.")
}

#' Compute tau score (A.1 helper)
#'
#' Tau = how concentrated a marker's expression is in a single type.
#'
#' @inheritParams compute_marker_AUC_batch
#' @export
compute_tau_score <- function(expr, labels, genes = NULL) {
  stop("Implementation provided by code/R/09_FUN_marker_specificity_AUC.R.")
}

#' Compute mutual-exclusivity score (A.1 helper)
#'
#' Mutex score = how clearly a marker is on in target vs off in non-target.
#' @inheritParams compute_marker_AUC_batch
#' @references Wang, F. et al. (2020). SCMarker. *PLOS Comp Biol* 16(3), e1007445.
#' @export
compute_mutex_score <- function(expr, labels, genes = NULL) {
  stop("Implementation provided by code/R/09_FUN_marker_specificity_AUC.R.")
}

#' Update specificity column empirically (A.1)
#'
#' Recomputes `markers_df$specificity_index` using AUC on a reference atlas.
#'
#' @param markers_df Marker data frame.
#' @param expr Reference expression matrix.
#' @param labels Reference cell-type labels.
#' @return Updated `markers_df`.
#' @export
update_specificity_empirical <- function(markers_df, expr, labels) {
  stop("Implementation provided by code/R/09_FUN_marker_specificity_AUC.R.")
}

# ----------------------- A.2: hierarchical reject ---------------------

#' Per-cell softmax confidence (A.2 helper)
#' @param score_matrix Numeric matrix (cells x types) of unnormalised scores.
#' @return Numeric matrix of softmax probabilities.
#' @export
compute_confidence <- function(score_matrix) {
  stop("Implementation provided by code/R/10_FUN_hierarchical_reject.R.")
}

#' Apply reject option to one cell (A.2 helper)
#' @export
hierarchical_reject_one <- function(probs, threshold) {
  stop("Implementation provided by code/R/10_FUN_hierarchical_reject.R.")
}

#' Apply hierarchical reject to a Seurat object (A.2)
#'
#' At each level k, defers cells with confidence < threshold_k to the parent
#' level. Adds `scONCO_L{k}_confidence` and `scONCO_L{k}_status` columns.
#'
#' @param seurat_obj A Seurat object previously annotated by [run_scONCO()].
#' @param thresholds Named numeric vector, e.g.
#'   `c(L1 = 0.7, L2 = 0.6, L3 = 0.55, L4 = 0.5)`.
#' @return Updated Seurat object.
#' @references de Lichtenberg, U.E. et al. (2024). Uncertainty-aware annotation
#'   with hierarchical reject. *Bioinformatics* 40(3), btae128.
#' @export
apply_hierarchical_reject <- function(seurat_obj,
                                      thresholds = c(L1 = 0.7, L2 = 0.6,
                                                     L3 = 0.55, L4 = 0.5)) {
  stop("Implementation provided by code/R/10_FUN_hierarchical_reject.R.")
}

#' Calibrate reject thresholds (A.2)
#'
#' Choose threshold_k for each level by sweeping a grid on a held-out
#' reference set to achieve a target reject rate or accuracy floor.
#' @export
calibrate_reject_thresholds <- function(seurat_obj, target = c("reject", "accuracy"),
                                         target_value = 0.1) {
  stop("Implementation provided by code/R/10_FUN_hierarchical_reject.R.")
}

# ----------------------- A.3: split conformal -------------------------

#' Softmax over a score matrix (helper)
#' @param x Numeric matrix.
#' @return Same-shape matrix of softmax probabilities.
#' @export
softmax_matrix <- function(x) {
  stop("Implementation provided by code/R/11_FUN_conformal_prediction.R.")
}

#' Calibrate a split conformal predictor (A.3)
#'
#' @param prob_calib Numeric matrix (cells x types) of softmax probabilities
#'   on the calibration set.
#' @param y_calib Integer vector of true label indices.
#' @param alpha Target miscoverage rate (default 0.1).
#' @return List with `q_hat` (nonconformity quantile), `alpha`, `n_calib`.
#' @references Boudaoud, R. et al. (2025). Conformal inference for reliable
#'   scRNA-seq annotation. *Bioinformatics* 41(10), btaf521.
#' @export
split_conformal_calibrate <- function(prob_calib, y_calib, alpha = 0.1) {
  stop("Implementation provided by code/R/11_FUN_conformal_prediction.R.")
}

#' Produce conformal prediction set (A.3)
#'
#' @param prob_test Numeric matrix (cells x types) of softmax probabilities.
#' @param calib Output of `split_conformal_calibrate()`.
#' @return List with `set` (list of integer vectors), `set_size` (integer).
#' @export
conformal_predict <- function(prob_test, calib) {
  stop("Implementation provided by code/R/11_FUN_conformal_prediction.R.")
}

#' Apply conformal prediction to a Seurat object (A.3)
#'
#' Splits cells into calibration and test subsets, computes conformal sets,
#' writes `conformal_set` and `conformal_set_size` to meta.data.
#'
#' @param seurat_obj A Seurat object previously annotated by [run_scONCO()].
#' @param alpha Target miscoverage rate.
#' @param calibration_split Fraction of cells to hold out for calibration.
#' @return Updated Seurat object.
#' @export
apply_conformal_to_seurat <- function(seurat_obj, alpha = 0.1,
                                      calibration_split = 0.2) {
  stop("Implementation provided by code/R/11_FUN_conformal_prediction.R.")
}

# ----------------------- A.4: empirical support -----------------------

#' Compute empirical AI support (A.4)
#'
#' Counts the fraction of AI source DBs that include each marker for each type.
#' @param markers_long_df Long-format DB with columns gene, celltype, source.
#' @return Data frame with `gene`, `celltype`, `ai_support`.
#' @export
compute_empirical_ai_support <- function(markers_long_df) {
  stop("Implementation provided by code/R/12_FUN_empirical_support.R.")
}

#' Compute empirical literature support (A.4)
#'
#' Queries PubMed (via easyPubMed) for co-occurrence of marker + cell type;
#' log-normalises to a 0–1 score.
#' @export
compute_empirical_lit_support <- function(genes, celltypes,
                                          max_hits = 1000, sleep = 0.34) {
  stop("Implementation provided by code/R/12_FUN_empirical_support.R.")
}

#' Refresh marker weights with empirical support (A.4)
#'
#' Replaces subjective ai_support / lit_support columns with empirically
#' derived ones.
#' @export
refresh_marker_weights_empirical <- function(markers_df,
                                              ai_source_dbs = NULL,
                                              lit_support_path = NULL) {
  stop("Implementation provided by code/R/12_FUN_empirical_support.R.")
}

#' Upgrade an existing scONCO pipeline to v1.1
#'
#' One-call convenience wrapper that applies A.1–A.4 in order.
#' @param markers_df DB to upgrade.
#' @param ref_expr Reference expression matrix for AUC computation.
#' @param ref_labels Reference cell-type labels.
#' @return Upgraded `markers_df` with `specificity_index`, `ai_support`,
#'   `lit_support` columns populated empirically.
#' @export
upgrade_to_scONCO_v1_1 <- function(markers_df, ref_expr = NULL,
                                     ref_labels = NULL) {
  stop("Implementation provided by code/R/12_FUN_empirical_support.R.")
}
