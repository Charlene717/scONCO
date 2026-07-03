###############################################################################
# 09_FUN_marker_specificity_AUC.R
#
# scONCO v1.1 升級項目 A.1 — 資訊論度量取代簡單 specificity
#
# 對應文獻:
#   - Pullin & McCarthy (2024) Genome Biology 25, 56.
#     "A comparison of marker gene selection methods..." — 59 方法 benchmark
#     發現 Wilcoxon AUC 為最佳單一指標
#   - Wang et al. (2020) PLOS Computational Biology — SCMarker
#     bimodal distribution + mutual exclusivity
#   - Dumitrascu et al. (2021) Nature Communications 12, 1186. — scGeneFit
#     linear programming for marker selection
#
# 取代 01_FUN_scONCO_core.R 內的 specificity = 1 - freq 設計
###############################################################################

suppressPackageStartupMessages({
  library(matrixStats)
  library(dplyr)
})

###############################################################################
## 1) Wilcoxon AUC per gene per target cell type
##    AUC ≡ concordance probability that a random cell of target type
##    expresses g higher than a random cell of any other type
##    Range [0, 1]; 0.5 = no signal; 1 = perfect marker
###############################################################################

#' Compute Wilcoxon AUC for one gene vs target cell type
#'
#' @param expr_vec numeric vector of expression (one gene across all cells)
#' @param cell_labels character vector of cell type labels
#' @param target_type the cell type to score against
#' @return AUC in [0, 1]; 0.5 means no separation
compute_marker_AUC_one <- function(expr_vec, cell_labels, target_type) {
  stopifnot(length(expr_vec) == length(cell_labels))
  in_target <- cell_labels == target_type
  n_pos <- sum(in_target)
  n_neg <- sum(!in_target)
  if (n_pos == 0 || n_neg == 0) return(NA_real_)
  # Wilcoxon U statistic via rank sum
  r <- rank(expr_vec, ties.method = "average")
  R_pos <- sum(r[in_target])
  U <- R_pos - n_pos * (n_pos + 1) / 2
  AUC <- U / (n_pos * n_neg)
  AUC
}

#' Batch Wilcoxon AUC for multiple genes × multiple cell types
#'
#' @param expr_mat genes (rows) × cells (cols) — Seurat data slot recommended
#' @param cell_labels character vector of length ncol(expr_mat)
#' @param markers_df data.frame with columns "gene" and "cell_type"
#'                   (one row per marker–cell type assignment)
#' @return markers_df with new column "auc_empirical"
compute_marker_AUC_batch <- function(expr_mat, cell_labels, markers_df) {
  stopifnot(all(c("gene", "cell_type") %in% colnames(markers_df)))
  stopifnot(ncol(expr_mat) == length(cell_labels))

  genes_avail <- intersect(unique(markers_df$gene), rownames(expr_mat))
  message(sprintf("[AUC] %d/%d markers available in expression matrix",
                  length(genes_avail), length(unique(markers_df$gene))))

  types_avail <- intersect(unique(markers_df$cell_type), unique(cell_labels))
  message(sprintf("[AUC] %d/%d cell types match labels",
                  length(types_avail), length(unique(markers_df$cell_type))))

  # Precompute rank matrix for efficiency (rank once per gene)
  expr_sub <- expr_mat[genes_avail, , drop = FALSE]
  if (inherits(expr_sub, "dgCMatrix")) expr_sub <- as.matrix(expr_sub)
  rank_mat <- t(matrixStats::colRanks(t(expr_sub), ties.method = "average"))
  rownames(rank_mat) <- rownames(expr_sub)

  markers_df$auc_empirical <- NA_real_
  for (i in seq_len(nrow(markers_df))) {
    g <- markers_df$gene[i]
    ct <- markers_df$cell_type[i]
    if (!(g %in% genes_avail) || !(ct %in% types_avail)) next

    in_target <- cell_labels == ct
    n_pos <- sum(in_target)
    n_neg <- sum(!in_target)
    if (n_pos == 0 || n_neg == 0) next

    r <- rank_mat[g, ]
    U <- sum(r[in_target]) - n_pos * (n_pos + 1) / 2
    markers_df$auc_empirical[i] <- U / (n_pos * n_neg)
  }

  markers_df
}

###############################################################################
## 2) Mutual information specificity (alternative to AUC)
##    Measures information gain between gene expression (binarized) and label
###############################################################################

#' Mutual information between binarized expression and cell type label
#'
#' @param expr_vec numeric vector of expression
#' @param cell_labels character/factor labels
#' @param threshold expression threshold (default = median > 0)
#' @return MI in nats (natural log)
compute_marker_MI_one <- function(expr_vec, cell_labels, threshold = NULL) {
  if (is.null(threshold)) {
    nonzero <- expr_vec[expr_vec > 0]
    threshold <- if (length(nonzero) > 0) median(nonzero) else 0
  }
  x_bin <- as.integer(expr_vec > threshold)
  y <- as.factor(cell_labels)

  joint <- table(x_bin, y) / length(x_bin)
  px <- rowSums(joint)
  py <- colSums(joint)

  mi <- 0
  for (i in seq_along(px)) {
    for (j in seq_along(py)) {
      pij <- joint[i, j]
      if (pij > 0) mi <- mi + pij * log(pij / (px[i] * py[j]))
    }
  }
  mi
}

###############################################################################
## 3) Tau-score (TS) — gene specificity across N cell types
##    TS = Σ_i (1 - x_i / x_max) / (N - 1)
##    Range [0, 1]; 1 = expressed in exactly one type
##    Reference: Yanai et al. (2005) Bioinformatics
###############################################################################

#' Compute tau-score specificity for each gene
#'
#' @param expr_mat genes × cells
#' @param cell_labels labels (length ncol)
#' @return named numeric vector (one tau per gene)
compute_tau_score <- function(expr_mat, cell_labels) {
  # Pseudobulk: mean expression per cell type per gene
  types <- unique(cell_labels)
  pseudobulk <- sapply(types, function(t) {
    idx <- which(cell_labels == t)
    if (length(idx) == 0) return(rep(NA, nrow(expr_mat)))
    if (inherits(expr_mat, "dgCMatrix")) {
      Matrix::rowMeans(expr_mat[, idx, drop = FALSE])
    } else {
      rowMeans(expr_mat[, idx, drop = FALSE])
    }
  })
  rownames(pseudobulk) <- rownames(expr_mat)
  colnames(pseudobulk) <- types

  # Tau per gene
  N <- ncol(pseudobulk)
  if (N < 2) stop("Need at least 2 cell types")

  tau <- apply(pseudobulk, 1, function(x) {
    x_max <- max(x, na.rm = TRUE)
    if (x_max <= 0 || is.na(x_max)) return(NA)
    sum(1 - x / x_max, na.rm = TRUE) / (N - 1)
  })
  tau
}

###############################################################################
## 4) Mutual exclusivity score (SCMarker-style)
##    For each marker g of type T, compute:
##      mutex(g, T) = fraction of OTHER markers of T whose expression is
##                    NOT co-expressed with g in non-T cells
##    Range [0, 1]; higher = more exclusive co-expression
###############################################################################

#' SCMarker-style mutual exclusivity score
#'
#' @param expr_mat genes × cells (counts or log-normalized)
#' @param cell_labels labels
#' @param markers_df with gene, cell_type
#' @return markers_df with new column mutex_score
compute_mutex_score <- function(expr_mat, cell_labels, markers_df) {
  stopifnot(all(c("gene", "cell_type") %in% colnames(markers_df)))

  # Binarize expression at non-zero
  bin_mat <- (expr_mat > 0) * 1
  if (inherits(bin_mat, "dgCMatrix")) bin_mat <- as.matrix(bin_mat)

  markers_df$mutex_score <- NA_real_
  for (ct in unique(markers_df$cell_type)) {
    ct_markers <- markers_df$gene[markers_df$cell_type == ct]
    ct_markers <- intersect(ct_markers, rownames(bin_mat))
    if (length(ct_markers) < 2) next

    out_cells <- which(cell_labels != ct)
    if (length(out_cells) == 0) next

    # Co-expression matrix in out-of-type cells
    sub_mat <- bin_mat[ct_markers, out_cells, drop = FALSE]
    # P(g_i AND g_j) in non-target cells
    coexpr <- (sub_mat %*% t(sub_mat)) / length(out_cells)
    diag(coexpr) <- NA

    # For each marker, mutex = 1 - mean P(coexpr with other markers)
    mutex_per_marker <- 1 - rowMeans(coexpr, na.rm = TRUE)
    idx <- which(markers_df$cell_type == ct & markers_df$gene %in% ct_markers)
    for (i in idx) {
      g <- markers_df$gene[i]
      if (g %in% names(mutex_per_marker)) {
        markers_df$mutex_score[i] <- mutex_per_marker[g]
      }
    }
  }
  markers_df
}

###############################################################################
## 5) Unified specificity update — combines AUC + tau + mutex
##    Replaces the simple `specificity = 1 - freq` in 01_FUN_scONCO_core.R
###############################################################################

#' Update markers_df with empirical specificity from reference dataset
#'
#' @param markers_df with gene, cell_type
#' @param expr_mat reference expression matrix (genes × cells)
#' @param cell_labels reference cell labels
#' @param weights how to combine AUC, tau, mutex (default 0.6/0.3/0.1)
#' @return markers_df with new columns: auc_empirical, tau_score, mutex_score,
#'                                       specificity_empirical
update_specificity_empirical <- function(markers_df,
                                         expr_mat,
                                         cell_labels,
                                         weights = c(auc = 0.6, tau = 0.3, mutex = 0.1)) {
  stopifnot(abs(sum(weights) - 1) < 1e-6)

  message("[1/3] Computing Wilcoxon AUC per marker-celltype pair...")
  markers_df <- compute_marker_AUC_batch(expr_mat, cell_labels, markers_df)

  message("[2/3] Computing tau-score per gene...")
  tau <- compute_tau_score(expr_mat, cell_labels)
  markers_df$tau_score <- tau[markers_df$gene]

  message("[3/3] Computing mutex score per marker-celltype pair...")
  markers_df <- compute_mutex_score(expr_mat, cell_labels, markers_df)

  # Composite: AUC dominates; missing values get 0.5 (no info)
  auc_term <- ifelse(is.na(markers_df$auc_empirical), 0.5, markers_df$auc_empirical)
  # AUC > 0.5 is informative; rescale so 0.5 → 0, 1 → 1
  auc_score <- pmax(0, 2 * (auc_term - 0.5))

  tau_term   <- ifelse(is.na(markers_df$tau_score), 0, markers_df$tau_score)
  mutex_term <- ifelse(is.na(markers_df$mutex_score), 0, markers_df$mutex_score)

  markers_df$specificity_empirical <-
    weights["auc"]   * auc_score +
    weights["tau"]   * tau_term +
    weights["mutex"] * mutex_term

  markers_df
}

###############################################################################
## Example usage
###############################################################################
# # Assuming you have a Seurat object with curated cell type labels
# library(Seurat)
# seu_reynolds <- readRDS("path/to/reynolds_2021_annotated.rds")
# expr_mat <- GetAssayData(seu_reynolds, assay = "RNA", layer = "data")
# cell_labels <- seu_reynolds$cell_type_level3
#
# # Update markers_df (must already exist from DB load)
# markers_df <- update_specificity_empirical(markers_df, expr_mat, cell_labels)
#
# # Then replace `specificity` with `specificity_empirical` in
# # sconco_compute_marker_weights() — or pass via custom alpha-weighted formula
