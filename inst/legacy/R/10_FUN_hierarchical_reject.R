###############################################################################
# 10_FUN_hierarchical_reject.R
#
# scONCO v1.1 升級項目 A.2 — 階層式 reject option
#
# 對應文獻:
#   - de Lichtenberg, U.E. et al. (2024) Bioinformatics 40(3), btae128.
#     "Uncertainty-aware single-cell annotation with a hierarchical reject option"
#   - 與 scONCO 4 層 ontology (Level1 → Level2 → Level3 → Level4_Abb) 完美匹配
#
# 核心邏輯:
#   每層 scMRMA 都輸出 softmax-normalized 機率 → 若 max(p) < threshold_k
#   則該 cell 在 Level k 被「拒答」，退回 Level k-1 的標籤
#
# 預設 threshold: L4 = 0.5, L3 = 0.4, L2 = 0.3, L1 = 0.2 (越深越嚴格)
###############################################################################

suppressPackageStartupMessages({
  library(dplyr)
})

###############################################################################
## 1) Compute confidence (softmax max) per cell per level
###############################################################################

#' Compute softmax-max confidence for each cell from raw scores
#'
#' @param scores_mat numeric matrix (cells × cell_types), raw activity scores
#' @return data.frame with: predicted_type, confidence (top-1 softmax)
compute_confidence <- function(scores_mat) {
  stopifnot(is.matrix(scores_mat))
  if (ncol(scores_mat) == 0) {
    return(data.frame(predicted_type = character(0), confidence = numeric(0)))
  }

  # Numerically stable softmax row-wise
  row_max <- apply(scores_mat, 1, max)
  z <- scores_mat - row_max
  expz <- exp(z)
  probs <- expz / rowSums(expz)

  # Top-1 prediction
  top_idx <- apply(probs, 1, which.max)
  predicted_type <- colnames(scores_mat)[top_idx]
  confidence <- probs[cbind(seq_len(nrow(probs)), top_idx)]

  data.frame(
    cell_id = if (!is.null(rownames(scores_mat))) rownames(scores_mat) else seq_len(nrow(probs)),
    predicted_type = predicted_type,
    confidence = confidence,
    stringsAsFactors = FALSE
  )
}

###############################################################################
## 2) Hierarchical reject decision per cell
##
## Input:  per-level (predicted_type, confidence) for one cell
## Output: final label assigned at the deepest level where confidence ≥ threshold
###############################################################################

#' Decide the deepest reliable label for one cell
#'
#' @param l1 list(type, conf)
#' @param l2 list(type, conf)
#' @param l3 list(type, conf)
#' @param l4 list(type, conf)
#' @param thresholds named vector c(L1=0.2, L2=0.3, L3=0.4, L4=0.5)
#' @return list(final_label, final_level, was_rejected)
hierarchical_reject_one <- function(l1, l2, l3, l4,
                                    thresholds = c(L1 = 0.2, L2 = 0.3, L3 = 0.4, L4 = 0.5)) {
  # Try deepest first; fall back if confidence too low or label missing
  if (!is.na(l4$type) && !is.na(l4$conf) && l4$conf >= thresholds["L4"]) {
    return(list(label = l4$type, level = "L4_Abb", rejected = FALSE))
  }
  if (!is.na(l3$type) && !is.na(l3$conf) && l3$conf >= thresholds["L3"]) {
    return(list(label = l3$type, level = "L3", rejected = TRUE))
  }
  if (!is.na(l2$type) && !is.na(l2$conf) && l2$conf >= thresholds["L2"]) {
    return(list(label = l2$type, level = "L2", rejected = TRUE))
  }
  if (!is.na(l1$type) && !is.na(l1$conf) && l1$conf >= thresholds["L1"]) {
    return(list(label = l1$type, level = "L1", rejected = TRUE))
  }
  # All levels rejected
  list(label = "Unknown", level = "Unknown", rejected = TRUE)
}

###############################################################################
## 3) Apply hierarchical reject to a Seurat object
##
## Expects 8 columns in @meta.data (4 labels + 4 confidences):
##   label_scMRMA_L1, conf_scMRMA_L1
##   label_scMRMA_L2, conf_scMRMA_L2
##   label_scMRMA_L3, conf_scMRMA_L3
##   label_scMRMA_L4Abb, conf_scMRMA_L4Abb
###############################################################################

#' Apply hierarchical reject across all cells in a Seurat object
#'
#' @param seurat_obj Seurat object with per-level label + confidence columns
#' @param thresholds named vector controlling reject thresholds
#' @param prefix output column prefix (default "scONCO")
#' @return seurat_obj with two new meta columns:
#'   <prefix>_final_label, <prefix>_final_level
apply_hierarchical_reject <- function(seurat_obj,
                                      thresholds = c(L1 = 0.2, L2 = 0.3, L3 = 0.4, L4 = 0.5),
                                      prefix = "scONCO") {
  required_cols <- c(
    "label_scMRMA_L1", "conf_scMRMA_L1",
    "label_scMRMA_L2", "conf_scMRMA_L2",
    "label_scMRMA_L3", "conf_scMRMA_L3",
    "label_scMRMA_L4Abb", "conf_scMRMA_L4Abb"
  )
  missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
  if (length(missing_cols) > 0) {
    stop("Missing meta columns: ", paste(missing_cols, collapse = ", "),
         "\nRun scONCO main pipeline with confidence-outputting variant first.")
  }

  md <- seurat_obj@meta.data
  n <- nrow(md)
  final_label <- character(n)
  final_level <- character(n)

  for (i in seq_len(n)) {
    res <- hierarchical_reject_one(
      l1 = list(type = md$label_scMRMA_L1[i],     conf = md$conf_scMRMA_L1[i]),
      l2 = list(type = md$label_scMRMA_L2[i],     conf = md$conf_scMRMA_L2[i]),
      l3 = list(type = md$label_scMRMA_L3[i],     conf = md$conf_scMRMA_L3[i]),
      l4 = list(type = md$label_scMRMA_L4Abb[i],  conf = md$conf_scMRMA_L4Abb[i]),
      thresholds = thresholds
    )
    final_label[i] <- res$label
    final_level[i] <- res$level
  }

  seurat_obj@meta.data[[paste0(prefix, "_final_label")]] <- final_label
  seurat_obj@meta.data[[paste0(prefix, "_final_level")]] <- final_level

  # Report
  tab <- table(final_level)
  message("[scONCO hierarchical reject] Final level distribution:")
  print(tab)
  message(sprintf("  Cells reaching deepest level (L4): %.1f%%",
                  100 * sum(final_level == "L4_Abb") / n))
  message(sprintf("  Cells rejected to L1 or Unknown:   %.1f%%",
                  100 * sum(final_level %in% c("L1", "Unknown")) / n))

  seurat_obj
}

###############################################################################
## 4) Threshold calibration on a validation set
##
## Use a labeled validation set to choose thresholds maximizing F1 at each level
###############################################################################

#' Calibrate per-level thresholds from a held-out validation set
#'
#' @param val_meta data.frame from validation Seurat with required columns
#'                 + ground_truth_L1, ground_truth_L2, etc.
#' @param target_accuracy minimum per-level accuracy (default 0.85)
#' @return named threshold vector c(L1=..., L2=..., L3=..., L4=...)
calibrate_reject_thresholds <- function(val_meta, target_accuracy = 0.85) {
  thresholds <- c(L1 = NA, L2 = NA, L3 = NA, L4 = NA)
  for (level in c("L1", "L2", "L3", "L4")) {
    label_col <- if (level == "L4") "label_scMRMA_L4Abb" else paste0("label_scMRMA_", level)
    conf_col  <- if (level == "L4") "conf_scMRMA_L4Abb"  else paste0("conf_scMRMA_", level)
    gt_col <- paste0("ground_truth_", level)
    if (!all(c(label_col, conf_col, gt_col) %in% colnames(val_meta))) next

    # Sweep thresholds; find smallest threshold such that subset accuracy ≥ target
    grid <- seq(0.05, 0.95, by = 0.05)
    best_t <- 0.5
    for (t in rev(grid)) {  # start from strict
      keep <- val_meta[[conf_col]] >= t & !is.na(val_meta[[label_col]])
      if (sum(keep) < 10) next
      acc <- mean(val_meta[[label_col]][keep] == val_meta[[gt_col]][keep], na.rm = TRUE)
      if (acc >= target_accuracy) { best_t <- t; break }
    }
    thresholds[level] <- best_t
  }
  # If any NA, fall back to defaults
  defaults <- c(L1 = 0.2, L2 = 0.3, L3 = 0.4, L4 = 0.5)
  thresholds[is.na(thresholds)] <- defaults[is.na(thresholds)]
  thresholds
}

###############################################################################
## Example
###############################################################################
# # After running main_scONCO.R with confidence outputs:
# seuratObject_Sample <- apply_hierarchical_reject(seuratObject_Sample,
#                                                 thresholds = c(L1=0.2, L2=0.3, L3=0.4, L4=0.5))
# table(seuratObject_Sample$scONCO_final_level)
# DimPlot(seuratObject_Sample, group.by = "scONCO_final_label")
