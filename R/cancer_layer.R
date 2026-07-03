# =============================================================
# scONCO cancer layer — CNV-aware malignant gate (A3)
# (Implementation lives in code/R/13_FUN_cnv_malignant_gate.R and is bridged
#  via zzz_source_legacy.R until devtools::document() is run.)
# =============================================================

#' Compute a per-cell CNV (aneuploidy) score in [0,1]
#'
#' Obtains a genome-wide copy-number signal per cell, used by the malignant
#' gate as the gold-standard malignancy cue. Either reads a precomputed CNV
#' call from `meta.data` (e.g. from CopyKAT / scATOMIC / inferCNV wrappers) or
#' runs the chosen CNV tool directly.
#'
#' @param seurat_obj A Seurat object with raw counts in the RNA assay.
#' @param method One of `"precomputed"`, `"copykat"`, `"infercnv"`, `"scevan"`.
#' @param cnv_col `meta.data` column with a precomputed CNV call/score
#'   (used when `method = "precomputed"`).
#' @param ... Passed to the underlying CNV tool.
#' @return Named numeric vector (cell barcode -> score in `[0,1]`), higher =
#'   more aneuploid / malignant.
#' @references Gao, R. et al. (2021) CopyKAT, *Nat Biotechnol* 39, 599-608;
#'   Tickle, T. et al. (2019) inferCNV, Broad Institute;
#'   De Falco, A. et al. (2023) SCEVAN, *Nat Commun* 14, 1074.
#' @export
compute_cnv_score <- function(seurat_obj, method = "precomputed", cnv_col = NULL, ...) {
  stop("Implementation provided by code/R/13_FUN_cnv_malignant_gate.R ",
       "(loaded on package attach via zzz_source_legacy.R).")
}

#' Apply the CNV-aware malignant gate (A3)
#'
#' Fuses scONCO's marker-based malignant evidence, a per-cell CNV score, and the
#' tumour-type module prior into a per-cell malignant probability
#' `P(malignant) = sigma(w_marker*marker + w_cnv*CNV + w_module*module)`, then
#' classifies each (epithelial-lineage) cell as Malignant / Non_malignant /
#' Uncertain and, optionally, reconciles the marker label with the gate. This is
#' scONCO's answer to the central tumour-scRNA-seq problem of distinguishing
#' malignant cells from their normal cell-of-origin.
#'
#' Adds to `meta.data`: `scONCO_cnv_score`, `scONCO_malignant_prob`,
#' `scONCO_malignant_call`, `scONCO_L1_gated`, `scONCO_L4_gated`,
#' `scONCO_gate_corrected`.
#'
#' @param seurat_obj A Seurat object previously annotated by [run_scONCO()].
#' @param cnv_score Optional named numeric vector (cell -> `[0,1]`). If `NULL`,
#'   computed via [compute_cnv_score()].
#' @param method,cnv_col Passed to [compute_cnv_score()] when `cnv_score` is `NULL`.
#' @param weights Named numeric `c(marker, cnv, module)` (auto-normalised).
#' @param fusion `"weighted"` or `"logistic"`.
#' @param threshold Decision threshold on `P(malignant)` (default 0.5).
#' @param uncertain_margin Half-width of the Uncertain band (0 disables).
#' @param restrict_to_epithelial Gate only epithelial-lineage cells (default TRUE).
#' @param correct_labels Reconcile marker label with the gate (default TRUE).
#' @param verbose Logical.
#' @return Updated Seurat object.
#' @export
apply_cnv_malignant_gate <- function(seurat_obj, cnv_score = NULL,
                                     method = "precomputed", cnv_col = NULL,
                                     weights = c(marker = 0.40, cnv = 0.45, module = 0.15),
                                     fusion = c("weighted", "logistic"),
                                     threshold = 0.5, uncertain_margin = 0.05,
                                     restrict_to_epithelial = TRUE,
                                     correct_labels = TRUE, verbose = TRUE) {
  stop("Implementation provided by code/R/13_FUN_cnv_malignant_gate.R.")
}

#' Evaluate gate malignant calls against a ground-truth flag
#'
#' @param seurat_obj A gated Seurat object (output of [apply_cnv_malignant_gate()]).
#' @param truth_col `meta.data` column with the ground-truth cell type or
#'   malignant flag (malignancy inferred by keyword for a cell-type column).
#' @return Named numeric `c(accuracy, sensitivity, specificity, f1)`.
#' @export
evaluate_malignant_gate <- function(seurat_obj, truth_col) {
  stop("Implementation provided by code/R/13_FUN_cnv_malignant_gate.R.")
}
