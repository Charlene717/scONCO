#' Run the scONCO hierarchical annotation pipeline
#'
#' Top-level entry point. Annotates a pre-clustered Seurat object using the
#' 4-level pan-cancer marker database (Level1 -> Level4_Abb) and writes one
#' meta.data column per level.
#'
#' @param seurat_obj A Seurat object that has been QC'd, normalised, scaled,
#'   PCA-projected, and clustered (so that `meta.data[[cluster_col]]` exists).
#' @param markers_df Marker database in long format with columns
#'   `gene`, `Level1`, `Level2`, `Level3`, `Level4_Abb`, optionally `sign`,
#'   `ai_support`, `lit_support`, `specificity_index`.
#' @param cluster_col Name of the cluster column in `meta.data`. Default
#'   `"seurat_clusters"`.
#' @param hierarchy Character vector of column names in `markers_df` defining
#'   the hierarchy from coarse to fine. Default
#'   `c("Level1", "Level2", "Level3", "Level4_Abb")`.
#' @param weights Named numeric vector of scoring weights:
#'   `c(alpha = 0.6, beta = 0.2, gamma = 0.2, delta = 0.2)`.
#' @param apply_confidence Logical. If `TRUE`, also runs `apply_hierarchical_reject()`
#'   and `apply_conformal_to_seurat()` with default parameters.
#' @param apply_cnv_gate Logical. If `TRUE`, also runs the CNV-aware malignant
#'   gate ([apply_cnv_malignant_gate()]), adding `scONCO_malignant_call`,
#'   `scONCO_malignant_prob`, `scONCO_cnv_score`, `scONCO_L1_gated`,
#'   `scONCO_L4_gated`, `scONCO_gate_corrected`.
#' @param cnv_method,cnv_col,cnv_score Passed to the CNV gate. With
#'   `cnv_method = "precomputed"` (default) supply `cnv_col` (a meta.data column
#'   holding a CopyKAT/scATOMIC/inferCNV call) or a precomputed `cnv_score`.
#' @param verbose Whether to print progress.
#'
#' @return Seurat object with extra meta.data columns:
#'   `scONCO_L1` ... `scONCO_L4` (or whatever `hierarchy` names),
#'   and (if `apply_confidence`) `scONCO_L{k}_confidence`,
#'   `scONCO_L{k}_status`, `conformal_set`, `conformal_set_size`.
#'
#' @examples
#' \dontrun{
#'   library(Seurat); library(scONCO)
#'   markers_df <- load_cancer_marker_db()
#'   seu <- run_scONCO(seu, markers_df)
#'   table(seu$scONCO_L4)
#' }
#'
#' @export
run_scONCO <- function(seurat_obj,
                        markers_df,
                        cluster_col      = "seurat_clusters",
                        hierarchy        = c("Level1","Level2","Level3","Level4_Abb"),
                        weights          = c(alpha = 0.6, beta = 0.2, gamma = 0.2, delta = 0.2),
                        apply_confidence = FALSE,
                        apply_cnv_gate   = FALSE,
                        cnv_method       = "precomputed",
                        cnv_col          = NULL,
                        cnv_score        = NULL,
                        verbose          = TRUE) {

  if (verbose) message("scONCO: hierarchical annotation across ",
                       length(hierarchy), " levels...")

  if (!exists("sconco_annotate_seurat_clusters", mode = "function")) {
    stop("scONCO core not loaded. Set SCONCO_LEGACY_DIR or run usethis::create_package().")
  }

  seurat_obj <- sconco_annotate_seurat_clusters(
    seurat_obj  = seurat_obj,
    markers_df  = markers_df,
    cluster_col = cluster_col,
    hierarchy   = hierarchy,
    weights     = weights
  )

  if (apply_confidence) {
    if (verbose) message("scONCO: applying v1.1 confidence layer...")
    seurat_obj <- apply_hierarchical_reject(seurat_obj)
    seurat_obj <- apply_conformal_to_seurat(seurat_obj, alpha = 0.1)
  }

  if (apply_cnv_gate) {
    if (verbose) message("scONCO: applying CNV-aware malignant gate...")
    if (!exists("apply_cnv_malignant_gate", mode = "function"))
      stop("CNV gate not loaded (code/R/13_FUN_cnv_malignant_gate.R).")
    seurat_obj <- apply_cnv_malignant_gate(
      seurat_obj, cnv_score = cnv_score,
      method = cnv_method, cnv_col = cnv_col, verbose = verbose)
  }

  if (verbose) message("scONCO: done.")
  invisible(seurat_obj)
}
