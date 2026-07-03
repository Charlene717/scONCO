# =============================================================
# v1.0 core layer — roxygen stubs.
# Implementations bridged from code/R/01-08_FUN_*.R via
# zzz_source_legacy.R until full devtools::document() refactor.
# =============================================================

#' Compute marker weights (v1.0 core)
#'
#' weight(g) = alpha * specificity(g) + beta * ai_support(g) + gamma * lit_support(g)
#'
#' @param markers_df Marker data frame.
#' @param alpha,beta,gamma Numeric weights. Defaults 0.6, 0.2, 0.2.
#' @return `markers_df` with a `weight` column added.
#' @export
sconco_compute_marker_weights <- function(markers_df,
                                           alpha = 0.6, beta = 0.2,
                                           gamma = 0.2) {
  stop("Implementation provided by code/R/01_FUN_scONCO_core.R.")
}

#' Cell-type activity scores
#'
#' score(cell, type) = sum_{g in M+(t)} weight(g) * expr(g, cell)
#'                   - delta * sum_{g in M-(t)} expr(g, cell)
#'
#' @export
sconco_activity_scores <- function(expr, markers_df, delta = 0.2) {
  stop("Implementation provided by code/R/01_FUN_scONCO_core.R.")
}

#' Annotate cells directly from a score matrix
#' @export
sconco_annotate <- function(score_matrix, hierarchy_table = NULL) {
  stop("Implementation provided by code/R/01_FUN_scONCO_core.R.")
}

#' Annotate cells in a Seurat object (single-pass)
#'
#' Convenience wrapper for cell-level annotation; for cluster-level annotation
#' (the default in run_scONCO), use [sconco_annotate_seurat_clusters()].
#'
#' @param seurat_obj A Seurat object.
#' @param markers_df Marker DB.
#' @param assay Assay to pull expression from. Default `"RNA"`.
#' @return Seurat object with `scONCO_celltype` meta.data column added.
#' @export
sconco_annotate_seurat <- function(seurat_obj, markers_df, assay = "RNA") {
  stop("Implementation provided by code/R/01_FUN_scONCO_core.R.")
}

#' Annotate cells in a Seurat object using existing clusters
#'
#' Hierarchical pipeline used internally by [run_scONCO()].
#'
#' @param seurat_obj A Seurat object with cluster IDs.
#' @param markers_df Marker DB.
#' @param cluster_col Cluster column name.
#' @param hierarchy Character vector of level columns in `markers_df`.
#' @param weights Named numeric vector.
#' @return Seurat object with one column per level in `hierarchy`.
#' @export
sconco_annotate_seurat_clusters <- function(seurat_obj, markers_df,
                                              cluster_col = "seurat_clusters",
                                              hierarchy   = c("Level1","Level2","Level3","Level4_Abb"),
                                              weights     = c(alpha=0.6, beta=0.2, gamma=0.2, delta=0.2)) {
  stop("Implementation provided by code/R/01_FUN_scONCO_core.R.")
}

#' Map cluster Level4 labels to broad cell types
#'
#' Useful for downstream plots that need only the 5 broad lineages.
#'
#' @param seurat_obj A Seurat object with `scONCO_L4` or equivalent column.
#' @param mapping_df Optional override for the L4 → broad map.
#' @return Seurat object with `broad_cell_type` meta.data column.
#' @export
annotate_broad_cell_clusters <- function(seurat_obj, mapping_df = NULL) {
  stop("Implementation provided by code/R/02_FUN_annotate_broad_clusters.R.")
}

#' Normalize a free-text cell type to scONCO canonical name
#' @export
normalize_celltype <- function(x) {
  stop("Implementation provided by code/R/03_FUN_celltype_alias_cancer.R.")
}

#' Map a cell type to its canonical scONCO name
#' @export
map_to_canonical <- function(x) {
  stop("Implementation provided by code/R/03_FUN_celltype_alias_cancer.R.")
}

#' Per-cell-type DotPlot panel
#' @export
plot_marker_dotplots <- function(seurat_obj, markers_df, group_by = "scONCO_L4") {
  stop("Implementation provided by code/R/04_FUN_marker_dotplots.R.")
}

#' scONCO-style marker bubble plot
#' @export
plot_marker_bubble <- function(seurat_obj, markers_df, group_by = "scONCO_L4",
                                    color_low = "#3B82F6", color_high = "#EF4444") {
  stop("Implementation provided by code/R/05_FUN_marker_bubble.R.")
}

#' Auto-reclustering by resolution gradient + ROGUE
#' @export
auto_recluster <- function(seurat_obj, resolutions = seq(0.2, 1.6, by = 0.2),
                            min_rogue = 0.85) {
  stop("Implementation provided by code/R/06_FUN_auto_reclustering.R.")
}

#' Merge similar numbered cell-type clusters by Spearman correlation
#' @export
merge_similar_numbered_celltypes_by_correlation <- function(seurat_obj,
                                                             threshold = 0.85,
                                                             top_genes = 50) {
  stop("Implementation provided by code/R/07_FUN_merge_clusters_correlation.R.")
}

#' Merge similar numbered cell-type clusters by marker Jaccard
#' @export
merge_similar_numbered_celltypes <- function(seurat_obj, markers_df,
                                              jaccard_threshold = 0.8) {
  stop("Implementation provided by code/R/08_FUN_merge_clusters_marker.R.")
}
