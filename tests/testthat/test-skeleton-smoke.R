# Smoke tests for the v1.1 skeleton.
# Real numerical tests are in ../code/tests/test_v1_1_modules.R; this file
# only checks that the exported function names resolve and that the
# legacy bridge has loaded the implementations.

test_that("scONCO package attaches without error", {
  expect_true("scONCO" %in% loadedNamespaces())
})

test_that("Core v1.0 functions are exported", {
  exports <- c(
    "run_scONCO",
    "sconco_compute_marker_weights",
    "sconco_activity_scores",
    "sconco_annotate",
    "sconco_annotate_seurat",
    "sconco_annotate_seurat_clusters",
    "annotate_broad_cell_clusters",
    "normalize_celltype",
    "map_to_canonical",
    "plot_marker_dotplots",
    "plot_marker_bubble",
    "auto_recluster",
    "merge_similar_numbered_celltypes",
    "merge_similar_numbered_celltypes_by_correlation"
  )
  for (e in exports) {
    expect_true(exists(e, mode = "function", envir = asNamespace("scONCO")),
                info = paste("missing export:", e))
  }
})

test_that("v1.1 confidence-layer functions are exported", {
  exports <- c(
    "compute_marker_AUC_batch", "compute_tau_score", "compute_mutex_score",
    "update_specificity_empirical",
    "compute_confidence", "hierarchical_reject_one",
    "apply_hierarchical_reject", "calibrate_reject_thresholds",
    "softmax_matrix", "split_conformal_calibrate",
    "conformal_predict", "apply_conformal_to_seurat",
    "compute_empirical_ai_support", "compute_empirical_lit_support",
    "refresh_marker_weights_empirical", "upgrade_to_scONCO_v1_1"
  )
  for (e in exports) {
    expect_true(exists(e, mode = "function", envir = asNamespace("scONCO")),
                info = paste("missing export:", e))
  }
})

test_that("load_cancer_marker_db returns a data.frame with required columns", {
  skip_on_cran()
  skip_if_not(file.exists(system.file("..", "..", "database", "current",
                                       "DB_pancancer_human_v1.0.R",
                                       package = "scONCO")),
              "DB v1.0 file not present in this skeleton install")
  m <- load_cancer_marker_db()
  expect_s3_class(m, "data.frame")
  for (col in c("gene", "Level1", "Level2", "Level3", "Level4_Abb")) {
    expect_true(col %in% colnames(m), info = paste("missing column:", col))
  }
})
