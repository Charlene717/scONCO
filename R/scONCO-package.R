#' scONCO: Pan-Cancer Single-Cell Annotation Framework
#'
#' scONCO combines a 4-level hierarchical pan-cancer marker database, a
#' specificity-weighted scoring algorithm extending scMRMA, and a statistical
#' confidence layer (hierarchical reject, split conformal, AUC specificity,
#' empirical AI/literature support) to annotate cancer single-cell RNA-seq data.
#'
#' @section Core functions:
#' * [run_scONCO()] — top-level hierarchical pipeline.
#' * [sconco_annotate_seurat()] — single-pass annotation.
#' * [load_cancer_marker_db()] — load the curated marker database.
#'
#' @section v1.1 statistical-confidence layer:
#' * [apply_hierarchical_reject()] — defer ambiguous predictions to parent level.
#' * [apply_conformal_to_seurat()] — split conformal prediction sets.
#' * [compute_marker_AUC_batch()] — Wilcoxon AUC marker specificity.
#' * [upgrade_to_scONCO_v1_1()] — one-call upgrade wrapper.
#'
#' @section Cancer layer (CNV-aware malignant gate):
#' * [apply_cnv_malignant_gate()] — fuse marker + CNV + module into a per-cell
#'   malignant call (the central malignant-vs-normal-origin problem).
#' * [compute_cnv_score()] — per-cell CNV signal from CopyKAT / inferCNV / SCEVAN.
#' * [evaluate_malignant_gate()] — malignant-calling accuracy vs ground-truth.
#'
#' @section Visualization:
#' * [plot_marker_bubble()]
#' * [plot_marker_dotplots()]
#'
#' @section Citation:
#' Chang C-J., Dai Y-H., et al. (2027). scONCO: a hierarchical, AI-augmented,
#' statistically calibrated marker-based annotation framework for cancer
#' single-cell transcriptomics. *Nucleic Acids Research* (in preparation).
#'
#' @section Imports:
#' The `@importFrom` directives below are the single source of truth for the
#' package NAMESPACE. They were moved here on 2026-09-14: the NAMESPACE file had
#' been maintained by hand ("Manual draft for v1.1 skeleton") and carried 21
#' `importFrom` entries that no roxygen tag reproduced, so the first
#' `devtools::document()` run would silently have dropped every one of them.
#' Add new imports here, not to NAMESPACE.
#'
#' `R CMD check` reports a NOTE that `Matrix`, `ggplot2`, `methods`, `stringr`
#' and `tibble` are declared in Imports but never imported from. This is an
#' artefact of the skeleton-phase legacy bridge, not a stale declaration: the
#' implementations that use them live in `inst/legacy/R/` and are sourced at
#' run time by `zzz_source_legacy.R`, so `R CMD check` — which scans only `R/` —
#' cannot see the usage. Verified 2026-09-14: `Matrix::rowMeans`, `ggplot2::aes`
#' and friends, `stringr::str_trim` and `tibble::tibble` are all called from the
#' legacy modules, and `methods` is required for S4 slot access on Seurat
#' objects (`seurat_obj@meta.data`). Do not remove them to silence the NOTE.
#' `tidyr` genuinely was unused and was dropped from Imports on 2026-09-14.
#' The NOTE resolves itself once the legacy modules are moved into `R/`.
#'
#' @importFrom Seurat AddMetaData FetchData FindClusters FindNeighbors
#' @importFrom Seurat FindVariableFeatures GetAssayData NormalizeData
#' @importFrom Seurat RunPCA RunUMAP ScaleData
#' @importFrom dplyr %>% filter group_by mutate summarise
#' @importFrom stats median quantile sd wilcox.test
#' @importFrom stats na.omit setNames
#' @importFrom utils head tail
#'
#' @docType package
#' @name scONCO-package
#' @aliases scONCO
"_PACKAGE"
