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
#' single-cell transcriptomics. *Nucleic Acids Research* [in preparation].
#'
#' @docType package
#' @name scONCO-package
#' @aliases scONCO
"_PACKAGE"
