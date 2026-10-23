# Generated from code/R/16_FUN_sconco.R by tools/sync_package.py - do not edit here:
# edit the source and re-run the script.
###############################################################################
# 16_FUN_sconco.R - scONCO end-to-end entry point.
#
#   sconco(counts)                 counts: genes x cells raw counts (dgCMatrix
#                                   or matrix) of ONE sample (the method is
#                                   per sample: the tumour test and the copy-
#                                   number reference are sample-level objects)
#   sconco_seurat(seu)             the same on a Seurat object, results written
#                                   to its meta.data (sconco_L1 ... sconco_status)
#   sconco_study(counts, sample)   every sample of a study, each
#                                   given the study's immune reference pool
#                                   (sconco_reference_pool) for use only where
#                                   its own immune reference is < 20 cells
#
# The revised version (docs/METHOD_DEVELOPMENT_PLAN.md, the revision): database v3.3 with
# the glioblastoma state modules used only to name cells the copy-number step
# has declared malignant; the cross-sample immune pool. Everything a single
# sample with its own reference gets is what the initial version gave it.
#
# Pipeline (docs/METHOD_DEVELOPMENT_PLAN.md):
#   C4 evidence-weighted marker scores (database v3.3, w = spec x log2(1+n_sources))
#   C3 multi-resolution hierarchical annotation (Level 1-4, cluster-wise, reject option)
#   C1 own copy-number signal on the epithelial-lineage node, immune/stromal reference
#   C2 sample-level tumour-presence test; no malignant call when the test is negative
#   C5 P(malignant) per epithelial-lineage cell from the C1 mixture posterior
#
# Defaults resolve relative to SCONCO_ROOT (or the working directory):
#   database/current_v3_3/DB_pancancer_human_v3.3.R and
#   database/reference_dbs/gene_positions.tsv
###############################################################################

.sconco_root <- function() Sys.getenv("SCONCO_ROOT", unset = getwd())

# package copy only: the shipped files first, the project folder (SCONCO_ROOT or the working
# directory) as before when the package is not installed
.sconco_extdata <- function(name, repo_dir, repo_name = name) {
  f <- system.file("extdata", name, package = "scONCO")
  if (nzchar(f)) f else file.path(.sconco_root(), repo_dir, repo_name) }

# The revised version ships database v3.3. The initial version used v3.2; to reproduce it pass
# sconco_load_db(system.file("extdata", "DB_pancancer_human_v3.2.R", package = "scONCO")).
#' Load a scONCO marker database
#'
#' Database v3.3 (the default) is the one the method ships; v3.2 reproduces the
#' initial version: `sconco_load_db(system.file("extdata",
#' "DB_pancancer_human_v3.2.R", package = "scONCO"))`. Every core marker has at least
#' three independent public sources; the per-marker evidence is in the project
#' repository (database/current_v3_3/).
#'
#' @param file Path to a database file (an R script defining `markers_df`).
#' @return A data frame, one row per marker entry (`gene`, `Level1`-`Level4_Abb`,
#'   specificity, number of sources, `part`, ...).
#' @export
sconco_load_db <- function(file = .sconco_extdata("DB_pancancer_human_v3.3.R", "database/current_v3_3")) {
  e <- new.env(); suppressMessages(source(file, local = e)); as.data.frame(e$markers_df) }

#' Load gene positions for the copy-number step
#'
#' @param file Path to a tab-separated file with columns `gene`, `chrom`, `start`,
#'   `end` (gzip allowed). The default is the Ensembl table shipped with the package.
#' @return A data frame.
#' @export
sconco_load_gene_positions <- function(file = .sconco_extdata("gene_positions.tsv.gz", "database/reference_dbs", "gene_positions.tsv"))
  utils::read.delim(file, stringsAsFactors = FALSE)

#' Annotate one tumour sample with scONCO
#'
#' Runs the whole method on one sample: evidence-weighted marker scores against the
#' marker database, multi-resolution hierarchical annotation (Levels 1-4), and, on the
#' epithelial-lineage node, a copy-number signal computed within the sample against
#' its own immune cells, preceded by a sample-level test of whether a tumour is
#' present. When the test is negative no cell is called malignant.
#'
#' The method is per sample: the tumour test and the copy-number reference are
#' sample-level objects. Run [sconco_study()] for several samples of one study.
#'
#' @param counts Raw counts of ONE sample, genes x cells (`dgCMatrix` or matrix),
#'   gene symbols as row names.
#' @param db Marker database (a data frame as returned by [sconco_load_db()]);
#'   `NULL` loads database v3.3, shipped with the package.
#' @param gene_pos Gene positions (columns `gene`, `chrom`, `start`, `end`); `NULL`
#'   loads the Ensembl positions shipped with the package.
#' @param resolutions Clustering resolutions for Levels 1, 2 and 3.
#' @param weight_evidence Weight each marker by `log2(1 + n_sources)` as well as its
#'   specificity (`FALSE` reproduces the unweighted ablation).
#' @param cnv Run the copy-number step and the tumour test. With `FALSE`, Level 2 of
#'   the epithelial lineage is not decided.
#' @param cnv_args List of further arguments for [sconco_cnv()], e.g. the two options
#'   of the revision that are off by default: `exclude_keratinocyte = TRUE` and
#'   `amp_min = 0.12`.
#' @param normalized `TRUE` if `counts` is already log-normalised.
#' @param reference_pool Optional cross-sample immune reference built by
#'   [sconco_reference_pool()]; used only when the sample has fewer than
#'   `pool_min_own` immune cells of its own.
#' @param pool_min_own Own immune cells below which the pool is used.
#' @param seed Random seed (clustering).
#' @param verbose Print progress.
#' @return An object of class `sconco`: a list with `annotation` (one row per cell:
#'   `L1`-`L4`, `conf`, `status`, and, where the copy-number step ran,
#'   `cnv_burden` and `p_malignant`), `tumour_test` (the sample-level test),
#'   `verdict` (one sentence), `probs`, `n_cells`, `minutes` and `settings`.
#' @seealso [sconco_study()], [sconco_seurat()], [sconco_load_db()]
#' @examples
#' \dontrun{
#' res <- sconco(counts)          # counts: one sample, genes x cells
#' res                            # verdict and cell counts per level
#' head(res$annotation)
#' }
#' @export
sconco <- function(counts, db = NULL, gene_pos = NULL, resolutions = c(0.2, 0.5, 1.2),
                    weight_evidence = TRUE, cnv = TRUE, cnv_args = list(), normalized = FALSE,
                    reference_pool = NULL, pool_min_own = 20L, seed = 1L, verbose = TRUE) {
  if (is.null(db)) db <- sconco_load_db()
  if (is.null(gene_pos) && cnv) gene_pos <- sconco_load_gene_positions()
  t0 <- Sys.time()
  h <- sconco_hier(counts, db, gene_pos, resolutions = resolutions, weight_evidence = weight_evidence,
                   cnv = cnv, cnv_args = cnv_args, normalized = normalized, reference_pool = reference_pool,
                   pool_min_own = pool_min_own, seed = seed, verbose = verbose)
  a <- h$annotation
  tt <- h$tumour_test
  verdict <- if (is.null(tt)) "not tested (no copy-number run: too few epithelial-lineage or reference cells)"
             else if (isTRUE(tt$below_amplitude_floor)) sprintf("no tumour detected (%.0f%% of candidates above the diploid null, but no cluster reaches the copy-number amplitude floor: largest arm deviation %.3f, floor %.2f) - no malignant calls made", 100 * tt$frac_above_null, tt$amplitude_max, tt$amp_min)
             else if (isTRUE(tt$present) && isTRUE(tt$n_pool > 0)) sprintf("tumour present (%.0f%% of epithelial-lineage cells above the diploid null; reference includes %d pooled immune cells from other samples)", 100 * tt$frac_above_null, tt$n_pool)
             else if (isTRUE(tt$present)) sprintf("tumour present (%.0f%% of epithelial-lineage cells above the diploid null)", 100 * tt$frac_above_null)
             else sprintf("no tumour detected (%.0f%% above the diploid null; threshold %.0f%%) - no malignant calls made", 100 * tt$frac_above_null, 100 * tt$frac_min)
  out <- list(annotation = a, tumour_test = tt, verdict = verdict, probs = h$probs,
              n_cells = nrow(a), minutes = round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2),
              settings = list(resolutions = resolutions, weight_evidence = weight_evidence, cnv = cnv, cnv_args = cnv_args, seed = seed))
  class(out) <- "sconco"; out }

#' @rdname sconco
#' @param x,object A `sconco` result.
#' @param ... Ignored.
#' @method print sconco
#' @export
print.sconco <- function(x, ...) {
  a <- x$annotation
  cat("scONCO annotation:", x$n_cells, "cells,", x$minutes, "min\n")
  cat("  sample verdict:", x$verdict, "\n")
  cat("  Level 1:", paste(names(table(a$L1)), table(a$L1), collapse = "; "), "\n")
  epi <- a[!is.na(a$L2) & a$L2 %in% c("Malignant", "Normal epithelium", "Epithelial (undetermined)"), ]
  if (nrow(epi)) cat("  epithelial lineage:", paste(names(table(epi$L2)), table(epi$L2), collapse = "; "), "\n")
  cat("  Level 4 status:", paste(names(table(a$status)), table(a$status), collapse = "; "), "\n")
  invisible(x) }

#' @rdname sconco
#' @method summary sconco
#' @export
summary.sconco <- function(object, ...) {
  a <- object$annotation
  s <- as.data.frame(table(L1 = a$L1, L2 = a$L2, L3 = a$L3, L4 = a$L4), stringsAsFactors = FALSE)
  s <- s[s$Freq > 0, ]; s <- s[order(s$L1, s$L2, s$L3, -s$Freq), ]; rownames(s) <- NULL
  names(s)[5] <- "n_cells"; s }

# ---- revision: cross-sample reference pool ----------------------------------
# sconco_reference_pool(): the immune-node cells of every sample of a study,
# annotated by the marker hierarchy alone (no copy number), capped per sample.
# The pool is the method's own annotation, not a user claim: a cell enters only
# because the same rule that would make it a reference in its own sample
# (Level 1 = Immune cell) applies. sconco_study() runs sconco() on every
# sample of a study, giving each sample the pool minus its own cells; the pool
# is used only where a sample's own immune reference is below pool_min_own.
#' Cross-sample immune reference pool
#'
#' Annotates every sample of a study by the marker hierarchy alone (no copy number)
#' and pools the cells placed at Level 1 = Immune cell, at most `per_sample_cap` per
#' sample. The pool is the method's own annotation, not a user claim. Used by
#' [sconco_study()] for samples whose own immune reference is too small.
#'
#' @param counts Raw counts of all samples, genes x cells.
#' @param sample Sample of each cell (length `ncol(counts)`).
#' @param db Marker database; `NULL` loads v3.3.
#' @param per_sample_cap Maximum immune cells taken from one sample.
#' @param min_cells Samples with fewer cells are skipped.
#' @param resolutions,weight_evidence,seed,verbose As in [sconco()].
#' @return A counts matrix of the pooled cells, with attribute `origin` (sample of
#'   each cell).
#' @export
sconco_reference_pool <- function(counts, sample, db = NULL, per_sample_cap = 300L, min_cells = 60L,
                                   resolutions = c(0.2, 0.5, 1.2), weight_evidence = TRUE, seed = 1L, verbose = TRUE) {
  if (is.null(db)) db <- sconco_load_db()
  stopifnot(length(sample) == ncol(counts))
  keep <- character(0); origin <- character(0)
  for (s in unique(sample)) {
    i <- which(sample == s); if (length(i) < min_cells) next
    h <- sconco_hier(counts[, i], db, NULL, resolutions = resolutions, weight_evidence = weight_evidence, cnv = FALSE, seed = seed, verbose = FALSE)
    imm <- rownames(h$annotation)[h$annotation$L1 == "Immune cell"]
    if (length(imm) > per_sample_cap) { set.seed(seed); imm <- sample(imm, per_sample_cap) }
    keep <- c(keep, imm); origin <- c(origin, rep(s, length(imm)))
    if (verbose) message("  pool: ", s, " contributes ", length(imm), " immune cells")
  }
  structure(counts[, keep, drop = FALSE], origin = setNames(origin, keep))
}

#' Run scONCO on every sample of a study
#'
#' Builds the cross-sample immune reference pool ([sconco_reference_pool()]) and runs
#' [sconco()] on each sample (60 cells or more), giving each sample the pool for use
#' only where its own immune reference has fewer than `pool_min_own` cells.
#'
#' @param counts Raw counts of all samples, genes x cells.
#' @param sample Sample of each cell.
#' @param db,gene_pos As in [sconco()].
#' @param per_sample_cap Passed to [sconco_reference_pool()].
#' @param pool_min_own,verbose As in [sconco()].
#' @param ... Passed to [sconco()].
#' @return A named list of `sconco` results, one per sample.
#' @export
sconco_study <- function(counts, sample, db = NULL, gene_pos = NULL, per_sample_cap = 300L, pool_min_own = 20L, verbose = TRUE, ...) {
  if (is.null(db)) db <- sconco_load_db()
  if (is.null(gene_pos)) gene_pos <- sconco_load_gene_positions()
  if (verbose) message("building the cross-sample immune reference pool")
  pool <- sconco_reference_pool(counts, sample, db, per_sample_cap = per_sample_cap, verbose = verbose)
  out <- list()
  for (s in unique(sample)) {
    i <- which(sample == s); if (length(i) < 60) next
    if (verbose) message("\n== ", s, ": ", length(i), " cells")
    out[[s]] <- sconco(counts[, i], db, gene_pos, reference_pool = pool, pool_min_own = pool_min_own, verbose = verbose, ...)
  }
  out
}

#' Run scONCO on a Seurat object
#'
#' Calls [sconco()] on the counts of a Seurat object holding ONE sample and writes
#' the result to its metadata (`sconco_L1` ... `sconco_L4`, `sconco_conf`,
#' `sconco_status`, `sconco_cnv_burden`, `sconco_p_malignant`); the sample-level
#' verdict goes to `seu@misc$sconco`.
#'
#' @param seu A Seurat object (one sample).
#' @param assay Assay to read; `NULL` for the default assay.
#' @param layer Layer holding raw counts.
#' @param ... Passed to [sconco()].
#' @return `seu` with the columns added.
#' @export
sconco_seurat <- function(seu, assay = NULL, layer = "counts", ...) {
  stopifnot(requireNamespace("Seurat", quietly = TRUE))
  counts <- Seurat::GetAssayData(seu, assay = assay, layer = layer)
  r <- sconco(counts, ...)
  a <- r$annotation[colnames(seu), ]
  for (col in c("L1", "L2", "L3", "L4", "conf", "status", "cnv_burden", "p_malignant"))
    if (col %in% names(a)) seu[[paste0("sconco_", col)]] <- a[[col]]
  seu@misc$sconco <- r[c("tumour_test", "verdict", "settings", "minutes")]
  seu }
