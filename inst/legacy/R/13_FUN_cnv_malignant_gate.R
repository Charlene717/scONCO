###############################################################################
## 13_FUN_cnv_malignant_gate.R — scONCO cancer layer (A3)
##
## CNV-aware malignant-vs-normal gate. The central difficulty of tumour scRNA-seq
## is telling a MALIGNANT epithelial cell apart from its NORMAL cell-of-origin —
## a problem marker-based annotators cannot solve from expression alone, because
## malignant cells keep expressing their lineage markers (EPCAM/KRT...). scONCO
## fuses three signals into a per-cell malignant probability:
##
##   P(malignant | cell) = sigma( w_marker * marker + w_cnv * CNV + w_module * module )
##
##   - marker : scONCO marker-based malignant evidence (Level1 == "Malignant cell"
##              weighted by its confidence)
##   - CNV    : genome-wide copy-number signal from CopyKAT / inferCNV / SCEVAN
##              (aneuploid -> malignant), the gold-standard malignancy cue
##   - module : tumour-type prior from the Level4 module (Mal_<type> high, normal
##              epithelial neutral, immune/stromal/endothelial low)
##
## The gate is applied only within the EPITHELIAL compartment (malignant + normal
## epithelial); immune / stromal / endothelial cells are never malignant epithelial
## tumour cells and are passed through as Non_malignant. When `correct_labels`
## is TRUE, disagreements between marker call and the gate are reconciled (a
## "Malignant" marker call without CNV support is demoted to normal epithelial,
## and a normal-epithelial call with strong CNV support is promoted to malignant).
##
## Outputs added to meta.data:
##   scONCO_cnv_score        per-cell CNV signal in [0,1]
##   scONCO_malignant_prob   fused P(malignant) in [0,1]
##   scONCO_malignant_call   "Malignant" / "Non_malignant" / "Uncertain"
##   scONCO_L1_gated         Level1 after gating
##   scONCO_L4_gated         Level4_Abb after gating
##   scONCO_gate_corrected   TRUE where the gate changed the marker call
###############################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    message("[scONCO cancer layer] Seurat not installed; functions will error if called.")
  }
})

# ---------------------------------------------------------------------------
# 0) helpers — detect scONCO output columns robustly
# ---------------------------------------------------------------------------
.sconco_detect_col <- function(md, candidates, pattern = NULL) {
  hit <- intersect(candidates, colnames(md))
  if (length(hit) > 0L) return(hit[1])
  if (!is.null(pattern)) {
    g <- grep(pattern, colnames(md), value = TRUE)
    if (length(g) > 0L) return(g[1])
  }
  NA_character_
}

.sconco_l1_col <- function(md)
  .sconco_detect_col(md, c("scONCO_L1", "scONCO_Level1", "label_scMRMA_L1"),
                     pattern = "scONCO_L1|_L1$")
.sconco_l4_col <- function(md)
  .sconco_detect_col(md, c("scONCO_L4_Abb", "scONCO_L4", "scONCO_L4Abb_refined",
                           "label_scMRMA_L4Abb"),
                     pattern = "scONCO_L4|_L4")
.sconco_conf_col <- function(md)
  .sconco_detect_col(md, c("scONCO_L4_confidence", "scONCO_L4_Abb_confidence",
                           "scONCO_confidence"),
                     pattern = "L4.*confidence|_confidence$")

.is_malignant_label <- function(x)
  grepl("malignant|tumou?r|carcinoma|melanoma|glioma|^mal($|_)", x, ignore.case = TRUE)
.is_epithelial_label <- function(x)
  grepl("epithel|hepatocyte|alveolar|goblet|^epi($|_)|^hep$", x, ignore.case = TRUE)

# ---------------------------------------------------------------------------
# 1) compute_cnv_score — obtain a per-cell CNV signal in [0,1]
# ---------------------------------------------------------------------------
#' Compute a per-cell CNV (aneuploidy) score in [0,1]
#'
#' @param seurat_obj Seurat object (raw counts in RNA assay).
#' @param method "precomputed" (read `cnv_col` from meta.data), "copykat",
#'   "infercnv", or "scevan".
#' @param cnv_col meta.data column holding a precomputed CNV call/score
#'   (e.g. CopyKAT "aneuploid"/"diploid", or a numeric burden). Used when
#'   method = "precomputed".
#' @param ... Passed to the underlying CNV tool.
#' @return Named numeric vector (cell barcode -> score in [0,1]); higher = more
#'   aneuploid / malignant.
compute_cnv_score <- function(seurat_obj,
                              method  = c("precomputed", "copykat", "infercnv", "scevan"),
                              cnv_col = NULL, ...) {
  method <- match.arg(method)
  cells  <- colnames(seurat_obj)

  .to01 <- function(v) {
    if (is.numeric(v)) {
      rng <- range(v, na.rm = TRUE)
      return(if (diff(rng) > 0) (v - rng[1]) / diff(rng) else rep(0.5, length(v)))
    }
    sv <- as.character(v)
    ifelse(grepl("aneuploid|malign|tumou?r|^mal", sv, ignore.case = TRUE), 1,
           ifelse(grepl("diploid|normal|non[_-]?malign", sv, ignore.case = TRUE), 0, NA_real_))
  }

  if (method == "precomputed") {
    if (is.null(cnv_col) || !cnv_col %in% colnames(seurat_obj@meta.data))
      stop("method='precomputed' needs `cnv_col` to be a meta.data column ",
           "(e.g. a CopyKAT/scATOMIC/inferCNV call). Available: ",
           paste(utils::head(colnames(seurat_obj@meta.data), 20), collapse = ", "))
    sc <- .to01(seurat_obj@meta.data[[cnv_col]])
    names(sc) <- rownames(seurat_obj@meta.data)
    return(sc[cells])
  }

  if (method == "copykat") {
    if (!requireNamespace("copykat", quietly = TRUE))
      stop("copykat not installed. remotes::install_github('navinlabcode/copykat')")
    counts <- as.matrix(Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "counts"))
    ck <- copykat::copykat(rawmat = counts, id.type = "S", sam.name = "scONCO_gate",
                           n.cores = 1, output.seg = FALSE, ...)
    pa <- ck$prediction
    sc <- .to01(pa$copykat.pred); names(sc) <- as.character(pa$cell.names)
    out <- sc[cells]; names(out) <- cells
    return(out)
  }

  if (method == "scevan") {
    if (!requireNamespace("SCEVAN", quietly = TRUE))
      stop("SCEVAN not installed. remotes::install_github('AntonioDeFalco/SCEVAN')")
    counts <- as.matrix(Seurat::GetAssayData(seurat_obj, assay = "RNA", slot = "counts"))
    res <- SCEVAN::pipelineCNA(counts, ...)
    sc <- .to01(res$class %||% res[["pred"]]); names(sc) <- rownames(res)
    out <- sc[cells]; names(out) <- cells
    return(out)
  }

  # infercnv: heavy; recommend running benchmarks/tools/run_inferCNV.R and
  # feeding the result back via method = "precomputed".
  stop("method='infercnv' is best run via benchmarks/tools/run_inferCNV.R; ",
       "then pass its call back with method='precomputed', cnv_col=<column>.")
}

`%||%` <- function(a, b) if (!is.null(a)) a else b

# ---------------------------------------------------------------------------
# 2) component signals
# ---------------------------------------------------------------------------
.marker_malignant_prob <- function(md, l1_col, conf_col) {
  l1 <- as.character(md[[l1_col]])
  conf <- if (!is.na(conf_col)) suppressWarnings(as.numeric(md[[conf_col]])) else rep(NA_real_, nrow(md))
  conf[is.na(conf)] <- 0.7  # neutral-ish default when confidence not available
  ifelse(.is_malignant_label(l1), conf, 1 - conf)
}

.module_prior <- function(md, l1_col, l4_col) {
  l1 <- as.character(md[[l1_col]]); l4 <- as.character(md[[l4_col]])
  prior <- rep(0.1, length(l1))                       # immune / stromal / endothelial / other
  prior[.is_epithelial_label(l1) | .is_epithelial_label(l4)] <- 0.4   # normal epithelial: neutral-low
  prior[.is_malignant_label(l1) | grepl("^Mal", l4)] <- 0.8           # malignant module: high
  prior
}

# ---------------------------------------------------------------------------
# 3) apply_cnv_malignant_gate — the gate
# ---------------------------------------------------------------------------
#' Apply the CNV-aware malignant gate to a scONCO-annotated Seurat object (A3)
#'
#' @param seurat_obj Seurat object previously annotated by run_scONCO().
#' @param cnv_score Optional named numeric vector (cell -> [0,1]); if NULL it is
#'   computed via compute_cnv_score(seurat_obj, method, cnv_col).
#' @param method,cnv_col Passed to compute_cnv_score() when cnv_score is NULL.
#' @param weights Named numeric c(marker, cnv, module); auto-normalised to sum 1.
#' @param fusion "weighted" (weighted mean) or "logistic" (sigmoid of weighted
#'   logit), per the manuscript P = sigma(.) form.
#' @param threshold Decision threshold on P(malignant). Default 0.5.
#' @param uncertain_margin Cells with |P - threshold| < margin -> "Uncertain"
#'   (0 disables). Default 0.05.
#' @param restrict_to_epithelial If TRUE (default) only epithelial-lineage cells
#'   are gated; non-epithelial are forced Non_malignant.
#' @param correct_labels If TRUE (default) reconcile marker call with the gate
#'   in scONCO_L1_gated / scONCO_L4_gated.
#' @param verbose Logical.
#' @return Updated Seurat object with the cancer-layer meta.data columns.
apply_cnv_malignant_gate <- function(seurat_obj,
                                     cnv_score        = NULL,
                                     method           = "precomputed",
                                     cnv_col          = NULL,
                                     weights          = c(marker = 0.40, cnv = 0.45, module = 0.15),
                                     fusion           = c("weighted", "logistic"),
                                     threshold        = 0.5,
                                     uncertain_margin = 0.05,
                                     restrict_to_epithelial = TRUE,
                                     correct_labels   = TRUE,
                                     verbose          = TRUE) {
  fusion <- match.arg(fusion)
  if (!inherits(seurat_obj, "Seurat")) stop("seurat_obj must be a Seurat object")
  md <- seurat_obj@meta.data
  cells <- rownames(md)

  l1_col <- .sconco_l1_col(md); l4_col <- .sconco_l4_col(md)
  conf_col <- .sconco_conf_col(md)
  if (is.na(l1_col) || is.na(l4_col))
    stop("Could not find scONCO Level1/Level4 columns. Run run_scONCO() first.")

  # --- CNV score ---
  if (is.null(cnv_score))
    cnv_score <- compute_cnv_score(seurat_obj, method = method, cnv_col = cnv_col)
  cnv <- as.numeric(cnv_score[cells]); cnv[is.na(cnv)] <- 0.5   # neutral when missing

  # --- component signals ---
  marker <- .marker_malignant_prob(md, l1_col, conf_col)
  module <- .module_prior(md, l1_col, l4_col)

  # --- fuse ---
  w <- weights / sum(weights)
  lin <- w["marker"] * marker + w["cnv"] * cnv + w["module"] * module
  P <- if (fusion == "logistic") 1 / (1 + exp(-4 * (lin - 0.5))) else as.numeric(lin)

  l1 <- as.character(md[[l1_col]]); l4 <- as.character(md[[l4_col]])
  epi_lineage <- .is_epithelial_label(l1) | .is_malignant_label(l1) |
                 .is_epithelial_label(l4) | grepl("^Mal", l4)

  call <- ifelse(P >= threshold, "Malignant", "Non_malignant")
  if (uncertain_margin > 0)
    call[abs(P - threshold) < uncertain_margin] <- "Uncertain"
  if (restrict_to_epithelial) {
    call[!epi_lineage] <- "Non_malignant"
    P[!epi_lineage]    <- pmin(P[!epi_lineage], 0.05)
  }

  # --- reconcile labels ---
  l1_gated <- l1; l4_gated <- l4; corrected <- rep(FALSE, length(l1))
  if (correct_labels) {
    # demotion: marker said malignant, gate says non-malignant -> normal epithelial
    demote <- .is_malignant_label(l1) & call == "Non_malignant"
    l1_gated[demote] <- "Epithelial cell"; l4_gated[demote] <- "Epi"; corrected[demote] <- TRUE
    # promotion: marker said normal epithelial, gate says malignant -> malignant
    promote <- .is_epithelial_label(l1) & !.is_malignant_label(l1) & call == "Malignant"
    l1_gated[promote] <- "Malignant cell"
    l4_gated[promote] <- ifelse(grepl("^Mal", l4[promote]), l4[promote], "Mal")
    corrected[promote] <- TRUE
  }

  seurat_obj@meta.data$scONCO_cnv_score      <- cnv
  seurat_obj@meta.data$scONCO_malignant_prob <- P
  seurat_obj@meta.data$scONCO_malignant_call <- call
  seurat_obj@meta.data$scONCO_L1_gated       <- l1_gated
  seurat_obj@meta.data$scONCO_L4_gated       <- l4_gated
  seurat_obj@meta.data$scONCO_gate_corrected <- corrected

  if (verbose) {
    message(sprintf("[CNV gate] %d cells | Malignant=%d Non_malignant=%d Uncertain=%d | corrected=%d (%.1f%%)",
                    length(call), sum(call == "Malignant"), sum(call == "Non_malignant"),
                    sum(call == "Uncertain"), sum(corrected),
                    100 * mean(corrected)))
  }
  invisible(seurat_obj)
}

# ---------------------------------------------------------------------------
# 4) evaluation helper — malignant-calling accuracy vs ground-truth
# ---------------------------------------------------------------------------
#' Evaluate gate malignant calls against a ground-truth malignant flag
#' @param seurat_obj Gated Seurat object.
#' @param truth_col meta.data column with ground-truth cell type or a malignant
#'   flag; malignant is inferred by keyword if it is a cell-type column.
#' @return Named numeric c(accuracy, sensitivity, specificity, f1).
evaluate_malignant_gate <- function(seurat_obj, truth_col) {
  md <- seurat_obj@meta.data
  if (!truth_col %in% colnames(md)) stop("truth_col not found: ", truth_col)
  true_mal <- .is_malignant_label(as.character(md[[truth_col]]))
  pred_mal <- md$scONCO_malignant_call == "Malignant"
  tp <- sum(pred_mal & true_mal); tn <- sum(!pred_mal & !true_mal)
  fp <- sum(pred_mal & !true_mal); fn <- sum(!pred_mal & true_mal)
  sens <- if ((tp+fn) > 0) tp/(tp+fn) else NA_real_
  spec <- if ((tn+fp) > 0) tn/(tn+fp) else NA_real_
  prec <- if ((tp+fp) > 0) tp/(tp+fp) else NA_real_
  f1   <- if (!is.na(prec) && !is.na(sens) && (prec+sens) > 0) 2*prec*sens/(prec+sens) else NA_real_
  c(accuracy = (tp+tn)/nrow(md), sensitivity = sens, specificity = spec, f1 = f1)
}
