# Generated from code/R/14_FUN_sconco_cnv.R by tools/sync_package.py - do not edit here:
# edit the source and re-run the script.
###############################################################################
# 14_FUN_sconco_cnv.R - scONCO's own copy-number signal and sample-level tumour
# test (components C1 and C2 of docs/METHOD_DEVELOPMENT_PLAN.md).
#
# No external caller. The scheme is the one Tirosh et al. 2016 introduced and
# inferCNV/CopyKAT/SCEVAN all build on - expression smoothed along the genome
# relative to a diploid reference - implemented here from scratch so that the
# whole pipeline, reference selection included, is scONCO's.
#
#   sconco_cnv_profile(counts, ref_cells, gene_pos, ...)
#       -> smoothed per-cell, per-window CNV matrix (cells x windows)
#   sconco_cnv_score(profile, ref_cells)
#       -> per cell: burden, coherence
#   sconco_tumour_test(score, ref_cells, epi_cells, ...)
#       -> sample-level verdict + null threshold (C2)
#   sconco_cnv_call(score, ref_cells, epi_cells, verdict, ...)
#       -> per cell: aneuploid / diploid / not.defined + P(aneuploid)
#   sconco_cnv(counts, ref_cells, epi_cells, gene_pos, ...)
#       -> runs all of the above
#
# Design choices, fixed before any result was seen (see the plan):
#   - genes ordered by genomic position; only genes expressed in >= min_frac of cells
#   - log2(1 + CP10k); per-gene centring on the REFERENCE mean; clipping at +/- clip
#   - moving average of `window` genes within each chromosome (no smoothing across
#     chromosome boundaries); a chromosome with fewer than `window` genes is skipped
#   - per-cell re-centring on the cell's own median window (removes library effects)
#   - burden    = mean of squared smoothed values
#   - coherence = Pearson correlation of the cell's profile with the mean profile of
#                 the top 5% burden cells among the CANDIDATE (epithelial) cells;
#                 when there is no tumour, that "top 5%" is noise and coherence is ~0
#   - null: the reference cells' burden; a candidate cell is aneuploid when its
#           burden exceeds the null's `null_q` quantile AND its coherence exceeds
#           `coh_min`; diploid when within the null; above the null but incoherent
#           -> not.defined
#   - sample-level test: the fraction of candidate cells above the null quantile;
#           tumour present when that fraction >= `frac_min` (with >= `n_min` cells)
###############################################################################

# Revision, change A2a: a second family of genomically clustered lineage
# programmes, excluded like the immune receptor loci. The epidermal
# differentiation complex on 1q21.3 (S100A, SPRR, LCE, IVL, FLG, LOR, TCHH, HRNR,
# RPTN, CRNN, CRCT1, PGLYRP3/4) and the keratin clusters on 12q13.13 (type II) and
# 17q21.2 (type I, with the KRTAP genes) are expressed at extreme levels by
# differentiating keratinocytes; in GSE144236 they dominated the burden of normal
# keratinocytes (P1: 12 of the windows that most separated normal from tumour
# keratinocytes were on 1q21) and raised two false alarms on normal skin.
SCONCO_EXCLUDE_KERATINOCYTE <- "^KRT[0-9]+[A-Z]?$|^KRTAP|^S100A[0-9]+[A-Z]?$|^SPRR|^LCE[0-9]|^IVL$|^FLG2?$|^LOR$|^TCHH(L1)?$|^HRNR$|^RPTN$|^CRNN$|^CRCT1$|^PGLYRP[34]$"

sconco_cnv_profile <- function(counts, ref_cells, gene_pos, window = 100L, min_frac = 0.05,
                               clip = 3, exclude_genes = NULL, ref_groups = NULL,
                               exclude_keratinocyte = FALSE, verbose = TRUE) {
  # ref_groups: optional named vector (cell -> group) for the reference cells. When
  # given, each reference cell is centred on ITS OWN group's mean, so the null
  # burden measures within-type noise only and is not inflated by the expression
  # differences between, say, T cells and macrophages. Candidates are always
  # centred on the pooled reference mean.
  stopifnot(is.character(ref_cells), length(ref_cells) >= 20)
  if (!inherits(counts, "dgCMatrix")) counts <- methods::as(counts, "CsparseMatrix")
  ref_cells <- intersect(ref_cells, colnames(counts))
  # genes with positions, expressed in enough cells
  gp <- gene_pos[gene_pos$gene %in% rownames(counts), ]
  # immunoglobulin, T-cell-receptor and HLA loci are genomically clustered and
  # expressed at extreme levels by plasma cells, T cells and antigen-presenting
  # cells: in a smoothed profile they look like focal copy-number events and
  # put the reference null out of reach (CID44991: 1,513 plasma cells in the
  # reference, 4,018 tumour cells declared "no tumour"). Excluded, as inferCNV
  # and CopyKAT do.
  gp <- gp[!grepl("^IG[HKL][VDJCAGME]|^TR[ABGD][VDJC]|^HLA-", gp$gene), ]
  if (isTRUE(exclude_keratinocyte)) gp <- gp[!grepl(SCONCO_EXCLUDE_KERATINOCYTE, gp$gene), ]   # revision, A2a
  if (!is.null(exclude_genes)) gp <- gp[!gp$gene %in% exclude_genes, ]
  keep <- Matrix::rowSums(counts[gp$gene, , drop = FALSE] > 0) >= min_frac * ncol(counts)
  gp <- gp[keep, ]
  gp <- gp[order(factor(gp$chrom, levels = c(1:22, "X")), gp$start), ]
  # log CP10k
  cs <- pmax(Matrix::colSums(counts), 1)
  X <- counts[gp$gene, , drop = FALSE] %*% Matrix::Diagonal(x = 1e4 / cs)
  X <- methods::as(X, "CsparseMatrix"); X@x <- log2(1 + X@x)
  dimnames(X) <- list(gp$gene, colnames(counts))      # the product drops names
  X <- as.matrix(X)                                   # genes x cells (dense; genes ~ 5-8k)
  # centre on the reference mean (reference cells on their own group's mean), clip
  mu <- rowMeans(X[, ref_cells, drop = FALSE])
  X <- X - mu
  if (!is.null(ref_groups)) {
    rg <- ref_groups[ref_cells]; rg[is.na(rg)] <- "ref"
    for (g in unique(rg)) { cg <- ref_cells[rg == g]; if (length(cg) >= 10)
      X[, cg] <- X[, cg, drop = FALSE] - rowMeans(X[, cg, drop = FALSE]) }
  }
  X[X >  clip] <-  clip; X[X < -clip] <- -clip
  # moving average within chromosome
  prof <- list(); wpos <- list()
  for (ch in unique(gp$chrom)) {
    idx <- which(gp$chrom == ch); n <- length(idx)
    if (n < window) next
    cx <- apply(X[idx, , drop = FALSE], 2, cumsum)   # genes x cells cumulative
    cx <- rbind(0, cx)
    sm <- (cx[(window + 1):(n + 1), , drop = FALSE] - cx[1:(n - window + 1), , drop = FALSE]) / window
    # step of window/4 to keep windows ~independent and the matrix small
    step <- max(1L, window %/% 4L); pick <- seq(1, nrow(sm), by = step)
    prof[[ch]] <- t(sm[pick, , drop = FALSE])         # cells x windows
    wpos[[ch]] <- data.frame(chrom = ch, start = gp$start[idx[pick]], end = gp$start[idx[pmin(pick + window - 1, n)]])
  }
  P <- do.call(cbind, prof); W <- do.call(rbind, wpos)
  # per-cell re-centring (each cell's own median window), then re-centre on the reference
  P <- P - apply(P, 1, stats::median)
  P <- sweep(P, 2, colMeans(P[ref_cells, , drop = FALSE]))
  if (verbose) message(sprintf("  cnv profile: %d genes, %d windows over %d chromosomes, %d reference cells",
                               nrow(gp), ncol(P), length(prof), length(ref_cells)))
  list(profile = P, windows = W, genes = gp$gene, ref_cells = ref_cells)
}

# Lineage correction with a lineage-matched anchor. Epithelial cells differ from
# an immune/stromal reference in expression programmes that are themselves
# clustered on the genome (the keratin clusters on 12q13 and 17q21 are the
# obvious case), so every epithelial cell - normal or malignant - carries a
# shared, CNV-like offset. Estimating that offset from "the least aberrant
# candidates" fails on tumours with subtle copy-number signal, because those
# cells ARE malignant (measured: 14/19 -> 9/19 tumours detected). Instead the
# offset is estimated from ANCHOR cells: candidates the marker step assigned to
# a terminally differentiated normal epithelial type that is not a carcinoma
# cell of origin (ciliated, goblet, islet, Paneth, tuft, enteroendocrine).
# Those cells are diploid with high probability, lineage-matched, and only
# exist because database v3.2 names them. With >= 20 anchors their median
# profile is subtracted from every candidate and their own (leave-one-out)
# burden becomes the null; with fewer than 20 nothing is subtracted and the
# immune/stromal null is used.
SCONCO_ANCHOR_TYPES <- c("Epi_Cil", "Epi_Gob", "Islet", "Epi_Paneth", "Epi_Tuft", "Epi_EEC")

sconco_lineage_correct <- function(prof, candidate_cells, anchor_cells, min_n = 10L, top_frac = 0.05, coh_max = 0.2) {
  P <- prof$profile; cand <- intersect(candidate_cells, rownames(P)); anc <- intersect(anchor_cells, cand)
  prof$anchor_cells <- anc; prof$lineage_corrected <- FALSE; prof$anchors_dropped <- 0L
  # guard: an anchor is a marker call, and in a sample that is almost entirely
  # tumour the "goblet cells" are tumour cells (CID4290A: 51 of 54 anchors were
  # malignant, the correction subtracted the tumour's own profile and 4,025
  # tumour cells were declared normal). Before correcting, the putative clone
  # template is taken from the top `top_frac` burden candidates; an anchor whose
  # uncorrected profile correlates with that template above `coh_max` is not
  # diploid-looking and is dropped. When there is no tumour the template is
  # noise and genuine anchors correlate with it at ~0.
  if (length(anc) >= min_n && length(cand) > length(anc)) {
    b <- rowMeans(P[cand, , drop = FALSE]^2); k <- max(5L, ceiling(top_frac * length(cand)))
    tmpl <- colMeans(P[cand[order(b, decreasing = TRUE)][seq_len(k)], , drop = FALSE])
    coh <- as.numeric(stats::cor(t(P[anc, , drop = FALSE]), tmpl))
    prof$anchors_dropped <- sum(coh > coh_max); anc <- anc[coh <= coh_max]; prof$anchor_cells <- anc
  }
  if (length(anc) < min_n) return(prof)
  L <- apply(P[anc, , drop = FALSE], 2, stats::median)
  Pc <- sweep(P[cand, , drop = FALSE], 2, L)
  # anchors are corrected leave-one-out (each on the median of the OTHER
  # anchors): an anchor corrected on a median it helped define has a deflated
  # burden, and with 20-40 anchors that deflation pushed the null below the
  # normal epithelial cells it is meant to cover (KIM_P20)
  A <- P[anc, , drop = FALSE]
  for (i in seq_along(anc)) Pc[anc[i], ] <- A[i, ] - apply(A[-i, , drop = FALSE], 2, stats::median)
  P[cand, ] <- Pc
  prof$profile <- P; prof$lineage_corrected <- TRUE; prof
}

sconco_cnv_score <- function(prof, candidate_cells, top_frac = 0.05, lowfreq = 8L) {
  P <- prof$profile
  # burden on the low-frequency component: the 100-gene windows are re-smoothed
  # over `lowfreq` consecutive windows within each chromosome (8 windows at a
  # step of 25 genes ~ 275 genes ~ an arm-scale event). A copy-number event
  # keeps its amplitude; the focal expression offsets that make a normal
  # epithelial cell differ from an immune reference (keratin clusters and the
  # like, one to three windows wide) are diluted. lowfreq = 1 gives the plain
  # burden.
  Pl <- P
  if (lowfreq > 1L && !is.null(prof$windows)) {
    for (ch in unique(prof$windows$chrom)) {
      j <- which(prof$windows$chrom == ch); if (length(j) < 2) next
      w <- min(lowfreq, length(j)); cs <- cbind(0, t(apply(P[, j, drop = FALSE], 1, cumsum)))
      k <- length(j) - w + 1L
      sm <- (cs[, (w + 1):(w + k), drop = FALSE] - cs[, 1:k, drop = FALSE]) / w
      # pad to the chromosome's width so windows keep their columns (edge windows repeat the end values)
      Pl[, j] <- sm[, pmin(pmax(seq_along(j) - (w %/% 2), 1L), k), drop = FALSE]
    }
  }
  burden <- rowMeans(Pl^2)
  cand <- intersect(candidate_cells, rownames(P))
  if (length(cand) < 10) cand <- rownames(P)
  k <- max(5L, ceiling(top_frac * length(cand)))
  top <- cand[order(burden[cand], decreasing = TRUE)][seq_len(k)]
  template <- colMeans(P[top, , drop = FALSE])
  coherence <- as.numeric(stats::cor(t(P), template))
  names(coherence) <- rownames(P)
  data.frame(cell = rownames(P), burden = burden, coherence = coherence, row.names = rownames(P))
}

sconco_tumour_test <- function(score, ref_cells, candidate_cells, null_q = 0.99, frac_min = 0.15,
                               n_min = 30L, null_cells = NULL, null_min_n = 10L, n_empirical = 30L) {
  # null_cells: the cells whose burden defines the null (anchors when the lineage
  # correction was applied, else the immune/stromal reference).
  # Anchor null: the empirical null_q quantile of the anchors' leave-one-out
  # burden (sconco_lineage_correct). Two alternatives were tried and rejected:
  # requiring >= 30 anchors (KUL_KUL21 fell to the reference null, MR 0.97 ->
  # 0.33) and a log-normal floor on the anchor burden (the deflated values make
  # the log-scale MAD huge; KIM_P20 / KIM_P34 / KUL_KUL21 all became "absent").
  use_anchor <- length(null_cells) >= null_min_n
  nul <- intersect(if (use_anchor) null_cells else ref_cells, score$cell)
  cand <- setdiff(intersect(candidate_cells, score$cell), if (use_anchor) nul else character(0))
  b_nul <- score[nul, "burden"]
  thr <- stats::quantile(b_nul, null_q, names = FALSE)
  # small anchor sets (10-29): the empirical quantile is the maximum of a few
  # values; the anchors are one normal type corrected leave-one-out, so their
  # burden is taken as log-normal and the threshold is the larger of the two.
  # (Never applied to the mixed-type reference, where it over-shoots.)
  if (use_anchor && length(nul) < n_empirical) { lb <- log(pmax(b_nul, 1e-8)); thr <- max(thr, exp(stats::median(lb) + stats::qnorm(null_q) * stats::mad(lb))) }
  med <- stats::median(b_nul)
  frac <- if (length(cand)) mean(score[cand, "burden"] > thr) else NA_real_
  present <- length(cand) >= n_min && !is.na(frac) && frac >= frac_min
  list(present = present, frac_above_null = frac, n_candidates = length(cand), n_null = length(nul),
       null_source = if (use_anchor) "anchor" else "reference",
       null_threshold = thr, null_median = med, null_q = null_q, frac_min = frac_min)
}

# Revision, change A2b: the amplitude of a cell group's copy-number profile
# against its own genome-wide baseline. The profile is centred on the reference,
# so a lineage programme shifts many windows a little and a copy-number event
# shifts one arm a lot: the statistic is the largest deviation of a chromosome
# arm (>= min_windows windows; arms split at the median window of each
# chromosome, no centromere table) from the median of all arm means, in the
# units of the profile (log2 CP10k, 100-gene windows). In GSE144236 P5 normal
# skin a cluster of differentiated keratinocytes was called aneuploid on a
# whole-X shift of 0.05; a single-copy arm loss is 0.15-0.5 in these units.
.sconco_arm_of <- function(windows) { chrom <- as.character(windows$chrom)
  paste0(chrom, ifelse(stats::ave(seq_along(chrom), chrom, FUN = function(i) i <= stats::median(i)) == 1, "p", "q")) }
sconco_amplitude <- function(profile, cells, windows, min_windows = 4L) {
  cells <- intersect(cells, rownames(profile)); if (length(cells) < 3 || is.null(windows)) return(NA_real_)
  d <- colMeans(profile[cells, , drop = FALSE]); arm <- .sconco_arm_of(windows)
  am <- tapply(d, arm, mean); nw <- tapply(d, arm, length); am <- am[nw >= min_windows]
  if (!length(am)) return(NA_real_); max(abs(am - stats::median(am)))
}

sconco_cnv_call <- function(score, ref_cells, candidate_cells, test, coh_min = 0.2, profile = NULL,
                            cluster_level = TRUE, max_k = 6L, min_cluster = 5L, secondary_cells = character(0), coh_sec = 0.5,
                            windows = NULL, amp_min = NULL) {
  # amp_min (revision, A2b): a cluster above the null is called aneuploid only if its
  # amplitude (sconco_amplitude) reaches amp_min; when no primary cluster does,
  # the sample is not tumour-present and every candidate is diploid. NULL = off.
  cand <- intersect(candidate_cells, score$cell)
  call <- setNames(rep("not.defined", nrow(score)), score$cell)
  p <- setNames(rep(NA_real_, nrow(score)), score$cell)
  if (!test$present) {                                # no tumour: every candidate is diploid
    call[cand] <- "diploid"; p[cand] <- 0
    return(data.frame(cell = score$cell, cnv = call, p_aneuploid = p, row.names = score$cell))
  }
  # secondary candidates (stromal / endothelial marker groups): never part of
  # the null, never in the tumour test; called aneuploid only when a cluster of
  # them is BOTH above the null and coherent with the primary candidates' clone
  # template - the mesenchymal-like melanoma that the marker step files under
  # "stromal" (GSE72056 T81). A fibroblast or endothelial cluster with its own
  # lineage offset is above the null and correlates with the template at
  # 0.2-0.35 (the shared "not immune" component); tumour cells filed as
  # stromal correlate at ~0.9 (T81: median 0.88). Hence the stricter `coh_sec`.
  sec <- setdiff(intersect(secondary_cells, score$cell), cand)
  if (length(sec) >= min_cluster && !is.null(profile)) {
    bs <- score[sec, "burden"]; cs <- score[sec, "coherence"]
    k <- max(2L, min(max_k, 2L + length(sec) %/% 300L))
    grp <- if (length(sec) >= 2 * min_cluster) stats::cutree(stats::hclust(stats::dist(profile[sec, , drop = FALSE]), method = "ward.D2"), k = k) else rep(1L, length(sec))
    for (g in unique(grp)) { i <- which(grp == g); if (length(i) < min_cluster) next
      if (stats::median(bs[i]) > test$null_threshold && stats::median(cs[i]) >= coh_sec) { call[sec[i]] <- "aneuploid"; p[sec[i]] <- mean(bs[i] > test$null_threshold) } }
  }
  b <- score[cand, "burden"]; ch <- score[cand, "coherence"]
  above <- b > test$null_threshold
  # per-cell rule: above the null and coherent -> aneuploid; within the null ->
  # diploid; above the null but incoherent (isolated noisy cell) -> not.defined
  cl <- ifelse(above, ifelse(ch >= coh_min, "aneuploid", "not.defined"), "diploid")
  z <- (b - test$null_threshold) / (test$null_threshold - test$null_median + 1e-9)
  pa <- stats::plogis(2 * z); pa[ch < coh_min] <- pmin(pa[ch < coh_min], 0.5)
  # cluster-level rule: a single cell's burden is noisy, a copy-number clone is
  # not. Candidates are clustered on their smoothed profiles and each cluster
  # is called as a unit - aneuploid when its median burden clears the null,
  # diploid when its median is within the null - and a cell inherits its
  # cluster's call with a probability equal to the cluster's share of cells
  # above the null. This is the clone-level consensus that CopyKAT and SCEVAN
  # also rely on; here it is done without an external caller.
  if (cluster_level && !is.null(profile) && length(cand) >= 2 * min_cluster) {
    Pc <- profile[cand, , drop = FALSE]
    k <- max(2L, min(max_k, 2L + length(cand) %/% 300L))
    hc <- stats::hclust(stats::dist(Pc), method = "ward.D2"); grp <- stats::cutree(hc, k = k)
    amp_by_cluster <- c(); an_clusters <- list()
    for (g in unique(grp)) {
      i <- which(grp == g); if (length(i) < min_cluster) next
      share <- mean(above[i]); medb <- stats::median(b[i])
      if (medb > test$null_threshold) { cl[i] <- "aneuploid"; pa[i] <- pmax(pa[i], share); an_clusters[[as.character(g)]] <- i
        if (!is.null(amp_min) && !is.null(windows)) amp_by_cluster[as.character(g)] <- sconco_amplitude(profile, cand[i], windows) }
      else if (share < 0.25) { cl[i] <- "diploid"; pa[i] <- pmin(pa[i], share) }
      # clusters with a median inside the null but many cells above it keep their per-cell calls
    }
    # revision (A2b): the amplitude floor is applied to the UNION of the clusters above
    # the null - the clone as a whole - not cluster by cluster, because ward
    # sub-clusters of one clone split its events between them and a tumour whose
    # union sits just above the floor (KIM_P28, 0.133) would lose sub-clusters.
    # A programme cluster (GSE144236 P5: 0.099) is the whole union of its sample.
    # The union is every primary candidate called aneuploid at this point, whether
    # by its cluster or by the per-cell rule (a tumour whose clusters straddle the
    # null threshold is called cell by cell: KUL01), so the gate sees the whole
    # aneuploid population; with fewer than 3 such cells there is nothing to judge.
    amp_union <- NA_real_; passed <- NA_integer_
    if (!is.null(amp_min) && !is.null(windows)) {
      iu <- which(cl == "aneuploid")
      if (length(iu) >= 3) { amp_union <- sconco_amplitude(profile, cand[iu], windows)
        if (is.finite(amp_union)) { passed <- as.integer(amp_union >= amp_min)
          if (passed == 0L) {   # above the null, below the floor: a programme, not a clone -> no tumour
            cl[] <- "diploid"; pa <- pmin(pa, 0.5); call[sec] <- "not.defined"; p[sec] <- NA_real_ } } }
    }
    attr(call, "amplitude") <- amp_by_cluster; attr(call, "amplitude_union") <- amp_union; attr(call, "amplitude_passed") <- passed
  }
  call[cand] <- cl; p[cand] <- pa
  out <- data.frame(cell = score$cell, cnv = call, p_aneuploid = p, row.names = score$cell)
  attr(out, "amplitude") <- attr(call, "amplitude"); attr(out, "amplitude_union") <- attr(call, "amplitude_union"); attr(out, "amplitude_passed") <- attr(call, "amplitude_passed"); out
}

#' Copy-number signal and sample-level tumour test
#'
#' Expression smoothed along the genome relative to a diploid reference (the scheme
#' of Tirosh et al. 2016, implemented here from scratch), a per-cell burden and
#' coherence, a sample-level test (the share of candidate cells above the reference
#' null; tumour present at `frac_min` or more) and a per-cell call that is made only
#' when the test is positive.
#'
#' @param counts Raw counts, genes x cells.
#' @param ref_cells Diploid reference cells (immune / stromal), at least 20.
#' @param candidate_cells Cells to be called (the epithelial lineage).
#' @param gene_pos Gene positions, see [sconco_load_gene_positions()].
#' @param anchor_cells Normal cells of the candidates' lineage used to correct the
#'   lineage programme, when enough are found.
#' @param ref_groups Optional group of each reference cell (each centred on its own
#'   group).
#' @param window Genes per moving-average window.
#' @param min_frac Minimum fraction of cells expressing a gene.
#' @param clip Clipping of centred expression.
#' @param null_q Quantile of the reference burden that defines the null.
#' @param frac_min Share of candidates above the null for a tumour to be present.
#' @param coh_min Minimum coherence with the clone for an aneuploid call.
#' @param lowfreq Low-frequency components removed from the profile.
#' @param secondary_cells Stromal / endothelial cells called only as coherent clusters.
#' @param exclude_keratinocyte Drop the epidermal-differentiation and keratin loci
#'   (revision option A2a; off).
#' @param amp_min Copy-number amplitude floor (revision option A2b; `NULL` = off).
#' @param verbose Print progress.
#' @return A list with `score`, `test`, `call`, `windows` and `profile`.
#' @export
sconco_cnv <- function(counts, ref_cells, candidate_cells, gene_pos, anchor_cells = character(0),
                       ref_groups = NULL, window = 100L, min_frac = 0.05, clip = 3, null_q = 0.99,
                       frac_min = 0.15, coh_min = 0.2, lowfreq = 8L, secondary_cells = character(0),
                       exclude_keratinocyte = FALSE, amp_min = NULL, verbose = TRUE) {
  prof <- sconco_cnv_profile(counts, ref_cells, gene_pos, window, min_frac, clip, ref_groups = ref_groups,
                             exclude_keratinocyte = exclude_keratinocyte, verbose = verbose)
  prof <- sconco_lineage_correct(prof, candidate_cells, anchor_cells)
  sc   <- sconco_cnv_score(prof, candidate_cells, lowfreq = lowfreq)
  tt   <- sconco_tumour_test(sc, ref_cells, candidate_cells, null_q, frac_min,
                             null_cells = if (isTRUE(prof$lineage_corrected)) prof$anchor_cells else NULL)
  cl   <- sconco_cnv_call(sc, ref_cells, candidate_cells, tt, coh_min, profile = prof$profile, secondary_cells = secondary_cells,
                          windows = prof$windows, amp_min = amp_min)
  if (!is.null(amp_min)) {   # revision (A2b): the verdict carries the amplitude, and "present" requires a cluster above the floor
    tt$amplitude_max <- attr(cl, "amplitude_union"); tt$amplitude_clusters <- attr(cl, "amplitude")
    tt$amplitude_passed <- attr(cl, "amplitude_passed"); tt$amp_min <- amp_min
    if (isTRUE(tt$present) && isTRUE(tt$amplitude_passed == 0L)) { tt$present <- FALSE; tt$below_amplitude_floor <- TRUE } }
  if (verbose) message(sprintf("  lineage anchor: %d cells (%d dropped as clone-like) -> correction %s; null from %s",
                               length(prof$anchor_cells), prof$anchors_dropped, if (isTRUE(prof$lineage_corrected)) "applied" else "skipped", tt$null_source))
  if (verbose) message(sprintf("  tumour test: %s (%.1f%% of %d candidates above the null 99th percentile); calls: %s",
                               if (tt$present) "PRESENT" else "absent", 100 * tt$frac_above_null, tt$n_candidates,
                               paste(names(table(cl$cnv)), table(cl$cnv), collapse = " ")))
  list(score = sc, test = tt, call = cl, windows = prof$windows, profile = prof$profile)
}
