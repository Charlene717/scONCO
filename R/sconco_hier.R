# Generated from code/R/15_FUN_sconco_hier.R by tools/sync_package.py - do not edit here:
# edit the source and re-run the script.
###############################################################################
# 15_FUN_sconco_hier.R - multi-resolution hierarchical annotation (component C3
# of docs/METHOD_DEVELOPMENT_PLAN.md), with the epithelial node decided by the
# copy-number signal (C1/C2) and evidence-weighted marker scores (C4).
#
# Scheme (after scMRMA, Li et al. 2022, re-implemented for our hierarchy):
#   Level 1  cells are clustered at low resolution; each cluster takes the
#            Level-1 node whose best member type has the highest mean per-cell
#            probability (best member, not the sum over members).
#            The Malignant and normal-Epithelial compartments of the database
#            form ONE node here, "Epithelial-lineage": markers are not asked to
#            separate them.
#   Level 2+ within each Level-1 group the cells are re-clustered at a higher
#            resolution and each cluster takes the best CHILD of its parent;
#            the same again for Level 3.
#   Level 4  per cell, among the children of its Level-3 node, with a reject
#            option (unassigned / ambiguous) as in sconco_annotate(), thresholds
#            scaled to chance level 1/k for nodes with many children.
#   Epithelial-lineage node: Level 2 is Malignant vs Normal epithelium and is
#            decided by sconco_cnv() (C1) after the sample-level tumour test
#            (C2), with the immune/stromal/endothelial Level-1 groups as the
#            diploid reference and the anchor types as lineage anchors. Cells
#            called malignant continue to the cancer-type modules; cells called
#            normal continue to the normal epithelial types. Cells of the
#            "Other" node and cells whose best type is a malignant module are
#            copy-number candidates too (non-epithelial tumours: melanoma,
#            glioma); called aneuploid they get Level 2 = Malignant under
#            their marker Level 1.
#
# Scoring (C4): marker weight w(g,t) = spec(g) * log2(1 + n_sources(g,t)) when the
# database carries n_sources, else spec(g) alone (`weight_evidence = FALSE` turns
# the source term off for the ablation).
#
# Dependencies: Seurat (clustering), Matrix. sconco_cnv() from 14_FUN_sconco_cnv.R.
###############################################################################

.sconco_hier_nodes <- function(db) {
  h <- unique(db[, c("Level1", "Level2", "Level3", "Level4_Abb")])
  h$node1 <- ifelse(h$Level1 %in% c("Malignant cell", "Epithelial cell"), "Epithelial-lineage", h$Level1)
  h$node2 <- ifelse(h$node1 == "Epithelial-lineage",
                    ifelse(h$Level1 == "Malignant cell", "Malignant", "Normal epithelium"), h$Level2)
  h
}

# per-cell weighted marker scores -> softmax probabilities over Level-4 types
.sconco_type_probs <- function(X, db, weight_evidence = TRUE) {
  types <- unique(db$Level4_Abb)
  n_types_using <- table(db$gene)
  S <- matrix(0, ncol(X), length(types), dimnames = list(colnames(X), types))
  for (ty in types) {
    d <- db[db$Level4_Abb == ty & db$gene %in% rownames(X), ]
    if (!nrow(d)) { S[, ty] <- NA; next }
    spec <- 1 - as.numeric(n_types_using[d$gene]) / length(types)
    w <- if (weight_evidence && "n_sources" %in% names(d)) spec * log2(1 + d$n_sources) else spec
    w <- w / sum(w)
    S[, ty] <- as.numeric(t(X[d$gene, , drop = FALSE]) %*% w)
  }
  S <- S[, !is.na(S[1, ]), drop = FALSE]
  Z <- S - apply(S, 1, max); P <- exp(Z); P <- P / rowSums(P)
  list(score = S, prob = P)
}

.sconco_cluster <- function(seu, cells, resolution, min_cells = 50L, npcs = 30L, seed = 1L) {
  if (length(cells) < min_cells) return(setNames(rep("c0", length(cells)), cells))
  s <- suppressWarnings(subset(seu, cells = cells))
  s <- suppressWarnings(Seurat::FindVariableFeatures(s, nfeatures = 2000, verbose = FALSE))
  s <- Seurat::ScaleData(s, features = Seurat::VariableFeatures(s), verbose = FALSE)
  np <- min(npcs, length(cells) - 1L, length(Seurat::VariableFeatures(s)) - 1L)
  s <- suppressWarnings(Seurat::RunPCA(s, npcs = np, verbose = FALSE, seed.use = seed))
  s <- Seurat::FindNeighbors(s, dims = seq_len(np), verbose = FALSE, k.param = min(20L, length(cells) - 1L))
  s <- Seurat::FindClusters(s, resolution = resolution, verbose = FALSE, random.seed = seed)
  setNames(paste0("c", as.character(Seurat::Idents(s))), cells)
}

.assign_clusters <- function(clu, P, candidates, mixed_frac = 0.25) {
  # clu: named cluster vector; P: prob matrix (cells x types); candidates: list(node -> types)
  # Each cluster takes the node whose best member type has the highest mean
  # per-cell probability (best member, not the sum: summing over a node's types
  # rewards the node with the most types, which is the compartment bias the
  # benchmark paper diagnosed). A cluster in which at least `mixed_frac` of the
  # cells individually prefer a different node is treated as mixed and its
  # cells are assigned one by one - the case of NK cells co-clustering with T.
  nodes <- names(candidates)
  node_of_cell <- function(cells) {
    sc <- sapply(nodes, function(nd) { ty <- intersect(candidates[[nd]], colnames(P))
      if (!length(ty)) rep(-Inf, length(cells)) else apply(P[cells, ty, drop = FALSE], 1, max) })
    sc <- matrix(sc, nrow = length(cells), dimnames = list(cells, nodes))
    list(best = nodes[max.col(sc, "first")], score = sc) }
  out <- setNames(rep(NA_character_, length(clu)), names(clu))
  for (k in unique(clu)) {
    cells <- names(clu)[clu == k]; nc <- node_of_cell(cells)
    cl_node <- nodes[which.max(colMeans(nc$score))]
    if (length(nodes) > 1 && mean(nc$best != cl_node) >= mixed_frac) out[cells] <- nc$best else out[cells] <- cl_node
  }
  out
}

#' Hierarchical annotation with the copy-number split of the epithelial lineage
#'
#' The engine behind [sconco()]; most users call that instead. Cells are clustered at
#' increasing resolution and each cluster takes the best child of its parent node
#' (after scMRMA); Level 4 is decided per cell with a reject option. Malignant and
#' normal epithelium form one node for the markers and are split by [sconco_cnv()].
#'
#' @param counts,db,gene_pos,resolutions,weight_evidence,cnv,cnv_args,normalized,reference_pool,pool_min_own,seed,verbose
#'   As in [sconco()]; here `db` is required.
#' @param reject_p,reject_gap Level-4 reject thresholds (top probability, gap to the
#'   second), scaled to chance level for nodes with many children.
#' @param group_reference Centre each reference cell on its own group (off; see the
#'   source for why).
#' @param mixed_frac Share of a cluster that must disagree with its label before the
#'   cluster is split per cell.
#' @param keep_cnv Keep the copy-number profile and calls in the result (diagnostics).
#' @return A list with `annotation`, `tumour_test`, `probs`, `clusters_L1` (and `cnv`
#'   when `keep_cnv = TRUE`).
#' @export
sconco_hier <- function(counts, db, gene_pos = NULL, resolutions = c(0.2, 0.5, 1.2),
                        weight_evidence = TRUE, reject_p = 0.3, reject_gap = 0.1,
                        cnv = TRUE, cnv_args = list(), group_reference = FALSE, mixed_frac = 0.25, normalized = FALSE,
                        reference_pool = NULL, pool_min_own = 20L, keep_cnv = FALSE, seed = 1L, verbose = TRUE) {
  # reference_pool (revision, change A1): a counts matrix (genes x cells) of
  # immune-node cells from OTHER samples of the same study and platform, as
  # built by sconco_reference_pool(). Used only when this sample's own immune
  # reference has fewer than pool_min_own cells - the initial version's floor below which no
  # copy-number step runs (GSE131928: 15 of 25 CD45-depleted tumours had no
  # reference at all). Giving the pool to samples with 20-49 own immune cells
  # made the null wider and lost tumours (revision log: MGH101 MR 0.99 -> 0.17), so
  # the pool fills a reference that would otherwise not exist, nothing more.
  # The pool cells are appended to the reference for the copy-number step
  # only; each is centred on the mean of ITS OWN sample of origin (ref_groups),
  # so the null measures within-sample noise and not the offsets between
  # samples; the candidates are centred on the pooled reference mean. Pool
  # cells are never annotated, never candidates, and any pool cell that is
  # also in `counts` is dropped.
  # group_reference: centre each reference cell on its own Level-3 group (a null of
  # within-type noise only). Off by default: that null is tighter than what a
  # normal epithelial cell - which belongs to none of the reference groups - can
  # meet, and it flipped whole normal-epithelium clusters to aneuploid (KIM_P20,
  # NFM 0.003 -> 0.99) and a tumour to diploid (KUL_KUL01). The pooled null keeps
  # the between-type differences and is the conservative one.
  stopifnot(requireNamespace("Seurat", quietly = TRUE))
  if (!inherits(counts, "dgCMatrix")) counts <- methods::as(methods::as(counts, "CsparseMatrix"), "dgCMatrix")
  # Revision: "state" modules (db$part == "state"; v3.3: the glioblastoma
  # AC-/OPC-/NPC-/MES-like states) are malignant sub-state programmes whose genes
  # are those of the normal lineage each state mimics. They take no part in any
  # node decision or in copy-number candidacy - the lineage probabilities P are
  # computed WITHOUT them, exactly as in the initial version - and are used only at Level 4,
  # inside a Level-3 node of the Malignant branch, to name cells the copy-number
  # step has already declared malignant.
  state_types <- if ("part" %in% names(db)) unique(db$Level4_Abb[db$part == "state"]) else character(0)
  db_all <- db[db$gene %in% rownames(counts), ]
  db <- db_all[!db_all$Level4_Abb %in% state_types, ]
  H <- .sconco_hier_nodes(db)
  # normalised matrix on the marker genes for scoring
  seu <- suppressWarnings(Seurat::CreateSeuratObject(counts, min.cells = 0, min.features = 0))
  seu <- if (normalized) Seurat::SetAssayData(seu, layer = "data", new.data = counts) else Seurat::NormalizeData(seu, verbose = FALSE)
  X <- as.matrix(Seurat::GetAssayData(seu, layer = "data")[unique(db$gene), , drop = FALSE])
  pr <- .sconco_type_probs(X, db, weight_evidence); P <- pr$prob
  # state scores, on the same normalised matrix, only if the database has any
  if (length(state_types)) { Xs <- as.matrix(Seurat::GetAssayData(seu, layer = "data")[unique(db_all$gene), , drop = FALSE])
    S_state <- .sconco_type_probs(Xs, db_all, weight_evidence)$score; Hs <- .sconco_hier_nodes(db_all) } else S_state <- NULL
  cells <- colnames(counts)
  L1 <- L2 <- L3 <- L4 <- setNames(rep(NA_character_, length(cells)), cells)
  # ---- Level 1 -----------------------------------------------------------------
  c1 <- .sconco_cluster(seu, cells, resolutions[1], seed = seed)
  cand1 <- split(H$Level4_Abb, H$node1)
  L1[] <- .assign_clusters(c1, P, cand1, mixed_frac)
  if (verbose) message("  L1: ", paste(names(table(L1)), table(L1), collapse = " "))
  # ---- Levels 2-3 within each Level-1 group (non-epithelial) --------------------
  for (nd in setdiff(unique(L1), "Epithelial-lineage")) {
    g <- cells[L1 == nd]; Hn <- H[H$node1 == nd, ]
    c2 <- .sconco_cluster(seu, g, resolutions[2], seed = seed)
    L2[g] <- .assign_clusters(c2, P, split(Hn$Level4_Abb, Hn$node2), mixed_frac)
    for (n2 in unique(L2[g])) {
      g2 <- g[L2[g] == n2]; H2 <- Hn[Hn$node2 == n2, ]
      c3 <- .sconco_cluster(seu, g2, resolutions[3], seed = seed)
      L3[g2] <- .assign_clusters(c3, P, split(H2$Level4_Abb, H2$Level3), mixed_frac)
    }
  }
  # ---- copy-number node ---------------------------------------------------------
  # Candidates are not only the epithelial-lineage cells: a melanoma or a glioma
  # is not epithelial, its cells land in the Stromal / Other nodes (GSE72056 T81:
  # 70 of 133 melanoma cells clustered as stromal, 25 as "Other", and became the
  # diploid REFERENCE, so the tumour cancelled itself out - 0 % above the null).
  # Candidates = epithelial-lineage node + the "Other" node (neural / glial /
  # cycling: where a glioma goes) + any cell whose best type is a malignant
  # module. Reference = immune / stromal / endothelial nodes minus those cells.
  cnv_res <- NULL
  epi <- cells[L1 == "Epithelial-lineage"]
  best_all <- colnames(P)[max.col(P, "first")]; names(best_all) <- cells
  mal_types <- H$Level4_Abb[H$Level1 == "Malignant cell"]
  cand <- union(union(epi, cells[L1 == "Other"]), cells[best_all %in% mal_types])
  # Reference = the immune node: the one lineage a solid tumour does not
  # imitate. Stromal and endothelial cells are SECONDARY candidates: scored,
  # not in the null, called malignant only as a cluster coherent with the
  # primary clone template (sconco_cnv_call). With too few immune cells the
  # endothelial, then the stromal, node is added to the reference.
  ref <- setdiff(cells[L1 == "Immune cell"], cand)
  # revision (A1): a cross-sample immune pool, when given, fills a short immune
  # reference BEFORE the sample's own endothelial / stromal cells do - immune
  # cells of another sample are still the lineage no solid tumour imitates,
  # whereas stromal cells are what a mesenchymal-like tumour resembles.
  counts_cnv <- counts; pool_used <- 0L; pool_groups <- NULL
  if (cnv && !is.null(reference_pool) && length(ref) < pool_min_own) {
    pc <- setdiff(colnames(reference_pool), cells)
    if (length(pc)) {
      g <- intersect(rownames(counts), rownames(reference_pool))
      P_ <- reference_pool[g, pc, drop = FALSE]; if (!inherits(P_, "dgCMatrix")) P_ <- methods::as(methods::as(P_, "CsparseMatrix"), "dgCMatrix")
      counts_cnv <- cbind(counts[g, , drop = FALSE], P_); n_own <- length(ref); ref <- union(ref, pc); pool_used <- length(pc)
      org <- attr(reference_pool, "origin"); pool_groups <- setNames(rep("own", length(ref)), ref)
      if (!is.null(org)) pool_groups[pc] <- paste0("pool_", org[pc]) else pool_groups[pc] <- "pool"
      if (verbose) message("  reference pool: ", pool_used, " immune cells from ", length(unique(pool_groups[pc])), " other samples added to a reference of ", n_own) }
  }
  if (length(ref) < 50) ref <- union(ref, setdiff(cells[L1 == "Endothelial cell"], cand))
  if (length(ref) < 50) ref <- union(ref, setdiff(cells[L1 == "Stromal cell"], cand))
  secondary <- setdiff(cells[L1 %in% c("Stromal cell", "Endothelial cell")], union(ref, cand))
  if (cnv && length(cand) >= 30 && length(ref) >= 20 && !is.null(gene_pos)) {
    # anchors: candidate cells whose best type WITHIN THE EPITHELIAL LINEAGE
    # (normal types and malignant modules together) is an anchor type. "Best
    # among the normal epithelial types" made anchors out of melanoma cells (a
    # melanoma cell's best normal type is arbitrary); "best over all types"
    # competes a goblet cell against 35 immune types it never resembles.
    # Clone-like anchors are additionally dropped inside sconco_lineage_correct().
    ty_epi <- intersect(H$Level4_Abb[H$node1 == "Epithelial-lineage"], colnames(P))
    best_epi <- ty_epi[max.col(P[epi, ty_epi, drop = FALSE], "first")]
    anchors <- epi[best_epi %in% SCONCO_ANCHOR_TYPES]
    cnv_res <- do.call(sconco_cnv, c(list(counts = counts_cnv, ref_cells = ref, candidate_cells = cand, gene_pos = gene_pos,
                                          anchor_cells = anchors, ref_groups = if (group_reference) L3[ref] else pool_groups,
                                          secondary_cells = secondary, verbose = verbose), cnv_args))
    cand_all <- union(cand, secondary)
    cl <- cnv_res$call[cand_all, "cnv"]; names(cl) <- cand_all
    # epithelial-lineage cells: the copy-number call IS the Level-2 decision
    L2[epi] <- ifelse(cl[epi] == "aneuploid", "Malignant", ifelse(cl[epi] == "diploid", "Normal epithelium", "Epithelial (undetermined)"))
    # candidates from other nodes: aneuploid -> Malignant (Level 1 kept as the
    # marker lineage, so the user sees "Stromal cell / Malignant"); otherwise
    # they keep the marker annotation they already have
    oth <- setdiff(cand_all, epi); oth_mal <- oth[cl[oth] == "aneuploid"]
    if (length(oth_mal)) { L2[oth_mal] <- "Malignant"; L3[oth_mal] <- NA }
  } else if (length(epi)) {
    L2[epi] <- "Epithelial (undetermined)"          # no CNV run: markers are NOT asked to decide
  }
  # Level 3 within Malignant (cancer modules + Mal) / Normal epithelium (normal types), by cluster
  for (n2 in c("Malignant", "Normal epithelium")) {
    g2 <- cells[!is.na(L2) & L2 == n2]; if (!length(g2)) next
    H2 <- H[H$node2 == n2, ]
    c3 <- .sconco_cluster(seu, g2, resolutions[3], seed = seed)
    L3[g2] <- .assign_clusters(c3, P, split(H2$Level4_Abb, H2$Level3), mixed_frac)
  }
  # ---- Level 4 per cell among the children of the Level-3 node, with reject ----
  conf <- setNames(rep(NA_real_, length(cells)), cells); status <- setNames(rep("assigned", length(cells)), cells)
  for (n3 in unique(na.omit(L3))) {
    g3 <- cells[!is.na(L3) & L3 == n3]; ty <- intersect(H$Level4_Abb[H$Level3 == n3], colnames(P))
    if (!length(ty)) next
    Pp <- P[g3, ty, drop = FALSE]; Pp <- Pp / rowSums(Pp)
    # state modules join the choice only inside a Malignant-branch Level-3 node
    # that has them (v3.3: "Glioma / glioblastoma cell"), for cells whose Level 2
    # is Malignant; the renormalised softmax of the scores over children + states
    # is the same quantity Pp is for the children alone
    sty <- if (!is.null(S_state)) intersect(Hs$Level4_Abb[Hs$Level3 == n3 & Hs$node2 == "Malignant"], state_types) else character(0)
    sty <- intersect(sty, colnames(S_state))
    if (length(sty) && all(!is.na(L2[g3]) & L2[g3] == "Malignant")) {
      ty <- c(ty, sty); Z <- S_state[g3, ty, drop = FALSE]; Z <- Z - apply(Z, 1, max); Pp <- exp(Z); Pp <- Pp / rowSums(Pp) }
    top <- max.col(Pp, "first"); p1 <- Pp[cbind(seq_along(g3), top)]
    p2 <- if (ncol(Pp) > 1) apply(Pp, 1, function(r) sort(r, decreasing = TRUE)[2]) else 0
    L4[g3] <- ty[top]; conf[g3] <- p1
    # reject thresholds are relative to chance (1/k) when a node has many children:
    # a T-cell node with 15 subtypes cannot give p1 >= 0.3 to most cells, so the
    # absolute rule of sconco_annotate() rejected two thirds of a melanoma atlas
    k <- ncol(Pp); tp <- min(reject_p, 1.5 / k); tg <- min(reject_gap, 0.3 / k)
    status[g3][p1 < tp] <- "unassigned"; status[g3][p1 - p2 < tg & p1 >= tp] <- "ambiguous"
  }
  undet <- cells[!is.na(L2) & L2 == "Epithelial (undetermined)"]
  L3[undet] <- "Epithelial (undetermined)"; L4[undet] <- "Epi_undetermined"; status[undet] <- "deferred"
  res <- data.frame(cell = cells, L1 = L1, L2 = L2, L3 = L3, L4 = L4, conf = conf, status = status, row.names = cells)
  if (!is.null(cnv_res)) { res$cnv_burden <- NA_real_; res$p_malignant <- NA_real_
    ca <- union(cand, secondary); res[ca, "cnv_burden"] <- cnv_res$score[ca, "burden"]; res[ca, "p_malignant"] <- cnv_res$call[ca, "p_aneuploid"] }
  tt <- if (!is.null(cnv_res)) cnv_res$test else NULL
  if (!is.null(tt)) tt$n_pool <- pool_used
  out <- list(annotation = res, tumour_test = tt, probs = P, clusters_L1 = c1)
  if (keep_cnv && !is.null(cnv_res)) out$cnv <- c(cnv_res, list(ref_cells = ref, candidate_cells = cand))   # profile, windows, score, call (development / diagnostics)
  out
}
