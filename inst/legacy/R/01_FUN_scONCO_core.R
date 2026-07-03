###############################################
# scONCO core algorithm (R implementation)
# - Gene x cell expression matrix
# - Marker table with pan-cancer-specific ontology
# - Specificity-weighted scoring
# - Hierarchical annotation with confidence
###############################################

#' Softmax helper (for numeric vector)
softmax_vec <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(x)
  z <- x - max(x)
  exp_z <- exp(z)
  exp_z / sum(exp_z)
}

#' Ensure a column exists, otherwise create default
ensure_col <- function(df, col, default = 0) {
  if (!col %in% colnames(df)) {
    df[[col]] <- default
  }
  df
}

#' Compute marker weights for scONCO
#'
#' @param marker_df data.frame with at least:
#'        gene, cell_type, lineage
#'        optional: negative (0/1), ai_support (0/1), lit_support (0/1)
#' @param alpha weight for specificity
#' @param beta  weight for AI support
#' @param gamma weight for literature support
#' @return marker_df with added columns: freq, specificity, weight
sconco_compute_marker_weights <- function(marker_df,
                                           alpha = 0.6,
                                           beta  = 0.2,
                                           gamma = 0.2) {
  stopifnot("gene"      %in% colnames(marker_df),
            "cell_type" %in% colnames(marker_df),
            "lineage"   %in% colnames(marker_df))
  
  # Ensure optional columns exist
  marker_df <- ensure_col(marker_df, "negative",   0)
  marker_df <- ensure_col(marker_df, "ai_support", 0)
  marker_df <- ensure_col(marker_df, "lit_support",0)
  
  # Frequency of each gene across cell types (positive markers only)
  pos_idx <- marker_df$negative == 0
  gene_cell_counts <- table(marker_df$gene[pos_idx])
  total_ctypes      <- length(unique(marker_df$cell_type))
  
  freq_vec <- as.numeric(gene_cell_counts) / total_ctypes
  names(freq_vec) <- names(gene_cell_counts)
  
  # Map frequency to each marker row
  marker_df$freq <- freq_vec[marker_df$gene]
  marker_df$freq[is.na(marker_df$freq)] <- 0
  
  # Specificity: higher for genes used in fewer cell types
  marker_df$specificity <- 1 - marker_df$freq
  
  # Normalize AI/Lit support to [0,1]
  marker_df$ai_support  <- pmin(pmax(marker_df$ai_support, 0), 1)
  marker_df$lit_support <- pmin(pmax(marker_df$lit_support, 0), 1)
  
  # Composite weight for positive markers
  marker_df$weight <- alpha * marker_df$specificity +
    beta  * marker_df$ai_support   +
    gamma * marker_df$lit_support
  
  # For negative markers, we still keep weight but will use as penalty
  marker_df$weight[marker_df$negative == 1] <- 
    marker_df$weight[marker_df$negative == 1] # could re-scale if needed
  
  return(marker_df)
}

#' Compute activity scores for each cell and each cell type
#'
#' @param expr_mat matrix genes x cells (rownames = genes, colnames = cells)
#' @param marker_df marker table with columns:
#'        gene, cell_type, weight, negative (0/1)
#' @param delta penalty factor for negative markers
#' @return list with:
#'         scores: matrix cells x cell_types (raw scores)
#'         probs : matrix cells x cell_types (softmax-normalized)
sconco_activity_scores <- function(expr_mat,
                                    marker_df,
                                    delta = 0.2) {
  # Check expression orientation
  if (is.null(rownames(expr_mat))) {
    stop("expr_mat must have rownames = gene symbols.")
  }
  
  # Ensure necessary columns
  stopifnot(all(c("gene", "cell_type", "negative", "weight") %in% colnames(marker_df)))
  
  # Filter markers that exist in expression
  valid_genes <- intersect(marker_df$gene, rownames(expr_mat))
  marker_df   <- marker_df[marker_df$gene %in% valid_genes, , drop = FALSE]
  if (nrow(marker_df) == 0L) {
    stop("No markers overlap with expression matrix rownames.")
  }
  
  cell_types <- sort(unique(marker_df$cell_type))
  n_cells    <- ncol(expr_mat)
  n_ctypes   <- length(cell_types)
  
  score_mat  <- matrix(NA_real_, nrow = n_cells, ncol = n_ctypes)
  colnames(score_mat) <- cell_types
  rownames(score_mat) <- colnames(expr_mat)
  
  # Pre-split marker_df by cell_type for efficiency
  split_by_ct <- split(marker_df, marker_df$cell_type)
  
  for (ct in cell_types) {
    msub <- split_by_ct[[ct]]
    if (is.null(msub)) next
    
    # Positive markers
    pos_genes <- msub$gene[msub$negative == 0]
    pos_w     <- msub$weight[msub$negative == 0]
    
    # Negative markers
    neg_genes <- msub$gene[msub$negative == 1]
    neg_w     <- msub$weight[msub$negative == 1]
    
    # Positive contribution
    pos_score <- rep(0, n_cells)
    if (length(pos_genes) > 0) {
      common_pos <- intersect(pos_genes, rownames(expr_mat))
      if (length(common_pos) > 0) {
        # Align weights
        w_vec <- pos_w[match(common_pos, pos_genes)]
        # expr_mat[genes x cells], multiply each gene row by its weight
        expr_sub <- expr_mat[common_pos, , drop = FALSE]
        # scale rows by weights
        expr_weighted <- expr_sub * matrix(w_vec, nrow = length(common_pos), ncol = n_cells, byrow = FALSE)
        pos_score <- colSums(expr_weighted)
      }
    }
    
    # Negative penalty
    neg_penalty <- rep(0, n_cells)
    if (length(neg_genes) > 0) {
      common_neg <- intersect(neg_genes, rownames(expr_mat))
      if (length(common_neg) > 0) {
        w_neg <- neg_w[match(common_neg, neg_genes)]
        expr_neg <- expr_mat[common_neg, , drop = FALSE]
        expr_weighted_neg <- expr_neg * matrix(w_neg, nrow = length(common_neg), ncol = n_cells, byrow = FALSE)
        neg_penalty <- delta * colSums(expr_weighted_neg)
      }
    }
    
    score_mat[, ct] <- pos_score - neg_penalty
  }
  
  # Convert raw scores to probabilities via softmax per cell
  prob_mat <- t(apply(score_mat, 1, softmax_vec))
  colnames(prob_mat) <- colnames(score_mat)
  rownames(prob_mat) <- rownames(score_mat)
  
  return(list(scores = score_mat, probs = prob_mat))
}

#' Hierarchical annotation with confidence-based decisions
#'
#' @param expr_mat  gene x cell expression matrix
#' @param marker_df marker table (with lineage, etc.)
#' @param alpha, beta, gamma specificity/AI/lit weights (passed to sconco_compute_marker_weights)
#' @param delta penalty factor for negative markers
#' @param conf_low minimum score gap to avoid ambiguous label
#' @param min_prob minimum probability to accept a label
#'
#' @return list with:
#'         cell_type   : final label per cell
#'         lineage     : major lineage per cell
#'         conf        : confidence (score gap) per cell
#'         prob_matrix : cell x cell_type probabilities
#'         score_matrix: cell x cell_type raw scores
sconco_annotate <- function(expr_mat,
                             marker_df,
                             alpha    = 0.6,
                             beta     = 0.2,
                             gamma    = 0.2,
                             delta    = 0.2,
                             conf_low = 0.1,
                             min_prob = 0.3) {
  stopifnot("lineage" %in% colnames(marker_df))
  
  # 1. Compute marker weights
  marker_df <- sconco_compute_marker_weights(marker_df,
                                              alpha = alpha,
                                              beta  = beta,
                                              gamma = gamma)
  # 2. Compute scores & probabilities
  as_res <- sconco_activity_scores(expr_mat,
                                    marker_df,
                                    delta = delta)
  score_mat <- as_res$scores  # cells x cell_types
  prob_mat  <- as_res$probs
  
  cell_ids  <- rownames(score_mat)
  cell_types <- colnames(score_mat)
  
  # 3. Derive lineage mapping for each cell type
  ct2lin <- marker_df[, c("cell_type", "lineage")]
  ct2lin <- unique(ct2lin)
  lineage_vec <- setNames(ct2lin$lineage, ct2lin$cell_type)
  
  # Some cell_types may share same name but different lineage (should be rare), handle by first occurrence
  lineage_vec <- lineage_vec[!duplicated(names(lineage_vec))]
  
  # 4. Stage 1: major lineage assignment (summed probability per lineage)
  lineages <- sort(unique(marker_df$lineage))
  lin_score <- matrix(0, nrow = length(cell_ids), ncol = length(lineages))
  rownames(lin_score) <- cell_ids
  colnames(lin_score) <- lineages
  
  for (lin in lineages) {
    # all cell types in this lineage
    ct_in_lin <- names(lineage_vec)[lineage_vec == lin]
    common_ct <- intersect(ct_in_lin, cell_types)
    if (length(common_ct) == 0) next
    lin_score[, lin] <- rowSums(prob_mat[, common_ct, drop = FALSE])
  }
  
  major_lineage <- apply(lin_score, 1, function(x) {
    if (all(x == 0 | is.na(x))) return(NA_character_)
    names(which.max(x))
  })
  
  # 5. Stage 2: subtype assignment within each lineage
  final_label <- rep("unassigned", length(cell_ids))
  names(final_label) <- cell_ids
  conf_vec <- rep(NA_real_, length(cell_ids))
  names(conf_vec) <- cell_ids
  
  for (cid in cell_ids) {
    lin <- major_lineage[[cid]]
    if (is.na(lin)) {
      final_label[[cid]] <- "unassigned"
      conf_vec[[cid]]    <- NA_real_
      next
    }
    # cell types in this lineage
    ct_in_lin <- names(lineage_vec)[lineage_vec == lin]
    common_ct <- intersect(ct_in_lin, cell_types)
    if (length(common_ct) == 0) {
      final_label[[cid]] <- paste0("lineage_", lin, "_noSubtype")
      conf_vec[[cid]]    <- NA_real_
      next
    }
    
    # Probabilities & scores for this cell within lineage
    probs_ct  <- prob_mat[cid, common_ct, drop = TRUE]
    scores_ct <- score_mat[cid, common_ct, drop = TRUE]
    
    # Order by score (raw scores, not probs)
    ord <- order(scores_ct, decreasing = TRUE, na.last = TRUE)
    scores_ord <- scores_ct[ord]
    ct_ord     <- names(scores_ord)
    
    if (length(ct_ord) == 0 || all(!is.finite(scores_ord))) {
      final_label[[cid]] <- "unassigned"
      conf_vec[[cid]]    <- NA_real_
      next
    }
    
    top1 <- ct_ord[1]
    top1_score <- scores_ord[1]
    top1_prob  <- probs_ct[top1]
    
    if (length(ct_ord) > 1) {
      top2_score <- scores_ord[2]
    } else {
      top2_score <- -Inf
    }
    
    conf <- top1_score - top2_score
    conf_vec[[cid]] <- conf
    
    if (!is.finite(top1_prob) || is.na(top1_prob) || top1_prob < min_prob) {
      final_label[[cid]] <- "unassigned"
    } else if (!is.finite(conf) || conf < conf_low) {
      # ambiguous between top2
      if (length(ct_ord) > 1) {
        final_label[[cid]] <- paste0("ambiguous_", top1, "_", ct_ord[2])
      } else {
        final_label[[cid]] <- top1
      }
    } else {
      final_label[[cid]] <- top1
    }
  }
  
  res <- list(
    cell_type    = final_label,
    lineage      = major_lineage,
    conf         = conf_vec,
    prob_matrix  = prob_mat,
    score_matrix = score_mat,
    lineage_probs= lin_score
  )
  class(res) <- c("scONCO_annotation", class(res))
  return(res)
}

#' Simple print method for scONCO_annotation
print.scONCO_annotation <- function(x, ...) {
  cat("scONCO annotation result\n")
  cat("Number of cells:", length(x$cell_type), "\n")
  cat("Number of unique labels:", length(unique(x$cell_type)), "\n")
  cat("Head of labels:\n")
  print(head(x$cell_type))
  invisible(x)
}

###############################################
# Seurat wrapper
###############################################

#' Run scONCO annotation on a Seurat object
#'
#' @param seurat_obj A Seurat object
#' @param marker_df  Marker table used by scONCO
#' @param assay      Assay name to use (default: DefaultAssay(seurat_obj))
#' @param slot       Which slot to use from the assay ("data" or "counts")
#' @param alpha,beta,gamma,delta,conf_low,min_prob Parameters passed to sconco_annotate()
#'
#' @return Seurat object with new metadata columns:
#'         scONCO_celltype, scONCO_lineage, scONCO_conf
#'
sconco_annotate_seurat <- function(seurat_obj,
                                    marker_df,
                                    assay    = NULL,
                                    slot     = "data",
                                    alpha    = 0.6,
                                    beta     = 0.2,
                                    gamma    = 0.2,
                                    delta    = 0.2,
                                    conf_low = 0.1,
                                    min_prob = 0.3) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required for sconco_annotate_seurat(). Please install it.")
  }
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seurat_obj)
  }
  
  # Get expression matrix: genes x cells
  expr_mat <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = slot)
  expr_mat <- as.matrix(expr_mat)  # ensure it's a matrix
  
  # Run core scONCO annotation
  annot <- sconco_annotate(expr_mat,
                            marker_df  = marker_df,
                            alpha      = alpha,
                            beta       = beta,
                            gamma      = gamma,
                            delta      = delta,
                            conf_low   = conf_low,
                            min_prob   = min_prob)
  
  # Align to cell order of Seurat object
  cells <- colnames(seurat_obj)
  
  seurat_obj$scONCO_celltype <- annot$cell_type[cells]
  seurat_obj$scONCO_lineage  <- annot$lineage[cells]
  seurat_obj$scONCO_conf     <- annot$conf[cells]
  
  return(seurat_obj)
}


#' Run scONCO annotation on Seurat clusters (cluster-based annotation)
#'
#' @param seurat_obj A Seurat object
#' @param marker_df  Marker table used by scONCO
#' @param assay      Assay name to use (default: DefaultAssay(seurat_obj))
#' @param slot       Which slot to use from the assay ("data" or "counts")
#' @param cluster_col Name of metadata column indicating clusters.
#'        If NULL (default), use Idents(seurat_obj).
#' @param alpha,beta,gamma,delta,conf_low,min_prob Parameters passed to sconco_annotate()
#'
#' @return Seurat object with new metadata columns:
#'         scONCO_cluster                  = cluster ID per cell
#'         scONCO_cluster_celltype         = scONCO label per cluster (mapped to cell)
#'         scONCO_cluster_lineage          = major lineage per cluster (mapped to cell)
#'         scONCO_cluster_conf             = confidence per cluster (mapped to cell)
#'
sconco_annotate_seurat_clusters <- function(seurat_obj,
                                             marker_df,
                                             assay      = NULL,
                                             slot       = "data",
                                             cluster_col= NULL,
                                             alpha      = 0.6,
                                             beta       = 0.2,
                                             gamma      = 0.2,
                                             delta      = 0.2,
                                             conf_low   = 0.1,
                                             min_prob   = 0.3) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required for sconco_annotate_seurat_clusters(). Please install it.")
  }
  
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(seurat_obj)
  }
  
  # Determine cluster labels for each cell
  if (is.null(cluster_col)) {
    clusters <- Seurat::Idents(seurat_obj)
  } else {
    if (!cluster_col %in% colnames(seurat_obj@meta.data)) {
      stop("cluster_col not found in seurat_obj@meta.data.")
    }
    clusters <- seurat_obj@meta.data[[cluster_col]]
  }
  clusters <- as.factor(clusters)
  
  # Save cluster ID into metadata (for convenience)
  seurat_obj$scONCO_cluster <- clusters
  
  # Compute average expression per cluster
  avg_list <- Seurat::AverageExpression(
    seurat_obj,
    assays  = assay,
    slot    = slot,
    group.by = "scONCO_cluster",
    verbose = FALSE
  )
  
  # AverageExpression returns a list per assay
  if (!assay %in% names(avg_list)) {
    stop("Assay not found in AverageExpression output.")
  }
  expr_mat_cluster <- avg_list[[assay]]  # genes x clusters
  
  expr_mat_cluster <- as.matrix(expr_mat_cluster)
  
  # Run scONCO on cluster-level expression
  annot_cluster <- sconco_annotate(
    expr_mat   = expr_mat_cluster,
    marker_df  = marker_df,
    alpha      = alpha,
    beta       = beta,
    gamma      = gamma,
    delta      = delta,
    conf_low   = conf_low,
    min_prob   = min_prob
  )
  
  # Map cluster labels back to cells
  cluster_ids           <- colnames(expr_mat_cluster)  # cluster levels
  cluster_celltype_map  <- annot_cluster$cell_type
  cluster_lineage_map   <- annot_cluster$lineage
  cluster_conf_map      <- annot_cluster$conf
  
  # Ensure names align
  # annot_cluster$cell_type 是以 cluster 名為 rowname
  # 我們依照 factor level mapping 回去
  ct_per_cell   <- cluster_celltype_map[as.character(clusters)]
  lin_per_cell  <- cluster_lineage_map[as.character(clusters)]
  conf_per_cell <- cluster_conf_map[as.character(clusters)]
  
  seurat_obj$scONCO_cluster_celltype <- ct_per_cell
  seurat_obj$scONCO_cluster_lineage  <- lin_per_cell
  seurat_obj$scONCO_cluster_conf     <- conf_per_cell
  
  return(seurat_obj)
}


###############################################
# Example usage with a Seurat object
###############################################

# library(Seurat)
#
# # 假設你已經有一個 Seurat 物件，名為 seu
#
# marker_df <- data.frame(
#   gene        = c("KRT14","KRT5","KRT10","KRT1","COL1A1","COL3A1"),
#   cell_type   = c("Basal_KC","Basal_KC","Spinous_KC","Spinous_KC",
#                   "Papillary_FB","Papillary_FB"),
#   lineage     = c("Keratinocyte","Keratinocyte","Keratinocyte","Keratinocyte",
#                   "Fibroblast","Fibroblast"),
#   negative    = c(0,0,0,0,0,0),
#   ai_support  = c(1,1,1,1,0,0),
#   lit_support = c(1,1,1,1,1,1),
#   stringsAsFactors = FALSE
# )
#
# seu <- sconco_annotate_seurat(
#   seurat_obj = seu,
#   marker_df  = marker_df,
#   assay      = NULL,   # use DefaultAssay(seu)
#   slot       = "data"  # use normalized data
# )
#
# head(seu$scONCO_celltype)
# head(seu$scONCO_lineage)
# head(seu$scONCO_conf)
###############################################

###############################################
# Example usage (Seurat)
###############################################
# library(Seurat)
#
# marker_df <- data.frame(
#   gene        = c("KRT14","KRT5","KRT10","KRT1","COL1A1","COL3A1"),
#   cell_type   = c("Basal_KC","Basal_KC","Spinous_KC","Spinous_KC",
#                   "Papillary_FB","Papillary_FB"),
#   lineage     = c("Keratinocyte","Keratinocyte","Keratinocyte","Keratinocyte",
#                   "Fibroblast","Fibroblast"),
#   negative    = c(0,0,0,0,0,0),
#   ai_support  = c(1,1,1,1,0,0),
#   lit_support = c(1,1,1,1,1,1),
#   stringsAsFactors = FALSE
# )
#
# # cell-based annotation
# seu <- sconco_annotate_seurat(seu, marker_df)
#
# # cluster-based annotation (using Idents)
# seu <- sconco_annotate_seurat_clusters(seu, marker_df)

# head(seu$scONCO_celltype)
# head(seu$scONCO_lineage)
# head(seu$scONCO_conf)
###############################################




###############################################
# Example usage (pseudo-code, not executed):
#
# expr_mat  <- log1p(counts_mat)  # genes x cells
#
# marker_df <- data.frame(
#   gene        = c("KRT14","KRT5","KRT10","KRT1","COL1A1","COL3A1"),
#   cell_type   = c("Basal_KC","Basal_KC","Spinous_KC","Spinous_KC",
#                   "Papillary_FB","Papillary_FB"),
#   lineage     = c("Keratinocyte","Keratinocyte","Keratinocyte","Keratinocyte",
#                   "Fibroblast","Fibroblast"),
#   negative    = c(0,0,0,0,0,0),
#   ai_support  = c(1,1,1,1,0,0),
#   lit_support = c(1,1,1,1,1,1),
#   stringsAsFactors = FALSE
# )
#
# annot <- sconco_annotate(expr_mat, marker_df)
# head(annot$cell_type)
###############################################



###############################################
# -[] Debug
# -[] 單一細胞標註和群集細胞標註函數合併改用參數設定
# -[] marker_df 的設計是否要改回多層
