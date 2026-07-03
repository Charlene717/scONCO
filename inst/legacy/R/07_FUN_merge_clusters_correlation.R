###############################################################################
# Script: MergeSimilar_NumberedCellTypes_ByAvgExprCorrelation_v4.R
# Fix v4:
#  - drop_suffix_if_singleton_prefix(): avoid dplyr::count(name=) incompatibility
#    -> use base R table() to count prefixes
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(tibble)
  library(Matrix)
})

merge_similar_numbered_celltypes_by_correlation <- function(
    obj,
    meta_col = "Cell_Type",
    assay = "RNA",
    layer_prefer = c("data", "counts"),
    
    # correlation settings
    cor_method = c("pearson", "spearman"),
    cor_cutoff = 0.9,
    use_features = c("variable", "all", "custom"),
    features_custom = NULL,
    n_features = 2000,
    
    # compare mode
    compare_mode = c("adjacent", "all_pairs"),
    
    # guards
    min_cells_per_cluster = 20,
    min_features_for_cor  = 200,
    
    # iteration / safety
    max_iter = 20,
    verbose = TRUE
) {
  stopifnot(inherits(obj, "Seurat"))
  stopifnot(meta_col %in% colnames(obj@meta.data))
  
  cor_method   <- match.arg(cor_method)
  use_features <- match.arg(use_features)
  compare_mode <- match.arg(compare_mode)
  layer_prefer <- match.arg(layer_prefer)
  
  # ---- helper: parse labels ending with _number
  parse_df <- function(labels) {
    labels <- trimws(as.character(labels))
    tibble::tibble(label = labels) %>%
      dplyr::mutate(
        has_num = stringr::str_detect(label, "_\\d+$"),
        prefix  = dplyr::if_else(has_num, stringr::str_replace(label, "_\\d+$", ""), NA_character_),
        idx     = dplyr::if_else(has_num, as.integer(stringr::str_extract(label, "\\d+$")), NA_integer_)
      )
  }
  
  # ---- helper: apply label mapping
  apply_label_map <- function(obj, map_vec_named) {
    v <- as.character(obj[[meta_col]][, 1])
    hit <- v %in% names(map_vec_named)
    v[hit] <- unname(map_vec_named[v[hit]])
    obj[[meta_col]] <- v
    obj
  }
  
  # ---- helper: renumber within prefix sequentially
  renumber_prefix <- function(obj, prefix) {
    labs <- unique(as.character(obj[[meta_col]][, 1]))
    df <- parse_df(labs) %>%
      dplyr::filter(has_num, prefix == !!prefix) %>%
      dplyr::arrange(idx)
    
    if (nrow(df) == 0) return(obj)
    
    new_labels <- paste0(prefix, "_", seq_len(nrow(df)))
    rename_map <- stats::setNames(new_labels, df$label)
    apply_label_map(obj, rename_map)
  }
  
  # ---- helper: if a prefix ends up with only "prefix_1", drop suffix -> "prefix"
  # (base R version to avoid dplyr::count(name=...) compatibility issues)
  drop_suffix_if_singleton_prefix <- function(obj) {
    labs <- unique(as.character(obj[[meta_col]][, 1]))
    df <- parse_df(labs)
    df <- df[df$has_num %in% TRUE, , drop = FALSE]
    if (nrow(df) == 0) return(obj)
    
    # count numbered labels per prefix (base R)
    k_tab <- table(df$prefix)
    one_prefix <- names(k_tab)[k_tab == 1]
    if (length(one_prefix) == 0) return(obj)
    
    map_vec <- c()
    for (px in one_prefix) {
      only_label <- df$label[df$prefix == px]
      if (length(only_label) == 1 && identical(only_label, paste0(px, "_1"))) {
        map_vec[only_label] <- px
      }
    }
    if (length(map_vec) == 0) return(obj)
    
    apply_label_map(obj, map_vec)
  }
  
  # ---- helper: robust get assay matrix (Seurat v4/v5)
  .get_assay_mat <- function(obj, assay, layer_try = c("data", "counts")) {
    Seurat::DefaultAssay(obj) <- assay
    for (ly in layer_try) {
      mat <- NULL
      mat <- tryCatch(Seurat::GetAssayData(obj, assay = assay, slot = ly), error = function(e) NULL)
      if (is.null(mat)) {
        mat <- tryCatch(Seurat::GetAssayData(obj, assay = assay, layer = ly), error = function(e) NULL)
      }
      if (!is.null(mat) && nrow(mat) > 0 && ncol(mat) > 0) {
        return(list(mat = mat, used = ly))
      }
    }
    stop("Cannot retrieve assay matrix for assay=", assay,
         " with layers/slots: ", paste(layer_try, collapse = ", "))
  }
  
  # ---- helper: pick features
  pick_features <- function(obj) {
    Seurat::DefaultAssay(obj) <- assay
    
    if (use_features == "custom") {
      feats <- unique(features_custom %||% character())
      feats <- feats[feats %in% rownames(obj)]
      if (length(feats) < 50) stop("features_custom too small after intersect with genes.")
      return(feats)
    }
    
    if (use_features == "all") {
      return(rownames(obj))
    }
    
    feats <- tryCatch(Seurat::VariableFeatures(obj), error = function(e) character())
    if (length(feats) == 0) {
      obj <- tryCatch(Seurat::FindVariableFeatures(obj, nfeatures = max(2000, n_features)),
                      error = function(e) obj)
      feats <- tryCatch(Seurat::VariableFeatures(obj), error = function(e) character())
    }
    if (length(feats) == 0) {
      g <- .get_assay_mat(obj, assay, layer_try = c("data", "counts"))$mat
      av <- Matrix::rowMeans(g)
      feats <- names(sort(av, decreasing = TRUE))[seq_len(min(length(av), n_features))]
    }
    feats <- feats[seq_len(min(length(feats), n_features))]
    feats
  }
  
  # ---- helper: pseudo-bulk mean expression per cluster (manual)
  compute_pseudobulk_means <- function(mat, clusters, cells_by_cluster) {
    out <- matrix(NA_real_, nrow = nrow(mat), ncol = length(clusters))
    rownames(out) <- rownames(mat)
    colnames(out) <- clusters
    for (k in seq_along(clusters)) {
      cl <- clusters[k]
      cells <- cells_by_cluster[[cl]]
      if (length(cells) == 0) next
      sub <- mat[, cells, drop = FALSE]
      out[, k] <- Matrix::rowMeans(sub)
    }
    out
  }
  
  # ---- helper: correlation for two vectors with guards
  cor_vec <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]; y <- y[ok]
    if (length(x) < min_features_for_cor) return(NA_real_)
    if (stats::sd(x) == 0 || stats::sd(y) == 0) return(NA_real_)
    suppressWarnings(stats::cor(x, y, method = cor_method, use = "pairwise.complete.obs"))
  }
  
  # ============================================================
  # Iterative merging loop
  # ============================================================
  history <- list()
  iter <- 0
  
  repeat {
    iter <- iter + 1
    if (iter > max_iter) {
      warning("Reached max_iter = ", max_iter, ". Stopping.")
      break
    }
    
    labs <- unique(as.character(obj[[meta_col]][, 1]))
    df_labs <- parse_df(labs)
    prefixes <- df_labs %>% dplyr::filter(has_num) %>% dplyr::distinct(prefix) %>% dplyr::pull(prefix)
    
    if (length(prefixes) == 0) {
      if (verbose) message("[Iter ", iter, "] No suffix _number labels detected. Done.")
      break
    }
    
    feats <- pick_features(obj)
    
    mat_info <- .get_assay_mat(obj, assay,
                               layer_try = c(layer_prefer, ifelse(layer_prefer == "data", "counts", "data")))
    mat_all <- mat_info$mat
    used_layer <- mat_info$used
    
    if (verbose) {
      message("\n==============================")
      message("[Iter ", iter, "] prefixes: ", paste(prefixes, collapse = ", "))
      message("[Iter ", iter, "] assay=", assay, " | used_layer/slot=", used_layer,
              " | cor_method=", cor_method, " | cor_cutoff=", cor_cutoff,
              " | features=", use_features, " (n=", length(feats), ")",
              " | compare_mode=", compare_mode)
    }
    
    any_merge <- FALSE
    Seurat::Idents(obj) <- obj[[meta_col]][, 1]
    
    for (px in prefixes) {
      cur_labels <- unique(as.character(obj[[meta_col]][, 1]))
      cur_df <- parse_df(cur_labels) %>% dplyr::filter(has_num, prefix == px) %>% dplyr::arrange(idx)
      if (nrow(cur_df) <= 1) next
      
      tab <- table(Seurat::Idents(obj))
      keep <- intersect(cur_df$label, names(tab)[tab >= min_cells_per_cluster])
      if (length(keep) <= 1) next
      
      merged_in_prefix <- TRUE
      while (merged_in_prefix) {
        merged_in_prefix <- FALSE
        
        Seurat::Idents(obj) <- obj[[meta_col]][, 1]
        cur_labels <- unique(as.character(obj[[meta_col]][, 1]))
        cur_df <- parse_df(cur_labels) %>% dplyr::filter(has_num, prefix == px) %>% dplyr::arrange(idx)
        tab <- table(Seurat::Idents(obj))
        keep <- intersect(cur_df$label, names(tab)[tab >= min_cells_per_cluster])
        if (length(keep) <= 1) break
        
        cells_by_cluster <- lapply(keep, function(cl) Seurat::WhichCells(obj, idents = cl))
        names(cells_by_cluster) <- keep
        
        feats_use <- feats[feats %in% rownames(mat_all)]
        if (length(feats_use) < min_features_for_cor) break
        
        cells_use <- unique(unlist(cells_by_cluster, use.names = FALSE))
        mat_sub <- mat_all[feats_use, cells_use, drop = FALSE]
        
        pb <- compute_pseudobulk_means(mat_sub, keep, cells_by_cluster)
        
        if (compare_mode == "adjacent") {
          df_keep <- parse_df(keep) %>% dplyr::arrange(idx)
          if (nrow(df_keep) <= 1) break
          pairs <- tibble::tibble(
            a = df_keep$label[-nrow(df_keep)],
            b = df_keep$label[-1]
          )
        } else {
          if (length(keep) <= 1) break
          comb <- t(combn(keep, 2))
          pairs <- tibble::tibble(a = comb[, 1], b = comb[, 2])
        }
        
        pairs$cor <- mapply(function(a, b) cor_vec(pb[, a], pb[, b]), pairs$a, pairs$b)
        
        if (all(is.na(pairs$cor))) {
          if (verbose) {
            message("  [", px, "] all correlations are NA. ",
                    "Try layer_prefer='counts' or use_features='all'.")
          }
          break
        }
        
        pairs <- pairs[order(-pairs$cor), , drop = FALSE]
        best <- pairs[1, , drop = FALSE]
        best_cor <- best$cor[1]
        a_best <- best$a[1]
        b_best <- best$b[1]
        
        if (verbose) {
          message(sprintf("  [%s] best pair %s vs %s : cor=%.3f (%s)",
                          px, a_best, b_best, best_cor, compare_mode))
        }
        
        if (is.finite(best_cor) && best_cor >= cor_cutoff) {
          idx_a <- parse_df(a_best)$idx[1]
          idx_b <- parse_df(b_best)$idx[1]
          into <- if (idx_a <= idx_b) a_best else b_best
          from <- if (idx_a <= idx_b) b_best else a_best
          
          obj <- apply_label_map(obj, stats::setNames(into, from))
          obj <- renumber_prefix(obj, px)
          
          any_merge <- TRUE
          merged_in_prefix <- TRUE
          
          history[[length(history) + 1]] <- tibble::tibble(
            iter = iter,
            prefix = px,
            merged_from = from,
            merged_into = into,
            similarity = best_cor,
            sim_method = paste0("cor_", cor_method),
            compare_mode = compare_mode,
            assay = assay,
            used_layer = used_layer,
            n_features = length(feats_use),
            min_cells_per_cluster = min_cells_per_cluster
          )
          
          if (verbose) message("    ✅ MERGE: ", from, " -> ", into,
                               " | cor=", sprintf("%.3f", best_cor),
                               " | renumber within ", px)
        } else {
          break
        }
      }
    }
    
    if (!any_merge) {
      if (verbose) message("[Iter ", iter, "] No merges. Converged ✅")
      break
    } else {
      if (verbose) message("[Iter ", iter, "] Merges happened; continue...")
    }
  }
  
  # ✅ post-process: drop "_1" when a prefix ends up singleton "prefix_1"
  obj <- drop_suffix_if_singleton_prefix(obj)
  
  history_df <- if (length(history) > 0) dplyr::bind_rows(history) else tibble::tibble()
  
  final_counts <- as.data.frame(table(as.character(obj[[meta_col]][, 1])), stringsAsFactors = FALSE)
  colnames(final_counts) <- c("Cell_Type", "n")
  final_counts <- final_counts[order(-final_counts$n, final_counts$Cell_Type), , drop = FALSE]
  rownames(final_counts) <- NULL
  
  list(
    obj = obj,
    merge_history = history_df,
    final_counts = final_counts
  )
}

###############################################################################
# ✅ Example usage
###############################################################################
res_merge <- merge_similar_numbered_celltypes_by_correlation(
  obj = seuratObject_Sample,
  meta_col = "Cell_Type",
  assay = "RNA",
  layer_prefer = "data",        # 若仍 NA，改成 "counts"
  cor_method = "pearson",
  cor_cutoff = 0.9,
  use_features = "variable",
  n_features = 2000,
  compare_mode = "all_pairs",
  min_cells_per_cluster = 20,
  min_features_for_cor = 200,
  max_iter = 20,
  verbose = TRUE
)
seuratObject_Sample <- res_merge$obj
res_merge$merge_history
res_merge$final_counts
sort(unique(seuratObject_Sample$Cell_Type))

