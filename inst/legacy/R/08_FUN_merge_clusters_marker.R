###############################################################################
# Script: MergeSimilar_NumberedCellTypes_RecursiveRenumber.R
# Purpose:
#   - Detect labels ending with _number (e.g., FB_1, vEC_2, KC_Spinous_6)
#   - Within each prefix, compute marker sets for each numbered cluster
#   - Compare adjacent numbered clusters by Jaccard(marker sets)
#   - If similar -> merge, renumber, and iterate recursively until stable
#   - Marker mode selectable: "global" vs "within_prefix"
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(stringr)
  library(tibble)
})

# ============================================================
# Helper: robust FindAllMarkers wrapper (Seurat v4/v5 friendly)
# ============================================================
.safe_FindAllMarkers <- function(
    obj,
    assay = NULL,
    slot = "data",
    only.pos = TRUE,
    logfc.threshold = 0.25,
    min.pct = 0.1,
    test.use = "wilcox",
    verbose = FALSE
) {
  if (!is.null(assay)) Seurat::DefaultAssay(obj) <- assay
  
  # SCT often needs PrepSCTFindMarkers
  if (Seurat::DefaultAssay(obj) == "SCT") {
    obj <- tryCatch(Seurat::PrepSCTFindMarkers(obj, verbose = verbose),
                    error = function(e) obj)
  }
  
  out <- tryCatch({
    Seurat::FindAllMarkers(
      object = obj,
      only.pos = only.pos,
      logfc.threshold = logfc.threshold,
      min.pct = min.pct,
      test.use = test.use,
      slot = slot,
      verbose = verbose
    )
  }, error = function(e1) {
    # fallback without slot (for some edge cases)
    tryCatch({
      Seurat::FindAllMarkers(
        object = obj,
        only.pos = only.pos,
        logfc.threshold = logfc.threshold,
        min.pct = min.pct,
        test.use = test.use,
        verbose = verbose
      )
    }, error = function(e2) {
      stop("FindAllMarkers failed. First error: ", conditionMessage(e1),
           " | Second error: ", conditionMessage(e2))
    })
  })
  
  out
}

# ============================================================
# Main: merge similar numbered labels within each prefix
# ============================================================
merge_similar_numbered_celltypes <- function(
    obj,
    meta_col = "Cell_Type",
    assay = NULL,
    slot = "data",
    
    # marker mode
    marker_mode = c("global", "within_prefix"),  # ✅ global = vs all other cells; within_prefix = vs same prefix only
    
    # marker settings
    only.pos = TRUE,
    logfc.threshold = 0.25,
    min.pct = 0.1,
    test.use = "wilcox",
    top_n_markers = 50,
    
    # similarity settings
    jaccard_cutoff = 0.6,
    min_intersection = 15,
    min_union = 30,
    
    # iteration / safety
    max_iter = 20,
    verbose = TRUE
) {
  stopifnot(inherits(obj, "Seurat"))
  stopifnot(meta_col %in% colnames(obj@meta.data))
  marker_mode <- match.arg(marker_mode)
  
  # --- parse labels like "FB_1", "KC_Spinous_2", "vEC_3"
  #     Only requires suffix _number
  parse_df <- function(labels) {
    labels <- trimws(as.character(labels))
    tibble(label = labels) %>%
      mutate(
        has_num = str_detect(label, "_\\d+$"),
        prefix  = if_else(has_num, str_replace(label, "_\\d+$", ""), NA_character_),
        idx     = if_else(has_num, as.integer(str_extract(label, "\\d+$")), NA_integer_)
      )
  }
  
  # --- Jaccard similarity
  jaccard <- function(a, b) {
    a <- unique(a); b <- unique(b)
    inter <- length(intersect(a, b))
    uni   <- length(union(a, b))
    if (uni == 0) return(list(j = 0, inter = 0, uni = 0))
    list(j = inter / uni, inter = inter, uni = uni)
  }
  
  # --- apply label mapping to meta_col
  apply_label_map <- function(obj, map_vec_named) {
    v <- as.character(obj[[meta_col]][, 1])
    hit <- v %in% names(map_vec_named)
    v[hit] <- unname(map_vec_named[v[hit]])
    obj[[meta_col]] <- v
    obj
  }
  
  # --- renumber within a prefix sequentially by idx order
  renumber_prefix <- function(obj, prefix) {
    labs <- unique(as.character(obj[[meta_col]][, 1]))
    df <- parse_df(labs) %>%
      filter(has_num, prefix == !!prefix) %>%
      arrange(idx)
    
    if (nrow(df) == 0) return(obj)
    
    new_labels <- paste0(prefix, "_", seq_len(nrow(df)))
    rename_map <- setNames(new_labels, df$label)
    apply_label_map(obj, rename_map)
  }
  
  # --- compute marker sets (top N genes) for clusters in groups_use
  get_marker_sets <- function(obj, groups_use) {
    # only keep existing and suffix _number
    all_labels_now <- unique(as.character(obj[[meta_col]][, 1]))
    groups_use <- groups_use[groups_use %in% all_labels_now & grepl("_\\d+$", groups_use)]
    if (length(groups_use) <= 1) {
      return(tibble(cluster = character(), marker_set = list()))
    }
    
    # set identities
    Idents(obj) <- obj[[meta_col]][, 1]
    
    # choose object used to compute markers
    if (marker_mode == "within_prefix") {
      # markers computed within prefix only
      obj_use <- subset(obj, idents = groups_use)
    } else {
      # markers computed globally (vs all other cells)
      obj_use <- obj
    }
    
    mk <- .safe_FindAllMarkers(
      obj_use,
      assay = assay,
      slot = slot,
      only.pos = only.pos,
      logfc.threshold = logfc.threshold,
      min.pct = min.pct,
      test.use = test.use,
      verbose = FALSE
    )
    
    # detect FC column
    fc_col <- intersect(c("avg_log2FC", "avg_logFC", "avg_diff"), colnames(mk))
    if (length(fc_col) == 0) stop("Marker table has no FC column (avg_log2FC/avg_logFC/avg_diff).")
    fc_col <- fc_col[1]
    
    mk_df <- as.data.frame(mk)
    
    # keep only our clusters of interest (important for global mode)
    mk_df <- mk_df[mk_df$cluster %in% groups_use, , drop = FALSE]
    if (nrow(mk_df) == 0) {
      return(tibble(cluster = character(), marker_set = list()))
    }
    
    # split + rank by FC using base R (avoid .data / dplyr mask issues)
    out_list <- split(mk_df, mk_df$cluster)
    out_list <- lapply(out_list, function(df) {
      df <- df[order(-df[[fc_col]]), , drop = FALSE]
      head(unique(df$gene), top_n_markers)
    })
    
    tibble::enframe(out_list, name = "cluster", value = "marker_set")
  }
  
  # ============================================================
  # Iterative merging loop
  # ============================================================
  all_history <- list()
  iter <- 0
  
  repeat {
    iter <- iter + 1
    if (iter > max_iter) {
      warning("Reached max_iter = ", max_iter, ". Stopping to avoid infinite loop.")
      break
    }
    
    labs <- unique(as.character(obj[[meta_col]][, 1]))
    df_labs <- parse_df(labs)
    
    numbered_prefixes <- df_labs %>%
      filter(has_num) %>%
      distinct(prefix) %>%
      pull(prefix)
    
    if (length(numbered_prefixes) == 0) {
      if (verbose) message("[Iter ", iter, "] No suffix _number labels detected. Done.")
      break
    }
    
    if (verbose) {
      message("\n==============================")
      message("[Iter ", iter, "] prefixes: ", paste(numbered_prefixes, collapse = ", "))
      message("[Iter ", iter, "] marker_mode: ", marker_mode)
    }
    
    any_merge_this_iter <- FALSE
    iter_log <- list()
    
    # process each prefix separately
    for (px in numbered_prefixes) {
      
      # ensure current labels for this prefix (sorted by idx)
      cur_labels <- unique(as.character(obj[[meta_col]][, 1]))
      cur_df <- parse_df(cur_labels) %>%
        filter(has_num, prefix == px) %>%
        arrange(idx)
      
      if (nrow(cur_df) <= 1) next
      
      merged_in_prefix <- TRUE
      while (merged_in_prefix) {
        merged_in_prefix <- FALSE
        
        # refresh
        cur_labels <- unique(as.character(obj[[meta_col]][, 1]))
        cur_df <- parse_df(cur_labels) %>%
          filter(has_num, prefix == px) %>%
          arrange(idx)
        if (nrow(cur_df) <= 1) break
        
        cur_labels_px <- cur_df$label
        
        # recompute marker sets (because merges change membership)
        mk_sets <- get_marker_sets(obj, cur_labels_px)
        mk_map <- setNames(mk_sets$marker_set, mk_sets$cluster)
        
        # compare adjacent by idx: 1-2, 2-3, ...
        for (i in seq_len(nrow(cur_df) - 1)) {
          a <- cur_df$label[i]
          b <- cur_df$label[i + 1]
          
          if (!(a %in% names(mk_map)) || !(b %in% names(mk_map))) next
          
          jc <- jaccard(mk_map[[a]], mk_map[[b]])
          
          if (verbose) {
            message(sprintf("  [%s] %s vs %s : J=%.3f (inter=%d, union=%d)",
                            px, a, b, jc$j, jc$inter, jc$uni))
          }
          
          if (is.finite(jc$j) &&
              jc$uni >= min_union &&
              jc$inter >= min_intersection &&
              jc$j >= jaccard_cutoff) {
            
            # merge rule: b -> a (e.g., FB_2 becomes FB_1)
            obj <- apply_label_map(obj, setNames(a, b))
            
            # renumber within prefix after merge (FB_3->FB_2 ...)
            obj <- renumber_prefix(obj, px)
            
            any_merge_this_iter <- TRUE
            merged_in_prefix <- TRUE
            
            iter_log[[length(iter_log) + 1]] <- tibble(
              iter = iter,
              prefix = px,
              merged_from = b,
              merged_into = a,
              jaccard = jc$j,
              intersection = jc$inter,
              union = jc$uni,
              marker_mode = marker_mode
            )
            
            if (verbose) message("    ✅ MERGE: ", b, " -> ", a, " | then renumber within ", px)
            
            # labels changed -> restart scanning
            break
          }
        }
      }
    }
    
    if (length(iter_log) > 0) {
      all_history[[iter]] <- bind_rows(iter_log)
    }
    
    if (!any_merge_this_iter) {
      if (verbose) message("[Iter ", iter, "] No merges. Converged ✅")
      break
    } else {
      if (verbose) message("[Iter ", iter, "] Merges happened; continue...")
    }
  }
  
  history_df <- if (length(all_history) > 0) bind_rows(all_history) else tibble()
  
  # final_counts: base R (avoid rename conflicts)
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
# Example usage
###############################################################################

# 1) 確認 meta 欄位存在
stopifnot("Cell_Type" %in% colnames(seuratObject_Sample@meta.data))

# 2) 執行：合併高度相似的 *_number 群集（遞迴到穩定）
res_merge <- merge_similar_numbered_celltypes(
  obj = seuratObject_Sample,
  meta_col = "Cell_Type",
  assay = "RNA",
  slot  = "data",

  marker_mode = "global",      # ✅ 你要的：vs 全體其他細胞
  # marker_mode = "within_prefix", # ✅ 若要同 lineage 內比對

  top_n_markers = 50,
  jaccard_cutoff = 0.6,
  min_intersection = 15,
  min_union = 30,

  logfc.threshold = 0.25,
  min.pct = 0.1,
  only.pos = TRUE,

  max_iter = 20,
  verbose = TRUE
)
#
# 3) 更新回你的物件
seuratObject_Sample <- res_merge$obj

# 4) 查看合併歷史
print(res_merge$merge_history)

# 5) 查看最後各 Cell_Type 數量
print(res_merge$final_counts)

# 6)（可選）看看最後的 unique labels
print(sort(unique(seuratObject_Sample$Cell_Type)))

