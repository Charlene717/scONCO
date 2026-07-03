###############################################################################
## ✅ FIXED Reclustering refinement block (auto-detect graph.name)
## - Works with Seurat v5 graph naming (reclust / reclust_snn / reclust.nn)
###############################################################################

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(scMRMA)
})

DefaultAssay(seuratObject_Sample) <- "RNA"

###############################################################################
## Settings (same as before)
###############################################################################
RECLUST_BASE_RES   <- 0.4
RECLUST_RES_FACTOR <- 1.5
RECLUST_RES_MAX    <- 6.0
RECLUST_MAX_ITERS  <- 8
RECLUST_MIN_CELLS_PER_SUBCLUSTER <- 15

PCA_NPCS <- 30
PCA_DIMS <- 1:30
NN_KPARAM <- 20
GRAPH_PREFIX <- "reclust"

Set_scMRMA_P   <- get0("Set_scMRMA_P", ifnotfound = 0.05)
Set_SubClust_k <- get0("Set_SubClust_k", ifnotfound = 20)
MIN_FEATURES_FOR_SCMRMA <- 200

COL_ORIG_CLU      <- "seurat_clusters_orig"
COL_RECLUST_TAG   <- "scONCO_reclust_tag"
COL_RECLUST_CLU   <- "scONCO_reclust_clusters"
COL_L4_FINAL      <- "scONCO_L4Abb_refined"

###############################################################################
## Helpers (minimal needed; assumes you already defined these earlier, but re-define safely)
###############################################################################
`%||%` <- function(a, b) if (!is.null(a)) a else b

.get_layer <- function(seu, assay, layer) {
  SeuratObject::GetAssayData(seu, assay = assay, layer = layer)
}

.pick_layer <- function(seu, assay = "RNA", min_features = 200) {
  m_data <- tryCatch(.get_layer(seu, assay, "data"), error = function(e) NULL)
  if (!is.null(m_data) && nrow(m_data) >= min_features && ncol(m_data) == ncol(seu)) return("data")
  m_counts <- tryCatch(.get_layer(seu, assay, "counts"), error = function(e) NULL)
  if (!is.null(m_counts) && nrow(m_counts) >= min_features && ncol(m_counts) == ncol(seu)) return("counts")
  "data"
}

.add_tmp_clusters_0based <- function(seu, from_col, to_col = "scMRMA_tmp_clusters") {
  v <- as.character(seu[[from_col]][, 1])
  u <- sort(unique(v))
  map <- setNames(as.character(seq_along(u) - 1L), u)
  tmp <- unname(map[v])
  tmp <- factor(tmp, levels = as.character(0:(length(u) - 1L)))
  seu[[to_col]] <- tmp
  seu
}

.get_clusters_aligned <- function(seu, cluster_col) {
  clu <- as.character(seu[[cluster_col]][, 1])
  names(clu) <- colnames(seu)
  clu <- clu[colnames(seu)]
  if (!identical(names(clu), colnames(seu))) stop("Cluster names not aligned to cellnames.")
  if (length(clu) != ncol(seu)) stop("Cluster length != nCells.")
  as.factor(clu)
}

.check_unique_by_reclust <- function(reclust_ids, labels) {
  df <- data.frame(rc = as.character(reclust_ids), lb = as.character(labels), stringsAsFactors = FALSE)
  ok <- TRUE
  maj_map <- setNames(rep(NA_character_, length(unique(df$rc))), unique(df$rc))
  for (cc in unique(df$rc)) {
    sub <- df[df$rc == cc & !is.na(df$lb) & df$lb != "", , drop=FALSE]
    if (nrow(sub) == 0) { ok <- FALSE; next }
    tt <- sort(table(sub$lb), decreasing = TRUE)
    if (length(tt) >= 2 && tt[1] == tt[2]) { ok <- FALSE; next }
    maj_map[cc] <- names(tt)[1]
  }
  list(ok = ok, majority_map = maj_map)
}

.run_scMRMA_cell_level <- function(seu_sub, selfDB, cluster_col, p, k, assay="RNA", layer=NULL, verbose=TRUE) {
  DefaultAssay(seu_sub) <- assay
  layer_use <- layer %||% .pick_layer(seu_sub, assay = assay, min_features = MIN_FEATURES_FOR_SCMRMA)
  m <- tryCatch(.get_layer(seu_sub, assay, layer_use), error=function(e) NULL)
  if (is.null(m) || nrow(m) < MIN_FEATURES_FOR_SCMRMA) return(NULL)
  
  genes_present <- intersect(selfDB$gene, rownames(seu_sub))
  selfDB2 <- selfDB[selfDB$gene %in% genes_present, , drop=FALSE]
  if (nrow(selfDB2) < 5 || length(unique(selfDB2$LevelS)) < 2) return(NULL)
  
  clu <- .get_clusters_aligned(seu_sub, cluster_col)
  
  res <- tryCatch({
    scMRMA(
      input = seu_sub,
      p = p,
      normalizedData = FALSE,
      selfDB = selfDB2,
      selfClusters = clu,
      k = k
    )
  }, error=function(e) {
    if (verbose) message("scMRMA ERROR: ", conditionMessage(e))
    return(NULL)
  })
  if (is.null(res)) return(NULL)
  
  ann <- res[["multiR"]][["annotationResult"]]
  if (length(ann) != ncol(seu_sub)) return(NULL)
  as.character(ann)
}

# ✅ NEW: pick graph name after FindNeighbors
.pick_graph_for_clustering <- function(seu, graph_prefix = "reclust") {
  gnames <- names(seu@graphs)
  if (length(gnames) == 0) stop("No graphs found in object after FindNeighbors().")
  
  # Prefer anything that starts with prefix and contains "snn"
  cand1 <- gnames[grepl(paste0("^", graph_prefix), gnames) & grepl("snn", gnames)]
  if (length(cand1) > 0) return(cand1[1])
  
  # Otherwise prefer exact prefix (your current case: "reclust")
  if (graph_prefix %in% gnames) return(graph_prefix)
  
  # Otherwise any graph starting with prefix
  cand2 <- gnames[grepl(paste0("^", graph_prefix), gnames)]
  if (length(cand2) > 0) return(cand2[1])
  
  # As last resort, use any *_snn (integrated_snn etc.)
  cand3 <- gnames[grepl("snn", gnames)]
  if (length(cand3) > 0) return(cand3[1])
  
  # Last: first available
  gnames[1]
}

###############################################################################
## Prepare DB (Level4_Abb markers)
###############################################################################
selfDB_L4 <- markers_df %>% dplyr::select(gene, Level4_Abb) %>% dplyr::distinct()
colnames(selfDB_L4)[2] <- "LevelS"

###############################################################################
## You must already have amb_targets$orig_clu from your ambiguous detection step.
## If not, define it here as a vector of clusters you want to refine, e.g.:
## amb_targets <- data.frame(orig_clu = c("1","5"))
###############################################################################
stopifnot(exists("amb_targets"))
stopifnot("orig_clu" %in% colnames(amb_targets))

# init refined labels if not exists
if (!(COL_L4_FINAL %in% colnames(seuratObject_Sample@meta.data))) {
  seuratObject_Sample[[COL_L4_FINAL]] <- seuratObject_Sample$scONCO_L4Abb
}
if (!(COL_RECLUST_TAG %in% colnames(seuratObject_Sample@meta.data))) {
  seuratObject_Sample[[COL_RECLUST_TAG]] <- NA_character_
}
if (!(COL_RECLUST_CLU %in% colnames(seuratObject_Sample@meta.data))) {
  seuratObject_Sample[[COL_RECLUST_CLU]] <- NA_character_
}
if (!(COL_ORIG_CLU %in% colnames(seuratObject_Sample@meta.data))) {
  seuratObject_Sample[[COL_ORIG_CLU]] <- seuratObject_Sample$seurat_clusters
}

###############################################################################
## Reclustering refinement loop
###############################################################################
for (orig in amb_targets$orig_clu) {
  message("\n=== Refining original cluster: ", orig, " ===")
  
  cells_orig <- rownames(seuratObject_Sample@meta.data)[as.character(seuratObject_Sample@meta.data[[COL_ORIG_CLU]]) == orig]
  if (length(cells_orig) < 2 * RECLUST_MIN_CELLS_PER_SUBCLUSTER) {
    message("Skip: too few cells in this cluster.")
    next
  }
  
  sub <- subset(seuratObject_Sample, cells = cells_orig)
  DefaultAssay(sub) <- "RNA"
  
  # ✅ avoid "Different features in new layer data than already exists for scale.data"
  # wipe scale.data and redo processing cleanly
  try({
    sub[["RNA"]]@scale.data <- matrix(numeric(0), nrow = 0, ncol = 0)
  }, silent = TRUE)
  
  sub <- NormalizeData(sub, verbose = FALSE)
  sub <- FindVariableFeatures(sub, verbose = FALSE)
  sub <- ScaleData(sub, verbose = FALSE)
  sub <- RunPCA(sub, npcs = PCA_NPCS, verbose = FALSE)
  
  res_now <- RECLUST_BASE_RES
  iter <- 0
  best_majority_map <- NULL
  best_subcluster_ids <- NULL
  best_iter_tag <- NULL
  
  while (iter < RECLUST_MAX_ITERS && res_now <= RECLUST_RES_MAX) {
    iter <- iter + 1
    message("  Iter ", iter, " | resolution=", signif(res_now, 3))
    
    sub <- FindNeighbors(
      sub,
      reduction = "pca",
      dims = PCA_DIMS,
      k.param = NN_KPARAM,
      graph.name = GRAPH_PREFIX,
      verbose = FALSE
    )
    
    graph_name_use <- .pick_graph_for_clustering(sub, graph_prefix = GRAPH_PREFIX)
    message("    Using graph.name = ", graph_name_use)
    
    sub <- FindClusters(
      sub,
      graph.name = graph_name_use,
      resolution = res_now,
      verbose = FALSE
    )
    
    tabc <- table(sub$seurat_clusters)
    if (any(tabc < RECLUST_MIN_CELLS_PER_SUBCLUSTER)) {
      message("    Warning: tiny subclusters present (min=", RECLUST_MIN_CELLS_PER_SUBCLUSTER, ").")
    }
    
    sub <- .add_tmp_clusters_0based(sub, from_col = "seurat_clusters", to_col = "scMRMA_tmp_clusters")
    
    lab_cell <- .run_scMRMA_cell_level(
      seu_sub = sub,
      selfDB = selfDB_L4,
      cluster_col = "scMRMA_tmp_clusters",
      p = Set_scMRMA_P,
      k = Set_SubClust_k,
      assay = "RNA",
      layer = NULL,
      verbose = FALSE
    )
    
    if (is.null(lab_cell)) {
      message("    scMRMA failed/NULL. Increase resolution...")
      res_now <- res_now * RECLUST_RES_FACTOR
      next
    }
    
    chk <- .check_unique_by_reclust(sub$seurat_clusters, lab_cell)
    
    best_majority_map <- chk$majority_map
    best_subcluster_ids <- as.character(sub$seurat_clusters)
    names(best_subcluster_ids) <- colnames(sub)
    best_iter_tag <- paste0("orig_", orig, "_iter", iter, "_res", signif(res_now, 3))
    
    if (chk$ok) {
      message("    ✅ Unique labeling achieved.")
      break
    } else {
      message("    Not unique yet. Increase resolution...")
      res_now <- res_now * RECLUST_RES_FACTOR
    }
  }
  
  if (is.null(best_majority_map) || is.null(best_subcluster_ids)) {
    message("  Failed to refine: no usable scMRMA result.")
    next
  }
  
  final_l4 <- unname(best_majority_map[best_subcluster_ids])
  main_idx <- match(names(best_subcluster_ids), colnames(seuratObject_Sample))
  
  fallback_prev <- seuratObject_Sample$scONCO_L4Abb[main_idx]
  final_l4[is.na(final_l4) | final_l4 == ""] <- fallback_prev[is.na(final_l4) | final_l4 == ""]
  
  seuratObject_Sample[[COL_L4_FINAL]][main_idx, 1] <- final_l4
  seuratObject_Sample[[COL_RECLUST_TAG]][main_idx, 1] <- best_iter_tag
  seuratObject_Sample[[COL_RECLUST_CLU]][main_idx, 1] <- paste0("orig", orig, "_sub", best_subcluster_ids)
  
  message("  Applied refined labels to ", length(main_idx), " cells.")
}

message("\n✅ Refinement finished. New label column: ", COL_L4_FINAL)

###############################################################################
## 4) Quick check plots
###############################################################################
message("\n✅ Refinement finished.")
message("New refined label column: ", COL_L4_FINAL)

# UMAP: before vs after
p_before <- DimPlot(seuratObject_Sample, reduction = "umap", group.by = "scONCO_L4Abb", label = TRUE) +
  ggplot2::ggtitle("Before refinement: scONCO_L4Abb")
p_after  <- DimPlot(seuratObject_Sample, reduction = "umap", group.by = COL_L4_FINAL, label = TRUE) +
  ggplot2::ggtitle("After refinement: scONCO_L4Abb_refined")

p_before + p_after

