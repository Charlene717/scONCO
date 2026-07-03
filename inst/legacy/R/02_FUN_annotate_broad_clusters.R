library(dplyr)
library(Seurat)
library(tibble)

# 定義函數
annotate_broad_cell_clusters <- function(seurat_object, broad_cell_type_col, seurat_cluster_col) {
  # 如果仍然需要 plyr，請在別處用 library(plyr)，或這裡加上 library(plyr)
  # 請注意載入順序以及衝突函數
  # library(plyr)
  
  # 這裡確保載入 dplyr / tibble
  library(dplyr)
  library(tibble)
  
  # 先取 meta.data，保留細胞名稱到欄位 cell
  meta <- seurat_object@meta.data %>%
    tibble::rownames_to_column(var = "cell")
  
  # 確保這兩個欄位都為 character（防止 factor / numeric 干擾）
  meta[[broad_cell_type_col]] <- as.character(meta[[broad_cell_type_col]])
  meta[[seurat_cluster_col]] <- as.character(meta[[seurat_cluster_col]])
  
  # 計算每個 cell type 出現的 seurat_clusters 數
  celltype_cluster_counts <- meta %>%
    dplyr::group_by(.data[[broad_cell_type_col]]) %>%
    dplyr::summarise(n_clusters = dplyr::n_distinct(.data[[seurat_cluster_col]])) %>%
    dplyr::ungroup()
  
  # 找到出現在多個 cluster 的細胞類型
  multi_cluster_celltypes <- celltype_cluster_counts %>%
    dplyr::filter(n_clusters > 1) %>%
    dplyr::pull(.data[[broad_cell_type_col]])
  
  # 建立 lookup 表：哪些 BroadCellType 對應哪些 cluster
  lookup <- meta %>%
    dplyr::filter(.data[[broad_cell_type_col]] %in% multi_cluster_celltypes) %>%
    dplyr::select(.data[[broad_cell_type_col]], .data[[seurat_cluster_col]]) %>%
    dplyr::distinct() %>%
    dplyr::group_by(.data[[broad_cell_type_col]]) %>%
    dplyr::arrange(.data[[seurat_cluster_col]]) %>%
    dplyr::mutate(suffix = dplyr::row_number()) %>%
    dplyr::ungroup()
  
  # 合併 lookup，並依需要加後綴
  meta_updated <- meta %>%
    dplyr::left_join(
      lookup,
      by = setNames(
        c(broad_cell_type_col, seurat_cluster_col),
        c(broad_cell_type_col, seurat_cluster_col)
      )
    ) %>%
    dplyr::mutate(
      BroadCellTypeAnnot_SeuratClusters = ifelse(
        .data[[broad_cell_type_col]] %in% multi_cluster_celltypes,
        paste0(.data[[broad_cell_type_col]], "_", suffix),
        .data[[broad_cell_type_col]]
      )
    ) %>%
    dplyr::select(-suffix)
  
  # 將 rownames 復原
  meta_updated <- meta_updated %>%
    tibble::column_to_rownames(var = "cell")
  
  # 更新 Seurat 物件並設定新的 identity
  seurat_object@meta.data <- meta_updated
  Idents(seurat_object) <- seurat_object$BroadCellTypeAnnot_SeuratClusters
  
  return(seurat_object)
}


# # 使用示例
# # 假設 `seuratObject_Sample` 是你的 Seurat 對象
# # broad_cell_type_col = "BroadCellTypeAnnot"
# # seurat_cluster_col = "seurat_clusters"
# 
# seuratObject_Sample <- annotate_broad_cell_clusters(
#   seurat_object = seuratObject_Sample,
#   broad_cell_type_col = "BroadCellTypeAnnot",
#   seurat_cluster_col = "seurat_clusters"
# )
# 
# # 驗證修改後的標註
# unique_annotations <- unique(seuratObject_Sample$BroadCellTypeAnnot_SeuratClusters)
# print(unique_annotations)
# 
# # 查看 identities 的分佈
# print(table(Idents(seuratObject_Sample)))
# 
# # 繪製 UMAP 圖
# DimPlot(seuratObject_Sample, group.by = "BroadCellTypeAnnot_SeuratClusters", label = TRUE, reduction = "umap") %>% print()
