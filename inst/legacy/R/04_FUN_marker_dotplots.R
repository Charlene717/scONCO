#### Visualization ####
# Idents(seuratObject_Sample) <- "seurat_clusters"

# Set_Idents <- "seurat_clusters"
Idents(seuratObject_Sample) <- Set_Idents

# 逐一產生每個細胞類型的 DotPlot
dotplots <- lapply(names(marker_sets), function(cell_type) {
  
  # 1. 去掉重複基因
  feats <- unique(marker_sets[[cell_type]])
  
  # 2. 保留實際存在於 Seurat 物件中的基因
  feats <- feats[feats %in% rownames(seuratObject_Sample)]
  
  # 3. 如果全部被過濾掉就跳過，避免空向量報錯
  if (length(feats) == 0) return(NULL)
  
  DotPlot(
    object    = seuratObject_Sample,
    features  = feats,
    cols      = c("white", "darkred"),
    scale     = FALSE,
    dot.scale = 8
  ) +
    RotatedAxis() +
    labs(title = cell_type) +
    theme(plot.title = element_text(hjust = 0.5, size = 24))
})

# 去掉因為 NULL 而產生的空元素（如果有）
dotplots <- dotplots[!vapply(dotplots, is.null, logical(1))]

if(!require('patchwork')) {install.packages('patchwork'); library(patchwork)}
# 以 patchwork 合併所有 DotPlot；可依資料量調整 ncol 或 nrow
# Plot_combined <- wrap_plots(dotplots, ncol = 3)
Plot_combined <- wrap_plots(dotplots, ncol = Set_Marker_DotPlot_ncol)


#### Export pdf ####
# 輸出成一張 PDF
pdf(paste0(Name_ExportFolder_CTAnnot, "/", Name_Export, "_Classical_markers_Combined_",Set_Idents,".pdf"), 
    width = Set_Marker_DotPlot_width, height = Set_Marker_DotPlot_height)
print(Plot_combined)
dev.off()

png(
  filename = paste0(Name_ExportFolder_CTAnnot, "/", Name_Export, "_Classical_markers_Combined_",Set_Idents,".png"),
  width = Set_Marker_DotPlot_width,        # 寬度 (英吋)
  height = Set_Marker_DotPlot_height,      # 高度 (英吋)
  units = "in",                 # 設定單位為英吋
  res = 300                     # 解析度 (dpi)，可依需求調整
)
print(Plot_combined)
dev.off()


# rm(dotplots, Plot_combined)
