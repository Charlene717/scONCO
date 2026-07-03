###############################################################################
#  Pan-cancer / TME marker list — source: ChatGPT-5.1 Instant V2
#  scONCO alternative AI-source curation (compact, valid reconstruction).
#  Defines markers_df via marker_sets_ChatGPT + level_info + build_scMRMA,
#  matching DB_pancancer_human_v1.0.R schema (gene|Level1|Level2|Level3|Level4_Abb).
#  NOTE: compact reconstruction of an interrupted earlier draft; extend as needed.
###############################################################################
if (!require("dplyr")) { install.packages("dplyr"); library(dplyr) }
uniq_keep_order <- function(x) x[!duplicated(x)]

marker_sets_ChatGPT <- list(
  "Malignant cell (generic)" = uniq_keep_order(c("EPCAM","KRT8","KRT18","KRT19","CDH1","MUC1")),
  "Cycling malignant cell"   = uniq_keep_order(c("MKI67","TOP2A","UBE2C","BIRC5","CENPF")),
  "Malignant melanoma cell"  = uniq_keep_order(c("MLANA","PMEL","TYR","MITF","SOX10")),
  "Epithelial cells"         = uniq_keep_order(c("EPCAM","CDH1","KRT8","KRT18","KRT19")),
  "T cells"                  = uniq_keep_order(c("CD3D","CD3E","TRAC","CD2","IL7R")),
  "CD8+ T cells"             = uniq_keep_order(c("CD8A","CD8B","GZMK","NKG7","CCL5")),
  "CD8+ exhausted T cells"   = uniq_keep_order(c("PDCD1","HAVCR2","LAG3","TIGIT","TOX")),
  "Regulatory T cells"       = uniq_keep_order(c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18")),
  "NK cells"                 = uniq_keep_order(c("NCAM1","NKG7","GNLY","KLRD1","NCR1")),
  "B cells"                  = uniq_keep_order(c("MS4A1","CD19","CD79A","CD79B","BANK1")),
  "Plasma cells"             = uniq_keep_order(c("MZB1","JCHAIN","XBP1","SDC1","DERL3")),
  "Monocytes"                = uniq_keep_order(c("CD14","FCN1","VCAN","S100A8","S100A9")),
  "Macrophages"              = uniq_keep_order(c("CD68","CD163","C1QA","C1QB","APOE")),
  "TAM SPP1 (angiogenic)"    = uniq_keep_order(c("SPP1","MARCO","FN1","MMP9","INHBA")),
  "Conventional dendritic cells type 1 (cDC1)" = uniq_keep_order(c("CLEC9A","XCR1","BATF3","IRF8")),
  "Plasmacytoid dendritic cells" = uniq_keep_order(c("LILRA4","CLEC4C","IL3RA","GZMB","IRF7")),
  "Mast cells"               = uniq_keep_order(c("TPSAB1","TPSB2","CPA3","KIT","MS4A2")),
  "Fibroblasts"              = uniq_keep_order(c("COL1A1","COL1A2","DCN","LUM","PDGFRA")),
  "Myofibroblastic CAF (myCAF)" = uniq_keep_order(c("ACTA2","TAGLN","POSTN","FAP","CTHRC1")),
  "Inflammatory CAF (iCAF)"  = uniq_keep_order(c("IL6","CXCL12","CXCL14","CFD","DPT")),
  "Pericytes"                = uniq_keep_order(c("RGS5","PDGFRB","NOTCH3","CSPG4","KCNJ8")),
  "Endothelial cells"        = uniq_keep_order(c("PECAM1","VWF","CDH5","CLDN5","CD34")),
  "Lymphatic endothelial cells" = uniq_keep_order(c("PROX1","LYVE1","PDPN","FLT4","CCL21"))
)

level_info <- list(
  "Malignant cell (generic)" = c("Malignant cell","Malignant (pan-cancer)","Malignant epithelial cell","Malignant cell (generic)"),
  "Cycling malignant cell"   = c("Malignant cell","Malignant (pan-cancer)","Cycling malignant cell","Cycling malignant cell"),
  "Malignant melanoma cell"  = c("Malignant cell","Malignant (melanoma)","Melanoma cell","Malignant melanoma cell"),
  "Epithelial cells"         = c("Epithelial cell","Normal epithelial cell","Epithelial cell","Epithelial cells"),
  "T cells"                  = c("Immune cell","Lymphoid cell","T cell","T cells"),
  "CD8+ T cells"             = c("Immune cell","Lymphoid cell","T cell","CD8+ T cells"),
  "CD8+ exhausted T cells"   = c("Immune cell","Lymphoid cell","T cell","CD8+ exhausted T cells"),
  "Regulatory T cells"       = c("Immune cell","Lymphoid cell","T cell","Regulatory T cells"),
  "NK cells"                 = c("Immune cell","Lymphoid cell","Natural killer cell","NK cells"),
  "B cells"                  = c("Immune cell","Lymphoid cell","B cell","B cells"),
  "Plasma cells"             = c("Immune cell","Lymphoid cell","Plasma cell","Plasma cells"),
  "Monocytes"                = c("Immune cell","Myeloid cell","Monocyte","Monocytes"),
  "Macrophages"              = c("Immune cell","Myeloid cell","Macrophage","Macrophages"),
  "TAM SPP1 (angiogenic)"    = c("Immune cell","Myeloid cell","Tumour-associated macrophage","TAM SPP1 (angiogenic)"),
  "Conventional dendritic cells type 1 (cDC1)" = c("Immune cell","Myeloid cell","Dendritic cell","Conventional dendritic cells type 1 (cDC1)"),
  "Plasmacytoid dendritic cells" = c("Immune cell","Myeloid cell","Dendritic cell","Plasmacytoid dendritic cells"),
  "Mast cells"               = c("Immune cell","Myeloid cell","Mast cell","Mast cells"),
  "Fibroblasts"              = c("Stromal cell","Fibroblast","Fibroblast","Fibroblasts"),
  "Myofibroblastic CAF (myCAF)" = c("Stromal cell","Fibroblast","Cancer-associated fibroblast","Myofibroblastic CAF (myCAF)"),
  "Inflammatory CAF (iCAF)"  = c("Stromal cell","Fibroblast","Cancer-associated fibroblast","Inflammatory CAF (iCAF)"),
  "Pericytes"                = c("Stromal cell","Perivascular mural cell","Pericyte","Pericytes"),
  "Endothelial cells"        = c("Endothelial cell","Vascular endothelial cell","Endothelial cell","Endothelial cells"),
  "Lymphatic endothelial cells" = c("Endothelial cell","Lymphatic endothelial cell","Lymphatic endothelial cell","Lymphatic endothelial cells")
)

build_scMRMA <- function(marker_sets, level_info) {
  do.call(rbind, lapply(names(marker_sets), function(ct) {
    genes <- uniq_keep_order(marker_sets[[ct]]); lv <- level_info[[ct]]
    data.frame(gene = genes, Level1 = lv[1], Level2 = lv[2],
               Level3 = lv[3], Level4_Abb = lv[4], stringsAsFactors = FALSE)
  }))
}
marker_df_scMRMA_ChatGPT <- build_scMRMA(marker_sets_ChatGPT, level_info)
markers_df <- marker_df_scMRMA_ChatGPT
cat("[ChatGPT-5.1 Instant V2] markers_df:", nrow(markers_df), "rows /",
    length(unique(markers_df$Level4_Abb)), "cell types\n")
