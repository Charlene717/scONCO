###############################################################################
#  Pan-cancer / TME marker list (MOUSE) — source: DeepResearch early 20250927
###############################################################################
if (!require("dplyr")) { install.packages("dplyr"); library(dplyr) }
uniq_keep_order <- function(x) x[!duplicated(x)]
marker_sets_ChatGPT <- list(
  "Malignant cell (generic)" = uniq_keep_order(c("Epcam","Krt8","Krt18","Krt19","Cdh1")),
  "T cells"                  = uniq_keep_order(c("Cd3d","Cd3e","Trac","Cd2","Il7r")),
  "CD8+ exhausted T cells"   = uniq_keep_order(c("Pdcd1","Havcr2","Lag3","Tox","Gzmb")),
  "NK cells"                 = uniq_keep_order(c("Ncr1","Klrb1c","Klrd1","Gzmb","Nkg7")),
  "B cells"                  = uniq_keep_order(c("Ms4a1","Cd19","Cd79a","Cd79b","Bank1")),
  "Macrophages"              = uniq_keep_order(c("Adgre1","Cd68","Mrc1","C1qa","Apoe")),
  "Fibroblasts"              = uniq_keep_order(c("Col1a1","Col1a2","Dcn","Lum","Pdgfra")),
  "Endothelial cells"        = uniq_keep_order(c("Pecam1","Cdh5","Cldn5","Egfl7","Flt1"))
)
level_info <- list(
  "Malignant cell (generic)" = c("Malignant cell","Malignant (pan-cancer)","Malignant epithelial cell","Malignant cell (generic)"),
  "T cells"                  = c("Immune cell","Lymphoid cell","T cell","T cells"),
  "CD8+ exhausted T cells"   = c("Immune cell","Lymphoid cell","T cell","CD8+ exhausted T cells"),
  "NK cells"                 = c("Immune cell","Lymphoid cell","Natural killer cell","NK cells"),
  "B cells"                  = c("Immune cell","Lymphoid cell","B cell","B cells"),
  "Macrophages"              = c("Immune cell","Myeloid cell","Macrophage","Macrophages"),
  "Fibroblasts"              = c("Stromal cell","Fibroblast","Fibroblast","Fibroblasts"),
  "Endothelial cells"        = c("Endothelial cell","Vascular endothelial cell","Endothelial cell","Endothelial cells")
)
build_scMRMA <- function(marker_sets, level_info) {
  do.call(rbind, lapply(names(marker_sets), function(ct) {
    g <- uniq_keep_order(marker_sets[[ct]]); lv <- level_info[[ct]]
    data.frame(gene=g, Level1=lv[1], Level2=lv[2], Level3=lv[3], Level4_Abb=lv[4], stringsAsFactors=FALSE)
  }))
}
markers_df <- build_scMRMA(marker_sets_ChatGPT, level_info)
cat("[DeepResearch early mouse] markers_df:", nrow(markers_df), "rows\n")
