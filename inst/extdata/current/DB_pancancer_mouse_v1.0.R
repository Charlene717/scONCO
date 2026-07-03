###############################################################################
#  Pan-Cancer Marker Reference (MOUSE)  ->  scMRMA format
#  (gene | Level1 | Level2 | Level3 | Level4_Abb)
#
#  scONCO v1.0 — mouse pan-cancer tumour-microenvironment (TME) marker DB.
#  Compact core-TME companion to DB_pancancer_human_v1.0.R, using MGI symbols.
#  Intended for mouse syngeneic / GEMM tumour models. Malignant cell types are
#  kept generic (tumour-model-specific signatures should be added per study).
###############################################################################

if (!require("dplyr")) { install.packages("dplyr"); library(dplyr) }

uniq_keep_order <- function(x) x[!duplicated(x)]

marker_sets_ChatGPT <- list(
  "Malignant cell (generic)" = uniq_keep_order(c(
    "Epcam","Krt8","Krt18","Krt19","Cdh1","Muc1","Krt7","Cd24a","Elf3","Spint2")),
  "Cycling malignant cell" = uniq_keep_order(c(
    "Mki67","Top2a","Cdk1","Ube2c","Ccnb1","Birc5","Cenpf","Pcna","Stmn1","Tyms")),
  "Epithelial cells" = uniq_keep_order(c(
    "Epcam","Cdh1","Krt8","Krt18","Krt19","Krt7","Cldn3","Cldn4","Cd24a")),

  "T cells" = uniq_keep_order(c(
    "Cd3d","Cd3e","Cd3g","Trac","Cd2","Cd5","Il7r","Cd4","Cd8a","Cd8b1","Lck")),
  "CD4+ T cells" = uniq_keep_order(c(
    "Cd4","Il7r","Cd40lg","Ccr7","Tcf7","Lef1","Sell","Foxp3","Rora","Ltb")),
  "CD8+ T cells" = uniq_keep_order(c(
    "Cd8a","Cd8b1","Gzmk","Gzma","Nkg7","Ccl5","Gzmb","Cd3d","Eomes","Cxcr3")),
  "CD8+ exhausted T cells" = uniq_keep_order(c(
    "Pdcd1","Havcr2","Lag3","Tigit","Ctla4","Tox","Entpd1","Cxcl13","Gzmb")),
  "Regulatory T cells" = uniq_keep_order(c(
    "Foxp3","Il2ra","Ctla4","Tnfrsf18","Ikzf2","Tigit","Tnfrsf4","Ccr8","Cd4")),
  "NK cells" = uniq_keep_order(c(
    "Ncr1","Klrb1c","Klrd1","Gzmb","Prf1","Nkg7","Klrk1","Eomes","Tyrobp")),

  "B cells" = uniq_keep_order(c(
    "Ms4a1","Cd19","Cd79a","Cd79b","Bank1","Ighd","Ighm","Cd22","Pax5")),
  "Plasma cells" = uniq_keep_order(c(
    "Mzb1","Jchain","Xbp1","Sdc1","Prdm1","Derl3","Tnfrsf17","Ighg1","Cd38")),

  "Monocytes" = uniq_keep_order(c(
    "Ly6c2","Fcn1","Vcan","S100a8","S100a9","Lyz2","Csf1r","Ccr2","Plac8")),
  "Macrophages" = uniq_keep_order(c(
    "Adgre1","Cd68","Mrc1","C1qa","C1qb","C1qc","Apoe","Csf1r","Fcgr1")),
  "TAM SPP1 (angiogenic)" = uniq_keep_order(c(
    "Spp1","Arg1","Fn1","Mmp12","Vegfa","Ctsb","Slc2a1","Cd68")),
  "TAM C1QC (immunoregulatory)" = uniq_keep_order(c(
    "C1qa","C1qb","C1qc","Apoe","Apoc1","Trem2","Gpnmb","Folr2","Cd163")),
  "Conventional dendritic cells type 1 (cDC1)" = uniq_keep_order(c(
    "Xcr1","Clec9a","Batf3","Irf8","Cadm1","Naaa","Itgae")),
  "Conventional dendritic cells type 2 (cDC2)" = uniq_keep_order(c(
    "Cd1c","Itgax","Clec10a","Sirpa","Cd209a","H2-Ab1")),
  "Mature regulatory DC (mregDC/LAMP3)" = uniq_keep_order(c(
    "Ccr7","Fscn1","Lamp3","Ccl22","Cd40","Relb","Il4i1")),
  "Plasmacytoid dendritic cells" = uniq_keep_order(c(
    "Siglech","Bst2","Irf7","Tcf4","Ly6d","Klk1","Spib")),
  "Neutrophils (TAN)" = uniq_keep_order(c(
    "S100a8","S100a9","Retnlg","Csf3r","Mmp9","Ly6g","Cxcr2","G0s2")),
  "Mast cells" = uniq_keep_order(c(
    "Cpa3","Mcpt4","Cma1","Kit","Ms4a2","Gata2","Tpsb2")),

  "Fibroblasts" = uniq_keep_order(c(
    "Col1a1","Col1a2","Col3a1","Dcn","Lum","Pdgfra","Mgp","Gsn","Sfrp2")),
  "Myofibroblastic CAF (myCAF)" = uniq_keep_order(c(
    "Acta2","Tagln","Myl9","Postn","Col12a1","Thbs2","Fap","Inhba")),
  "Inflammatory CAF (iCAF)" = uniq_keep_order(c(
    "Il6","Cxcl12","Cxcl1","Has1","Pdgfra","Dpt","C3","Clec3b")),
  "Pericytes" = uniq_keep_order(c(
    "Rgs5","Pdgfrb","Notch3","Cspg4","Acta2","Kcnj8","Higd1b","Cox4i2")),

  "Endothelial cells" = uniq_keep_order(c(
    "Pecam1","Cdh5","Cldn5","Egfl7","Ramp2","Flt1","Kdr","Aqp1","Emcn")),
  "Lymphatic endothelial cells" = uniq_keep_order(c(
    "Prox1","Lyve1","Pdpn","Flt4","Ccl21a","Mmrn1","Reln")),

  "Erythrocytes" = uniq_keep_order(c(
    "Hbb-bs","Hba-a1","Alas2","Slc4a1","Gypa","Bnip3l"))
)

level_info <- list(
  "Malignant cell (generic)" = c("Malignant cell","Malignant (pan-cancer)","Malignant epithelial / tumour cell","Malignant cell (generic)"),
  "Cycling malignant cell"   = c("Malignant cell","Malignant (pan-cancer)","Cycling malignant cell","Cycling malignant cell"),
  "Epithelial cells"         = c("Epithelial cell","Normal epithelial cell","Epithelial cell","Epithelial cells"),
  "T cells"                  = c("Immune cell","Lymphoid cell","T cell","T cells"),
  "CD4+ T cells"             = c("Immune cell","Lymphoid cell","T cell","CD4+ T cells"),
  "CD8+ T cells"             = c("Immune cell","Lymphoid cell","T cell","CD8+ T cells"),
  "CD8+ exhausted T cells"   = c("Immune cell","Lymphoid cell","T cell","CD8+ exhausted T cells"),
  "Regulatory T cells"       = c("Immune cell","Lymphoid cell","T cell","Regulatory T cells"),
  "NK cells"                 = c("Immune cell","Lymphoid cell","Natural killer cell","NK cells"),
  "B cells"                  = c("Immune cell","Lymphoid cell","B cell","B cells"),
  "Plasma cells"             = c("Immune cell","Lymphoid cell","Plasma cell","Plasma cells"),
  "Monocytes"                = c("Immune cell","Myeloid cell","Monocyte","Monocytes"),
  "Macrophages"              = c("Immune cell","Myeloid cell","Macrophage","Macrophages"),
  "TAM SPP1 (angiogenic)"    = c("Immune cell","Myeloid cell","Tumour-associated macrophage","TAM SPP1 (angiogenic)"),
  "TAM C1QC (immunoregulatory)" = c("Immune cell","Myeloid cell","Tumour-associated macrophage","TAM C1QC (immunoregulatory)"),
  "Conventional dendritic cells type 1 (cDC1)" = c("Immune cell","Myeloid cell","Dendritic cell","Conventional dendritic cells type 1 (cDC1)"),
  "Conventional dendritic cells type 2 (cDC2)" = c("Immune cell","Myeloid cell","Dendritic cell","Conventional dendritic cells type 2 (cDC2)"),
  "Mature regulatory DC (mregDC/LAMP3)" = c("Immune cell","Myeloid cell","Dendritic cell","Mature regulatory DC (mregDC/LAMP3)"),
  "Plasmacytoid dendritic cells" = c("Immune cell","Myeloid cell","Dendritic cell","Plasmacytoid dendritic cells"),
  "Neutrophils (TAN)"        = c("Immune cell","Myeloid cell","Granulocyte","Neutrophils (TAN)"),
  "Mast cells"               = c("Immune cell","Myeloid cell","Mast cell","Mast cells"),
  "Fibroblasts"              = c("Stromal cell","Fibroblast","Fibroblast","Fibroblasts"),
  "Myofibroblastic CAF (myCAF)" = c("Stromal cell","Fibroblast","Cancer-associated fibroblast","Myofibroblastic CAF (myCAF)"),
  "Inflammatory CAF (iCAF)"  = c("Stromal cell","Fibroblast","Cancer-associated fibroblast","Inflammatory CAF (iCAF)"),
  "Pericytes"                = c("Stromal cell","Perivascular mural cell","Pericyte","Pericytes"),
  "Endothelial cells"        = c("Endothelial cell","Vascular endothelial cell","Endothelial cell","Endothelial cells"),
  "Lymphatic endothelial cells" = c("Endothelial cell","Lymphatic endothelial cell","Lymphatic endothelial cell","Lymphatic endothelial cells"),
  "Erythrocytes"             = c("Other","Erythroid cell","Erythrocyte","Erythrocytes")
)

build_scMRMA <- function(marker_sets, level_info) {
  missing <- setdiff(names(marker_sets), names(level_info))
  if (length(missing)) stop("Missing level_info for: ", paste(missing, collapse = ", "))
  do.call(rbind, lapply(names(marker_sets), function(ct) {
    genes <- uniq_keep_order(marker_sets[[ct]]); lv <- level_info[[ct]]
    data.frame(gene = genes, Level1 = lv[1], Level2 = lv[2],
               Level3 = lv[3], Level4_Abb = lv[4], stringsAsFactors = FALSE)
  }))
}

marker_df_scMRMA_ChatGPT <- build_scMRMA(marker_sets_ChatGPT, level_info)

level4_abbrev_map <- c(
  "Malignant cell (generic)" = "Mal", "Cycling malignant cell" = "Mal_Cyc",
  "Epithelial cells" = "Epi",
  "T cells" = "Tcell", "CD4+ T cells" = "CD4T", "CD8+ T cells" = "CD8T",
  "CD8+ exhausted T cells" = "CD8_Tex", "Regulatory T cells" = "Treg", "NK cells" = "NK",
  "B cells" = "Bcell", "Plasma cells" = "Plasma",
  "Monocytes" = "Mono", "Macrophages" = "Mac",
  "TAM SPP1 (angiogenic)" = "TAM_SPP1", "TAM C1QC (immunoregulatory)" = "TAM_C1QC",
  "Conventional dendritic cells type 1 (cDC1)" = "cDC1",
  "Conventional dendritic cells type 2 (cDC2)" = "cDC2",
  "Mature regulatory DC (mregDC/LAMP3)" = "mregDC",
  "Plasmacytoid dendritic cells" = "pDC", "Neutrophils (TAN)" = "Neu", "Mast cells" = "Mast",
  "Fibroblasts" = "FB", "Myofibroblastic CAF (myCAF)" = "myCAF",
  "Inflammatory CAF (iCAF)" = "iCAF", "Pericytes" = "Peri",
  "Endothelial cells" = "EC", "Lymphatic endothelial cells" = "lEC",
  "Erythrocytes" = "RBC"
)

old_abb <- marker_df_scMRMA_ChatGPT$Level4_Abb
marker_df_scMRMA_ChatGPT$Level4_Abb <- ifelse(
  old_abb %in% names(level4_abbrev_map),
  unname(level4_abbrev_map[old_abb]), old_abb)

cat("scONCO pan-cancer DB v1.0 (mouse)\n")
cat("Cell types: ", length(unique(marker_df_scMRMA_ChatGPT$Level4_Abb)),
    " | rows: ", nrow(marker_df_scMRMA_ChatGPT), "\n", sep = "")

markers_df <- marker_df_scMRMA_ChatGPT
