#### Expert-curated Marker DB — Pan-Cancer / TME (Human) ####
## Dr. 戴揚紘 (Yang-Hong Dai) — radiation oncology / bioinformatics curator
## 泛癌腫瘤微環境 (TME) expert-curated marker DB.
## Style C: create_marker_rows(genes, level1, level2, level3, level4) + bind_rows -> DBMarker_Expert

library(dplyr)
library(tibble)

create_marker_rows <- function(genes, level1, level2, level3, level4) {
  tibble(gene = genes, Level1 = level1, Level2 = level2, Level3 = level3, Level4 = level4)
}

DBMarker_Expert <- bind_rows(
  create_marker_rows(c("EPCAM","KRT8","KRT18","KRT19","CDH1","MUC1"),
    "Malignant cell","Malignant (pan-cancer)","Malignant epithelial cell","Malignant cell (generic)"),
  create_marker_rows(c("MKI67","TOP2A","UBE2C","BIRC5","CENPF"),
    "Malignant cell","Malignant (pan-cancer)","Cycling malignant cell","Cycling malignant cell"),
  create_marker_rows(c("VIM","ZEB1","SNAI2","TWIST1","FN1"),
    "Malignant cell","Malignant (pan-cancer)","EMT-like malignant cell","EMT-like malignant cell"),
  create_marker_rows(c("MLANA","PMEL","TYR","DCT","MITF","SOX10"),
    "Malignant cell","Malignant (melanoma)","Melanoma cell","Malignant melanoma cell"),
  create_marker_rows(c("NKX2-1","NAPSA","SFTPC","SFTPB","CEACAM5"),
    "Malignant cell","Malignant (LUAD)","Lung adenocarcinoma cell","Malignant lung adenocarcinoma cell"),
  create_marker_rows(c("ESR1","ERBB2","GATA3","FOXA1","TFF1"),
    "Malignant cell","Malignant (BRCA)","Breast carcinoma cell","Malignant breast cell"),
  create_marker_rows(c("ALB","APOA1","TTR","SERPINA1","GPC3"),
    "Malignant cell","Malignant (HCC)","Hepatocellular carcinoma cell","Malignant hepatocellular cell"),
  create_marker_rows(c("CDX2","KRT20","CEACAM5","TFF3","OLFM4"),
    "Malignant cell","Malignant (CRC)","Colorectal carcinoma cell","Malignant colorectal cell"),
  create_marker_rows(c("GFAP","SOX2","OLIG2","PTPRZ1","EGFR"),
    "Malignant cell","Malignant (GBM/glioma)","Glioma / glioblastoma cell","Malignant glioma cell"),
  create_marker_rows(c("EPCAM","CDH1","KRT8","KRT18","KRT19"),
    "Epithelial cell","Normal epithelial cell","Epithelial cell","Epithelial cells"),
  create_marker_rows(c("CD3D","CD3E","TRAC","CD2","IL7R"),
    "Immune cell","Lymphoid cell","T cell","T cells"),
  create_marker_rows(c("CD4","IL7R","CD40LG","CCR7","TCF7"),
    "Immune cell","Lymphoid cell","T cell","CD4+ T cells"),
  create_marker_rows(c("CD8A","CD8B","GZMK","NKG7","CCL5"),
    "Immune cell","Lymphoid cell","T cell","CD8+ T cells"),
  create_marker_rows(c("PDCD1","HAVCR2","LAG3","TIGIT","TOX","CXCL13"),
    "Immune cell","Lymphoid cell","T cell","CD8+ exhausted T cells"),
  create_marker_rows(c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF18"),
    "Immune cell","Lymphoid cell","T cell","Regulatory T cells"),
  create_marker_rows(c("NCAM1","NKG7","GNLY","KLRD1","NCR1"),
    "Immune cell","Lymphoid cell","Natural killer cell","NK cells"),
  create_marker_rows(c("MS4A1","CD19","CD79A","CD79B","BANK1"),
    "Immune cell","Lymphoid cell","B cell","B cells"),
  create_marker_rows(c("MZB1","JCHAIN","XBP1","SDC1","DERL3"),
    "Immune cell","Lymphoid cell","Plasma cell","Plasma cells"),
  create_marker_rows(c("CD14","FCN1","VCAN","S100A8","S100A9"),
    "Immune cell","Myeloid cell","Monocyte","Monocytes"),
  create_marker_rows(c("CD68","CD163","C1QA","C1QB","APOE"),
    "Immune cell","Myeloid cell","Macrophage","Macrophages"),
  create_marker_rows(c("C1QC","APOC1","TREM2","GPNMB","FOLR2"),
    "Immune cell","Myeloid cell","Tumour-associated macrophage","TAM C1QC (immunoregulatory)"),
  create_marker_rows(c("SPP1","MARCO","FN1","MMP9","INHBA"),
    "Immune cell","Myeloid cell","Tumour-associated macrophage","TAM SPP1 (angiogenic)"),
  create_marker_rows(c("CLEC9A","XCR1","BATF3","IRF8","CADM1"),
    "Immune cell","Myeloid cell","Dendritic cell","Conventional dendritic cells type 1 (cDC1)"),
  create_marker_rows(c("CD1C","FCER1A","CLEC10A","CD1E"),
    "Immune cell","Myeloid cell","Dendritic cell","Conventional dendritic cells type 2 (cDC2)"),
  create_marker_rows(c("LILRA4","CLEC4C","IL3RA","GZMB","IRF7"),
    "Immune cell","Myeloid cell","Dendritic cell","Plasmacytoid dendritic cells"),
  create_marker_rows(c("FCGR3B","CSF3R","S100A8","CXCR2","G0S2"),
    "Immune cell","Myeloid cell","Granulocyte","Neutrophils (TAN)"),
  create_marker_rows(c("TPSAB1","TPSB2","CPA3","KIT","MS4A2"),
    "Immune cell","Myeloid cell","Mast cell","Mast cells"),
  create_marker_rows(c("COL1A1","COL1A2","DCN","LUM","PDGFRA"),
    "Stromal cell","Fibroblast","Fibroblast","Fibroblasts"),
  create_marker_rows(c("ACTA2","TAGLN","POSTN","FAP","CTHRC1"),
    "Stromal cell","Fibroblast","Cancer-associated fibroblast","Myofibroblastic CAF (myCAF)"),
  create_marker_rows(c("IL6","CXCL12","CXCL14","CFD","DPT"),
    "Stromal cell","Fibroblast","Cancer-associated fibroblast","Inflammatory CAF (iCAF)"),
  create_marker_rows(c("RGS5","PDGFRB","NOTCH3","CSPG4","KCNJ8"),
    "Stromal cell","Perivascular mural cell","Pericyte","Pericytes"),
  create_marker_rows(c("PECAM1","VWF","CDH5","CLDN5","CD34"),
    "Endothelial cell","Vascular endothelial cell","Endothelial cell","Endothelial cells"),
  create_marker_rows(c("PROX1","LYVE1","PDPN","FLT4","CCL21"),
    "Endothelial cell","Lymphatic endothelial cell","Lymphatic endothelial cell","Lymphatic endothelial cells")
)

cat("Expert pan-cancer DB (human): rows =", nrow(DBMarker_Expert),
    "| cell types =", length(unique(DBMarker_Expert$Level4)), "\n")
