#### Expert-curated Marker DB — Pan-cancer (Mouse) ####
## Dr. 戴揚紘 (Yang-Hong Dai) — radiation oncology / bioinformatics curator
## Expert-curated mouse pan-cancer / syngeneic-tumour TME marker list.
## Companion to DBMarker_PanCancer_Expert.R (human). MGI symbols.

library(dplyr)
library(tibble)

create_marker_rows <- function(genes, level1, level2, level3, level4) {
  tibble(gene = genes, Level1 = level1, Level2 = level2, Level3 = level3, Level4 = level4)
}

DBMarker_Expert_Mouse <- bind_rows(
  create_marker_rows(c("Epcam","Krt8","Krt18","Krt19","Cdh1","Muc1"),
    "Malignant cell","Malignant (pan-cancer)","Malignant epithelial cell","Malignant cell (generic)"),
  create_marker_rows(c("Mki67","Top2a","Ube2c","Birc5","Cenpf"),
    "Malignant cell","Malignant (pan-cancer)","Cycling malignant cell","Cycling malignant cell"),
  create_marker_rows(c("Cd3d","Cd3e","Trac","Cd2","Il7r"),
    "Immune cell","Lymphoid cell","T cell","T cells"),
  create_marker_rows(c("Cd8a","Cd8b1","Gzmk","Nkg7","Pdcd1","Havcr2","Lag3","Tox"),
    "Immune cell","Lymphoid cell","T cell","CD8+ exhausted T cells"),
  create_marker_rows(c("Foxp3","Il2ra","Ctla4","Ikzf2","Tnfrsf18"),
    "Immune cell","Lymphoid cell","T cell","Regulatory T cells"),
  create_marker_rows(c("Ncr1","Klrb1c","Klrd1","Gzmb","Nkg7"),
    "Immune cell","Lymphoid cell","Natural killer cell","NK cells"),
  create_marker_rows(c("Ms4a1","Cd19","Cd79a","Cd79b","Bank1"),
    "Immune cell","Lymphoid cell","B cell","B cells"),
  create_marker_rows(c("Mzb1","Jchain","Xbp1","Sdc1","Derl3"),
    "Immune cell","Lymphoid cell","Plasma cell","Plasma cells"),
  create_marker_rows(c("Ly6c2","Fcn1","Vcan","S100a8","Lyz2"),
    "Immune cell","Myeloid cell","Monocyte","Monocytes"),
  create_marker_rows(c("Adgre1","Cd68","Mrc1","C1qa","Apoe"),
    "Immune cell","Myeloid cell","Macrophage","Macrophages"),
  create_marker_rows(c("Spp1","Arg1","Fn1","Mmp12","Vegfa"),
    "Immune cell","Myeloid cell","Tumour-associated macrophage","TAM SPP1 (angiogenic)"),
  create_marker_rows(c("Xcr1","Clec9a","Batf3","Irf8","Cadm1"),
    "Immune cell","Myeloid cell","Dendritic cell","Conventional dendritic cells type 1 (cDC1)"),
  create_marker_rows(c("S100a8","S100a9","Retnlg","Csf3r","Mmp9"),
    "Immune cell","Myeloid cell","Granulocyte","Neutrophils (TAN)"),
  create_marker_rows(c("Cpa3","Mcpt4","Cma1","Kit","Ms4a2"),
    "Immune cell","Myeloid cell","Mast cell","Mast cells"),
  create_marker_rows(c("Col1a1","Col1a2","Dcn","Lum","Pdgfra"),
    "Stromal cell","Fibroblast","Fibroblast","Fibroblasts"),
  create_marker_rows(c("Acta2","Tagln","Postn","Fap","Inhba"),
    "Stromal cell","Fibroblast","Cancer-associated fibroblast","Myofibroblastic CAF (myCAF)"),
  create_marker_rows(c("Rgs5","Pdgfrb","Notch3","Cspg4","Kcnj8"),
    "Stromal cell","Perivascular mural cell","Pericyte","Pericytes"),
  create_marker_rows(c("Pecam1","Cdh5","Cldn5","Egfl7","Flt1"),
    "Endothelial cell","Vascular endothelial cell","Endothelial cell","Endothelial cells"),
  create_marker_rows(c("Prox1","Lyve1","Pdpn","Flt4","Ccl21a"),
    "Endothelial cell","Lymphatic endothelial cell","Lymphatic endothelial cell","Lymphatic endothelial cells")
)

cat("Expert pan-cancer DB (mouse): rows =", nrow(DBMarker_Expert_Mouse),
    "| cell types =", length(unique(DBMarker_Expert_Mouse$Level4)), "\n")
