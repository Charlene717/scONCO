###############################################################################
#  Pan-Cancer Marker Reference  ->  scMRMA format
#  (gene | Level1 | Level2 | Level3 | Level4_Abb)
#
#  scONCO v1.0 — pan-cancer tumour-microenvironment (TME) marker database.
#
#  Design uses the same 4-list structure:
#    - marker_sets_ChatGPT : one list (cell type -> gene vector)
#    - level_info              : one list (cell type -> 4-level hierarchy)
#    - level4_abbrev_map       : full name -> abbreviation
#    - build_scMRMA()          : assemble long-format markers_df
#
#  Scope (modular):
#    * PAN-CANCER CORE  — malignant (generic), immune, stromal, endothelial,
#                         normal-epithelial compartments shared across cancers.
#    * CANCER-TYPE MODULES — tissue-of-origin malignant signatures
#      (melanoma, LUAD, LUSC, BRCA, HCC, CRC, PDAC, STAD, OV, PRAD, GBM, RCC,
#       HNSC, BLCA). Tag = Level2 "Malignant (<type>)". A user can subset by
#       module via load_cancer_marker_db(module = "melanoma") (see loader).
#
#  Primary references for marker choices (see literature/ + endnote_library):
#    - scATOMIC pan-cancer reference (Nofech-Mozes 2023, Nat Commun)
#    - scCancer2 TME annotation (Chen 2024, Bioinformatics)
#    - Pan-cancer T-cell atlas (Zheng L. 2021, Science)
#    - Pan-cancer myeloid atlas (Cheng S. 2021, Cell)
#    - Pan-cancer CAF atlas (Galbo 2021; Lavie 2022, Nat Cancer)
#    - Tumour endothelial atlas (Goveia 2020, Cancer Cell)
###############################################################################

if (!require("dplyr")) { install.packages("dplyr"); library(dplyr) }

###############################################################################
## 0) Helpers
###############################################################################
uniq_keep_order <- function(x) x[!duplicated(x)]

###############################################################################
## 1) marker_sets_ChatGPT (defined completely here)
###############################################################################
marker_sets_ChatGPT <- list(

  # =========================================================================
  # MALIGNANT — generic / programme-level (pan-cancer)
  # =========================================================================
  "Malignant cell (generic)" = uniq_keep_order(c(
    "EPCAM","KRT8","KRT18","KRT19","KRT7","CDH1","MUC1","KRT17",
    "S100A4","SOX4","CLDN3","CLDN4","CD24","ELF3","KRT5","SPINT2"
  )),
  "Cycling malignant cell" = uniq_keep_order(c(
    "MKI67","TOP2A","CDK1","UBE2C","CCNB1","CCNB2","CENPF","BIRC5",
    "TYMS","PCNA","RRM2","NUSAP1","ASPM","CDC20","STMN1","HMGB2"
  )),
  "EMT-like malignant cell" = uniq_keep_order(c(
    "VIM","ZEB1","ZEB2","SNAI1","SNAI2","TWIST1","FN1","CDH2",
    "SPARC","TGFBI","LOXL2","COL1A1","SERPINE1","MMP2","ITGB1","S100A4"
  )),

  # =========================================================================
  # MALIGNANT — cancer-type modules (tissue-of-origin signatures)
  # =========================================================================
  "Malignant melanoma cell" = uniq_keep_order(c(
    "MLANA","PMEL","TYR","TYRP1","DCT","MITF","SOX10","S100B",
    "MIA","GPNMB","SLC45A2","PRAME","AXL","NGFR","EDNRB"
  )),
  "Malignant lung adenocarcinoma cell" = uniq_keep_order(c(
    "NKX2-1","NAPSA","SFTPC","SFTPB","SFTPA1","SFTPA2","MUC1","SLPI",
    "SCGB3A1","KRT7","EPCAM","CEACAM5","SFTA2","LAMP3"
  )),
  "Malignant lung squamous cell" = uniq_keep_order(c(
    "TP63","KRT5","KRT6A","KRT14","KRT17","SOX2","DSG3","PKP1",
    "S100A2","KRT15","SERPINB5","NTRK2","PERP"
  )),
  "Malignant breast cell" = uniq_keep_order(c(
    "EPCAM","KRT8","KRT18","KRT19","ESR1","ERBB2","GATA3","FOXA1",
    "KRT7","ANKRD30A","TFF1","XBP1","AGR2","MUCL1","KRT5","KRT14"
  )),
  "Malignant hepatocellular cell" = uniq_keep_order(c(
    "ALB","APOA1","APOA2","APOC3","TTR","SERPINA1","FGB","FGA",
    "APOB","AHSG","TF","GPC3","AFP","HP","CYP2E1"
  )),
  "Malignant colorectal cell" = uniq_keep_order(c(
    "CDX2","KRT20","CEACAM5","CEACAM6","EPCAM","TFF3","FABP1","KRT8",
    "MUC2","LGR5","OLFM4","REG4","PHGR1","GUCA2B"
  )),
  "Malignant pancreatic cell" = uniq_keep_order(c(
    "KRT19","KRT8","KRT18","MUC1","CEACAM5","CEACAM6","TFF1","TFF2",
    "S100P","LGALS4","SPINK1","AGR2","CLDN18","KRT7","SOX9"
  )),
  "Malignant gastric cell" = uniq_keep_order(c(
    "MUC5AC","TFF1","TFF2","GKN1","GKN2","PGC","MUC6","LIPF",
    "CLDN18","KRT8","KRT18","EPCAM","CEACAM5","AGR2"
  )),
  "Malignant ovarian cell" = uniq_keep_order(c(
    "PAX8","WT1","MUC16","KRT8","KRT18","KRT7","EPCAM","CD24",
    "MSLN","FOLR1","SOX17","CLDN3","WFDC2","STAR"
  )),
  "Malignant prostate cell" = uniq_keep_order(c(
    "KLK3","KLK2","AR","NKX3-1","ACPP","FOLH1","STEAP2","TMPRSS2",
    "AMACR","HOXB13","MSMB","KRT8","KRT18","EPCAM"
  )),
  "Malignant glioma cell" = uniq_keep_order(c(
    "GFAP","SOX2","OLIG1","OLIG2","S100B","VIM","EGFR","PTPRZ1",
    "BCAN","SOX9","NES","CD44","CHI3L1","APOE","MBP"
  )),
  "Malignant renal cell" = uniq_keep_order(c(
    "CA9","NDUFA4L2","VEGFA","ANGPTL4","SLC17A3","KRT18","KRT8","PAX8",
    "VIM","NNMT","SPP1","ALDOB","UMOD","CP"
  )),
  "Malignant head and neck squamous cell" = uniq_keep_order(c(
    "TP63","KRT5","KRT6A","KRT14","KRT17","SFN","S100A2","PERP",
    "KRT15","DSG3","PI3","SPRR1B","CDH1","EPCAM"
  )),
  "Malignant bladder urothelial cell" = uniq_keep_order(c(
    "UPK1B","UPK2","UPK3A","KRT20","KRT13","KRT5","GATA3","FOXA1",
    "PPARG","S100P","KRT8","KRT18","EPCAM","CD24"
  )),

  # =========================================================================
  # EPITHELIAL — normal / tissue-resident (context compartments)
  # =========================================================================
  "Epithelial cells" = uniq_keep_order(c(
    "EPCAM","CDH1","KRT8","KRT18","KRT19","KRT7","CLDN3","CLDN4",
    "MUC1","CD24","ELF3","SPINT2","KRT5","CDH3"
  )),
  "Alveolar epithelial cells" = uniq_keep_order(c(
    "SFTPC","SFTPB","SFTPA1","SFTPA2","AGER","PDPN","NAPSA","SCGB1A1",
    "SCGB3A1","AQP5","HOPX","NKX2-1"
  )),
  "Hepatocytes" = uniq_keep_order(c(
    "ALB","TTR","APOA1","SERPINA1","CYP2E1","HP","TF","FGB",
    "APOC3","ASGR1","HNF4A","CPS1"
  )),
  "Secretory/Goblet epithelial cells" = uniq_keep_order(c(
    "MUC2","TFF3","SPINK4","FCGBP","ZG16","AGR2","REG4","CLCA1",
    "MUC5AC","TFF1","LGALS4","KRT8"
  )),

  # =========================================================================
  # IMMUNE — T / NK lymphoid
  # =========================================================================
  "T cells" = uniq_keep_order(c(
    "CD3D","CD3E","CD3G","TRAC","TRBC1","TRBC2","CD2","CD5",
    "IL7R","CD4","CD8A","CD8B","CCR7","LCK"
  )),
  "CD4+ T cells" = uniq_keep_order(c(
    "CD4","IL7R","CD40LG","CCR7","TCF7","LEF1","SELL","FOXP3",
    "RORA","CXCR4","ANXA1","LTB"
  )),
  "CD8+ T cells" = uniq_keep_order(c(
    "CD8A","CD8B","GZMK","GZMA","NKG7","CCL5","GZMH","CD3D",
    "EOMES","CST7","CXCR3","KLRG1"
  )),
  "CD8+ cytotoxic T cells" = uniq_keep_order(c(
    "CD8A","GZMB","PRF1","GNLY","NKG7","FGFBP2","KLRD1","GZMH",
    "FCGR3A","CX3CR1","FASLG","KLRG1"
  )),
  "CD8+ exhausted T cells" = uniq_keep_order(c(
    "PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","TOX","ENTPD1","LAYN",
    "CXCL13","GZMB","ITGAE","TNFRSF9"
  )),
  "Naive/memory T cells" = uniq_keep_order(c(
    "CCR7","SELL","TCF7","LEF1","IL7R","CD27","CD28","S1PR1",
    "MAL","NOSIP","TXK","ACTN1"
  )),
  "Regulatory T cells" = uniq_keep_order(c(
    "FOXP3","IL2RA","CTLA4","TNFRSF18","IKZF2","TIGIT","TNFRSF4","BATF",
    "CCR8","LAYN","ENTPD1","CD4"
  )),
  "T follicular helper cells" = uniq_keep_order(c(
    "CXCL13","BCL6","PDCD1","ICOS","CD200","TOX2","IL21","CXCR5",
    "BTLA","CD4","MAF","TOX"
  )),
  "MAIT cells" = uniq_keep_order(c(
    "SLC4A10","KLRB1","RORC","RORA","IL7R","CCR6","ZBTB16","NCR3",
    "DPP4","CXCR6","GZMK","TRAV1-2"
  )),
  "γδ T cells" = uniq_keep_order(c(
    "TRDC","TRGC1","TRGC2","TRDV2","TRGV9","KLRD1","NKG7","GZMB",
    "CD3D","KLRC1","CD7","FCGR3A"
  )),
  "Proliferating T cells" = uniq_keep_order(c(
    "MKI67","TOP2A","STMN1","TUBB","HMGB2","TYMS","CD3D","UBE2C",
    "BIRC5","CKS1B","RRM2","CENPF"
  )),
  "NK cells" = uniq_keep_order(c(
    "NCAM1","NKG7","GNLY","KLRD1","KLRF1","NCR1","FCGR3A","GZMB",
    "PRF1","KLRC1","XCL1","FGFBP2","TYROBP","KLRB1"
  )),

  # =========================================================================
  # IMMUNE — B / Plasma
  # =========================================================================
  "B cells" = uniq_keep_order(c(
    "MS4A1","CD19","CD79A","CD79B","CD20","BANK1","HLA-DRA","TCL1A",
    "IGHD","IGHM","CD22","VPREB3","PAX5"
  )),
  "Germinal center B cells" = uniq_keep_order(c(
    "BCL6","AICDA","RGS13","MEF2B","LRMP","SUGCT","STMN1","MKI67",
    "CD83","NEIL1","MS4A1","CD79A"
  )),
  "Plasma cells" = uniq_keep_order(c(
    "MZB1","JCHAIN","XBP1","SDC1","PRDM1","DERL3","TNFRSF17","IGHG1",
    "IGHA1","IGKC","CD38","SLAMF7","FKBP11"
  )),

  # =========================================================================
  # IMMUNE — Myeloid (mono / macro / TAM / DC / neutrophil / mast)
  # =========================================================================
  "Monocytes" = uniq_keep_order(c(
    "CD14","FCN1","VCAN","S100A8","S100A9","S100A12","LYZ","FCGR3A",
    "CSF1R","CD300E","SERPINA1","CLEC12A","EREG"
  )),
  "Macrophages" = uniq_keep_order(c(
    "CD68","CD163","C1QA","C1QB","C1QC","MRC1","APOE","FCGR3A",
    "LYZ","CSF1R","MARCO","SPP1","TREM2"
  )),
  "TAM C1QC (immunoregulatory)" = uniq_keep_order(c(
    "C1QA","C1QB","C1QC","APOE","APOC1","TREM2","GPNMB","SLC40A1",
    "FOLR2","SELENOP","CD163","MRC1"
  )),
  "TAM SPP1 (angiogenic)" = uniq_keep_order(c(
    "SPP1","MARCO","FN1","MMP9","MMP12","CCL2","VCAN","INHBA",
    "SLC2A1","CTSB","MIF","CD68"
  )),
  "TAM FOLR2 (tissue-resident-like)" = uniq_keep_order(c(
    "FOLR2","LYVE1","SELENOP","SLC40A1","F13A1","MRC1","CD163L1","STAB1",
    "PLTP","C1QA","CD163","RNASE1"
  )),
  "M1-like macrophages" = uniq_keep_order(c(
    "CXCL9","CXCL10","CXCL11","GBP1","STAT1","IDO1","IL1B","TNF",
    "CCL5","ISG15","CD68","FCGR1A"
  )),
  "M2-like macrophages" = uniq_keep_order(c(
    "MRC1","CD163","MSR1","MERTK","STAB1","MAF","CCL18","SLC40A1",
    "F13A1","LGMN","CD68","TGFB1"
  )),
  "Conventional dendritic cells type 1 (cDC1)" = uniq_keep_order(c(
    "CLEC9A","XCR1","BATF3","IRF8","CADM1","C1ORF54","SNX22","WDFY4",
    "THBD","IDO1","HLA-DRA"
  )),
  "Conventional dendritic cells type 2 (cDC2)" = uniq_keep_order(c(
    "CD1C","FCER1A","CLEC10A","CD1E","FCGR2B","HLA-DQA1","CD1A","ITGAX",
    "SIRPA","CLEC4A","HLA-DRA"
  )),
  "Mature regulatory DC (mregDC/LAMP3)" = uniq_keep_order(c(
    "LAMP3","CCR7","FSCN1","CCL19","CCL22","BIRC3","IL7R","CD40",
    "RELB","IDO1","MARCKSL1","CD274"
  )),
  "Plasmacytoid dendritic cells" = uniq_keep_order(c(
    "LILRA4","CLEC4C","IL3RA","GZMB","JCHAIN","IRF7","TCF4","SPIB",
    "PLD4","SERPINF1","ITM2C","BCL11A"
  )),
  "Neutrophils (TAN)" = uniq_keep_order(c(
    "FCGR3B","CSF3R","S100A8","S100A9","CXCR2","FFAR2","G0S2","CMTM2",
    "PROK2","BASP1","SOD2","IL1B"
  )),
  "Mast cells" = uniq_keep_order(c(
    "TPSAB1","TPSB2","CPA3","KIT","MS4A2","HPGDS","GATA2","CMA1",
    "IL1RL1","SLC18A2","HDC","RGS13"
  )),

  # =========================================================================
  # STROMAL — Fibroblast / CAF / mural
  # =========================================================================
  "Fibroblasts" = uniq_keep_order(c(
    "COL1A1","COL1A2","COL3A1","DCN","LUM","PDGFRA","FBLN1","MGP",
    "COL6A1","COL6A2","SFRP2","CFD","GSN"
  )),
  "Myofibroblastic CAF (myCAF)" = uniq_keep_order(c(
    "ACTA2","TAGLN","MYH11","MYL9","POSTN","COL10A1","COL11A1","FAP",
    "INHBA","CTHRC1","TPM2","MMP11","LRRC15"
  )),
  "Inflammatory CAF (iCAF)" = uniq_keep_order(c(
    "IL6","CXCL12","CXCL14","CFD","PDGFRA","DPT","CXCL1","CXCL2",
    "CCL2","C3","C7","HAS1","PLA2G2A"
  )),
  "Antigen-presenting CAF (apCAF)" = uniq_keep_order(c(
    "CD74","HLA-DRA","HLA-DRB1","HLA-DPA1","HLA-DPB1","SLPI","SAA1","CD9",
    "PIGR","CLU","HLA-DQA1","KRT8"
  )),
  "Pericytes" = uniq_keep_order(c(
    "RGS5","PDGFRB","NOTCH3","CSPG4","ACTA2","MCAM","KCNJ8","ABCC9",
    "HIGD1B","CALD1","NDUFA4L2","COX4I2"
  )),
  "Smooth muscle cells" = uniq_keep_order(c(
    "ACTA2","MYH11","TAGLN","CNN1","DES","MYL9","LMOD1","MYLK",
    "PLN","ACTG2","SMTN","CASQ2"
  )),

  # =========================================================================
  # ENDOTHELIAL
  # =========================================================================
  "Endothelial cells" = uniq_keep_order(c(
    "PECAM1","VWF","CDH5","CLDN5","CD34","ENG","FLT1","KDR",
    "EGFL7","RAMP2","CLEC14A","ECSCR","AQP1"
  )),
  "Tip / angiogenic endothelial cells" = uniq_keep_order(c(
    "ESM1","NID2","PGF","APLN","ANGPT2","INSR","COL4A1","COL4A2",
    "PXDN","NOTCH4","KDR","RGCC"
  )),
  "Venous endothelial cells" = uniq_keep_order(c(
    "ACKR1","SELP","SELE","VWF","VCAM1","CCL14","CPE","IL1R1",
    "NR2F2","PLVAP","SOCS3"
  )),
  "Lymphatic endothelial cells" = uniq_keep_order(c(
    "PROX1","LYVE1","PDPN","FLT4","CCL21","MMRN1","TFF3","RELN",
    "NTS","CLDN5","PECAM1","NRP2"
  )),

  # =========================================================================
  # OTHER / housekeeping
  # =========================================================================
  "Cycling cells (lineage-agnostic)" = uniq_keep_order(c(
    "MKI67","TOP2A","CDK1","UBE2C","BIRC5","CCNB1","CENPF","NUSAP1",
    "PCNA","STMN1","HMGB2","TYMS"
  )),
  "Erythrocytes" = uniq_keep_order(c(
    "HBB","HBA1","HBA2","ALAS2","SLC4A1","GYPA","AHSP","CA1",
    "HBM","SNCA","BNIP3L"
  ))
)

###############################################################################
## 2) level_info  (Level1, Level2, Level3_full, Level4_Abbrev_source)
###############################################################################
level_info <- list(
  # ---- Malignant generic ----
  "Malignant cell (generic)" = c("Malignant cell","Malignant (pan-cancer)","Malignant epithelial / tumour cell","Malignant cell (generic)"),
  "Cycling malignant cell"   = c("Malignant cell","Malignant (pan-cancer)","Cycling malignant cell","Cycling malignant cell"),
  "EMT-like malignant cell"  = c("Malignant cell","Malignant (pan-cancer)","EMT-like malignant cell","EMT-like malignant cell"),

  # ---- Malignant cancer-type modules ----
  "Malignant melanoma cell"                 = c("Malignant cell","Malignant (melanoma)","Melanoma cell","Malignant melanoma cell"),
  "Malignant lung adenocarcinoma cell"      = c("Malignant cell","Malignant (LUAD)","Lung adenocarcinoma cell","Malignant lung adenocarcinoma cell"),
  "Malignant lung squamous cell"            = c("Malignant cell","Malignant (LUSC)","Lung squamous carcinoma cell","Malignant lung squamous cell"),
  "Malignant breast cell"                   = c("Malignant cell","Malignant (BRCA)","Breast carcinoma cell","Malignant breast cell"),
  "Malignant hepatocellular cell"           = c("Malignant cell","Malignant (HCC)","Hepatocellular carcinoma cell","Malignant hepatocellular cell"),
  "Malignant colorectal cell"               = c("Malignant cell","Malignant (CRC)","Colorectal carcinoma cell","Malignant colorectal cell"),
  "Malignant pancreatic cell"               = c("Malignant cell","Malignant (PDAC)","Pancreatic ductal carcinoma cell","Malignant pancreatic cell"),
  "Malignant gastric cell"                  = c("Malignant cell","Malignant (STAD)","Gastric carcinoma cell","Malignant gastric cell"),
  "Malignant ovarian cell"                  = c("Malignant cell","Malignant (OV)","Ovarian carcinoma cell","Malignant ovarian cell"),
  "Malignant prostate cell"                 = c("Malignant cell","Malignant (PRAD)","Prostate carcinoma cell","Malignant prostate cell"),
  "Malignant glioma cell"                   = c("Malignant cell","Malignant (GBM/glioma)","Glioma / glioblastoma cell","Malignant glioma cell"),
  "Malignant renal cell"                    = c("Malignant cell","Malignant (RCC)","Renal cell carcinoma cell","Malignant renal cell"),
  "Malignant head and neck squamous cell"   = c("Malignant cell","Malignant (HNSC)","Head & neck squamous carcinoma cell","Malignant head and neck squamous cell"),
  "Malignant bladder urothelial cell"       = c("Malignant cell","Malignant (BLCA)","Bladder urothelial carcinoma cell","Malignant bladder urothelial cell"),

  # ---- Normal / tissue-resident epithelial ----
  "Epithelial cells"                  = c("Epithelial cell","Normal epithelial cell","Epithelial cell","Epithelial cells"),
  "Alveolar epithelial cells"         = c("Epithelial cell","Normal epithelial cell","Alveolar epithelial cell","Alveolar epithelial cells"),
  "Hepatocytes"                       = c("Epithelial cell","Normal epithelial cell","Hepatocyte","Hepatocytes"),
  "Secretory/Goblet epithelial cells" = c("Epithelial cell","Normal epithelial cell","Secretory/Goblet cell","Secretory/Goblet epithelial cells"),

  # ---- Immune: T / NK ----
  "T cells"                  = c("Immune cell","Lymphoid cell","T cell","T cells"),
  "CD4+ T cells"             = c("Immune cell","Lymphoid cell","T cell","CD4+ T cells"),
  "CD8+ T cells"             = c("Immune cell","Lymphoid cell","T cell","CD8+ T cells"),
  "CD8+ cytotoxic T cells"   = c("Immune cell","Lymphoid cell","T cell","CD8+ cytotoxic T cells"),
  "CD8+ exhausted T cells"   = c("Immune cell","Lymphoid cell","T cell","CD8+ exhausted T cells"),
  "Naive/memory T cells"     = c("Immune cell","Lymphoid cell","T cell","Naive/memory T cells"),
  "Regulatory T cells"       = c("Immune cell","Lymphoid cell","T cell","Regulatory T cells"),
  "T follicular helper cells"= c("Immune cell","Lymphoid cell","T cell","T follicular helper cells"),
  "MAIT cells"               = c("Immune cell","Lymphoid cell","T cell","MAIT cells"),
  "γδ T cells"               = c("Immune cell","Lymphoid cell","T cell","γδ T cells"),
  "Proliferating T cells"    = c("Immune cell","Lymphoid cell","T cell","Proliferating T cells"),
  "NK cells"                 = c("Immune cell","Lymphoid cell","Natural killer cell","NK cells"),

  # ---- Immune: B / Plasma ----
  "B cells"                  = c("Immune cell","Lymphoid cell","B cell","B cells"),
  "Germinal center B cells"  = c("Immune cell","Lymphoid cell","B cell","Germinal center B cells"),
  "Plasma cells"             = c("Immune cell","Lymphoid cell","Plasma cell","Plasma cells"),

  # ---- Immune: Myeloid ----
  "Monocytes"                              = c("Immune cell","Myeloid cell","Monocyte","Monocytes"),
  "Macrophages"                            = c("Immune cell","Myeloid cell","Macrophage","Macrophages"),
  "TAM C1QC (immunoregulatory)"            = c("Immune cell","Myeloid cell","Tumour-associated macrophage","TAM C1QC (immunoregulatory)"),
  "TAM SPP1 (angiogenic)"                  = c("Immune cell","Myeloid cell","Tumour-associated macrophage","TAM SPP1 (angiogenic)"),
  "TAM FOLR2 (tissue-resident-like)"       = c("Immune cell","Myeloid cell","Tumour-associated macrophage","TAM FOLR2 (tissue-resident-like)"),
  "M1-like macrophages"                    = c("Immune cell","Myeloid cell","Macrophage","M1-like macrophages"),
  "M2-like macrophages"                    = c("Immune cell","Myeloid cell","Macrophage","M2-like macrophages"),
  "Conventional dendritic cells type 1 (cDC1)" = c("Immune cell","Myeloid cell","Dendritic cell","Conventional dendritic cells type 1 (cDC1)"),
  "Conventional dendritic cells type 2 (cDC2)" = c("Immune cell","Myeloid cell","Dendritic cell","Conventional dendritic cells type 2 (cDC2)"),
  "Mature regulatory DC (mregDC/LAMP3)"    = c("Immune cell","Myeloid cell","Dendritic cell","Mature regulatory DC (mregDC/LAMP3)"),
  "Plasmacytoid dendritic cells"           = c("Immune cell","Myeloid cell","Dendritic cell","Plasmacytoid dendritic cells"),
  "Neutrophils (TAN)"                      = c("Immune cell","Myeloid cell","Granulocyte","Neutrophils (TAN)"),
  "Mast cells"                             = c("Immune cell","Myeloid cell","Mast cell","Mast cells"),

  # ---- Stromal: Fibroblast / CAF / mural ----
  "Fibroblasts"                            = c("Stromal cell","Fibroblast","Fibroblast","Fibroblasts"),
  "Myofibroblastic CAF (myCAF)"            = c("Stromal cell","Fibroblast","Cancer-associated fibroblast","Myofibroblastic CAF (myCAF)"),
  "Inflammatory CAF (iCAF)"                = c("Stromal cell","Fibroblast","Cancer-associated fibroblast","Inflammatory CAF (iCAF)"),
  "Antigen-presenting CAF (apCAF)"         = c("Stromal cell","Fibroblast","Cancer-associated fibroblast","Antigen-presenting CAF (apCAF)"),
  "Pericytes"                              = c("Stromal cell","Perivascular mural cell","Pericyte","Pericytes"),
  "Smooth muscle cells"                    = c("Stromal cell","Perivascular mural cell","Smooth muscle cell","Smooth muscle cells"),

  # ---- Endothelial ----
  "Endothelial cells"                      = c("Endothelial cell","Vascular endothelial cell","Endothelial cell","Endothelial cells"),
  "Tip / angiogenic endothelial cells"     = c("Endothelial cell","Vascular endothelial cell","Tip endothelial cell","Tip / angiogenic endothelial cells"),
  "Venous endothelial cells"               = c("Endothelial cell","Vascular endothelial cell","Venous endothelial cell","Venous endothelial cells"),
  "Lymphatic endothelial cells"            = c("Endothelial cell","Lymphatic endothelial cell","Lymphatic endothelial cell","Lymphatic endothelial cells"),

  # ---- Other ----
  "Cycling cells (lineage-agnostic)"       = c("Other","Proliferating cell","Cycling cell","Cycling cells (lineage-agnostic)"),
  "Erythrocytes"                           = c("Other","Erythroid cell","Erythrocyte","Erythrocytes")
)

###############################################################################
## 3) Build scMRMA long table  (outputs Level4_Abb directly)
###############################################################################
build_scMRMA <- function(marker_sets, level_info) {
  missing <- setdiff(names(marker_sets), names(level_info))
  if (length(missing)) {
    stop("The following cell types exist in marker_sets but are missing in level_info (names must match):\n",
         paste(missing, collapse = ", "))
  }

  do.call(rbind, lapply(names(marker_sets), function(ct) {
    genes <- uniq_keep_order(marker_sets[[ct]])
    lv <- level_info[[ct]]
    data.frame(
      gene       = genes,
      Level1     = lv[1],
      Level2     = lv[2],
      Level3     = lv[3],
      Level4_Abb = lv[4],
      stringsAsFactors = FALSE
    )
  }))
}

marker_df_scMRMA_ChatGPT <- build_scMRMA(marker_sets_ChatGPT, level_info)

###############################################################################
## 4) Level4_Abb abbreviation mapping (full name -> abbreviation)
###############################################################################
level4_abbrev_map <- c(
  # ---- Malignant generic ----
  "Malignant cell (generic)" = "Mal",
  "Cycling malignant cell"   = "Mal_Cyc",
  "EMT-like malignant cell"  = "Mal_EMT",

  # ---- Malignant cancer-type modules ----
  "Malignant melanoma cell"               = "Mal_Mel",
  "Malignant lung adenocarcinoma cell"    = "Mal_LUAD",
  "Malignant lung squamous cell"          = "Mal_LUSC",
  "Malignant breast cell"                 = "Mal_BRCA",
  "Malignant hepatocellular cell"         = "Mal_HCC",
  "Malignant colorectal cell"             = "Mal_CRC",
  "Malignant pancreatic cell"             = "Mal_PDAC",
  "Malignant gastric cell"                = "Mal_STAD",
  "Malignant ovarian cell"                = "Mal_OV",
  "Malignant prostate cell"               = "Mal_PRAD",
  "Malignant glioma cell"                 = "Mal_GBM",
  "Malignant renal cell"                  = "Mal_RCC",
  "Malignant head and neck squamous cell" = "Mal_HNSC",
  "Malignant bladder urothelial cell"     = "Mal_BLCA",

  # ---- Normal epithelial ----
  "Epithelial cells"                  = "Epi",
  "Alveolar epithelial cells"         = "Epi_AT",
  "Hepatocytes"                       = "Hep",
  "Secretory/Goblet epithelial cells" = "Epi_Gob",

  # ---- Immune: T / NK ----
  "T cells"                   = "Tcell",
  "CD4+ T cells"              = "CD4T",
  "CD8+ T cells"              = "CD8T",
  "CD8+ cytotoxic T cells"    = "CD8_Tcyt",
  "CD8+ exhausted T cells"    = "CD8_Tex",
  "Naive/memory T cells"      = "T_naive_mem",
  "Regulatory T cells"        = "Treg",
  "T follicular helper cells" = "Tfh",
  "MAIT cells"                = "MAIT",
  "γδ T cells"                = "gdT",
  "Proliferating T cells"     = "T_prolif",
  "NK cells"                  = "NK",

  # ---- Immune: B / Plasma ----
  "B cells"                  = "Bcell",
  "Germinal center B cells"  = "B_GC",
  "Plasma cells"             = "Plasma",

  # ---- Immune: Myeloid ----
  "Monocytes"                              = "Mono",
  "Macrophages"                            = "Mac",
  "TAM C1QC (immunoregulatory)"            = "TAM_C1QC",
  "TAM SPP1 (angiogenic)"                  = "TAM_SPP1",
  "TAM FOLR2 (tissue-resident-like)"       = "TAM_FOLR2",
  "M1-like macrophages"                    = "Mac_M1",
  "M2-like macrophages"                    = "Mac_M2",
  "Conventional dendritic cells type 1 (cDC1)" = "cDC1",
  "Conventional dendritic cells type 2 (cDC2)" = "cDC2",
  "Mature regulatory DC (mregDC/LAMP3)"    = "mregDC",
  "Plasmacytoid dendritic cells"           = "pDC",
  "Neutrophils (TAN)"                      = "Neu",
  "Mast cells"                             = "Mast",

  # ---- Stromal ----
  "Fibroblasts"                            = "FB",
  "Myofibroblastic CAF (myCAF)"            = "myCAF",
  "Inflammatory CAF (iCAF)"                = "iCAF",
  "Antigen-presenting CAF (apCAF)"         = "apCAF",
  "Pericytes"                              = "Peri",
  "Smooth muscle cells"                    = "SMC",

  # ---- Endothelial ----
  "Endothelial cells"                      = "EC",
  "Tip / angiogenic endothelial cells"     = "EC_tip",
  "Venous endothelial cells"               = "vEC",
  "Lymphatic endothelial cells"            = "lEC",

  # ---- Other ----
  "Cycling cells (lineage-agnostic)"       = "Cycling",
  "Erythrocytes"                           = "RBC"
)

old_abb <- marker_df_scMRMA_ChatGPT$Level4_Abb
marker_df_scMRMA_ChatGPT$Level4_Abb <- ifelse(
  old_abb %in% names(level4_abbrev_map),
  unname(level4_abbrev_map[old_abb]),
  old_abb
)

###############################################################################
## 5) Sanity check + Export
###############################################################################
cat("scONCO pan-cancer DB v1.0\n")
cat("Cell types (Level4_Abb): ", length(unique(marker_df_scMRMA_ChatGPT$Level4_Abb)), "\n", sep = "")
cat("Total marker rows       : ", nrow(marker_df_scMRMA_ChatGPT), "\n", sep = "")
cat("Unique genes            : ", length(unique(marker_df_scMRMA_ChatGPT$gene)), "\n", sep = "")
cat("Level1 lineages         : ", paste(unique(marker_df_scMRMA_ChatGPT$Level1), collapse = ", "), "\n\n", sep = "")

markers_df <- marker_df_scMRMA_ChatGPT
