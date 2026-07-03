## Ref: https://chatgpt.com/c/692ff182-9fd8-8323-8f61-e1e8b798f00a

## ============================================================
##  Pan-cancer / tumour microenvironment (TME) scMRMA Marker Database
##  泛癌 / 腫瘤微環境 (TME) scMRMA 標記基因資料庫 (human gene symbol)
##  - pancancer_markers: Level1~4 + Gene + Source 的 tibble
##  - 每個 cell type 一個 tibble，最後以 bind_rows() 合併
## ============================================================

library(tibble)
library(dplyr)

# ============================================================
# A. Malignant compartment 惡性腫瘤細胞
# ============================================================

# Generic malignant cell 泛癌惡性細胞
malignant_generic <- tibble(
  Level1 = "Malignant cell",
  Level2 = "Malignant (pan-cancer)",
  Level3 = "Malignant cell (generic)",
  Level4 = "Malignant epithelial (generic)",
  Gene   = c("EPCAM","KRT8","KRT18","KRT19","KRT7","CDH1",
             "MUC1","CLDN3","CLDN4","CD24"),
  Source = "Pan-cancer TME atlas"
)

# Cycling malignant cell 增殖中惡性細胞
malignant_cycling <- tibble(
  Level1 = "Malignant cell",
  Level2 = "Malignant (pan-cancer)",
  Level3 = "Cycling malignant cell",
  Level4 = "Proliferating tumour cell",
  Gene   = c("MKI67","TOP2A","CDK1","UBE2C","CCNB1","BIRC5",
             "CENPF","PCNA","STMN1","TYMS"),
  Source = "Pan-cancer TME atlas"
)

# EMT-like malignant cell 上皮間質轉化惡性細胞
malignant_emt <- tibble(
  Level1 = "Malignant cell",
  Level2 = "Malignant (pan-cancer)",
  Level3 = "EMT-like malignant cell",
  Level4 = "Mesenchymal-like tumour cell",
  Gene   = c("VIM","ZEB1","ZEB2","SNAI1","SNAI2","TWIST1",
             "FN1","CDH2","SPARC","MMP2"),
  Source = "Pan-cancer TME atlas"
)

# Malignant melanoma cell 惡性黑色素瘤細胞
malignant_melanoma <- tibble(
  Level1 = "Malignant cell",
  Level2 = "Malignant (melanoma)",
  Level3 = "Malignant melanoma cell",
  Level4 = "Melanoma tumour cell",
  Gene   = c("MLANA","PMEL","TYR","TYRP1","DCT","MITF",
             "SOX10","S100B","GPNMB","PRAME"),
  Source = "Pan-cancer TME atlas"
)

# Malignant lung adenocarcinoma cell 肺腺癌惡性細胞
malignant_luad <- tibble(
  Level1 = "Malignant cell",
  Level2 = "Malignant (LUAD)",
  Level3 = "Malignant lung adenocarcinoma cell",
  Level4 = "LUAD tumour cell",
  Gene   = c("NKX2-1","NAPSA","SFTPC","SFTPB","SFTPA1","MUC1",
             "KRT7","CEACAM5","LAMP3"),
  Source = "Pan-cancer TME atlas"
)

# Malignant lung squamous cell 肺鱗癌惡性細胞
malignant_lusc <- tibble(
  Level1 = "Malignant cell",
  Level2 = "Malignant (LUSC)",
  Level3 = "Malignant lung squamous cell",
  Level4 = "LUSC tumour cell",
  Gene   = c("TP63","KRT5","KRT6A","KRT14","KRT17","SOX2",
             "DSG3","PKP1","S100A2"),
  Source = "Pan-cancer TME atlas"
)

# Malignant breast cell 乳癌惡性細胞
malignant_brca <- tibble(
  Level1 = "Malignant cell",
  Level2 = "Malignant (BRCA)",
  Level3 = "Malignant breast cell",
  Level4 = "Breast tumour cell",
  Gene   = c("EPCAM","KRT8","KRT18","KRT19","ESR1","ERBB2",
             "GATA3","FOXA1","TFF1","XBP1"),
  Source = "Pan-cancer TME atlas"
)

# Malignant hepatocellular cell 肝細胞癌惡性細胞
malignant_hcc <- tibble(
  Level1 = "Malignant cell",
  Level2 = "Malignant (HCC)",
  Level3 = "Malignant hepatocellular cell",
  Level4 = "HCC tumour cell",
  Gene   = c("ALB","APOA1","APOA2","TTR","SERPINA1","FGB",
             "FGA","GPC3","AFP","CYP2E1"),
  Source = "Pan-cancer TME atlas"
)

# Malignant colorectal cell 大腸直腸癌惡性細胞
malignant_crc <- tibble(
  Level1 = "Malignant cell",
  Level2 = "Malignant (CRC)",
  Level3 = "Malignant colorectal cell",
  Level4 = "CRC tumour cell",
  Gene   = c("CDX2","KRT20","CEACAM5","CEACAM6","TFF3","FABP1",
             "MUC2","LGR5","OLFM4","REG4"),
  Source = "Pan-cancer TME atlas"
)

# ============================================================
# B. Normal / non-malignant epithelium 正常上皮
# ============================================================

# Epithelial cells 上皮細胞
epithelial <- tibble(
  Level1 = "Epithelial cell",
  Level2 = "Normal epithelial cell",
  Level3 = "Epithelial cells",
  Level4 = "Epithelial cell",
  Gene   = c("EPCAM","CDH1","KRT8","KRT18","KRT19","KRT7",
             "CLDN3","CLDN4","MUC1","CD24"),
  Source = "Pan-cancer TME atlas"
)

# Alveolar epithelial cells 肺泡上皮細胞
alveolar_epi <- tibble(
  Level1 = "Epithelial cell",
  Level2 = "Normal epithelial cell",
  Level3 = "Alveolar epithelial cells",
  Level4 = "Alveolar epithelial cell",
  Gene   = c("SFTPC","SFTPB","SFTPA1","SFTPA2","AGER","PDPN",
             "NAPSA","SCGB1A1","HOPX","NKX2-1"),
  Source = "Pan-cancer TME atlas"
)

# Hepatocytes 肝細胞
hepatocytes <- tibble(
  Level1 = "Epithelial cell",
  Level2 = "Normal epithelial cell",
  Level3 = "Hepatocytes",
  Level4 = "Hepatocyte",
  Gene   = c("ALB","TTR","APOA1","SERPINA1","CYP2E1","HP",
             "TF","FGB","ASGR1","HNF4A"),
  Source = "Pan-cancer TME atlas"
)

# Secretory / Goblet epithelial cells 分泌型 / 杯狀上皮細胞
goblet_epi <- tibble(
  Level1 = "Epithelial cell",
  Level2 = "Normal epithelial cell",
  Level3 = "Secretory/Goblet epithelial cells",
  Level4 = "Goblet cell",
  Gene   = c("MUC2","TFF3","SPINK4","FCGBP","ZG16","AGR2",
             "REG4","CLCA1"),
  Source = "Pan-cancer TME atlas"
)

# ============================================================
# C. Immune compartment - Lymphoid 免疫細胞 - 淋巴系
# ============================================================

# T cells T 細胞
t_cells <- tibble(
  Level1 = "Immune cell",
  Level2 = "Lymphoid cell",
  Level3 = "T cell",
  Level4 = "T cell",
  Gene   = c("CD3D","CD3E","CD3G","TRAC","CD2","CD5",
             "IL7R","CD4","CD8A","CD8B"),
  Source = "Pan-cancer TME atlas"
)

# CD4+ T cells CD4+ T 細胞
cd4_t <- tibble(
  Level1 = "Immune cell",
  Level2 = "Lymphoid cell",
  Level3 = "T cell",
  Level4 = "CD4+ T cell",
  Gene   = c("CD4","IL7R","CD40LG","CCR7","TCF7","LEF1","SELL"),
  Source = "Pan-cancer TME atlas"
)

# CD8+ T cells CD8+ T 細胞
cd8_t <- tibble(
  Level1 = "Immune cell",
  Level2 = "Lymphoid cell",
  Level3 = "T cell",
  Level4 = "CD8+ T cell",
  Gene   = c("CD8A","CD8B","GZMK","GZMA","NKG7","CCL5",
             "GZMH","EOMES"),
  Source = "Pan-cancer TME atlas"
)

# CD8+ cytotoxic T cells CD8+ 細胞毒性 T 細胞
cd8_cytotoxic <- tibble(
  Level1 = "Immune cell",
  Level2 = "Lymphoid cell",
  Level3 = "CD8+ cytotoxic T cells",
  Level4 = "CD8+ cytotoxic T cell",
  Gene   = c("CD8A","GZMB","PRF1","GNLY","NKG7","FGFBP2",
             "KLRD1","FCGR3A"),
  Source = "Pan-cancer TME atlas"
)

# CD8+ exhausted T cells CD8+ 耗竭 T 細胞
cd8_exhausted <- tibble(
  Level1 = "Immune cell",
  Level2 = "Lymphoid cell",
  Level3 = "CD8+ exhausted T cells",
  Level4 = "CD8+ exhausted T cell",
  Gene   = c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","TOX",
             "ENTPD1","LAYN","CXCL13"),
  Source = "Pan-cancer TME atlas"
)

# Naive/memory T cells 初始 / 記憶 T 細胞
naive_memory_t <- tibble(
  Level1 = "Immune cell",
  Level2 = "Lymphoid cell",
  Level3 = "Naive/memory T cells",
  Level4 = "Naive/memory T cell",
  Gene   = c("CCR7","SELL","TCF7","LEF1","IL7R","CD27",
             "CD28","MAL"),
  Source = "Pan-cancer TME atlas"
)

# Regulatory T cells 調節型 T 細胞
treg <- tibble(
  Level1 = "Immune cell",
  Level2 = "Lymphoid cell",
  Level3 = "Regulatory T cells",
  Level4 = "Regulatory T cell",
  Gene   = c("FOXP3","IL2RA","CTLA4","TNFRSF18","IKZF2","TIGIT",
             "CCR8","LAYN"),
  Source = "Pan-cancer TME atlas"
)

# T follicular helper cells 濾泡輔助型 T 細胞
tfh <- tibble(
  Level1 = "Immune cell",
  Level2 = "Lymphoid cell",
  Level3 = "T follicular helper cells",
  Level4 = "T follicular helper cell",
  Gene   = c("CXCL13","BCL6","PDCD1","ICOS","CD200","IL21",
             "CXCR5","MAF"),
  Source = "Pan-cancer TME atlas"
)

# MAIT cells 黏膜相關恆定 T 細胞
mait <- tibble(
  Level1 = "Immune cell",
  Level2 = "Lymphoid cell",
  Level3 = "MAIT cells",
  Level4 = "MAIT cell",
  Gene   = c("SLC4A10","KLRB1","RORC","RORA","IL7R","CCR6",
             "ZBTB16","TRAV1-2"),
  Source = "Pan-cancer TME atlas"
)

# Gamma-delta T cells γδ T 細胞
gd_t <- tibble(
  Level1 = "Immune cell",
  Level2 = "Lymphoid cell",
  Level3 = "Gamma-delta T cells",
  Level4 = "Gamma-delta T cell",
  Gene   = c("TRDC","TRGC1","TRGC2","TRDV2","TRGV9","KLRD1","NKG7"),
  Source = "Pan-cancer TME atlas"
)

# Proliferating T cells 