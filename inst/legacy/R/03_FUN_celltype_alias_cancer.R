## ===============================================================
## KGDLab — Pan-cancer cell type normalization by regex (specific→general)
##
## Maps free-text / author-supplied cell-type labels onto scONCO canonical
## cell types. Rules are ordered SPECIFIC -> GENERIC: the first matching
## pattern wins, so malignant subtypes and TME subsets resolve before their
## broad parents.
## ===============================================================

suppressPackageStartupMessages({
  if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr")
  if (!requireNamespace("dplyr",   quietly = TRUE)) install.packages("dplyr")
})
library(stringr)
library(dplyr)


# ---------- 0) Pre-processing: normalise the surface form ----------
.normalize_label <- function(x) {
  x |>
    str_to_lower() |>
    str_replace_all("[_/\\-]+", " ") |>
    str_squish() |>
    str_replace_all("\\bcells\\b", "cell") |>
    str_replace_all("\\s*\\+\\s*", "+")
}

# ---------- 1) Rule table: specific -> generic (ORDER MATTERS!) ----------
regex_rules <- c(
  # ----- Malignant (tumour) cells: type-specific first -----
  "Malignant melanoma cell"     = "\\b(malignant|tumou?r|cancer|neoplastic)\\b.*\\bmelanoma\\b|\\bmelanoma\\s*(cell|tumou?r)\\b",
  "Malignant hepatocellular cell"="\\b(hcc|hepatocellular)\\b.*\\b(cell|carcinoma|tumou?r)\\b|\\bmalignant\\s*hepatocyte\\b",
  "Malignant lung adenocarcinoma cell" = "\\b(luad|lung\\s*adeno)\\w*\\b|\\badenocarcinoma\\b.*\\blung\\b",
  "Malignant lung squamous cell"= "\\b(lusc)\\b|\\blung\\b.*\\bsquamous\\b",
  "Malignant colorectal cell"   = "\\b(crc|colorectal|colon|rectal)\\b.*\\b(cell|carcinoma|tumou?r|epithel)\\b",
  "Malignant pancreatic cell"   = "\\b(pdac|pancrea\\w+)\\b.*\\b(cell|ductal|carcinoma|tumou?r)\\b",
  "Malignant gastric cell"      = "\\bgastric\\b.*\\b(cell|carcinoma|tumou?r)\\b",
  "Malignant ovarian cell"      = "\\bovarian\\b.*\\b(cell|carcinoma|tumou?r)\\b",
  "Malignant prostate cell"     = "\\bprostate\\b.*\\b(cell|carcinoma|tumou?r)\\b|\\bprad\\b",
  "Malignant glioma cell"       = "\\b(glioma|glioblastoma|gbm|astrocytoma)\\b",
  "Malignant renal cell"        = "\\b(rcc|renal\\s*cell)\\b|\\bclear\\s*cell\\s*renal\\b",
  "Malignant head and neck squamous cell" = "\\b(hnsc|hnscc|head\\s*and\\s*neck)\\b",
  "Malignant bladder urothelial cell" = "\\b(blca|bladder|urothelial)\\b.*\\b(cell|carcinoma|tumou?r)\\b",
  "Malignant breast cell"       = "\\b(brca|breast)\\b.*\\b(cell|carcinoma|tumou?r|epithel)\\b",
  "Cycling malignant cell"      = "\\b(cycling|proliferat\\w+|mitotic)\\b.*\\b(malignant|tumou?r|cancer)\\b",
  "EMT-like malignant cell"     = "\\bemt\\b.*\\b(malignant|tumou?r|cell)\\b|\\bmesenchymal[- ]like\\s*tumou?r\\b",
  "Malignant cell (generic)"    = "\\b(malignant|tumou?r|cancer|neoplastic|epithelial\\s*tumou?r)\\b.*\\bcell\\b|\\bmalignant\\b|\\btumou?r\\s*cell\\b",

  # ----- Normal / tissue-resident epithelial -----
  "Alveolar epithelial cells"   = "\\b(at1|at2|alveolar)\\b.*\\bepi(thel\\w+)?\\b|\\bat[12]\\b",
  "Hepatocytes"                 = "\\bhepatocytes?\\b",
  "Secretory/Goblet epithelial cells" = "\\b(goblet|secretory)\\b.*\\bcell\\b",
  "Epithelial cells"            = "\\bepithelial\\s*cell\\b|\\bepi\\b|\\bnormal\\s*epithel\\w+\\b",

  # ----- T / NK lymphoid (specific -> generic) -----
  "Regulatory T cells"          = "\\bt\\s*-?\\s*regs?\\b|\\bregulatory\\s*t\\b|\\bfoxp3\\+?\\s*t\\b",
  "CD8+ exhausted T cells"      = "\\bexhaust\\w+\\s*(cd8|t)\\b|\\b(cd8|t)\\b.*\\bexhaust\\w+\\b|\\bdysfunctional\\s*cd8\\b",
  "CD8+ cytotoxic T cells"      = "\\bcytotoxic\\s*(cd8|t)\\b|\\bctl\\b|\\beffector\\s*cd8\\b|\\bcd8\\b.*\\bgzmb\\b",
  "CD8+ T cells"                = "\\bcd8(?:\\+|\\s*positive)?\\s*t(?:\\s*cell)?\\b",
  "T follicular helper cells"   = "\\b(tfh|t\\s*follicular\\s*helper)\\b|\\bcxcl13\\+?\\s*t\\b",
  "Naive/memory T cells"        = "\\b(naive|memory|central\\s*memory)\\s*t\\b",
  "CD4+ T cells"                = "\\bcd4(?:\\+|\\s*positive)?\\s*t(?:\\s*cell)?\\b|\\bhelper\\s*t\\b",
  "MAIT cells"                  = "\\bmait\\b|\\bmucosal[- ]associated\\s*invariant\\b",
  "γδ T cells"                  = "\\b(gamma\\s*delta|gd|γδ)\\s*t\\b|\\btgd\\b",
  "Proliferating T cells"       = "\\b(cycling|proliferat\\w+)\\s*t\\b|\\bt\\b.*\\bmki67\\b",
  "NK cells"                    = "\\bnk\\s*cell\\b|\\bnatural\\s*killer\\b|\\bnk\\b",
  "T cells"                     = "\\bt\\s*cell\\b|\\bt\\s*lymphocyte\\b|\\btil\\b",

  # ----- B / Plasma -----
  "Plasma cells"                = "\\bplasma\\s*cell\\b|\\bantibody\\s*secreting\\b|\\bplasmablast\\b",
  "Germinal center B cells"     = "\\bgerminal\\s*cent\\w+\\s*b\\b|\\bgc\\s*b\\b",
  "B cells"                     = "\\bb\\s*cell\\b|\\bb\\s*lymphocyte\\b",

  # ----- Myeloid (specific -> generic) -----
  "TAM C1QC (immunoregulatory)" = "\\bc1qc\\b.*\\b(tam|macrophage)\\b|\\btam[_ ]c1qc\\b",
  "TAM SPP1 (angiogenic)"       = "\\bspp1\\b.*\\b(tam|macrophage)\\b|\\btam[_ ]spp1\\b",
  "TAM FOLR2 (tissue-resident-like)" = "\\bfolr2\\b.*\\b(tam|macrophage)\\b|\\blyve1\\+?\\s*macrophage\\b",
  "M1-like macrophages"         = "\\bm1[- ]?like\\b|\\bm1\\s*macrophage\\b",
  "M2-like macrophages"         = "\\bm2[- ]?like\\b|\\bm2\\s*macrophage\\b",
  "Mature regulatory DC (mregDC/LAMP3)" = "\\b(mregdc|lamp3\\+?\\s*dc)\\b|\\bmature\\s*regulatory\\s*dc\\b",
  "Conventional dendritic cells type 1 (cDC1)" = "\\bcdc1\\b|\\bclec9a\\+?\\s*dc\\b|\\bxcr1\\+?\\s*dc\\b",
  "Conventional dendritic cells type 2 (cDC2)" = "\\bcdc2\\b|\\bcd1c\\+?\\s*dc\\b",
  "Plasmacytoid dendritic cells"= "\\bpdc\\b|\\bplasmacytoid\\s*dendritic\\b",
  "Tumour-associated macrophage" = "\\btam\\b|\\btumou?r[- ]associated\\s*macrophage\\b",
  "Macrophages"                 = "\\bmacrophage\\b|\\bm(?:φ|phi)\\b|\\bkupffer\\b|\\bmicroglia\\b",
  "Monocytes"                   = "\\bmonocytes?\\b|\\bmono\\b|\\bcd14\\+?\\b",
  "Neutrophils (TAN)"           = "\\b(neutrophil|tan|pmn)\\b",
  "Mast cells"                  = "\\bmast\\s*cell\\b",

  # ----- Stromal: CAF / fibroblast / mural -----
  "Myofibroblastic CAF (myCAF)" = "\\bmycaf\\b|\\bmyofibroblastic\\s*caf\\b|\\bacta2\\+?\\s*caf\\b",
  "Inflammatory CAF (iCAF)"     = "\\bicaf\\b|\\binflammatory\\s*caf\\b|\\bil6\\+?\\s*caf\\b",
  "Antigen-presenting CAF (apCAF)" = "\\bapcaf\\b|\\bantigen[- ]presenting\\s*caf\\b",
  "Fibroblasts"                 = "\\b(caf|cancer[- ]associated\\s*fibroblast)\\b|\\bfibroblast\\b|\\bfb\\b",
  "Pericytes"                   = "\\bpericytes?\\b|\\bmural\\s*cell\\b|\\brgs5\\+?\\b",
  "Smooth muscle cells"         = "\\bsmooth\\s*muscle\\b|\\bsmc\\b|\\bvsmc\\b",

  # ----- Endothelial -----
  "Tip / angiogenic endothelial cells" = "\\btip\\s*(cell|ec)\\b|\\bangiogenic\\s*endothel\\w+\\b|\\besm1\\+?\\s*ec\\b",
  "Venous endothelial cells"    = "\\bvenous\\s*endothel\\w+\\b|\\backr1\\+?\\s*ec\\b",
  "Lymphatic endothelial cells" = "\\blymphatic\\s*endothel\\w+\\b|\\blec\\b|\\bprox1\\+?\\b",
  "Endothelial cells"           = "\\bendothelial\\s*cell\\b|\\bec\\b|\\bvascular\\s*endothel\\w+\\b",

  # ----- Dendritic generic (after specific DC subsets) -----
  "Conventional dendritic cells type 2 (cDC2)_generic" = "\\bdendritic\\s*cell\\b|\\bdc\\b",

  # ----- Other -----
  "Cycling cells (lineage-agnostic)" = "\\b(cycling|proliferat\\w+|mki67\\+?)\\s*cell\\b",
  "Erythrocytes"                = "\\b(erythrocyte|red\\s*blood\\s*cell|rbc)\\b"
)

# ---------- 2) Single-label normalisation ----------
normalize_celltype <- function(raw_label) {
  nm <- .normalize_label(raw_label)
  for (canonical in names(regex_rules)) {
    if (str_detect(nm, regex(regex_rules[[canonical]], ignore_case = TRUE))) {
      return(sub("_generic$", "", canonical))
    }
  }
  return(NA_character_)
}

# ---------- 3) Batch ----------
map_to_canonical <- function(labels) {
  vapply(labels, normalize_celltype, FUN.VALUE = character(1))
}

# ---------- 4) Self-test (pan-cancer example labels) ----------
.test <- function() {
  test_labels <- c(
    "Malignant melanoma cells", "Melanoma_tumor", "LUAD epithelial",
    "Tumor cells", "CD8+ exhausted T cells", "CD8-Tex", "Tregs",
    "myCAF", "iCAF", "SPP1+ TAM", "C1QC macrophage", "LAMP3+ DC",
    "pDC", "Tip cells", "Plasma_cells", "NK Cells", "Pericytes",
    "cancer-associated fibroblast"
  )
  data.frame(raw_label = test_labels,
             canonical = map_to_canonical(test_labels),
             stringsAsFactors = FALSE)
}
# print(.test())
