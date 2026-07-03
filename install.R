# scONCO install.R
# ---------------------------------------------------------------
# Installs scONCO and its dependencies in a reproducible way.
# Run interactively:
#   source("install.R")
# Or from terminal:
#   Rscript install.R
# ---------------------------------------------------------------

cat("scONCO installer — v1.1\n")
cat("==========================\n\n")

# ---- helpers ----
ensure_pkg <- function(pkg, source = c("CRAN", "Bioc", "GitHub"), repo = NULL) {
  source <- match.arg(source)
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  v %s already installed (%s)\n", pkg,
                as.character(packageVersion(pkg))))
    return(invisible(TRUE))
  }
  cat(sprintf("  -> installing %s from %s ...\n", pkg, source))
  if (source == "CRAN") {
    install.packages(pkg, dependencies = TRUE)
  } else if (source == "Bioc") {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  } else if (source == "GitHub") {
    if (!requireNamespace("remotes", quietly = TRUE))
      install.packages("remotes")
    remotes::install_github(repo)
  }
}

# ---- R version check ----
if (getRversion() < "4.4.0") {
  stop("scONCO requires R >= 4.4.0; you have ", getRversion())
}

# ---- 1. CRAN deps ----
cran_deps <- c(
  "Seurat", "dplyr", "tibble", "tidyr", "stringr", "ggplot2",
  "Matrix", "testthat", "knitr", "rmarkdown", "pheatmap", "pROC",
  "easyPubMed", "devtools", "remotes", "BiocManager",
  # benchmark deps
  "optparse", "scCATCH", "SCINA", "aricode", "caret", "ggrepel",
  "future", "future.apply"
)
cat("[1/4] CRAN dependencies\n")
invisible(lapply(cran_deps, ensure_pkg, source = "CRAN"))

# ---- 2. Bioconductor deps ----
bioc_deps <- c("ComplexHeatmap", "SingleR", "SummarizedExperiment",
               "SingleCellExperiment", "celldex", "CelliD", "CHETAH",
               "infercnv")
cat("\n[2/4] Bioconductor dependencies\n")
invisible(lapply(bioc_deps, ensure_pkg, source = "Bioc"))

# ---- 3. GitHub deps ----
cat("\n[3/4] GitHub dependencies\n")
ensure_pkg("scMRMA",       source = "GitHub", repo = "JiaLiVUMC/scMRMA")
ensure_pkg("GPTCelltype",  source = "GitHub", repo = "Winnie09/GPTCelltype")
ensure_pkg("scATOMIC",     source = "GitHub", repo = "abelson-lab/scATOMIC")
ensure_pkg("copykat",      source = "GitHub", repo = "navinlabcode/copykat")

# ---- 4. scONCO itself ----
cat("\n[4/4] scONCO\n")
if (!requireNamespace("scONCO", quietly = TRUE)) {
  if (file.exists("DESCRIPTION")) {
    cat("  -> installing scONCO from local source ...\n")
    remotes::install_local(".", upgrade = "never")
  } else {
    cat("  -> installing scONCO from GitHub ...\n")
    remotes::install_github("scONCO-tool/scONCO", upgrade = "never")
  }
} else {
  cat(sprintf("  v scONCO already installed (%s)\n",
              as.character(packageVersion("scONCO"))))
}

cat("\nDone. Try:\n")
cat("  library(scONCO)\n")
cat("  ?scONCO::run_scONCO\n")
