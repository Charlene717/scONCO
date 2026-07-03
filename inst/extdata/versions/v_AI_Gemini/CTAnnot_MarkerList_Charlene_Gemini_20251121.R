###############################################################################
# scONCO — Cell-Type Annotation Speed Benchmark (Full: Cell + Cluster)
# Author: scONCO + GPT-5.1 Thinking
# Date: 2025-10-07
#
# Methods supported
#   - SingleR: HPCA/BPE (cell-level & cluster-level)
#   - scMRMA: single-level / multi-level (cell-level & cluster-level)
#   - CelliD: PanglaoDB gene sets (cell-level & cluster-level)
#   - scCATCH: cluster-level
#
# Usage
#   1) Load your Seurat object as `seuratObject_Sample` (DefaultAssay = "RNA").
#   2) Source your scMRMA marker DB script to create `marker_df_scMRMA_ChatGPT`.
#   3) source() this file.
#
# Outputs
#   ./scONCO_CTAnnotBench_5/
#     - <timestamp>_CTAnnot_Bench_results.csv
#     - <timestamp>_CTAnnot_Bench_summary.csv
#     - <timestamp>_BenchPlots.pdf / _BenchPlots.png
#     - <timestamp>_sessionInfo.txt
###############################################################################

suppressPackageStartupMessages({
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  # Base CRAN
  for (p in c("Seurat","tidyverse","ggplot2","patchwork","scales","R.utils","devtools")) {
    if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos="https://cloud.r-project.org")
  }
  # Bioc
  for (p in c("SingleR","celldex","BiocParallel")) {
    if (!requireNamespace(p, quietly = TRUE)) BiocManager::install(p, ask=FALSE, update=FALSE)
  }
  # Optional methods
  if (!requireNamespace("scMRMA", quietly = TRUE)) devtools::install_github("JiaLiVUMC/scMRMA", upgrade="never", quiet=TRUE)
  if (!requireNamespace("CelliD", quietly = TRUE)) BiocManager::install("CelliD", ask=FALSE, update=FALSE)
  if (!requireNamespace("scCATCH", quietly = TRUE)) install.packages("scCATCH", repos="https://cloud.r-project.org")
})

suppressPackageStartupMessages({
  library(Seurat)
  library(tidyverse)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(R.utils)
  library(devtools)
  library(SingleR)
  library(celldex)
  library(BiocParallel)
  library(scMRMA)
})

# Optional libraries (CelliD, scCATCH) are loaded on demand in wrappers.

###############################################################################
# USER CONFIG
###############################################################################
if (!exists("Name_ExportFolder")) Name_ExportFolder <- getwd()
Name_ExportFolder_CTAnnot <- file.path(Name_ExportFolder, "scONCO_CTAnnotBench_5")
if (!dir.exists(Name_ExportFolder_CTAnnot)) dir.create(Name_ExportFolder_CTAnnot, recursive = TRUE)
Name_Export <- format(Sys.time(), "%Y%m%d%H%M")

# Methods to run
methods_to_run <- c(
  "SingleR_HPCA","SingleR_BPE",
  "SingleR_HPCA_cluster","SingleR_BPE_cluster",
  "scMRMA_singleL",           # "scMRMA_multiL",
  "scMRMA_singleL_cluster",   # "scMRMA_multiL_cluster",
  "CelliD_panglao",          # CelliD cell-level (Panglao)
  "CelliD_panglao_cluster",  # CelliD cluster-level (Panglao)
  "scCATCH_cluster"          # scCATCH cluster-level
)

# Toggles
enable_CelliD  <- TRUE
enable_scCATCH <- TRUE

# Subsample sizes & reps
subsample_sizes <-   c(400, 1000, 10000, 20000, 30000, 50000, 80000)# c(400,1000, 5000, 8000, 10000) # c(100, 400,600) #  c(400,1000, 5000, 6000, 8000, 10000,12000)
n_reps <- 2
random_seed <- 42

# Timeouts
per_run_timeout_sec <- 36000  # 1h per run; set Inf to disable

# Parallel
Set_BPPARAM <- BiocParallel::bpparam()

# scMRMA params
Set_scMRMA_P   <- 0.05
Set_SubClust_k <- 20

# Ensure marker DB present
if (!exists("marker_df_scMRMA_ChatGPT")) {
  warning("`marker_df_scMRMA_ChatGPT` not found. Please source your marker DB file before running scMRMA.")
}
marker_df_scMRMA_ChatGPT_M <- NULL
if (exists("marker_df_scMRMA_ChatGPT")) {
  if ("Level4" %in% colnames(marker_df_scMRMA_ChatGPT)) {
    marker_df_scMRMA_ChatGPT_M <- marker_df_scMRMA_ChatGPT[,c(1,5)]
    colnames(marker_df_scMRMA_ChatGPT_M)[2] <- "LevelS"
  } else {
    warning("`marker_df_scMRMA_ChatGPT` does not have Level1–Level4 columns. scMRMA_multiL will fail.")
  }
}

###############################################################################
# PREP: Seurat object checks
###############################################################################
if (!exists("seuratObject_Sample")) stop("Please load `seuratObject_Sample` before sourcing this script.")
DefaultAssay(seuratObject_Sample) <- "RNA"
try({ seuratObject_Sample <- JoinLayers(seuratObject_Sample) }, silent = TRUE)
if (!("seurat_clusters" %in% colnames(seuratObject_Sample@meta.data))) {
  message("No `seurat_clusters` column found; creating a dummy single cluster (OK for timing).")
  seuratObject_Sample$seurat_clusters <- factor("0")
} else {
  seuratObject_Sample$seurat_clusters <- as.factor(seuratObject_Sample$seurat_clusters)
}

###############################################################################
# Preload references OUTSIDE timing
###############################################################################
message("Preloading celldex references (HPCA/BPE) outside timings...")
Ref_hpca.se <- celldex::HumanPrimaryCellAtlasData()
Ref_bpe.se  <- celldex::BlueprintEncodeData()

###############################################################################
# Helpers
###############################################################################
subsample_seurat <- function(obj, n_cells) {
  n_cells <- min(n_cells, ncol(obj))
  set.seed(random_seed)
  sel <- sample(colnames(obj), n_cells)
  sub_obj <- subset(obj, cells = sel)
  DefaultAssay(sub_obj) <- "RNA"
  try({ sub_obj <- JoinLayers(sub_obj) }, silent = TRUE)
  if (!("seurat_clusters" %in% colnames(sub_obj@meta.data))) {
    sub_obj$seurat_clusters <- factor("0")
  } else {
    sub_obj$seurat_clusters <- as.factor(sub_obj$seurat_clusters)
  }
  sub_obj
}

time_with_timeout <- function(expr, timeout_sec = per_run_timeout_sec) {
  gc()
  elapsed <- NA_real_
  res <- NULL
  timed <- try({
    elapsed <- system.time({
      res <- R.utils::withTimeout({
        force(expr)
      }, timeout = timeout_sec, onTimeout = "error")
    })[["elapsed"]]
  }, silent = TRUE)
  if (inherits(timed, "try-error")) {
    list(result = NULL, elapsed = NA_real_, error = TRUE)
  } else list(result = res, elapsed = as.numeric(elapsed), error = FALSE)
}

###############################################################################
# Benchmark wrappers
###############################################################################
bench_SingleR <- function(obj, which_ref = c("HPCA","BPE")) {
  which_ref <- match.arg(which_ref)
  ref_se <- if (which_ref == "HPCA") Ref_hpca.se else Ref_bpe.se
  counts <- Seurat::GetAssayData(obj)
  time_with_timeout({
    SingleR::SingleR(test = counts, ref = ref_se, assay.type.test = 1,
                     labels = ref_se$label.main, BPPARAM = Set_BPPARAM)
  })
}

bench_SingleR_cluster <- function(obj, which_ref = c("HPCA","BPE")) {
  which_ref <- match.arg(which_ref)
  ref_se <- if (which_ref == "HPCA") Ref_hpca.se else Ref_bpe.se
  Idents(obj) <- "seurat_clusters"
  counts <- Seurat::GetAssayData(obj)
  time_with_timeout({
    SingleR::SingleR(test = counts, ref = ref_se, assay.type.test = 1,
                     labels = ref_se$label.main, clusters = Idents(obj),
                     BPPARAM = Set_BPPARAM)
  })
}

# scMRMA: single/multi x (cell/cluster)
bench_scMRMA <- function(obj, mode = c("single","multi"), cluster_based = FALSE) {
  mode <- match.arg(mode)
  DefaultAssay(obj) <- "RNA"
  obj$seurat_clusters <- as.factor(obj$seurat_clusters)
  
  # choose DB
  db_use <- if (mode == "single") marker_df_scMRMA_ChatGPT_M else marker_df_scMRMA_ChatGPT
  if (is.null(db_use)) stop("scMRMA marker DB not prepared.")
  
  if (cluster_based) {
    # cluster-based: per-cluster label
    Idents(obj) <- "seurat_clusters"
    time_with_timeout({
      scMRMA::scMRMA(
        input = obj, p = Set_scMRMA_P, normalizedData = FALSE,
        selfDB = db_use, selfClusters = Idents(obj), k = Set_SubClust_k
      )
    })
  } else {
    # cell-level
    time_with_timeout({
      scMRMA::scMRMA(
        input = obj, p = Set_scMRMA_P, normalizedData = FALSE,
        selfDB = db_use, selfClusters = NULL, k = Set_SubClust_k
      )
    })
  }
}

# CelliD (cell-level): Panglao gene sets + RunMCA + RunCellHGT
bench_CelliD <- function(obj) {
  if (!enable_CelliD) return(list(result=NULL, elapsed=NA_real_, error=TRUE))
  DefaultAssay(obj) <- "RNA"
  # Load on demand
  if (!requireNamespace("CelliD", quietly = TRUE)) stop("CelliD not installed.")
  suppressPackageStartupMessages(library(CelliD))
  
  # Pre-download PanglaoDB (robust to retries)
  panglao_try <- try({
    Ref_panglao <- readr::read_tsv(
      "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz",
      progress = FALSE, show_col_types = FALSE
    )
  }, silent = TRUE)
  if (inherits(panglao_try, "try-error")) {
    return(list(result=NULL, elapsed=NA_real_, error=TRUE))
  } else {
    Ref_panglao <- panglao_try
  }
  Ref_panglao_all <- Ref_panglao %>% dplyr::filter(stringr::str_detect(species, "Hs"))
  Ref_panglao_all_gs <- Ref_panglao_all %>%
    dplyr::group_by(`cell type`) %>%
    dplyr::summarise(`cell type` = dplyr::first(`cell type`),
                     geneset = list(unique(`official gene symbol`)),
                     .groups = "drop") %>%
    dplyr::relocate(`cell type`)
  Ref_panglao_all_gs_Flt <- setNames(Ref_panglao_all_gs$geneset, Ref_panglao_all_gs$`cell type`)
  Ref_panglao_all_gs_Flt <- Ref_panglao_all_gs_Flt[sapply(Ref_panglao_all_gs_Flt, length) >= 10]
  
  time_with_timeout({
    obj <- CelliD::RunMCA(obj, nmcs = 50)
    mat <- CelliD::RunCellHGT(obj, pathways = Ref_panglao_all_gs_Flt, dims = 1:50)
    pred <- rownames(mat)[apply(mat, 2, which.max)]
    list(pred = pred, max_score = apply(mat, 2, max))
  })
}

# CelliD (cluster-level): Panglao gene sets + RunMCA + RunCellHGT + cluster-average
bench_CelliD_cluster <- function(obj) {
  if (!enable_CelliD) return(list(result = NULL, elapsed = NA_real_, error = TRUE))
  DefaultAssay(obj) <- "RNA"
  
  if (!requireNamespace("CelliD", quietly = TRUE)) stop("CelliD not installed.")
  suppressPackageStartupMessages(library(CelliD))
  
  # Pre-download PanglaoDB
  panglao_try <- try({
    Ref_panglao <- readr::read_tsv(
      "https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz",
      progress = FALSE, show_col_types = FALSE
    )
  }, silent = TRUE)
  if (inherits(panglao_try, "try-error")) {
    return(list(result = NULL, elapsed = NA_real_, error = TRUE))
  } else {
    Ref_panglao <- panglao_try
  }
  
  Ref_panglao_all <- Ref_panglao %>%
    dplyr::filter(stringr::str_detect(species, "Hs"))
  Ref_panglao_all_gs <- Ref_panglao_all %>%
    dplyr::group_by(`cell type`) %>%
    dplyr::summarise(
      `cell type` = dplyr::first(`cell type`),
      geneset     = list(unique(`official gene symbol`)),
      .groups     = "drop"
    ) %>%
    dplyr::relocate(`cell type`)
  
  Ref_panglao_all_gs_Flt <- setNames(
    Ref_panglao_all_gs$geneset,
    Ref_panglao_all_gs$`cell type`
  )
  Ref_panglao_all_gs_Flt <- Ref_panglao_all_gs_Flt[
    sapply(Ref_panglao_all_gs_Flt, length) >= 10
  ]
  
  # ensure clusters
  if (!("seurat_clusters" %in% colnames(obj@meta.data))) {
    obj$seurat_clusters <- factor("0")
  } else {
    obj$seurat_clusters <- as.factor(obj$seurat_clusters)
  }
  Seurat::Idents(obj) <- "seurat_clusters"
  
  time_with_timeout({
    # 1) MCA
    obj <- CelliD::RunMCA(obj, nmcs = 50)
    # 2) pathways x cells matrix
    mat <- CelliD::RunCellHGT(obj, pathways = Ref_panglao_all_gs_Flt, dims = 1:50)
    
    cluster_ids     <- Seurat::Idents(obj)
    cluster_levels  <- levels(cluster_ids)
    cluster_ids_vec <- as.character(cluster_ids)
    cells_all       <- colnames(mat)
    
    if (length(cluster_ids_vec) != length(cells_all)) {
      stop("Length of cluster_ids does not match number of cells in CelliD matrix.")
    }
    
    cluster_scores <- matrix(
      NA_real_,
      nrow = nrow(mat),
      ncol = length(cluster_levels),
      dimnames = list(rownames(mat), cluster_levels)
    )
    
    for (cl in cluster_levels) {
      idx <- which(cluster_ids_vec == cl)
      if (length(idx) == 0) next
      if (length(idx) == 1L) {
        cluster_scores[, cl] <- mat[, idx, drop = FALSE]
      } else {
        cluster_scores[, cl] <- rowMeans(mat[, idx, drop = FALSE])
      }
    }
    
    pred_cluster <- rownames(cluster_scores)[apply(cluster_scores, 2, which.max)]
    
    list(
      pred      = pred_cluster,
      max_score = apply(cluster_scores, 2, max)
    )
  })
}

# scCATCH (cluster-level)
bench_scCATCH <- function(obj) {
  if (!enable_scCATCH) return(list(result=NULL, elapsed=NA_real_, error=TRUE))
  if (!requireNamespace("scCATCH", quietly = TRUE)) stop("scCATCH not installed.")
  suppressPackageStartupMessages(library(scCATCH))
  DefaultAssay(obj) <- "RNA"
  Idents(obj) <- "seurat_clusters"
  time_with_timeout({
    data.input <- Seurat::GetAssayData(obj, assay = "RNA", slot = "data")
    data.input <- scCATCH::rev_gene(data = data.input, data_type = "data",
                                    species = "Human", geneinfo = scCATCH::geneinfo)
    sobj <- scCATCH::createscCATCH(data = data.input, cluster = as.character(Idents(obj)))
    sobj <- scCATCH::findmarkergene(
      object = sobj, species = "Human",
      marker = scCATCH::cellmatch,
      tissue = c('Tumour','Tumour microenvironment','Lymph node','Peripheral blood','Bone marrow','Spleen'),
      use_method = "1"
    )
    sobj <- scCATCH::findcelltype(object = sobj)
    sobj
  })
}

###############################################################################
# Warm-up (optional; not recorded)
###############################################################################
invisible(try({
  small_obj <- subsample_seurat(seuratObject_Sample, min(100, ncol(seuratObject_Sample)))
  bench_SingleR(small_obj,"HPCA")
  bench_SingleR(small_obj,"BPE")
  if (!is.null(marker_df_scMRMA_ChatGPT_M)) bench_scMRMA(small_obj,"single",FALSE)
  if (exists("marker_df_scMRMA_ChatGPT"))    bench_scMRMA(small_obj,"multi",FALSE)
  if (enable_CelliD) {
    bench_CelliD(small_obj)
    bench_CelliD_cluster(small_obj)
  }
  if (enable_scCATCH) bench_scCATCH(small_obj)
}, silent = TRUE))

###############################################################################
# Benchmark loop
###############################################################################
sizes <- sort(unique(pmin(subsample_sizes, ncol(seuratObject_Sample))))
if (length(sizes) == 0) stop("No valid subsample sizes.")

results <- list(); row_id <- 1L
total_runs <- length(sizes) * n_reps * length(methods_to_run)
pb <- txtProgressBar(min = 0, max = total_runs, style = 3)
progress_count <- 0

for (n in sizes) {
  for (rep_i in seq_len(n_reps)) {
    sub_obj <- subsample_seurat(seuratObject_Sample, n)
    for (m in methods_to_run) {
      progress_count <- progress_count + 1; setTxtProgressBar(pb, progress_count)
      bench_out <- switch(m,
                          "SingleR_HPCA"           = bench_SingleR(sub_obj,"HPCA"),
                          "SingleR_BPE"            = bench_SingleR(sub_obj,"BPE"),
                          "SingleR_HPCA_cluster"   = bench_SingleR_cluster(sub_obj,"HPCA"),
                          "SingleR_BPE_cluster"    = bench_SingleR_cluster(sub_obj,"BPE"),
                          "scMRMA_singleL"          = bench_scMRMA(sub_obj,"single",FALSE),
                          # "scMRMA_multiL"           = bench_scMRMA(sub_obj,"multi",FALSE),
                          "scMRMA_singleL_cluster"  = bench_scMRMA(sub_obj,"single",TRUE),
                          # "scMRMA_multiL_cluster"   = bench_scMRMA(sub_obj,"multi",TRUE),
                          "CelliD_panglao"         = bench_CelliD(sub_obj),
                          "CelliD_panglao_cluster" = bench_CelliD_cluster(sub_obj),
                          "scCATCH_cluster"        = bench_scCATCH(sub_obj),
                          { warning(sprintf("Unknown method: %s", m)); list(result=NULL, elapsed=NA_real_, error=TRUE) }
      )
      results[[row_id]] <- data.frame(
        method      = m,
        n_cells     = n,
        replicate   = rep_i,
        elapsed_sec = bench_out$elapsed,
        error       = isTRUE(bench_out$error),
        stringsAsFactors = FALSE
      )
      row_id <- row_id + 1L
      rm(bench_out); gc()
    }
    rm(sub_obj); gc()
  }
}
close(pb)

results_df <- dplyr::bind_rows(results)

###############################################################################
# Summaries & Plots （時間改成分鐘，細胞與時間都不用 log）
###############################################################################

# 秒數轉成分鐘
results_df <- results_df %>%
  dplyr::mutate(
    elapsed_min = elapsed_sec / 60
  )

summary_df <- results_df %>%
  dplyr::group_by(method, n_cells) %>%
  dplyr::summarise(
    n_ok       = sum(!error, na.rm = TRUE),
    n_failed   = sum(error, na.rm = TRUE),
    mean_min   = mean(elapsed_min, na.rm = TRUE),
    median_min = median(elapsed_min, na.rm = TRUE),
    sd_min     = sd(elapsed_min, na.rm = TRUE),
    .groups    = "drop"
  ) %>%
  dplyr::mutate(
    Type = ifelse(grepl("cluster", method), "Cluster-based", "Cell-level")
  )

# Save CSVs
csv_path <- file.path(Name_ExportFolder_CTAnnot, paste0(Name_Export, "_CTAnnot_Bench_results.csv"))
readr::write_csv(results_df, csv_path)

summary_csv_path <- file.path(Name_ExportFolder_CTAnnot, paste0(Name_Export, "_CTAnnot_Bench_summary.csv"))
readr::write_csv(summary_df, summary_csv_path)

# Plots（全部用線性刻度，不取 log）
p_median <- ggplot(summary_df, aes(x = n_cells, y = median_min, color = method, linetype = Type)) +
  geom_line() +
  geom_point() +
  labs(
    x = "Number of cells",
    y = "Median elapsed time (minutes)",
    title = "Cell-type Annotation Speed Benchmark (Cell vs Cluster)"
  ) +
  theme_bw()

p_all <- ggplot(results_df, aes(x = n_cells, y = elapsed_min, color = method)) +
  geom_point(alpha = 0.6, position = position_jitter(width = 0.05, height = 0)) +
  labs(
    x = "Number of cells",
    y = "Elapsed time per run (minutes)",
    title = "Per-run timing (all replicates)"
  ) +
  theme_bw()

plot_pdf <- file.path(Name_ExportFolder_CTAnnot, paste0(Name_Export, "_BenchPlots.pdf"))
plot_png <- file.path(Name_ExportFolder_CTAnnot, paste0(Name_Export, "_BenchPlots.png"))

ggsave(plot_pdf, p_median + patchwork::plot_spacer() + p_all, width = 12, height = 8)
ggsave(plot_png, p_median + patchwork::plot_spacer() + p_all, width = 12, height = 8, dpi = 300)

# Print quick view
print(summary_df %>% arrange(method, n_cells))

# Session info
sink(file.path(Name_ExportFolder_CTAnnot, paste0(Name_Export, "_sessionInfo.txt")))
print(sessionInfo())
sink()

message("✅ Done. Results saved in: ", Name_ExportFolder_CTAnnot)


##################################################################################
##################################################################################
##################################################################################
library(dplyr)
library(ggplot2)

# ---------- 0) 小工具：淺/深色轉換（不需額外套件） ----------
lighten_hex <- function(col, factor = 0.18){   # ⭐ 調小亮化幅度
  rgb <- grDevices::col2rgb(col) / 255
  rgb <- pmin(1, rgb + (1 - rgb) * factor)
  grDevices::rgb(rgb[1], rgb[2], rgb[3])
}
darken_hex <- function(col, factor = 0.12){    # ⭐ 調小變暗幅度
  rgb <- grDevices::col2rgb(col) / 255
  rgb <- pmax(0, rgb * (1 - factor))
  grDevices::rgb(rgb[1], rgb[2], rgb[3])
}

# ---------- 1) 清理與改名 ----------
summary_df2 <- summary_df %>%
  # filter(method != "scMRMA_singleL") %>%
  mutate(method = ifelse(method == "scMRMA_singleL_cluster", "scONCO_cluster", method)) %>%
  mutate(method = ifelse(method == "scMRMA_singleL", "scONCO", method)) %>%
  mutate(family = sub("_cluster$", "", method))

# ---------- 2) 指定家族底色 ----------
okabe_ito <- c(
  "#009E73", "#56B4E9", "#E69F00", "#939294",
  "#fc5362", "#8ac96f"
)
families <- unique(summary_df2$family)

base_map <- setNames(rep(okabe_ito, length.out = length(families)), families)
if ("CelliD_panglao" %in% names(base_map)) base_map["CelliD_panglao"] <- "#009E73"

# ---------- 3) 生成 Cell/Cluster 淺深色（已反轉 + 深淺差縮小） ----------
method_levels <- unique(summary_df2$method)
method_colors <- sapply(method_levels, function(m){
  fam <- sub("_cluster$", "", m)
  base <- base_map[[fam]]
  
  if (grepl("_cluster$", m)) {
    lighten_hex(base, 0.18)   # ⭐ Cluster-based → 淺色（更接近）
  } else {
    darken_hex(base, 0.12)    # ⭐ Cell-level → 深色（更接近）
  }
})
method_colors <- setNames(method_colors, method_levels)

# ---------- 3b) 設定線型（反轉） ----------
linetype_map <- c(
  "Cell-level"     = "dashed",
  "Cluster-based"  = "solid"
)

# ---------- 4) 作圖 ----------
p_median <- ggplot(
  summary_df2,
  aes(x = n_cells, y = median_min, color = method, linetype = Type)
) +
  geom_line(size = 1.6) +
  geom_point(size = 4.2) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = linetype_map) +
  labs(
    x = "Number of cells",
    y = "Elapsed time (minutes)",
    title = "Cell-type Annotation Speed Benchmark"
  ) +
  theme_bw(base_size = 20) +
  theme(
    plot.title = element_text(size = 26, face = "bold"),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 22),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
    panel.border = element_rect(color = "black", fill = NA, size = 2),
    legend.key.size = unit(1.5, "lines"),
    legend.spacing.y = unit(0.4, "cm")
  )

p_median


################################################################################


##################################################################################
library(dplyr)
library(ggplot2)

# ---------- 0) 小工具：淺/深色轉換（不需額外套件） ----------
lighten_hex <- function(col, factor = 0.18){
  rgb <- grDevices::col2rgb(col) / 255
  rgb <- pmin(1, rgb + (1 - rgb) * factor)
  grDevices::rgb(rgb[1], rgb[2], rgb[3])
}
darken_hex <- function(col, factor = 0.12){
  rgb <- grDevices::col2rgb(col) / 255
  rgb <- pmax(0, rgb * (1 - factor))
  grDevices::rgb(rgb[1], rgb[2], rgb[3])
}

# ---------- 1) 清理與改名 ----------
summary_df2 <- summary_df %>%
  mutate(method = ifelse(method == "scMRMA_singleL_cluster", "scONCO_cluster", method)) %>%
  mutate(method = ifelse(method == "scMRMA_singleL", "scONCO", method)) %>%
  mutate(family = sub("_cluster$", "", method))

# ---------- 2) 指定家族底色 ----------
okabe_ito <- c(
  "#009E73", "#56B4E9", "#E69F00", "#939294",
  "#fc5362", "#8ac96f"
)
families <- unique(summary_df2$family)

base_map <- setNames(rep(okabe_ito, length.out = length(families)), families)
if ("CelliD_panglao" %in% names(base_map)) base_map["CelliD_panglao"] <- "#009E73"

# ---------- 3) 生成 Cell/Cluster 顏色 ----------
method_levels <- unique(summary_df2$method)
method_colors <- sapply(method_levels, function(m){
  fam <- sub("_cluster$", "", m)
  base <- base_map[[fam]]
  
  if (grepl("_cluster$", m)) {
    lighten_hex(base, 0.18)
  } else {
    darken_hex(base, 0.12)
  }
})
method_colors <- setNames(method_colors, method_levels)

# ---------- 4) 分開資料 ----------
summary_cell <- summary_df2 %>% filter(Type == "Cell-level")
summary_cluster <- summary_df2 %>% filter(Type == "Cluster-based")

# ---------- 5) 統一 Y 軸 ----------
y_limits <- c(0, 9)

# ---------- 6) 作圖：Cell-level ----------
p_cell <- ggplot(
  summary_cell,
  aes(x = n_cells, y = median_min, color = method, linetype = Type)
) +
  geom_line(size = 1.6) +
  geom_point(size = 4.2) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = c(
    "Cell-level" = "21", # "Cell-level" = "dashed",
    "Cluster-based" = "solid"
  )) +                                    # ⭐ 強制指定虛線/實線
  guides(
    color = guide_legend(order = 1),      # ⭐ Method 在上方
    linetype = guide_legend(order = 2)    # ⭐ Type 在下方
  ) +
  labs(
    x = "Number of cells",
    y = "Elapsed time (minutes)",
    title = "Cell-type Annotation Speed Benchmark (Cell-level)"
  ) +
  ylim(y_limits) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 22),
    panel.border = element_rect(color = "black", fill = NA, size = 2)
  )

# ---------- 7) 作圖：Cluster-based ----------
p_cluster <- ggplot(
  summary_cluster,
  aes(x = n_cells, y = median_min, color = method, linetype = Type)
) +
  geom_line(size = 1.6) +
  geom_point(size = 4.2) +
  scale_color_manual(values = method_colors) +
  scale_linetype_manual(values = c(
    "Cell-level" = "dashed",
    "Cluster-based" = "solid"
  )) +
  guides(
    color = guide_legend(order = 1),
    linetype = guide_legend(order = 2)
  ) +
  labs(
    x = "Number of cells",
    y = "Elapsed time (minutes)",
    title = "Cell-type Annotation Speed Benchmark (Cluster-based)"
  ) +
  ylim(y_limits) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 16),
    plot.title = element_text(size = 22, face = "bold"),
    axis.title = element_text(size = 22, face = "bold"),
    axis.text = element_text(size = 22),
    panel.border = element_rect(color = "black", fill = NA, size = 2)
  )


# 顯示

p_clust