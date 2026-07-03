#### To-do list ####
# -[]顏色參數改放到函數中
# -[]客製化設定Marker
# -[]加入自動化流程
# -[]整理範例並Comment out

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(stringr)
})

`%||%` <- function(a, b) if (!is.null(a)) a else b

## =========================================================
## ✅ 0) PARAMS：只改這區就好
## =========================================================

PARAM_THEME <- "brown"       # "brown" or "diverging"
PARAM_pt_size_range <- c(0.1, 4)

PARAM_brown_color_low  <- "#f2f2f2"
PARAM_brown_color_high <- "#8b0000"
PARAM_brown_panel_bg   <- "#e9e2d8"
PARAM_brown_outer_bg   <- "white"
PARAM_brown_na         <- "#f2f2f2"

PARAM_div_color_low      <- "#2166ac"
PARAM_div_color_mid      <- "white"
PARAM_div_color_high     <- "#b2182b"
PARAM_div_color_midpoint <- 0
PARAM_div_na             <- "white"
PARAM_div_panel_bg       <- "#e9e2d8"
PARAM_div_outer_bg       <- "white"

PARAM_add_gene_vlines   <- TRUE
PARAM_gene_vline_color  <- "grey70"
PARAM_gene_vline_size   <- 0.3

PARAM_add_marker_ct_separators <- TRUE
PARAM_marker_ct_sep_color <- "black"
PARAM_marker_ct_sep_size  <- 0.9

PARAM_x_text_angle     <- 90
PARAM_top_label_size   <- 3.2
PARAM_top_label_vjust  <- -1.2
PARAM_top_margin_pt    <- 60
PARAM_show_legend      <- TRUE


## =========================================================
## ✅ 1) 主 function（修正版：自動檢查 gene 是否存在、提供 debug、可選 strict）
## =========================================================
## =========================================================
## ✅ Safe version: no `.data`, force dplyr namespace
## =========================================================

`%||%` <- function(a, b) if (!is.null(a)) a else b

plot_marker_bubble_topN_allow_dup_gene <- function(
    obj,
    markers_df,
    group.by = "Cell_Type",
    assay = NULL,
    slot = "data",
    level5_col = "Level4",
    top_n = 5,
    
    # ordering
    y_order = c("sort", "as_is", "custom"),
    y_custom = NULL,
    marker_ct_order = c("as_is", "sort", "custom"),
    marker_ct_custom = NULL,
    
    # debug
    strict = TRUE,
    verbose = TRUE
) {
  y_order <- match.arg(y_order)
  marker_ct_order <- match.arg(marker_ct_order)
  
  ## ---- safety message for masking
  if ("package:plyr" %in% search()) {
    message("[WARN] package:plyr is attached. This function forces dplyr:: calls to avoid masking.")
  }
  
  ## ---- deps (use namespaces)
  if (!requireNamespace("Seurat", quietly = TRUE)) stop("Need Seurat.")
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Need dplyr.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Need ggplot2.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Need stringr.")
  
  ## ---- assay handling ----
  if (!is.null(assay)) Seurat::DefaultAssay(obj) <- assay
  assay_use <- assay %||% Seurat::DefaultAssay(obj)
  
  if (!assay_use %in% names(obj@assays)) {
    stop("Assay '", assay_use, "' not found. Available assays: ",
         paste(names(obj@assays), collapse = ", "))
  }
  
  if (!all(c("gene", level5_col) %in% colnames(markers_df))) {
    stop("markers_df must contain columns: gene and ", level5_col)
  }
  
  ## ---- normalize gene text (NO `.data`) ----
  markers_df <- dplyr::mutate(
    markers_df,
    gene = stringr::str_trim(as.character(.data$gene))  # 這裡 .data$gene 是安全的（固定欄位）
  )
  markers_df$.marker_ct <- as.character(markers_df[[level5_col]])
  
  ## 1) mapping
  map_raw <- markers_df |>
    dplyr::transmute(
      marker_ct     = .data$.marker_ct,
      gene          = .data$gene,
      priority_rank = dplyr::row_number()
    ) |>
    dplyr::filter(!is.na(.data$marker_ct), .data$marker_ct != "",
                  !is.na(.data$gene), .data$gene != "") |>
    dplyr::distinct(.data$marker_ct, .data$gene, .keep_all = TRUE)
  
  ## ---- genes in object ----
  genes_in_obj <- rownames(obj[[assay_use]])
  genes_in_obj_trim <- stringr::str_trim(genes_in_obj)
  
  map_raw$gene_match <- map_raw$gene
  present <- map_raw$gene_match %in% genes_in_obj_trim
  missing_genes <- unique(map_raw$gene_match[!present])
  
  if (isTRUE(verbose)) {
    message("[INFO] assay_use = ", assay_use)
    message("[INFO] markers total unique genes = ", length(unique(map_raw$gene_match)))
    message("[INFO] genes found in assay = ", sum(present), " (",
            round(100 * sum(present) / max(1, length(present)), 1), "%)")
    if (length(missing_genes) > 0) {
      message("[INFO] example missing genes (up to 20): ",
              paste(head(missing_genes, 20), collapse = ", "))
    }
  }
  
  map_raw <- map_raw |>
    dplyr::filter(.data$gene_match %in% genes_in_obj_trim) |>
    dplyr::mutate(gene = .data$gene_match) |>
    dplyr::select(-.data$gene_match)
  
  if (nrow(map_raw) == 0) {
    msg <- paste0(
      "No marker genes found in the Seurat object assay.\n",
      "Assay used: ", assay_use, "\n",
      "Likely reasons:\n",
      "  (1) markers_df gene symbols do not match rownames(obj[[", assay_use, "]]).\n",
      "  (2) Wrong assay.\n",
      "  (3) Ensembl vs symbol / species mismatch.\n"
    )
    if (isTRUE(strict)) stop(msg)
    return(list(
      plot = NULL, data = NULL, top_genes = NULL, ct_sep = NULL,
      error = msg, missing_genes = missing_genes, assay_used = assay_use
    ))
  }
  
  ## 2) DotPlot once
  dotplot_has_slot <- "slot" %in% names(formals(Seurat::DotPlot))
  dp_args <- list(
    object   = obj,
    features = unique(map_raw$gene),
    group.by = group.by,
    assay    = assay_use
  )
  if (dotplot_has_slot) dp_args$slot <- slot
  
  dp  <- do.call(Seurat::DotPlot, dp_args)
  df0 <- dp$data
  
  required_cols <- c("id", "features.plot", "avg.exp.scaled", "pct.exp")
  if (!all(required_cols %in% colnames(df0))) {
    stop("DotPlot output missing expected columns: ",
         paste(setdiff(required_cols, colnames(df0)), collapse = ", "))
  }
  df0$gene <- df0[["features.plot"]]
  
  ## 3) rank
  gene_rank <- df0 |>
    dplyr::inner_join(map_raw, by = "gene", relationship = "many-to-many") |>
    dplyr::group_by(.data$marker_ct, .data$gene, .data$priority_rank) |>
    dplyr::summarise(score = mean(.data$avg.exp.scaled, na.rm = TRUE), .groups = "drop") |>
    dplyr::arrange(.data$marker_ct, .data$priority_rank, dplyr::desc(.data$score))
  
  ## 4) top_n per marker_ct
  top_genes_df <- gene_rank |>
    dplyr::group_by(.data$marker_ct) |>
    dplyr::slice_head(n = top_n) |>
    dplyr::ungroup()
  
  if (nrow(top_genes_df) == 0) stop("No top genes selected; check expression or markers_df.")
  
  ## marker_ct order
  marker_ct_levels <- unique(top_genes_df$marker_ct)
  if (marker_ct_order == "sort") {
    marker_ct_levels <- rev(sort(marker_ct_levels))
  } else if (marker_ct_order == "custom") {
    if (is.null(marker_ct_custom)) stop("marker_ct_order='custom' requires marker_ct_custom.")
    marker_ct_levels <- c(marker_ct_custom, setdiff(marker_ct_levels, marker_ct_custom))
  }
  
  top_genes_df <- dplyr::mutate(
    top_genes_df,
    gene_key   = paste0(.data$marker_ct, "||", .data$gene),
    gene_label = .data$gene
  )
  
  gene_key_levels <- unlist(lapply(marker_ct_levels, function(ct){
    top_genes_df |>
      dplyr::filter(.data$marker_ct == ct) |>
      dplyr::arrange(.data$priority_rank, dplyr::desc(.data$score)) |>
      dplyr::pull(.data$gene_key)
  }), use.names = FALSE)
  
  df <- df0 |>
    dplyr::inner_join(
      top_genes_df |> dplyr::select(.data$marker_ct, .data$gene, .data$gene_key, .data$gene_label) |> dplyr::distinct(),
      by = "gene",
      relationship = "many-to-many"
    ) |>
    dplyr::filter(.data$gene_key %in% gene_key_levels) |>
    dplyr::mutate(
      gene_key = factor(.data$gene_key, levels = gene_key_levels),
      pct.exp  = pmax(0, pmin(100, .data$pct.exp))
    )
  
  ## Y order
  y_levels <- unique(as.character(df$id))
  if (y_order == "sort") {
    y_levels <- rev(sort(y_levels))
  } else if (y_order == "custom") {
    if (is.null(y_custom)) stop("y_order='custom' requires y_custom.")
    y_levels <- c(y_custom, setdiff(y_levels, y_custom))
  }
  df$id <- factor(as.character(df$id), levels = y_levels)
  
  ## label map
  x_labels <- top_genes_df |>
    dplyr::select(.data$gene_key, .data$gene_label) |>
    dplyr::distinct()
  x_label_map <- setNames(as.character(x_labels$gene_label), as.character(x_labels$gene_key))
  
  ## positions for ct label & separators
  x_pos <- data.frame(
    gene_key = factor(gene_key_levels, levels = gene_key_levels),
    xi = seq_along(gene_key_levels),
    stringsAsFactors = FALSE
  ) |>
    dplyr::left_join(
      top_genes_df |> dplyr::select(.data$marker_ct, .data$gene_key) |> dplyr::distinct(),
      by = "gene_key"
    )
  x_pos$marker_ct <- as.character(x_pos$marker_ct)
  
  ct_mid <- x_pos |>
    dplyr::group_by(.data$marker_ct) |>
    dplyr::summarise(
      x_mid = (min(.data$xi) + max(.data$xi)) / 2,
      x_end = max(.data$xi),
      .groups = "drop"
    ) |>
    dplyr::arrange(match(.data$marker_ct, marker_ct_levels))
  
  ## vlines data (這些 PARAM_* 仍沿用你原本全域參數)
  gene_vline_df <- NULL
  if (isTRUE(PARAM_add_gene_vlines)) {
    x_n <- length(gene_key_levels)
    if (x_n >= 2) gene_vline_df <- data.frame(xintercept = seq(1.5, x_n - 0.5, by = 1))
  }
  
  ct_sep_df <- NULL
  if (isTRUE(PARAM_add_marker_ct_separators)) {
    x_ends <- ct_mid$x_end
    if (length(x_ends) >= 2) ct_sep_df <- data.frame(xintercept = x_ends[-length(x_ends)] + 0.5)
  }
  
  ## theme params
  if (PARAM_THEME == "brown") {
    panel_bg <- PARAM_brown_panel_bg
    outer_bg <- PARAM_brown_outer_bg
    na_col   <- PARAM_brown_na
  } else if (PARAM_THEME == "diverging") {
    panel_bg <- PARAM_div_panel_bg
    outer_bg <- PARAM_div_outer_bg
    na_col   <- PARAM_div_na
  } else {
    stop("PARAM_THEME must be 'brown' or 'diverging'.")
  }
  
  ## base plot
  p <- ggplot2::ggplot(df, ggplot2::aes(x = gene_key, y = id)) +
    ggplot2::geom_point(ggplot2::aes(size = pct.exp, color = avg.exp.scaled), alpha = 0.95) +
    ggplot2::scale_size(range = PARAM_pt_size_range, limits = c(0, 100)) +
    ggplot2::scale_x_discrete(labels = x_label_map) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      panel.background   = ggplot2::element_rect(fill = panel_bg, color = NA),
      plot.background    = ggplot2::element_rect(fill = outer_bg, color = NA),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.title         = ggplot2::element_blank(),
      axis.text.x        = ggplot2::element_text(angle = PARAM_x_text_angle, hjust = 1, vjust = 0.5, color = "black"),
      axis.text.y        = ggplot2::element_text(color = "black"),
      plot.margin        = ggplot2::margin(t = PARAM_top_margin_pt, r = 10, b = 10, l = 10),
      legend.position    = if (isTRUE(PARAM_show_legend)) "bottom" else "none"
    ) +
    ggplot2::coord_cartesian(clip = "off")
  
  ## color scale
  if (PARAM_THEME == "brown") {
    p <- p + ggplot2::scale_color_gradient(
      low = PARAM_brown_color_low,
      high = PARAM_brown_color_high,
      na.value = na_col
    )
  } else {
    p <- p + ggplot2::scale_color_gradient2(
      low = PARAM_div_color_low,
      mid = PARAM_div_color_mid,
      high = PARAM_div_color_high,
      midpoint = PARAM_div_color_midpoint,
      na.value = na_col
    )
  }
  
  ## vlines
  if (isTRUE(PARAM_add_gene_vlines) && !is.null(gene_vline_df)) {
    p <- p + ggplot2::geom_vline(
      data = gene_vline_df,
      ggplot2::aes(xintercept = xintercept),
      color = PARAM_gene_vline_color,
      linewidth = PARAM_gene_vline_size,
      inherit.aes = FALSE
    )
  }
  if (isTRUE(PARAM_add_marker_ct_separators) && !is.null(ct_sep_df)) {
    p <- p + ggplot2::geom_vline(
      data = ct_sep_df,
      ggplot2::aes(xintercept = xintercept),
      color = PARAM_marker_ct_sep_color,
      linewidth = PARAM_marker_ct_sep_size,
      inherit.aes = FALSE
    )
  }
  
  ## ✅ top marker_ct labels
  p <- p + ggplot2::geom_text(
    data = ct_mid,
    ggplot2::aes(x = x_mid, y = Inf, label = marker_ct),
    vjust = PARAM_top_label_vjust,
    size  = PARAM_top_label_size,
    color = "black",
    inherit.aes = FALSE
  )
  
  list(
    plot = p, data = df, top_genes = top_genes_df, ct_sep = ct_sep_df,
    assay_used = assay_use, missing_genes = missing_genes
  )
}




## =========================================================
## ✅ 2) Run
## =========================================================
source("CTAnnot_MarkerList_Wilson_ChatGPT5_20251219.R")
source("CTAnnot_MarkerList_Charlene_ChatGPT5_2025092701.R")
source("CTAnnot_MarkerList_Charlene_ChatGPT5_20251209_20251228M.R")
markers_df <- marker_df_scMRMA_ChatGPT

## ---- base key from seurat ----
ct_key <- str_remove(as.character(seuratObject_Sample$Cell_Type), "_\\d+$")
map_ct_to_lv4 <- setNames(unique(ct_key), unique(ct_key))

level4_keep <- unique(map_ct_to_lv4[ct_key])
level4_keep <- level4_keep[!is.na(level4_keep) & level4_keep != ""]
# if ("SMC" %in% unique(as.character(markers_df$Level4))) {
#   level4_keep <- unique(c(level4_keep, "SMC","M1_like","M2_like","moMac"))
# }

# level4_keep <- unique(c(level4_keep, "SMC","Peri","M1_like","M2_like","MyoFB","Plasma","moMac"))
# level4_keep <- unique(c(level4_keep, "SMC","Peri","MyoFB","Plasma"))
level4_keep <- unique(c(level4_keep, "SMC","Peri","MyoFB", "KC_Basal", "KC_Spinous","KC_Granular","Plasma"))


markers_df <- markers_df %>% filter(Level4 %in% level4_keep)


## quick check
unique(markers_df$Level4)


## 使用棕色系：
PARAM_THEME <- "brown"
res_brown <- plot_marker_bubble_topN_allow_dup_gene(
  obj = seuratObject_Sample,
  markers_df = markers_df,
  group.by = "Cell_Type",
  level5_col = "Level4",
  top_n = 10,
  y_order = "sort"
)
res_brown$plot

## 使用正負色（0白）：
PARAM_THEME <- "diverging"
res_div <- plot_marker_bubble_topN_allow_dup_gene(
  obj = seuratObject_Sample,
  markers_df = markers_df,
  group.by = "Cell_Type",
  level5_col = "Level4",
  top_n = 10,
  y_order = "sort"
)
res_div$plot


## =========================================================
## ✅ 3) 免疫細胞版本（外觀同樣只改最上面 PARAM_*）
## =========================================================


## ---- Immune subset ----
immune_types <- c("Tcell_res","NK/Tcyt_res","M1_like","M2_like","moMac","DC_p","Mastcell")
markers_df_immune <- subset(markers_df, Level4 %in% immune_types)

## =========================================================
## ✅ 3) Plot
## =========================================================

PARAM_THEME <- "brown"
res_immune <- plot_marker_bubble_topN_allow_dup_gene(
  obj        = seuratObject_Sample,
  markers_df = markers_df_immune,
  group.by   = "Cell_Type",
  level5_col = "Level4",
  top_n      = 15,
  y_order    = "sort",
  marker_ct_order = "as_is",
  strict     = FALSE,   # ✅ 不要直接停掉，先把 debug 資訊吐出來
  verbose    = TRUE
)

## 如果真的沒有任何 gene match，這裡會顯示原因與 missing gene
if (is.null(res_immune$plot)) {
  cat(res_immune$error, "\n")
  cat("Example missing genes:\n")
  print(head(res_immune$missing_genes, 30))
} else {
  res_immune$plot
}
