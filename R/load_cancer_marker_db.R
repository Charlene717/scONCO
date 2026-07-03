#' Load the pan-cancer marker database (optionally a single cancer-type module)
#'
#' Loads the curated pan-cancer marker database shipped with scONCO. By default
#' it returns the full pan-cancer `markers_df` (malignant + immune + stromal +
#' endothelial + normal-epithelial compartments). A `module` can be supplied to
#' restrict the malignant compartment to a single cancer type while keeping the
#' shared tumour-microenvironment (TME) cell types.
#'
#' @param version Character. One of `"v1.0"` (default), or any of
#'   `"v_AI_ChatGPT5_DeepResearch"`, `"v_AI_ChatGPT5_Pro"`,
#'   `"v_AI_ChatGPT5_Thinking"`, `"v_AI_ChatGPT5_Instant"`,
#'   `"v_AI_Gemini"`, `"v_Expert"` to load a single source.
#' @param organism `"human"` (default) or `"mouse"`.
#' @param module Optional character. `"pan"` (default; everything) or one of the
#'   cancer-type modules: `"melanoma"`, `"LUAD"`, `"LUSC"`, `"BRCA"`, `"HCC"`,
#'   `"CRC"`, `"PDAC"`, `"STAD"`, `"OV"`, `"PRAD"`, `"GBM"`, `"RCC"`, `"HNSC"`,
#'   `"BLCA"`. When a specific module is given, malignant cell types from the
#'   other cancers are dropped, but the generic malignant programmes and all
#'   TME / normal compartments are kept.
#'
#' @return A data.frame with columns:
#'   * `gene` — HUGO symbol (or MGI for mouse)
#'   * `Level1` … `Level4_Abb` — hierarchical labels
#'   * `sign` — `"positive"` (default) or `"negative"`, where annotated
#'   * `ai_support` — fraction of AI sources agreeing (0–1)
#'   * `lit_support` — literature support score (0–1)
#'   * `specificity_index` — Wilcoxon AUC on a pan-cancer reference (0–1) if computed
#'
#' @examples
#' \dontrun{
#'   m   <- load_cancer_marker_db()                       # full pan-cancer, human
#'   mm  <- load_cancer_marker_db(organism = "mouse")     # mouse pan-cancer
#'   mel <- load_cancer_marker_db(module = "melanoma")    # melanoma + TME only
#' }
#'
#' @export
load_cancer_marker_db <- function(version  = "v1.0",
                                  organism = c("human", "mouse"),
                                  module   = "pan") {
  organism <- match.arg(organism)

  pkg_root  <- find.package("scONCO", quiet = TRUE)
  proj_root <- normalizePath(file.path(pkg_root, "..", ".."), mustWork = FALSE)
  extdata   <- system.file("extdata", package = "scONCO")  # DB bundled inside the package

  if (identical(version, "v1.0")) {
    fname <- if (organism == "human") "DB_pancancer_human_v1.0.R" else "DB_pancancer_mouse_v1.0.R"
    path  <- file.path(extdata, "current", fname)                                    # 1) installed package
    if (!file.exists(path))
      path <- file.path(proj_root, "database", "current", fname)                     # 2) dev-mode fallback
  } else {
    vers_dir <- file.path(extdata, "versions", version)                              # 1) installed package
    if (!dir.exists(vers_dir))
      vers_dir <- file.path(proj_root, "database", "versions", version)              # 2) dev-mode fallback
    if (!dir.exists(vers_dir))
      stop("Unknown DB version: ", version)
    files <- list.files(vers_dir, pattern = "\\.R$", full.names = TRUE)
    if (length(files) == 0L) stop("No .R files found in ", vers_dir)
    path <- files[which.max(file.info(files)$mtime)]
  }

  if (!file.exists(path))
    stop("Marker DB file not found: ", path)

  env <- new.env(parent = emptyenv())
  source(path, local = env)

  if (!exists("markers_df", envir = env))
    stop("Sourced DB file did not define `markers_df`: ", path)

  markers_df <- env$markers_df

  # --- Optional cancer-type module subsetting -------------------------------
  if (!identical(module, "pan")) {
    keep_tag       <- .sconco_module_tag(module)
    is_mal         <- markers_df$Level1 == "Malignant cell"
    is_generic_mal <- is_mal & grepl("pan-cancer", markers_df$Level2, fixed = TRUE)
    is_this_module <- is_mal & markers_df$Level2 == keep_tag
    # Keep: all non-malignant + generic malignant programmes + this module
    markers_df <- markers_df[!is_mal | is_generic_mal | is_this_module, , drop = FALSE]
    if (sum(is_this_module) == 0L)
      warning("Module '", module, "' matched no malignant cell types; ",
              "returning TME + generic malignant only.")
  }

  invisible(markers_df)
}

#' @keywords internal
.sconco_module_tag <- function(module) {
  map <- c(
    melanoma = "Malignant (melanoma)",
    LUAD     = "Malignant (LUAD)",
    LUSC     = "Malignant (LUSC)",
    BRCA     = "Malignant (BRCA)",
    HCC      = "Malignant (HCC)",
    CRC      = "Malignant (CRC)",
    PDAC     = "Malignant (PDAC)",
    STAD     = "Malignant (STAD)",
    OV       = "Malignant (OV)",
    PRAD     = "Malignant (PRAD)",
    GBM      = "Malignant (GBM/glioma)",
    RCC      = "Malignant (RCC)",
    HNSC     = "Malignant (HNSC)",
    BLCA     = "Malignant (BLCA)"
  )
  if (!module %in% names(map))
    stop("Unknown module '", module, "'. Valid modules: ",
         paste(names(map), collapse = ", "), ", or 'pan'.")
  unname(map[module])
}
