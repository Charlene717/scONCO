#!/usr/bin/env Rscript
###############################################################################
# sconco — command-line interface for scONCO v1.1
#
# Usage:
#   sconco annotate <seurat.rds> [--db v1.0] [--out result.rds]
#                    [--cluster-col seurat_clusters]
#                    [--with-confidence]
#   sconco db-info  [--version v1.0]
#   sconco version
#   sconco --help
#
# Examples:
#   sconco annotate path/to/preclustered.rds --with-confidence --out annot.rds
#   sconco db-info --version v_AI_Gemini
#   sconco version
#
# Installation:
#   # After installing the scONCO package
#   ln -s $(Rscript -e 'cat(system.file("cli", "sconco.R", package="scONCO"))') \
#         /usr/local/bin/sconco
#   chmod +x /usr/local/bin/sconco
###############################################################################

suppressPackageStartupMessages({
  library(optparse)
})

args <- commandArgs(trailingOnly = TRUE)

usage <- function() {
  cat("sconco — scONCO command-line interface\n\n")
  cat("Usage:\n")
  cat("  sconco annotate <seurat.rds> [options]\n")
  cat("  sconco db-info [options]\n")
  cat("  sconco version\n\n")
  cat("Subcommand options:\n")
  cat("  annotate:\n")
  cat("    --db VERSION           DB version [default v1.0]\n")
  cat("    --organism human|mouse [default human]\n")
  cat("    --out PATH             Output RDS file (annotated Seurat)\n")
  cat("    --cluster-col COL      Cluster column name [default seurat_clusters]\n")
  cat("    --with-confidence      Apply v1.1 reject + conformal\n")
  cat("  db-info:\n")
  cat("    --version VERSION      DB version to inspect [default v1.0]\n")
  cat("    --organism human|mouse [default human]\n")
}

if (length(args) == 0L) { usage(); quit(save = "no", status = 0L) }
sub <- args[1]; rest <- args[-1]

#------------------------------------------------------------------ version
if (sub == "version" || sub == "--version" || sub == "-v") {
  cat("scONCO v", as.character(utils::packageVersion("scONCO")), "\n", sep = "")
  quit(save = "no", status = 0L)
}

#------------------------------------------------------------------ db-info
if (sub == "db-info") {
  opt <- parse_args(OptionParser(option_list = list(
    make_option(c("--version"),  type = "character", default = "v1.0"),
    make_option(c("--organism"), type = "character", default = "human")
  )), args = rest)

  m <- scONCO::load_cancer_marker_db(version = opt$version, organism = opt$organism)
  cat(sprintf("scONCO DB: %s (%s)\n", opt$version, opt$organism))
  cat(sprintf("  unique genes:      %d\n", length(unique(m$gene))))
  cat(sprintf("  unique Level1:     %d\n", length(unique(m$Level1))))
  cat(sprintf("  unique Level4_Abb: %d\n", length(unique(m$Level4_Abb))))
  if ("sign" %in% colnames(m))
    cat(sprintf("  sign breakdown:    %s\n",
                paste(names(table(m$sign)), table(m$sign), sep = "=",
                      collapse = " "))
        )
  quit(save = "no", status = 0L)
}

#------------------------------------------------------------------ annotate
if (sub == "annotate") {
  if (length(rest) == 0L) {
    cat("Error: sconco annotate requires <seurat.rds> as first argument.\n\n")
    usage(); quit(save = "no", status = 1L)
  }
  rds_in <- rest[1]
  rest_opts <- rest[-1]
  opt <- parse_args(OptionParser(option_list = list(
    make_option(c("--db"),              type = "character", default = "v1.0"),
    make_option(c("--organism"),        type = "character", default = "human"),
    make_option(c("--out"),             type = "character", default = NULL),
    make_option(c("--cluster-col"),     type = "character", default = "seurat_clusters",
                dest = "cluster_col"),
    make_option(c("--with-confidence"), action = "store_true", default = FALSE,
                dest = "with_confidence")
  )), args = rest_opts)

  if (!file.exists(rds_in))
    stop("Input Seurat RDS not found: ", rds_in)

  cat("Loading Seurat object:", rds_in, "\n")
  seu <- readRDS(rds_in)

  cat("Loading marker DB:", opt$db, "(", opt$organism, ")\n")
  m <- scONCO::load_cancer_marker_db(version = opt$db, organism = opt$organism)

  cat("Running scONCO annotation...\n")
  seu <- scONCO::run_scONCO(
    seurat_obj       = seu,
    markers_df       = m,
    cluster_col      = opt$cluster_col,
    apply_confidence = opt$with_confidence
  )

  out_path <- if (is.null(opt$out))
    sub("\\.rds$", "_scONCO.rds", rds_in, ignore.case = TRUE) else opt$out
  saveRDS(seu, out_path)
  cat("Wrote:", out_path, "\n")
  quit(save = "no", status = 0L)
}

cat("Unknown subcommand:", sub, "\n\n"); usage(); quit(save = "no", status = 1L)
