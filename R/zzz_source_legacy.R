# =============================================================
# scONCO — legacy bridge during skeleton phase
# =============================================================
# During the v1.1 skeleton phase, the canonical implementations
# live in ../code/R/. This file sources them on package attach so
# the exported function names resolve immediately. Once
# `devtools::document()` is run and each function has its own .R
# file with roxygen2 documentation, this file should be deleted.
# =============================================================

.sconco_legacy_dir <- function() {
  # Resolve the parent project's code/R/ directory.
  # Order of attempts:
  #   1. Environment variable SCONCO_LEGACY_DIR.
  #   2. ../code/R/ relative to the installed package.
  #   3. inst/legacy/R/ inside the installed package.
  env <- Sys.getenv("SCONCO_LEGACY_DIR", unset = "")
  if (nzchar(env) && dir.exists(env)) return(env)

  pkg <- find.package("scONCO", quiet = TRUE)
  cand <- normalizePath(file.path(pkg, "..", "..", "code", "R"),
                        mustWork = FALSE)
  if (dir.exists(cand)) return(cand)

  cand2 <- system.file("legacy", "R", package = "scONCO")
  if (nzchar(cand2)) return(cand2)

  return(NULL)
}

.onAttach <- function(libname, pkgname) {
  d <- .sconco_legacy_dir()
  if (is.null(d)) {
    packageStartupMessage(
      "scONCO v1.1 skeleton: legacy code/R/ directory not found.\n",
      "Set SCONCO_LEGACY_DIR=/path/to/code/R or copy modules to inst/legacy/R/."
    )
    return(invisible(NULL))
  }
  files <- list.files(d, pattern = "^[0-9]+_FUN_.*\\.R$",
                      full.names = TRUE)
  for (f in files) {
    tryCatch(
      source(f, local = topenv()),
      error = function(e) packageStartupMessage(
        "scONCO: failed to source ", basename(f), " — ", conditionMessage(e)
      )
    )
  }
  packageStartupMessage(
    "scONCO v1.1 loaded (", length(files), " legacy modules)."
  )
}
