# scONCO <img src="man/figures/logo.png" align="right" height="120"/>

**Single-cell ONcology Reference-guided Marker-based Annotation**

[![R-CMD-check](https://github.com/scONCO-tool/scONCO/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/scONCO-tool/scONCO/actions/workflows/R-CMD-check.yaml)
[![Codecov test coverage](https://codecov.io/gh/scONCO-tool/scONCO/branch/main/graph/badge.svg)](https://app.codecov.io/gh/scONCO-tool/scONCO)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.XXXXXXX.svg)](https://doi.org/10.5281/zenodo.XXXXXXX)

scONCO is a cancer cell-type annotation framework for scRNA-seq, built on the [scMRMA](https://github.com/JiaLiVUMC/scMRMA) hierarchical idea and augmented with an AI × Expert co-curated marker database and a specificity-weighted scoring algorithm with statistical confidence guarantees.

> ⚠️ **Status**: pre-release skeleton (v1.1). The functions in `R/` are placeholders that source the working scripts in `../code/R/`. Full `usethis::create_package()` + `devtools::document()` integration is the next milestone.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("scONCO-tool/scONCO")
```

Or from source:

```r
remotes::install_local("path/to/scONCO_package")
```

---

## Quick Start

```r
library(scONCO)
library(Seurat)

# 1. Load pan-cancer marker DB v1.0
markers_df <- load_cancer_marker_db(version = "v1.0", organism = "human")

# 2. Annotate a pre-clustered Seurat object
seu <- run_scONCO(
  seurat_obj   = seu,
  markers_df   = markers_df,
  cluster_col  = "seurat_clusters",
  hierarchy    = c("Level1", "Level2", "Level3", "Level4_Abb")
)

# 3. Apply v1.1 statistical confidence layer
seu <- apply_hierarchical_reject(seu, thresholds = c(L1 = 0.7, L2 = 0.6, L3 = 0.55, L4 = 0.5))
seu <- apply_conformal_to_seurat(seu, alpha = 0.1, calibration_split = 0.2)

# 4. Visualise
plot_marker_bubble(seu, markers_df, group_by = "scONCO_L4")
```

Results land in `seu@meta.data`:

| Column                         | Meaning                                  |
|--------------------------------|------------------------------------------|
| `scONCO_L1` … `scONCO_L4`    | Hierarchical labels at each level        |
| `scONCO_L{k}_confidence`      | Softmax confidence at level k            |
| `scONCO_L{k}_status`          | `assigned` / `deferred` per cell         |
| `conformal_set`                | Prediction set (list of plausible types) |
| `conformal_set_size`           | Size of the set                          |

---

## What's in this version

### v1.0 (modules 01–08)
- Specificity-weighted hierarchical scoring
- AI-curated marker DB (5 LLM modes + Expert curation)
- Iterative cluster merging by correlation and Jaccard
- Auto-reclustering (resolution gradient + ROGUE)

### v1.1 (modules 09–12, NEW)
- **A.1** AUC-based marker specificity (Pullin & McCarthy 2024)
- **A.2** Hierarchical reject option (de Lichtenberg et al. 2024)
- **A.3** Split conformal prediction (Boudaoud et al. 2025)
- **A.4** Empirical AI/literature support (replaces subjective scores)

### Coming in v2.0 (Advanced Methods companion paper)
- M1 Tier system, M2 Spatial regularization, M3 Trajectory constraint, …, M10 DB drift monitor

---

## Citation

If you use scONCO in published work, please cite:

> Chang, C-J., Dai, Y-H., et al. (2027). scONCO: a hierarchical, AI-augmented, statistically calibrated marker-based annotation framework for cancer single-cell transcriptomics. *Nucleic Acids Research* [in preparation].

BibTeX in `inst/CITATION`.

---

## Documentation

- `vignettes/scONCO-quickstart.Rmd` — quick start (10-minute tour).
- `vignettes/scONCO-v1_1-confidence.Rmd` — calibrated confidence with A.1–A.4.
- `vignettes/scONCO-custom-DB.Rmd` — how to extend the marker DB for new cell types.
- `docs/` (parent repository): algorithm pseudocode, DB schema, benchmark plan.

---

## Repository layout

```
scONCO_package/
├── DESCRIPTION
├── NAMESPACE
├── LICENSE
├── README.md
├── R/                         # All exported R functions (roxygen2-documented)
├── man/                       # Auto-generated Rd files
├── tests/testthat/            # Unit tests for v1.0 + v1.1 modules
├── inst/
│   ├── extdata/               # Reference DB snapshots (Zenodo-frozen versions)
│   └── CITATION
├── vignettes/                 # Long-form tutorials
├── data/                      # `markers_df` Rda for fast load
├── data-raw/                  # Scripts to regenerate `data/`
└── .github/workflows/         # CI (R-CMD-check, lintr, pkgdown)
```

---

## Development

```bash
# Inside the package directory
R -e "devtools::document()"
R -e "devtools::test()"
R -e "devtools::check()"
R -e "pkgdown::build_site()"
```

Continuous integration is configured in `.github/workflows/R-CMD-check.yaml`. Push triggers a multi-platform R-CMD-check + lintr + coverage upload to codecov.io.

---

## License

MIT (provisional; final license to be confirmed by PI before public release).

---

