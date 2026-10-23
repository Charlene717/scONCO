# scONCO

**Single-cell ONcology Reference-guided Marker-based Annotation**

[![R-CMD-check](https://github.com/Charlene717/scONCO/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Charlene717/scONCO/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

scONCO annotates the cells of a tumour sample along a four-level cell-type hierarchy, scored against a pan-cancer marker database in which every core marker is supported by at least three independent public sources. Markers are not asked to decide malignancy: the epithelial lineage is split into malignant and normal epithelium by a copy-number signal computed inside the sample against the sample's own immune cells, after a sample-level test of whether a tumour is present at all; when the test is negative, no cell is called malignant.

> **Status**: version 1.2.0 accompanies the manuscript (in preparation). `sconco()`, `sconco_study()` and `sconco_seurat()` are the method of the paper; they are generated from the project's `code/R/14-16` by `tools/sync_package.py`, and `tools/check_package_identity.R` confirms the package computes exactly what those files compute. The functions of the earlier flat marker-scoring framework (`run_scONCO()` and its confidence layer) are kept to reproduce the benchmark.

---

## Installation

```r
# install.packages("remotes")
remotes::install_github("Charlene717/scONCO")   # private until the paper is submitted
```

Or from source:

```r
remotes::install_local("path/to/scONCO_tool")
```

---

## Quick Start

```r
library(scONCO)

# one sample: raw counts, genes x cells (dgCMatrix or matrix)
res <- sconco(counts)
res                       # the sample verdict (tumour present or not) and cell counts per level
head(res$annotation)      # L1-L4, conf, status, cnv_burden, p_malignant per cell

# every sample of a study (cross-sample immune reference for samples with few immune cells)
all <- sconco_study(counts, sample = sample_id)

# on a Seurat object holding one sample: results go to its meta.data (sconco_L1 ... sconco_p_malignant)
seu <- sconco_seurat(seu)
```

The default database is version 3.3, shipped in `inst/extdata`; version 3.2 reproduces the initial version of the method: `sconco(counts, db = sconco_load_db(system.file("extdata", "DB_pancancer_human_v3.2.R", package = "scONCO")))`.

### The earlier flat marker-scoring framework (benchmark)

```r
markers_df <- load_cancer_marker_db(version = "v1.0", organism = "human")
seu <- run_scONCO(seurat_obj = seu, markers_df = markers_df, cluster_col = "seurat_clusters")
seu <- apply_hierarchical_reject(seu, thresholds = c(L1 = 0.7, L2 = 0.6, L3 = 0.55, L4 = 0.5))
seu <- apply_conformal_to_seurat(seu, alpha = 0.1, calibration_split = 0.2)
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

### 1.2.0 — the method of the paper
- `sconco()`, `sconco_study()`, `sconco_seurat()`, `sconco_reference_pool()`; building blocks `sconco_hier()` and `sconco_cnv()`
- Marker database v3.3 (default) and v3.2, Ensembl gene positions (`inst/extdata`)

### 1.0 / 1.1 — the earlier flat marker-scoring framework
- Specificity-weighted hierarchical scoring, cluster merging, auto-reclustering (modules 01–08)
- AUC-based marker specificity, hierarchical reject option, split conformal prediction (modules 09–12)

---

## Citation

If you use scONCO in published work, please cite:

> Chang C-J, Dai Y-H. scONCO: hierarchical cell-type annotation with a sample-level tumour test for single-cell RNA sequencing of solid tumours. In preparation.

BibTeX in `inst/CITATION`.

---

## Documentation

- `vignettes/scONCO-quickstart.Rmd` — quick start (10-minute tour).
- `docs/` (parent repository): algorithm pseudocode, DB schema, benchmark plan.

---

## Repository layout

```
scONCO_tool/
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

MIT. See `LICENSE`.

---

## Acknowledgements

- Dr. Yang-Hong Dai (NCKU) — radiation oncology / bioinformatics co-curation.
- Sharon — EPPK cohort sharing.
- National Cheng Kung University.
- Funded in part by NSTC Taiwan [grant ID].

We thank the developers of [scMRMA](https://github.com/JiaLiVUMC/scMRMA), [Seurat](https://satijalab.org/seurat/), [SingleR](https://bioconductor.org/packages/SingleR/), and [scCATCH](https://github.com/ZJUFanLab/scCATCH) — scONCO stands on their shoulders.
