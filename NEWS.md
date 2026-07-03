# scONCO — NEWS

Reverse-chronological changelog. Date format: YYYY-MM-DD.

---

## v1.1.0 (2026-05-29, in development)

### New features

* **A.1 — AUC-based marker specificity** (module `09_FUN_marker_specificity_AUC.R`)
  - `compute_marker_AUC_batch()`: Wilcoxon AUC per (gene, cell-type) on a reference atlas.
  - `compute_tau_score()`, `compute_mutex_score()`: concentration + mutual-exclusivity metrics (Wang et al. 2020 SCMarker).
  - `update_specificity_empirical()`: refresh `markers_df$specificity_index` from new reference data without retraining.
  - References: Pullin & McCarthy (2024) *Genome Biology* 25, 56; Wang et al. (2020) *PLOS Comp Biol*.

* **A.2 — Hierarchical reject option** (module `10_FUN_hierarchical_reject.R`)
  - `compute_confidence()` / `hierarchical_reject_one()` per-cell softmax confidence.
  - `apply_hierarchical_reject()`: defer ambiguous cells (confidence < threshold_k) to parent level.
  - `calibrate_reject_thresholds()`: sweep thresholds on a held-out set to hit target reject rate.
  - Adds `scONCO_L{k}_confidence` and `scONCO_L{k}_status` ∈ {assigned, deferred} columns.
  - Reference: de Lichtenberg et al. (2024) *Bioinformatics* 40(3), btae128.

* **A.3 — Split conformal prediction** (module `11_FUN_conformal_prediction.R`)
  - `split_conformal_calibrate()`, `conformal_predict()`: marginal coverage guarantee at user-chosen α.
  - `apply_conformal_to_seurat()`: writes `conformal_set` + `conformal_set_size` to meta.data.
  - Marginal coverage guarantee: `P(true_label ∈ C(cell)) ≥ 1 − α`, distribution-free.
  - Reference: Boudaoud et al. (2025) *Bioinformatics* 41(10), btaf521.

* **A.4 — Empirical AI/literature support** (module `12_FUN_empirical_support.R`)
  - `compute_empirical_ai_support()`: counts fraction of AI source DBs that include each marker.
  - `compute_empirical_lit_support()`: PubMed co-occurrence score via `easyPubMed`.
  - `refresh_marker_weights_empirical()`: replace subjective `ai_support` / `lit_support` columns.
  - `upgrade_to_scONCO_v1_1()`: one-call wrapper that applies A.1–A.4 in sequence.

* **`run_scONCO()` top-level entry**
  - Single function for the full pipeline; new `apply_confidence` argument toggles the A.1–A.4 layer.
  - CLI counterpart: `sconco annotate <seurat.rds> --with-confidence --out annotated.rds`.

* **DB schema additions**
  - `sign` column (positive / negative) in `markers_df` — enables the δ negative-marker penalty.
  - Negative markers filled via `database/scripts/01_add_negative_marker_column.R` from a curator template (Dr. Dai).
  - OBO + JSON exports via `database/scripts/03_export_obo_json.R` for external Zenodo / HCA citation.
  - Cross-version consensus (Jaccard) quantification via `database/scripts/02_compute_consensus_overlap.R`.

* **R package skeleton**
  - DESCRIPTION / NAMESPACE / roxygen2 docs / testthat smoke tests / vignettes / inst/CITATION.
  - install.R one-call dependency setup (CRAN + Bioconductor + GitHub).
  - renv.lock for fully reproducible installs.
  - Dockerfile + GitHub Actions docker-publish workflow.
  - `inst/cli/sconco.R` shell-callable CLI (`annotate` / `db-info` / `version` subcommands).

### Benchmark infrastructure (new)

* `benchmarks/tools/` — 9 unified tool wrappers conforming to `00_wrapper_interface.R` contract:
  scONCO, scMRMA, scCATCH, CellID, SCINA (marker); SingleR, Seurat_LT, CHETAH (reference); GPTCelltype (LLM).

* `benchmarks/references/` — 4 reference loaders (Reynolds 2021, celldex HPCA, Hao 2021 PBMC, Joost 2020 mouse) sharing `00_reference_interface.R`.

* `benchmarks/datasets/` — added 5 tier-1 loaders (Reynolds 2021, Hughes 2020 psoriasis, Rojahn 2020 AD, Tabib 2021 SSc, Joost 2020 mouse) plus the existing 3.

* `benchmarks/run_all.R` — master execution with 3 sub-plans union (marker / reference / llm by tool category), CLI for `--datasets / --tools / --db / --reference / --parallel / --dry-run / --results`. Auto-filters cross-species mismatches.

* `benchmarks/figures/` — Fig 2 (DB axis), Fig 3 (Tool axis), Fig 4 (Ablation), Fig 5 (Coverage), Fig 6 (Reference axis) ggplot skeletons + `make_all.R` orchestrator. Stub data fallback so layouts preview even before benchmarks run.

* `benchmarks/00_self_check.R` — pre-flight environment checklist (9 sections).
* `benchmarks/smoke_test_synthetic.R` — 19 sandbox checks (PASS), no Seurat / no real data.
* `benchmarks/smoke_test_GSE163973.R` — first end-to-end real-dataset smoke test.

### Bug fixes

* `run_all.R::.resolve_script_dir()` — replaced fragile `sys.frame(1)$ofile` with a 4-fallback resolver
  that works under `Rscript`, `R -f`, `source()`, env-var override, and interactive.

### Verified in agent sandbox (R 4.1.2)

* 74 / 74 .R files parse clean.
* 9 / 9 tool wrappers' signatures match the interface contract.
* 19 / 19 checks pass in `smoke_test_synthetic.R`.
* 5 dry-run plans verified, including kitchen-sink 354 combinations (280 marker + 66 reference + 8 llm).
* Cross-species filter drops 30 organism-mismatches correctly.

---

## v1.0.0 (2026-05-10)

Initial public release of the scONCO framework.

### Core algorithm (modules 01–08)

* `01_FUN_scONCO_core.R`: specificity-weighted hierarchical scoring
  `weight(g) = α · spec + β · ai_support + γ · lit_support`
  `score(cell, type) = Σ weight(g) · expr(g, cell) − δ · neg_penalty(cell, type)`
* `02_FUN_annotate_broad_clusters.R`: L4 → broad lineage mapping (5 lineages).
* `03_FUN_celltype_alias_cancer.R`: regex normalisation of cancer cell-type names.
* `04_FUN_marker_dotplots.R` + `05_FUN_marker_bubble.R`: visualisation.
* `06_FUN_auto_reclustering.R`: ROGUE-guided resolution search.
* `07_FUN_merge_clusters_correlation.R`: merge clusters by Spearman ρ on top markers.
* `08_FUN_merge_clusters_marker.R`: merge clusters by marker Jaccard.

### Marker database

* `database/current/DB_pancancer_human_v1.0.R` — ~539 markers × 61 cell types × 4 hierarchy levels (pan-cancer).
* `database/versions/` — 9 source variants (5 AI modes × 2 versions + Expert).
* `scONCO_PanCancer_Marker_Comprehensive.xlsx` — 10-sheet curated DB overview.

### Documentation

* Project README + `code/README.md` + `database/README.md` + `literature/README.md`.
* `manuscript/scONCO_NSTC_Proposal_v3_EndNote.docx` — NSTC proposal (3-year program).
* `manuscript/scONCO_Paper_Outline.docx` — initial NAR paper outline (superseded by v2 main + advanced in 2026-05-29).
* `presentations/scONCO_Progress_PaperPlanning_v1.pptx` — 22-slide milestone deck.

---

## Roadmap

### v1.2 (planned, ~2026-08)
* Mouse pan-cancer marker DB (`DB_pancancer_mouse_v1.0.R`) anchored on Joost 2020 GSE129218.
* Bioconductor first submission.

### v1.5 (planned, ~2026-10)
* popV-style ensemble vote across scONCO + scMRMA + SingleR + GPTCelltype.
* HCE loss + SOCAM calibration on softmax confidence.

### v2.0 (planned, ~2027 Q2; see Advanced Methods companion paper)
* M1 Marker Importance Tier System
* M2 Spatial Context Regularization (Visium / Xenium)
* M3 Trajectory-aware Hierarchical Constraint (scVelo / CellRank)
* M4 Disease-state Conditional Markers
* M5 Counterfactual Explanation Layer
* M6 Patient-specific Active Conformal Calibration
* M7 Negative Marker Auto-discovery
* M8 Cell-Cell Communication Aware Marker Weighting
* M9 Cross-omics Validation (RNA + ATAC + Protein)
* M10 DB Drift Monitoring

See `manuscript/scONCO_Paper_Outline_v2_Advanced.docx` for full M1–M10 specification.
