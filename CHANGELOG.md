# Changelog

All notable changes to msProteomiX will be documented in this file.

Format follows [Keep a Changelog](https://keepachangelog.com/).

## [0.3.0] - 2026-06-17

### Added
- **Plasma QC Module** — blood contamination index evaluation
  - `calc_contamination_index()` with 4 marker panels (PLT/RBC/PBMC)
  - `plot_ci_barplot()`, `plot_ci_summary()`, `plot_ci_heatmap()`
  - `generate_plasma_report()` (EN/ZH, HTML + Chrome PDF)
  - Script `17_血浆质控评估.R`
- **QC Evaluation Panel** — 14 method evaluation plots in one call
  - `plot_qc_panel()` (peptide length, charge, missed cleavage, modifications, GRAVY, pI, M/Z, cumulative intensity, Cys/alkylation, RT, mass error, missing values, intensity boxplot, protein rank)
  - Scripts `10_QC方法评估.R`
- **Summary Table** — per-sample QC metrics CSV
  - `generate_summary_table()` (PSM, Peptide, ProteinGroups, 0-miss%, Cys%, Alk%)
  - Script `18_数据总结表.R`
- **Column Comparison Module** — chromatography column sensitivity comparison
  - `import_alphatims_results()`, `plot_lod_curve()`
  - Scripts `30–32`
- **Additional analysis modules** — scripts `11–16`
  - KEGG/Reactome enrichment, differential heatmap, PPI network, GSEA, marker boxplot, HTML report generation
- **Windows Deployment** — one-click deployment for non-technical users
  - `deploy/build_windows_bundle.sh` (Mac-side bundle builder)
  - `deploy/windows_install.R` (Windows-side installer with Tsinghua mirrors)
- **Automated Testing** — testthat framework with 9 smoke tests
- **Branch Strategy** — `release/v0.3.0` frozen for market team

### Fixed
- PCA crash on zero-variance proteins (`prcomp(scale.=TRUE)` error)
- Chrome path detection: empty Windows env vars producing phantom paths
- Windows installer: `readline()` → `utils::menu()` for non-interactive mode
- Windows installer: removed `grid` (base package) from install list

### Changed
- Repository cleanup: moved 8.5GB test data outside repo (13MB → clean)
- Untracked `project_wiki/` from git (kept locally)
- Cleaned up `.gitignore` (removed stale entries, added vignette artifacts)
- Updated README with all 23 scripts and deployment instructions

## [0.2.0] - 2026-05-29

### Added
- `update_scripts()` — one-command sync of project scripts with latest package version
- Auto version check in `setup_workdir()` — warns when scripts are outdated
- `R/imports.R` + `R/globals.R` — proper NAMESPACE imports and NSE variable declarations
- Median labels on CV boxplot and violin plots
- Correlation heatmap: R² values with white-to-red gradient (0.8 cutoff)
- Coverage calculation: per-group peptide filtering using Start/End columns from FragPipe
- Barplot: 3-panel chart (Protein + Peptide + PSM counts)
- PCA: 95% confidence ellipses + sample labels
- Venn: graceful skip for >4 groups (auto-fallback to UpSet)
- CHANGELOG.md, design_specification.md, CONTRIBUTING.md

### Changed
- **parsers.R split** into 5 focused files (parsers.R, parser_fragpipe.R, parser_spectronaut.R, parser_maxquant.R, parser_pd_diann.R)
- `create_project()` now always overwrites scripts (ensures latest version)
- Replaced `readr` with base R `read.delim()` globally (fixes vroom memory bug on Apple Silicon)
- All runtime strings use ASCII only (no emoji/em-dash) to prevent UTF-8 locale warnings
- `.msProteomiX_env` properly exported in NAMESPACE

### Fixed
- `%>%` not found: added `magrittr` + `rlang` to Imports
- `.msProteomiX_env` not found after clean R session
- UpSet PDF blank page: `onefile=FALSE` + `print(p)`
- Coverage inflated: was using global peptides instead of per-group
- Limma `topTable` row misalignment: `sort.by = "none"`
- FragPipe quantification column regex: anchored to `MaxLFQ Intensity$`
- Contrast selection: accepts group name or number
- `stat_ellipse` crash with <3 samples: conditional guard
- Spectronaut: tibble subsetting, trailing semicolons, double-log2 protection

### Removed
- `readr` from Imports (replaced by base R)
- Tracked artifacts from git (output/, test_output/, project_wiki/)

## [0.1.0] - 2026-05-28

### Added
- Initial release: FragPipe + Spectronaut parsers
- 10 analysis scripts (00-09) + AP-MS module (20)
- MsDataSet S3 class with unified data structure
- Interactive grouping wizard
- Full visualization suite: barplot, Venn/UpSet, PCA, CV, correlation, volcano
- Differential analysis (Limma + t-test)
- Sequence coverage analysis
- GO enrichment analysis
