# Changelog

All notable changes to msProteomiX will be documented in this file.

Format follows [Keep a Changelog](https://keepachangelog.com/).

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
