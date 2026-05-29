# msProteomiX

**Multi-Engine Mass Spectrometry Proteomics Data Processing**

A comprehensive R package for processing and visualizing mass spectrometry proteomics data from multiple search engines. Designed for both interactive use (RStudio Source mode) and programmatic workflows.

## Supported Search Engines

| Engine | Status | Detection | Input Files |
|--------|--------|-----------|-------------|
| **FragPipe** | ✅ Complete | Auto | `combined_protein.tsv`, `combined_peptide.tsv` |
| **Spectronaut** | ✅ Complete | Auto | `*_Report.tsv` (Run Pivot export) |
| **MaxQuant** | 🔲 Planned | Auto | `proteinGroups.txt` |
| **Proteome Discoverer** | 🔲 Planned | Auto | `*_Proteins.txt` |
| **DIA-NN** | 🔲 Planned | Auto | `report.pg_matrix.tsv` |

## Analysis Features

| Module | Scripts | Functions |
|--------|---------|-----------|
| Data Import | `01` | `read_ms_data()`, `interactive_grouping()` |
| ID Barplot | `02` | `plot_id_barplot()` (Protein + Peptide + PSM) |
| Venn / UpSet | `03` | `plot_venn()`, `plot_upset()` |
| PCA | `04` | `plot_pca()` (95% ellipses + labels) |
| CV Assessment | `05` | `plot_cv_boxplot()`, `plot_cv_violin()` |
| Correlation Heatmap | `06` | `plot_corr_heatmap()` (R² values) |
| Differential Analysis | `07` | `run_diff_analysis()`, `plot_volcano()` |
| Sequence Coverage | `08` | `calc_coverage()`, `plot_coverage()` |
| GO Enrichment | `09` | `run_go_enrichment()`, `plot_go_bubble()` |
| AP-MS Analysis | `20` | `normalize_bait()`, `run_anova_timecourse()`, `run_mfuzz_cluster()` |

## Quick Start

### Installation

```r
install.packages("devtools")
devtools::install_github("BayOmics/msProteomiX")
```

### Interactive Mode (Recommended for non-programmers)

```r
library(msProteomiX)

# 1. Create project structure
create_project("~/Desktop/MyProject")

# 2. Place search engine output in wkdir/
# 3. Open scripts/01 in RStudio → Source
# 4. Run scripts 02-09 as needed
```

### Programmatic Mode

```r
library(msProteomiX)

# Read data (auto-detects engine)
ms <- read_ms_data("path/to/results/")

# Set groups
groups <- set_groups(ms$sample_names, c("Ctrl", "Ctrl", "Treatment", "Treatment"))

# Differential analysis
result <- run_diff_analysis(ms, groups, method = "limma",
                            ref_group = "Ctrl", test_group = "Treatment")

# Volcano plot
plot_volcano(result, label_top = 20)
```

### Updating Scripts

After package update, sync your project scripts:

```r
library(msProteomiX)
update_scripts()  # auto-detects project, overwrites outdated scripts
```

## Documentation

- [Design Specification](vignettes/design_specification.md) — Architecture, data structures, extension guide
- [CHANGELOG](CHANGELOG.md) — Version history
- [Contributing](.github/CONTRIBUTING.md) — Development standards and pitfalls

## License

MIT
