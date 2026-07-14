# msProteomiX

**Multi-Engine Mass Spectrometry Proteomics Data Processing**

A comprehensive R package for processing and visualizing mass spectrometry proteomics data from multiple search engines. Designed for both interactive use (RStudio Source mode) and programmatic workflows. All analysis runs **100% offline**.

## Supported Search Engines

| Engine | Status | Detection | Input Files |
|--------|--------|-----------|-------------|
| **FragPipe** | ✅ Complete | Auto | `combined_protein.tsv`, `combined_peptide.tsv` |
| **Spectronaut** | ✅ Complete | Auto | `*_Report.tsv` (Run Pivot export) |
| **MaxQuant** | 🔲 Planned | Auto | `proteinGroups.txt` |
| **Proteome Discoverer** | 🔲 Planned | Auto | `*_Proteins.txt` |
| **DIA-NN** | 🔲 Planned | Auto | `report.pg_matrix.tsv` |

## Analysis Modules

### Core Analysis (Scripts 01–09)

| Module | Script | Key Functions |
|--------|--------|---------------|
| Data Import & Grouping | `01` | `read_ms_data()`, `interactive_grouping()` |
| ID Barplot | `02` | `plot_id_barplot()` |
| Venn / UpSet | `03` | `plot_venn()`, `plot_upset()` |
| PCA | `04` | `plot_pca()` |
| CV Assessment | `05` | `plot_cv_boxplot()`, `plot_cv_violin()` |
| Correlation Heatmap | `06` | `plot_corr_heatmap()` |
| Differential Analysis | `07` | `run_diff_analysis()`, `plot_volcano()` |
| Sequence Coverage | `08` | `calc_coverage()`, `plot_coverage()` |
| GO Enrichment | `09` | `run_go_enrichment()`, `plot_go_bubble()` |

### Advanced Analysis (Scripts 10–18)

| Module | Script | Key Functions |
|--------|--------|---------------|
| QC Evaluation (14 plots) | `10` | `plot_qc_panel()` |
| KEGG & Reactome | `11` | `run_kegg_enrichment()`, `run_reactome_enrichment()` |
| Differential Heatmap | `12` | `plot_diff_heatmap()` |
| PPI Network | `13` | `run_ppi_network()` |
| GSEA | `14` | `run_gsea_analysis()` |
| Marker Expression | `15` | `plot_marker_boxplot()` |
| Report Generation | `16` | `generate_report()` |
| Plasma QC | `17` | `calc_contamination_index()`, `generate_plasma_report()` |
| Summary Table | `18` | `generate_summary_table()` |

### Specialized Modules

| Module | Script | Key Functions |
|--------|--------|---------------|
| AP-MS Time Course | `20` | `normalize_bait()`, `run_anova_timecourse()` |
| Column Comparison | `30–32` | `import_alphatims_results()`, `plot_lod_curve()` |

## Quick Start

### Installation

```r
# One-liner install from GitHub (recommended)
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("BayOmics/msProteomiX")
```

> **Note**: Some optional features (GO/KEGG enrichment, GSEA, plasma QC) require
> Bioconductor packages. Install them when needed:
> ```r
> if (!requireNamespace("BiocManager")) install.packages("BiocManager")
> BiocManager::install(c("limma", "clusterProfiler", "org.Hs.eg.db"))
> ```

<details>
<summary>Alternative: install from source tarball</summary>

```r
install.packages("msProteomiX_0.3.0.tar.gz", repos = NULL, type = "source")
```
</details>

### Interactive Mode (Recommended)

```r
library(msProteomiX)

# 1. Create project structure
create_project("~/Desktop/MyProject")

# 2. Place search engine output in wkdir/
# 3. Open scripts/01 in RStudio → click Source
# 4. Run scripts 02–18 as needed
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

## Windows Deployment

A one-click deployment bundle is available for non-technical users:

```bash
# On Mac (developer): build the bundle
bash deploy/build_windows_bundle.sh

# On Windows (user): open windows_install.R in RStudio → Source
```

See `deploy/` directory for details.

## Repository Structure

```
msProteomiX/
├── R/                    # 33 source files
├── inst/
│   ├── extdata/          # Built-in data (marker panels)
│   └── scripts/          # 23 user scripts (00–18, 20, 30–32)
├── man/                  # roxygen2 documentation
├── tests/testthat/       # Automated tests
├── deploy/               # Windows deployment tools
├── vignettes/            # Design specification
├── DESCRIPTION           # Package metadata
└── NAMESPACE             # Auto-generated exports
```

## Documentation

- [Design Specification](vignettes/design_specification.md) — Architecture, data structures, extension guide
- [CHANGELOG](CHANGELOG.md) — Version history
- [Contributing](.github/CONTRIBUTING.md) — Development standards and pitfalls

## License

MIT
