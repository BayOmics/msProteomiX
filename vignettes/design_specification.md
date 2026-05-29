# msProteomiX — Design Specification

> Architecture, data structures, and extension guide

## 1. Overview

msProteomiX is an R package for multi-engine mass spectrometry proteomics data processing. It provides:

- **Unified parsing** of 5 search engines into a common data structure
- **Interactive analysis scripts** for users with no coding experience
- **Modular visualization** and statistical analysis functions

### Target Users

| User Type | Interface | Description |
|-----------|-----------|-------------|
| Lab technician / Marketing | `inst/scripts/` + RStudio Source | No coding required |
| Bioinformatician | `library(msProteomiX)` API | Full programmatic access |

---

## 2. Architecture

### 2.1 Core Data Structure: `MsDataSet`

An S3 class (defined in `R/utils.R`) that standardizes output from all search engines:

```
MsDataSet (list)
├── proteins        # numeric matrix [n_proteins × n_samples] — quantification values
├── protein_info    # data.frame — protein metadata (ID, Gene, Description, etc.)
├── sample_names    # character vector — clean sample names
├── peptides        # numeric matrix [n_peptides × n_samples] (optional)
├── peptide_info    # data.frame — peptide metadata (Sequence, Start, End, etc.)
├── ions            # numeric matrix (optional, FragPipe only)
├── ion_info        # data.frame (optional)
├── psm_list        # list of data.frames (optional, FragPipe only)
├── engine          # character — "fragpipe", "spectronaut", etc.
└── is_log2         # logical — TRUE if data is already log2-transformed
```

**Key invariant**: `ncol(proteins) == length(sample_names)`, column order matches.

### 2.2 Cross-Script Communication: `.msProteomiX_env`

A package-level environment exported for user scripts in RStudio Source mode:

```r
# Script 01 stores data:
.msProteomiX_env$ms_data    <- ms
.msProteomiX_env$group_info <- gi

# Script 02+ retrieves data:
ms <- .msProteomiX_env$ms_data
gi <- .msProteomiX_env$group_info
```

**Why not global variables?** Package environments persist across `source()` calls within the same R session, while avoiding global namespace pollution.

### 2.3 Project Structure

```
user_project/
├── scripts/     ← copied from inst/scripts/ via create_project()
│   ├── 00_安装指南.R
│   ├── 01_数据导入与分组.R
│   ├── ...
│   └── 09_GO富集分析.R
├── wkdir/       ← user places search engine output here
│   ├── combined_protein.tsv    (FragPipe)
│   ├── combined_peptide.tsv    (FragPipe)
│   └── ...
└── output/      ← analysis results auto-saved here
    ├── Project_IDBarplot.pdf
    ├── Project_PCA.pdf
    └── ...
```

### 2.4 Script Versioning

`setup_workdir()` auto-checks script MD5 against installed package version.
`update_scripts()` overwrites outdated scripts.

---

## 3. Search Engine Adapters

### 3.1 Supported Engines

| Engine | File | Status | Detection Pattern | Quantification Column |
|--------|------|--------|-------------------|-----------------------|
| FragPipe | `parser_fragpipe.R` | ✅ Complete | `combined_protein.tsv` | `<sample> MaxLFQ Intensity` |
| Spectronaut | `parser_spectronaut.R` | ✅ Complete | `*_Report.tsv` + PG. headers | Format A: clean / Format B: `[N] file.d.PG.Quantity` |
| MaxQuant | `parser_maxquant.R` | 🔲 Skeleton | `proteinGroups.txt` | `LFQ intensity <sample>` |
| Proteome Discoverer | `parser_pd_diann.R` | 🔲 Skeleton | `*_Proteins.txt` | `Abundance: <sample>` |
| DIA-NN | `parser_pd_diann.R` | 🔲 Skeleton | `report.pg_matrix.tsv` | Each column = one sample |

### 3.2 Adding a New Engine

1. Create `R/parser_<engine>.R`
2. Implement `parse_<engine>(path)` returning `new_MsDataSet(...)`
3. Add detection logic in `detect_engine()` (in `R/parsers.R`)
4. Add case in `read_ms_data()` switch statement
5. Export in roxygen: `#' @export`
6. Run `devtools::document()`

**Template**:

```r
#' @export
parse_newengine <- function(path) {
  if (.is_directory(path)) {
    path <- .find_file(path, "expected_file.tsv")
  }
  
  df <- .read_omics_file(path)
  
  # Extract quantification columns
  quant_cols <- grep("your_pattern", colnames(df), value = TRUE)
  quant_mat  <- as.matrix(df[, quant_cols])
  
  # Clean sample names
  sample_names <- sub("your_prefix", "", quant_cols)
  colnames(quant_mat) <- sample_names
  
  # Build protein info
  meta_cols <- setdiff(colnames(df), quant_cols)
  protein_info <- df[, meta_cols, drop = FALSE]
  protein_info <- .build_label_name(protein_info)
  
  new_MsDataSet(
    proteins     = quant_mat,
    protein_info = protein_info,
    sample_names = sample_names,
    engine       = "newengine",
    is_log2      = FALSE
  )
}
```

---

## 4. Module Inventory

### R/ Directory (22 files)

| File | Lines | Purpose |
|------|-------|---------|
| `parsers.R` | 192 | `read_ms_data()`, `detect_engine()`, shared helpers |
| `parser_fragpipe.R` | 188 | FragPipe parser |
| `parser_spectronaut.R` | 246 | Spectronaut parser (dual format) |
| `parser_maxquant.R` | 209 | MaxQuant parser (skeleton) |
| `parser_pd_diann.R` | 373 | PD + DIA-NN parsers (skeleton) |
| `utils.R` | 328 | `MsDataSet`, `create_project()`, `update_scripts()`, `setup_workdir()` |
| `imports.R` | 10 | Package-level imports (magrittr, rlang, stats, utils) |
| `globals.R` | 25 | `globalVariables()` for NSE columns |
| `zzz.R` | 22 | `.onLoad` / `.onAttach` hooks |
| `grouping.R` | 236 | Interactive grouping wizard |
| `diff_analysis.R` | 222 | Limma / t-test differential analysis |
| `coverage.R` | 187 | Sequence coverage (Biostrings + FASTA) |
| `enrichment.R` | 133 | GO enrichment (clusterProfiler) |
| `apms.R` | 487 | AP-MS: bait normalization, ANOVA, Mfuzz clustering |
| `qc_identification.R` | 146 | ID counts, zero-miss cleavage, alkylation efficiency |
| `qc_evaluation.R` | 83 | pI, GRAVY, missed cleavage distributions |
| `visualization_barplot.R` | 95 | ID barplot (3-panel) |
| `visualization_venn.R` | 135 | Venn + UpSet plots |
| `visualization_pca.R` | 119 | PCA with ellipses |
| `visualization_cv.R` | 125 | CV boxplot + violin |
| `visualization_correlation.R` | 94 | R² correlation heatmap |
| `visualization_volcano.R` | 101 | Volcano plot |
| `visualization_apms.R` | 332 | AP-MS dotplot + heatmap + Mfuzz plots |

---

## 5. Known Pitfalls

### Critical

| # | Issue | Mitigation |
|---|-------|------------|
| 1 | **readr/vroom memory bug** on Apple Silicon RStudio | Use `base::read.delim()` exclusively |
| 2 | **Limma `topTable` default sort** misaligns rows | Always use `sort.by = "none"` |
| 3 | **`.msProteomiX_env` must be exported** | `@export` tag in utils.R |
| 4 | **`stat_ellipse` crashes** with <3 samples per group | Conditional guard in `plot_pca()` |
| 5 | **FragPipe regex** must anchor `MaxLFQ Intensity$` | Suffix match, not prefix |

### Engine-Specific

| Engine | Issue | Fix |
|--------|-------|-----|
| FragPipe | Old versions have `Intensity=0` | Fallback to Spectral Count |
| FragPipe | `combined_protein.tsv` regex too broad | Exact basename match |
| Spectronaut | `protein_info` returned as tibble | Force `as.data.frame()` in `new_MsDataSet()` |
| Spectronaut | Gene column trailing semicolons | `str_remove(";+$")` |
| Spectronaut | `PG.Log2Quantity` → double log2 | `is_log2 = TRUE` flag |
| Spectronaut | NaN values (not NA) | `is.nan()` check before conversion |

### Encoding

| Issue | Fix |
|-------|-----|
| Emoji in `message()` → UTF-8 warning | Use ASCII symbols only (`>>>`, not `ℹ️`) |
| Em-dash `—` in strings | Replace with ASCII `-` |
| Non-ASCII in comments | Allowed (R ignores comments at runtime) |

---

## 6. Dependencies

### Required (Imports)

| Package | Purpose |
|---------|---------|
| `dplyr` | Data manipulation |
| `magrittr` | Pipe operator `%>%` |
| `rlang` | `.data` pronoun for NSE |
| `stringr` | String operations |
| `ggplot2` | Plotting |
| `ggrepel` | Non-overlapping labels |
| `reshape2` | `melt()` / `dcast()` |
| `scales` | Axis formatting |
| `tibble` | Data structure utilities |

### Optional (Suggests)

| Package | Purpose |
|---------|---------|
| `limma` | Differential analysis |
| `Biostrings` | FASTA parsing for coverage |
| `clusterProfiler` | GO enrichment |
| `org.Hs.eg.db` | Human gene annotation |
| `ggvenn` | Venn diagrams |
| `UpSetR` | UpSet plots |
| `Mfuzz` | Fuzzy clustering (AP-MS) |
