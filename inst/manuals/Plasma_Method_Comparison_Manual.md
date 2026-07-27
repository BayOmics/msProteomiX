# BayOmics msProteomiX -- Plasma Method Comparison User Manual

> Target audience: BayOmics marketing team
> Last updated: 2026-07

---

## 1. Overview

This script compares different bead-based plasma proteomics sample preparation
methods. It generates a **comprehensive evaluation report** with 14+ quality
metrics, exported as a self-contained HTML file that can be printed to PDF.

**What you get**: A professional comparison report showing which method identifies
more proteins, which has better quantification precision, depletion efficiency, etc.

---

## 2. Prerequisites

### 2.1 Install RStudio

Download and install RStudio Desktop from: https://posit.co/download/rstudio-desktop/

### 2.2 Install msProteomiX

Open RStudio, and run this in the Console (bottom-left panel):

```r
# Only needed once
install.packages("devtools")
devtools::install_github("BayOmics/msProteomiX")
install.packages("base64enc")   # For HTML report generation
```

### 2.3 Create a Project

In the RStudio Console, run:

```r
library(msProteomiX)
create_project("~/Desktop/PlasmaComparison")
```

This creates a project folder with all scripts and a `manuals/` folder.

---

## 3. Data Preparation

### 3.1 What to Request from the Bioinformatics Team

**For each condition** (e.g., "BayOmics_3uL", "Nanomics_5uL"), request:

#### Spectronaut results (3 files per condition):

| File | Purpose | Example |
|------|---------|---------|
| `*_Report_ProtienGroup_MSC (Pivot).tsv` | Protein quantification | Required |
| `*_IdentificationsOverview.tsv` | Per-sample ID counts | Required |
| `*_Peptide_Report.tsv` | Peptide-level analysis | Optional but recommended |

#### DIA-NN results (1 file per condition):

| File | Purpose |
|------|---------|
| `*_report.pg_matrix.tsv` | Protein quantification |

### 3.2 Organize Files

Place each condition's files in its own subfolder under `wkdir/`:

```
PlasmaComparison/
  wkdir/
    Spectronaut/
      Condition_A/        <-- Put all 3 files here
      Condition_B/        <-- Put all 3 files here
    DIA_NN/
      pg_matrix_A.tsv     <-- One file per condition
      pg_matrix_B.tsv
  scripts/
    19_...R               <-- The script you will edit
  output/                 <-- Results will appear here
```

---

## 4. Edit the Configuration

Open `scripts/19_血浆前处理方案比较.R` in RStudio.

You only need to edit the **CONFIG** section (lines 22-42). Everything below
the line `# EXECUTION -- Do NOT edit below this line` should not be touched.

### 4.1 Set Conditions

Edit the `conditions` list. Each condition needs:
- `name`: A short label (use English, no spaces)
- `path`: Relative path to the data folder (starting with `./wkdir/`)
- `engine`: `"spectronaut"` or `"diann"`

**Example** (comparing 2 methods, each with 2 volume conditions):

```r
conditions <- list(
  list(name = "BayOmics_3",   path = "./wkdir/Spectronaut/BM_3",   engine = "spectronaut"),
  list(name = "BayOmics_5",   path = "./wkdir/Spectronaut/BM_5",   engine = "spectronaut"),
  list(name = "Competitor_3", path = "./wkdir/Spectronaut/NM_3",   engine = "spectronaut"),
  list(name = "Competitor_5", path = "./wkdir/Spectronaut/NM_5",   engine = "spectronaut")
)
```

### 4.2 Set Groups

Edit `compare_groups` to define which conditions belong to which group.
Groups determine the color coding in charts.

```r
compare_groups <- list(
  BayOmics   = c("BayOmics_3", "BayOmics_5"),
  Competitor = c("Competitor_3", "Competitor_5")
)
```

### 4.3 Set Project Name

```r
project_name <- "Customer_PlasmaG"   # Appears in report header
```

---

## 5. Run the Script

1. In RStudio, open the script file
2. Press **Ctrl+A** (select all)
3. Press **Ctrl+Enter** (run)
4. Wait for the message: `All evaluations complete!`

This typically takes 1-3 minutes depending on the number of conditions.

---

## 6. Find Your Results

After the script finishes, check the `output/` directory:

```
output/
  BayOmics_3_evaluation/      <-- Per-condition CI + QC results
  Nanomics_3_evaluation/
  ...
  comparison/                  <-- Cross-condition comparison
    10_comparison_report_en.html   <-- English report (main deliverable)
    10_comparison_report_en.pdf    <-- English PDF (if Chrome installed)
    10_comparison_report_zh.html   <-- Chinese report
    10_comparison_report_zh.pdf    <-- Chinese PDF
    10_pg_comparison.pdf           <-- Individual charts (for PPT use)
    11_precursor.pdf
    ...
```

### 6.1 The Main Report

Open `10_comparison_report_en.html` (or `_zh.html`) in your web browser.
This is a **self-contained** file with all charts embedded. You can:

- **Email it directly** to clients
- **Print to PDF**: File > Print > Save as PDF (if auto-PDF failed)

### 6.2 Individual Charts

All comparison charts are also saved as individual PDF files in
`output/comparison/` for use in PowerPoint presentations.

---

## 7. Troubleshooting

| Error | Solution |
|-------|----------|
| `Error: package 'msProteomiX' not found` | Run `devtools::install_github("BayOmics/msProteomiX")` |
| `Error: package 'base64enc' not found` | Run `install.packages("base64enc")` |
| `Error: cannot open file '...'` | Check that the file path is correct and the file exists |
| `Chrome not found` | Install Google Chrome, or open the HTML file in a browser and use Print > Save as PDF |
| Chinese characters in path cause errors | Move data to a path with only English characters |

---

## 8. Appendix: What Each Metric Means

| Metric | What It Measures | Good Result |
|--------|-----------------|-------------|
| **Protein Groups** | How many proteins were identified | More is better |
| **Precursors** | How many peptide precursors identified | More is better |
| **Data Completeness** | % proteins detected in all replicates | Higher is better |
| **Quant CV** | Reproducibility of protein quantities | Lower is better (<20% ideal) |
| **Correlation** | Agreement between replicates | Higher is better (>0.95 ideal) |
| **PCA** | Sample clustering pattern | Replicates should cluster together |
| **CI (Contamination Index)** | Blood cell contamination level | Pass = clean sample |
| **Top-N Fraction** | How much signal comes from top proteins | Lower = better depletion |
| **High-Abundance Rank** | Where ALB/IgG rank in intensity | Higher rank = better depletion |
| **Missed Cleavage** | Digestion efficiency (0-miss rate) | Higher 0-miss = better digestion |
| **Peptide Length** | Length of identified peptides | 7-25 aa is typical |
| **Dynamic Range** | Intensity span (P1-P99) | Wider = better sensitivity |
| **Rank Curve** | Protein abundance distribution | Smooth curve is expected |

---

*Manual by BayOmics Bioinformatics Team. For support, contact the bioinformatics team.*
