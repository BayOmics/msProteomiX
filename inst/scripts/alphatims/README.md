# AlphaTims Scripts for msProteomiX

Python scripts for processing timsTOF Pro raw data using [AlphaTims](https://github.com/MannLabs/alphatims).

## Prerequisites

```bash
pip install alphatims
```

## Scripts

### extract_precursors.py

Extract MS1 precursors from timsTOF `.d` / `.hdf` files and compute
**Contamination Ratio (CR)** per RT slice.

#### Quick Start

```bash
# Single file (Windows — reads .d directly with Bruker library)
python extract_precursors.py "C:\data\sample.d" --output_dir ./cr_output

# Single file (Mac/Linux — use .hdf for accurate mz/mobility)
python extract_precursors.py sample.hdf --output_dir ./cr_output

# Batch: process all .d/.hdf in a directory
python extract_precursors.py /path/to/data/ --output_dir ./cr_output

# Export per-slice CSVs for heatmap visualization (at RT=32 min)
python extract_precursors.py sample.d --export_slices --heatmap_rts 32
```

#### Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `input` | (required) | `.d` folder, `.hdf` file, or directory |
| `--output_dir` | `./cr_output` | Output directory |
| `--n_slices` | `10` | Number of equally-spaced RT slices |
| `--line_points` | `350,0.8,950,1.3` | CR dividing line: `mz1,im1,mz2,im2` |
| `--export_slices` | `False` | Export per-slice CSV files |
| `--heatmap_rts` | `None` | RT values for heatmap export (e.g., `32,48`) |

#### Output

```
cr_output/
├── {sample}_cr_summary.csv         # Summary: rt_min, tot_inten, contamination_ratio
├── {sample}/                       # Per-slice CSVs (with --export_slices)
│   ├── {sample}_slice_rt00.0.csv
│   ├── {sample}_slice_rt08.0.csv
│   └── ...
└── all_cr_summary.csv              # Combined summary (batch mode)
```

## .d → .hdf Conversion (for macOS users)

> **Important**: On macOS, alphatims cannot use Bruker's native library.
> The `mz_values` and `mobility_values` are **estimated** from metadata,
> with errors potentially up to 6 Th. For accurate values, convert `.d`
> to `.hdf` on **Windows** first.

### On Windows (PowerShell)

```powershell
# Single file
alphatims export hdf "C:\path\to\sample.d"
# -> Generates C:\path\to\sample.hdf

# Batch conversion
Get-ChildItem -Path "C:\data\" -Filter "*.d" -Directory | ForEach-Object {
    alphatims export hdf $_.FullName
}
```

### On Windows (Python)

```python
import alphatims.bruker

data = alphatims.bruker.TimsTOF("sample.d")
data.save_as_hdf(directory="./", file_name="sample.hdf", overwrite=True)
```

Transfer the `.hdf` file to Mac, then use it directly:

```bash
python extract_precursors.py sample.hdf --output_dir ./cr_output
```

## Downstream R Analysis

After running the Python script, use the R functions in msProteomiX:

```r
library(msProteomiX)

# Import CR results
cr_data <- import_cr_results("cr_output", group_map = list(
  "HeLa"          = c("Hela"),
  "W/o Ext. Wash" = c("DDM"),
  "W/ Ext. Wash"  = c("FFHE")
))

# Plot CR gradient (Fig. 2f)
plot_cr_gradient(cr_data)

# Plot IM-m/z heatmap (Fig. 2g)
plot_cr_heatmap("cr_output/sample_name/sample_slice_rt32.0.csv")
```
