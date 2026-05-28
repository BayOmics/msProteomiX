# msProteomiX <img src="man/figures/logo.png" align="right" height="138" />

**Multi-Engine Mass Spectrometry Proteomics Data Processing**

一个面向质谱蛋白组学的 R 包，支持多种搜库引擎的统一数据处理和可视化。

## 支持的搜库引擎

| Engine | Status | Input Files |
|--------|--------|-------------|
| **FragPipe** | ✅ Supported | `combined_protein.tsv`, `combined_peptide.tsv` |
| **MaxQuant** | 🔜 Coming | `proteinGroups.txt` |
| **Proteome Discoverer** | 🔜 Coming | `_Proteins.txt` |
| **DIA-NN** | 🔜 Coming | `report.pg_matrix.tsv` |

## 分析功能

- 📊 定性分析柱状图 (Protein/Peptide/PSM counts)
- 🔵 Venn/UpSet 蛋白重叠分析
- 📈 PCA 主成分分析
- 📉 CV 变异系数评估
- 🗺️ 相关性热图
- 🌋 差异分析 & 火山图 (Limma / t-test)
- 📏 序列覆盖度分析
- 🧬 GO 富集分析
- 🔬 AP-MS 亲和纯化分析 (Coming Soon)

## 快速上手

### 安装

```r
# 在 RStudio 中运行
install.packages("devtools")
devtools::install_github("BayOmics/msProteomiX")
```

### 使用 (RStudio 交互模式)

1. 打开 `inst/scripts/00_安装指南.R`，点击 Source
2. 打开 `inst/scripts/01_数据导入与分组.R`，点击 Source
3. 按需运行 02-09 分析脚本

### 使用 (编程模式)

```r
library(msProteomiX)

# 读取数据
ms <- read_ms_data("path/to/fragpipe_output/", engine = "auto")

# 设置分组
groups <- set_groups(
  get_sample_names(ms),
  c("Control", "Control", "Treatment", "Treatment")
)

# 差异分析
result <- run_diff_analysis(ms, groups, method = "limma",
                             ref_group = "Control", test_group = "Treatment")

# 火山图
plot_volcano(result, label_top = 20)
```

## License

MIT
