# msProteomiX 使用说明

## 快速上手 — 3 步完成分析

### 步骤 1: 安装（仅需一次）

在 RStudio 中打开 `00_安装指南.R`，点击右上角 **Source** 按钮，等待安装完成。

### 步骤 2: 创建分析项目

安装完成后，在 RStudio 中新建一个 R Project（或打开已有的），然后在 Console 中输入：

```r
library(msProteomiX)
create_project(getwd())
```

运行后会在你的项目中自动创建以下目录结构：

```
你的项目/
├── scripts/   ← 分析脚本（自动从包中复制）
├── wkdir/     ← 搜库结果放这里
└── output/    ← 分析结果自动保存在这里
```

### 步骤 3: 运行分析

1. **放入数据**：将搜库软件的结果文件放入 `wkdir/` 目录
   - FragPipe — 复制整个输出文件夹
   - Spectronaut — 复制 `*_Report.tsv`（Run Pivot 导出）
   - MaxQuant — 复制 `txt/` 文件夹（含 `proteinGroups.txt`）
   - PD — 复制 `*_Proteins.txt`
   - DIA-NN — 复制 `report.pg_matrix.tsv`

2. **运行脚本**：在 RStudio 中打开 `scripts/01_数据导入与分组.R` → 点击 Source → 按提示完成分组

3. **后续分析**：打开 `scripts/02` ~ `09` 脚本，Source 运行（01 必须首先运行，02-09 顺序任意）

4. **查看结果**：`output/` 文件夹中包含 PDF 图表和 CSV 源数据

---

## 脚本列表

| 脚本 | 功能 |
|------|------|
| `00_安装指南.R` | 一键安装所有依赖（仅需运行一次） |
| `01_数据导入与分组.R` | 读取数据 + 交互式分组（**必须首先运行**） |
| `02_定性分析_柱状图.R` | 蛋白/肽段/PSM 鉴定数量柱状图 |
| `03_Venn维恩图.R` | 蛋白重叠分析 Venn / UpSet 图 |
| `04_PCA主成分分析.R` | PCA 主成分分析 |
| `05_CV变异系数.R` | 定量重复性 CV 箱线图 / 小提琴图 |
| `06_相关性热图.R` | 样品相关性热图 |
| `07_差异分析_火山图.R` | Limma / t-test 差异分析 + 火山图 |
| `08_序列覆盖度.R` | 蛋白序列覆盖度分析（需要 FASTA） |
| `09_GO富集分析.R` | GO 细胞组分富集气泡图 |
| `20_AP-MS_时序分析.R` | AP-MS 亲和纯化时序分析（Bait 归一化 + ANOVA + DotPlot） |

---

## 常见问题

**Q: "Error: 请先运行 01\_数据导入与分组.R"**
> 说明你还没有运行第一步。请先打开并 Source 01 脚本。

**Q: "Error: 未找到 combined\_protein.tsv"**
> 请确认你的搜库结果已放入 `wkdir/` 文件夹中。

**Q: 安装 Bioconductor 包失败**
> 尝试手动安装：
> ```r
> if (!require("BiocManager")) install.packages("BiocManager")
> BiocManager::install("limma")
> ```

**Q: 想修改分组怎么办？**
> 重新运行 `01_数据导入与分组.R`，选择"重新分组"即可。

**Q: 火山图参数怎么调？**
> 打开 `07_差异分析_火山图.R`，修改顶部"用户设置"区域的参数。

**Q: "vector memory exhausted" 内存不足**
> 在 Console 中运行以下命令，然后重启 RStudio：
> ```r
> writeLines("R_MAX_VSIZE=32Gb", "~/.Renviron")
> ```

---

*技术支持：如遇到问题，请联系 BayOmics 技术团队。*
