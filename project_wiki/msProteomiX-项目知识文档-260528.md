# msProteomiX — 项目知识文档

> 供新对话/新 workspace 快速恢复上下文

## 项目基本信息

- **路径**: `/Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/`
- **类型**: 独立 R 包 (非 Shiny app)
- **用户画像**: 不懂代码的客户/市场推广人员，通过 RStudio "Source" 按钮运行
- **来源**: 从 `/Users/wangxi/Documents/WPSCloud/DT/Fragpipe_Data_Processing_marketing_v20260115/` 的 23 个散乱脚本重构而来
- **Git**: 4 commits on `main`, 未推送到 remote

## 架构决策

| 决策 | 选择 | 原因 |
|------|------|------|
| 包形式 | 独立 R 包 | 用户要求"保持独立R包" |
| 用户交互 | `inst/scripts/` + Source 按钮 | 用户要求"完全不懂代码的人可用" |
| 数据持久化 | `.msProteomiX_env` (package env) | Source 模式下跨脚本传递数据 |
| 分组机制 | `interactive_grouping()` 统一函数 | 原项目 6 处重复代码 |
| 差异分析 | 统一 `run_diff_analysis()` | 原项目 6 个火山图脚本变体 |
| 定量列解析 | 后缀匹配 `MaxLFQ Intensity$` | FragPipe 列名 = "SampleName MaxLFQ Intensity" |

## 模块清单 (15 个 R 文件)

### 核心
- [utils.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/R/utils.R) — `%||%`, `.msProteomiX_env`, `MsDataSet` S3 class, `save_plot_and_data()`, `mspx_colors()`
- [parsers.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/R/parsers.R) — `read_ms_data()`, `detect_engine()`, `parse_fragpipe()` (MQ/PD/DIA-NN 占位)
- [grouping.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/R/grouping.R) — `interactive_grouping()`, `set_groups()`, `check_prerequisites()`
- [zzz.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/R/zzz.R) — `.onLoad` / `.onAttach`

### QC
- [qc_identification.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/R/qc_identification.R) — `count_ids()`, `calc_zero_miss()`, `calc_cys_percent()`, `calc_alk_efficiency()`
- [qc_evaluation.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/R/qc_evaluation.R) — `calc_pI()`, `calc_gravy()`, `calc_missed_cleavage()`

### 可视化
- [visualization_barplot.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/R/visualization_barplot.R) — `plot_id_barplot()` (TODO: target 参数未实现)
- [visualization_venn.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/R/visualization_venn.R) — `plot_venn()`, `plot_upset()`
- [visualization_pca.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/R/visualization_pca.R) — `plot_pca()` (含 <3 样本组椭圆保护)
- [visualization_cv.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/R/visualization_cv.R) — `plot_cv_boxplot()`, `plot_cv_violin()`
- [visualization_correlation.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/R/visualization_correlation.R) — `plot_corr_heatmap()`
- [visualization_volcano.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/R/visualization_volcano.R) — `plot_volcano()` (用 `%||%` from utils.R)

### 分析
- [diff_analysis.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/R/diff_analysis.R) — `run_diff_analysis()` (Limma + t-test, **topTable sort.by="none" 已修复**)
- [coverage.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/R/coverage.R) — `calc_coverage()`, `plot_coverage()` (需 Biostrings + FASTA)
- [enrichment.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/R/enrichment.R) — `run_go_enrichment()`, `plot_go_bubble()` (需 clusterProfiler)

## 客户脚本 (inst/scripts/)

```
00_安装指南.R           → 一键安装所有依赖
01_数据导入与分组.R     → read_ms_data + interactive_grouping (必须首先运行)
02_定性分析_柱状图.R    → plot_id_barplot
03_Venn维恩图.R         → plot_venn + plot_upset
04_PCA主成分分析.R      → plot_pca
05_CV变异系数.R         → plot_cv_boxplot + plot_cv_violin
06_相关性热图.R         → plot_corr_heatmap
07_差异分析_火山图.R    → run_diff_analysis + plot_volcano
08_序列覆盖度.R         → calc_coverage + plot_coverage
09_GO富集分析.R         → run_go_enrichment + plot_go_bubble
```

## 验证状态

**3 个真实 FragPipe 数据集，14 项测试全部通过:**

| 数据集 | 蛋白 | 样本 | 状态 |
|--------|------|------|------|
| FISAP HeLa (202512001) | 453 | 14 | ✅ 14/14 |
| 17min Gradient | 321 | 6 | ✅ ALL PASS |
| HTProX (R96_8) | 5,317 | 24 | ✅ ALL PASS |

数据位置: `/Users/wangxi/Documents/WPSCloud/DT/Fragpipe_Data_Processing_marketing_v20260115/data/`

## 深度审核 — 已修复的 Bug

| # | 严重度 | 问题 | 修复 |
|---|--------|------|------|
| C1 | 🔴 | 定量列正则匹配前缀而非后缀 | `MaxLFQ Intensity$` |
| C2 | 🔴 | `^sp\|^tr\|` 无效正则 | `^(sp|tr)\|` |
| C3+C4 | 🔴 | Limma topTable 默认排序导致行错位 | `sort.by = "none"` |
| H1 | 🟠 | `[Cc]` 应为 `C` | 已修 |
| H3 | 🟠 | zzz.R 影子环境 | 删除重复定义 |
| H4 | 🟠 | `ylim()` 裁剪数据 | `coord_cartesian()` |
| M2 | 🟡 | coverage str_locate_all regex风险 | `fixed()` |
| M3 | 🟡 | stat_ellipse <3样本崩溃 | 条件保护 |
| L1 | 🟢 | `%||%` 是 R 4.4.0+ 才有 | 恢复自定义 |

## 未完成工作 (Phase 2+3)

### Phase 2: 扩展引擎 + AP-MS
- [ ] `parse_maxquant()` — 解析 `proteinGroups.txt`
- [ ] `parse_pd()` — 解析 Proteome Discoverer `_Proteins.txt`
- [ ] `parse_diann()` — 解析 `report.pg_matrix.tsv`
- [ ] AP-MS 模块 — SAINTexpress/CRAPome 集成
- [ ] `devtools::document()` 生成 man/ 页面
- [ ] `R CMD check` 通过

### Phase 3: 磷酸化蛋白组学
- [ ] 磷酸化位点解析
- [ ] Motif 分析
- [ ] Kinase 活性推断

## Git 历史

```
5f41785 chore: add test artifacts to .gitignore
d381d8a fix: restore %||% operator
82db997 fix: 14 bugs from deep audit
b53dd05 feat: initial msProteomiX R package scaffold
```

## 关键陷阱提醒

1. **`.msProteomiX_env` 只在 `utils.R` 中定义**，`zzz.R` 不再重新创建
2. **Limma `topTable` 必须用 `sort.by = "none"`**，否则行错位
3. **FragPipe 列名是 `SampleName MaxLFQ Intensity`**（样本名在前，后缀固定），正则必须锚定 `$`
4. **R 4.1.0 没有 base `%||%`**，需要自定义（4.4.0 才有）
5. **`stat_ellipse` 要求每组 ≥3 样本**，必须检查
