# msProteomiX — 仓库结构优化计划

> 基于完整审计: R CMD check 结果 + 文件结构 + 代码扫描

## 审计发现总结

| 类别 | 发现 | 严重度 |
|------|------|--------|
| R CMD check | 1 error, 3 warnings, 4 notes | 🔴 |
| 缺少 imports | `aov`, `p.adjust`, `complete.cases` 等 stats/utils 函数未导入 | 🔴 |
| NSE 变量绑定 | 26 个 `no visible binding` (dplyr/ggplot2 的 NSE 列名) | 🟠 |
| 仓库杂物 | `output/`, `test_output/`, `test_validation.R`, `搜库软件输出结果demo/` 在仓库中 | 🟠 |
| parsers.R 过大 | 1160 行单文件，5 个引擎解析器混合 | 🟡 |
| 缺少文档 | 无 CHANGELOG, 无 design spec, README 过时 | 🟡 |
| 编码问题 | 残留 emoji/非 ASCII 字符在 runtime strings 中 | 🟡 |
| project_wiki 位置 | 内部文档放在仓库根目录，不属于包 | 🟢 |

---

## Phase 1: R CMD check 合规 (必须)

> [!IMPORTANT]
> 这是发布到 GitHub 的最低要求。当前 1 error + 3 warnings 无法通过 CI。

### [MODIFY] [imports.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/R/imports.R)
- 添加缺失的 stats/utils imports:
  ```r
  #' @importFrom stats aov complete.cases p.adjust rnorm setNames cor
  #' @importFrom utils flush.console
  ```

### [MODIFY] 多个 visualization_*.R + diff_analysis.R + apms.R
- 在每个使用 dplyr/ggplot2 NSE 的文件中添加 `utils::globalVariables()` 消除 NOTE:
  ```r
  utils::globalVariables(c("PC1", "PC2", "Group", "ShortLabel", ...))
  ```
- 涉及文件:
  - `visualization_barplot.R` → `mean_val`, `sd_val`
  - `visualization_pca.R` → `PC1`, `PC2`, `Group`, `ShortLabel`
  - `visualization_volcano.R` → `logFC`
  - `visualization_apms.R` → `Group`, `Value`, `DisplayLabel`, `Label`, etc.
  - `visualization_cv.R` → `CV`, `MedianCV`, `group`
  - `visualization_correlation.R` → `Sample1`, `Sample2`, `R2`
  - `coverage.R` → `Coverage`, `MedianCov`
  - `enrichment.R` → `p.adjust`, `Description`, `Count`

### 验证
```bash
Rscript -e 'devtools::check(".", args="--no-examples --no-tests --no-vignettes")'
# 目标: 0 errors, 0 warnings, ≤1 note
```

---

## Phase 2: 仓库卫生清理

> [!WARNING]
> 以下文件/目录不应出现在 Git 仓库中。清理后 `git push` 前请确认。

### 删除已追踪的杂物

| 路径 | 原因 | 操作 |
|------|------|------|
| `output/` | 运行时产物 | `git rm -r --cached` |
| `test_output/` | 测试产物 | `git rm -r --cached` |
| `test_validation.R` | 开发调试脚本 | `git rm --cached` |
| `搜库软件输出结果demo/` | 用户数据 | `git rm -r --cached` |
| `project_wiki/` | 内部文档，移入 `vignettes/` 或 `.github/` | `git rm -r --cached` |
| `.Rhistory` | RStudio 历史 | `git rm --cached` |
| `.DS_Store` (多处) | macOS 系统文件 | `git rm -r --cached` |

### [MODIFY] [.gitignore](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/.gitignore)
```diff
+ # Development artifacts
+ project_wiki/
+ .Rhistory
```
(确保以上所有路径已在 .gitignore 中)

### [MODIFY] [.Rbuildignore](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/.Rbuildignore)
```diff
+ ^project_wiki$
+ ^test_validation\.R$
+ ^test_output$
+ ^output$
+ ^\.github$
```

---

## Phase 3: parsers.R 拆分

> [!NOTE]
> `parsers.R` 有 1160 行，包含 5 个引擎的解析器 + 公共函数。拆分后更容易维护。

### 拆分方案

| 新文件 | 内容 | 预计行数 |
|--------|------|----------|
| `R/parsers.R` | `read_ms_data()`, `detect_engine()`, `new_MsDataSet()`, `.read_omics_file()` | ~120 |
| `R/parser_fragpipe.R` | `parse_fragpipe()`, `.extract_fp_quant_cols()` | ~250 |
| `R/parser_spectronaut.R` | `parse_spectronaut()`, `.detect_spectronaut_header()`, `.extract_sn_quant_cols()`, `.build_sn_label_name()` | ~350 |
| `R/parser_maxquant.R` | `parse_maxquant()`, `.extract_mq_quant_cols()`, `.build_mq_label_name()` | ~150 |
| `R/parser_pd_diann.R` | `parse_pd()`, `parse_diann()`, 相关辅助函数 | ~200 |

### 操作方式
- 纯文件拆分，不改变任何函数签名或逻辑
- 每个文件保留原有 roxygen 注释

---

## Phase 4: 文档补全

### [NEW] `CHANGELOG.md`
- 记录所有版本变更 (从 v0.1.0 开始)
- 格式: [Keep a Changelog](https://keepachangelog.com/)
- 内容摘自 git log

### [NEW] `vignettes/design_specification.Rmd`
- **架构设计**: MsDataSet S3 对象结构、.msProteomiX_env 跨脚本通信机制
- **引擎适配**: 各搜库引擎的列名约定、定量列检测策略
- **用户工作流**: create_project → setup_workdir → scripts/ 模式说明
- **扩展指南**: 如何添加新引擎 (parse_xxx 模板)

### [MODIFY] [README.md](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/README.md)
- 更新引擎状态表 (Spectronaut 已完成)
- 添加 badge: R CMD check, version
- 添加 `update_scripts()` 使用说明
- 添加功能截图/示例输出
- 英文为主 + 中文补充

### [NEW] `.github/CONTRIBUTING.md`
- 开发环境搭建
- 代码规范 (base R I/O, 无 readr, ASCII-only runtime strings)
- 已知陷阱清单 (从 project_wiki 迁移)

### [MODIFY] `project_wiki/` → 迁移后删除
- 内容整合到 `vignettes/design_specification.Rmd` 和 `.github/CONTRIBUTING.md`

---

## Phase 5: 编码硬化

### 全面 ASCII 扫描
- 扫描所有 `R/*.R` 文件中的非 ASCII 字符 (在 `message()`, `stop()`, `cat()`, 字符串字面量中)
- **注释中的中文保留** (不影响运行)
- **runtime 字符串全部 ASCII 化** (防止不同 locale 出 warning)

### [MODIFY] [00_安装指南.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/inst/scripts/00_安装指南.R)
- `readr` 仍在 CRAN 安装列表中但已从 Imports 移除 → 删除
- 添加 `magrittr`, `rlang` 到安装列表

---

## 优先级与工时估计

| Phase | 优先级 | 预计时间 | 依赖 |
|-------|--------|----------|------|
| Phase 1: R CMD check | 🔴 必须 | 15 min | 无 |
| Phase 2: 仓库卫生 | 🔴 必须 | 10 min | 无 |
| Phase 5: 编码硬化 | 🟠 推荐 | 10 min | 无 |
| Phase 3: parsers 拆分 | 🟡 建议 | 20 min | Phase 1 |
| Phase 4: 文档补全 | 🟡 建议 | 30 min | Phase 2, 3 |

> **总计: ~85 min**

## Open Questions

> [!IMPORTANT]
> 1. **parsers.R 拆分**: 是否同意拆分为 5 个文件？还是保持单文件？
> 2. **project_wiki**: 迁移到 `vignettes/` 还是 `.github/`？还是保留在根目录？
> 3. **design_specification**: 用 Rmd (可生成 HTML) 还是纯 Markdown？
> 4. **Phase 3+4 是否现在执行**：还是先确保 01-09 脚本全部跑通后再做？
