# msProteomiX — 仓库结构优化计划 (v2 — 已确认)

> 所有 5 个 Phase 已批准，立即执行

## Phase 1: R CMD check 合规 🔴

### [MODIFY] [imports.R](file:///Users/wangxi/Developer/BayOmics/GitHub/msProteomiX/R/imports.R)
- 添加 `stats::aov`, `stats::complete.cases`, `stats::p.adjust`, `stats::rnorm`, `stats::setNames`, `stats::cor`
- 添加 `utils::flush.console`

### [MODIFY] 8 个 R 文件 — 添加 `utils::globalVariables()` 消除 NSE NOTE
- `visualization_barplot.R` → `mean_val`, `sd_val`
- `visualization_pca.R` → `PC1`, `PC2`, `Group`, `ShortLabel`
- `visualization_volcano.R` → `logFC`
- `visualization_apms.R` → `Group`, `Value`, `DisplayLabel`, `Label`, `MaxGroup`, `Mean`, `SD`, `RelAbundance`
- `visualization_cv.R` → `CV`, `MedianCV`, `group`, `label`
- `visualization_correlation.R` → `Sample1`, `Sample2`, `R2`
- `coverage.R` → `Coverage`, `MedianCov`
- `enrichment.R` → `p.adjust`, `Description`, `Count`

### 验证: `R CMD check` → 0 errors, 0 warnings

---

## Phase 2: 仓库卫生 🔴

### Step 1: 备份到 WPSCloud
```
/Users/wangxi/Documents/WPSCloud/msProteomiX_instant/
├── dev_backup/              ← 新建
│   ├── test_validation.R
│   ├── test_output/
│   ├── output/
│   ├── project_wiki/
│   └── 搜库软件输出结果demo/
```

### Step 2: `git rm --cached` 从仓库移除 (本地已有备份)
- `output/`, `test_output/`, `test_validation.R`, `搜库软件输出结果demo/`, `project_wiki/`, `.Rhistory`, `.DS_Store`

### Step 3: 更新 `.gitignore` + `.Rbuildignore`

---

## Phase 3: parsers.R 拆分 ✅

| 新文件 | 内容 | 预计行数 |
|--------|------|----------|
| `R/parsers.R` | `read_ms_data()`, `detect_engine()`, `new_MsDataSet()`, 公共辅助 | ~120 |
| `R/parser_fragpipe.R` | `parse_fragpipe()` + 辅助函数 | ~250 |
| `R/parser_spectronaut.R` | `parse_spectronaut()` + 辅助函数 | ~350 |
| `R/parser_maxquant.R` | `parse_maxquant()` + 辅助函数 | ~150 |
| `R/parser_pd_diann.R` | `parse_pd()` + `parse_diann()` + 辅助函数 | ~200 |

纯文件拆分，不改变任何函数签名或逻辑。

---

## Phase 4: 文档补全 ✅

### [NEW] `CHANGELOG.md` — 从 git log 生成版本记录
### [NEW] `vignettes/design_specification.md` — 纯 Markdown
- 架构设计 (MsDataSet, .msProteomiX_env)
- 引擎适配策略
- 用户工作流
- 扩展指南 (添加新引擎模板)
- 已知陷阱 (从 project_wiki 迁移)

### [MODIFY] `README.md` — 更新引擎状态 + badge + update_scripts 说明
### [NEW] `.github/CONTRIBUTING.md` — 开发规范 + 陷阱清单

---

## Phase 5: 编码硬化 ✅

- 全面扫描 `R/*.R` runtime strings 中的非 ASCII 字符
- 更新 `00_安装指南.R`: 删除 readr, 添加 magrittr/rlang

---

## 执行顺序

```
Phase 1 (R CMD check) → Phase 5 (编码) → Phase 2 (清理备份) → Phase 3 (拆分) → Phase 4 (文档)
```

预计总时间: ~85 min
