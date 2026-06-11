# ╔══════════════════════════════════════════════════════════════╗
# ║  msProteomiX — 步骤 31: 色谱柱对比汇总报告                 ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【前置条件】已运行 30_色谱柱灵敏度对比.R                   ║
# ║                                                              ║
# ║  【功能说明】                                                ║
# ║   整合所有分析结果，生成 HTML 综合报告，包含：               ║
# ║   • 实验概述（样本 + 分析方法）                              ║
# ║   • BSA LOD 曲线（ColA vs ColB）                             ║
# ║   • BSA 强度线性（R² 对比）                                  ║
# ║   • iRT 3-panel（峰面积/FWHM/CV%）                           ║
# ║   • HeLa 蛋白数 vs 上样量曲线                               ║
# ║   • KPI 汇总表                                               ║
# ║   • 结论（自动判断推荐色谱柱）                               ║
# ║                                                              ║
# ║  【输出】output/<PROJECT>_ColComparison_Report.html          ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)
setup_workdir()

# ════════════════════════════════════════════════════════════════
# 用户配置区
# ════════════════════════════════════════════════════════════════
PROJECT_NAME <- "ColComparison_20260608"
OUTPUT       <- "output"
REPORT_FILE  <- file.path(OUTPUT, paste0(PROJECT_NAME, "_ColComparison_Report.html"))
# ════════════════════════════════════════════════════════════════


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 检查前置脚本是否已运行
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
if (!exists("alphatims_res", envir = .msProteomiX_env)) {
  stop("Please run 30_\u8272\u8c31\u67f1\u7075\u654f\u5ea6\u5bf9\u6bd4.R first.")
}

alphatims_res <- .msProteomiX_env$alphatims_res
fragpipe_res  <- .msProteomiX_env$fragpipe_res
kpi_table     <- .msProteomiX_env$kpi_table

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 生成 HTML 报告
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
message("\n=== Generating Column Comparison HTML Report ===")

# 收集已生成的图片路径
png_files <- list.files(OUTPUT, pattern = "\\.png$", full.names = TRUE)
png_files <- png_files[grepl(PROJECT_NAME, png_files)]

# 自动结论：比较 KPI
conclusion <- ""
if (!is.null(kpi_table) && nrow(kpi_table) > 0 &&
    all(c("ColA", "ColB") %in% colnames(kpi_table))) {
  # 判断 HeLa 蛋白数
  hela_row <- kpi_table[kpi_table$metric == "HeLa_mean_protein_IDs", ]
  if (nrow(hela_row) > 0) {
    n_a <- hela_row$ColA
    n_b <- hela_row$ColB
    winner <- if (!is.na(n_a) && !is.na(n_b)) {
      if (n_a > n_b * 1.05) "ColA" else if (n_b > n_a * 1.05) "ColB" else "comparable"
    } else "insufficient data"
    conclusion <- sprintf(
      "HeLa mean protein IDs: ColA = %s, ColB = %s. Recommendation: %s.",
      ifelse(is.na(n_a), "N/A", format(n_a, big.mark = ",")),
      ifelse(is.na(n_b), "N/A", format(n_b, big.mark = ",")),
      winner
    )
  }
}

# 构建 HTML
.write_col_report <- function(report_file, project_name,
                               kpi_table, png_files, conclusion) {
  # KPI 表格 HTML
  kpi_html <- if (!is.null(kpi_table) && nrow(kpi_table) > 0) {
    rows_html <- paste(apply(kpi_table, 1, function(r) {
      sprintf("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
              r["metric"], r["unit"],
              ifelse(is.na(r["ColA"]), "—", r["ColA"]),
              ifelse(is.na(r["ColB"]), "—", r["ColB"]))
    }), collapse = "\n")
    sprintf(
      "<table class='kpi-table'><thead><tr>
       <th>Metric</th><th>Unit</th><th>ColA</th><th>ColB</th>
       </tr></thead><tbody>%s</tbody></table>", rows_html
    )
  } else "<p><em>KPI table not available (run Step 30 first).</em></p>"

  # 图片 HTML
  img_html <- if (length(png_files) > 0) {
    paste(vapply(sort(png_files), function(f) {
      fname <- basename(f)
      # 使用相对路径（仅文件名），确保 HTML 在 output/ 目录旁可移植
      sprintf(
        "<div class='figure'><h3>%s</h3>\n         <img src='%s' alt='%s' style='max-width:100%%;'></div>",
        tools::file_path_sans_ext(fname), fname, fname
      )
    }, character(1)), collapse = "\n")
  } else "<p><em>No figures found. Run 30_\u8272\u8c31\u67f1\u7075\u654f\u5ea6\u5bf9\u6bd4.R to generate plots.</em></p>"

  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

  html <- sprintf('<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>%s — Column Comparison Report</title>
<style>
  body { font-family: "Helvetica Neue", Arial, sans-serif; max-width: 1200px;
         margin: 0 auto; padding: 20px; background: #f8f9fa; color: #333; }
  h1   { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
  h2   { color: #34495e; margin-top: 40px; }
  h3   { color: #555; }
  .summary-box { background: #fff; border-left: 4px solid #3498db;
                 padding: 16px 24px; margin: 20px 0; border-radius: 4px;
                 box-shadow: 0 2px 4px rgba(0,0,0,.08); }
  .kpi-table { border-collapse: collapse; width: 100%%; margin: 20px 0; }
  .kpi-table th { background: #3498db; color: #fff; padding: 10px 16px; text-align: left; }
  .kpi-table td { padding: 8px 16px; border-bottom: 1px solid #e0e0e0; }
  .kpi-table tr:nth-child(even) { background: #f5f5f5; }
  .figure { background: #fff; padding: 20px; margin: 20px 0;
            border-radius: 8px; box-shadow: 0 2px 8px rgba(0,0,0,.1); }
  .conclusion { background: #eafaf1; border-left: 4px solid #27ae60;
                padding: 16px 24px; margin: 20px 0; border-radius: 4px; }
  .footer { color: #aaa; font-size: 12px; margin-top: 60px;
            border-top: 1px solid #ddd; padding-top: 10px; }
</style>
</head>
<body>
<h1>%s — Nano-LC Column Sensitivity Comparison Report</h1>

<div class="summary-box">
<h2>Experiment Overview</h2>
<p><strong>Instrument:</strong> TimsTOF Pro, DDA-PASEF</p>
<p><strong>BSA+iRT:</strong> 4 concentrations (0.02, 0.2, 2, 20 fmol) &times; 2 columns &times; 2 replicates = 16 raw files</p>
<p><strong>HeLa Digest:</strong> 3 loading amounts (0.25, 0.5, 1 ng) &times; 2 columns &times; 2 replicates = 12 raw files</p>
<p><strong>iRT Kit:</strong> Biognosys iRT-Kit (11 synthetic peptides, irtfusion.fasta)</p>
<p><strong>Analysis Tools:</strong> AlphaTims (EIC extraction) + FragPipe (database search) + msProteomiX (statistics)</p>
<p><strong>FASTA:</strong> BSA P02769 + irtfusion + cRAP (117 entries) | Human SwissProt + cRAP (20,546 entries)</p>
<p><strong>Report generated:</strong> %s</p>
</div>

<h2>KPI Summary</h2>
%s

%s

<div class="conclusion">
<h2>Conclusion</h2>
<p>%s</p>
</div>

<div class="footer">Generated by msProteomiX v0.2.0 | %s</div>
</body></html>',
    project_name, project_name, timestamp,
    kpi_html, img_html,
    if (nchar(conclusion) > 0) conclusion else "Run Step 30 to generate KPI data for automatic conclusion.",
    timestamp
  )

  writeLines(html, report_file)
  message(sprintf("\u2705 HTML report saved: %s", report_file))
}

.write_col_report(REPORT_FILE, PROJECT_NAME, kpi_table, png_files, conclusion)

message(sprintf("\n=== Done! Open the report: %s ===", REPORT_FILE))
