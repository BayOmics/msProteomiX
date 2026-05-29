# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 14: GSEA 预排序富集分析            ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【前提条件】                                                ║
# ║   • 必须先运行 脚本 07 (差异分析)                           ║
# ║                                                              ║
# ║  【说明】                                                    ║
# ║   • GSEA 使用所有蛋白的 logFC 排序 (不限于显著蛋白)        ║
# ║   • 基于 KEGG 通路数据库 (gseKEGG)                         ║
# ║   • NES > 0: 通路在实验组中上调/激活                       ║
# ║   • NES < 0: 通路在实验组中下调/抑制                       ║
# ║                                                              ║
# ║  【输出文件】                                                ║
# ║   • output/GSEA_KEGG_*        — GSEA 气泡图                ║
# ║   • output/GSEA_KEGG_*.csv    — 富集结果表                 ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

# --- 自动定位工作目录 ---
library(msProteomiX)
setup_workdir()

# ━━━━━━━━━━━━━━━━━ 用户设置 (可修改) ━━━━━━━━━━━━━━━━━
org_db   <- "org.Hs.eg.db"   # 物种注释包: org.Hs.eg.db(人), org.Mm.eg.db(鼠)
organism <- "hsa"             # KEGG 物种: hsa(人), mmu(鼠)
top_n    <- 20                # 展示前 N 个通路
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# --- 主程序 ---
diff_result <- .msProteomiX_env$diff_result
if (is.null(diff_result)) {
  stop("Please run script 07 (Diff Analysis) first.")
}

project_name <- .msProteomiX_env$project_name %||% "Project"

# 1. GSEA 分析
message("\n>>> Running GSEA analysis...")
gsea_result <- run_gsea_analysis(
  diff_result,
  org_db   = org_db,
  organism = organism
)

# 2. 可视化
if (!is.null(gsea_result) && nrow(gsea_result) > 0) {
  message("\n>>> Plotting GSEA results...")
  p <- plot_gsea_result(
    gsea_result,
    top_n = top_n,
    output_dir = "output",
    project_name = project_name
  )

  # 保存结果到共享环境
  .msProteomiX_env$gsea_result <- gsea_result
} else {
  message(">>> No significant GSEA results. Try adjusting parameters.")
}

message("\n========================================")
message("  GSEA analysis complete!")
message("========================================")
