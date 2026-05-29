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
# ║   • output/GSEA_KEGG_*        — GSEA NES 气泡图             ║
# ║   • output/GSEA_ES_*_top*     — 经典富集曲线图 (合并)       ║
# ║   • output/GSEA_ES_*_pathway  — 单通路富集曲线图            ║
# ║   • output/GSEA_*.csv         — 富集结果表                  ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

# --- 自动定位工作目录 ---
library(msProteomiX)
setup_workdir()

# ━━━━━━━━━━━━━━━━━ 用户设置 (可修改) ━━━━━━━━━━━━━━━━━
gene_sets <- "go_bp"          # 基因集: go_bp(默认,完全离线), go_cc, go_mf, kegg, hallmark, reactome, kegg_online(联网)
organism <- "hsa"             # 物种: hsa(人), mmu(鼠)
top_n    <- 20                # NES 气泡图展示前 N 个通路
es_top_n <- 3                 # 富集曲线图: 各展示前 N 个激活/抑制通路
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# --- 主程序 ---
diff_result <- .msProteomiX_env$diff_result
if (is.null(diff_result)) {
  stop("Please run script 07 (Diff Analysis) first.")
}

project_name <- .msProteomiX_env$project_name %||% "Project"

# 1. GSEA 分析 (默认离线, 无需联网)
message("\n>>> Running GSEA analysis...")
gsea_result <- run_gsea_analysis(
  diff_result,
  organism  = organism,
  gene_sets = gene_sets
)

# 2. 可视化
if (!is.null(gsea_result) && nrow(gsea_result) > 0) {
  # 2a. NES 气泡图 (总览)
  message("\n>>> Plotting GSEA NES bubble chart...")
  p <- plot_gsea_result(
    gsea_result,
    top_n = top_n,
    output_dir = "output",
    project_name = project_name
  )

  # 2b. 经典富集曲线图 (running enrichment score)
  message("\n>>> Plotting GSEA enrichment curves...")
  plot_gsea_enrichment(
    gsea_result,
    top_n = es_top_n,
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
