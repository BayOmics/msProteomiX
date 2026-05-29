# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 15: Marker 基因表达分析            ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【前提条件】                                                ║
# ║   • 必须先运行 脚本 01 (数据导入与分组)                    ║
# ║                                                              ║
# ║  【说明】                                                    ║
# ║   • 展示指定 Marker 基因在各组中的表达量分布               ║
# ║   • 以箱线图 + 散点的形式展示                              ║
# ║   • 可自定义 Marker 基因列表                               ║
# ║                                                              ║
# ║  【输出文件】                                                ║
# ║   • output/Marker_Expression.*   — Marker 基因箱线图       ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

# --- 自动定位工作目录 ---
library(msProteomiX)
setup_workdir()

# ━━━━━━━━━━━━━━━━━ 用户设置 (必须修改!) ━━━━━━━━━━━━━━
# 在此列出你关心的 Marker 基因 (Gene Symbol)
# 示例: 常见癌症相关蛋白
marker_genes <- c(
  "TP53",      # 肿瘤抑制蛋白
  "EGFR",      # 表皮生长因子受体
  "MYC",       # 原癌基因
  "CDH1",      # E-cadherin
  "VIM",       # Vimentin
  "ACTB"       # Beta-actin (内参)
)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# --- 主程序 ---
ms_data    <- .msProteomiX_env$ms_data
group_info <- .msProteomiX_env$group_info
if (is.null(ms_data)) {
  stop("Please run script 01 (Data Import) first.")
}

project_name <- .msProteomiX_env$project_name %||% "Project"

# 绘制 Marker 表达箱线图
message("\n>>> Plotting marker gene expression...")
p <- plot_marker_boxplot(
  ms_data,
  group_info,
  markers      = marker_genes,
  output_dir   = "output",
  project_name = project_name
)

if (!is.null(p)) {
  message("\n========================================")
  message("  Marker gene analysis complete!")
  message("========================================")
} else {
  message("\n>>> No markers were found. Please check gene symbols.")
}
