# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 6: 相关性热图                      ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【前置条件】已运行 01_数据导入与分组.R                     ║
# ║  【操作】点击 Source 按钮即可                                ║
# ║                                                              ║
# ║  【可调参数】                                                ║
# ║   • corr_method — 相关系数方法                              ║
# ║     "pearson" (默认，线性相关) 或 "spearman" (秩相关)       ║
# ║                                                              ║
# ║  【输出】output/CorrHeatmap_*.pdf + .csv                    ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)
check_prerequisites("01")

# ━━━━━━━━━━━━━━━━━ 用户设置 (可修改) ━━━━━━━━━━━━━━━━━
corr_method <- "pearson"    # "pearson" 或 "spearman"
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

plot_corr_heatmap(
  ms_data      = .msProteomiX_env$ms_data,
  group_info   = .msProteomiX_env$group_info,
  method       = corr_method,
  project_name = "Project"
)

message("\n  Done! Check the output/ folder.")
