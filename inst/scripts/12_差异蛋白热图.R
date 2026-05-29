# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 12: 差异蛋白聚类热图               ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【前置条件】已运行 01 和 07 (需要差异分析结果)             ║
# ║                                                              ║
# ║  【可调参数】                                                ║
# ║   • top_n — 显示 Top N 差异蛋白 (默认 50)                  ║
# ║   • scale — 行标准化方法: "row" (Z-score) 或 "none"        ║
# ║                                                              ║
# ║  【输出】output/Heatmap_<对比组名>.pdf + .csv               ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)
setup_workdir()
check_prerequisites("07")

# ━━━━━━━━━━━━━━━━━ 用户设置 (可修改) ━━━━━━━━━━━━━━━━━
top_n  <- 50       # 显示 Top N 差异蛋白
scale  <- "row"    # "row" = Z-score 行标准化, "none" = 不标准化
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# --- 差异蛋白聚类热图 ---
plot_diff_heatmap(
  diff_result  = result,
  ms_data      = .msProteomiX_env$ms_data,
  group_info   = .msProteomiX_env$group_info,
  top_n        = top_n,
  scale        = scale,
  project_name = "Project"
)

message("\n  Done! Check the output/ folder.")
message("  To analyze another comparison, run 07 first then re-Source this script.")
