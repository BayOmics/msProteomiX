# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 5: CV 变异系数分析                 ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【前置条件】已运行 01_数据导入与分组.R                     ║
# ║  【操作】点击 Source 按钮即可                                ║
# ║  【输出】output/CV_Boxplot.pdf + CV_Violin.pdf              ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)
check_prerequisites("01")

# 箱线图
plot_cv_boxplot(
  ms_data      = .msProteomiX_env$ms_data,
  group_info   = .msProteomiX_env$group_info,
  project_name = "Project"
)

# 小提琴图
plot_cv_violin(
  ms_data      = .msProteomiX_env$ms_data,
  group_info   = .msProteomiX_env$group_info,
  project_name = "Project"
)

message("\n  Done! Check the output/ folder.")
