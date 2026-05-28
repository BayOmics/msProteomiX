# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 3: Venn 维恩图                     ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【前置条件】已运行 01_数据导入与分组.R                     ║
# ║  【操作】点击 Source 按钮即可                                ║
# ║  【输出】output/Venn_*.pdf + .csv                           ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)
setup_workdir()
check_prerequisites("01")

# Venn 图
plot_venn(
  ms_data      = .msProteomiX_env$ms_data,
  group_info   = .msProteomiX_env$group_info,
  project_name = "Project"
)

# UpSet 图
plot_upset(
  ms_data      = .msProteomiX_env$ms_data,
  group_info   = .msProteomiX_env$group_info,
  project_name = "Project"
)

message("\n  Done! Check the output/ folder.")
