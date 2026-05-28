# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 2: 定性分析柱状图                  ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【前置条件】已运行 01_数据导入与分组.R                     ║
# ║                                                              ║
# ║  【操作】点击 Source 按钮即可                                ║
# ║                                                              ║
# ║  【输出】output/Identification_ProteinGroups.pdf + .csv     ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)
setup_workdir()
check_prerequisites("01")

plot_id_barplot(
  ms_data      = .msProteomiX_env$ms_data,
  group_info   = .msProteomiX_env$group_info,
  project_name = "Project"
)

message("\n  Done! Check the output/ folder.")
