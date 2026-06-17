# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 18: 数据总结表                     ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【前置条件】已运行 01_数据导入与分组.R                     ║
# ║  【操作】点击 Source 按钮即可                                ║
# ║  【输出】output/Project_Summary_Table.csv                   ║
# ║                                                              ║
# ║  【包含指标】                                                ║
# ║   - PSM: 谱图-肽段匹配数量                                  ║
# ║   - Peptide: 鉴定到的肽段数                                 ║
# ║   - Protein_Groups: 鉴定到的蛋白组数                        ║
# ║   - Zero_Miss_Pct: 0漏切率 (%)                              ║
# ║   - Cys_Peptide_Pct: 含半胱氨酸肽段比例 (%)                 ║
# ║   - Alkylation_Pct: 烷基化效率 (%)                          ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)
setup_workdir()
check_prerequisites("01")

summary_df <- generate_summary_table(
  ms_data      = .msProteomiX_env$ms_data,
  group_info   = .msProteomiX_env$group_info,
  project_name = "Project"
)

# Print to console for quick review
print(summary_df)

message("\n  Done! Check the output/ folder for Summary_Table.csv")
