# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 10: QC 方法评估                    ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【前置条件】已运行 01_数据导入与分组.R                     ║
# ║                                                              ║
# ║  【功能说明】                                                ║
# ║   生成 9 类 QC 评估图 (按组分面):                           ║
# ║   • 肽段长度分布 (Peptide Length)                            ║
# ║   • 前体离子电荷分布 (Charge Distribution)                  ║
# ║   • 漏切分布 (Missed Cleavage)                              ║
# ║   • 修饰类型分布 (Modification Types)                       ║
# ║   • 疏水性 GRAVY 分布                                       ║
# ║   • 等电点 pI 分布                                          ║
# ║   • M/Z 分布                                                ║
# ║   • 累积强度曲线 (动态范围)                                 ║
# ║   • 巯基肽占比 + 烷基化效率                                 ║
# ║                                                              ║
# ║  【数据要求】需要 PSM 水平数据 (FragPipe: psm.tsv)          ║
# ║   其他引擎暂不支持, 缺少PSM数据的图将自动跳过              ║
# ║                                                              ║
# ║  【输出】output/QC_*.pdf + .csv                             ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)
setup_workdir()
check_prerequisites("01")

# ━━━━━━━━━━━━━━━━━ 一键生成所有 QC 图 (按组) ━━━━━━━━━━━━━━━━━
plot_qc_panel(
  ms_data      = .msProteomiX_env$ms_data,
  group_info   = .msProteomiX_env$group_info,
  project_name = "Project"
)

message("\n  Done! Check the output/ folder.")
