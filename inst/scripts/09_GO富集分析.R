# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 9: GO 富集分析                     ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【前置条件】已运行 01_数据导入与分组.R                     ║
# ║                                                              ║
# ║  【可调参数】                                                ║
# ║   • species_db — 物种注释数据库                             ║
# ║     "org.Hs.eg.db" (人), "org.Mm.eg.db" (鼠)              ║
# ║   • go_category — GO 分类                                   ║
# ║     "CC" (细胞组分), "BP" (生物过程), "MF" (分子功能)      ║
# ║                                                              ║
# ║  【输出】output/GO_Enrichment.pdf + .csv                    ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)
check_prerequisites("01")

# ━━━━━━━━━━━━━━━━━ 用户设置 (可修改) ━━━━━━━━━━━━━━━━━
species_db  <- "org.Hs.eg.db"   # 人: "org.Hs.eg.db", 鼠: "org.Mm.eg.db"
go_category <- "CC"              # "CC", "BP", "MF"
top_n       <- 10                # 每组保留前N个结果
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

go_result <- run_go_enrichment(
  ms_data    = .msProteomiX_env$ms_data,
  group_info = .msProteomiX_env$group_info,
  org_db     = species_db,
  ont        = go_category,
  top_n      = top_n
)

plot_go_bubble(go_result, project_name = "Project")

message("\n  Done! Check the output/ folder.")
