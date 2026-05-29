# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 11: KEGG / Reactome 通路富集       ║
# ╠══════════════════════════════════════════════════════════════╣
# ║  ⚠ 注意: KEGG 富集需要联网! Reactome 使用本地数据库(离线)  ║
# ║  离线替代: 使用脚本 10 (GO富集) 或 14 (GSEA, 默认离线)     ║
# ║  【前置条件】已运行 01_数据导入与分组.R                     ║
# ║   可选: 已运行 07_差异分析 (用于差异蛋白富集)               ║
# ║                                                              ║
# ║  【可调参数】                                                ║
# ║   • species_db   — 物种注释数据库                           ║
# ║     "org.Hs.eg.db" (人), "org.Mm.eg.db" (鼠)              ║
# ║   • species_kegg — KEGG 物种代码                            ║
# ║     "hsa" (人), "mmu" (鼠)                                 ║
# ║   • use_diff     — 是否使用差异蛋白 (需先运行07)           ║
# ║                                                              ║
# ║  【输出】output/KEGG_Enrichment.pdf + .csv                  ║
# ║         output/Reactome_Enrichment.pdf + .csv               ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)
setup_workdir()
check_prerequisites("01")

# ━━━━━━━━━━━━━━━━━ 用户设置 (可修改) ━━━━━━━━━━━━━━━━━
species_db   <- "org.Hs.eg.db"   # 人: "org.Hs.eg.db", 鼠: "org.Mm.eg.db"
species_kegg <- "hsa"             # 人: "hsa", 鼠: "mmu"
top_n_kegg   <- 15                # KEGG 每组保留前N个结果
top_n_react  <- 20                # Reactome 保留前N个结果
use_diff     <- FALSE             # TRUE = 仅对差异蛋白富集 (需先运行07)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# --- KEGG Pathway ---
message("\n>>> Starting KEGG Pathway Enrichment...")

if (use_diff) {
  if (!exists("result")) {
    stop("Diff analysis result not found. Run 07 first, or set use_diff <- FALSE.")
  }
  kegg_result <- run_kegg_enrichment(
    diff_result = result,
    org_db      = species_db,
    organism    = species_kegg,
    top_n       = top_n_kegg
  )
} else {
  kegg_result <- run_kegg_enrichment(
    ms_data    = .msProteomiX_env$ms_data,
    group_info = .msProteomiX_env$group_info,
    org_db     = species_db,
    organism   = species_kegg,
    top_n      = top_n_kegg
  )
}

plot_kegg_bar(kegg_result, project_name = "Project")

# --- Reactome Pathway ---
message("\n>>> Starting Reactome Pathway Enrichment...")

tryCatch({
  if (use_diff) {
    reactome_result <- run_reactome_enrichment(
      diff_result = result,
      org_db      = species_db,
      top_n       = top_n_react
    )
  } else {
    reactome_result <- run_reactome_enrichment(
      ms_data    = .msProteomiX_env$ms_data,
      group_info = .msProteomiX_env$group_info,
      org_db     = species_db,
      top_n      = top_n_react
    )
  }

  plot_reactome_bar(reactome_result, project_name = "Project")
}, error = function(e) {
  message(">>> Reactome skipped: ", e$message)
  message("    Install with: BiocManager::install('ReactomePA')")
})

message("\n  Done! Check the output/ folder.")
