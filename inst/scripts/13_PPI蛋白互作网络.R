# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 13: PPI 蛋白互作网络               ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【前置条件】已运行 01 和 07 (需要差异分析结果)             ║
# ║  【网络要求】需要联网访问 STRING 数据库                     ║
# ║                                                              ║
# ║  【可调参数】                                                ║
# ║   • species_id      — STRING 物种 NCBI ID                   ║
# ║     9606 (人), 10090 (鼠)                                   ║
# ║   • score_threshold — Combined Score 阈值 (0-1000)          ║
# ║     700 = 高可信度, 400 = 中等可信度                        ║
# ║                                                              ║
# ║  【输出】                                                    ║
# ║   output/PPI_Network_<对比组名>.pdf — 网络图                ║
# ║   output/PPI_Network_<对比组名>_edges.csv — 互作关系表      ║
# ║   output/PPI_Network_<对比组名>_nodes.csv — 节点信息表      ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)
setup_workdir()
check_prerequisites("07")

# ━━━━━━━━━━━━━━━━━ 用户设置 (可修改) ━━━━━━━━━━━━━━━━━
species_id      <- 9606    # 9606 = 人, 10090 = 鼠
score_threshold <- 700     # 700 = 高可信度, 400 = 中等可信度
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# --- 查询 STRING 数据库 ---
tryCatch({
  ppi_result <- run_ppi_network(
    diff_result     = result,
    species         = species_id,
    score_threshold = score_threshold
  )

  # --- 可视化 PPI 网络 ---
  if (!is.null(ppi_result)) {
    plot_ppi_network(
      ppi_result   = ppi_result,
      diff_result  = result,
      project_name = "Project"
    )
  }
}, error = function(e) {
  message(">>> PPI network skipped: ", e$message)
  message("    Install with: BiocManager::install('STRINGdb')")
  message("    Also needs:   install.packages('igraph')")
})

message("\n  Done! Check the output/ folder.")
message("  To analyze another comparison, run 07 first then re-Source this script.")
