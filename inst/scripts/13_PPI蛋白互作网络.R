# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 13: 蛋白互作/共表达网络             ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【前置条件】已运行 01 和 07 (需要差异分析结果)             ║
# ║                                                              ║
# ║  【分析模式】                                                ║
# ║   • "correlation" (默认, 完全离线)                           ║
# ║      基于蛋白丰度 Pearson 相关性构建共表达网络               ║
# ║   • "string" (需联网, 访问 STRING 数据库)                    ║
# ║      查询蛋白-蛋白相互作用 (PPI)                             ║
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
mode            <- "correlation" # "correlation"(离线) 或 "string"(联网)
cor_threshold   <- 0.8           # 相关系数阈值 (仅 correlation 模式)
species_id      <- 9606          # 9606 = 人, 10090 = 鼠 (仅 string 模式)
score_threshold <- 700           # STRING 阈值 (仅 string 模式)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

ms_data <- .msProteomiX_env$ms_data
result  <- .msProteomiX_env$diff_result
project_name <- .msProteomiX_env$project_name %||% "Project"

# --- 构建网络 ---
tryCatch({
  ppi_result <- run_ppi_network(
    diff_result     = result,
    ms_data         = ms_data,
    mode            = mode,
    cor_threshold   = cor_threshold,
    species         = species_id,
    score_threshold = score_threshold
  )

  # --- 可视化网络 ---
  if (!is.null(ppi_result)) {
    plot_ppi_network(
      ppi_result   = ppi_result,
      diff_result  = result,
      project_name = project_name
    )
  }
}, error = function(e) {
  message(">>> Network analysis skipped: ", e$message)
  if (mode == "string") {
    message("    Install with: BiocManager::install('STRINGdb')")
    message("    Also needs:   install.packages('igraph')")
  }
})

message("\n  Done! Check the output/ folder.")
