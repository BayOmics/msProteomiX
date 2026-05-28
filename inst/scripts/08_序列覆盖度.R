# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 8: 序列覆盖度分析                  ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【前置条件】已运行 01_数据导入与分组.R                     ║
# ║  【额外需要】FASTA 数据库文件 (放入 wkdir/ 目录)           ║
# ║  【操作】修改下方 fasta_path 路径后，点击 Source            ║
# ║  【输出】output/Sequence_Coverage.pdf + .csv                ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)
check_prerequisites("01")

# ━━━━━━━━━━━━━━━━━ 用户设置 (可修改) ━━━━━━━━━━━━━━━━━
fasta_path <- ""  # 填入 FASTA 文件路径，如 "wkdir/FP_xxx/UP000005640.fasta"
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

if (fasta_path == "") {
  # 自动搜索
  fasta_files <- list.files("wkdir", pattern = "\\.(fasta|fa)$",
                             recursive = TRUE, full.names = TRUE)
  if (length(fasta_files) > 0) {
    fasta_path <- fasta_files[1]
    message(paste("  Auto-detected FASTA:", fasta_path))
  } else {
    stop("Please specify fasta_path or place a .fasta file in wkdir/")
  }
}

cov_df <- calc_coverage(
  ms_data    = .msProteomiX_env$ms_data,
  fasta_path = fasta_path,
  group_info = .msProteomiX_env$group_info
)

plot_coverage(cov_df, project_name = "Project")

message("\n  Done! Check the output/ folder.")
