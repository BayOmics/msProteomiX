# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 8: 序列覆盖度分析                  ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【前置条件】已运行 01_数据导入与分组.R                     ║
# ║  【额外需要】FASTA 数据库文件 (放入 wkdir/ 目录)           ║
# ║  【操作】修改下方 fasta_path 路径后，点击 Source            ║
# ║  【输出】output/Sequence_Coverage.pdf + .csv                ║
# ║                                                              ║
# ╠══════════════════════════════════════════════════════════════╣
# ║  【计算方法】                                                ║
# ║                                                              ║
# ║  序列覆盖度 = 被鉴定肽段覆盖的氨基酸残基数 / 蛋白总长 × 100%  ║
# ║                                                              ║
# ║  1. 从 combined_peptide.tsv 获取每个蛋白鉴定到的肽段        ║
# ║  2. 利用肽段的 Start/End 位置 (搜库引擎给出的精确位置)      ║
# ║     将肽段映射到蛋白全长序列上                               ║
# ║     (若无 Start/End 列, 则在 FASTA 序列中做字符串匹配)      ║
# ║  3. 用 boolean mask 标记每个被覆盖的残基位置                ║
# ║     (重叠区域只计算一次)                                     ║
# ║  4. 按组统计: 每组中任意样品检出的肽段均纳入计算            ║
# ║                                                              ║
# ║  输出: 每组每蛋白一个覆盖度值, 箱线图展示各组分布           ║
# ║  中位数标注在图上, 便于直观比较                              ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)
setup_workdir()
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
