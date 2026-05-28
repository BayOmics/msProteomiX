# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 1: 数据导入与样品分组              ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【操作步骤】                                                ║
# ║   1. 将搜库软件的结果文件夹放入 wkdir/ 目录                 ║
# ║      • FragPipe:    复制整个输出文件夹                      ║
# ║      • Spectronaut: 复制 *_Report.tsv (Run Pivot 导出)      ║
# ║      • MaxQuant:    复制 txt/ 文件夹                        ║
# ║      • DIA-NN:      复制 report.pg_matrix.tsv               ║
# ║                                                              ║
# ║   2. 在 RStudio 中点击右上角 "Source" 按钮运行              ║
# ║                                                              ║
# ║   3. 按照屏幕提示进行分组:                                  ║
# ║      • 输入样品编号范围 (如 1-3)，按回车                    ║
# ║      • 输入组名 (如 Control)，按回车                        ║
# ║      • 重复上述步骤直到所有样品分组完毕                     ║
# ║      • 输入 y 确认完成                                      ║
# ║                                                              ║
# ║  【输出文件】                                                ║
# ║   • wkdir/group_info.csv  — 分组信息 (后续脚本自动读取)    ║
# ║                                                              ║
# ║  【注意事项】                                                ║
# ║   • 这一步只需运行一次，后续脚本会自动读取分组文件         ║
# ║   • 如需修改分组，重新运行此脚本即可                       ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)

# ━━━━━━━━━━━━━━━━━ 用户设置 (可修改) ━━━━━━━━━━━━━━━━━
data_dir <- "wkdir"            # 数据文件夹路径
engine   <- "auto"             # 搜库引擎: "auto", "fragpipe", "spectronaut", "maxquant", "pd", "diann"
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# --- 主程序 (不需要修改) ---
if (!dir.exists(data_dir)) dir.create(data_dir)
if (!dir.exists("output")) dir.create("output")

# 1. 读取数据
message(">>> Reading data...")
ms_data <- read_ms_data(data_dir, engine = engine)
print(ms_data)

# 2. 交互式分组
message("\n>>> Starting grouping wizard...")
group_info <- interactive_grouping(
  sample_names = get_sample_names(ms_data),
  group_file   = file.path(data_dir, "group_info.csv"),
  context      = "Data Import"
)

# 3. 保存到共享环境 (供后续脚本使用)
.msProteomiX_env$ms_data    <- ms_data
.msProteomiX_env$group_info <- group_info

message("\n========================================")
message("  Data import and grouping complete!")
message(paste("  Samples:", length(get_sample_names(ms_data))))
message(paste("  Groups:", paste(unique(group_info$user_group), collapse = ", ")))
message("  Next: Open 02 - 10 scripts to run analyses.")
message("========================================")
