# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 16: 生成分析报告                   ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【前提条件】                                                ║
# ║   • 运行完所有需要的分析脚本 (01-15)                        ║
# ║                                                              ║
# ║  【说明】                                                    ║
# ║   • 自动扫描 output/ 目录中的所有分析结果                   ║
# ║   • 按模块分类组织图表和数据表                              ║
# ║   • 生成 HTML 格式的综合报告                                ║
# ║                                                              ║
# ║  【输出文件】                                                ║
# ║   • output/{project}_Report.html — 综合分析报告             ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

# --- 自动定位工作目录 ---
library(msProteomiX)
setup_workdir()

# ━━━━━━━━━━━━━━━━━ 用户设置 (可修改) ━━━━━━━━━━━━━━━━━
output_dir   <- "output"          # 分析结果目录
format       <- "html"            # 报告格式: "html" 或 "pdf"
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# --- 主程序 ---
project_name <- .msProteomiX_env$project_name %||% "Project"

message("\n>>> Generating analysis report...")
report_path <- generate_report(
  output_dir   = output_dir,
  project_name = project_name,
  format       = format
)

message("\n========================================")
message("  Report generation complete!")
message(paste("  Report:", report_path))
message("========================================")
