# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 7: 差异分析与火山图                ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【前置条件】已运行 01_数据导入与分组.R                     ║
# ║                                                              ║
# ║  【操作步骤】                                                ║
# ║   1. 修改下方 "用户设置" 中的参数 (可选)                    ║
# ║   2. 点击 Source 按钮运行                                   ║
# ║   3. 按提示选择对照组和实验组                               ║
# ║   4. 查看 output/ 文件夹中的结果                            ║
# ║                                                              ║
# ║  【可调参数】                                                ║
# ║   • analysis_method — 差异分析方法                          ║
# ║     "limma" (推荐，更稳健) 或 "ttest" (经典t检验)           ║
# ║   • p_cutoff — P值阈值 (默认 0.05)                         ║
# ║   • fc_cutoff — Fold Change阈值 (默认 1，即2倍变化)        ║
# ║   • use_adj_p — 是否使用校正P值 (默认 FALSE)               ║
# ║   • label_top_n — 标注前N个差异蛋白名 (默认 20)            ║
# ║                                                              ║
# ║  【输出文件】                                                ║
# ║   • output/Volcano_<组名>.pdf — 火山图                      ║
# ║   • output/DiffExpr_<组名>.csv — 完整差异分析结果表         ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)
setup_workdir()
check_prerequisites("01")

# ━━━━━━━━━━━━━━━━━ 用户设置 (可修改) ━━━━━━━━━━━━━━━━━
analysis_method <- "limma"      # "limma" 或 "ttest"
p_cutoff        <- 0.05         # P值阈值
fc_cutoff       <- 1            # log2 Fold Change 阈值
use_adj_p       <- FALSE        # TRUE = 使用BH校正P值
label_top_n     <- 20           # 标注前N个差异蛋白
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# --- 主程序 ---
result <- run_diff_analysis(
  ms_data    = .msProteomiX_env$ms_data,
  group_info = .msProteomiX_env$group_info,
  method     = analysis_method,
  p_cutoff   = p_cutoff,
  fc_cutoff  = fc_cutoff,
  p_type     = ifelse(use_adj_p, "adj", "raw")
)

plot_volcano(
  result,
  label_top    = label_top_n,
  p_type       = ifelse(use_adj_p, "adj", "raw"),
  project_name = "Project"
)

message("\n  Done! Check the output/ folder.")
message("  To analyze another comparison, click Source again.")
