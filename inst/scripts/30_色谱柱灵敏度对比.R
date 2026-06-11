# ╔══════════════════════════════════════════════════════════════╗
# ║  msProteomiX — 步骤 30: 色谱柱灵敏度对比                   ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【实验设计】                                                ║
# ║   BSA+iRT：4 浓度 × 2 柱 × 2 重复 = 16 个 .d 文件         ║
# ║   HeLa：   3 上样量 × 2 柱 × 2 重复 = 12 个 .d 文件        ║
# ║                                                              ║
# ║  【前置条件】                                                ║
# ║   1. AlphaTims 脚本已运行：                                  ║
# ║      python wkdir/02_alphatims_bsa_irt_analysis.py           ║
# ║      -> 生成 wkdir/irt_metrics.csv + wkdir/bsa_metrics.csv  ║
# ║   2. FragPipe 搜库已完成：                                   ║
# ║      bash wkdir/01_fragpipe_batch_search.sh                  ║
# ║      -> fragpipe_out/BSA_all/ + fragpipe_out/HeLa_all/       ║
# ║                                                              ║
# ║  【输出】output/ 目录下的图表 + CSV + 汇总 HTML 报告         ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)
setup_workdir()

# ════════════════════════════════════════════════════════════════
# 用户配置区（根据实际路径修改）
# ════════════════════════════════════════════════════════════════
PROJECT_NAME <- "ColComparison_20260608"
IRT_CSV  <- "wkdir/irt_metrics.csv"   # AlphaTims iRT 输出
BSA_CSV  <- "wkdir/bsa_metrics.csv"   # AlphaTims BSA 输出
BSA_DIR  <- "wkdir/fragpipe_out/BSA_all"   # FragPipe BSA 搜库目录
HELA_DIR <- "wkdir/fragpipe_out/HeLa_all"  # FragPipe HeLa 搜库目录
OUTPUT   <- "output"
# ════════════════════════════════════════════════════════════════


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 1: 加载 AlphaTims 分析结果（BSA 特征肽 + iRT）
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
message("\n=== Step 1: Loading AlphaTims results ===")
alphatims_res <- import_alphatims_results(
  bsa_csv = BSA_CSV,
  irt_csv = IRT_CSV
)

# 保存到共享环境（供后续脚本使用）
.msProteomiX_env$alphatims_res <- alphatims_res


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 2: 加载 FragPipe 搜库结果（BSA + HeLa）
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
message("\n=== Step 2: Loading FragPipe results ===")
fragpipe_res <- merge_fragpipe_runs(
  bsa_dir  = BSA_DIR,
  hela_dir = HELA_DIR
)
.msProteomiX_env$fragpipe_res <- fragpipe_res


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 3: BSA LOD 曲线（检出肽段数 vs 浓度）
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
message("\n=== Step 3: BSA LOD Curve ===")
p_lod <- plot_lod_curve(
  bsa_df       = alphatims_res$bsa,
  output_dir   = OUTPUT,
  project_name = PROJECT_NAME
)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 4: BSA 强度线性（log-log 回归 R²）
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
message("\n=== Step 4: BSA Linearity ===")
p_lin <- plot_linearity(
  bsa_df       = alphatims_res$bsa,
  output_dir   = OUTPUT,
  project_name = PROJECT_NAME
)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 5: iRT 色谱柱对比（峰面积/FWHM/CV%）
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
message("\n=== Step 5: iRT Column Comparison ===")
p_irt <- plot_irt_column_compare(
  irt_df       = alphatims_res$irt,
  output_dir   = OUTPUT,
  project_name = PROJECT_NAME
)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 6: HeLa 上样量 vs 蛋白数曲线
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
message("\n=== Step 6: HeLa Loading Curve ===")
p_hela <- plot_hela_loading_curve(
  hela_ms      = fragpipe_res$hela,
  output_dir   = OUTPUT,
  project_name = PROJECT_NAME
)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 7: 计算所有 KPI 并输出汇总表
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
message("\n=== Step 7: KPI Summary Table ===")
kpi_table <- calc_column_metrics(
  alphatims_res = alphatims_res,
  fragpipe_res  = fragpipe_res
)
.msProteomiX_env$kpi_table <- kpi_table
if (nrow(kpi_table) > 0) {
  utils::write.csv(kpi_table,
    file.path(OUTPUT, paste0(PROJECT_NAME, "_KPI_summary.csv")),
    row.names = FALSE
  )
  message("  KPI table saved -> output/KPI_summary.csv")
}


message("\n=== Step 30 Done! Check output/ folder ===")
message("  Next: run 31_色谱柱对比_汇总报告.R")
