# ==============================================================================
# msProteomiX — 17. 血浆样本质量评估 (Contamination Index)
# ==============================================================================
# 功能: 评估血浆蛋白组学样本的细胞污染程度
# 指标: 血小板(PLT)、红细胞(RBC)、PBMC 的污染指数(CI)
# 原理: CI = SUM(标记蛋白强度) / SUM(全部蛋白强度)
# 参考: Geyer 2019, Gao 2026, Korff 2025
# ==============================================================================
#
# 使用说明:
#   1. 将此脚本放在 Spectronaut/FragPipe 输出目录中
#   2. 修改下方配置参数
#   3. 在 RStudio 中逐段运行
#
# 适用范围: 仅适用于人源血浆蛋白组学数据
# ==============================================================================

library(msProteomiX)

# ==============================================================================
# 配置参数 — 请根据实际情况修改
# ==============================================================================

data_dir   <- "."                    # Spectronaut/FragPipe 输出目录
fbs_dir    <- NULL                   # FBS 对照目录 (可选, BOVIN FASTA 搜库)
                                     # 例如: fbs_dir <- "../spectronaut_FBS_BOVIN"

project    <- "My_Plasma_Project"    # 项目名
instrument <- "Thermo Orbitrap"      # 仪器型号
search_mode <- "directDIA"           # 搜库模式 ("directDIA", "DDA", 等)


# ==============================================================================
# 1. 数据导入
# ==============================================================================

ms <- read_ms_data(data_dir)
message(sprintf(">>> Loaded: %d proteins x %d samples",
                nrow(ms$proteins), ncol(ms$proteins)))


# ==============================================================================
# 2. 污染指数 (CI) 计算
# ==============================================================================

ci <- calc_contamination_index(ms)

# 查看结果
print(ci[, c("Sample", "PGs", "CI_PLT", "CI_RBC", "CI_PBMC",
             "J_PLT", "J_RBC", "J_PBMC", "Overall")])

# 保存 CI 结果
dir.create("output", showWarnings = FALSE)
write.csv(ci, file.path("output", "contamination_index.csv"), row.names = FALSE)
message(">>> CI results saved to output/contamination_index.csv")


# ==============================================================================
# 3. CI 可视化
# ==============================================================================

# CI 柱状图 (3-面板: PLT / RBC / PBMC)
plot_ci_barplot(ci, output_dir = "output")

# CI 摘要图 (结果表格可视化)
plot_ci_summary(ci, output_dir = "output")

# Marker 蛋白强度热图 (需要 pheatmap 包)
plot_ci_heatmap(ms, ci, output_dir = "output")


# ==============================================================================
# 4. 生成评估报告 (可选)
# ==============================================================================

# 加载 FBS 对照数据 (如有)
ms_fbs <- if (!is.null(fbs_dir)) read_ms_data(fbs_dir) else NULL

# 生成中文报告
generate_plasma_report(
  ms, ci,
  ms_fbs      = ms_fbs,
  project     = project,
  instrument  = instrument,
  search_mode = search_mode,
  lang        = "zh",               # "en" 英文 | "zh" 中文
  output_dir  = "output"
)

# 如需英文版本:
# generate_plasma_report(ms, ci, ms_fbs = ms_fbs, project = project, lang = "en")

message(">>> Done! Check output/ directory for results.")
