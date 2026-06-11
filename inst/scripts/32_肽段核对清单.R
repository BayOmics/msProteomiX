#!/usr/bin/env Rscript
# ==============================================================================
# 32_肽段核对清单.R
# 生成三份核对表，用于在 AlphaTims 中手动验证肽段检出情况
# ==============================================================================
# 输出文件（在 output/ 目录下）：
#   <Project>_BSA_verification_table.csv / .xlsx
#     — BSA 肽段：FragPipe PSM × AlphaTims EIC
#     — 列：序列 | m/z(obs/theo/calc) | ppm误差 | RT(两来源) | 峰强度 | 峰面积 | SNR | FWHM
#
#   <Project>_IRT_verification_table.csv / .xlsx
#     — iRT 肽段：FragPipe PSM × AlphaTims EIC（同上）
#
#   <Project>_HeLa_verification_table.csv / .xlsx
#     — HeLa 肽段：FragPipe combined_peptide.tsv（含各样本 MaxLFQ 强度列）
#     — 列：序列 | m/z | charge | RT | 各样本强度列 | 蛋白/基因
# ==============================================================================

library(msProteomiX)

# ====== 工作目录设置 ======
# 本脚本所有路径均相对于项目根目录（msProteomiX/），
# 请确保工作目录正确，或取消下行注释手动指定：
# setwd("/path/to/your/project")   # ← 按实际路径修改
#
# 自动检测：如果在 RStudio 中 Source 此脚本，工作目录可能不对
# 检查当前工作目录：getwd()
# ===========================

# ====== 用户配置区 ======
PROJECT_NAME <- "ColumnComparison"

# FragPipe 搜库输出目录（包含 psm.tsv 和 combined_peptide.tsv）
BSA_FRAGPIPE_DIR  <- "wkdir/fragpipe_out/BSA_all"
HELA_FRAGPIPE_DIR <- "wkdir/fragpipe_out/HeLa_all"

# AlphaTims 输出 CSV
BSA_AT_CSV  <- "wkdir/bsa_metrics.csv"
IRT_AT_CSV  <- "wkdir/irt_metrics.csv"

# 输出目录
OUTPUT_DIR  <- "output"

# PeptideProphet 最低概率阈值（0.9 = 90%）
MIN_PP_PROB <- 0.9
# =======================

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

cat("==============================================\n")
cat("  肽段核对清单生成\n")
cat("==============================================\n\n")

# --- [1] BSA 核对表（PSM + AlphaTims）---
cat("--- [1] BSA 核对表 ---\n")
if (dir.exists(BSA_FRAGPIPE_DIR)) {
  bsa_tbl <- make_verification_table(
    fragpipe_dir = BSA_FRAGPIPE_DIR,
    at_csv       = BSA_AT_CSV,
    table_type   = "bsa",
    output_dir   = OUTPUT_DIR,
    project_name = PROJECT_NAME,
    min_pp_prob  = MIN_PP_PROB
  )
  cat(sprintf("  BSA 肽段数: %d\n\n", nrow(bsa_tbl)))
} else {
  cat("  [跳过] BSA FragPipe 目录不存在:", BSA_FRAGPIPE_DIR, "\n\n")
}

# --- [2] iRT 核对表（PSM + AlphaTims）---
cat("--- [2] iRT 核对表 ---\n")
if (dir.exists(BSA_FRAGPIPE_DIR)) {
  irt_tbl <- make_verification_table(
    fragpipe_dir = BSA_FRAGPIPE_DIR,   # iRT 与 BSA 在同一次搜库中
    at_csv       = IRT_AT_CSV,
    table_type   = "irt",
    output_dir   = OUTPUT_DIR,
    project_name = PROJECT_NAME,
    min_pp_prob  = MIN_PP_PROB
  )
  cat(sprintf("  iRT 肽段数: %d\n\n", nrow(irt_tbl)))
} else {
  cat("  [跳过] BSA/iRT FragPipe 目录不存在:", BSA_FRAGPIPE_DIR, "\n\n")
}

# --- [3] HeLa 核对表（combined_peptide.tsv）---
cat("--- [3] HeLa 核对表 ---\n")
if (dir.exists(HELA_FRAGPIPE_DIR)) {
  hela_tbl <- make_verification_table(
    fragpipe_dir = HELA_FRAGPIPE_DIR,
    table_type   = "hela",
    output_dir   = OUTPUT_DIR,
    project_name = PROJECT_NAME,
    min_pp_prob  = MIN_PP_PROB
  )
  cat(sprintf("  HeLa 肽段数: %d\n\n", nrow(hela_tbl)))
} else {
  cat("  [跳过] HeLa FragPipe 目录不存在:", HELA_FRAGPIPE_DIR, "\n\n")
}

cat("==============================================\n")
cat("  输出文件位于:", OUTPUT_DIR, "\n")
cat("  BSA_verification_table.csv/.xlsx\n")
cat("  IRT_verification_table.csv/.xlsx\n")
cat("  HeLa_verification_table.csv/.xlsx\n")
cat("==============================================\n")
