# ╔══════════════════════════════════════════════════════════════╗
# ║       msProteomiX — 步骤 20: AP-MS 时序分析                ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【适用场景】                                                ║
# ║   亲和纯化质谱 (AP-MS) 实验，如:                            ║
# ║   • Co-IP / Flag-tag / HA-tag pulldown                      ║
# ║   • APEX2 proximity labeling (时序实验)                     ║
# ║   • BioID proximity labeling                                ║
# ║                                                              ║
# ║  【操作步骤】                                                ║
# ║   1. 先运行 01_数据导入与分组.R 完成数据导入                ║
# ║   2. 修改下方 bait 基因名                                   ║
# ║   3. 点击 "Source" 运行                                     ║
# ║                                                              ║
# ║  【输出文件】                                                ║
# ║   output/APMS_DotPlot.*           — Dot Plot 气泡图         ║
# ║   output/APMS_Heatmap.*           — 相对丰度热图            ║
# ║   output/APMS_Abundance_Curves.*  — 蛋白丰度曲线            ║
# ║   output/APMS_Mfuzz_*             — 时序聚类结果            ║
# ║   output/APMS_ANOVA_results.csv   — ANOVA 统计结果          ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)
setup_workdir()

# ━━━━━━━━━━━━━━━━━ 用户设置 (请修改) ━━━━━━━━━━━━━━━━━
bait       <- "UBASH3B"     # ★ Bait 蛋白的基因名 (如 "CBL", "UBASH3B")
p_cutoff   <- 0.05          # ANOVA 显著性阈值
n_clusters <- 4             # Mfuzz 聚类数
do_crapome <- TRUE          # 是否过滤 CRAPome 污染蛋白
do_mfuzz   <- TRUE          # 是否做 Mfuzz 时序聚类
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# --- 前置检查 ---
check_prerequisites("Data Import")
ms_data    <- .msProteomiX_env$ms_data
group_info <- .msProteomiX_env$group_info

if (!dir.exists("output")) dir.create("output")

# =============================================
# 1. CRAPome 污染蛋白标记
# =============================================
if (do_crapome) {
  message("\n>>> Step 1: CRAPome filtering...")
  ms_data <- filter_crapome(ms_data, action = "remove")
}

# =============================================
# 2. Bait 归一化
# =============================================
message("\n>>> Step 2: Bait normalization...")
ms_norm <- normalize_bait(ms_data, bait_gene = bait)

# =============================================
# 3. 计算相对丰度
# =============================================
message("\n>>> Step 3: Relative abundance...")
rel_abund <- calc_relative_abundance(ms_norm, group_info)

# =============================================
# 4. ANOVA 时序检验
# =============================================
message("\n>>> Step 4: ANOVA test...")
anova_result <- run_anova_timecourse(ms_norm, group_info, p_cutoff = p_cutoff)

# 保存 ANOVA 结果
write.csv(anova_result, "output/APMS_ANOVA_results.csv", row.names = FALSE)

# 提取显著蛋白
sig_ids <- anova_result$Protein[anova_result$significant == "YES"]
message(sprintf("  Significant proteins: %d / %d", length(sig_ids), nrow(anova_result)))

# =============================================
# 5. Dot Plot
# =============================================
message("\n>>> Step 5: Dot Plot...")
plot_apms_dotplot(rel_abund, sig_proteins = sig_ids, project_name = "APMS")

# =============================================
# 6. 热图
# =============================================
message("\n>>> Step 6: Heatmap...")
plot_apms_heatmap(rel_abund, sig_proteins = sig_ids, project_name = "APMS")

# =============================================
# 7. 丰度曲线 (显著蛋白)
# =============================================
message("\n>>> Step 7: Abundance curves...")
plot_abundance_curve(rel_abund, proteins = sig_ids, project_name = "APMS")

# =============================================
# 8. Mfuzz 时序聚类 (可选)
# =============================================
if (do_mfuzz && length(sig_ids) >= n_clusters) {
  message("\n>>> Step 8: Mfuzz clustering...")
  mfuzz_res <- run_mfuzz_cluster(ms_norm, group_info,
                                  sig_proteins = sig_ids,
                                  n_clusters = n_clusters)
  plot_mfuzz(mfuzz_res, project_name = "APMS")
} else if (do_mfuzz) {
  message("\n>>> Step 8: Skipped (too few significant proteins for clustering)")
}

# =============================================
# 完成
# =============================================
message("\n========================================")
message("  AP-MS analysis complete!")
message(sprintf("  Bait: %s", bait))
message(sprintf("  Significant: %d / %d proteins", length(sig_ids), nrow(anova_result)))
message("  Check output/ folder for results.")
message("========================================")
