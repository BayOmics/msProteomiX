# ╔══════════════════════════════════════════════════════════════╗
# ║  msProteomiX — Step 33: Contamination Ratio (CR) Evaluation ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  Evaluate single-charge contamination in timsTOF Pro data    ║
# ║  by analyzing IM-MS heatmaps across the LC gradient.         ║
# ║                                                              ║
# ║  Reference:                                                  ║
# ║    "Multimodal single cell-resolved spatial proteomics       ║
# ║     reveal pancreatic tumor heterogeneity" (Fig. 2f, 2g)     ║
# ║                                                              ║
# ║  【操作流程】                                                ║
# ║   Step 0: (可选) 从 .d/.hdf 提取前体离子并计算 CR            ║
# ║           - 需要 Python + alphatims                           ║
# ║           - 如已有 cr_output，跳过此步                        ║
# ║   Step 1: 自动发现 CR 结果文件                                ║
# ║   Step 2: 交互式分组 — 支持任意组数和命名                    ║
# ║   Step 3: 绘制 CR Gradient 曲线                              ║
# ║   Step 4: 绘制 TIC 强度曲线                                  ║
# ║   Step 5: 绘制 IM-m/z Density Map                            ║
# ║   Step 6: 生成 HTML 报告                                     ║
# ║                                                              ║
# ║  【输出】output/ 目录                                        ║
# ║   - CR_Gradient.pdf/png  — CR vs RT 曲线                    ║
# ║   - CR_TIC.pdf/png       — TIC 强度曲线                     ║
# ║   - CR_Heatmap_*.pdf/png — IM-m/z 密度图                    ║
# ║   - CR_Report.html       — 汇总 HTML 报告                   ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

library(msProteomiX)
setup_workdir()

# ════════════════════════════════════════════════════════════════
# 用户配置区（根据实际情况修改）
# ════════════════════════════════════════════════════════════════
PROJECT_NAME <- "CR_Analysis"
CR_DIR       <- "wkdir/cr_output"     # extract_precursors.py 输出目录
OUTPUT       <- "output"

# CR 分界线参数（论文默认值，通常无需修改）
LINE_POINT1 <- c(350, 0.8)   # (m/z, 1/K0)
LINE_POINT2 <- c(950, 1.3)   # (m/z, 1/K0)
# ════════════════════════════════════════════════════════════════


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 0 (Optional): Extract CR data from .d/.hdf files
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# If you already have CR results in CR_DIR, skip this step.
# Uncomment and modify the path below to run extraction from R:

# extract_cr_data(
#   input_path   = "D:/data/",          # .d folder, .hdf file, or directory
#   output_dir   = CR_DIR,
#   n_slices     = 10,
#   export_slices = TRUE
# )


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 1: Auto-discover CR results
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
message("\n=== Step 1: Auto-discovering CR results ===")

if (!dir.exists(CR_DIR)) {
  stop("CR output directory not found: ", CR_DIR,
       "\n  Either:\n",
       "  (a) Run extract_cr_data() in Step 0 above, or\n",
       "  (b) Run extract_precursors.py from command line first.")
}

# Find all *_cr_summary.csv files
cr_csvs <- list.files(CR_DIR, pattern = "_cr_summary\\.csv$",
                      full.names = TRUE)
if (length(cr_csvs) == 0) {
  stop("No *_cr_summary.csv files found in: ", CR_DIR)
}

# Extract sample names from filenames
sample_names <- sub("_cr_summary\\.csv$", "", basename(cr_csvs))
message(sprintf(">>> Found %d sample(s):", length(sample_names)))
for (i in seq_along(sample_names)) {
  message(sprintf("    [%d] %s", i, sample_names[i]))
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 2: Interactive grouping
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
message("\n=== Step 2: Sample grouping ===")
message(">>> Assign samples to groups (e.g., HeLa, DDM, FFHE).")
message("    Supports any number of groups.")

group_file <- file.path(CR_DIR, "cr_group_info.csv")

group_info <- interactive_grouping(
  sample_names = sample_names,
  group_file   = group_file,
  context      = "CR Analysis"
)

# Build group_map from interactive result
group_map <- split(group_info$sample_name, group_info$user_group)
# Convert to list of exact-match patterns
group_map <- lapply(group_map, function(samps) {
  paste0("^", gsub("([.|()\\^{}+$*?]|\\[|\\])", "\\\\\\1", samps), "$")
})

message(sprintf(">>> Groups defined: %s",
                paste(names(group_map), collapse = ", ")))


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 3: Import CR data with group assignment
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
message("\n=== Step 3: Importing CR results with groups ===")
cr_data <- import_cr_results(
  cr_dir    = CR_DIR,
  group_map = group_map
)

message("\n>>> CR data preview:")
print(head(cr_data))


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 4: Plot CR Gradient (Fig. 2f upper)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
message("\n=== Step 4: CR Gradient ===")
plot_cr_gradient(cr_data,
                 sd_factor = 0.5,
                 output_dir = OUTPUT,
                 project_name = PROJECT_NAME)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 5: Plot TIC Intensity (Fig. 2f lower)
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
message("\n=== Step 5: TIC Intensity ===")
plot_cr_tic(cr_data,
            output_dir = OUTPUT,
            project_name = PROJECT_NAME)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 6: Auto-discover and plot IM-m/z Density Maps
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
message("\n=== Step 6: IM-m/z Density Maps ===")

# Auto-discover ALL slice CSVs from all sample subdirectories
sample_dirs <- list.dirs(CR_DIR, recursive = FALSE, full.names = TRUE)
heatmap_csvs <- character(0)

for (sample_dir in sample_dirs) {
  found <- list.files(sample_dir, pattern = "_slice_rt.*\\.csv$",
                      full.names = TRUE)
  heatmap_csvs <- c(heatmap_csvs, found)
}

if (length(heatmap_csvs) > 0) {
  message(sprintf(">>> Found %d heatmap slice(s) across %d sample(s)",
                  length(heatmap_csvs), length(sample_dirs)))

  for (csv_path in heatmap_csvs) {
    message(sprintf("    Plotting: %s", basename(csv_path)))
    tryCatch({
      plot_cr_heatmap(csv_path,
                      line_point1 = LINE_POINT1,
                      line_point2 = LINE_POINT2,
                      output_dir = OUTPUT,
                      project_name = PROJECT_NAME)
    }, error = function(e) {
      message(sprintf("    [SKIP] %s: %s", basename(csv_path), e$message))
    })
  }
} else {
  message(">>> No heatmap slice CSVs found.")
  message("    To generate them, set export_slices = TRUE in extract_cr_data().")
}


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# Step 7: Generate HTML Report
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
message("\n=== Step 7: Generating HTML Report ===")
generate_cr_report(
  cr_data,
  heatmap_csvs = if (length(heatmap_csvs) > 0) heatmap_csvs else NULL,
  line_point1 = LINE_POINT1,
  line_point2 = LINE_POINT2,
  output_dir  = OUTPUT,
  project_name = PROJECT_NAME
)


# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
message("\n========================================")
message("  CR Analysis complete!")
message(sprintf("  Samples: %d | Groups: %s",
                length(sample_names),
                paste(names(group_map), collapse = ", ")))
message(sprintf("  Heatmaps: %d slice(s) plotted", length(heatmap_csvs)))
message("  Results saved to: output/")
message("========================================")
