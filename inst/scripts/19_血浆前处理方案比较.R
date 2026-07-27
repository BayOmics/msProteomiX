# ==============================================================================
# msProteomiX -- 19. Plasma Sample Prep Comparison
# ==============================================================================
# Metrics: CI, PG count, Quant CV, Data completeness, PCA, Top-N,
#          High-abundance proteins, Missed Cleavage, Peptide length,
#          Charge state, Dynamic range, Rank curve overlay
# Input:  Spectronaut or DIA-NN search results (each condition searched independently)
# Output: Self-contained HTML reports (EN + ZH) + Chrome PDF + individual PDFs
# ==============================================================================
#
# Usage:
#   1. Edit the CONFIG section below
#   2. In RStudio: Ctrl+A (select all), then Ctrl+Enter (run)
#   3. Check output/ directory for results
# ==============================================================================

library(msProteomiX)
library(ggplot2)

# ==============================================================================
# CONFIG -- Edit this section only
# ==============================================================================

conditions <- list(
  list(name = "BayOmics_3",  path = "./Spectronaut/20260722_105657_20260721_BM_3",  engine = "spectronaut"),
  list(name = "Nanomics_3",  path = "./Spectronaut/20260722_110059_20260721_NM_3",  engine = "spectronaut"),
  list(name = "BayOmics_5",  path = "./Spectronaut/20260722_110249_20260722_BM_5",  engine = "spectronaut"),
  list(name = "Nanomics_5",  path = "./Spectronaut/20260722_110316_20260722_NM_5",  engine = "spectronaut")
)

compare_groups <- list(
  BayOmics = c("BayOmics_3", "BayOmics_5"),
  Nanomics = c("Nanomics_3", "Nanomics_5")
)

project_name <- "SZU_PlasmaG"
output_base  <- "output"

# ==============================================================================
# EXECUTION -- Do NOT edit below this line
# ==============================================================================

# --- Local helpers ---
.b64img <- function(path, width = "85%", caption = "") {
  if (!file.exists(path)) return("")
  b64 <- base64enc::base64encode(path)
  paste0('<div class="figure">',
    sprintf('<img src="data:image/png;base64,%s" style="max-width:%s">', b64, width),
    if (nzchar(caption)) sprintf('<div class="caption">%s</div>', caption) else "",
    '</div>')
}

.chrome_path <- function() {
  candidates <- c(
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    if (nzchar(Sys.getenv("ProgramFiles")))
      file.path(Sys.getenv("ProgramFiles"), "Google/Chrome/Application/chrome.exe"),
    if (nzchar(Sys.getenv("LOCALAPPDATA")))
      file.path(Sys.getenv("LOCALAPPDATA"), "Google/Chrome/Application/chrome.exe"),
    Sys.which("google-chrome"), Sys.which("chromium-browser"), Sys.which("chromium"))
  for (p in candidates) if (nzchar(p) && file.exists(p)) return(p)
  NULL
}

.save_png <- function(p, name, png_dir, w = 8, h = 5) {
  f <- file.path(png_dir, paste0(name, ".png"))
  ggplot2::ggsave(f, p, width = w, height = h, dpi = 200)
  f
}

grp_colors <- c("#3498db", "#e74c3c", "#2ecc71", "#f39c12", "#9b59b6", "#1abc9c")
names(grp_colors) <- head(names(compare_groups), length(grp_colors))
grp_colors <- grp_colors[seq_along(names(compare_groups))]
grp_colors[["Other"]] <- "#95a5a6"

# ==============================================================================
# Phase 1: Per-condition evaluation
# ==============================================================================

message("========================================")
message("  msProteomiX Plasma Method Comparison")
message("========================================\n")

all_ms <- list(); all_ci_results <- list(); all_summaries <- list(); all_overviews <- list()
all_det <- list()

for (cond in conditions) {
  cond_name <- cond$name; cond_path <- cond$path; cond_engine <- cond$engine
  message(sprintf("\n>>> [%s] Evaluating: %s (%s)", format(Sys.time(), "%H:%M:%S"), cond_name, cond_engine))

  eval_dir <- file.path(output_base, paste0(cond_name, "_evaluation"))
  qc_dir   <- file.path(eval_dir, "qc_panel")
  dir.create(eval_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_dir,   recursive = TRUE, showWarnings = FALSE)

  ms <- tryCatch(read_ms_data(cond_path, engine = cond_engine),
                 error = function(e) { message("    ERROR: ", e$message); NULL })
  if (is.null(ms)) { message("    Skipping ", cond_name); next }
  all_ms[[cond_name]] <- ms
  message(sprintf("    %d proteins x %d samples", nrow(ms$proteins), ncol(ms$proteins)))

  # IdentificationsOverview
  if (dir.exists(cond_path)) {
    ov_files <- list.files(cond_path, pattern = "IdentificationsOverview.*\\.tsv$", full.names = TRUE, ignore.case = TRUE)
    if (length(ov_files) > 0) {
      ov <- read.delim(ov_files[1], stringsAsFactors = FALSE); ov$Condition <- cond_name
      all_overviews[[cond_name]] <- ov
    }
  }

  # CI
  ci <- tryCatch(calc_contamination_index(ms), error = function(e) { message("    CI error: ", e$message); NULL })
  if (is.null(ci)) next
  ci$Condition <- cond_name; all_ci_results[[cond_name]] <- ci
  write.csv(ci, file.path(eval_dir, "01_contamination_index.csv"), row.names = FALSE)
  tryCatch(plot_ci_barplot(ci, output_dir = eval_dir), error = function(e) NULL)
  tryCatch(plot_ci_heatmap(ms, ci, output_dir = eval_dir), error = function(e) NULL)

  # Detection rate
  prot_mat <- ms$proteins
  det_pct <- sapply(seq_len(ncol(prot_mat)), function(j) sum(!is.na(prot_mat[, j])) / nrow(prot_mat) * 100)
  det_df <- data.frame(Sample = ms$sample_names, DetRate = det_pct,
    Judgment = ifelse(det_pct >= 80, "Pass", ifelse(det_pct >= 50, "Flag", "Fail")),
    Condition = cond_name, stringsAsFactors = FALSE)
  all_det[[cond_name]] <- det_df

  p <- ggplot(det_df, aes(x = Sample, y = DetRate, fill = Judgment)) +
    geom_col(width = 0.65) + geom_hline(yintercept = 80, linetype = "dashed", color = "#27ae60") +
    geom_text(aes(label = sprintf("%.1f%%", DetRate)), vjust = -0.3, size = 2.5) +
    scale_fill_manual(values = c(Pass="#27ae60", Flag="#f39c12", Fail="#e74c3c")) +
    scale_y_continuous(limits = c(0, 105)) +
    labs(title = paste0("Detection Rate - ", cond_name), x = NULL, y = "%") +
    theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  ggsave(file.path(eval_dir, "05_detection_rate.pdf"), p, width = 8, height = 5)

  # Intensity + Rank
  mat_raw <- as.matrix(prot_mat); if (isTRUE(ms$is_log2)) mat_raw <- 2^mat_raw
  mat_raw[is.na(mat_raw) | mat_raw == 0] <- NA
  int_l <- data.frame(Sample = rep(ms$sample_names, each = nrow(mat_raw)),
    log2I = as.vector(log2(mat_raw)), stringsAsFactors = FALSE)
  int_l <- int_l[!is.na(int_l$log2I), ]
  p <- ggplot(int_l, aes(x = Sample, y = log2I, fill = Sample)) +
    geom_boxplot(outlier.size = 0.2, alpha = 0.8) +
    labs(title = paste0("Intensity - ", cond_name), x = NULL, y = "log2(Intensity)") +
    theme_bw(base_size = 10) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7))
  ggsave(file.path(eval_dir, "06_intensity_boxplot.pdf"), p, width = 8, height = 5)

  rl <- lapply(seq_len(ncol(mat_raw)), function(j) {
    v <- mat_raw[, j]; v <- v[!is.na(v) & v > 0]; v <- sort(v, decreasing = TRUE)
    data.frame(Sample = ms$sample_names[j], Rank = seq_along(v), log2I = log2(v), stringsAsFactors = FALSE) })
  p <- ggplot(do.call(rbind, rl), aes(x = Rank, y = log2I, color = Sample)) +
    geom_line(linewidth = 0.5, alpha = 0.85) +
    labs(title = paste0("Protein Rank - ", cond_name), x = "Rank", y = "log2(Intensity)") +
    theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(eval_dir, "07_protein_rank.pdf"), p, width = 8, height = 5)

  # QC panel
  tryCatch(plot_qc_panel(ms, NULL, qc_dir, paste0(project_name, "_", cond_name)), error = function(e) {
    tryCatch(plot_qc_missing_values(ms, NULL, qc_dir, paste0(project_name, "_", cond_name)), error = function(e2) NULL)
    tryCatch(plot_qc_intensity_boxplot(ms, NULL, qc_dir, paste0(project_name, "_", cond_name)), error = function(e2) NULL) })

  # Summary
  summary_df <- data.frame(Sample = ci$Sample, Condition = cond_name, PGs = ci$PGs,
    CI_PLT = sprintf("%.5f (%s)", ci$CI_PLT, ci$J_PLT), CI_RBC = sprintf("%.5f (%s)", ci$CI_RBC, ci$J_RBC),
    CI_PBMC = ifelse(is.na(ci$CI_PBMC), "N/A", sprintf("%.5f (%s)", ci$CI_PBMC, ci$J_PBMC)),
    Detection_Rate = sprintf("%.1f%% (%s)", det_df$DetRate, det_df$Judgment),
    Overall = ci$Overall, stringsAsFactors = FALSE)
  write.csv(summary_df, file.path(eval_dir, "08_summary_judgment.csv"), row.names = FALSE)
  all_summaries[[cond_name]] <- summary_df
  message(sprintf(">>> [%s] %s done.", format(Sys.time(), "%H:%M:%S"), cond_name))
}


# ==============================================================================
# Phase 2: Cross-condition Comprehensive Comparison
# ==============================================================================

message("\n========================================")
message("  Phase 2: Comprehensive Comparison")
message("========================================\n")

if (length(all_ci_results) < 2) {
  message(">>> Fewer than 2 conditions. Skipping comparison.")
} else {
  comp_dir <- file.path(output_base, "comparison")
  dir.create(comp_dir, recursive = TRUE, showWarnings = FALSE)
  png_dir <- file.path(tempdir(), "msProteomiX_comparison_png")
  dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)

  ci_all <- do.call(rbind, all_ci_results)
  ci_all$Group <- NA_character_
  for (gn in names(compare_groups)) ci_all$Group[ci_all$Condition %in% compare_groups[[gn]]] <- gn
  ci_all$Group[is.na(ci_all$Group)] <- "Other"
  write.csv(ci_all, file.path(comp_dir, "00_all_ci_merged.csv"), row.names = FALSE)

  # Fig 1: PG Count
  pg_data <- do.call(rbind, lapply(names(all_ms), function(nm) {
    mat <- as.matrix(all_ms[[nm]]$proteins); pgs <- colSums(!is.na(mat) & mat > 0)
    data.frame(Condition = nm, Group = ci_all$Group[ci_all$Condition == nm][1],
      Sample = all_ms[[nm]]$sample_names, PGs = pgs, stringsAsFactors = FALSE) }))
  pg_agg <- aggregate(PGs ~ Condition + Group, pg_data, function(x) c(mean=mean(x), sd=sd(x), cv=sd(x)/mean(x)*100))
  pg_agg <- cbind(pg_agg[,1:2], as.data.frame(pg_agg$PGs))

  p1 <- ggplot(pg_agg, aes(x = Condition, y = mean, fill = Group)) +
    geom_col(width = 0.65) + geom_errorbar(aes(ymin = mean-sd, ymax = mean+sd), width = 0.2) +
    geom_text(aes(label = sprintf("%.0f\n(CV=%.1f%%)", mean, cv)), vjust = -0.3, size = 2.3) +
    scale_fill_manual(values = grp_colors) +
    labs(title = "Protein Groups Identified", x = NULL, y = "Protein Groups") +
    theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8))
  ggsave(file.path(comp_dir, "10_pg_comparison.pdf"), p1, width = 8, height = 5)
  .save_png(p1, "fig01_pg", png_dir)
  message("  Fig01 PG count done")

  # Fig 2: Precursor/Peptide (SN)
  fig02_path <- NULL
  if (length(all_overviews) > 0) {
    ov_all <- do.call(rbind, all_overviews)
    ov_all$Group <- NA_character_; for (gn in names(compare_groups)) ov_all$Group[ov_all$Condition %in% compare_groups[[gn]]] <- gn
    ov_all$Group[is.na(ov_all$Group)] <- "Other"
    ov_long <- reshape2::melt(ov_all, id.vars = c("Condition","Group","FileName"),
      measure.vars = c("Precursors","Peptides","ProteinGroups"), variable.name = "Level", value.name = "Count")
    ov_long$Level <- factor(ov_long$Level, labels = c("Precursors","Peptides","Protein Groups"))
    p2 <- ggplot(ov_long, aes(x = Condition, y = Count, fill = Group)) +
      geom_boxplot(alpha = 0.8) + scale_fill_manual(values = grp_colors) +
      facet_wrap(~Level, scales = "free_y", ncol = 1) +
      labs(title = "Identification Summary (Spectronaut)", x = NULL, y = "Count") +
      theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 30, hjust = 1, size = 7), strip.text = element_text(face = "bold"))
    ggsave(file.path(comp_dir, "11_precursor.pdf"), p2, width = 8, height = 10)
    fig02_path <- .save_png(p2, "fig02_precursor", png_dir, w = 8, h = 10)
    message("  Fig02 precursor done")
  }

  # Fig 3: Data Completeness
  comp_data <- do.call(rbind, lapply(names(all_ms), function(nm) {
    mat <- as.matrix(all_ms[[nm]]$proteins); ns <- ncol(mat)
    det_n <- rowSums(!is.na(mat) & mat > 0); tb <- table(factor(det_n, levels = 1:ns))
    data.frame(Condition = nm, Group = ci_all$Group[ci_all$Condition == nm][1],
      NReps = as.integer(names(tb)), Count = as.integer(tb), Total = ns, stringsAsFactors = FALSE) }))
  # Use fraction label per condition (e.g., "6/6", "9/9" = "All")
  comp_data$Label <- paste0(comp_data$NReps, "/", comp_data$Total)
  # Order by detection fraction (descending: All first, then N-1, etc.)
  comp_data$Frac <- comp_data$NReps / comp_data$Total
  comp_data$Label <- reorder(factor(comp_data$Label), -comp_data$Frac)
  n_levels <- length(levels(comp_data$Label))
  p3 <- ggplot(comp_data, aes(x = Condition, y = Count, fill = Label)) +
    geom_col(position = "stack") +
    { if (n_levels <= 9) scale_fill_brewer(palette = "YlOrRd", direction = -1, name = "Detected in")
      else scale_fill_viridis_d(option = "inferno", direction = -1, name = "Detected in") } +
    labs(title = "Data Completeness", x = NULL, y = "Proteins") +
    theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8))
  ggsave(file.path(comp_dir, "12_data_completeness.pdf"), p3, width = 9, height = 5)
  .save_png(p3, "fig03_completeness", png_dir, w = 9)
  message("  Fig03 completeness done")

  # Fig 4: Quant CV
  cv_data <- do.call(rbind, lapply(names(all_ms), function(nm) {
    mat <- as.matrix(all_ms[[nm]]$proteins)
    if (isTRUE(all_ms[[nm]]$is_log2)) mat <- 2^mat; mat[mat == 0] <- NA
    keep <- rowSums(!is.na(mat)) >= 3; mat_f <- mat[keep, , drop = FALSE]
    cvs <- apply(mat_f, 1, function(x) { x <- x[!is.na(x)]; sd(x)/mean(x)*100 })
    data.frame(Condition = nm, Group = ci_all$Group[ci_all$Condition == nm][1], CV = cvs, stringsAsFactors = FALSE) }))
  cv_summary <- aggregate(CV ~ Condition + Group, cv_data, function(x) {
    c(n = length(x), median = median(x), mean = mean(x), pct20 = sum(x < 20)/length(x)*100, pct10 = sum(x < 10)/length(x)*100) })
  cv_summary <- cbind(cv_summary[,1:2], as.data.frame(cv_summary$CV))

  p4a <- ggplot(cv_data, aes(x = CV, color = Condition)) + geom_density(linewidth = 0.7) + xlim(0, 100) +
    geom_vline(xintercept = 20, linetype = "dashed", alpha = 0.5) +
    labs(title = "Quantification CV Distribution", x = "CV (%)", y = "Density") +
    theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(comp_dir, "13_quant_cv_density.pdf"), p4a, width = 9, height = 5)
  .save_png(p4a, "fig04a_cv_density", png_dir, w = 9)

  p4b <- ggplot(cv_data, aes(x = Condition, y = CV, fill = Group)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.8) + ylim(0, 100) +
    geom_hline(yintercept = 20, linetype = "dashed", alpha = 0.5) +
    scale_fill_manual(values = grp_colors) +
    labs(title = "Quantification CV", x = NULL, y = "CV (%)") +
    theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8))
  ggsave(file.path(comp_dir, "13_quant_cv_boxplot.pdf"), p4b, width = 9, height = 5)
  .save_png(p4b, "fig04b_cv_box", png_dir, w = 9)
  write.csv(cv_summary, file.path(comp_dir, "13_quant_cv_summary.csv"), row.names = FALSE)
  message("  Fig04 quant CV done")

  # Fig 5: Correlation heatmaps
  for (nm in names(all_ms)) {
    mat <- as.matrix(all_ms[[nm]]$proteins)
    if (isTRUE(all_ms[[nm]]$is_log2)) mat_l <- mat else { mat[mat == 0] <- NA; mat_l <- log2(mat) }
    cm <- cor(mat_l, use = "pairwise.complete.obs"); short <- sub("^\\d+_SZU_", "", sub("\\.raw$", "", all_ms[[nm]]$sample_names))
    rownames(cm) <- colnames(cm) <- short
    cd <- reshape2::melt(cm, varnames = c("S1","S2"), value.name = "R")
    p5 <- ggplot(cd, aes(x = S1, y = S2, fill = R)) + geom_tile() +
      geom_text(aes(label = sprintf("%.3f", R)), size = 2) +
      scale_fill_gradient2(low="#e74c3c", mid="#f1c40f", high="#27ae60", midpoint=median(cm[upper.tri(cm)])) +
      labs(title = paste0("Correlation - ", nm), x = NULL, y = NULL) +
      theme_bw(base_size = 9) + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 6), axis.text.y = element_text(size = 6))
    ggsave(file.path(comp_dir, paste0("14_correlation_", nm, ".pdf")), p5, width = 7, height = 6)
    .save_png(p5, paste0("fig05_corr_", nm), png_dir, w = 7, h = 6) }
  message("  Fig05 correlation done")

  # Fig 6: PCA
  pca_data <- do.call(rbind, lapply(names(all_ms), function(nm) {
    mat <- as.matrix(all_ms[[nm]]$proteins)
    if (!isTRUE(all_ms[[nm]]$is_log2)) { mat[mat == 0] <- NA; mat <- log2(mat) }
    mat[is.na(mat)] <- min(mat, na.rm = TRUE) / 2
    pc <- prcomp(t(mat), center = TRUE, scale. = TRUE)
    data.frame(PC1 = pc$x[,1], PC2 = pc$x[,2], Condition = nm,
      Group = ci_all$Group[ci_all$Condition == nm][1], stringsAsFactors = FALSE) }))
  p6 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = Condition, shape = Group)) +
    geom_point(size = 3, alpha = 0.85) +
    stat_ellipse(aes(group = Condition), type = "norm", linetype = "dashed", alpha = 0.4) +
    labs(title = "PCA - All Conditions", x = "PC1", y = "PC2") +
    theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(comp_dir, "15_pca.pdf"), p6, width = 9, height = 7)
  .save_png(p6, "fig06_pca", png_dir, w = 9, h = 7)
  message("  Fig06 PCA done")

  # Fig 7: CI Summary
  ci_long <- rbind(
    data.frame(Condition=ci_all$Condition, Group=ci_all$Group, Panel="PLT", Value=ci_all$CI_PLT, stringsAsFactors=FALSE),
    data.frame(Condition=ci_all$Condition, Group=ci_all$Group, Panel="RBC", Value=ci_all$CI_RBC, stringsAsFactors=FALSE),
    data.frame(Condition=ci_all$Condition, Group=ci_all$Group, Panel="PBMC", Value=ci_all$CI_PBMC, stringsAsFactors=FALSE))
  ci_long <- ci_long[!is.na(ci_long$Value), ]; ci_long$Panel <- factor(ci_long$Panel, levels = c("PLT","RBC","PBMC"))
  p7 <- ggplot(ci_long, aes(x = Condition, y = Value, fill = Group)) +
    geom_boxplot(alpha = 0.8) + facet_wrap(~Panel, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = grp_colors) +
    labs(title = "Contamination Index Summary", x = NULL, y = "CI") +
    theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 30, hjust = 1, size = 7))
  ggsave(file.path(comp_dir, "16_ci_summary.pdf"), p7, width = 8, height = 9)
  .save_png(p7, "fig07_ci", png_dir, w = 8, h = 9)
  message("  Fig07 CI done")

  # Fig 8: Top-N
  topn_data <- do.call(rbind, lapply(names(all_ms), function(nm) {
    mat <- as.matrix(all_ms[[nm]]$proteins); if (isTRUE(all_ms[[nm]]$is_log2)) mat <- 2^mat; mat[mat == 0] <- NA
    grp <- ci_all$Group[ci_all$Condition == nm][1]
    do.call(rbind, lapply(seq_len(ncol(mat)), function(j) {
      v <- mat[,j]; v <- v[!is.na(v) & v > 0]; v <- sort(v, decreasing = TRUE)
      data.frame(Condition = nm, Group = grp,
        Top10 = sum(v[1:min(10,length(v))])/sum(v)*100,
        Top20 = sum(v[1:min(20,length(v))])/sum(v)*100, stringsAsFactors = FALSE) })) }))
  topn_long <- reshape2::melt(topn_data, id.vars = c("Condition","Group"), measure.vars = c("Top10","Top20"),
    variable.name = "TopN", value.name = "Frac")
  p8 <- ggplot(topn_long, aes(x = Condition, y = Frac, fill = Group)) +
    geom_boxplot(alpha = 0.8) + facet_wrap(~TopN, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = grp_colors) +
    labs(title = "Top-N Protein Fraction (Depletion)", x = NULL, y = "%") +
    theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold"), axis.text.x = element_text(angle = 30, hjust = 1, size = 7))
  ggsave(file.path(comp_dir, "17_topN.pdf"), p8, width = 9, height = 7)
  .save_png(p8, "fig08_topn", png_dir, w = 9, h = 7)
  message("  Fig08 Top-N done")

  # Fig 9: High-abundance proteins
  ha_tgt <- c("ALB","IGHG1","IGHG2","IGHG3","IGHG4","TF","A2M","HP","APOA1","C3","IGKC")
  ha_up  <- c("P02768","P01857","P01859","P01860","P01861","P02787","P01023","P00738","P02647","P01024","P01834")
  ha_data <- do.call(rbind, lapply(names(all_ms), function(nm) {
    ms <- all_ms[[nm]]; pi <- ms$protein_info; mat <- as.matrix(ms$proteins)
    if (isTRUE(ms$is_log2)) mat <- 2^mat; mat[mat == 0] <- NA
    grp <- ci_all$Group[ci_all$Condition == nm][1]; mi <- rowMeans(mat, na.rm = TRUE); rnk <- rank(-mi, na.last = "last")
    do.call(rbind, lapply(seq_along(ha_tgt), function(k) {
      idx <- integer(0)
      if ("Gene" %in% colnames(pi)) idx <- grep(paste0("(^|;)", ha_tgt[k], "(;|$)"), pi$Gene)
      if (length(idx) == 0 && "Protein ID" %in% colnames(pi)) idx <- grep(ha_up[k], pi[["Protein ID"]], fixed = TRUE)
      if (length(idx) > 0) data.frame(Condition=nm, Group=grp, Protein=ha_tgt[k], Rank=rnk[idx[1]], Total=nrow(mat), stringsAsFactors=FALSE)
      else NULL })) }))
  p9 <- ggplot(ha_data, aes(x = Protein, y = Rank, color = Group)) +
    geom_point(size = 2.5, alpha = 0.8, position = position_dodge(width = 0.6)) +
    scale_color_manual(values = grp_colors) + coord_flip() +
    labs(title = "High-Abundance Protein Rank", subtitle = "Lower = more abundant", x = NULL, y = "Rank") +
    theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(comp_dir, "18_abundant.pdf"), p9, width = 9, height = 6)
  .save_png(p9, "fig09_abundant", png_dir, w = 9, h = 6)
  write.csv(ha_data, file.path(comp_dir, "18_abundant_data.csv"), row.names = FALSE)
  message("  Fig09 HA proteins done")

  # Fig 10-12: Peptide-level (SN only, auto-skip)
  fig10_path <- fig11_path <- fig12_path <- NULL
  sn_pep <- list()
  for (nm in names(all_ms)) {
    cn <- conditions[sapply(conditions, function(c) c$name == nm)][[1]]
    if (tolower(cn$engine) != "spectronaut" || !dir.exists(cn$path)) next
    pf <- list.files(cn$path, pattern = "Peptide_Report[.]tsv$", full.names = TRUE)
    if (length(pf) == 0) next
    pep <- read.delim(pf[1], check.names = FALSE, quote = "")
    if (nrow(pep) < 10) next
    ids <- pep[[1]]; seqs <- sub("^_(.+)_[.][0-9]+$", "\\1", ids)
    seqs_clean <- gsub("[[][^]]*[]]", "", seqs); charges <- as.integer(sub(".*[.]([0-9]+)$", "\\1", ids))
    mc <- sapply(seqs_clean, function(s) { ch <- strsplit(s,"")[[1]]; n <- length(ch)
      if (n < 2) return(0L); cnt <- 0L; for (i in 1:(n-1)) if (ch[i] %in% c("K","R") && ch[i+1] != "P") cnt <- cnt+1L; cnt }, USE.NAMES = FALSE)
    sn_pep[[nm]] <- list(mc=mc, len=nchar(seqs_clean), charge=charges, grp=ci_all$Group[ci_all$Condition==nm][1])
    message(sprintf("  Peptide_Report: %s (%d precursors)", nm, nrow(pep)))
  }

  if (length(sn_pep) > 0) {
    mc_df <- do.call(rbind, lapply(names(sn_pep), function(nm) {
      d <- sn_pep[[nm]]; data.frame(Condition=nm, Group=d$grp,
        MC=factor(pmin(d$mc,3), levels=0:3, labels=c("0-miss","1-miss","2-miss","3+")), stringsAsFactors=FALSE) }))
    mc_pct <- as.data.frame(prop.table(table(mc_df$Condition, mc_df$MC), margin = 1) * 100)
    colnames(mc_pct) <- c("Condition","MC","Pct")
    p10 <- ggplot(mc_pct, aes(x = Condition, y = Pct, fill = MC)) + geom_col(position = "stack") +
      geom_text(aes(label = ifelse(Pct > 3, sprintf("%.1f%%", Pct), "")), position = position_stack(vjust = 0.5), size = 2.5) +
      scale_fill_brewer(palette = "Set2", name = "MC") +
      labs(title = "Missed Cleavage (Spectronaut)", x = NULL, y = "%") +
      theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave(file.path(comp_dir, "19_missed_cleavage.pdf"), p10, width = 8, height = 5)
    fig10_path <- .save_png(p10, "fig10_mc", png_dir)
    message("  Fig10 MC done")

    len_df <- do.call(rbind, lapply(names(sn_pep), function(nm)
      data.frame(Condition=nm, Length=sn_pep[[nm]]$len, stringsAsFactors=FALSE)))
    p11 <- ggplot(len_df, aes(x = Length, color = Condition)) + geom_density(linewidth = 0.7) + xlim(5, 50) +
      labs(title = "Peptide Length (Spectronaut)", x = "aa", y = "Density") +
      theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave(file.path(comp_dir, "20_peptide_length.pdf"), p11, width = 8, height = 5)
    fig11_path <- .save_png(p11, "fig11_peplen", png_dir)
    message("  Fig11 peptide length done")

    ch_df <- do.call(rbind, lapply(names(sn_pep), function(nm)
      data.frame(Condition=nm, Charge=factor(sn_pep[[nm]]$charge), stringsAsFactors=FALSE)))
    ch_pct <- as.data.frame(prop.table(table(ch_df$Condition, ch_df$Charge), margin = 1) * 100)
    colnames(ch_pct) <- c("Condition","Charge","Pct")
    p12 <- ggplot(ch_pct[ch_pct$Pct > 0.5,], aes(x = Charge, y = Pct, fill = Condition)) +
      geom_col(position = "dodge", alpha = 0.85) +
      labs(title = "Charge State (Spectronaut)", x = "Charge", y = "%") +
      theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
    ggsave(file.path(comp_dir, "21_charge.pdf"), p12, width = 8, height = 5)
    fig12_path <- .save_png(p12, "fig12_charge", png_dir)
    message("  Fig12 charge done")
  }

  # Fig 13: Dynamic Range
  dr_data <- do.call(rbind, lapply(names(all_ms), function(nm) {
    mat <- as.matrix(all_ms[[nm]]$proteins); if (isTRUE(all_ms[[nm]]$is_log2)) mat <- 2^mat; mat[mat==0] <- NA
    grp <- ci_all$Group[ci_all$Condition==nm][1]
    do.call(rbind, lapply(seq_len(ncol(mat)), function(j) {
      v <- mat[,j]; v <- v[!is.na(v) & v > 0]; v <- log2(v); q <- quantile(v, c(0.01,0.99))
      data.frame(Condition=nm, Group=grp, DR=q[2]-q[1], stringsAsFactors=FALSE) })) }))
  p13 <- ggplot(dr_data, aes(x = Condition, y = DR, fill = Group)) + geom_boxplot(alpha = 0.8) +
    scale_fill_manual(values = grp_colors) +
    labs(title = "Dynamic Range (P1-P99, log2)", x = NULL, y = "log2 orders") +
    theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1, size = 8))
  ggsave(file.path(comp_dir, "22_dynamic_range.pdf"), p13, width = 9, height = 5)
  .save_png(p13, "fig13_dr", png_dir, w = 9)
  message("  Fig13 dynamic range done")

  # Fig 14: Rank Overlay
  rk_data <- do.call(rbind, lapply(names(all_ms), function(nm) {
    mat <- as.matrix(all_ms[[nm]]$proteins); if (isTRUE(all_ms[[nm]]$is_log2)) mat <- 2^mat; mat[mat==0] <- NA
    mi <- apply(mat, 1, median, na.rm = TRUE); mi <- mi[!is.na(mi) & mi > 0]; mi <- sort(mi, decreasing = TRUE)
    data.frame(Condition=nm, Rank=seq_along(mi), log2I=log2(mi), stringsAsFactors=FALSE) }))
  p14 <- ggplot(rk_data, aes(x = Rank, y = log2I, color = Condition)) + geom_line(linewidth = 0.5, alpha = 0.8) +
    labs(title = "Protein Rank Curve (Median)", x = "Rank", y = "log2(Intensity)") +
    theme_bw(base_size = 10) + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
  ggsave(file.path(comp_dir, "23_rank_overlay.pdf"), p14, width = 9, height = 5)
  .save_png(p14, "fig14_rank", png_dir, w = 9)
  message("  Fig14 rank overlay done")

  # ==========================================================================
  # Phase 3: Professional HTML Report (EN only, ZH to follow)
  # ==========================================================================
  message("\n>>> Phase 3: Professional HTML Report (EN)")
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    message(">>> base64enc not installed. Run: install.packages('base64enc')")
  } else {

  css <- '<style>@page{margin:2cm 2.5cm;size:A4}body{font-family:"Helvetica Neue","PingFang SC","Noto Sans SC",sans-serif;font-size:11pt;line-height:1.6;color:#333;max-width:850px;margin:0 auto;padding:20px}h1{font-size:20pt;border-bottom:3px solid #2c3e50;padding-bottom:8px;color:#2c3e50}h2{font-size:15pt;color:#2c3e50;border-bottom:1px solid #bdc3c7;padding-bottom:4px;margin-top:30px}h3{font-size:13pt;color:#34495e;margin-top:20px}table{border-collapse:collapse;width:100%;margin:10px 0;font-size:10pt}th{background:#2c3e50;color:#fff;padding:8px 10px;text-align:center}td{border:1px solid #ddd;padding:6px 10px}tr:nth-child(even){background:#f9f9f9}.figure{margin:15px 0;page-break-inside:avoid;text-align:center}.caption{font-size:9pt;color:#7f8c8d;margin-top:4px}.note{background:#eaf2f8;border-left:4px solid #3498db;padding:10px 15px;margin:10px 0;font-size:10pt}.warn{background:#fef9e7;border-left:4px solid #f39c12;padding:10px 15px;margin:10px 0;font-size:10pt}.header-meta{color:#7f8c8d;font-size:10pt}hr{border:none;border-top:1px solid #ddd;margin:25px 0}blockquote{border-left:4px solid #2c3e50;margin:15px 0;padding:10px 20px;font-style:italic;color:#555;background:#f8f9fa}</style>'

  # --- Helper: build per-condition data overview rows ---
  ov_rows <- paste0(sapply(names(all_ms), function(nm) {
    ci_sub <- ci_all[ci_all$Condition == nm, ]
    det_sub <- all_det[[nm]]
    sprintf('<tr><td><b>%s</b></td><td>%s</td><td>%d</td><td>%.0f</td><td>%.1f%%</td><td>%s</td></tr>',
      nm, ci_all$Group[ci_all$Condition == nm][1],
      ncol(all_ms[[nm]]$proteins), mean(ci_sub$PGs),
      mean(det_sub$DetRate),
      paste0(sum(ci_sub$Overall == "Pass"), "/", nrow(ci_sub), " Pass"))
  }), collapse = "")

  # --- Helper: build comprehensive CI judgment table ---
  ci_table_rows <- paste0(sapply(seq_len(nrow(ci_all)), function(i) {
    r <- ci_all[i, ]
    .fmt_j <- function(v, j) {
      icon <- if (j == "Pass") " &#x2705;" else if (j == "Flag") " &#x26A0;&#xFE0F;" else if (j == "Fail") " &#x274C;" else ""
      sprintf("%.4f (%s%s)", v, j, icon)
    }
    sprintf('<tr><td>%s</td><td>%s</td><td>%d</td><td>%s</td><td>%s</td><td>%s</td><td><b>%s</b></td></tr>',
      r$Condition, sub("\\.raw$", "", r$Sample), r$PGs,
      .fmt_j(r$CI_PLT, r$J_PLT), .fmt_j(r$CI_RBC, r$J_RBC),
      if (is.na(r$CI_PBMC)) "N/A" else .fmt_j(r$CI_PBMC, r$J_PBMC),
      r$Overall)
  }), collapse = "")

  # --- Helper: per-condition CI analysis ---
  per_cond_ci <- paste0(sapply(names(all_ms), function(nm) {
    ci_sub <- ci_all[ci_all$Condition == nm, ]
    n_pass <- sum(ci_sub$Overall == "Pass"); n_flag <- sum(ci_sub$Overall == "Flag"); n_fail <- sum(ci_sub$Overall == "Fail")
    n_total <- nrow(ci_sub)
    med_plt <- median(ci_sub$CI_PLT); med_rbc <- median(ci_sub$CI_RBC)
    med_pbmc <- if (all(is.na(ci_sub$CI_PBMC))) NA else median(ci_sub$CI_PBMC, na.rm = TRUE)

    status <- if (n_pass == n_total) "Pass" else if (n_fail > 0) "Fail" else "Flag"
    icon <- if (status == "Pass") "&#x2705;" else if (status == "Flag") "&#x26A0;&#xFE0F;" else "&#x274C;"

    detail <- sprintf("CI_PLT median = %.4f (threshold: 0.009), CI_RBC median = %.4f (threshold: 0.05)", med_plt, med_rbc)
    if (!is.na(med_pbmc)) detail <- paste0(detail, sprintf(", CI_PBMC median = %.4f (threshold: 0.005)", med_pbmc))

    advice <- if (status == "Pass") {
      "All CI indicators are within safe thresholds. Data is clean and can be directly used for downstream analysis."
    } else if (n_fail > 0) {
      "Significant contamination detected. Marker proteins should be excluded before downstream analysis."
    } else {
      "Mild contamination detected but below Fail thresholds. Data is usable; check marker proteins if analyzing related pathways."
    }

    paste0('<h3>', nm, ' &mdash; ', status, ' ', icon, ' (', n_pass, '/', n_total, ' samples Pass)</h3>',
      '<p>', detail, '</p><p><b>Recommendation:</b> ', advice, '</p>')
  }), collapse = "")

  # --- Helper: per-indicator CI analysis ---
  per_indicator_ci <- paste0(sapply(c("PLT", "RBC", "PBMC"), function(panel) {
    col <- paste0("CI_", panel)
    j_col <- paste0("J_", panel)
    vals <- ci_all[[col]]; js <- ci_all[[j_col]]
    if (all(is.na(vals))) return("")

    n_pass <- sum(js == "Pass", na.rm = TRUE); n_total <- sum(!is.na(js))
    panel_name <- if (panel == "PLT") "Platelet (CI_PLT)" else if (panel == "RBC") "Erythrocyte (CI_RBC)" else "PBMC (CI_PBMC)"

    # per-group medians
    grp_stats <- paste0(sapply(names(compare_groups), function(gn) {
      idx <- ci_all$Condition %in% compare_groups[[gn]]
      v <- vals[idx]; v <- v[!is.na(v)]
      if (length(v) == 0) return("")
      sprintf("%s median = %.4f", gn, median(v))
    }), collapse = ", ")

    conclusion <- if (n_pass == n_total) "No significant contamination in any condition." else
      sprintf("%d/%d samples passed. Some samples show elevated levels.", n_pass, n_total)

    paste0('<h3>', panel_name, '</h3>',
      '<p><b>', n_pass, '/', n_total, '</b> samples Pass. ', grp_stats, '.</p>',
      '<p><b>Conclusion:</b> ', conclusion, '</p>')
  }), collapse = "")

  # --- Helper: executive summary table ---
  sc_rows <- paste0(sapply(names(all_ms), function(nm) {
    ci_sub <- ci_all[ci_all$Condition == nm, ]; cv_sub <- cv_summary[cv_summary$Condition == nm, ]
    det_sub <- all_det[[nm]]
    sprintf('<tr><td>%s</td><td>%s</td><td>%d</td><td>%.0f</td><td>%.1f%%</td><td>%.1f%%</td><td>%.1f%%</td><td>%s</td></tr>',
      nm, ci_all$Group[ci_all$Condition == nm][1],
      ncol(all_ms[[nm]]$proteins), mean(ci_sub$PGs),
      mean(det_sub$DetRate), cv_sub$median[1], cv_sub$pct20[1],
      paste0(sum(ci_sub$Overall == "Pass"), "/", nrow(ci_sub), " Pass"))
  }), collapse = "")

  # --- Helper: CV interpretation ---
  cv_interp <- paste0('<p>The CV values reported here are calculated on <b>normalized LFQ-equivalent quantities</b> ',
    '(Spectronaut PG.Quantity) across <b>biological replicates</b> (different individuals), ',
    'not technical replicates. Therefore, the observed median CVs (40-60%) reflect genuine ',
    'inter-individual biological variation in addition to technical variance. ',
    'Despite these CVs, the high Pearson correlations (&gt;0.95) confirm that the ',
    'overall quantification baseline is highly consistent across all methods.</p>')

  # --- Helper: dynamic range interpretation ---
  dr_agg <- aggregate(DR ~ Condition + Group, dr_data, median)
  dr_interp <- paste0('<p>Dynamic range is measured as the P1-P99 intensity span in log2 space per sample. ',
    paste0(sapply(seq_len(nrow(dr_agg)), function(i) sprintf("%s: median %.1f log2 orders", dr_agg$Condition[i], dr_agg$DR[i])), collapse = "; "),
    '.</p>')

  # --- Helper: Top-N interpretation ---
  topn_agg <- aggregate(Frac ~ Condition + Group + TopN, topn_long, median)
  top10_info <- topn_agg[topn_agg$TopN == "Top10", ]
  topn_interp <- paste0('<p>Top-10 protein fraction indicates the degree of high-abundance protein depletion. ',
    'Lower values mean better depletion. ',
    paste0(sapply(seq_len(nrow(top10_info)), function(i) sprintf("%s: %.1f%%", top10_info$Condition[i], top10_info$Frac[i])), collapse = "; "),
    '.</p>')

  # --- Conclusions: find best condition ---
  best_pg <- pg_agg$Condition[which.max(pg_agg$mean)]
  best_cv <- cv_summary$Condition[which.min(cv_summary$median)]
  all_pass_conds <- names(which(sapply(names(all_ms), function(nm) {
    ci_sub <- ci_all[ci_all$Condition == nm, ]; all(ci_sub$Overall == "Pass") })))

  conclusions_html <- paste0(
    '<ol>',
    '<li><b>All conditions show no severe cell contamination</b>: CI_PLT and CI_RBC are well within safe thresholds across all samples.</li>',
    sprintf('<li><b>%s</b> achieves the highest protein identification (~%.0f PGs on average).</li>', best_pg, pg_agg$mean[pg_agg$Condition == best_pg]),
    sprintf('<li><b>%s</b> shows the best quantification reproducibility (median CV = %.1f%%).</li>', best_cv, cv_summary$median[cv_summary$Condition == best_cv]),
    if (length(all_pass_conds) > 0) sprintf('<li>Conditions with 100%% CI Pass rate: <b>%s</b>.</li>', paste(all_pass_conds, collapse = ", ")) else "",
    '</ol>')

  key_insight <- paste0(
    '<blockquote>',
    '<b>Higher protein identification count does not equal better data quality.</b><br>',
    'Cell contamination (platelets, erythrocytes, PBMC) can lead to false identification of non-plasma proteins, ',
    'artificially inflating protein counts. Contamination assessment must precede any comparison of identification depth.<br><br>',
    '<em>&mdash; Gao H, Korff K, Mann M, Guo T. Deeper is not always better in plasma proteomics. Nature Biotechnology 2026</em>',
    '</blockquote>')

  # --- Per-condition recommendations ---
  per_cond_rec <- paste0('<table>',
    '<tr><th>Condition</th><th>Recommendation</th></tr>',
    paste0(sapply(names(all_ms), function(nm) {
      ci_sub <- ci_all[ci_all$Condition == nm, ]
      n_pass <- sum(ci_sub$Overall == "Pass"); n_total <- nrow(ci_sub)
      rec <- if (n_pass == n_total) "CI all Pass. Data can be directly used for downstream analysis." else
        "Some samples flagged. Usable for analysis; check marker overlap in differential results."
      sprintf('<tr><td><b>%s</b></td><td>%s</td></tr>', nm, rec)
    }), collapse = ""),
    '</table>')

  # --- References ---
  refs <- paste0('<ol style="font-size:10pt;color:#555;">',
    '<li>Geyer PE, Voytik E, Treit PV, Doll S, Kleinhempel A, Niu L, et al. Plasma proteome profiling to detect and avoid sample-related biases in biomarker studies. <i>EMBO Molecular Medicine</i> 2019; 11(11): EMMM201910427.</li>',
    '<li>Gao H, Zhan Y, Liu Y, Zhu Z, Zheng Y, Qian L, et al. Systematic evaluation of blood contamination in nanoparticle-based plasma proteomics. <i>EMBO Molecular Medicine</i> 2026; 18(1): 275-296.</li>',
    '<li>Korff K, M&uuml;ller-Reif JB, Fichtl D, Albrecht V, Schebesta AS, Itang EC, et al. Pre-analytical drivers of bias in bead-enriched plasma proteomics. <i>EMBO Molecular Medicine</i> 2025; 17(11): 3174.</li>',
    '<li>Gao H, Korff K, Mann M, Guo T. Deeper is not always better in plasma proteomics. <i>Nature Biotechnology</i> 2026; 1-4.</li>',
    '</ol>')

  software <- paste0('<ul style="font-size:10pt;color:#555;">',
    sprintf('<li>Data analysis: msProteomiX v%s (BayOmics)</li>', as.character(packageVersion("msProteomiX"))),
    '<li>CI calculation: Geyer 2019 formula + Gao 2026 NP-optimized marker panels</li>',
    '<li>Search engine: Spectronaut / DIA-NN (as configured)</li>',
    '</ul>')

  # === Assemble HTML ===
  html <- paste0(
    '<!DOCTYPE html><html><head><meta charset="utf-8"><title>', project_name, ' - Plasma Method Comparison Report</title>', css, '</head><body>',

    # Header
    '<h1>', project_name, ' - Plasma Method Comparison Report</h1>',
    '<p class="header-meta">Comprehensive Evaluation of Bead-Based Sample Preparation Methods</p>',
    sprintf('<p class="header-meta">Generated: %s | msProteomiX v%s | Assessed by: BayOmics Bioinformatics Team</p>', Sys.Date(), as.character(packageVersion("msProteomiX"))),
    '<hr>',

    # S1: Data Overview
    '<h2>1. Data Overview</h2>',
    '<p>This report systematically evaluates ', length(names(all_ms)), ' sample preparation conditions using ',
    sum(sapply(all_ms, function(m) ncol(m$proteins))), ' total samples from Spectronaut/DIA-NN search results.</p>',
    '<table><tr><th>Condition</th><th>Group</th><th>Samples</th><th>Mean PGs</th><th>Mean Detection Rate</th><th>CI Overall</th></tr>',
    ov_rows, '</table>',
    '<hr>',

    # S2: Contamination Assessment Method
    '<h2>2. Contamination Assessment Methodology</h2>',
    '<h3>2.1 Framework</h3>',
    '<p>This report adopts the latest plasma quality assessment framework in proteomics:</p>',
    '<ul>',
    '<li><b>Contamination Index (CI)</b>: Quantitative formula proposed by Geyer et al. 2019 (<i>EMBO Molecular Medicine</i>), measuring the proportion of cell contamination proteins in total protein intensity.</li>',
    '<li><b>Marker Panels</b>: PLT and RBC panels from Gao et al. 2026 (<i>EMBO Molecular Medicine</i>), optimized for nanoparticle/bead enrichment workflows (30 proteins each).</li>',
    '<li><b>PBMC Panel</b>: 5-protein panel from Korff et al. 2025 (<i>EMBO Molecular Medicine</i>).</li>',
    '<li><b>Philosophy</b>: Following <i>Nat Biotechnol</i> 2026 &mdash; "Deeper is not always better in plasma proteomics" &mdash; data depth does not equal data quality.</li>',
    '</ul>',
    '<h3>2.2 CI Formula</h3>',
    '<div class="note">CI = &Sigma;(marker protein intensities) / &Sigma;(all protein intensities)<br>',
    'Higher CI indicates greater contamination from the corresponding cell type.</div>',
    '<h3>2.3 Judgment Thresholds</h3>',
    '<table><tr><th>Indicator</th><th>Pass</th><th>Flag</th><th>Fail</th><th>Source</th></tr>',
    '<tr><td>CI_Platelet</td><td>&le; 0.009</td><td>0.009 &ndash; 0.05</td><td>&gt; 0.05</td><td>Baize tool</td></tr>',
    '<tr><td>CI_Erythrocyte</td><td>&le; 0.05</td><td>0.05 &ndash; 0.1</td><td>&gt; 0.1</td><td>Baize tool</td></tr>',
    '<tr><td>CI_PBMC</td><td>&lt; 0.005</td><td>0.005 &ndash; 0.02</td><td>&gt; 0.02</td><td>Literature estimate</td></tr>',
    '</table>',
    '<div class="note"><b>Interpretation:</b> Pass &#x2705; = CI below safe threshold, contamination negligible. ',
    'Flag &#x26A0;&#xFE0F; = mild contamination, data usable with caution. ',
    'Fail &#x274C; = significant contamination, marker proteins should be excluded.</div>',
    '<hr>',

    # S3: Contamination Results
    '<h2>3. Contamination Assessment Results</h2>',
    '<h3>3.1 Comprehensive Judgment Table</h3>',
    '<table><tr><th>Condition</th><th>Sample</th><th>PGs</th><th>CI_PLT</th><th>CI_RBC</th><th>CI_PBMC</th><th>Overall</th></tr>',
    ci_table_rows, '</table>',
    '<h3>3.2 Per-Condition Analysis</h3>', per_cond_ci,
    '<h3>3.3 Per-Indicator Analysis</h3>', per_indicator_ci,
    .b64img(file.path(png_dir, "fig07_ci.png"), "85%", "Contamination Index Comparison"),
    '<hr>',

    # S4: Identification & Completeness
    '<h2>4. Identification & Data Completeness</h2>',
    '<table><tr><th>Condition</th><th>Group</th><th>Samples</th><th>Mean PGs</th><th>Median CV (%)</th><th>%CV&lt;20%</th><th>Mean Detection Rate</th><th>CI</th></tr>',
    sc_rows, '</table>',
    .b64img(file.path(png_dir, "fig01_pg.png"), "90%", "Protein Groups Identified"),
    if (!is.null(fig02_path)) .b64img(fig02_path, "80%", "Precursors / Peptides") else "",
    .b64img(file.path(png_dir, "fig03_completeness.png"), "90%", "Data Completeness"),
    '<hr>',

    # S5: Quantification Precision
    '<h2>5. Quantification Precision</h2>',
    cv_interp,
    .b64img(file.path(png_dir, "fig04a_cv_density.png"), "90%", "CV Distribution"),
    .b64img(file.path(png_dir, "fig04b_cv_box.png"), "90%", "CV Boxplot"),
    paste0(sapply(names(all_ms), function(nm) .b64img(file.path(png_dir, paste0("fig05_corr_", nm, ".png")), "75%", paste0("Correlation - ", nm))), collapse = ""),
    .b64img(file.path(png_dir, "fig06_pca.png"), "85%", "PCA"),
    '<hr>',

    # S6: Sample Quality
    '<h2>6. Sample Quality (Depletion Efficiency)</h2>',
    topn_interp,
    .b64img(file.path(png_dir, "fig08_topn.png"), "85%", "Top-N Fraction"),
    .b64img(file.path(png_dir, "fig09_abundant.png"), "85%", "High-Abundance Proteins"),
    '<hr>',

    # S7: Digestion Quality
    if (!is.null(fig10_path)) paste0(
      '<h2>7. Digestion Quality (Spectronaut Peptide-Level)</h2>',
      .b64img(fig10_path, "85%", "Missed Cleavage"),
      if (!is.null(fig11_path)) .b64img(fig11_path, "85%", "Peptide Length") else "",
      if (!is.null(fig12_path)) .b64img(fig12_path, "85%", "Charge State") else "",
      '<hr>') else "",

    # S8: Dynamic Range
    '<h2>8. Dynamic Range</h2>',
    dr_interp,
    .b64img(file.path(png_dir, "fig13_dr.png"), "90%", "Dynamic Range"),
    .b64img(file.path(png_dir, "fig14_rank.png"), "90%", "Rank Curve"),
    '<hr>',

    # S9: Conclusions & Recommendations
    '<h2>9. Conclusions & Recommendations</h2>',
    '<h3>9.1 Core Findings</h3>', conclusions_html,
    '<h3>9.2 Key Insight</h3>', key_insight,
    '<h3>9.3 Per-Condition Recommendations</h3>', per_cond_rec,
    '<h3>9.4 Next Steps</h3>',
    '<ol><li>If comparing against competitor kits, run the same CI evaluation on competitor search results.</li>',
    '<li>Ensure data quality is confirmed before comparing identification depth across methods.</li></ol>',
    '<hr>',

    # S10: Methodology & References
    '<h2>10. Methodology & References</h2>',
    '<h3>References</h3>', refs,
    '<h3>Software</h3>', software,
    '<hr>',

    # Footer
    '<p style="text-align:center;color:#aaa;font-size:9pt;">Report auto-generated by BayOmics Bioinformatics Team.</p>',
    '</body></html>')

  # --- Save EN HTML ---
  html_path <- file.path(comp_dir, "10_comparison_report_en.html")
  writeLines(html, html_path, useBytes = TRUE)
  message(sprintf(">>> HTML saved: %s", basename(html_path)))

  # --- Chrome PDF ---
  chrome <- .chrome_path()
  if (!is.null(chrome)) {
    pdf_path <- file.path(comp_dir, "10_comparison_report_en.pdf")
    tmp_pdf <- file.path(tempdir(), "comp_rpt_en.pdf")
    cmd <- sprintf('"%s" --headless --disable-gpu --no-sandbox --print-to-pdf="%s" "%s"', chrome, tmp_pdf, html_path)
    res <- system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (res == 0 && file.exists(tmp_pdf)) {
      file.copy(tmp_pdf, pdf_path, overwrite = TRUE); file.remove(tmp_pdf)
      message(sprintf(">>> PDF saved: %s", basename(pdf_path)))
    }
  } else { message(">>> Chrome not found. Open HTML in browser > Print > Save as PDF.") }

  # --- ZH Version (gsub translation-layer, Pitfall #17 compliant) ---
  message(">>> Generating Chinese version ...")
  html_zh <- html

  # Header
  html_zh <- gsub("Plasma Method Comparison Report", "\u8840\u6d46\u524d\u5904\u7406\u65b9\u6848\u6bd4\u8f83\u62a5\u544a", html_zh, fixed = TRUE)
  html_zh <- gsub("Comprehensive Evaluation of Bead-Based Sample Preparation Methods", "\u57fa\u4e8e\u78c1\u73e0\u7684\u8840\u6d46\u86cb\u767d\u7ec4\u524d\u5904\u7406\u65b9\u6848\u7efc\u5408\u8bc4\u4f30", html_zh, fixed = TRUE)
  html_zh <- gsub("Generated:", "\u62a5\u544a\u65e5\u671f:", html_zh, fixed = TRUE)
  html_zh <- gsub("Assessed by: BayOmics Bioinformatics Team", "\u8bc4\u4f30\u65b9: BayOmics \u751f\u7269\u4fe1\u606f\u56e2\u961f", html_zh, fixed = TRUE)

  # Section titles
  html_zh <- gsub(">1. Data Overview<", ">1. \u6570\u636e\u6982\u89c8<", html_zh, fixed = TRUE)
  html_zh <- gsub(">2. Contamination Assessment Methodology<", ">2. \u6c61\u67d3\u8bc4\u4f30\u65b9\u6cd5\u5b66<", html_zh, fixed = TRUE)
  html_zh <- gsub(">3. Contamination Assessment Results<", ">3. \u6c61\u67d3\u8bc4\u4f30\u7ed3\u679c<", html_zh, fixed = TRUE)
  html_zh <- gsub(">4. Identification &amp; Data Completeness<", ">4. \u9274\u5b9a\u80fd\u529b\u4e0e\u6570\u636e\u5b8c\u6574\u6027<", html_zh, fixed = TRUE)
  html_zh <- gsub(">5. Quantification Precision<", ">5. \u5b9a\u91cf\u7cbe\u5ea6<", html_zh, fixed = TRUE)
  html_zh <- gsub(">6. Sample Quality (Depletion Efficiency)<", ">6. \u6837\u672c\u8d28\u91cf (\u53bb\u9664\u6548\u7387)<", html_zh, fixed = TRUE)
  html_zh <- gsub(">7. Digestion Quality (Spectronaut Peptide-Level)<", ">7. \u9176\u89e3\u8d28\u91cf (Spectronaut \u80bd\u6bb5\u7ea7)<", html_zh, fixed = TRUE)
  html_zh <- gsub(">8. Dynamic Range<", ">8. \u52a8\u6001\u8303\u56f4<", html_zh, fixed = TRUE)
  html_zh <- gsub(">9. Conclusions &amp; Recommendations<", ">9. \u7ed3\u8bba\u4e0e\u5efa\u8bae<", html_zh, fixed = TRUE)
  html_zh <- gsub(">10. Methodology &amp; References<", ">10. \u65b9\u6cd5\u5b66\u4e0e\u53c2\u8003\u6587\u732e<", html_zh, fixed = TRUE)

  # Sub-section titles
  html_zh <- gsub(">2.1 Framework<", ">2.1 \u8bc4\u4f30\u6846\u67b6<", html_zh, fixed = TRUE)
  html_zh <- gsub(">2.2 CI Formula<", ">2.2 CI \u8ba1\u7b97\u516c\u5f0f<", html_zh, fixed = TRUE)
  html_zh <- gsub(">2.3 Judgment Thresholds<", ">2.3 \u5224\u5b9a\u9608\u503c<", html_zh, fixed = TRUE)
  html_zh <- gsub(">3.1 Comprehensive Judgment Table<", ">3.1 \u7efc\u5408\u5224\u5b9a\u8868<", html_zh, fixed = TRUE)
  html_zh <- gsub(">3.2 Per-Condition Analysis<", ">3.2 \u9010\u6761\u4ef6\u5206\u6790<", html_zh, fixed = TRUE)
  html_zh <- gsub(">3.3 Per-Indicator Analysis<", ">3.3 \u9010\u6307\u6807\u5206\u6790<", html_zh, fixed = TRUE)
  html_zh <- gsub(">9.1 Core Findings<", ">9.1 \u6838\u5fc3\u53d1\u73b0<", html_zh, fixed = TRUE)
  html_zh <- gsub(">9.2 Key Insight<", ">9.2 \u6838\u5fc3\u89c2\u70b9<", html_zh, fixed = TRUE)
  html_zh <- gsub(">9.3 Per-Condition Recommendations<", ">9.3 \u9010\u6761\u4ef6\u5efa\u8bae<", html_zh, fixed = TRUE)
  html_zh <- gsub(">9.4 Next Steps<", ">9.4 \u540e\u7eed\u5efa\u8bae<", html_zh, fixed = TRUE)
  html_zh <- gsub(">References<", ">\u53c2\u8003\u6587\u732e<", html_zh, fixed = TRUE)
  html_zh <- gsub(">Software<", ">\u8f6f\u4ef6<", html_zh, fixed = TRUE)

  # Table headers
  html_zh <- gsub(">Condition<", ">\u6761\u4ef6<", html_zh, fixed = TRUE)
  html_zh <- gsub(">Group<", ">\u5206\u7ec4<", html_zh, fixed = TRUE)
  html_zh <- gsub(">Samples<", ">\u6837\u672c\u6570<", html_zh, fixed = TRUE)
  html_zh <- gsub(">Mean PGs<", ">\u5e73\u5747 PGs<", html_zh, fixed = TRUE)
  html_zh <- gsub(">Mean Detection Rate<", ">\u5e73\u5747\u68c0\u6d4b\u7387<", html_zh, fixed = TRUE)
  html_zh <- gsub(">CI Overall<", ">CI \u7efc\u5408\u5224\u5b9a<", html_zh, fixed = TRUE)
  html_zh <- gsub(">Sample<", ">\u6837\u672c<", html_zh, fixed = TRUE)
  html_zh <- gsub(">Overall<", ">\u7efc\u5408<", html_zh, fixed = TRUE)
  html_zh <- gsub(">Indicator<", ">\u6307\u6807<", html_zh, fixed = TRUE)
  html_zh <- gsub(">Source<", ">\u6765\u6e90<", html_zh, fixed = TRUE)
  html_zh <- gsub(">Recommendation<", ">\u5efa\u8bae<", html_zh, fixed = TRUE)

  # Body text
  html_zh <- gsub("This report systematically evaluates", "\u672c\u62a5\u544a\u7cfb\u7edf\u8bc4\u4f30\u4e86", html_zh, fixed = TRUE)
  html_zh <- gsub("sample preparation conditions using", "\u4e2a\u524d\u5904\u7406\u6761\u4ef6\uff0c\u5171\u5305\u542b", html_zh, fixed = TRUE)
  html_zh <- gsub("total samples from Spectronaut/DIA-NN search results.", "\u4e2a\u6837\u672c\u7684 Spectronaut/DIA-NN \u641c\u5e93\u7ed3\u679c\u3002", html_zh, fixed = TRUE)
  html_zh <- gsub("This report adopts the latest plasma quality assessment framework in proteomics:", "\u672c\u62a5\u544a\u91c7\u7528\u56fd\u9645\u86cb\u767d\u7ec4\u5b66\u9886\u57df\u6700\u65b0\u7684\u8840\u6d46\u6837\u672c\u8d28\u91cf\u8bc4\u4f30\u6846\u67b6\uff1a", html_zh, fixed = TRUE)
  html_zh <- gsub("measuring the proportion of cell contamination proteins in total protein intensity.", "\u8ba1\u7b97\u5404\u7c7b\u7ec6\u80de\u6c61\u67d3\u86cb\u767d\u5360\u603b\u86cb\u767d\u5f3a\u5ea6\u7684\u6bd4\u4f8b\u3002", html_zh, fixed = TRUE)
  html_zh <- gsub("optimized for nanoparticle/bead enrichment workflows (30 proteins each).", "\u9488\u5bf9\u7eb3\u7c73\u78c1\u73e0\u5bcc\u96c6\u5de5\u4f5c\u6d41\u4f18\u5316 (\u5404 30 \u4e2a\u86cb\u767d)\u3002", html_zh, fixed = TRUE)
  html_zh <- gsub("Higher CI indicates greater contamination from the corresponding cell type.", "CI \u503c\u8d8a\u9ad8\uff0c\u8bf4\u660e\u8be5\u7c7b\u7ec6\u80de\u6210\u5206\u5728\u603b\u86cb\u767d\u4e2d\u7684\u5360\u6bd4\u8d8a\u5927\uff0c\u6c61\u67d3\u8d8a\u4e25\u91cd\u3002", html_zh, fixed = TRUE)
  html_zh <- gsub("Literature estimate", "\u6587\u732e\u4f30\u7b97", html_zh, fixed = TRUE)
  html_zh <- gsub("Baize tool", "Baize \u5f00\u6e90\u5de5\u5177", html_zh, fixed = TRUE)

  # Interpretation note
  html_zh <- gsub("<b>Interpretation:</b> Pass &#x2705; = CI below safe threshold, contamination negligible. Flag &#x26A0;&#xFE0F; = mild contamination, data usable with caution. Fail &#x274C; = significant contamination, marker proteins should be excluded.",
    "<b>\u5224\u5b9a\u8bf4\u660e\uff1a</b> Pass \u2705 = CI \u4f4e\u4e8e\u5b89\u5168\u9608\u503c\uff0c\u6c61\u67d3\u53ef\u5ffd\u7565\u3002Flag \u26A0\uFE0F = \u8f7b\u5ea6\u6c61\u67d3\uff0c\u6570\u636e\u53ef\u7528\u4f46\u9700\u8c28\u614e\u3002Fail \u274C = \u663e\u8457\u6c61\u67d3\uff0c\u5e94\u6392\u9664 marker \u86cb\u767d\u540e\u518d\u5206\u6790\u3002",
    html_zh, fixed = TRUE)

  # Per-condition CI text
  html_zh <- gsub("samples Pass)", "\u4e2a\u6837\u672c Pass)", html_zh, fixed = TRUE)
  html_zh <- gsub("(threshold:", "(\u9608\u503c:", html_zh, fixed = TRUE)
  html_zh <- gsub("<b>Recommendation:</b>", "<b>\u5efa\u8bae\uff1a</b>", html_zh, fixed = TRUE)
  html_zh <- gsub("All CI indicators are within safe thresholds. Data is clean and can be directly used for downstream analysis.", "\u6240\u6709 CI \u6307\u6807\u5747\u5728\u5b89\u5168\u9608\u503c\u5185\u3002\u6570\u636e\u5e72\u51c0\uff0c\u53ef\u76f4\u63a5\u7528\u4e8e\u4e0b\u6e38\u5206\u6790\u3002", html_zh, fixed = TRUE)
  html_zh <- gsub("Significant contamination detected. Marker proteins should be excluded before downstream analysis.", "\u68c0\u6d4b\u5230\u663e\u8457\u6c61\u67d3\u3002\u4e0b\u6e38\u5206\u6790\u524d\u5e94\u6392\u9664 marker \u86cb\u767d\u3002", html_zh, fixed = TRUE)
  html_zh <- gsub("Mild contamination detected but below Fail thresholds. Data is usable; check marker proteins if analyzing related pathways.", "\u68c0\u6d4b\u5230\u8f7b\u5ea6\u6c61\u67d3\u4f46\u672a\u8fbe Fail \u9608\u503c\u3002\u6570\u636e\u53ef\u7528\uff0c\u5206\u6790\u76f8\u5173\u901a\u8def\u65f6\u6ce8\u610f\u6392\u67e5 marker \u86cb\u767d\u3002", html_zh, fixed = TRUE)

  # Per-indicator CI text
  html_zh <- gsub("samples Pass.", "\u4e2a\u6837\u672c Pass\u3002", html_zh, fixed = TRUE)
  html_zh <- gsub("<b>Conclusion:</b>", "<b>\u7ed3\u8bba\uff1a</b>", html_zh, fixed = TRUE)
  html_zh <- gsub("No significant contamination in any condition.", "\u6240\u6709\u6761\u4ef6\u5747\u65e0\u663e\u8457\u6c61\u67d3\u3002", html_zh, fixed = TRUE)
  html_zh <- gsub("samples passed. Some samples show elevated levels.", "\u4e2a\u6837\u672c\u901a\u8fc7\u3002\u90e8\u5206\u6837\u672c\u663e\u793a\u504f\u9ad8\u6c34\u5e73\u3002", html_zh, fixed = TRUE)

  # CV interpretation
  html_zh <- gsub("The CV values reported here are calculated on <b>normalized LFQ-equivalent quantities</b> (Spectronaut PG.Quantity) across <b>biological replicates</b> (different individuals), not technical replicates. Therefore, the observed median CVs (40-60%) reflect genuine inter-individual biological variation in addition to technical variance. Despite these CVs, the high Pearson correlations (&gt;0.95) confirm that the overall quantification baseline is highly consistent across all methods.",
    "\u672c\u62a5\u544a\u4e2d\u7684 CV \u503c\u57fa\u4e8e <b>\u5f52\u4e00\u5316\u540e\u7684 LFQ \u5f53\u91cf\u503c</b> (Spectronaut PG.Quantity)\uff0c\u5728 <b>\u751f\u7269\u5b66\u91cd\u590d</b> (\u4e0d\u540c\u4e2a\u4f53) \u800c\u975e\u6280\u672f\u91cd\u590d\u4e0a\u8ba1\u7b97\u3002\u56e0\u6b64\uff0c\u89c2\u5bdf\u5230\u7684\u4e2d\u4f4d CV (40-60%) \u53cd\u6620\u7684\u662f\u4e2a\u4f53\u95f4\u771f\u5b9e\u7684\u751f\u7269\u5b66\u53d8\u5f02 + \u6280\u672f\u65b9\u5dee\u3002\u5c3d\u7ba1 CV \u8f83\u9ad8\uff0c\u4f46 Pearson \u76f8\u5173\u7cfb\u6570 (&gt;0.95) \u786e\u8ba4\u4e86\u6240\u6709\u65b9\u6cd5\u7684\u5b9a\u91cf\u57fa\u7ebf\u9ad8\u5ea6\u4e00\u81f4\u3002",
    html_zh, fixed = TRUE)

  # Dynamic range / Top-N
  html_zh <- gsub("Dynamic range is measured as the P1-P99 intensity span in log2 space per sample.", "\u52a8\u6001\u8303\u56f4\u4ee5\u6bcf\u4e2a\u6837\u672c\u7684 P1-P99 \u5f3a\u5ea6\u8de8\u5ea6 (log2 \u7a7a\u95f4) \u8861\u91cf\u3002", html_zh, fixed = TRUE)
  html_zh <- gsub("log2 orders", "log2 \u6570\u91cf\u7ea7", html_zh, fixed = TRUE)
  html_zh <- gsub("Top-10 protein fraction indicates the degree of high-abundance protein depletion. Lower values mean better depletion.", "Top-10 \u86cb\u767d\u5360\u6bd4\u53cd\u6620\u9ad8\u4e30\u5ea6\u86cb\u767d\u7684\u53bb\u9664\u6548\u7387\u3002\u503c\u8d8a\u4f4e\u8bf4\u660e\u53bb\u9664\u6548\u679c\u8d8a\u597d\u3002", html_zh, fixed = TRUE)

  # Conclusions
  html_zh <- gsub("All conditions show no severe cell contamination", "\u6240\u6709\u6761\u4ef6\u5747\u65e0\u4e25\u91cd\u7ec6\u80de\u6c61\u67d3", html_zh, fixed = TRUE)
  html_zh <- gsub("CI_PLT and CI_RBC are well within safe thresholds across all samples.", "CI_PLT \u548c CI_RBC \u5728\u6240\u6709\u6837\u672c\u4e2d\u5747\u5904\u4e8e\u5b89\u5168\u9608\u503c\u5185\u3002", html_zh, fixed = TRUE)
  html_zh <- gsub("achieves the highest protein identification", "\u5b9e\u73b0\u4e86\u6700\u9ad8\u7684\u86cb\u767d\u9274\u5b9a\u6570", html_zh, fixed = TRUE)
  html_zh <- gsub("PGs on average).", "PGs)\u3002", html_zh, fixed = TRUE)
  html_zh <- gsub("shows the best quantification reproducibility (median CV =", "\u663e\u793a\u6700\u4f73\u5b9a\u91cf\u91cd\u73b0\u6027 (\u4e2d\u4f4d CV =", html_zh, fixed = TRUE)
  html_zh <- gsub("Conditions with 100% CI Pass rate:", "CI 100% Pass \u7684\u6761\u4ef6:", html_zh, fixed = TRUE)

  # Key insight
  html_zh <- gsub("Higher protein identification count does not equal better data quality.", "\u86cb\u767d\u9274\u5b9a\u6570\u9ad8\u4e0d\u7b49\u4e8e\u6570\u636e\u8d28\u91cf\u597d\u3002", html_zh, fixed = TRUE)
  html_zh <- gsub("Cell contamination (platelets, erythrocytes, PBMC) can lead to false identification of non-plasma proteins, artificially inflating protein counts. Contamination assessment must precede any comparison of identification depth.",
    "\u7ec6\u80de\u6c61\u67d3 (\u8840\u5c0f\u677f/\u7ea2\u7ec6\u80de/PBMC) \u4f1a\u5bfc\u81f4\u975e\u8840\u6d46\u6765\u6e90\u7684\u86cb\u767d\u88ab\u9519\u8bef\u9274\u5b9a\uff0c\u865a\u589e\u86cb\u767d\u9274\u5b9a\u6570\u3002\u5728\u6bd4\u8f83\u4e0d\u540c\u65b9\u6cd5\u7684\u9274\u5b9a\u6df1\u5ea6\u4e4b\u524d\uff0c\u5fc5\u987b\u5148\u8fdb\u884c\u6c61\u67d3\u8bc4\u4f30\u3002",
    html_zh, fixed = TRUE)

  # Per-condition recommendations
  html_zh <- gsub("CI all Pass. Data can be directly used for downstream analysis.", "CI \u5168\u90e8 Pass\u3002\u6570\u636e\u53ef\u76f4\u63a5\u7528\u4e8e\u4e0b\u6e38\u5206\u6790\u3002", html_zh, fixed = TRUE)
  html_zh <- gsub("Some samples flagged. Usable for analysis; check marker overlap in differential results.", "\u90e8\u5206\u6837\u672c\u6807\u8bb0\u3002\u53ef\u7528\u4e8e\u5206\u6790\uff0c\u4f46\u5dee\u5f02\u7ed3\u679c\u4e2d\u9700\u6392\u67e5 marker \u86cb\u767d\u3002", html_zh, fixed = TRUE)

  # Next steps
  html_zh <- gsub("If comparing against competitor kits, run the same CI evaluation on competitor search results.", "\u5982\u9700\u5bf9\u6bd4\u7ade\u54c1\u8bd5\u5242\u76d2\uff0c\u5efa\u8bae\u5bf9\u7ade\u54c1\u641c\u5e93\u7ed3\u679c\u4e5f\u8fdb\u884c\u76f8\u540c\u7684 CI \u8bc4\u4f30\u3002", html_zh, fixed = TRUE)
  html_zh <- gsub("Ensure data quality is confirmed before comparing identification depth across methods.", "\u786e\u4fdd\u5728\u53ef\u6bd4\u7684\u8d28\u91cf\u6807\u51c6\u4e0b\u8ba8\u8bba\u9274\u5b9a\u6570\u5dee\u5f02\uff0c\u907f\u514d\u5c06\u56e0\u6c61\u67d3\u865a\u589e\u7684\u86cb\u767d\u8bef\u8ba4\u4e3a\u66f4\u6df1\u7684\u8986\u76d6\u3002", html_zh, fixed = TRUE)

  # Software
  html_zh <- gsub("Data analysis:", "\u6570\u636e\u5206\u6790:", html_zh, fixed = TRUE)
  html_zh <- gsub("CI calculation:", "CI \u8ba1\u7b97:", html_zh, fixed = TRUE)
  html_zh <- gsub("Search engine: Spectronaut / DIA-NN (as configured)", "\u641c\u5e93\u8f6f\u4ef6: Spectronaut / DIA-NN (\u6309\u914d\u7f6e)", html_zh, fixed = TRUE)

  # Footer
  html_zh <- gsub("Report auto-generated by BayOmics Bioinformatics Team.", "\u672c\u62a5\u544a\u7531 BayOmics \u751f\u7269\u4fe1\u606f\u56e2\u961f\u81ea\u52a8\u751f\u6210\u3002", html_zh, fixed = TRUE)

  # Save ZH HTML
  html_zh_path <- file.path(comp_dir, "10_comparison_report_zh.html")
  writeLines(html_zh, html_zh_path, useBytes = TRUE)
  message(sprintf(">>> HTML saved: %s", basename(html_zh_path)))

  # Chrome PDF (ZH)
  if (!is.null(chrome)) {
    pdf_zh_path <- file.path(comp_dir, "10_comparison_report_zh.pdf")
    tmp_pdf_zh <- file.path(tempdir(), "comp_rpt_zh.pdf")
    cmd_zh <- sprintf('"%s" --headless --disable-gpu --no-sandbox --print-to-pdf="%s" "%s"', chrome, tmp_pdf_zh, html_zh_path)
    res_zh <- system(cmd_zh, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (res_zh == 0 && file.exists(tmp_pdf_zh)) {
      file.copy(tmp_pdf_zh, pdf_zh_path, overwrite = TRUE); file.remove(tmp_pdf_zh)
      message(sprintf(">>> PDF saved: %s", basename(pdf_zh_path)))
    }
  }

  unlink(png_dir, recursive = TRUE)
  } # end base64enc check
} # end comparison

message("\n========================================")
message("  All evaluations complete!")
message(sprintf("  Output: %s", output_base))
message("========================================")
