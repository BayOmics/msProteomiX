# ==============================================================================
# msProteomiX — Contamination Index Visualization
# ==============================================================================


#' Plot Contamination Index barplot
#'
#' Creates a 3-panel faceted barplot showing CI values for PLT, RBC, and PBMC
#' panels with Pass/Flag/Fail threshold lines.
#'
#' @param ci_result data.frame from \code{calc_contamination_index()}
#' @param output_dir Output directory for saving the plot
#' @param prefix Filename prefix (default "CI")
#' @return Invisible ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' ci <- calc_contamination_index(ms)
#' plot_ci_barplot(ci, output_dir = "output")
#' }
plot_ci_barplot <- function(ci_result, output_dir = "output", prefix = "CI") {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required for CI visualization.")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  jcol <- c(Pass = "#27ae60", Flag = "#f39c12", Fail = "#e74c3c",
            "N/A" = "#bdc3c7")

  # Build long format for faceted plot
  ci_long <- rbind(
    data.frame(Sample = ci_result$Sample, Panel = "CI_Platelet (Gao2026)",
               Value = ci_result$CI_PLT, Judgment = ci_result$J_PLT,
               stringsAsFactors = FALSE),
    data.frame(Sample = ci_result$Sample, Panel = "CI_Erythrocyte (Gao2026)",
               Value = ci_result$CI_RBC, Judgment = ci_result$J_RBC,
               stringsAsFactors = FALSE),
    data.frame(Sample = ci_result$Sample, Panel = "CI_PBMC (Korff2025)",
               Value = ci_result$CI_PBMC, Judgment = ci_result$J_PBMC,
               stringsAsFactors = FALSE)
  )
  ci_long$Panel <- factor(ci_long$Panel,
    levels = c("CI_Platelet (Gao2026)", "CI_Erythrocyte (Gao2026)", "CI_PBMC (Korff2025)"))

  # Clean sample names (remove .raw suffix)
  ci_long$Sample <- gsub("\\.raw$", "", ci_long$Sample)

  # Threshold lines
  thresh_df <- data.frame(
    Panel = factor(c(
      rep("CI_Platelet (Gao2026)", 2),
      "CI_Erythrocyte (Gao2026)",
      rep("CI_PBMC (Korff2025)", 2)
    ), levels = levels(ci_long$Panel)),
    yval = c(0.009, 0.05, 0.1, 0.005, 0.02),
    Label = c("Pass 0.009", "Fail 0.05", "Fail 0.1", "Pass 0.005", "Fail 0.02")
  )

  p <- ggplot2::ggplot(ci_long, ggplot2::aes(x = .data$Sample, y = .data$Value, fill = .data$Judgment)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::facet_wrap(~ .data$Panel, scales = "free_y", ncol = 1) +
    ggplot2::geom_hline(data = thresh_df, ggplot2::aes(yintercept = .data$yval),
                        linetype = "dashed", color = "grey40", linewidth = 0.5) +
    ggplot2::geom_text(data = thresh_df,
                       ggplot2::aes(x = Inf, y = .data$yval, label = .data$Label),
                       hjust = 1.1, vjust = -0.5, size = 2.8, color = "grey40",
                       inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(values = jcol, name = "Judgment") +
    ggplot2::labs(title = "Contamination Index (CI) Evaluation",
                  x = NULL, y = "Contamination Index") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
      strip.text = ggplot2::element_text(face = "bold")
    )

  # Save PDF + PNG
  pdf_path <- file.path(output_dir, paste0(prefix, "_barplot.pdf"))
  png_path <- file.path(output_dir, paste0(prefix, "_barplot.png"))
  ggplot2::ggsave(pdf_path, p, width = 8, height = 8, dpi = 300)
  ggplot2::ggsave(png_path, p, width = 8, height = 8, dpi = 300)
  message(sprintf(">>> CI barplot saved: %s", pdf_path))

  invisible(p)
}


#' Plot Contamination Index overview table as a heatmap-style plot
#'
#' Displays a styled summary of CI values and judgments per sample as a
#' tile-based heatmap with color-coded Pass/Flag/Fail status.
#'
#' @param ci_result data.frame from \code{calc_contamination_index()}
#' @param output_dir Output directory for saving the plot
#' @param prefix Filename prefix (default "CI")
#' @return Invisible ggplot object
#' @export
plot_ci_summary <- function(ci_result, output_dir = "output", prefix = "CI") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  jcol <- c(Pass = "#27ae60", Flag = "#f39c12", Fail = "#e74c3c",
            "N/A" = "#bdc3c7")

  # Build long-format judgment table
  samples_clean <- gsub("\\.raw$", "", ci_result$Sample)

  judge_long <- rbind(
    data.frame(Sample = samples_clean, Panel = "PLT",
               Judgment = ci_result$J_PLT,
               Label = ifelse(is.na(ci_result$CI_PLT), "N/A",
                              sprintf("%.4f", ci_result$CI_PLT)),
               stringsAsFactors = FALSE),
    data.frame(Sample = samples_clean, Panel = "RBC",
               Judgment = ci_result$J_RBC,
               Label = ifelse(is.na(ci_result$CI_RBC), "N/A",
                              sprintf("%.4f", ci_result$CI_RBC)),
               stringsAsFactors = FALSE),
    data.frame(Sample = samples_clean, Panel = "PBMC",
               Judgment = ci_result$J_PBMC,
               Label = ifelse(is.na(ci_result$CI_PBMC), "N/A",
                              sprintf("%.4f", ci_result$CI_PBMC)),
               stringsAsFactors = FALSE),
    data.frame(Sample = samples_clean, Panel = "Overall",
               Judgment = ci_result$Overall,
               Label = ci_result$Overall,
               stringsAsFactors = FALSE)
  )
  judge_long$Panel <- factor(judge_long$Panel,
    levels = c("PLT", "RBC", "PBMC", "Overall"))
  judge_long$Sample <- factor(judge_long$Sample, levels = rev(samples_clean))

  p <- ggplot2::ggplot(judge_long,
    ggplot2::aes(x = .data$Panel, y = .data$Sample, fill = .data$Judgment)) +
    ggplot2::geom_tile(color = "white", linewidth = 1.5) +
    ggplot2::geom_text(ggplot2::aes(label = .data$Label), size = 3.2, color = "white",
                       fontface = "bold") +
    ggplot2::scale_fill_manual(values = jcol, name = "Status") +
    ggplot2::labs(title = "Contamination Index Summary",
                  x = "Marker Panel", y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      panel.grid = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(face = "bold")
    )

  pdf_path <- file.path(output_dir, paste0(prefix, "_summary.pdf"))
  png_path <- file.path(output_dir, paste0(prefix, "_summary.png"))
  ggplot2::ggsave(pdf_path, p, width = 6, height = max(4, nrow(ci_result) * 0.6 + 1), dpi = 300)
  ggplot2::ggsave(png_path, p, width = 6, height = max(4, nrow(ci_result) * 0.6 + 1), dpi = 300)
  message(sprintf(">>> CI summary saved: %s", pdf_path))

  invisible(p)
}


#' Plot marker protein intensity heatmap
#'
#' Creates a heatmap showing log2 intensities of matched marker proteins
#' across samples, with panel grouping on the y-axis.
#'
#' @param ms MsDataSet object
#' @param ci_result data.frame from \code{calc_contamination_index()}.
#'   Used to retrieve match info (optional, marker panels are re-loaded if NULL).
#' @param output_dir Output directory for saving the plot
#' @param prefix Filename prefix (default "CI")
#' @return Invisible NULL
#' @export
plot_ci_heatmap <- function(ms, ci_result = NULL, output_dir = "output", prefix = "CI") {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Load panels and match
  panel_df <- .load_marker_panels("default")
  prot_info <- ms$protein_info

  # Collect matched markers from default panels
  panels_to_show <- c("PLT_Gao2026", "RBC_Gao2026", "PBMC_Korff2025")
  all_rows <- list()
  all_labels <- list()

  for (pname in panels_to_show) {
    sub_df <- panel_df[panel_df$Panel == pname, ]
    if (nrow(sub_df) == 0) next
    if (sub_df$Match_By[1] == "uniprot") {
      rows <- .match_markers_uniprot(sub_df$UniProt_ID, prot_info)
    } else {
      rows <- .match_markers_gene(sub_df$Gene, prot_info)
    }
    if (length(rows) > 0) {
      genes <- as.character(prot_info[["Gene"]][rows])
      genes[is.na(genes) | !nzchar(genes)] <- as.character(prot_info[["Protein ID"]][rows[is.na(genes) | !nzchar(genes)]])
      # Take first gene before semicolon
      genes <- sub(";.*", "", genes)
      # Add panel prefix
      labels <- paste0("[", sub("_.*", "", pname), "] ", genes)
      all_rows[[pname]] <- rows
      all_labels[[pname]] <- labels
    }
  }

  if (length(unlist(all_rows)) == 0) {
    message(">>> No marker proteins matched. Skipping heatmap.")
    return(invisible(NULL))
  }

  # Build matrix
  mat <- as.matrix(ms$proteins)
  if (!isTRUE(ms$is_log2)) mat <- log2(mat + 1)
  combined_rows <- unlist(all_rows)
  combined_labels <- unlist(all_labels)
  mat_sub <- mat[combined_rows, , drop = FALSE]
  rownames(mat_sub) <- combined_labels
  colnames(mat_sub) <- gsub("\\.raw$", "", colnames(mat_sub))

  # Use pheatmap if available
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    # Panel annotation
    panel_annot <- data.frame(
      Panel = rep(names(all_rows), sapply(all_rows, length)),
      row.names = combined_labels
    )
    panel_annot$Panel <- sub("_.*", "", panel_annot$Panel)

    annot_colors <- list(Panel = c(PLT = "#e74c3c", RBC = "#3498db", PBMC = "#2ecc71"))

    pdf_path <- file.path(output_dir, paste0(prefix, "_marker_heatmap.pdf"))
    png_path <- file.path(output_dir, paste0(prefix, "_marker_heatmap.png"))

    grDevices::pdf(pdf_path, width = 8, height = max(6, length(combined_rows) * 0.18 + 2))
    pheatmap::pheatmap(mat_sub,
      cluster_rows = FALSE, cluster_cols = FALSE,
      annotation_row = panel_annot,
      annotation_colors = annot_colors,
      main = "Marker Protein Intensities (log2)",
      fontsize_row = 6, fontsize_col = 9,
      color = grDevices::colorRampPalette(c("#f7fbff", "#3498db", "#2c3e50"))(50),
      na_col = "grey90"
    )
    grDevices::dev.off()

    grDevices::png(png_path, width = 2400, height = max(1800, length(combined_rows) * 55 + 600), res = 300)
    pheatmap::pheatmap(mat_sub,
      cluster_rows = FALSE, cluster_cols = FALSE,
      annotation_row = panel_annot,
      annotation_colors = annot_colors,
      main = "Marker Protein Intensities (log2)",
      fontsize_row = 6, fontsize_col = 9,
      color = grDevices::colorRampPalette(c("#f7fbff", "#3498db", "#2c3e50"))(50),
      na_col = "grey90"
    )
    grDevices::dev.off()

    message(sprintf(">>> CI marker heatmap saved: %s", pdf_path))
  } else {
    message(">>> Package 'pheatmap' not available. Skipping marker heatmap.")
    message("    Install with: install.packages('pheatmap')")
  }

  invisible(NULL)
}
