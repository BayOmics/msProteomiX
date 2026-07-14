# ==============================================================================
# msProteomiX -- Contamination Ratio (CR) Analysis
# timsTOF Pro IM-MS heatmap contamination evaluation
# ==============================================================================
# Reproduces the analysis from:
#   "Multimodal single cell-resolved spatial proteomics reveal
#    pancreatic tumor heterogeneity" (Fig. 2f, 2g)
#
# Upstream:
#   AlphaTims Python script (extract_precursors.py) -> *_cr_summary.csv
#
# Functions:
#   [1] import_cr_results()   -- Import CR summary CSVs
#   [2] plot_cr_gradient()    -- Fig. 2f upper: CR vs RT curve
#   [3] plot_cr_tic()         -- Fig. 2f lower: TIC intensity curve
#   [4] plot_cr_heatmap()     -- Fig. 2g: IM-m/z heatmap with dividing line
#   [5] generate_cr_report()  -- Combined HTML report
# ==============================================================================


# ==============================================================================
# [1] import_cr_results()
# ==============================================================================

#' Import Contamination Ratio results from AlphaTims Python script
#'
#' Reads one or more \code{*_cr_summary.csv} files produced by
#' \code{extract_precursors.py} and returns a combined data.frame
#' with group assignments.
#'
#' @param cr_dir Directory containing \code{*_cr_summary.csv} files.
#' @param pattern File name pattern to match (default: \code{"_cr_summary\\\\.csv$"}).
#' @param group_map Named list for group assignment. Each element name is a
#'   group label, and each value is a character vector of patterns.
#'   Samples matching any pattern are assigned to that group.
#'   Example: \code{list("HeLa" = "Hela", "W/o Wash" = "DDM")}.
#'   If \code{NULL}, all samples are assigned group "default".
#' @return data.frame with columns: sample, file_name, rt_min, tot_inten,
#'   contamination_ratio, group
#' @export
import_cr_results <- function(cr_dir = "wkdir/cr_output",
                               pattern = "_cr_summary\\.csv$",
                               group_map = NULL) {
  if (!dir.exists(cr_dir)) {
    stop("CR output directory not found: ", cr_dir,
         "\n  Run extract_precursors.py first.")
  }

  csv_files <- list.files(cr_dir, pattern = pattern, full.names = TRUE)
  if (length(csv_files) == 0) {
    stop("No *_cr_summary.csv files found in: ", cr_dir,
         "\n  Run extract_precursors.py first.")
  }

  message(sprintf(">>> Loading CR results: %d files from %s",
                  length(csv_files), cr_dir))

  # Read and combine all summary CSVs
  all_data <- lapply(csv_files, function(f) {
    df <- read.csv(f, stringsAsFactors = FALSE)
    # Extract sample name from filename
    df$sample <- sub("_cr_summary\\.csv$", "", basename(f))
    df
  })
  combined <- do.call(rbind, all_data)
  rownames(combined) <- NULL

  # Assign groups
  if (!is.null(group_map)) {
    combined$group <- "Other"
    for (grp_name in names(group_map)) {
      pats <- group_map[[grp_name]]
      for (pat in pats) {
        matched <- grepl(pat, combined$sample, ignore.case = TRUE)
        combined$group[matched] <- grp_name
      }
    }
  } else {
    combined$group <- "default"
  }

  # Ensure correct column types
  combined$rt_min <- as.numeric(combined$rt_min)
  combined$tot_inten <- as.numeric(combined$tot_inten)
  combined$contamination_ratio <- as.numeric(combined$contamination_ratio)

  n_samples <- length(unique(combined$sample))
  n_groups <- length(unique(combined$group))
  message(sprintf("    %d samples, %d groups, %d data points",
                  n_samples, n_groups, nrow(combined)))

  combined
}


# ==============================================================================
# [2] plot_cr_gradient()
# ==============================================================================

#' Plot Contamination Ratio along LC gradient
#'
#' Reproduces the upper panel of Fig. 2f: CR vs RT with mean line
#' and \eqn{\pm} SD ribbon per group.
#'
#' @param cr_data data.frame from \code{\link{import_cr_results}}.
#' @param group_col Column name for grouping (default: \code{"group"}).
#' @param sd_factor SD multiplier for ribbon (default: 0.5, matching paper).
#' @param xlim RT range in minutes (default: \code{NULL}, auto).
#' @param ylim CR range (default: \code{NULL}, auto).
#' @param output_dir Output directory for saving plots (default: \code{"output"}).
#' @param project_name Project name prefix (default: \code{"Project"}).
#' @return ggplot object (invisibly)
#' @export
plot_cr_gradient <- function(cr_data, group_col = "group",
                              sd_factor = 0.5,
                              xlim = NULL, ylim = NULL,
                              output_dir = "output",
                              project_name = "Project") {
  ensure_output_dir(output_dir)

  # Aggregate: mean + sd per group per RT bin
  # Round RT to nearest integer for grouping across replicates
  cr_data$rt_bin <- round(cr_data$rt_min)

  groups <- unique(cr_data[[group_col]])
  n_groups <- length(groups)

  agg_list <- lapply(groups, function(g) {
    sub_df <- cr_data[cr_data[[group_col]] == g, ]
    agg <- stats::aggregate(
      contamination_ratio ~ rt_bin,
      data = sub_df,
      FUN = function(x) c(mean = mean(x, na.rm = TRUE),
                          sd = stats::sd(x, na.rm = TRUE))
    )
    result <- data.frame(
      rt_bin = agg$rt_bin,
      mean_cr = agg$contamination_ratio[, "mean"],
      sd_cr = agg$contamination_ratio[, "sd"],
      group = g,
      stringsAsFactors = FALSE
    )
    # Replace NA sd with 0 (single replicate)
    result$sd_cr[is.na(result$sd_cr)] <- 0
    result
  })
  agg_df <- do.call(rbind, agg_list)
  agg_df$ymin <- agg_df$mean_cr - sd_factor * agg_df$sd_cr
  agg_df$ymax <- agg_df$mean_cr + sd_factor * agg_df$sd_cr
  agg_df$ymin[agg_df$ymin < 0] <- 0

  p <- ggplot2::ggplot(agg_df,
                       ggplot2::aes(x = .data$rt_bin, y = .data$mean_cr,
                                    color = .data$group, fill = .data$group)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$ymin, ymax = .data$ymax),
                         alpha = 0.15, color = NA) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 1.5) +
    ggplot2::scale_color_manual(values = mspx_colors(n_groups)) +
    ggplot2::scale_fill_manual(values = mspx_colors(n_groups)) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Contamination Ratio along LC Gradient",
      subtitle = sprintf("Ribbon = mean %s %.1f SD", "\u00b1", sd_factor),
      x = "Retention Time (min)",
      y = "Contamination Ratio (%)",
      color = "Group",
      fill = "Group"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, size = 9,
                                             color = "grey50"),
      legend.position = "bottom"
    )

  if (!is.null(xlim)) p <- p + ggplot2::xlim(xlim)
  if (!is.null(ylim)) p <- p + ggplot2::ylim(ylim)

  save_plot_and_data(p, agg_df, project_name, "CR_Gradient",
                     output_dir = output_dir, width = 8, height = 5)
  invisible(p)
}


# ==============================================================================
# [3] plot_cr_tic()
# ==============================================================================

#' Plot TIC intensity along LC gradient
#'
#' Reproduces part of the lower panel of Fig. 2f: total ion intensity
#' (log2-transformed) at each RT slice.
#'
#' @param cr_data data.frame from \code{\link{import_cr_results}}.
#' @param group_col Column name for grouping (default: \code{"group"}).
#' @param output_dir Output directory (default: \code{"output"}).
#' @param project_name Project name (default: \code{"Project"}).
#' @return ggplot object (invisibly)
#' @export
plot_cr_tic <- function(cr_data, group_col = "group",
                         output_dir = "output",
                         project_name = "Project") {
  ensure_output_dir(output_dir)

  groups <- unique(cr_data[[group_col]])
  n_groups <- length(groups)

  p <- ggplot2::ggplot(cr_data,
                       ggplot2::aes(x = .data$rt_min, y = .data$tot_inten,
                                    color = .data[[group_col]],
                                    group = .data$sample)) +
    ggplot2::geom_line(alpha = 0.6, linewidth = 0.5) +
    ggplot2::geom_point(size = 1) +
    ggplot2::scale_color_manual(values = mspx_colors(n_groups)) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "Total Ion Intensity along LC Gradient",
      x = "Retention Time (min)",
      y = "log2(Total Intensity)",
      color = "Group"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )

  save_plot_and_data(p, cr_data[, c("sample", "rt_min", "tot_inten",
                                     group_col)],
                     project_name, "CR_TIC",
                     output_dir = output_dir, width = 8, height = 4)
  invisible(p)
}


# ==============================================================================
# [4] plot_cr_heatmap()
# ==============================================================================

#' Plot IM-m/z heatmap with CR dividing line
#'
#' Reproduces Fig. 2g: scatter/density heatmap of precursors in the
#' trapped ion mobility (IM) vs m/z space at a specific RT, with the
#' CR dividing line overlaid as a white line.
#'
#' @param slice_csv Path to a single RT-slice CSV file exported by
#'   \code{extract_precursors.py} (with \code{--export_slices}).
#'   Must contain columns: \code{mz_values}, \code{mobility_values},
#'   \code{corrected_intensity_values}.
#' @param line_point1 Numeric vector \code{c(mz, im)} for the first point
#'   of the dividing line (default: \code{c(350, 0.8)}).
#' @param line_point2 Numeric vector \code{c(mz, im)} for the second point
#'   of the dividing line (default: \code{c(950, 1.3)}).
#' @param log_intensity If \code{TRUE}, use log10 transform for color scale
#'   (default: \code{TRUE}).
#' @param title Plot title (default: auto-generated from file name).
#' @param output_dir Output directory (default: \code{"output"}).
#' @param project_name Project name (default: \code{"Project"}).
#' @return ggplot object (invisibly)
#' @export
plot_cr_heatmap <- function(slice_csv,
                             line_point1 = c(350, 0.8),
                             line_point2 = c(950, 1.3),
                             log_intensity = TRUE,
                             title = NULL,
                             output_dir = "output",
                             project_name = "Project") {
  ensure_output_dir(output_dir)

  if (!file.exists(slice_csv)) {
    stop("Slice CSV not found: ", slice_csv,
         "\n  Run extract_precursors.py with --export_slices first.")
  }

  df <- read.csv(slice_csv, stringsAsFactors = FALSE)

  required_cols <- c("mz_values", "mobility_values", "corrected_intensity_values")
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) {
    stop("Missing columns in CSV: ", paste(missing, collapse = ", "))
  }

  # Remove zero-intensity rows
  df <- df[df$corrected_intensity_values > 0, ]
  if (nrow(df) == 0) {
    message(">>> No non-zero intensity precursors in: ", slice_csv)
    return(invisible(NULL))
  }

  # Compute dividing line
  k <- (line_point2[2] - line_point1[2]) / (line_point2[1] - line_point1[1])
  b <- line_point1[2] - k * line_point1[1]

  # Intensity for color mapping
  if (log_intensity) {
    df$plot_inten <- log10(df$corrected_intensity_values + 1)
    inten_label <- "log10(Intensity)"
  } else {
    df$plot_inten <- df$corrected_intensity_values
    inten_label <- "Intensity"
  }

  # Auto title
  if (is.null(title)) {
    title <- sub("\\.csv$", "", basename(slice_csv))
  }

  # Extract RT from the data for subtitle
  rt_val <- NULL
  if ("rt_values_min" %in% colnames(df)) {
    rt_val <- round(unique(df$rt_values_min)[1], 1)
  }

  # Dividing line segment (clip to data range)
  mz_range <- range(df$mz_values)
  line_df <- data.frame(
    x = mz_range,
    y = k * mz_range + b
  )

  # Build plot (black background, viridis-like color scale)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$mz_values,
                                         y = .data$mobility_values,
                                         color = .data$plot_inten)) +
    ggplot2::geom_point(size = 0.3, alpha = 0.6) +
    ggplot2::geom_line(data = line_df,
                       ggplot2::aes(x = .data$x, y = .data$y),
                       color = "white", linewidth = 0.6,
                       linetype = "dashed", inherit.aes = FALSE) +
    ggplot2::scale_color_viridis_c(option = "inferno", name = inten_label) +
    ggplot2::theme_dark() +
    ggplot2::theme(
      panel.background = ggplot2::element_rect(fill = "black"),
      plot.background = ggplot2::element_rect(fill = "grey10"),
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold",
                                          color = "white"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey70",
                                             size = 9),
      axis.title = ggplot2::element_text(color = "white"),
      axis.text = ggplot2::element_text(color = "grey80"),
      legend.background = ggplot2::element_rect(fill = "grey20"),
      legend.text = ggplot2::element_text(color = "white"),
      legend.title = ggplot2::element_text(color = "white")
    ) +
    ggplot2::labs(
      title = title,
      subtitle = if (!is.null(rt_val)) sprintf("RT = %.1f min", rt_val) else NULL,
      x = "m/z",
      y = expression("1/K"[0] ~ "(V" %.% "s" %.% "cm"^{-2} * ")")
    )

  # Determine suffix for file naming
  suffix <- paste0("CR_Heatmap_", sub("\\.csv$", "", basename(slice_csv)))

  # Save with black background
  base_name <- file.path(output_dir,
                         paste0(project_name, "_", suffix))
  ggplot2::ggsave(paste0(base_name, ".pdf"), p,
                  width = 8, height = 6, bg = "grey10")
  ggplot2::ggsave(paste0(base_name, ".png"), p,
                  width = 8, height = 6, dpi = 300, bg = "grey10")
  utils::write.csv(
    df[, c("mz_values", "mobility_values", "corrected_intensity_values")],
    paste0(base_name, ".csv"), row.names = FALSE
  )
  message(sprintf("\u2705 Generated: %s", suffix))
  print(p)

  invisible(p)
}


# ==============================================================================
# [5] generate_cr_report()
# ==============================================================================

#' Generate Contamination Ratio HTML report
#'
#' Combines CR gradient, TIC, and heatmap plots into a single HTML report.
#'
#' @param cr_data data.frame from \code{\link{import_cr_results}}.
#' @param heatmap_csvs Character vector of per-slice CSV paths for heatmap
#'   (from \code{extract_precursors.py --export_slices}).
#'   If \code{NULL}, heatmap section is skipped.
#' @param group_col Column name for grouping (default: \code{"group"}).
#' @param sd_factor SD multiplier for CR ribbon (default: 0.5).
#' @param line_point1 Dividing line point 1 (default: \code{c(350, 0.8)}).
#' @param line_point2 Dividing line point 2 (default: \code{c(950, 1.3)}).
#' @param output_dir Output directory (default: \code{"output"}).
#' @param project_name Project name (default: \code{"Project"}).
#' @return Path to generated HTML report (invisibly)
#' @export
generate_cr_report <- function(cr_data,
                                heatmap_csvs = NULL,
                                group_col = "group",
                                sd_factor = 0.5,
                                line_point1 = c(350, 0.8),
                                line_point2 = c(950, 1.3),
                                output_dir = "output",
                                project_name = "Project") {
  ensure_output_dir(output_dir)
  message("\n=== Generating CR Report ===\n")

  # Collect PNG paths for embedding
  png_files <- character(0)

  # 1. CR Gradient
  message(">>> [1/3] CR Gradient plot...")
  plot_cr_gradient(cr_data, group_col = group_col, sd_factor = sd_factor,
                   output_dir = output_dir, project_name = project_name)
  png_files <- c(png_files,
                 file.path(output_dir, paste0(project_name, "_CR_Gradient.png")))

  # 2. TIC
  message(">>> [2/3] TIC intensity plot...")
  plot_cr_tic(cr_data, group_col = group_col,
              output_dir = output_dir, project_name = project_name)
  png_files <- c(png_files,
                 file.path(output_dir, paste0(project_name, "_CR_TIC.png")))

  # 3. Heatmaps
  if (!is.null(heatmap_csvs) && length(heatmap_csvs) > 0) {
    message(sprintf(">>> [3/3] IM-m/z heatmaps (%d)...", length(heatmap_csvs)))
    for (csv_path in heatmap_csvs) {
      if (!file.exists(csv_path)) {
        message("    Skipping (not found): ", csv_path)
        next
      }
      plot_cr_heatmap(csv_path,
                      line_point1 = line_point1,
                      line_point2 = line_point2,
                      output_dir = output_dir,
                      project_name = project_name)
      suffix <- paste0("CR_Heatmap_", sub("\\.csv$", "", basename(csv_path)))
      png_files <- c(png_files,
                     file.path(output_dir,
                               paste0(project_name, "_", suffix, ".png")))
    }
  } else {
    message(">>> [3/3] No heatmap CSVs provided. Skipping.")
  }

  # Build HTML
  html_path <- file.path(output_dir,
                         paste0(project_name, "_CR_Report.html"))
  .build_cr_html(html_path, png_files, cr_data, project_name)

  message(sprintf("\n>>> CR Report saved: %s", html_path))
  invisible(html_path)
}


#' @keywords internal
.build_cr_html <- function(html_path, png_files, cr_data, project_name) {
  # Encode PNGs as base64
  img_tags <- vapply(png_files, function(f) {
    if (!file.exists(f)) return("")
    raw <- readBin(f, "raw", file.info(f)$size)
    b64 <- base64enc::base64encode(raw)
    paste0(
      '<div style="text-align:center;margin:20px 0;">',
      '<img src="data:image/png;base64,', b64,
      '" style="max-width:100%;border:1px solid #ddd;border-radius:8px;">',
      '<p style="color:#666;font-size:12px;">', basename(f), '</p></div>'
    )
  }, character(1))

  # Summary table
  n_samples <- length(unique(cr_data$sample))
  n_groups <- length(unique(cr_data$group))
  mean_cr <- mean(cr_data$contamination_ratio, na.rm = TRUE)

  html_content <- paste0(
    '<!DOCTYPE html><html><head><meta charset="utf-8">',
    '<title>', project_name, ' - CR Report</title>',
    '<style>body{font-family:Arial,sans-serif;max-width:900px;margin:0 auto;',
    'padding:20px;background:#fafafa;}h1{color:#333;border-bottom:2px solid #4CAF50;',
    'padding-bottom:10px;}h2{color:#555;margin-top:30px;}',
    'table{border-collapse:collapse;margin:15px 0;}',
    'td,th{border:1px solid #ddd;padding:8px 12px;text-align:left;}',
    'th{background:#4CAF50;color:white;}</style></head><body>',
    '<h1>', project_name, ' - Contamination Ratio Report</h1>',
    '<h2>Summary</h2>',
    '<table><tr><th>Metric</th><th>Value</th></tr>',
    '<tr><td>Samples</td><td>', n_samples, '</td></tr>',
    '<tr><td>Groups</td><td>', n_groups, '</td></tr>',
    '<tr><td>Mean CR</td><td>', sprintf("%.1f%%", mean_cr), '</td></tr>',
    '</table>',
    '<h2>CR along LC Gradient</h2>',
    if (length(img_tags) >= 1) img_tags[1] else "",
    '<h2>TIC Intensity</h2>',
    if (length(img_tags) >= 2) img_tags[2] else "",
    if (length(img_tags) >= 3) paste0('<h2>IM-m/z Heatmaps</h2>',
                                       paste(img_tags[3:length(img_tags)],
                                             collapse = "\n")) else "",
    '</body></html>'
  )

  writeLines(html_content, html_path)
}
