# ==============================================================================
# msProteomiX — Differential Protein Heatmap Visualization
# ==============================================================================

#' Plot differential protein clustering heatmap
#'
#' Generate a pheatmap-style clustering heatmap for top differentially expressed
#' proteins. Rows are Z-score normalized.
#'
#' @param diff_result Output from run_diff_analysis()
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame
#' @param top_n Number of top DE proteins to show (default 50)
#' @param scale Row scaling method: "row" (Z-score, default) or "none"
#' @param cluster_rows Cluster rows (default TRUE)
#' @param cluster_cols Cluster columns (default TRUE)
#' @param show_rownames Show row names (default TRUE, auto-off if top_n > 80)
#' @param color_palette Color palette for heatmap values
#' @param output_dir Output directory
#' @param project_name Project name for file naming
#' @return pheatmap object (invisible)
#' @export
plot_diff_heatmap <- function(diff_result,
                              ms_data,
                              group_info,
                              top_n = 50,
                              scale = "row",
                              cluster_rows = TRUE,
                              cluster_cols = TRUE,
                              show_rownames = TRUE,
                              color_palette = NULL,
                              output_dir = "output",
                              project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Please install pheatmap: install.packages('pheatmap')")
  }

  # --- Extract top DE proteins ---
  diff_df <- diff_result[diff_result$diff %in% c("UP", "DOWN"), ]
  if (nrow(diff_df) == 0) {
    message(">>> No differentially expressed proteins to plot.")
    return(invisible(NULL))
  }

  # Sort by |logFC| and take top_n
  diff_df <- diff_df[order(-abs(diff_df$logFC)), ]
  diff_df <- utils::head(diff_df, top_n)

  # Get row indices in original data
  row_indices <- as.integer(rownames(diff_df))

  # --- Build expression matrix ---
  contrast <- attr(diff_result, "contrast") %||% "Contrast"
  ref_grp  <- attr(diff_result, "ref_group")
  test_grp <- attr(diff_result, "test_group")

  # Only use samples from the contrast groups (if available)
  if (!is.null(ref_grp) && !is.null(test_grp)) {
    contrast_samples <- group_info$sample_name[
      group_info$user_group %in% c(ref_grp, test_grp)]
    valid_samples <- intersect(contrast_samples, ms_data$sample_names)
  } else {
    valid_samples <- intersect(ms_data$sample_names, group_info$sample_name)
  }
  expr_mat <- as.matrix(ms_data$proteins[row_indices, valid_samples, drop = FALSE])

  # Log2 transform if needed
  expr_mat[expr_mat == 0] <- NA
  if (!isTRUE(ms_data$is_log2)) {
    expr_mat <- log2(expr_mat)
  }

  # Handle Inf/NA for hclust compatibility
  expr_mat[is.infinite(expr_mat)] <- NA
  na_pct <- rowMeans(is.na(expr_mat))
  keep_rows <- na_pct < 0.8
  n_removed <- sum(!keep_rows)
  expr_mat <- expr_mat[keep_rows, , drop = FALSE]
  diff_df <- diff_df[keep_rows, , drop = FALSE]
  if (n_removed > 0) {
    message(sprintf(">>> Removed %d proteins with >80%% missing values.", n_removed))
  }
  if (nrow(expr_mat) == 0) {
    message(">>> All rows have >80% missing values. Cannot plot heatmap.")
    return(invisible(NULL))
  }
  for (i in seq_len(nrow(expr_mat))) {
    row_na <- is.na(expr_mat[i, ])
    if (any(row_na)) {
      row_min <- min(expr_mat[i, !row_na], na.rm = TRUE)
      expr_mat[i, row_na] <- row_min - 1
    }
  }

  # --- Row labels ---
  if ("Label_Name" %in% colnames(diff_df)) {
    row_labels <- diff_df$Label_Name
  } else if ("Gene" %in% colnames(diff_df)) {
    row_labels <- diff_df$Gene
  } else {
    row_labels <- rownames(expr_mat)
  }
  # Handle duplicates
  row_labels <- make.unique(as.character(row_labels), sep = "_")
  rownames(expr_mat) <- row_labels

  # Auto-hide row names if too many
  if (top_n > 80) show_rownames <- FALSE

  # --- Column annotation ---
  grp_vec <- group_info$user_group[match(valid_samples, group_info$sample_name)]
  col_anno <- data.frame(Group = factor(grp_vec), row.names = valid_samples)

  # Annotation colors
  n_groups <- length(unique(grp_vec))
  group_colors <- mspx_colors(n_groups)
  names(group_colors) <- levels(col_anno$Group)
  anno_colors <- list(Group = group_colors)

  # --- Color palette ---
  if (is.null(color_palette)) {
    color_palette <- grDevices::colorRampPalette(
      c("#3C5488", "#4DBBD5", "white", "#E64B35", "#8B0000")
    )(100)
  }

  # --- Plot ---
  fname_base <- file.path(output_dir,
                          paste0(project_name, "_Heatmap_", contrast))

  # Save to PDF (pheatmap is not ggplot, use pdf/dev.off)
  grDevices::pdf(paste0(fname_base, ".pdf"),
                 width = 8,
                 height = max(6, top_n * 0.15 + 2))

  ph <- pheatmap::pheatmap(
    expr_mat,
    scale            = scale,
    cluster_rows     = cluster_rows,
    cluster_cols     = cluster_cols,
    show_rownames    = show_rownames,
    show_colnames    = TRUE,
    annotation_col   = col_anno,
    annotation_colors = anno_colors,
    color            = color_palette,
    border_color     = NA,
    fontsize_row     = if (top_n <= 30) 8 else if (top_n <= 60) 6 else 5,
    fontsize_col     = 9,
    main             = paste("Differential Proteins:", contrast),
    silent           = TRUE
  )

  # Must explicitly draw in PDF device
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)

  grDevices::dev.off()

  # Save source data as CSV
  out_df <- as.data.frame(expr_mat)
  out_df$Gene <- rownames(expr_mat)
  out_df$diff <- as.character(diff_df$diff)
  out_df$logFC <- diff_df$logFC
  utils::write.csv(out_df, paste0(fname_base, ".csv"), row.names = FALSE)

  # Print to screen
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)

  message(sprintf(">>> Heatmap saved: %s (.pdf + .csv)", basename(fname_base)))
  message(sprintf("    Showing %d DE proteins (%d UP, %d DOWN)",
                  nrow(diff_df),
                  sum(diff_df$diff == "UP"),
                  sum(diff_df$diff == "DOWN")))

  invisible(ph)
}


#' Plot top N protein heatmap (by p-value or fold change)
#'
#' A simpler version that shows top N proteins ranked by significance,
#' regardless of UP/DOWN status.
#'
#' @param diff_result Output from run_diff_analysis()
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame
#' @param top_n Number of top proteins to show
#' @param rank_by Ranking criterion: "pvalue" or "fc" (default "pvalue")
#' @param output_dir Output directory
#' @param project_name Project name
#' @return pheatmap object (invisible)
#' @export
plot_top_heatmap <- function(diff_result,
                             ms_data,
                             group_info,
                             top_n = 30,
                             rank_by = "pvalue",
                             output_dir = "output",
                             project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Please install pheatmap: install.packages('pheatmap')")
  }

  # Rank proteins
  if (rank_by == "pvalue") {
    p_col <- if ("adj.P.Val" %in% colnames(diff_result)) "adj.P.Val" else "P.Value"
    diff_sorted <- diff_result[order(diff_result[[p_col]]), ]
  } else {
    diff_sorted <- diff_result[order(-abs(diff_result$logFC)), ]
  }

  diff_top <- utils::head(diff_sorted, top_n)
  row_indices <- as.integer(rownames(diff_top))

  # Build expression matrix
  valid_samples <- intersect(ms_data$sample_names, group_info$sample_name)
  expr_mat <- as.matrix(ms_data$proteins[row_indices, valid_samples, drop = FALSE])
  expr_mat[expr_mat == 0] <- NA
  if (!isTRUE(ms_data$is_log2)) {
    expr_mat <- log2(expr_mat)
  }

  # Handle Inf/NA for hclust compatibility
  expr_mat[is.infinite(expr_mat)] <- NA
  na_pct <- rowMeans(is.na(expr_mat))
  keep_rows <- na_pct <= 0.5
  expr_mat <- expr_mat[keep_rows, , drop = FALSE]
  diff_top <- diff_top[keep_rows, , drop = FALSE]
  if (nrow(expr_mat) == 0) {
    message(">>> All rows have >50% missing values. Cannot plot heatmap.")
    return(invisible(NULL))
  }
  for (i in seq_len(nrow(expr_mat))) {
    row_na <- is.na(expr_mat[i, ])
    if (any(row_na)) {
      row_min <- min(expr_mat[i, !row_na], na.rm = TRUE)
      expr_mat[i, row_na] <- row_min - 1
    }
  }

  # Row labels
  if ("Label_Name" %in% colnames(diff_top)) {
    rownames(expr_mat) <- make.unique(as.character(diff_top$Label_Name), sep = "_")
  }

  # Column annotation
  grp_vec <- group_info$user_group[match(valid_samples, group_info$sample_name)]
  col_anno <- data.frame(Group = factor(grp_vec), row.names = valid_samples)
  n_groups <- length(unique(grp_vec))
  group_colors <- mspx_colors(n_groups)
  names(group_colors) <- levels(col_anno$Group)

  contrast <- attr(diff_result, "contrast") %||% "Contrast"
  fname_base <- file.path(output_dir,
                          paste0(project_name, "_TopHeatmap_", contrast))

  grDevices::pdf(paste0(fname_base, ".pdf"), width = 8,
                 height = max(5, top_n * 0.18 + 2))

  ph <- pheatmap::pheatmap(
    expr_mat,
    scale            = "row",
    cluster_rows     = TRUE,
    cluster_cols     = TRUE,
    show_rownames    = TRUE,
    annotation_col   = col_anno,
    annotation_colors = list(Group = group_colors),
    color            = grDevices::colorRampPalette(
      c("#3C5488", "#4DBBD5", "white", "#E64B35", "#8B0000")
    )(100),
    border_color     = NA,
    fontsize_row     = if (top_n <= 30) 8 else 6,
    main             = paste("Top", top_n, "Proteins:", contrast),
    silent           = TRUE
  )

  grid::grid.newpage()
  grid::grid.draw(ph$gtable)
  grDevices::dev.off()

  # Save CSV
  out_df <- as.data.frame(expr_mat)
  out_df$Gene <- rownames(expr_mat)
  utils::write.csv(out_df, paste0(fname_base, ".csv"), row.names = FALSE)

  # Print to screen
  grid::grid.newpage()
  grid::grid.draw(ph$gtable)

  message(sprintf(">>> Top heatmap saved: %s", basename(fname_base)))
  invisible(ph)
}
