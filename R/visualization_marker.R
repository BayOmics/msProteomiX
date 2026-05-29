# ==============================================================================
# msProteomiX — Marker Gene Expression Visualization
# ==============================================================================

#' Plot marker gene boxplot
#'
#' Displays expression levels of user-specified marker genes across groups
#' as boxplots with individual data points.
#'
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame
#' @param markers Character vector of gene symbols to plot
#' @param output_dir Output directory
#' @param project_name Project name
#' @param ncol Number of columns in facet layout (default auto)
#' @return ggplot object or NULL
#' @export
plot_marker_boxplot <- function(ms_data, group_info,
                                 markers,
                                 output_dir = "output",
                                 project_name = "Project",
                                 ncol = NULL) {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  prot_mat <- ms_data$proteins
  pinfo <- ms_data$protein_info

  # Find gene column
  gene_col <- if ("Gene" %in% colnames(pinfo)) "Gene"
              else if ("Label_Name" %in% colnames(pinfo)) "Label_Name"
              else NULL
  if (is.null(gene_col)) {
    message(">>> Cannot find Gene column in protein_info.")
    return(invisible(NULL))
  }

  # Map gene to first symbol (handle multi-gene like "GENE1;GENE2")
  gene_first <- sub(";.*$", "", as.character(pinfo[[gene_col]]))

  # Match markers to protein rows
  found_markers <- c()
  plot_data_list <- list()

  for (marker in markers) {
    idx <- which(gene_first == marker)
    if (length(idx) == 0) {
      # Try case-insensitive
      idx <- which(tolower(gene_first) == tolower(marker))
    }
    if (length(idx) == 0) {
      message(sprintf("  >>> Marker '%s' not found in data.", marker))
      next
    }
    # If multiple rows match, use first
    row_idx <- idx[1]
    found_markers <- c(found_markers, marker)

    expr_vals <- as.numeric(prot_mat[row_idx, ])
    sample_names <- colnames(prot_mat)

    # Match samples to groups
    matched <- match(sample_names, group_info$sample_name)
    groups <- group_info$user_group[matched]

    # Handle is_log2 flag
    if (isTRUE(ms_data$is_log2)) {
      ylab <- "Expression (log2)"
    } else {
      # Apply log2 for visualization
      expr_vals <- log2(expr_vals + 1)
      ylab <- "Expression (log2)"
    }

    df <- data.frame(
      Sample = sample_names,
      Group = groups,
      Expression = expr_vals,
      Gene = marker,
      stringsAsFactors = FALSE
    )
    df <- df[!is.na(df$Group), ]
    plot_data_list[[marker]] <- df
  }

  if (length(plot_data_list) == 0) {
    message(">>> No markers found in the dataset.")
    return(invisible(NULL))
  }

  plot_df <- do.call(rbind, plot_data_list)
  plot_df$Gene <- factor(plot_df$Gene, levels = found_markers)

  # Auto ncol
  n_markers <- length(found_markers)
  if (is.null(ncol)) {
    ncol <- min(4, n_markers)
  }

  # Color palette
  n_groups <- length(unique(plot_df$Group))
  colors <- .marker_palette(n_groups)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Group, y = Expression, fill = Group)) +
    ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA, width = 0.6) +
    ggplot2::geom_jitter(width = 0.15, size = 1.5, alpha = 0.8,
                          ggplot2::aes(color = Group), show.legend = FALSE) +
    ggplot2::facet_wrap(~ Gene, scales = "free_y", ncol = ncol) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(title = "Marker Gene Expression",
                  x = NULL, y = ylab) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      strip.text = ggplot2::element_text(face = "bold", size = 11),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )

  # Adaptive figure size
  nrow_facet <- ceiling(n_markers / ncol)
  fig_w <- max(6, ncol * 2.5)
  fig_h <- max(4, nrow_facet * 3 + 1)

  save_plot_and_data(p, plot_df, project_name, "Marker_Expression",
                     output_dir = output_dir, width = fig_w, height = fig_h)
  p
}


#' Marker boxplot color palette
#' @keywords internal
.marker_palette <- function(n) {
  base_colors <- c("#E64B35", "#4DBBD5", "#00A087", "#F39B7F",
                    "#8491B4", "#91D1C2", "#DC6B5A", "#7E6148")
  if (n <= length(base_colors)) {
    return(base_colors[seq_len(n)])
  }
  grDevices::colorRampPalette(base_colors)(n)
}
