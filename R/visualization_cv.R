# ==============================================================================
# msProteomiX — CV 变异系数可视化
# ==============================================================================

#' 绘制 CV 变异系数箱线图
#'
#' @param ms_data MsDataSet 对象
#' @param group_info 分组信息 data.frame
#' @param output_dir 输出目录
#' @param project_name 项目名
#' @return ggplot 对象
#' @export
plot_cv_boxplot <- function(ms_data, group_info,
                             output_dir = "output",
                             project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  cv_df <- .calc_cv_by_group(ms_data, group_info)
  if (is.null(cv_df) || nrow(cv_df) == 0) {
    message("  Not enough replicates for CV calculation.")
    return(invisible(NULL))
  }

  n_groups <- length(unique(cv_df$Group))

  # 中位数标签
  medians <- cv_df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(MedianCV = stats::median(CV, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Label = paste0(round(MedianCV * 100, 1), "%"))

  p <- ggplot2::ggplot(cv_df, ggplot2::aes(x = Group, y = CV, fill = Group)) +
    ggplot2::geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    ggplot2::geom_text(data = medians,
                       ggplot2::aes(x = Group, y = MedianCV, label = Label),
                       vjust = -0.8, fontface = "bold", size = 4, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = mspx_colors(n_groups)) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::coord_cartesian(ylim = c(0, min(max(cv_df$CV, na.rm = TRUE) * 1.1, 1.5))) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Coefficient of Variation (CV)", y = "CV", x = "Group") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  save_plot_and_data(p, cv_df, project_name, "CV_Boxplot", output_dir = output_dir)
  p
}


#' 绘制 CV 变异系数小提琴图
#'
#' @param ms_data MsDataSet 对象
#' @param group_info 分组信息 data.frame
#' @param output_dir 输出目录
#' @param project_name 项目名
#' @return ggplot 对象
#' @export
plot_cv_violin <- function(ms_data, group_info,
                            output_dir = "output",
                            project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  cv_df <- .calc_cv_by_group(ms_data, group_info)
  if (is.null(cv_df) || nrow(cv_df) == 0) return(invisible(NULL))

  n_groups <- length(unique(cv_df$Group))

  # 中位数标签
  medians <- cv_df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(MedianCV = stats::median(CV, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Label = paste0(round(MedianCV * 100, 1), "%"))

  p <- ggplot2::ggplot(cv_df, ggplot2::aes(x = Group, y = CV, fill = Group)) +
    ggplot2::geom_violin(alpha = 0.6, trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    ggplot2::geom_text(data = medians,
                       ggplot2::aes(x = Group, y = MedianCV, label = Label),
                       vjust = -0.8, fontface = "bold", size = 4, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = mspx_colors(n_groups)) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::coord_cartesian(ylim = c(0, min(max(cv_df$CV, na.rm = TRUE) * 1.1, 1.5))) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "CV Distribution (Violin)", y = "CV", x = "Group") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )

  save_plot_and_data(p, cv_df, project_name, "CV_Violin", output_dir = output_dir)
  p
}


#' @keywords internal
.calc_cv_by_group <- function(ms_data, group_info) {
  prot_mat <- as.matrix(ms_data$proteins)
  unique_grps <- unique(group_info$user_group)
  cv_list <- list()

  for (grp in unique_grps) {
    samps <- group_info$sample_name[group_info$user_group == grp]
    col_idx <- match(samps, ms_data$sample_names)
    col_idx <- col_idx[!is.na(col_idx)]

    if (length(col_idx) < 2) next

    sub_dat <- prot_mat[, col_idx, drop = FALSE]
    sub_dat[sub_dat == 0] <- NA
    means <- rowMeans(sub_dat, na.rm = TRUE)
    sds <- apply(sub_dat, 1, stats::sd, na.rm = TRUE)
    cvs <- sds / means
    valid_cv <- cvs[!is.na(cvs) & cvs < 1.5]
    if (length(valid_cv) > 10) {
      cv_list[[grp]] <- data.frame(Group = grp, CV = valid_cv)
    }
  }

  if (length(cv_list) == 0) return(NULL)
  do.call(rbind, cv_list)
}
