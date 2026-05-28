# ==============================================================================
# msProteomiX — PCA 可视化
# ==============================================================================

#' 绘制 PCA 主成分分析图
#'
#' @param ms_data MsDataSet 对象
#' @param group_info 分组信息 data.frame
#' @param output_dir 输出目录
#' @param project_name 项目名
#' @return ggplot 对象
#' @export
plot_pca <- function(ms_data, group_info,
                     output_dir = "output",
                     project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  prot_mat <- as.matrix(ms_data$proteins)
  colnames(prot_mat) <- ms_data$sample_names

  # 匹配分组
  valid_samples <- intersect(ms_data$sample_names, group_info$sample_name)
  grp_vec <- group_info$user_group[match(valid_samples, group_info$sample_name)]
  prot_sub <- prot_mat[, valid_samples, drop = FALSE]

  # Log2 转换 + 移除全 NA 行
  prot_sub[prot_sub == 0] <- NA
  prot_log <- log2(prot_sub)

  # 移除缺失太多的行 (>50%)
  keep <- rowSums(!is.na(prot_log)) >= ncol(prot_log) * 0.5
  prot_log <- prot_log[keep, , drop = FALSE]

  # 插补剩余 NA (行均值)
  for (i in seq_len(nrow(prot_log))) {
    row_mean <- mean(prot_log[i, ], na.rm = TRUE)
    prot_log[i, is.na(prot_log[i, ])] <- row_mean
  }

  # PCA
  pca_result <- stats::prcomp(t(prot_log), center = TRUE, scale. = TRUE)

  # 提取方差解释比例
  var_pct <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 1)

  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Sample = valid_samples,
    Group = grp_vec
  )

  n_groups <- length(unique(grp_vec))
  colors <- mspx_colors(n_groups)

  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
    ggplot2::geom_point(size = 4, alpha = 0.8) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = "PCA - Protein Quantification",
      x = sprintf("PC1 (%s%%)", var_pct[1]),
      y = sprintf("PC2 (%s%%)", var_pct[2])
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )

  # stat_ellipse requires >=3 points per group; skip if any group is too small
  grp_counts <- table(grp_vec)
  if (all(grp_counts >= 3)) {
    p <- p + ggplot2::stat_ellipse(level = 0.68, linetype = "dashed", show.legend = FALSE)
  }

  grp_str <- paste(unique(grp_vec), collapse = "_")
  save_plot_and_data(p, pca_df, project_name, paste0("PCA_", grp_str),
                     output_dir = output_dir)
  p
}
