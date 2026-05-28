# ==============================================================================
# msProteomiX — 相关性热图可视化
# ==============================================================================

#' 绘制样品相关性热图
#'
#' @param ms_data MsDataSet 对象
#' @param group_info 分组信息 data.frame
#' @param method 相关系数方法: "pearson" 或 "spearman"
#' @param output_dir 输出目录
#' @param project_name 项目名
#' @return ggplot 对象
#' @export
plot_corr_heatmap <- function(ms_data, group_info,
                               method = "pearson",
                               output_dir = "output",
                               project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  prot_mat <- as.matrix(ms_data$proteins)
  colnames(prot_mat) <- ms_data$sample_names

  # 匹配分组
  valid_samples <- intersect(ms_data$sample_names, group_info$sample_name)
  prot_sub <- prot_mat[, valid_samples, drop = FALSE]

  # Log2 转换
  prot_sub[prot_sub == 0] <- NA
  prot_log <- log2(prot_sub)

  # 计算相关性矩阵
  cor_mat <- stats::cor(prot_log, use = "pairwise.complete.obs", method = method)

  # 转为长格式
  cor_df <- reshape2::melt(cor_mat)
  colnames(cor_df) <- c("Sample1", "Sample2", "Correlation")

  p <- ggplot2::ggplot(cor_df, ggplot2::aes(x = Sample1, y = Sample2, fill = Correlation)) +
    ggplot2::geom_tile(color = "white") +
    ggplot2::geom_text(ggplot2::aes(label = round(Correlation, 3)), size = 3) +
    ggplot2::scale_fill_gradient2(low = "#0073C2FF", mid = "white", high = "#CD534CFF",
                                   midpoint = median(cor_df$Correlation, na.rm = TRUE),
                                   limits = c(min(cor_df$Correlation, na.rm = TRUE), 1)) +
    ggplot2::theme_minimal() +
    ggplot2::labs(title = paste("Correlation Heatmap (", method, ")", sep = ""),
                  x = NULL, y = NULL) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = ggplot2::element_text(size = 10)
    )

  grp_str <- paste(unique(group_info$user_group), collapse = "_")
  save_plot_and_data(p, as.data.frame(cor_mat), project_name,
                     paste0("CorrHeatmap_", grp_str), output_dir = output_dir)
  p
}
