# ==============================================================================
# msProteomiX — 相关性热图可视化
# ==============================================================================

#' 绘制样品相关性热图 (R²)
#'
#' 使用 R² (决定系数) 绘制样品间相关性热图。
#' 颜色: 白色 (低 R²) → 红色 (高 R², ≥ 0.8)。
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

  # 计算相关性矩阵 → R²
  cor_mat <- stats::cor(prot_log, use = "pairwise.complete.obs", method = method)
  r2_mat <- cor_mat^2

  # 转为长格式
  r2_df <- reshape2::melt(r2_mat)
  colnames(r2_df) <- c("Sample1", "Sample2", "R2")

  # 保持样本顺序
  r2_df$Sample1 <- factor(r2_df$Sample1, levels = valid_samples)
  r2_df$Sample2 <- factor(r2_df$Sample2, levels = rev(valid_samples))

  # 自适应字体大小 (样本多时缩小)
  n_samples <- length(valid_samples)
  text_size <- if (n_samples <= 10) 3.5 else if (n_samples <= 20) 2.5 else 1.8
  axis_size <- if (n_samples <= 10) 10 else if (n_samples <= 20) 8 else 6

  # 只在下三角 (含对角线) 显示数字, 避免重叠
  r2_df$label <- ifelse(
    as.integer(r2_df$Sample1) <= (n_samples - as.integer(r2_df$Sample2) + 1),
    sprintf("%.2f", r2_df$R2),
    ""
  )

  p <- ggplot2::ggplot(r2_df, ggplot2::aes(x = Sample1, y = Sample2, fill = R2)) +
    ggplot2::geom_tile(color = "grey90", linewidth = 0.3) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = text_size, color = "black") +
    ggplot2::scale_fill_gradientn(
      colors = c("white", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"),
      values = scales::rescale(c(0.5, 0.7, 0.8, 0.9, 1.0)),
      limits = c(min(r2_df$R2, na.rm = TRUE), 1),
      name = expression(R^2)
    ) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = bquote("Sample Correlation Heatmap (" ~ R^2 ~ "," ~ .(method) ~ ")"),
      x = NULL, y = NULL
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = axis_size),
      axis.text.y = ggplot2::element_text(size = axis_size),
      panel.grid = ggplot2::element_blank(),
      legend.title = ggplot2::element_text(face = "bold")
    )

  # 自适应 PDF 尺寸
  pdf_w <- max(8, n_samples * 0.4 + 2)
  pdf_h <- max(7, n_samples * 0.35 + 2)

  grp_str <- paste(unique(group_info$user_group), collapse = "_")
  fname <- paste0("CorrHeatmap_", grp_str)
  base_name <- file.path(output_dir, paste0(project_name, "_", fname))
  utils::write.csv(as.data.frame(r2_mat), paste0(base_name, ".csv"))
  ggplot2::ggsave(paste0(base_name, ".pdf"), p, width = pdf_w, height = pdf_h)
  print(p)
  message(sprintf("  \u2705 Generated: %s", fname))

  p
}
