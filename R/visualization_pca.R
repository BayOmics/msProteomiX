# ==============================================================================
# msProteomiX — PCA 可视化
# ==============================================================================

#' 绘制 PCA 主成分分析图
#'
#' 包含 95% 置信椭圆、简短样本标签 (Group-1, Group-2)，
#' 使用 ggrepel 避免标签重叠。
#'
#' @param ms_data MsDataSet 对象
#' @param group_info 分组信息 data.frame
#' @param output_dir 输出目录
#' @param project_name 项目名
#' @param ellipse_level 置信椭圆水平 (默认 0.95)
#' @return ggplot 对象
#' @export
plot_pca <- function(ms_data, group_info,
                     output_dir = "output",
                     project_name = "Project",
                     ellipse_level = 0.95) {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  prot_mat <- as.matrix(ms_data$proteins)
  colnames(prot_mat) <- ms_data$sample_names

  # 匹配分组
  valid_samples <- intersect(ms_data$sample_names, group_info$sample_name)
  grp_vec <- group_info$user_group[match(valid_samples, group_info$sample_name)]
  prot_sub <- prot_mat[, valid_samples, drop = FALSE]

  # Log2 转换 + 缺失值处理 (参考实现: min_val * 0.9)
  prot_sub[prot_sub == 0] <- NA
  keep <- rowSums(!is.na(prot_sub)) > 0
  prot_log <- log2(prot_sub[keep, , drop = FALSE])

  # 用最小值 * 0.9 填补缺失值 (与参考脚本一致)
  min_val <- min(prot_log[!is.infinite(as.matrix(prot_log))], na.rm = TRUE)
  impute_val <- min_val * 0.9
  prot_log[is.na(prot_log)] <- impute_val
  prot_log[is.infinite(as.matrix(prot_log))] <- impute_val

  # Remove zero-variance proteins (constant across all samples after imputation)
  # prcomp(scale.=TRUE) cannot handle these
  row_vars <- apply(prot_log, 1, stats::var, na.rm = TRUE)
  prot_log <- prot_log[row_vars > 0 & !is.na(row_vars), , drop = FALSE]

  # PCA
  pca_result <- stats::prcomp(t(prot_log), center = TRUE, scale. = TRUE)

  # 提取方差解释比例
  var_pct <- round(pca_result$sdev^2 / sum(pca_result$sdev^2) * 100, 1)

  pca_df <- data.frame(
    PC1 = pca_result$x[, 1],
    PC2 = pca_result$x[, 2],
    Sample = valid_samples,
    Group = grp_vec,
    stringsAsFactors = FALSE
  )

  # 生成简短标签 (Group-1, Group-2, ...)
  pca_df <- do.call(rbind, lapply(split(pca_df, pca_df$Group), function(sub) {
    sub$ShortLabel <- paste0(sub$Group, "-", seq_len(nrow(sub)))
    sub
  }))
  rownames(pca_df) <- NULL

  # 检查每组样本数 (椭圆需要 >= 3)
  grp_counts <- table(pca_df$Group)

  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = Group)) +
    ggplot2::geom_point(size = 4, alpha = 0.9) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = ShortLabel),
      size = 3.5, show.legend = FALSE,
      max.overlaps = 50, box.padding = 0.5
    ) +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.title = ggplot2::element_text(size = 12, face = "bold"),
      legend.position = "right"
    ) +
    ggplot2::labs(
      title = paste0("PCA Analysis (with ", round(ellipse_level * 100), "% CI Ellipses)"),
      x = sprintf("PC1 (%s%%)", var_pct[1]),
      y = sprintf("PC2 (%s%%)", var_pct[2])
    )

  # 95% 置信椭圆 (填充色, 需要 >= 3 per group)
  if (all(grp_counts >= 3)) {
    p <- p +
      ggplot2::stat_ellipse(
        ggplot2::aes(fill = Group),
        geom = "polygon", level = ellipse_level,
        alpha = 0.15, show.legend = FALSE
      ) +
      ggplot2::scale_fill_brewer(palette = "Set1")
  } else {
    # 部分组 < 3, 只给 >= 3 的组画椭圆
    eligible <- names(grp_counts[grp_counts >= 3])
    if (length(eligible) > 0) {
      ellipse_data <- pca_df[pca_df$Group %in% eligible, ]
      p <- p +
        ggplot2::stat_ellipse(
          data = ellipse_data,
          ggplot2::aes(fill = Group),
          geom = "polygon", level = ellipse_level,
          alpha = 0.15, show.legend = FALSE
        ) +
        ggplot2::scale_fill_brewer(palette = "Set1")
    }
  }

  grp_str <- paste(unique(grp_vec), collapse = "_")
  if (nchar(grp_str) > 40) grp_str <- paste0(length(unique(grp_vec)), "Groups")
  suppressWarnings(
    save_plot_and_data(p, pca_df, project_name, paste0("PCA_", grp_str),
                       output_dir = output_dir)
  )
  p
}
