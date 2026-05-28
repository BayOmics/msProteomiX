# ==============================================================================
# msProteomiX — 火山图可视化
# ==============================================================================

#' 绘制火山图
#'
#' 统一的火山图函数，替代原来的 6 个火山图脚本变体。
#'
#' @param diff_result run_diff_analysis() 返回的 data.frame
#' @param label_top 标注前 N 个差异蛋白 (默认 20)
#' @param p_type P 值类型: "raw" 或 "adj"
#' @param output_dir 输出目录
#' @param project_name 项目名
#' @return ggplot 对象
#' @export
plot_volcano <- function(diff_result,
                          label_top = 20,
                          p_type = NULL,
                          output_dir = "output",
                          project_name = "Project") {
  ensure_output_dir(output_dir)

  # 从属性获取元数据
  contrast <- attr(diff_result, "contrast") %||% "Contrast"
  if (is.null(p_type)) p_type <- attr(diff_result, "p_type") %||% "raw"
  fc_cutoff <- attr(diff_result, "fc_cutoff") %||% 1
  p_cutoff <- attr(diff_result, "p_cutoff") %||% 0.05

  # 选择 P 值列
  p_col <- if (p_type == "adj" && "adj.P.Val" %in% colnames(diff_result)) {
    "adj.P.Val"
  } else {
    "P.Value"
  }

  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("Please install ggrepel: install.packages('ggrepel')")
  }

  # 提取标注标签
  label_col <- if ("Label_Name" %in% colnames(diff_result)) "Label_Name" else NULL

  # Top labels
  top_labels <- data.frame()
  if (!is.null(label_col) && label_top > 0) {
    # 按 |logFC| 降序排列 (标注变化最大的蛋白)
    top_up <- diff_result[diff_result$diff == "UP", ]
    top_up <- top_up[order(-abs(top_up$logFC)), ]
    top_up <- utils::head(top_up, label_top)

    top_down <- diff_result[diff_result$diff == "DOWN", ]
    top_down <- top_down[order(-abs(top_down$logFC)), ]
    top_down <- utils::head(top_down, label_top)

    top_labels <- rbind(top_up, top_down)
  }

  # 统计
  n_up <- sum(diff_result$diff == "UP", na.rm = TRUE)
  n_down <- sum(diff_result$diff == "DOWN", na.rm = TRUE)

  # 绘图
  p <- ggplot2::ggplot(diff_result,
                       ggplot2::aes(x = logFC, y = -log10(.data[[p_col]]), col = diff)) +
    ggplot2::geom_point(alpha = 0.6, size = 1.5) +
    ggplot2::scale_color_manual(values = c("UP" = "#CD534CFF", "DOWN" = "#0073C2FF", "NO" = "grey70")) +
    ggplot2::geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), lty = 2, lwd = 0.4, color = "grey40") +
    ggplot2::geom_hline(yintercept = -log10(p_cutoff), lty = 2, lwd = 0.4, color = "grey40") +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste("Volcano:", contrast),
      subtitle = sprintf("Up: %d | Down: %d", n_up, n_down),
      x = "log2 Fold Change",
      y = sprintf("-log10(%s)", ifelse(p_type == "adj", "adj.P", "P.Value"))
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = ggplot2::element_text(hjust = 0.5)
    )

  # 添加标签
  if (nrow(top_labels) > 0 && !is.null(label_col)) {
    p <- p + ggrepel::geom_text_repel(
      data = top_labels,
      ggplot2::aes(label = .data[[label_col]]),
      size = 3, max.overlaps = 50, show.legend = FALSE, color = "black"
    )
  }

  print(p)

  # 保存
  fname_base <- file.path(output_dir, paste0(project_name, "_Volcano_", contrast))
  ggplot2::ggsave(paste0(fname_base, ".pdf"), p, width = 8, height = 7)
  utils::write.csv(diff_result, paste0(fname_base, ".csv"), row.names = FALSE)
  message(paste("  Results saved:", fname_base))

  p
}


