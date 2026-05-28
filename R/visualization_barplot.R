# ==============================================================================
# msProteomiX — 柱状图可视化
# ==============================================================================

#' 绘制蛋白/肽段/PSM 鉴定数量柱状图
#'
#' 根据 target 参数绘制 Protein Groups、Peptide、PSM 的鉴定数量柱状图。
#' target = "all" 时同时输出三张图。
#'
#' @param ms_data MsDataSet 对象
#' @param group_info 分组信息 data.frame
#' @param target 绘图目标: "protein_groups", "peptide", "psm", 或 "all"
#' @param output_dir 输出目录
#' @param project_name 项目名 (用于文件命名)
#' @return ggplot 对象 (target="all" 时返回 list)
#' @export
plot_id_barplot <- function(ms_data, group_info,
                             target = "all",
                             output_dir = "output",
                             project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  id_df <- count_ids(ms_data, group_info)
  if (nrow(id_df) == 0) stop("No valid data for barplot.")

  # 定义绘图目标
  all_targets <- list(
    protein_groups = list(col = "protein_groups", title = "Protein Groups", y_lab = "Number of Protein Groups"),
    peptide        = list(col = "peptide",        title = "Peptides",       y_lab = "Number of Peptides"),
    psm            = list(col = "psm",            title = "PSMs",           y_lab = "Number of PSMs")
  )

  if (target == "all") {
    targets <- names(all_targets)
  } else {
    targets <- target
  }

  plots <- list()

  for (tgt in targets) {
    if (!tgt %in% names(all_targets)) next
    info <- all_targets[[tgt]]

    # 跳过全为 0 的指标
    if (all(id_df[[info$col]] == 0, na.rm = TRUE)) {
      message(sprintf(">>> %s: all zero, skipping.", info$title))
      next
    }

    # 统计每组均值和标准差
    stats_df <- id_df %>%
      dplyr::group_by(group) %>%
      dplyr::summarise(
        mean_val = mean(.data[[info$col]], na.rm = TRUE),
        sd_val   = stats::sd(.data[[info$col]], na.rm = TRUE),
        .groups  = "drop"
      )
    # 单样本组的 sd 为 NA，设为 0
    stats_df$sd_val[is.na(stats_df$sd_val)] <- 0

    y_max <- max(stats_df$mean_val + stats_df$sd_val, na.rm = TRUE) * 1.15

    p <- ggplot2::ggplot(stats_df, ggplot2::aes(x = group, y = mean_val)) +
      ggplot2::geom_bar(stat = "identity", position = ggplot2::position_dodge(),
                        fill = "#87CEFA", color = "black", width = 0.7) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin = mean_val - sd_val,
                                           ymax = mean_val + sd_val),
                             width = 0.2, position = ggplot2::position_dodge(0.9)) +
      ggplot2::geom_text(ggplot2::aes(label = round(mean_val, 0)),
                         vjust = -1.5, size = 4) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid.major = ggplot2::element_blank(),
        axis.line = ggplot2::element_line(colour = "black"),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12),
        axis.title = ggplot2::element_text(size = 14, face = "bold"),
        plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
      ) +
      ggplot2::labs(title = paste(info$title, "Identification"),
                    x = "Group", y = info$y_lab) +
      ggplot2::coord_cartesian(ylim = c(0, y_max))

    fname <- paste0("Identification_", info$title)
    save_plot_and_data(p, stats_df, project_name, fname, output_dir = output_dir)

    message(sprintf(">>> %s barplot saved.", info$title))
    print(p)
    plots[[tgt]] <- p
  }

  if (length(plots) == 1) return(plots[[1]])
  invisible(plots)
}
