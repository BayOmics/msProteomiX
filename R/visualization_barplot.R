# ==============================================================================
# msProteomiX — 柱状图可视化
# ==============================================================================

#' 绘制蛋白/肽段/PSM 鉴定数量柱状图
#'
#' @param ms_data MsDataSet 对象
#' @param group_info 分组信息 data.frame
#' @param target 绘图目标: "protein_groups", "peptide", "psm", 或 "all"
#' @param output_dir 输出目录
#' @param project_name 项目名 (用于文件命名)
#' @return ggplot 对象 (或 list)
#' @export
plot_id_barplot <- function(ms_data, group_info,
                             target = "all",
                             output_dir = "output",
                             project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  # TODO: support target = "peptide", "psm", "all" once peptide/PSM counts are

  # available in count_ids()
  id_df <- count_ids(ms_data, group_info)
  if (nrow(id_df) == 0) stop("No valid data for barplot.")

  # 统计每组均值和标准差
  stats_df <- id_df %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(
      mean_val = mean(protein_groups, na.rm = TRUE),
      sd_val   = stats::sd(protein_groups, na.rm = TRUE),
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
    ggplot2::labs(title = "Protein Groups Identification",
                  x = "Group", y = "Number of Protein Groups") +
    ggplot2::coord_cartesian(ylim = c(0, y_max))

  save_plot_and_data(p, stats_df, project_name, "Identification_ProteinGroups",
                     output_dir = output_dir)
  p
}
