# ==============================================================================
# msProteomiX — AP-MS 可视化模块
# ==============================================================================

#' AP-MS 相对丰度曲线图
#'
#' 绘制蛋白在不同时间点/条件下的相对丰度变化曲线。
#' 可指定特定蛋白，或绘制所有显著蛋白（生成多页 PDF）。
#'
#' @param rel_abundance calc_relative_abundance() 返回的 data.frame
#' @param proteins 要绘制的蛋白 ID 向量 (NULL = 全部)
#' @param max_plots 最多绘制多少个蛋白 (默认 50, 防止 PDF 过大)
#' @param project_name 项目名 (用于输出文件名)
#' @return ggplot 对象 (最后一个)
#' @export
plot_abundance_curve <- function(rel_abundance, proteins = NULL,
                                  max_plots = 50, project_name = "APMS") {
  if (!is.null(proteins)) {
    rel_abundance <- rel_abundance[rel_abundance$Protein %in% proteins, ]
  }

  unique_proteins <- unique(rel_abundance$Protein)
  if (length(unique_proteins) > max_plots) {
    message(sprintf("  \u26a0\ufe0f \u622a\u53d6\u524d %d \u4e2a\u86cb\u767d (\u5171 %d \u4e2a)",
                    max_plots, length(unique_proteins)))
    unique_proteins <- unique_proteins[seq_len(max_plots)]
    rel_abundance <- rel_abundance[rel_abundance$Protein %in% unique_proteins, ]
  }

  # 保存多页 PDF
  out_file <- file.path("output", paste0(project_name, "_Abundance_Curves"))
  if (!dir.exists("output")) dir.create("output")

  pdf_file <- paste0(out_file, ".pdf")
  grDevices::pdf(pdf_file, width = 8, height = 6)
  on.exit(grDevices::dev.off(), add = TRUE)

  last_plot <- NULL
  for (prot in unique_proteins) {
    df_single <- rel_abundance[rel_abundance$Protein == prot, ]
    label <- df_single$Label_Name[1]
    gene <- df_single$Gene[1]
    title <- if (!is.na(gene) && gene != label) {
      sprintf("%s (%s)", gene, prot)
    } else {
      prot
    }

    # 保持组顺序
    df_single$Group <- factor(df_single$Group, levels = unique(df_single$Group))

    p <- ggplot2::ggplot(df_single, ggplot2::aes(x = Group, y = RelAbundance, group = 1)) +
      ggplot2::geom_line(linewidth = 1.2, color = "#2563EB") +
      ggplot2::geom_point(size = 3, color = "#2563EB") +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = RelAbundance, ymax = RelAbundance + SD / max(Mean, 1e-10)),
        width = 0.3, linewidth = 0.8, color = "#2563EB"
      ) +
      ggplot2::coord_cartesian(ylim = c(0, 1.3)) +
      ggplot2::labs(
        title = title,
        x = NULL,
        y = "Relative Abundance"
      ) +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 12),
        plot.title = ggplot2::element_text(face = "bold", size = 14)
      )

    print(p)
    last_plot <- p
  }

  message(sprintf("  \u2705 \u5df2\u751f\u6210: %s (%d \u4e2a\u86cb\u767d)", basename(pdf_file), length(unique_proteins)))
  invisible(last_plot)
}


#' AP-MS Dot Plot (核心可视化)
#'
#' 经典 AP-MS 气泡图：X 轴为时间点/条件，Y 轴为基因名，
#' 气泡大小表示相对丰度，颜色表示最大丰度所在的组。
#'
#' @param rel_abundance calc_relative_abundance() 返回的 data.frame
#' @param sig_proteins 显著蛋白 ID 向量 (NULL = 使用全部)
#' @param project_name 项目名
#' @return ggplot 对象
#' @export
plot_apms_dotplot <- function(rel_abundance, sig_proteins = NULL,
                               project_name = "APMS") {
  if (!is.null(sig_proteins)) {
    rel_abundance <- rel_abundance[rel_abundance$Protein %in% sig_proteins, ]
  }

  if (nrow(rel_abundance) == 0) {
    warning("\u65e0\u6570\u636e\u53ef\u7ed8\u5236")
    return(invisible(NULL))
  }

  # 按最大丰度组排序
  rel_abundance$Group <- factor(rel_abundance$Group, levels = unique(rel_abundance$Group))

  # 用 Label_Name 作为 Y 轴
  # 按 MaxGroup 排序蛋白
  prot_order <- rel_abundance[!duplicated(rel_abundance$Protein), ]
  prot_order <- prot_order[order(prot_order$MaxGroup), ]
  # 用 make.unique 处理重复的 Label_Name
  unique_labels <- make.unique(prot_order$Label_Name)
  label_map <- setNames(unique_labels, prot_order$Protein)
  rel_abundance$DisplayLabel <- label_map[rel_abundance$Protein]
  rel_abundance$DisplayLabel <- factor(
    rel_abundance$DisplayLabel,
    levels = rev(unique_labels)
  )

  n_proteins <- length(unique(rel_abundance$Protein))
  plot_height <- max(8, n_proteins * 0.3)  # 动态高度

  # 配色
  groups <- levels(rel_abundance$Group)
  n_colors <- length(groups)
  dot_colors <- mspx_colors(n_colors)

  p <- ggplot2::ggplot(
    rel_abundance,
    ggplot2::aes(x = Group, y = DisplayLabel,
                 size = RelAbundance, color = MaxGroup)
  ) +
    ggplot2::geom_point() +
    ggplot2::scale_size_continuous(
      range = c(1, 8),
      name = "Relative\nAbundance"
    ) +
    ggplot2::scale_color_manual(
      values = dot_colors,
      name = "Peak Group"
    ) +
    ggplot2::labs(
      title = paste0("AP-MS Dot Plot - ", project_name),
      x = NULL, y = NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 9),
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      legend.position = "right"
    )

  save_plot_and_data(
    plot_obj     = p,
    data_df      = rel_abundance,
    project_name = project_name,
    suffix       = "DotPlot",
    width        = 10,
    height       = min(plot_height, 40)
  )

  p
}


#' AP-MS 热图
#'
#' 显著蛋白的相对丰度热图 (用 pheatmap)。
#' 行按峰值时间点分组排序。
#'
#' @param rel_abundance calc_relative_abundance() 返回的 data.frame
#' @param sig_proteins 显著蛋白 ID 向量 (NULL = 使用全部)
#' @param project_name 项目名
#' @return pheatmap 对象
#' @export
plot_apms_heatmap <- function(rel_abundance, sig_proteins = NULL,
                               project_name = "APMS") {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("\u274c \u9700\u8981\u5b89\u88c5 pheatmap \u5305: install.packages('pheatmap')")
  }

  if (!is.null(sig_proteins)) {
    rel_abundance <- rel_abundance[rel_abundance$Protein %in% sig_proteins, ]
  }

  if (nrow(rel_abundance) == 0) {
    warning("\u65e0\u6570\u636e\u53ef\u7ed8\u5236")
    return(invisible(NULL))
  }

  groups <- unique(rel_abundance$Group)

  # 构建宽格式矩阵
  mat <- matrix(NA, nrow = length(unique(rel_abundance$Protein)),
                ncol = length(groups))
  colnames(mat) <- groups
  unique_prots <- unique(rel_abundance$Protein)

  for (j in seq_along(groups)) {
    sub <- rel_abundance[rel_abundance$Group == groups[j], ]
    mat[, j] <- sub$RelAbundance[match(unique_prots, sub$Protein)]
  }

  # 行名用基因名
  label_names <- rel_abundance$Label_Name[match(unique_prots, rel_abundance$Protein)]
  rownames(mat) <- make.unique(label_names)

  # 按最大丰度组排序
  max_group <- rel_abundance$MaxGroup[match(unique_prots, rel_abundance$Protein)]
  mat <- mat[order(max_group), , drop = FALSE]

  mat[is.na(mat)] <- 0

  n_proteins <- nrow(mat)
  plot_height <- max(8, n_proteins * 0.2)

  out_file <- file.path("output", paste0(project_name, "_Heatmap"))
  if (!dir.exists("output")) dir.create("output")

  # 保存数据
  utils::write.csv(mat, paste0(out_file, ".csv"))

  grDevices::pdf(paste0(out_file, ".pdf"), width = 8, height = min(plot_height, 40))
  ph <- pheatmap::pheatmap(
    mat,
    cluster_rows = FALSE,
    cluster_cols = FALSE,
    color = grDevices::colorRampPalette(c("#F0F4FF", "#3B82F6", "#1E3A8A"))(100),
    border_color = NA,
    main = paste0("AP-MS Heatmap - ", project_name),
    fontsize_row = max(6, 12 - n_proteins * 0.05),
    fontsize_col = 12,
    show_rownames = n_proteins <= 100
  )
  grDevices::dev.off()

  message(sprintf("  \u2705 \u5df2\u751f\u6210: %s", basename(out_file)))
  ph
}


#' Mfuzz 聚类结果可视化
#'
#' 绘制 Mfuzz 软聚类的中心折线图和成员热图。
#'
#' @param mfuzz_result run_mfuzz_cluster() 返回的 list
#' @param project_name 项目名
#' @return ggplot 对象
#' @export
plot_mfuzz <- function(mfuzz_result, project_name = "APMS") {
  if (!requireNamespace("Mfuzz", quietly = TRUE)) {
    stop("\u274c \u9700\u8981\u5b89\u88c5 Mfuzz \u5305")
  }

  cl   <- mfuzz_result$cl
  eset <- mfuzz_result$eset
  n_cl <- mfuzz_result$n_clusters

  # 使用 Mfuzz 自带的绘图
  out_file <- file.path("output", paste0(project_name, "_Mfuzz_Clusters"))
  if (!dir.exists("output")) dir.create("output")

  grDevices::pdf(paste0(out_file, ".pdf"), width = 12, height = 3 * ceiling(n_cl / 2))
  Mfuzz::mfuzz.plot2(eset, cl = cl,
                     mfrow = c(ceiling(n_cl / 2), 2),
                     time.labels = colnames(Biobase::exprs(eset)),
                     centre = TRUE)
  grDevices::dev.off()

  # 同时用 ggplot 画中心曲线
  centers <- cl$centers
  groups <- colnames(centers)
  n_groups <- ncol(centers)

  center_df <- data.frame(
    Cluster = rep(seq_len(n_cl), each = n_groups),
    Group   = rep(groups, n_cl),
    Value   = as.vector(t(centers)),
    stringsAsFactors = FALSE
  )
  center_df$Group <- factor(center_df$Group, levels = groups)
  center_df$Cluster <- factor(center_df$Cluster)

  # 计算每个 cluster 的蛋白数
  cluster_assign <- cl$cluster
  cluster_sizes <- table(cluster_assign)
  center_df$ClusterLabel <- paste0("Cluster ", center_df$Cluster,
                                    " (n=", cluster_sizes[as.character(center_df$Cluster)], ")")

  p <- ggplot2::ggplot(center_df,
    ggplot2::aes(x = Group, y = Value, group = 1)) +
    ggplot2::geom_line(linewidth = 1.5, color = "#2563EB") +
    ggplot2::geom_point(size = 3, color = "#2563EB") +
    ggplot2::facet_wrap(~ ClusterLabel, scales = "free_y") +
    ggplot2::labs(
      title = paste0("Mfuzz Clusters - ", project_name),
      x = NULL, y = "Z-score"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 10),
      strip.text  = ggplot2::element_text(face = "bold", size = 11),
      plot.title  = ggplot2::element_text(face = "bold", size = 14)
    )

  save_plot_and_data(
    plot_obj     = p,
    data_df      = center_df,
    project_name = project_name,
    suffix       = "Mfuzz_Centers",
    width        = 12,
    height       = 3 * ceiling(n_cl / 2)
  )

  # 保存 cluster 成员列表
  membership <- mfuzz_result$membership
  cluster_assign <- cl$cluster
  member_df <- data.frame(
    Gene    = rownames(membership),
    Cluster = cluster_assign,
    membership,
    stringsAsFactors = FALSE
  )
  utils::write.csv(member_df,
    file.path("output", paste0(project_name, "_Mfuzz_Membership.csv")),
    row.names = FALSE)

  message(sprintf("  \u2705 Mfuzz \u53ef\u89c6\u5316\u5b8c\u6210: %d \u7c7b, %d \u86cb\u767d",
                  n_cl, nrow(membership)))
  p
}
