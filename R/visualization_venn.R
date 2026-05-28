# ==============================================================================
# msProteomiX — Venn / UpSet 可视化
# ==============================================================================

#' 绘制 Venn 维恩图
#'
#' @param ms_data MsDataSet 对象
#' @param group_info 分组信息 data.frame
#' @param selected_groups 选择绘制的组名 (NULL 表示全部，最多 4 组)
#' @param output_dir 输出目录
#' @param project_name 项目名
#' @return ggplot 对象
#' @export
plot_venn <- function(ms_data, group_info,
                      selected_groups = NULL,
                      output_dir = "output",
                      project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  if (!requireNamespace("ggvenn", quietly = TRUE)) {
    stop("Please install ggvenn: install.packages('ggvenn')")
  }
  ensure_output_dir(output_dir)

  # 构建蛋白 ID 列表
  venn_list <- .build_venn_list(ms_data, group_info)

  if (is.null(selected_groups)) selected_groups <- names(venn_list)
  if (length(selected_groups) < 2) stop("At least 2 groups needed for Venn plot.")
  if (length(selected_groups) > 4) stop("ggvenn supports at most 4 groups.")

  plot_list <- venn_list[selected_groups]

  p <- ggvenn::ggvenn(
    plot_list,
    fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF")[1:length(plot_list)],
    stroke_size = 0.5,
    set_name_size = 4,
    text_size = 4.5,
    show_percentage = TRUE
  ) +
    ggplot2::labs(title = "Protein Overlap Analysis") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14))

  grp_str <- paste(selected_groups, collapse = "_")
  suffix <- paste0("Venn_", grp_str)
  save_plot_and_data(p, .venn_to_df(plot_list), project_name, suffix,
                     output_dir = output_dir)

  p
}


#' 绘制 UpSet 图
#'
#' @param ms_data MsDataSet 对象
#' @param group_info 分组信息 data.frame
#' @param output_dir 输出目录
#' @param project_name 项目名
#' @export
plot_upset <- function(ms_data, group_info,
                       output_dir = "output",
                       project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  if (!requireNamespace("UpSetR", quietly = TRUE)) {
    stop("Please install UpSetR: install.packages('UpSetR')")
  }
  ensure_output_dir(output_dir)

  venn_list <- .build_venn_list(ms_data, group_info)

  if (length(venn_list) < 2) {
    message("  At least 2 groups needed for UpSet plot. Skipping.")
    return(invisible(NULL))
  }

  pdf_name <- file.path(output_dir, paste0(project_name, "_UpSet.pdf"))
  grDevices::pdf(pdf_name, width = 8, height = 6, onefile = FALSE)
  tryCatch({
    p <- UpSetR::upset(UpSetR::fromList(venn_list), order.by = "freq",
                        mainbar.y.label = "Protein Intersections",
                        sets.x.label = "Proteins Per Group")
    p
    message("  Generated: UpSet plot")
  }, error = function(e) message(paste("  UpSet plot failed:", e$message)))
  grDevices::dev.off()
}


# --- 内部辅助 ---

#' @keywords internal
.build_venn_list <- function(ms_data, group_info) {
  prot_mat <- ms_data$proteins
  prot_ids <- if ("Protein ID" %in% colnames(ms_data$protein_info)) {
    ms_data$protein_info$`Protein ID`
  } else if ("Label_Name" %in% colnames(ms_data$protein_info)) {
    ms_data$protein_info$Label_Name
  } else {
    paste0("Prot_", seq_len(nrow(prot_mat)))
  }

  unique_grps <- unique(group_info$user_group)
  venn_list <- list()

  for (grp in unique_grps) {
    samps <- group_info$sample_name[group_info$user_group == grp]
    col_idx <- match(samps, ms_data$sample_names)
    col_idx <- col_idx[!is.na(col_idx)]
    if (length(col_idx) == 0) next

    sub_mat <- prot_mat[, col_idx, drop = FALSE]
    means <- rowMeans(as.matrix(sub_mat), na.rm = TRUE)
    detected <- !is.na(means) & means > 0
    venn_list[[grp]] <- as.character(prot_ids[detected])
  }

  venn_list
}

#' @keywords internal
.venn_to_df <- function(venn_list) {
  all_prot <- unique(unlist(venn_list))
  df <- data.frame(Protein = all_prot, stringsAsFactors = FALSE)
  for (grp in names(venn_list)) {
    df[[grp]] <- df$Protein %in% venn_list[[grp]]
  }
  df$Partition <- apply(df[, -1, drop = FALSE], 1, function(row) {
    paste(names(row)[row], collapse = " & ")
  })
  df
}
