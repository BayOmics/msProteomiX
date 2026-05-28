# ==============================================================================
# msProteomiX — 序列覆盖度分析
# ==============================================================================

#' 计算蛋白序列覆盖度
#'
#' 基于 Peptide 表和 FASTA 数据库计算每个蛋白的序列覆盖度。
#'
#' @param ms_data MsDataSet 对象 (需含 peptides 数据)
#' @param fasta_path FASTA 数据库文件路径
#' @param group_info 分组信息 data.frame
#' @return data.frame，包含 ProteinID, Coverage, Group
#' @export
calc_coverage <- function(ms_data, fasta_path, group_info) {
  stopifnot(inherits(ms_data, "MsDataSet"))
  if (!requireNamespace("Biostrings", quietly = TRUE)) {
    stop("Please install Biostrings: BiocManager::install('Biostrings')")
  }
  if (nrow(ms_data$peptides) == 0) stop("No peptide data in MsDataSet.")

  pep_df <- ms_data$peptides
  fasta_data <- Biostrings::readAAStringSet(fasta_path)

  # 处理 FASTA ID
  raw_names <- names(fasta_data)
  fasta_ids <- stringr::str_extract(raw_names, "(?<=\\|)[^|]+(?=\\|)")
  idx_na <- is.na(fasta_ids)
  if (any(idx_na)) fasta_ids[idx_na] <- sub("\\s.*", "", raw_names[idx_na])
  names(fasta_data) <- fasta_ids

  # 处理 Peptide ID
  prot_col <- grep("Protein ID|Protein", colnames(pep_df), value = TRUE)[1]
  if (is.na(prot_col)) stop("Cannot find Protein ID column in peptide data.")

  pep_df$CleanID <- sapply(as.character(pep_df[[prot_col]]), function(x) {
    if (grepl("\\|", x)) {
      parts <- unlist(strsplit(x, "\\|"))
      if (length(parts) >= 2) return(parts[2])
    }
    x
  })

  common_ids <- intersect(pep_df$CleanID, fasta_ids)
  if (length(common_ids) == 0) stop("No matching protein IDs between peptides and FASTA.")
  message(sprintf("  Matched proteins: %d", length(common_ids)))

  # 序列列
  seq_col <- grep("Peptide Sequence|Sequence", colnames(pep_df), value = TRUE)[1]
  pep_raw_cols <- colnames(pep_df)

  coverage_list <- list()
  unique_grps <- unique(group_info$user_group)

  for (grp in unique_grps) {
    grp_samples <- group_info$sample_name[group_info$user_group == grp]
    grp_cols <- c()
    for (s in grp_samples) {
      s_safe <- escape_regex(s)
      pat <- paste0("^", s_safe, ".*(Intensity|Spectral Count)$")
      hits <- grep(pat, pep_raw_cols, ignore.case = TRUE, value = TRUE)
      hits <- hits[!grepl("(Unique|Total)", hits, ignore.case = TRUE)]
      grp_cols <- c(grp_cols, hits)
    }
    if (length(grp_cols) == 0) next

    is_detected <- rowSums(pep_df[, grp_cols, drop = FALSE] > 0, na.rm = TRUE) > 0
    sub_pep <- pep_df[is_detected, ]
    valid_peps <- sub_pep[sub_pep$CleanID %in% common_ids, ]
    if (nrow(valid_peps) == 0) next

    unique_prots <- unique(valid_peps$CleanID)
    grp_res <- list()

    for (pid in unique_prots) {
      fasta_idx <- match(pid, names(fasta_data))
      if (is.na(fasta_idx)) next

      curr_peps <- valid_peps[valid_peps$CleanID == pid, ]
      prot_len <- Biostrings::width(fasta_data[fasta_idx])
      if (prot_len == 0) next

      cover_mask <- logical(prot_len)
      pep_seqs <- unique(curr_peps[[seq_col]])
      prot_str <- as.character(fasta_data[[fasta_idx]])

      for (ps in pep_seqs) {
        matches <- stringr::str_locate_all(prot_str, stringr::fixed(ps))[[1]]
        if (nrow(matches) > 0) {
          for (m in seq_len(nrow(matches))) {
            cover_mask[matches[m, 1]:matches[m, 2]] <- TRUE
          }
        }
      }

      grp_res[[pid]] <- data.frame(
        ProteinID = pid,
        Coverage  = sum(cover_mask) / prot_len * 100,
        Group     = grp,
        stringsAsFactors = FALSE
      )
    }

    if (length(grp_res) > 0) coverage_list[[grp]] <- do.call(rbind, grp_res)
  }

  if (length(coverage_list) == 0) return(data.frame())
  do.call(rbind, coverage_list)
}


#' 绘制序列覆盖度箱线图
#'
#' @param coverage_df calc_coverage() 返回的 data.frame
#' @param output_dir 输出目录
#' @param project_name 项目名
#' @return ggplot 对象
#' @export
plot_coverage <- function(coverage_df,
                           output_dir = "output",
                           project_name = "Project") {
  ensure_output_dir(output_dir)

  if (is.null(coverage_df) || nrow(coverage_df) == 0) {
    message("  No coverage data to plot.")
    return(invisible(NULL))
  }

  n_groups <- length(unique(coverage_df$Group))

  # 中位数标签
  medians <- coverage_df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(MedianCov = stats::median(Coverage, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(Label = paste0(round(MedianCov, 1), "%"))

  p <- ggplot2::ggplot(coverage_df, ggplot2::aes(x = Group, y = Coverage, fill = Group)) +
    ggplot2::geom_boxplot(alpha = 0.6, outlier.shape = NA) +
    ggplot2::geom_text(data = medians,
                       ggplot2::aes(x = Group, y = MedianCov, label = Label),
                       vjust = -0.5, fontface = "bold", size = 4) +
    ggplot2::scale_fill_manual(values = mspx_colors(n_groups)) +
    ggplot2::labs(title = "Protein Sequence Coverage", y = "Coverage (%)") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  save_plot_and_data(p, coverage_df, project_name, "Sequence_Coverage",
                     output_dir = output_dir)
  p
}
