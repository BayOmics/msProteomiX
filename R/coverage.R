# ==============================================================================
# msProteomiX — 序列覆盖度分析
# ==============================================================================

#' 计算蛋白序列覆盖度
#'
#' 基于 Peptide 表和 FASTA 数据库计算每个蛋白的序列覆盖度。
#' 优化版: 预建肽段→蛋白映射, 用 base R gregexpr 替代 stringr。
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

  # Choose data source: peptides first, fallback to psms
  if (nrow(ms_data$peptides) > 0) {
    pep_df <- ms_data$peptides
    data_source <- "peptides"
  } else if (nrow(ms_data$psms) > 0) {
    pep_df <- ms_data$psms
    data_source <- "psms"
  } else {
    stop("No peptide or PSM data in MsDataSet.")
  }
  message(sprintf("  Using %s data (%d rows) for coverage calculation.", data_source, nrow(pep_df)))

  fasta_data <- Biostrings::readAAStringSet(fasta_path)

  # 处理 FASTA ID
  raw_names <- names(fasta_data)
  fasta_ids <- stringr::str_extract(raw_names, "(?<=\\|)[^|]+(?=\\|)")
  idx_na <- is.na(fasta_ids)
  if (any(idx_na)) fasta_ids[idx_na] <- sub("\\s.*", "", raw_names[idx_na])
  names(fasta_data) <- fasta_ids

  # 预提取所有蛋白序列为字符向量 (一次性, 避免重复转换)
  fasta_strings <- as.character(fasta_data)
  fasta_widths <- Biostrings::width(fasta_data)

  # 处理 Peptide ID
  prot_col <- grep("Protein ID|^Protein$", colnames(pep_df), value = TRUE)[1]
  if (is.na(prot_col)) stop("Cannot find Protein ID column in peptide data.")

  pep_df$CleanID <- sapply(as.character(pep_df[[prot_col]]), function(x) {
    if (grepl("\\|", x)) {
      parts <- unlist(strsplit(x, "\\|"))
      if (length(parts) >= 2) return(parts[2])
    }
    x
  }, USE.NAMES = FALSE)

  # 序列列
  seq_col <- grep("^Peptide Sequence$|^Peptide$|^Sequence$", colnames(pep_df), value = TRUE)[1]
  pep_raw_cols <- colnames(pep_df)

  common_ids <- intersect(pep_df$CleanID, fasta_ids)
  if (length(common_ids) == 0) stop("No matching protein IDs between peptides and FASTA.")
  message(sprintf("  Matched proteins: %d", length(common_ids)))

  # 预过滤: 只保留匹配上 FASTA 的肽段
  pep_df_valid <- pep_df[pep_df$CleanID %in% common_ids, ]
  has_pos <- all(c("Start", "End") %in% colnames(pep_df_valid))
  keep_cols <- c("CleanID", seq_col)
  if (has_pos) keep_cols <- c(keep_cols, "Start", "End")

  # Detect if data is long-format (has Spectrum File column, e.g. Spectronaut)
  is_long_format <- "Spectrum File" %in% colnames(pep_df_valid)

  # 核心覆盖度计算
  .calc_one_prot <- function(pid, pep_by_prot_grp) {
    fasta_idx <- match(pid, fasta_ids)
    if (is.na(fasta_idx)) return(NA_real_)
    prot_len <- fasta_widths[fasta_idx]
    if (prot_len == 0) return(NA_real_)

    prot_str <- fasta_strings[fasta_idx]
    pep_rows <- pep_by_prot_grp[[pid]]
    if (is.null(pep_rows) || nrow(pep_rows) == 0) return(0)

    cover_mask <- logical(prot_len)

    # 优先使用 Start/End 精确位置 (与参考脚本一致)
    if (has_pos) {
      starts <- pep_rows$Start
      ends <- pep_rows$End
      valid_idx <- !is.na(starts) & !is.na(ends) & starts >= 1 & ends <= prot_len
      if (any(valid_idx)) {
        for (k in which(valid_idx)) {
          cover_mask[starts[k]:ends[k]] <- TRUE
        }
      }
    } else {
      # Fallback: 字符串匹配
      pep_seqs <- unique(pep_rows[[seq_col]])
      for (ps in pep_seqs) {
        hits <- gregexpr(ps, prot_str, fixed = TRUE)[[1]]
        if (hits[1] > 0) {
          pep_len <- nchar(ps)
          for (h in hits) {
            cover_mask[h:(h + pep_len - 1L)] <- TRUE
          }
        }
      }
    }
    sum(cover_mask) / prot_len * 100
  }

  coverage_list <- list()
  unique_grps <- unique(group_info$user_group)
  n_grps <- length(unique_grps)

  for (g_idx in seq_along(unique_grps)) {
    grp <- unique_grps[g_idx]
    grp_samples <- group_info$sample_name[group_info$user_group == grp]

    if (is_long_format) {
      # Long-format (Spectronaut psms): filter by Spectrum File column
      # Spectrum File values may partially match sample names
      spec_files <- unique(pep_df_valid$`Spectrum File`)
      matched_files <- c()
      for (s in grp_samples) {
        matched_files <- c(matched_files, spec_files[grepl(s, spec_files, fixed = TRUE)])
      }
      if (length(matched_files) == 0) next
      grp_pep <- pep_df_valid[pep_df_valid$`Spectrum File` %in% matched_files, keep_cols, drop = FALSE]
    } else {
      # Wide-format (FragPipe peptides): filter by intensity columns
      grp_cols <- c()
      for (s in grp_samples) {
        s_safe <- escape_regex(s)
        pat <- paste0("^", s_safe, ".*(Intensity|Spectral Count)$")
        hits <- grep(pat, pep_raw_cols, ignore.case = TRUE, value = TRUE)
        hits <- hits[!grepl("(Unique|Total)", hits, ignore.case = TRUE)]
        grp_cols <- c(grp_cols, hits)
      }
      if (length(grp_cols) == 0) next

      # 关键: 只用当前组检出的肽段行 (与参考脚本一致)
      is_detected <- rowSums(pep_df_valid[, grp_cols, drop = FALSE] > 0, na.rm = TRUE) > 0
      grp_pep <- pep_df_valid[is_detected, keep_cols, drop = FALSE]
    }

    detected_prots <- unique(grp_pep$CleanID)
    if (length(detected_prots) == 0) next

    # 每组重建蛋白→肽段映射
    pep_by_prot_grp <- split(grp_pep, grp_pep$CleanID)

    message(sprintf("  [%d/%d] %s: calculating %d proteins...", g_idx, n_grps, grp, length(detected_prots)))

    covs <- vapply(detected_prots, function(pid) .calc_one_prot(pid, pep_by_prot_grp), numeric(1))
    valid <- !is.na(covs)

    if (any(valid)) {
      coverage_list[[grp]] <- data.frame(
        ProteinID = detected_prots[valid],
        Coverage  = covs[valid],
        Group     = grp,
        stringsAsFactors = FALSE
      )
    }
  }

  if (length(coverage_list) == 0) return(data.frame())
  result <- do.call(rbind, coverage_list)
  rownames(result) <- NULL
  result
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
