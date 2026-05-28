# ==============================================================================
# msProteomiX — 定性统计 QC 模块
# ==============================================================================

#' 计算 PSM/Peptide/ProteinGroup 鉴定数量
#'
#' 从 MsDataSet 中统计每个样本的 PSM、肽段和蛋白组鉴定数量。
#' FragPipe: 利用 psm.tsv 的 Spectrum File 列计算 PSM,
#'           利用 combined_peptide.tsv 的 Spectral Count 列计算 Peptide。
#' 其他引擎: 从蛋白矩阵非零值估算 Protein Groups。
#'
#' @param ms_data MsDataSet 对象
#' @param group_info 分组信息 data.frame
#' @return data.frame，包含 sample_name, group, protein_groups, peptide, psm 列
#' @export
count_ids <- function(ms_data, group_info) {
  stopifnot(inherits(ms_data, "MsDataSet"))

  prot_mat <- ms_data$proteins
  samples <- ms_data$sample_names

  # --- 统计 Protein Groups (所有引擎通用: 非零/非NA 计数) ---
  prot_counts <- sapply(samples, function(s) {
    if (s %in% colnames(prot_mat)) {
      sum(prot_mat[[s]] > 0, na.rm = TRUE)
    } else 0
  })


  # --- 统计 PSM (FragPipe: 从 psm.tsv 按 Spectrum File 计算) ---
  psm_counts <- rep(NA_integer_, length(samples))

  if (!is.null(ms_data$psms) && nrow(ms_data$psms) > 0 &&
      "Spectrum File" %in% colnames(ms_data$psms)) {
    # psm.tsv 的 Spectrum File 包含样本文件名 (如 "20241021_293T_96_1_A1.mzML")
    psm_df <- ms_data$psms
    sf_col <- psm_df[["Spectrum File"]]

    for (i in seq_along(samples)) {
      # 尝试匹配: 样本名出现在文件名中
      matched <- grepl(samples[i], sf_col, fixed = TRUE)
      psm_counts[i] <- sum(matched)
    }
  }

  # --- 统计 Peptide (FragPipe: 从 combined_peptide.tsv 的 Spectral Count 列) ---
  pep_counts <- rep(NA_integer_, length(samples))

  if (!is.null(ms_data$peptides) && nrow(ms_data$peptides) > 0) {
    pep_df <- ms_data$peptides

    # 方法1: 查找 per-sample Spectral Count 列
    sc_cols <- grep("Spectral Count$", colnames(pep_df), value = TRUE)
    sc_cols <- sc_cols[!grepl("(Combined|Total|Unique)", sc_cols)]

    if (length(sc_cols) > 0) {
      # 提取样本名 (去掉 " Spectral Count" 后缀)
      sc_samples <- trimws(sub("\\s*Spectral Count$", "", sc_cols))

      for (i in seq_along(samples)) {
        idx <- match(samples[i], sc_samples)
        if (!is.na(idx)) {
          vals <- suppressWarnings(as.numeric(pep_df[[sc_cols[idx]]]))
          pep_counts[i] <- sum(vals > 0, na.rm = TRUE)
        }
      }
    }
  }

  # --- 组装结果 ---
  row_list <- lapply(seq_along(samples), function(i) {
    samp <- samples[i]
    grp <- group_info$user_group[match(samp, group_info$sample_name)]
    if (is.na(grp)) return(NULL)

    data.frame(
      sample_name    = samp,
      group          = grp,
      protein_groups = prot_counts[i],
      peptide        = ifelse(is.na(pep_counts[i]), 0L, pep_counts[i]),
      psm            = ifelse(is.na(psm_counts[i]), 0L, psm_counts[i]),
      stringsAsFactors = FALSE
    )
  })

  do.call(rbind, Filter(Negate(is.null), row_list))
}


#' 计算 0-Miss Cleavage 比例
#'
#' 从 PSM 数据中计算去冗余后 0 漏切比例。
#'
#' @param psm_df PSM 水平 data.frame (需含 Peptide 和 Number of Missed Cleavages 列)
#' @return 百分比数值
#' @export
calc_zero_miss <- function(psm_df) {
  if (nrow(psm_df) == 0) return(NA_real_)

  # 去冗余
  if ("Peptide" %in% colnames(psm_df)) {
    psm_df <- psm_df[!duplicated(psm_df$Peptide), ]
  }

  mc_col <- grep("Missed.Cleavage|Number of Missed Cleavages",
                  colnames(psm_df), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(mc_col)) return(NA_real_)

  round((sum(psm_df[[mc_col]] == 0, na.rm = TRUE) / nrow(psm_df)) * 100, digits = 2)
}


#' 计算含 Cys 的 PSM 百分比
#'
#' @param psm_df PSM 水平 data.frame (需含 Peptide 列)
#' @return 百分比数值
#' @export
calc_cys_percent <- function(psm_df) {
  if (nrow(psm_df) == 0) return(NA_real_)
  if (!"Peptide" %in% colnames(psm_df)) return(NA_real_)

  cys_count <- sum(grepl("C", psm_df$Peptide))
  round((cys_count / nrow(psm_df)) * 100, digits = 2)
}


#' 计算烷基化效率
#'
#' @param psm_df PSM 水平 data.frame (需含 Peptide 和 Assigned Modifications 列)
#' @return 百分比数值
#' @export
calc_alk_efficiency <- function(psm_df) {
  if (nrow(psm_df) == 0) return(NA_real_)

  pep_col <- "Peptide"
  mod_col <- grep("Assigned.Modification", colnames(psm_df),
                   value = TRUE, ignore.case = TRUE)[1]

  if (is.na(mod_col) || !pep_col %in% colnames(psm_df)) return(NA_real_)

  cys_psm <- psm_df[grepl("C", psm_df[[pep_col]]), ]
  if (nrow(cys_psm) == 0) return(NA_real_)

  alk_count <- sum(grepl("C", cys_psm[[mod_col]]), na.rm = TRUE)
  round((alk_count / nrow(cys_psm)) * 100, digits = 2)
}
