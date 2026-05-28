# ==============================================================================
# msProteomiX — 定性统计 QC 模块
# ==============================================================================

#' 计算 PSM/Peptide/ProteinGroup 鉴定数量
#'
#' 从 MsDataSet 中统计每个样本的 PSM、肽段和蛋白组鉴定数量。
#'
#' @param ms_data MsDataSet 对象
#' @param group_info 分组信息 data.frame
#' @return data.frame，包含 Sample, Group, protein_groups, peptide, psm 等列
#' @export
count_ids <- function(ms_data, group_info) {
  stopifnot(inherits(ms_data, "MsDataSet"))

  prot_mat <- ms_data$proteins
  samples <- ms_data$sample_names

  row_list <- lapply(seq_along(samples), function(i) {
    samp <- samples[i]
    grp <- group_info$user_group[match(samp, group_info$sample_name)]
    if (is.na(grp)) return(NULL)

    prot_count <- sum(prot_mat[, i] > 0, na.rm = TRUE)

    data.frame(
      sample_name    = samp,
      group          = grp,
      protein_groups = prot_count,
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
