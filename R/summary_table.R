# ==============================================================================
# msProteomiX — Summary Table Generator
# ==============================================================================

#' Generate per-sample summary table
#'
#' Computes key QC metrics for each sample and returns a data.frame.
#' Metrics include: PSM count, Peptide count, Protein Groups,
#' 0-Miss Cleavage rate, Cysteine peptide percentage, and
#' Alkylation efficiency.
#'
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame
#' @param output_dir Output directory (default "output")
#' @param project_name Project name for output file naming
#' @return data.frame with summary metrics per sample (invisible)
#' @export
generate_summary_table <- function(ms_data, group_info,
                                   output_dir = "output",
                                   project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  samples <- ms_data$sample_names
  prot_mat <- ms_data$proteins
  psm_df  <- ms_data$psms

  # --- count_ids gives PSM, Peptide, Protein Groups ---
  id_df <- count_ids(ms_data, group_info)

  # --- Per-sample PSM-level metrics (0-miss, Cys%, Alk%) ---
  has_psm <- !is.null(psm_df) && nrow(psm_df) > 0
  has_sf  <- has_psm && "Spectrum File" %in% colnames(psm_df)

  # Prepare per-sample PSM metric containers
  zero_miss_vals <- rep(NA_real_, length(samples))
  cys_pct_vals   <- rep(NA_real_, length(samples))
  alk_eff_vals   <- rep(NA_real_, length(samples))

  if (has_psm && has_sf) {
    sf_col <- as.character(psm_df[["Spectrum File"]])
    for (i in seq_along(samples)) {
      matched <- grepl(samples[i], sf_col, fixed = TRUE)
      sub <- psm_df[matched, ]
      if (nrow(sub) == 0) next
      zero_miss_vals[i] <- calc_zero_miss(sub)
      cys_pct_vals[i]   <- calc_cys_percent(sub)
      alk_eff_vals[i]    <- calc_alk_efficiency(sub)
    }
  } else if (has_psm) {
    # No per-sample split possible, compute overall
    overall_zm  <- calc_zero_miss(psm_df)
    overall_cys <- calc_cys_percent(psm_df)
    overall_alk <- calc_alk_efficiency(psm_df)
    zero_miss_vals <- rep(overall_zm, length(samples))
    cys_pct_vals   <- rep(overall_cys, length(samples))
    alk_eff_vals   <- rep(overall_alk, length(samples))
  }

  # --- Assemble summary ---
  result <- data.frame(
    Sample         = id_df$sample_name,
    Group          = id_df$group,
    PSM            = id_df$psm,
    Peptide        = id_df$peptide,
    Protein_Groups = id_df$protein_groups,
    Zero_Miss_Pct  = zero_miss_vals[match(id_df$sample_name, samples)],
    Cys_Peptide_Pct = cys_pct_vals[match(id_df$sample_name, samples)],
    Alkylation_Pct = alk_eff_vals[match(id_df$sample_name, samples)],
    stringsAsFactors = FALSE
  )

  # --- Save CSV ---
  out_file <- file.path(output_dir,
                        paste0(project_name, "_Summary_Table.csv"))
  utils::write.csv(result, out_file, row.names = FALSE)
  message(sprintf("  Generated: %s (%d samples x %d metrics)",
                  basename(out_file), nrow(result), ncol(result) - 2))

  invisible(result)
}
