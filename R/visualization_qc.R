# ==============================================================================
# msProteomiX — QC Method Evaluation Visualization
# ==============================================================================

# ------------------------------------------------------------------------------
# Internal helper: match PSM rows to samples via Spectrum File column
# Returns a data.frame with original PSM columns + 'qc_sample' + 'qc_group'
# ------------------------------------------------------------------------------
.psm_with_groups <- function(psm_df, group_info, sample_names) {
  if (is.null(group_info) || !"Spectrum File" %in% colnames(psm_df)) {
    psm_df$qc_sample <- "All"
    psm_df$qc_group  <- "All"
    return(psm_df)
  }

  psm_df$qc_sample <- NA_character_
  psm_df$qc_group  <- NA_character_
  sf_col <- as.character(psm_df[["Spectrum File"]])

  for (i in seq_len(nrow(group_info))) {
    sname <- group_info$sample_name[i]
    matched <- grepl(sname, sf_col, fixed = TRUE)
    psm_df$qc_sample[matched] <- sname
    psm_df$qc_group[matched]  <- group_info$user_group[i]
  }

  psm_df <- psm_df[!is.na(psm_df$qc_group), ]
  psm_df
}


# ==============================================================================
# 1. Peptide length distribution (by group)
# ==============================================================================

#' Plot peptide length distribution
#'
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame (NULL = overall)
#' @param output_dir Output directory
#' @param project_name Project name
#' @return ggplot object or NULL
#' @export
plot_qc_peptide_length <- function(ms_data, group_info = NULL,
                                   output_dir = "output",
                                   project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  psm_df <- ms_data$psms
  if (is.null(psm_df) || nrow(psm_df) == 0) {
    message(">>> No PSM data available. Skipping peptide length plot.")
    return(invisible(NULL))
  }

  pep_col <- .find_column(psm_df, c("Peptide", "Modified Sequence",
                                     "Stripped.Sequence", "Sequence"))
  if (is.null(pep_col)) {
    message(">>> Cannot find peptide sequence column. Skipping.")
    return(invisible(NULL))
  }

  psm_df <- .psm_with_groups(psm_df, group_info, ms_data$sample_names)

  seqs <- gsub("[^A-Za-z]", "", as.character(psm_df[[pep_col]]))
  pep_lengths <- nchar(seqs)
  plot_df <- data.frame(peptide_length = pep_lengths,
                        Group = psm_df$qc_group,
                        stringsAsFactors = FALSE)
  plot_df <- plot_df[plot_df$peptide_length > 0 & !is.na(plot_df$peptide_length), ]

  if (nrow(plot_df) == 0) return(invisible(NULL))

  n_groups <- length(unique(plot_df$Group))
  colors <- mspx_colors(n_groups)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = peptide_length, color = Group)) +
    ggplot2::geom_density(linewidth = 1) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_x_continuous(limits = c(5, 55), breaks = seq(5, 55, 5)) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Peptide Length Distribution",
                  x = "Peptide Length (amino acids)", y = "Density") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  save_plot_and_data(p, plot_df, project_name, "QC_Peptide_Length",
                     output_dir = output_dir, width = 8, height = 5)
  p
}


# ==============================================================================
# 2. Precursor charge distribution (by group)
# ==============================================================================

#' Plot precursor charge distribution
#'
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame (NULL = overall)
#' @param output_dir Output directory
#' @param project_name Project name
#' @return ggplot object or NULL
#' @export
plot_qc_charge_distribution <- function(ms_data, group_info = NULL,
                                        output_dir = "output",
                                        project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  psm_df <- ms_data$psms
  if (is.null(psm_df) || nrow(psm_df) == 0) {
    message(">>> No PSM data available. Skipping charge distribution plot.")
    return(invisible(NULL))
  }

  charge_col <- .find_column(psm_df, c("Charge", "Precursor Charge",
                                        "charge", "PrecursorCharge"))
  if (is.null(charge_col)) {
    message(">>> Cannot find charge column. Skipping.")
    return(invisible(NULL))
  }

  psm_df <- .psm_with_groups(psm_df, group_info, ms_data$sample_names)

  charges <- suppressWarnings(as.integer(psm_df[[charge_col]]))
  plot_df <- data.frame(charge = factor(charges), Group = psm_df$qc_group,
                        stringsAsFactors = FALSE)
  plot_df <- plot_df[!is.na(plot_df$charge), ]

  if (nrow(plot_df) == 0) return(invisible(NULL))

  n_groups <- length(unique(plot_df$Group))
  colors <- mspx_colors(n_groups)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = charge, fill = Group)) +
    ggplot2::geom_bar(position = "dodge", color = "white") +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Precursor Charge Distribution",
                  x = "Charge State", y = "Number of PSMs") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  save_plot_and_data(p, plot_df, project_name, "QC_Charge_Distribution",
                     output_dir = output_dir, width = 7, height = 5)
  p
}


# ==============================================================================
# 3. Missed cleavage distribution (stacked by group)
# ==============================================================================

#' Plot missed cleavage distribution
#'
#' Stacked bar chart showing missed cleavage proportions per group.
#'
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame (NULL = overall)
#' @param output_dir Output directory
#' @param project_name Project name
#' @return ggplot object or NULL
#' @export
plot_qc_missed_cleavage <- function(ms_data, group_info = NULL,
                                    output_dir = "output",
                                    project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  psm_df <- ms_data$psms
  if (is.null(psm_df) || nrow(psm_df) == 0) {
    message(">>> No PSM data available. Skipping missed cleavage plot.")
    return(invisible(NULL))
  }

  psm_df <- .psm_with_groups(psm_df, group_info, ms_data$sample_names)

  mc_col <- .find_column(psm_df, c("Number of Missed Cleavages",
                                    "Missed.Cleavage", "Missed Cleavages"))
  if (!is.null(mc_col)) {
    mc_values <- suppressWarnings(as.integer(psm_df[[mc_col]]))
  } else {
    pep_col <- .find_column(psm_df, c("Peptide", "Stripped.Sequence", "Sequence"))
    if (is.null(pep_col)) {
      message(">>> Cannot determine missed cleavages. Skipping.")
      return(invisible(NULL))
    }
    mc_values <- calc_missed_cleavage(as.character(psm_df[[pep_col]]))
  }

  mc_str <- as.character(mc_values)
  mc_str[mc_values >= 2] <- "2+"
  plot_df <- data.frame(Group = psm_df$qc_group, MC = mc_str,
                        stringsAsFactors = FALSE)
  plot_df <- plot_df[!is.na(plot_df$MC), ]

  if (nrow(plot_df) == 0) return(invisible(NULL))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Group, fill = MC)) +
    ggplot2::geom_bar(position = "fill") +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::scale_fill_brewer(palette = "Blues") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Missed Cleavage Rate",
                  x = NULL, y = "Proportion", fill = "Missed\nCleavages") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  # Summary table
  summ_df <- as.data.frame(table(plot_df$Group, plot_df$MC))
  colnames(summ_df) <- c("Group", "MC", "Count")

  save_plot_and_data(p, summ_df, project_name, "QC_Missed_Cleavage",
                     output_dir = output_dir, width = 7, height = 5)
  p
}


# ==============================================================================
# 4. Modification type distribution
# ==============================================================================

#' Plot modification type distribution
#'
#' @param ms_data MsDataSet object
#' @param top_n Number of top modifications to show
#' @param output_dir Output directory
#' @param project_name Project name
#' @return ggplot object or NULL
#' @export
plot_qc_modification <- function(ms_data, top_n = 10,
                                 output_dir = "output",
                                 project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  psm_df <- ms_data$psms
  if (is.null(psm_df) || nrow(psm_df) == 0) {
    message(">>> No PSM data available. Skipping modification plot.")
    return(invisible(NULL))
  }

  mod_col <- .find_column(psm_df, c("Assigned Modifications", "Modifications",
                                     "Modified Sequence", "Modified Peptide",
                                     "assigned_modifications"))
  if (is.null(mod_col)) {
    message(">>> Cannot find modification column. Skipping.")
    return(invisible(NULL))
  }

  mod_strings <- as.character(psm_df[[mod_col]])
  mod_strings <- mod_strings[!is.na(mod_strings) & mod_strings != ""]
  if (length(mod_strings) == 0) return(invisible(NULL))

  # 检测格式: Spectronaut = [Mod (X)] inline; FragPipe = "N-term(x)" comma/semicolon-sep
  if (any(grepl("\\[", mod_strings))) {
    # Spectronaut format: extract [Modification (X)] patterns
    all_mods <- unlist(regmatches(mod_strings, gregexpr("\\[([^]]+)\\]", mod_strings)))
    all_mods <- gsub("^\\[|\\]$", "", all_mods)  # remove brackets
  } else {
    # FragPipe / generic format: split by , or ;
    all_mods <- unlist(strsplit(mod_strings, "[,;]"))
    all_mods <- trimws(all_mods)
    all_mods <- gsub("^[0-9]+", "", all_mods)
    all_mods <- trimws(all_mods)
  }
  all_mods <- all_mods[all_mods != ""]
  if (length(all_mods) == 0) return(invisible(NULL))

  mod_table <- sort(table(all_mods), decreasing = TRUE)
  mod_df <- data.frame(mod_type = names(mod_table),
                       mod_count = as.integer(mod_table),
                       stringsAsFactors = FALSE)
  mod_df <- utils::head(mod_df, top_n)
  mod_df$mod_type <- factor(mod_df$mod_type, levels = rev(mod_df$mod_type))

  p <- ggplot2::ggplot(mod_df, ggplot2::aes(x = mod_type, y = mod_count)) +
    ggplot2::geom_bar(stat = "identity", fill = "#8491B4", color = "white") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw() +
    ggplot2::labs(title = paste("Top", min(top_n, nrow(mod_df)),
                                "Modification Types"),
                  x = NULL, y = "Count") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  save_plot_and_data(p, mod_df, project_name, "QC_Modification_Types",
                     output_dir = output_dir, width = 9, height = 6)
  p
}


# ==============================================================================
# 5. GRAVY hydrophobicity distribution (by group)
# ==============================================================================

#' Plot GRAVY hydrophobicity distribution
#'
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame (NULL = overall)
#' @param output_dir Output directory
#' @param project_name Project name
#' @return ggplot object or NULL
#' @export
plot_qc_gravy <- function(ms_data, group_info = NULL,
                          output_dir = "output",
                          project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  psm_df <- ms_data$psms
  if (is.null(psm_df) || nrow(psm_df) == 0) {
    message(">>> No PSM data available. Skipping GRAVY plot.")
    return(invisible(NULL))
  }

  pep_col <- .find_column(psm_df, c("Peptide", "Stripped.Sequence", "Sequence"))
  if (is.null(pep_col)) {
    message(">>> Cannot find peptide column. Skipping GRAVY plot.")
    return(invisible(NULL))
  }

  psm_df <- .psm_with_groups(psm_df, group_info, ms_data$sample_names)

  # Subsample if large
  if (nrow(psm_df) > 20000) psm_df <- psm_df[sample(nrow(psm_df), 20000), ]

  gravy <- calc_gravy(as.character(psm_df[[pep_col]]))
  plot_df <- data.frame(GRAVY = gravy, Group = psm_df$qc_group,
                        stringsAsFactors = FALSE)
  plot_df <- plot_df[!is.na(plot_df$GRAVY), ]
  if (nrow(plot_df) == 0) return(invisible(NULL))

  n_groups <- length(unique(plot_df$Group))
  colors <- mspx_colors(n_groups)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = GRAVY, color = Group)) +
    ggplot2::geom_density(linewidth = 1) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Hydrophobicity (GRAVY) Distribution",
                  x = "GRAVY Score", y = "Density") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  save_plot_and_data(p, plot_df, project_name, "QC_GRAVY",
                     output_dir = output_dir, width = 8, height = 5)
  p
}


# ==============================================================================
# 6. pI isoelectric point distribution (by group)
# ==============================================================================

#' Plot peptide pI distribution
#'
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame (NULL = overall)
#' @param output_dir Output directory
#' @param project_name Project name
#' @return ggplot object or NULL
#' @export
plot_qc_pi <- function(ms_data, group_info = NULL,
                       output_dir = "output",
                       project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  psm_df <- ms_data$psms
  if (is.null(psm_df) || nrow(psm_df) == 0) {
    message(">>> No PSM data available. Skipping pI plot.")
    return(invisible(NULL))
  }

  pep_col <- .find_column(psm_df, c("Peptide", "Stripped.Sequence", "Sequence"))
  if (is.null(pep_col)) {
    message(">>> Cannot find peptide column. Skipping pI plot.")
    return(invisible(NULL))
  }

  psm_df <- .psm_with_groups(psm_df, group_info, ms_data$sample_names)

  if (nrow(psm_df) > 20000) psm_df <- psm_df[sample(nrow(psm_df), 20000), ]

  pI_vals <- calc_pI_seq(as.character(psm_df[[pep_col]]))
  plot_df <- data.frame(pI = pI_vals, Group = psm_df$qc_group,
                        stringsAsFactors = FALSE)
  plot_df <- plot_df[!is.na(plot_df$pI), ]
  if (nrow(plot_df) == 0) return(invisible(NULL))

  n_groups <- length(unique(plot_df$Group))
  colors <- mspx_colors(n_groups)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = pI, color = Group)) +
    ggplot2::geom_density(linewidth = 1) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Peptide pI Distribution",
                  x = "Isoelectric Point (pI)", y = "Density") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  save_plot_and_data(p, plot_df, project_name, "QC_pI",
                     output_dir = output_dir, width = 8, height = 5)
  p
}


# ==============================================================================
# 7. M/Z distribution (by group)
# ==============================================================================

#' Plot M/Z distribution
#'
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame (NULL = overall)
#' @param output_dir Output directory
#' @param project_name Project name
#' @return ggplot object or NULL
#' @export
plot_qc_mz <- function(ms_data, group_info = NULL,
                       output_dir = "output",
                       project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  psm_df <- ms_data$psms
  if (is.null(psm_df) || nrow(psm_df) == 0) {
    message(">>> No PSM data available. Skipping M/Z plot.")
    return(invisible(NULL))
  }

  mz_col <- .find_column(psm_df, c("Calibrated Observed M/Z", "Observed M/Z",
                                     "m/z", "MZ", "PrecursorMZ"))
  if (is.null(mz_col)) {
    message(">>> Cannot find M/Z column. Skipping.")
    return(invisible(NULL))
  }

  psm_df <- .psm_with_groups(psm_df, group_info, ms_data$sample_names)

  mz_vals <- suppressWarnings(as.numeric(psm_df[[mz_col]]))
  plot_df <- data.frame(MZ = mz_vals, Group = psm_df$qc_group,
                        stringsAsFactors = FALSE)
  plot_df <- plot_df[!is.na(plot_df$MZ) & plot_df$MZ > 0, ]
  if (nrow(plot_df) == 0) return(invisible(NULL))

  n_groups <- length(unique(plot_df$Group))
  colors <- mspx_colors(n_groups)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = MZ, color = Group)) +
    ggplot2::geom_density(linewidth = 1) +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "M/Z Distribution",
                  x = "M/Z", y = "Density") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  save_plot_and_data(p, plot_df, project_name, "QC_MZ",
                     output_dir = output_dir, width = 8, height = 5)
  p
}


# ==============================================================================
# 8. Cumulative intensity curve (dynamic range)
# ==============================================================================

#' Plot cumulative intensity curve
#'
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame
#' @param output_dir Output directory
#' @param project_name Project name
#' @return ggplot object or NULL
#' @export
plot_qc_cumulative_intensity <- function(ms_data, group_info = NULL,
                                         output_dir = "output",
                                         project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  prot_mat <- ms_data$proteins
  if (is.null(prot_mat) || ncol(prot_mat) == 0) {
    message(">>> No protein data. Skipping cumulative intensity plot.")
    return(invisible(NULL))
  }
  if (is.null(group_info)) {
    message(">>> group_info required for cumulative intensity. Skipping.")
    return(invisible(NULL))
  }

  unique_groups <- unique(group_info$user_group)
  cum_list <- list()

  for (grp in unique_groups) {
    grp_samples <- group_info$sample_name[group_info$user_group == grp]
    grp_samples <- intersect(grp_samples, colnames(prot_mat))
    if (length(grp_samples) == 0) next

    sub_mat <- as.matrix(prot_mat[, grp_samples, drop = FALSE])
    sub_mat[sub_mat == 0] <- NA
    means <- rowMeans(sub_mat, na.rm = TRUE)
    valid <- means[!is.na(means) & means > 0]
    if (length(valid) == 0) next

    sorted <- sort(valid, decreasing = TRUE)
    cum_prop <- cumsum(sorted) / sum(sorted)
    cum_list[[grp]] <- data.frame(Group = grp,
                                  Rank = seq_along(cum_prop),
                                  CumulativeProp = cum_prop,
                                  stringsAsFactors = FALSE)
  }

  if (length(cum_list) == 0) return(invisible(NULL))
  plot_df <- do.call(rbind, cum_list)

  n_groups <- length(unique(plot_df$Group))
  colors <- mspx_colors(n_groups)

  p <- ggplot2::ggplot(plot_df,
                        ggplot2::aes(x = Rank, y = CumulativeProp,
                                    color = Group)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey") +
    ggplot2::scale_color_manual(values = colors) +
    ggplot2::scale_y_continuous(labels = scales::percent) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Dynamic Range (Cumulative Intensity)",
                  x = "Protein Rank", y = "Cumulative %") +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))

  save_plot_and_data(p, plot_df, project_name, "QC_Cumulative_Intensity",
                     output_dir = output_dir, width = 8, height = 5)
  p
}


# ==============================================================================
# 9. Cysteine peptide ratio + Alkylation efficiency (per sample)
# ==============================================================================

#' Plot Cys-containing peptide ratio and alkylation efficiency
#'
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame
#' @param output_dir Output directory
#' @param project_name Project name
#' @return ggplot object or NULL
#' @export
plot_qc_cys_alkylation <- function(ms_data, group_info = NULL,
                                   output_dir = "output",
                                   project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  psm_df <- ms_data$psms
  if (is.null(psm_df) || nrow(psm_df) == 0) {
    message(">>> No PSM data available. Skipping Cys/Alkylation plot.")
    return(invisible(NULL))
  }

  psm_df <- .psm_with_groups(psm_df, group_info, ms_data$sample_names)
  unique_samples <- unique(psm_df$qc_sample)

  result_list <- list()
  for (samp in unique_samples) {
    sub <- psm_df[psm_df$qc_sample == samp, ]
    grp <- sub$qc_group[1]
    cys_pct <- calc_cys_percent(sub)
    alk_eff <- calc_alk_efficiency(sub)

    result_list[[samp]] <- data.frame(
      Sample = samp, Group = grp,
      Cys_Percent = ifelse(is.na(cys_pct), 0, cys_pct),
      Alk_Efficiency = ifelse(is.na(alk_eff), 0, alk_eff),
      stringsAsFactors = FALSE
    )
  }

  if (length(result_list) == 0) return(invisible(NULL))
  plot_df <- do.call(rbind, result_list)

  # Reshape to long format
  long_df <- rbind(
    data.frame(Sample = plot_df$Sample, Group = plot_df$Group,
               Metric = "Cys Peptide %", Value = plot_df$Cys_Percent,
               stringsAsFactors = FALSE),
    data.frame(Sample = plot_df$Sample, Group = plot_df$Group,
               Metric = "Alkylation Eff. %", Value = plot_df$Alk_Efficiency,
               stringsAsFactors = FALSE)
  )

  n_groups <- length(unique(long_df$Group))
  colors <- mspx_colors(n_groups)

  p <- ggplot2::ggplot(long_df,
                        ggplot2::aes(x = Sample, y = Value, fill = Group)) +
    ggplot2::geom_bar(stat = "identity", color = "white", width = 0.7) +
    ggplot2::facet_wrap(~ Metric, scales = "free_y", ncol = 1) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Cysteine Peptide Ratio & Alkylation Efficiency",
                  x = NULL, y = "%") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 7)
    )

  save_plot_and_data(p, plot_df, project_name, "QC_Cys_Alkylation",
                     output_dir = output_dir, width = 10, height = 7)
  p
}


# ==============================================================================
# 10. PSM over Retention Time (per sample)
# ==============================================================================

#' Plot PSM count over Retention Time
#'
#' Histogram of PSM identifications across retention time bins.
#' Requires "Retention" column in PSM data.
#'
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame
#' @param output_dir Output directory
#' @param project_name Project name
#' @param n_bins Number of RT bins (default 50)
#' @return ggplot object or NULL
#' @export
plot_qc_rt_distribution <- function(ms_data, group_info = NULL,
                                     output_dir = "output",
                                     project_name = "Project",
                                     n_bins = 50) {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  psm_df <- ms_data$psms
  if (is.null(psm_df) || nrow(psm_df) == 0) {
    message(">>> No PSM data. Skipping RT distribution.")
    return(invisible(NULL))
  }

  rt_col <- .find_column(psm_df, c("Retention", "RT", "EG.ApexRT"))
  if (is.null(rt_col)) {
    message(">>> No RT column found. Skipping RT distribution.")
    return(invisible(NULL))
  }

  psm_df <- .psm_with_groups(psm_df, group_info, ms_data$sample_names)
  psm_df$RT_val <- as.numeric(psm_df[[rt_col]])
  psm_df <- psm_df[!is.na(psm_df$RT_val), ]
  if (nrow(psm_df) == 0) return(invisible(NULL))

  n_groups <- length(unique(psm_df$qc_group))

  p <- ggplot2::ggplot(psm_df, ggplot2::aes(x = RT_val, fill = qc_group)) +
    ggplot2::geom_histogram(bins = n_bins, alpha = 0.7,
                            position = "identity", color = "white", linewidth = 0.1) +
    ggplot2::facet_wrap(~ qc_sample, scales = "free_y") +
    ggplot2::scale_fill_manual(values = mspx_colors(n_groups)) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "PSM Identifications over Retention Time",
                  x = "Retention Time (min)", y = "PSM Count", fill = "Group") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 7),
      strip.text = ggplot2::element_text(size = 7)
    )

  save_plot_and_data(p, data.frame(Sample = psm_df$qc_sample,
                                    Group = psm_df$qc_group,
                                    RT = psm_df$RT_val),
                     project_name, "QC_RT_Distribution",
                     output_dir = output_dir, width = 12, height = 8)
  p
}


# ==============================================================================
# 11. Mass Error (ppm) distribution
# ==============================================================================

#' Plot mass error (ppm) distribution
#'
#' Density plot of mass accuracy in ppm.
#' Requires "Calibrated Observed Mass" and "Calculated Peptide Mass" columns.
#' Currently available for FragPipe only.
#'
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame
#' @param output_dir Output directory
#' @param project_name Project name
#' @return ggplot object or NULL
#' @export
plot_qc_mass_error <- function(ms_data, group_info = NULL,
                                output_dir = "output",
                                project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  psm_df <- ms_data$psms
  if (is.null(psm_df) || nrow(psm_df) == 0) {
    message(">>> No PSM data. Skipping mass error plot.")
    return(invisible(NULL))
  }

  obs_col <- .find_column(psm_df, c("Calibrated Observed Mass", "Observed Mass"))
  theo_col <- .find_column(psm_df, c("Calculated Peptide Mass", "Calculated Mass"))

  if (is.null(obs_col) || is.null(theo_col)) {
    message(">>> No mass columns found. Skipping mass error plot.")
    return(invisible(NULL))
  }

  psm_df <- .psm_with_groups(psm_df, group_info, ms_data$sample_names)
  obs_mass <- as.numeric(psm_df[[obs_col]])
  theo_mass <- as.numeric(psm_df[[theo_col]])
  ppm <- (obs_mass - theo_mass) / theo_mass * 1e6
  psm_df$ppm <- ppm
  psm_df <- psm_df[!is.na(psm_df$ppm) & is.finite(psm_df$ppm), ]
  if (nrow(psm_df) == 0) return(invisible(NULL))

  med_ppm <- stats::median(psm_df$ppm, na.rm = TRUE)
  n_groups <- length(unique(psm_df$qc_group))

  p <- ggplot2::ggplot(psm_df, ggplot2::aes(x = ppm, fill = qc_group)) +
    ggplot2::geom_histogram(bins = 80, alpha = 0.7, position = "identity",
                            color = "white", linewidth = 0.1) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
    ggplot2::geom_vline(xintercept = med_ppm, linetype = "solid",
                        color = "#E64B35", linewidth = 0.8) +
    ggplot2::annotate("text", x = med_ppm, y = Inf, vjust = 2, hjust = -0.1,
                      label = sprintf("Median = %.2f ppm", med_ppm),
                      color = "#E64B35", fontface = "bold", size = 3.5) +
    ggplot2::facet_wrap(~ qc_sample, scales = "free_y") +
    ggplot2::scale_fill_manual(values = mspx_colors(n_groups)) +
    ggplot2::coord_cartesian(xlim = c(-20, 20)) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Mass Error Distribution",
                  x = "Mass Error (ppm)", y = "PSM Count", fill = "Group") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 7),
      strip.text = ggplot2::element_text(size = 7)
    )

  save_plot_and_data(p, data.frame(Sample = psm_df$qc_sample,
                                    Group = psm_df$qc_group,
                                    ppm = psm_df$ppm),
                     project_name, "QC_Mass_Error",
                     output_dir = output_dir, width = 12, height = 8)
  p
}


# ==============================================================================
# 12. Missing Value Pattern (protein level)
# ==============================================================================

#' Plot missing value pattern
#'
#' Bar chart of missing value percentage per sample + heatmap of missing pattern.
#' Uses protein-level abundance matrix (works for all engines).
#'
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame
#' @param output_dir Output directory
#' @param project_name Project name
#' @return ggplot object or NULL
#' @export
plot_qc_missing_values <- function(ms_data, group_info = NULL,
                                    output_dir = "output",
                                    project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  prot_mat <- ms_data$proteins
  if (is.null(prot_mat) || nrow(prot_mat) == 0) {
    message(">>> No protein data. Skipping missing value plot.")
    return(invisible(NULL))
  }

  # Compute per-sample missing %
  total_prots <- nrow(prot_mat)
  is_log2 <- isTRUE(ms_data$is_log2)

  missing_pct <- sapply(colnames(prot_mat), function(s) {
    vals <- prot_mat[, s]
    if (is_log2) {
      n_missing <- sum(is.na(vals) | is.nan(vals))
    } else {
      n_missing <- sum(is.na(vals) | is.nan(vals) | vals == 0)
    }
    round(n_missing / total_prots * 100, 2)
  })

  plot_df <- data.frame(
    Sample = names(missing_pct),
    Missing_Pct = as.numeric(missing_pct),
    stringsAsFactors = FALSE
  )

  # Add group info
  if (!is.null(group_info)) {
    plot_df$Group <- group_info$user_group[match(plot_df$Sample,
                                                  group_info$sample_name)]
  } else {
    plot_df$Group <- "All"
  }
  plot_df$Group[is.na(plot_df$Group)] <- "Unknown"

  n_groups <- length(unique(plot_df$Group))
  plot_df$Sample <- factor(plot_df$Sample, levels = plot_df$Sample)

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = Sample, y = Missing_Pct,
                                               fill = Group)) +
    ggplot2::geom_bar(stat = "identity", color = "white", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(round(Missing_Pct, 1), "%")),
                       vjust = -0.3, size = 3) +
    ggplot2::scale_fill_manual(values = mspx_colors(n_groups)) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Missing Value Percentage per Sample",
                  x = NULL, y = "Missing Values (%)",
                  subtitle = sprintf("Total proteins: %d", total_prots)) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8)
    )

  save_plot_and_data(p, plot_df, project_name, "QC_Missing_Values",
                     output_dir = output_dir, width = 10, height = 6)
  p
}


# ==============================================================================
# 13. Sample Intensity Distribution (protein level)
# ==============================================================================

#' Plot sample intensity distribution
#'
#' Box + violin plot of log2 protein intensities per sample.
#' Useful for assessing normalization effects and outlier samples.
#'
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame
#' @param output_dir Output directory
#' @param project_name Project name
#' @return ggplot object or NULL
#' @export
plot_qc_intensity_boxplot <- function(ms_data, group_info = NULL,
                                       output_dir = "output",
                                       project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  prot_mat <- ms_data$proteins
  if (is.null(prot_mat) || nrow(prot_mat) == 0) {
    message(">>> No protein data. Skipping intensity boxplot.")
    return(invisible(NULL))
  }

  is_log2 <- isTRUE(ms_data$is_log2)

  # Reshape to long format
  long_list <- list()
  for (s in colnames(prot_mat)) {
    vals <- prot_mat[, s]
    if (!is_log2) {
      vals[vals == 0] <- NA
      vals <- log2(vals)
    }
    vals <- vals[!is.na(vals) & is.finite(vals)]
    if (length(vals) == 0) next
    long_list[[s]] <- data.frame(
      Sample = s,
      Log2_Intensity = vals,
      stringsAsFactors = FALSE
    )
  }

  if (length(long_list) == 0) return(invisible(NULL))
  long_df <- do.call(rbind, long_list)

  # Add group info
  if (!is.null(group_info)) {
    long_df$Group <- group_info$user_group[match(long_df$Sample,
                                                  group_info$sample_name)]
  } else {
    long_df$Group <- "All"
  }
  long_df$Group[is.na(long_df$Group)] <- "Unknown"

  n_groups <- length(unique(long_df$Group))
  long_df$Sample <- factor(long_df$Sample,
                           levels = unique(long_df$Sample))

  p <- ggplot2::ggplot(long_df,
                        ggplot2::aes(x = Sample, y = Log2_Intensity,
                                     fill = Group)) +
    ggplot2::geom_violin(alpha = 0.3, scale = "width", width = 0.8) +
    ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.3, alpha = 0.8) +
    ggplot2::scale_fill_manual(values = mspx_colors(n_groups)) +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Sample Intensity Distribution",
                  x = NULL, y = "log2(Intensity)") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8)
    )

  save_plot_and_data(p,
                     data.frame(Sample = long_df$Sample,
                                Group = long_df$Group,
                                Log2_Intensity = long_df$Log2_Intensity),
                     project_name, "QC_Intensity_Distribution",
                     output_dir = output_dir, width = 10, height = 6)
  p
}


# ==============================================================================
# 14. QC Panel — all-in-one
# ==============================================================================

#' Generate QC evaluation panel (all QC plots)
#'
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame (NULL = overall)
#' @param output_dir Output directory
#' @param project_name Project name
#' @return List of generated ggplot objects (invisible)
#' @export
plot_qc_panel <- function(ms_data, group_info = NULL,
                          output_dir = "output",
                          project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  message(">>> Generating QC evaluation plots...")
  plots <- list()
  .try <- function(name, expr) {
    tryCatch(expr, error = function(e) {
      message(paste("  !", name, "failed:", e$message)); NULL
    })
  }

  plots$peptide_length <- .try("Peptide length",
    plot_qc_peptide_length(ms_data, group_info, output_dir, project_name))

  plots$charge <- .try("Charge distribution",
    plot_qc_charge_distribution(ms_data, group_info, output_dir, project_name))

  plots$missed_cleavage <- .try("Missed cleavage",
    plot_qc_missed_cleavage(ms_data, group_info, output_dir, project_name))

  plots$modification <- .try("Modification",
    plot_qc_modification(ms_data, output_dir = output_dir,
                         project_name = project_name))

  plots$gravy <- .try("GRAVY",
    plot_qc_gravy(ms_data, group_info, output_dir, project_name))

  plots$pi <- .try("pI",
    plot_qc_pi(ms_data, group_info, output_dir, project_name))

  plots$mz <- .try("M/Z",
    plot_qc_mz(ms_data, group_info, output_dir, project_name))

  plots$cumulative <- .try("Cumulative intensity",
    plot_qc_cumulative_intensity(ms_data, group_info, output_dir, project_name))

  plots$cys_alk <- .try("Cys/Alkylation",
    plot_qc_cys_alkylation(ms_data, group_info, output_dir, project_name))

  plots$rt_dist <- .try("RT distribution",
    plot_qc_rt_distribution(ms_data, group_info, output_dir, project_name))

  plots$mass_error <- .try("Mass error",
    plot_qc_mass_error(ms_data, group_info, output_dir, project_name))

  plots$missing_val <- .try("Missing values",
    plot_qc_missing_values(ms_data, group_info, output_dir, project_name))

  plots$intensity <- .try("Intensity distribution",
    plot_qc_intensity_boxplot(ms_data, group_info, output_dir, project_name))

  n_ok <- sum(!sapply(plots, is.null))
  message(sprintf(">>> QC panel complete: %d/%d plots generated.", n_ok, length(plots)))
  invisible(plots)
}


# ==============================================================================
# 10. Protein abundance rank plot (by group)
# ==============================================================================

#' Plot protein abundance rank
#'
#' Scatter plot of proteins ranked by mean intensity per group.
#' Useful for assessing dynamic range.
#'
#' @param ms_data MsDataSet object
#' @param group_info Group info data.frame
#' @param output_dir Output directory
#' @param project_name Project name
#' @return ggplot object or NULL
#' @export
plot_qc_protein_rank <- function(ms_data, group_info = NULL,
                                  output_dir = "output",
                                  project_name = "Project") {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ensure_output_dir(output_dir)

  prot_mat <- ms_data$proteins
  if (is.null(prot_mat) || nrow(prot_mat) == 0) {
    message(">>> No protein data. Skipping rank plot.")
    return(invisible(NULL))
  }

  if (is.null(group_info)) {
    mean_vals <- rowMeans(prot_mat, na.rm = TRUE)
    mean_vals <- mean_vals[!is.na(mean_vals) & mean_vals > 0]
    mean_vals <- sort(mean_vals, decreasing = TRUE)
    plot_df <- data.frame(
      Rank = seq_along(mean_vals),
      Intensity = mean_vals,
      Group = "All",
      stringsAsFactors = FALSE
    )
  } else {
    plot_list <- list()
    for (grp in unique(group_info$user_group)) {
      grp_samps <- group_info$sample_name[group_info$user_group == grp]
      grp_samps <- intersect(grp_samps, colnames(prot_mat))
      if (length(grp_samps) == 0) next
      mean_vals <- rowMeans(prot_mat[, grp_samps, drop = FALSE], na.rm = TRUE)
      mean_vals <- mean_vals[!is.na(mean_vals) & mean_vals > 0]
      mean_vals <- sort(mean_vals, decreasing = TRUE)
      plot_list[[grp]] <- data.frame(
        Rank = seq_along(mean_vals),
        Intensity = mean_vals,
        Group = grp,
        stringsAsFactors = FALSE
      )
    }
    plot_df <- do.call(rbind, plot_list)
  }

  if (nrow(plot_df) == 0) return(invisible(NULL))

  is_log2 <- isTRUE(ms_data$is_log2)
  if (is_log2) {
    plot_df$Intensity_plot <- 2^plot_df$Intensity
  } else {
    plot_df$Intensity_plot <- plot_df$Intensity
  }

  n_groups <- length(unique(plot_df$Group))
  colors <- c("#E64B35", "#4DBBD5", "#00A087", "#F39B7F",
              "#8491B4", "#91D1C2", "#DC6B5A", "#7E6148")
  if (n_groups > length(colors)) {
    colors <- grDevices::colorRampPalette(colors)(n_groups)
  }

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(
    x = Rank, y = Intensity_plot, color = Group
  )) +
    ggplot2::geom_point(size = 0.5, alpha = 0.6) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_color_manual(values = colors[seq_len(n_groups)]) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = "Protein Abundance Rank",
      x = "Protein Rank",
      y = "Intensity (log10 scale)"
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )

  save_plot_and_data(p, plot_df[, c("Rank", "Intensity", "Group")],
                     project_name, "QC_Protein_Rank",
                     output_dir = output_dir, width = 9, height = 6)
  p
}

# ==============================================================================
# Internal helpers
# ==============================================================================

#' Find a column by trying multiple candidate names
#' @keywords internal
.find_column <- function(df, candidates) {
  cols <- colnames(df)
  for (cand in candidates) {
    if (cand %in% cols) return(cand)
    # Case-insensitive fallback
    idx <- which(tolower(cols) == tolower(cand))
    if (length(idx) > 0) return(cols[idx[1]])
  }
  NULL
}
