# ==============================================================================
# msProteomiX — 色谱柱灵敏度对比分析
# TimsTOF Pro DDA-PASEF 纳升色谱柱比较专用函数
# ==============================================================================
# 实验设计：
#   BSA+iRT：4 浓度 × 2 柱 × 2 重复 = 16 个 .d 文件
#   HeLa：   3 上样量 × 2 柱 × 2 重复 = 12 个 .d 文件
#
# 依赖的上游工具：
#   AlphaTims (Python) -> irt_metrics.csv + bsa_metrics.csv
#   FragPipe (CLI)     -> fragpipe_out/BSA_all/ + fragpipe_out/HeLa_all/
# ==============================================================================


# ==============================================================================
# [1] import_alphatims_results()
# ==============================================================================

#' \u5bfc\u5165 AlphaTims \u5206\u6790\u7ed3\u679c
#'
#' \u8bfb\u53d6 AlphaTims Python \u811a\u672c\u8f93\u51fa\u7684 BSA \u7279\u5f81\u80bd\u548c iRT \u6307\u6807 CSV\uff0c
#' \u8fd4\u56de\u5e26\u5e26\u5206\u7ec4\u5143\u4fe1\u606f\u7684\u6807\u51c6\u5316 data.frame\u3002
#'
#' @param bsa_csv BSA \u7279\u5f81\u80bd CSV \u8def\u5f84 (bsa_metrics.csv)
#' @param irt_csv iRT \u7387\u5b9a CSV \u8def\u5f84 (irt_metrics.csv)
#' @return \u5217\u8868\uff0c\u542b\u4e24\u4e2a data.frame\uff1a$bsa \u548c $irt
#' @export
import_alphatims_results <- function(bsa_csv = "bsa_metrics.csv",
                                      irt_csv = "irt_metrics.csv") {
  result <- list(bsa = data.frame(), irt = data.frame())

  if (file.exists(bsa_csv)) {
    df <- read.csv(bsa_csv, stringsAsFactors = FALSE)
    df <- .alphatims_add_groups(df)
    result$bsa <- df
    message(sprintf(">>> BSA metrics: %d rows, %d samples",
                    nrow(df), length(unique(df$sample))))
  } else {
    message(sprintf(">>> BSA CSV not found: %s (run AlphaTims script first)", bsa_csv))
  }

  if (file.exists(irt_csv)) {
    df <- read.csv(irt_csv, stringsAsFactors = FALSE)
    df <- .alphatims_add_groups(df)
    result$irt <- df
    message(sprintf(">>> iRT metrics: %d rows, %d samples",
                    nrow(df), length(unique(df$sample))))
  } else {
    message(sprintf(">>> iRT CSV not found: %s (run AlphaTims script first)", irt_csv))
  }

  result
}


#' @keywords internal
.alphatims_add_groups <- function(df) {
  # amount \u8f6c\u6570\u5b57
  if ("amount" %in% colnames(df)) {
    df$amount <- suppressWarnings(as.numeric(df$amount))
  }
  # \u786e\u4fdd column \u5b57\u6bb5\u5b58\u5728
  if (!"column" %in% colnames(df) && "sample" %in% colnames(df)) {
    df$column <- ifelse(grepl("ColA", df$sample, ignore.case = TRUE), "ColA",
                 ifelse(grepl("ColB", df$sample, ignore.case = TRUE), "ColB", "Unknown"))
  }
  df
}


# ==============================================================================
# [2] merge_fragpipe_runs()
# ==============================================================================

#' \u5408\u5e76\u591a\u7ec4 FragPipe \u641c\u5e93\u7ed3\u679c
#'
#' \u5c06 BSA \u548c HeLa \u4e24\u4e2a\u72ec\u7acb\u7684 FragPipe \u641c\u5e93\u76ee\u5f55\u5206\u522b\u89e3\u6790\u5e76\u5b58\u5165\u5217\u8868\u3002
#'
#' @param bsa_dir BSA FragPipe \u8f93\u51fa\u76ee\u5f55 (fragpipe_out/BSA_all)
#' @param hela_dir HeLa FragPipe \u8f93\u51fa\u76ee\u5f55 (fragpipe_out/HeLa_all)
#' @return \u5217\u8868\uff0c\u542b $bsa (MsDataSet) \u548c $hela (MsDataSet)
#' @export
merge_fragpipe_runs <- function(bsa_dir = "fragpipe_out/BSA_all",
                                 hela_dir = "fragpipe_out/HeLa_all") {
  result <- list(bsa = NULL, hela = NULL)

  if (dir.exists(bsa_dir)) {
    message(">>> Loading BSA FragPipe results...")
    result$bsa <- parse_fragpipe(bsa_dir)
    message(sprintf("    Proteins: %d, Samples: %d",
                    nrow(result$bsa$proteins), length(result$bsa$sample_names)))
  } else {
    message(sprintf(">>> BSA dir not found: %s", bsa_dir))
  }

  if (dir.exists(hela_dir)) {
    message(">>> Loading HeLa FragPipe results...")
    result$hela <- parse_fragpipe(hela_dir)
    message(sprintf("    Proteins: %d, Samples: %d",
                    nrow(result$hela$proteins), length(result$hela$sample_names)))
  } else {
    message(sprintf(">>> HeLa dir not found: %s", hela_dir))
  }

  result
}


# ==============================================================================
# [3] calc_lod()
# ==============================================================================

#' \u8ba1\u7b97 LOD\uff08\u6700\u4f4e\u53ef\u68c0\u51fa\u6d53\u5ea6\uff09
#'
#' \u57fa\u4e8e AlphaTims \u7684 BSA \u7279\u5f81\u80bd\u68c0\u51fa\u6570\u636e\uff0c\u8ba1\u7b97\u6bcf\u6839\u8272\u8c31\u67f1
#' \u7684 LOD\uff08\u5b9a\u4e49\u4e3a\uff1a\u5728 2 \u4e2a\u91cd\u590d\u4e2d\u5747\u68c0\u51fa\u81f3\u5c11\u4e00\u4e2a\u7279\u5f81\u80bd\u7684\u6700\u4f4e\u6d53\u5ea6\uff09\u3002
#'
#' @param bsa_df AlphaTims BSA metrics data.frame
#' @param min_peptides \u5b9a\u4e49\u68c0\u51fa\u6240\u9700\u6700\u5c11\u80bd\u6570 (default: 1)
#' @return data.frame\uff0c\u542b\u5206\u7ec4\u548c LOD \u5c45\u5ea6
#' @export
calc_lod <- function(bsa_df, min_peptides = 1) {
  if (is.null(bsa_df) || nrow(bsa_df) == 0) {
    message(">>> No BSA data available.")
    return(data.frame())
  }

  # 每个样本的检出肽数
  detected_rows <- bsa_df[bsa_df$detected, ]
  if (nrow(detected_rows) == 0) {
    message(">>> No BSA peptides detected in any sample. LOD cannot be determined.")
    return(data.frame(column = character(), lod = numeric(),
                      lod_unit = character(), n_amounts = integer(),
                      stringsAsFactors = FALSE))
  }

  per_sample <- stats::aggregate(
    detected ~ sample + column + amount + unit,
    data = detected_rows,
    FUN  = sum
  )

  # \u6bcf\u4e2a\u5206\u7ec4\uff08column + amount\uff09\u7684\u5e73\u5747\u68c0\u51fa\u80bd\u6570
  per_group <- stats::aggregate(
    detected ~ column + amount + unit,
    data = per_sample,
    FUN  = mean
  )

  # LOD = \u5e73\u5747\u68c0\u51fa\u80bd\u6570 >= min_peptides \u7684\u6700\u4f4e\u6d53\u5ea6
  lod_list <- lapply(unique(per_group$column), function(col) {
    sub <- per_group[per_group$column == col, ]
    sub <- sub[order(sub$amount), ]
    detected_amts <- sub$amount[sub$detected >= min_peptides]
    lod_val <- if (length(detected_amts) > 0) min(detected_amts) else NA_real_
    data.frame(
      column    = col,
      lod       = lod_val,
      lod_unit  = if (nrow(sub) > 0) sub$unit[1] else "",
      n_amounts = nrow(sub),
      stringsAsFactors = FALSE
    )
  })

  lod_df <- do.call(rbind, lod_list)
  message(sprintf(">>> LOD (>= %d peptide):", min_peptides))
  for (i in seq_len(nrow(lod_df))) {
    message(sprintf("    %s: %.3g %s", lod_df$column[i],
                    lod_df$lod[i], lod_df$lod_unit[i]))
  }
  lod_df
}


# ==============================================================================
# [4] plot_lod_curve()
# ==============================================================================

#' \u7ed8\u5236 BSA LOD \u66f2\u7ebf
#'
#' \u4ee5\u6d53\u5ea6\u4e3a X \u8f74\uff0c\u68c0\u51fa\u80bd\u6bb5\u6570\u4e3a Y \u8f74\uff0c\u5bf9\u6bd4 ColA vs ColB\u3002
#'
#' @param bsa_df AlphaTims BSA metrics data.frame
#' @param output_dir \u8f93\u51fa\u76ee\u5f55
#' @param project_name \u9879\u76ee\u540d
#' @param peptide_type \u5206\u6790\u5bf9\u8c61 (\u53ef\u9009 "bsa" \u6216 "irt")
#' @return ggplot \u5bf9\u8c61
#' @export
plot_lod_curve <- function(bsa_df,
                            output_dir   = "output",
                            project_name = "Project",
                            peptide_type = "bsa") {
  if (is.null(bsa_df) || nrow(bsa_df) == 0) {
    message(">>> No BSA data available for LOD curve.")
    return(invisible(NULL))
  }
  ensure_output_dir(output_dir)

  # \u6bcf\u4e2a\u6837\u672c\u7684\u68c0\u51fa\u80bd\u6570
  per_sample <- stats::aggregate(
    detected ~ sample + column + amount + unit,
    data = bsa_df[!is.na(bsa_df$detected) & bsa_df$detected, ],
    FUN  = sum
  )
  colnames(per_sample)[colnames(per_sample) == "detected"] <- "n_detected"

  # 每组均値 + SD
  agg_mean <- stats::aggregate(n_detected ~ column + amount + unit,
                                data = per_sample, FUN = mean)
  agg_sd   <- stats::aggregate(n_detected ~ column + amount + unit,
                                data = per_sample, FUN = stats::sd)
  # 使用列名而非位置索引，避免 aggregate 输出列顺序变化
  names(agg_mean)[names(agg_mean) == "n_detected"] <- "mean_n"
  names(agg_sd)[names(agg_sd)     == "n_detected"] <- "sd_n"
  plot_df <- merge(agg_mean, agg_sd[, c("column", "amount", "sd_n")],
                   by = c("column", "amount"))
  plot_df$sd_n[is.na(plot_df$sd_n)] <- 0

  unit_label <- if (nrow(plot_df) > 0) plot_df$unit[1] else "fmol"
  n_peps     <- length(unique(bsa_df$sequence))
  cols       <- c("ColA" = "#377EB8", "ColB" = "#E41A1C")

  p <- ggplot2::ggplot(plot_df,
         ggplot2::aes(x = amount, y = mean_n, color = column, group = column)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = mean_n - sd_n, ymax = mean_n + sd_n),
      width = 0.05, linewidth = 0.8
    ) +
    ggplot2::scale_x_log10(
      breaks = sort(unique(plot_df$amount)),
      labels = scales::label_number(drop0trailing = TRUE)
    ) +
    ggplot2::scale_y_continuous(limits = c(0, n_peps + 0.5),
                                 breaks = seq(0, n_peps, 2)) +
    ggplot2::scale_color_manual(values = cols, name = "Column") +
    ggplot2::labs(
      title    = paste0("BSA LOD Curve \u2014 ", project_name),
      subtitle = paste0(n_peps, " signature peptides | error bars = SD (n=2)"),
      x        = paste0("Amount (", unit_label, ", log scale)"),
      y        = "Detected Peptides (mean)"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      legend.position   = "top",
      panel.grid.minor  = ggplot2::element_blank()
    )

  save_plot_and_data(p, plot_df, project_name, "30_BSA_LOD_Curve",
                     output_dir, width = 7, height = 5)
  p
}


# ==============================================================================
# [5] plot_linearity()
# ==============================================================================

#' \u7ed8\u5236 BSA \u5f3a\u5ea6\u7ebf\u6027\u56fe\uff08log-log\uff09
#'
#' \u5bf9\u6bcf\u6839\u8272\u8c71\u67f1\u548c\u6bcf\u6761 BSA \u80bd\u6bb5\uff0c\u7ed8\u5236 log2(\u5f3a\u5ea6) ~ log2(\u6d53\u5ea6)
#' \u7684\u56de\u5f52\u76f4\u7ebf\u548c R\u00b2\u3002
#'
#' @param bsa_df AlphaTims BSA metrics data.frame
#' @param output_dir \u8f93\u51fa\u76ee\u5f55
#' @param project_name \u9879\u76ee\u540d
#' @return ggplot \u5bf9\u8c61
#' @export
plot_linearity <- function(bsa_df,
                            output_dir   = "output",
                            project_name = "Project") {
  if (is.null(bsa_df) || nrow(bsa_df) == 0) {
    message(">>> No BSA data for linearity plot.")
    return(invisible(NULL))
  }
  ensure_output_dir(output_dir)

  # 同时需要 peak_area 和 apex_intensity 都 > 0
  has_intensity <- "apex_intensity" %in% colnames(bsa_df)
  df <- bsa_df[bsa_df$detected & bsa_df$peak_area > 0, ]
  if (nrow(df) == 0) {
    message(">>> No detected BSA peptides for linearity plot.")
    return(invisible(NULL))
  }

  cols <- c("ColA" = "#377EB8", "ColB" = "#E41A1C")

  # --- Panel A: Peak Area 线性 ---
  agg_area <- stats::aggregate(peak_area ~ column + amount, data = df, FUN = mean)
  agg_area$log2_amount <- log2(agg_area$amount)
  agg_area$log2_area   <- log2(agg_area$peak_area)

  r2_area <- do.call(rbind, lapply(unique(agg_area$column), function(col) {
    sub <- agg_area[agg_area$column == col, ]
    if (nrow(sub) < 3) return(NULL)
    fit <- stats::lm(log2_area ~ log2_amount, data = sub)
    r2  <- summary(fit)$r.squared
    data.frame(column = col, r2 = round(r2, 4),
               slope  = round(stats::coef(fit)[2], 3),
               label  = sprintf("R\u00b2=%.4f slope=%.3f", r2, stats::coef(fit)[2]),
               stringsAsFactors = FALSE)
  }))

  p_area <- ggplot2::ggplot(agg_area,
           ggplot2::aes(x = log2_amount, y = log2_area, color = column, group = column)) +
    ggplot2::geom_point(size = 3) +
    ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
    ggplot2::scale_color_manual(values = cols, name = "Column") +
    ggplot2::labs(
      title    = paste0("BSA Linearity (Peak Area) \u2014 ", project_name),
      subtitle = paste0("Mean peak area across ", length(unique(df$sequence)), " BSA signature peptides"),
      x = "log2(Amount / fmol)", y = "log2(Mean Peak Area)"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(legend.position = "top")

  if (!is.null(r2_area) && nrow(r2_area) > 0) {
    x_pos <- min(agg_area$log2_amount, na.rm = TRUE)
    y_pos <- max(agg_area$log2_area,   na.rm = TRUE)
    for (i in seq_len(nrow(r2_area))) {
      y_off <- (i - 1) * diff(range(agg_area$log2_area, na.rm = TRUE)) * 0.12
      p_area <- p_area + ggplot2::annotate(
        "text", x = x_pos, y = y_pos - y_off,
        label = r2_area$label[i], color = cols[r2_area$column[i]],
        size = 3.5, hjust = 0
      )
    }
  }
  save_plot_and_data(p_area, agg_area, project_name, "31a_BSA_Linearity_Area",
                     output_dir, width = 6, height = 5)

  # --- Panel B: Apex Intensity 线性 ---
  p_int <- NULL
  if (has_intensity) {
    df_int <- df[df$apex_intensity > 0, ]
    if (nrow(df_int) > 0) {
      agg_int <- stats::aggregate(apex_intensity ~ column + amount, data = df_int, FUN = mean)
      agg_int$log2_amount    <- log2(agg_int$amount)
      agg_int$log2_intensity <- log2(agg_int$apex_intensity)

      r2_int <- do.call(rbind, lapply(unique(agg_int$column), function(col) {
        sub <- agg_int[agg_int$column == col, ]
        if (nrow(sub) < 3) return(NULL)
        fit <- stats::lm(log2_intensity ~ log2_amount, data = sub)
        r2  <- summary(fit)$r.squared
        data.frame(column = col, r2 = round(r2, 4),
                   slope  = round(stats::coef(fit)[2], 3),
                   label  = sprintf("R\u00b2=%.4f slope=%.3f", r2, stats::coef(fit)[2]),
                   stringsAsFactors = FALSE)
      }))

      p_int <- ggplot2::ggplot(agg_int,
               ggplot2::aes(x = log2_amount, y = log2_intensity, color = column, group = column)) +
        ggplot2::geom_point(size = 3) +
        ggplot2::geom_smooth(method = "lm", se = FALSE, linewidth = 1.2) +
        ggplot2::scale_color_manual(values = cols, name = "Column") +
        ggplot2::labs(
          title    = paste0("BSA Linearity (Apex Intensity) \u2014 ", project_name),
          subtitle = paste0("Mean apex intensity across ", length(unique(df_int$sequence)), " BSA signature peptides"),
          x = "log2(Amount / fmol)", y = "log2(Mean Apex Intensity)"
        ) +
        ggplot2::theme_bw(base_size = 13) +
        ggplot2::theme(legend.position = "top")

      if (!is.null(r2_int) && nrow(r2_int) > 0) {
        x_pos <- min(agg_int$log2_amount, na.rm = TRUE)
        y_pos <- max(agg_int$log2_intensity, na.rm = TRUE)
        for (i in seq_len(nrow(r2_int))) {
          y_off <- (i - 1) * diff(range(agg_int$log2_intensity, na.rm = TRUE)) * 0.12
          p_int <- p_int + ggplot2::annotate(
            "text", x = x_pos, y = y_pos - y_off,
            label = r2_int$label[i], color = cols[r2_int$column[i]],
            size = 3.5, hjust = 0
          )
        }
      }
      save_plot_and_data(p_int, agg_int, project_name, "31b_BSA_Linearity_Intensity",
                         output_dir, width = 6, height = 5)
      message("\u2705 Apex intensity linearity plot saved.")
    }
  } else {
    message(">>> apex_intensity column not found; skipping intensity linearity plot.")
  }

  invisible(list(area = p_area, intensity = p_int))
}



# ==============================================================================
# [6] plot_irt_column_compare()
# ==============================================================================

#' \u7ed8\u5236 iRT \u8272\u8c71\u67f1\u5bf9\u6bd4\u56fe\uff083-panel\uff09
#'
#' \u9762\u677f 1\uff1a\u5404\u6d53\u5ea6\u4e0b iRT \u5e73\u5747\u5cf0\u9762\u79ef\uff08ColA vs ColB\uff09
#' \u9762\u677f 2\uff1a\u5404 iRT \u80bd\u6bb5\u7684\u5e73\u5747 FWHM\uff08\u5cf0\u5bbd\uff09
#' \u9762\u677f 3\uff1a\u91cd\u590d\u95f4 CV%\uff08\u5cf0\u9762\u79ef\uff09
#'
#' @param irt_df AlphaTims iRT metrics data.frame
#' @param output_dir \u8f93\u51fa\u76ee\u5f55
#' @param project_name \u9879\u76ee\u540d
#' @return ggplot \u5bf9\u8c61 (\u7ec4\u5408\u56fe)
#' @export
plot_irt_column_compare <- function(irt_df,
                                     output_dir   = "output",
                                     project_name = "Project") {
  if (is.null(irt_df) || nrow(irt_df) == 0) {
    message(">>> No iRT data available.")
    return(invisible(NULL))
  }
  ensure_output_dir(output_dir)

  cols     <- c("ColA" = "#377EB8", "ColB" = "#E41A1C")
  df_det   <- irt_df[!is.na(irt_df$detected) & irt_df$detected, ]

  # Panel 1: \u5404\u6d53\u5ea6\u5cf0\u9762\u79ef\uff08\u4e2d\u4f4d\uff09
  agg1 <- stats::aggregate(peak_area ~ column + amount,
                            data = df_det, FUN = stats::median)
  p1 <- ggplot2::ggplot(agg1,
          ggplot2::aes(x = factor(amount), y = peak_area, fill = column)) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    ggplot2::scale_fill_manual(values = cols, name = "Column") +
    ggplot2::scale_y_continuous(labels = scales::label_scientific()) +
    ggplot2::labs(title = "iRT Median Peak Area",
                  x = "Amount (fmol)", y = "Median Peak Area") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "top")

  # Panel 2: FWHM \u5206\u80fd\u5bbd\u5ea6\uff08\u6bcf\u6761\u80bd\u6bb5\uff09
  agg2 <- stats::aggregate(fwhm ~ column + sequence,
                            data = df_det[!is.na(df_det$fwhm), ], FUN = mean)
  # \u6309\u5e8f\u5217\u6539\u7b80\u5199\u5e8f\u5217\u540d
  agg2$seq_short <- substr(agg2$sequence, 1, 8)
  p2 <- ggplot2::ggplot(agg2,
          ggplot2::aes(x = stats::reorder(seq_short, fwhm),
                       y = fwhm, fill = column)) +
    ggplot2::geom_col(position = "dodge", width = 0.7) +
    ggplot2::scale_fill_manual(values = cols, name = "Column") +
    ggplot2::geom_hline(yintercept = 0.3, linetype = "dashed",
                        color = "grey50", linewidth = 0.8) +
    ggplot2::annotate("text", x = 1.5, y = 0.32, label = "target 0.3 min",
                      size = 3, color = "grey40", hjust = 0) +
    ggplot2::labs(title = "iRT Peak Width (FWHM)",
                  x = "Peptide (truncated)", y = "FWHM (min)") +
    ggplot2::coord_flip() +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "top")

  # Panel 3: \u91cd\u590d\u95f4 CV%\uff08\u6240\u6709\u6d53\u5ea6\u5408\u5e76\uff09
  cv_df <- stats::aggregate(peak_area ~ column + sequence,
                             data = df_det, FUN = function(x) {
                               if (length(x) < 2 || mean(x) == 0) return(NA_real_)
                               stats::sd(x) / mean(x) * 100
                             })
  colnames(cv_df)[3] <- "cv_pct"
  cv_df <- cv_df[!is.na(cv_df$cv_pct), ]

  p3 <- ggplot2::ggplot(cv_df,
          ggplot2::aes(x = column, y = cv_pct, fill = column)) +
    ggplot2::geom_boxplot(alpha = 0.8, outlier.shape = 16, outlier.size = 1.5) +
    ggplot2::geom_hline(yintercept = 5, linetype = "dashed",
                        color = "grey50", linewidth = 0.8) +
    ggplot2::annotate("text", x = 0.6, y = 6.5, label = "5% threshold",
                      size = 3, color = "grey40", hjust = 0) +
    ggplot2::scale_fill_manual(values = cols, guide = "none") +
    ggplot2::labs(title = "iRT Inter-Rep CV%",
                  subtitle = "Peak area reproducibility",
                  x = "Column", y = "CV (%)") +
    ggplot2::theme_bw(base_size = 11)

  # \u7ec4\u5408 3 \u4e2a\u9762\u677f
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    combined <- gridExtra::grid.arrange(p1, p2, p3, ncol = 3)
    ggplot2::ggsave(
      file.path(output_dir, paste0(project_name, "_32_iRT_Column_Compare.pdf")),
      gridExtra::arrangeGrob(p1, p2, p3, ncol = 3),
      width = 18, height = 6
    )
    ggplot2::ggsave(
      file.path(output_dir, paste0(project_name, "_32_iRT_Column_Compare.png")),
      gridExtra::arrangeGrob(p1, p2, p3, ncol = 3),
      width = 18, height = 6, dpi = 300, bg = "white"
    )
    message("\u2705 \u5df2\u751f\u6210: iRT_Column_Compare")
    invisible(combined)
  } else {
    message(">>> gridExtra not installed. Saving panels separately...")
    save_plot_and_data(p1, agg1, project_name, "32a_iRT_PeakArea", output_dir, 7, 5)
    save_plot_and_data(p2, agg2, project_name, "32b_iRT_FWHM",     output_dir, 7, 6)
    save_plot_and_data(p3, cv_df, project_name, "32c_iRT_CV",       output_dir, 5, 5)
    invisible(p1)
  }
}


# ==============================================================================
# [7] plot_hela_loading_curve()
# ==============================================================================

#' \u7ed8\u5236 HeLa \u4e0a\u6837\u91cf vs \u86cb\u767d\u8d28\u6570\u66f2\u7ebf
#'
#' \u4ee5\u4e0a\u6837\u91cf\u4e3a X \u8f74\uff0c\u68c0\u5230\u7684\u86cb\u767d\u8d28\u6570/\u80bd\u6bb5\u6570\u4e3a Y \u8f74\uff0c\u5bf9\u6bd4 ColA vs ColB\u3002
#'
#' @param hela_ms MsDataSet \u5bf9\u8c61 (FragPipe HeLa \u641c\u5e93\u7ed3\u679c)
#' @param output_dir \u8f93\u51fa\u76ee\u5f55
#' @param project_name \u9879\u76ee\u540d
#' @return ggplot \u5bf9\u8c61
#' @export
plot_hela_loading_curve <- function(hela_ms,
                                     output_dir   = "output",
                                     project_name = "Project") {
  if (is.null(hela_ms) || !inherits(hela_ms, "MsDataSet")) {
    message(">>> No HeLa MsDataSet available.")
    return(invisible(NULL))
  }
  ensure_output_dir(output_dir)

  prot_mat <- hela_ms$proteins
  samples  <- hela_ms$sample_names

  # \u89e3\u6790\u6837\u672c\u540d\uff1a HeLa_ColA_1ng_rep1 -> column + amount
  meta <- do.call(rbind, lapply(samples, function(s) {
    m <- regmatches(s, regexpr(
      "(ColA|ColB)_([\\.0-9]+)(ng|fmol)", s, ignore.case = TRUE
    ))
    if (length(m) == 0) {
      return(data.frame(sample = s, column = "Unknown",
                        amount = NA_real_, unit = "", stringsAsFactors = FALSE))
    }
    parts <- regmatches(m, regexec(
      "(ColA|ColB)_([\\.0-9]+)(ng|fmol)", m, ignore.case = TRUE
    ))[[1]]
    data.frame(sample = s, column = parts[2],
               amount = as.numeric(parts[3]),
               unit   = parts[4], stringsAsFactors = FALSE)
  }))

  # \u6bcf\u6837\u672c\u7684\u86cb\u767d\u8d28\u6570 (intensity > 0)
  n_prot <- colSums(prot_mat > 0, na.rm = TRUE)
  meta$n_proteins <- n_prot[match(meta$sample, names(n_prot))]

  # \u5747\u5024
  agg <- stats::aggregate(n_proteins ~ column + amount + unit,
                           data = meta, FUN = mean)
  agg_sd <- stats::aggregate(n_proteins ~ column + amount,
                              data = meta, FUN = stats::sd)
  colnames(agg_sd)[3] <- "sd_n"
  plot_df <- merge(agg, agg_sd[, c("column", "amount", "sd_n")],
                   by = c("column", "amount"))
  plot_df$sd_n[is.na(plot_df$sd_n)] <- 0

  unit_label <- if (nrow(plot_df) > 0) plot_df$unit[1] else "ng"
  cols <- c("ColA" = "#377EB8", "ColB" = "#E41A1C")

  p <- ggplot2::ggplot(plot_df,
         ggplot2::aes(x = amount, y = n_proteins, color = column, group = column)) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 4) +
    ggplot2::geom_errorbar(
      ggplot2::aes(ymin = n_proteins - sd_n, ymax = n_proteins + sd_n),
      width = 0.02, linewidth = 0.8
    ) +
    ggplot2::scale_x_log10(
      breaks = sort(unique(plot_df$amount)),
      labels = scales::label_number(drop0trailing = TRUE)
    ) +
    ggplot2::scale_color_manual(values = cols, name = "Column") +
    ggplot2::scale_y_continuous(labels = scales::label_comma()) +
    ggplot2::labs(
      title    = paste0("HeLa Proteome Depth \u2014 ", project_name),
      subtitle = "Protein IDs (1% FDR) | error bars = SD (n=2)",
      x        = paste0("Loading Amount (", unit_label, ", log scale)"),
      y        = "Protein IDs"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(legend.position = "top",
                   panel.grid.minor = ggplot2::element_blank())

  save_plot_and_data(p, plot_df, project_name, "33_HeLa_Loading_Curve",
                     output_dir, width = 7, height = 5)
  p
}


# ==============================================================================
# [8] calc_column_metrics()
# ==============================================================================

#' \u8ba1\u7b97\u8272\u8c71\u67f1\u5bf9\u6bd4 KPI \u6c47\u603b\u8868
#'
#' \u6574\u5408 AlphaTims + FragPipe \u7ed3\u679c\uff0c\u751f\u6210 8 \u4e2a\u6838\u5fc3 KPI \u7684\u6c47\u603b\u8868\u683c\u3002
#'
#' @param alphatims_res import_alphatims_results() \u7684\u8fd4\u56de\u5217\u8868
#' @param fragpipe_res merge_fragpipe_runs() \u7684\u8fd4\u56de\u5217\u8868
#' @return data.frame KPI \u6c47\u603b\u8868
#' @export
calc_column_metrics <- function(alphatims_res = NULL, fragpipe_res = NULL) {
  rows <- list()

  # --- iRT KPIs ---
  if (!is.null(alphatims_res$irt) && nrow(alphatims_res$irt) > 0) {
    irt_df <- alphatims_res$irt
    for (col in unique(irt_df$column)) {
      sub <- irt_df[irt_df$column == col, ]
      # \u5728\u6700\u9ad8\u6d53\u5ea6\u4e0b\u7684\u68c0\u51fa\u7387
      if ("amount" %in% colnames(sub)) {
        max_amt <- max(sub$amount, na.rm = TRUE)
        sub_max <- sub[sub$amount == max_amt, ]
        irt_rate <- mean(sub_max$detected, na.rm = TRUE) * 11
      } else {
        irt_rate <- sum(sub$detected, na.rm = TRUE)
      }
      median_area      <- stats::median(sub$peak_area[sub$detected], na.rm = TRUE)
      median_intensity <- if ("apex_intensity" %in% colnames(sub))
        stats::median(sub$apex_intensity[sub$detected], na.rm = TRUE) else NA_real_
      median_fwhm <- stats::median(sub$fwhm[sub$detected & !is.na(sub$fwhm)], na.rm = TRUE)
      # CV%：在最高浓度下跨重复计算（避免混入浓度差异导致 CV 虚高）
      max_amt_cv <- max(sub$amount, na.rm = TRUE)
      sub_repro  <- sub[!is.na(sub$amount) & sub$amount == max_amt_cv &
                        !is.na(sub$detected) & sub$detected, ]
      cv_pct <- if (nrow(sub_repro) >= 2 && mean(sub_repro$peak_area) > 0) {
        stats::sd(sub_repro$peak_area) /
          mean(sub_repro$peak_area) * 100
      } else NA_real_
      rows[[length(rows) + 1]] <- data.frame(
        column = col, metric = "iRT_detection_rate", value = irt_rate,
        unit = "peptides/11", stringsAsFactors = FALSE)
      rows[[length(rows) + 1]] <- data.frame(
        column = col, metric = "iRT_median_peak_area", value = median_area,
        unit = "AU", stringsAsFactors = FALSE)
      rows[[length(rows) + 1]] <- data.frame(
        column = col, metric = "iRT_median_apex_intensity", value = median_intensity,
        unit = "AU", stringsAsFactors = FALSE)
      rows[[length(rows) + 1]] <- data.frame(
        column = col, metric = "iRT_median_FWHM_min", value = median_fwhm,
        unit = "min", stringsAsFactors = FALSE)
      rows[[length(rows) + 1]] <- data.frame(
        column = col, metric = "iRT_CV_pct", value = cv_pct,
        unit = "%", stringsAsFactors = FALSE)
    }
  }

  # --- BSA KPIs ---
  if (!is.null(alphatims_res$bsa) && nrow(alphatims_res$bsa) > 0) {
    bsa_df <- alphatims_res$bsa
    lod_df <- calc_lod(bsa_df)
    for (i in seq_len(nrow(lod_df))) {
      rows[[length(rows) + 1]] <- data.frame(
        column = lod_df$column[i], metric = "BSA_LOD",
        value = lod_df$lod[i],
        unit = lod_df$lod_unit[i], stringsAsFactors = FALSE)
    }
  }

  # --- HeLa KPIs ---
  if (!is.null(fragpipe_res$hela) && inherits(fragpipe_res$hela, "MsDataSet")) {
    hela_ms  <- fragpipe_res$hela
    prot_mat <- hela_ms$proteins
    samples  <- hela_ms$sample_names
    for (col in c("ColA", "ColB")) {
      subs <- samples[grepl(col, samples, ignore.case = TRUE)]
      if (length(subs) == 0) next
      sub_mat  <- prot_mat[, subs, drop = FALSE]
      n_prot   <- mean(colSums(sub_mat > 0, na.rm = TRUE))
      cv_vals  <- apply(sub_mat, 1, function(x) {
        x <- x[x > 0]
        if (length(x) < 2) return(NA_real_)
        stats::sd(x) / mean(x) * 100
      })
      rows[[length(rows) + 1]] <- data.frame(
        column = col, metric = "HeLa_mean_protein_IDs",
        value = round(n_prot), unit = "proteins",
        stringsAsFactors = FALSE)
      rows[[length(rows) + 1]] <- data.frame(
        column = col, metric = "HeLa_median_CV_pct",
        value = round(stats::median(cv_vals, na.rm = TRUE), 1),
        unit = "%", stringsAsFactors = FALSE)
    }
  }

  if (length(rows) == 0) {
    message(">>> No data available to calculate metrics.")
    return(data.frame())
  }
  kpi_df <- do.call(rbind, rows)
  kpi_df$value <- round(kpi_df$value, 4)

  # 宽格式，方便对比
  wide <- reshape2::dcast(kpi_df, metric + unit ~ column, value.var = "value")
  message(">>> KPI Summary:")
  print(wide)
  wide
}


# ==============================================================================
# [9] make_verification_table()
# ==============================================================================

#' 生成肽段核对清单
#'
#' 合并 FragPipe PSM / combined_peptide 结果与 AlphaTims EIC 指标，
#' 输出含序列、m/z、charge、RT、峰强度、峰面积的核对表。
#' 支持三种模式：
#'   "bsa"  — BSA 肽段：PSM.tsv × bsa_metrics.csv
#'   "irt"  — iRT 肽段：PSM.tsv × irt_metrics.csv
#'   "hela" — HeLa 肽段：combined_peptide.tsv（无 AlphaTims）
#'
#' @param fragpipe_dir FragPipe 输出目录（含 psm.tsv 或 combined_peptide.tsv）
#' @param at_csv       AlphaTims CSV 路径（bsa_metrics.csv 或 irt_metrics.csv）
#' @param table_type   "bsa" | "irt" | "hela"
#' @param output_dir   输出目录
#' @param project_name 项目名前缀
#' @param min_pp_prob  PeptideProphet 最低概率（默认 0.9）
#' @return 合并后的核对 data.frame（不可见返回，同时保存 CSV/Excel）
#' @export
make_verification_table <- function(fragpipe_dir,
                                     at_csv       = NULL,
                                     table_type   = c("bsa", "irt", "hela"),
                                     output_dir   = "output",
                                     project_name = "Project",
                                     min_pp_prob  = 0.9) {
  table_type <- match.arg(table_type)
  ensure_output_dir(output_dir)

  # ============================================================
  # 内部辅助：FragPipe psm.tsv 列名标准化
  # ============================================================
  .pick_col <- function(df, candidates) {
    found <- intersect(candidates, colnames(df))
    if (length(found) == 0) return(rep(NA_character_, nrow(df)))
    as.character(df[[found[1]]])
  }

  .read_psm_dir <- function(dir) {
    fs <- list.files(dir, pattern = "^psm\\.tsv$",
                     full.names = TRUE, recursive = TRUE)
    if (length(fs) == 0) stop("No psm.tsv found under: ", dir)
    message(sprintf(">>> Found %d psm.tsv file(s)", length(fs)))
    raw <- do.call(rbind, lapply(fs, function(f) {
      df <- read.delim(f, stringsAsFactors = FALSE, check.names = FALSE)
      df$._sample <- basename(dirname(f))
      df
    }))
    data.frame(
      sample        = raw$`._sample`,
      sequence      = .pick_col(raw, c("Peptide","peptide","Sequence")),
      mod_sequence  = .pick_col(raw, c("Modified Peptide","modified_peptide")),
      charge        = suppressWarnings(as.integer(
                        .pick_col(raw, c("Charge","charge")))),
      mz_observed   = suppressWarnings(as.numeric(
                        .pick_col(raw, c("Calibrated Observed M/Z",
                                         "Observed M/Z","observed_mz")))),
      mz_calculated = suppressWarnings(as.numeric(
                        .pick_col(raw, c("Calculated M/Z","calculated_mz")))),
      rt_min        = suppressWarnings(as.numeric(
                        .pick_col(raw, c("Retention","Retention Time",
                                         "retention_time_min","RT")))),
      intensity_ms1 = suppressWarnings(as.numeric(
                        .pick_col(raw, c("Intensity","intensity")))),
      pp_prob       = suppressWarnings(as.numeric(
                        .pick_col(raw, c("PeptideProphet Probability",
                                         "pep_prob","Probability")))),
      hyperscore    = suppressWarnings(as.numeric(
                        .pick_col(raw, c("Hyperscore","hyperscore")))),
      delta_ppm     = suppressWarnings(as.numeric(
                        .pick_col(raw, c("Delta Mass","delta_mass")))),
      protein       = .pick_col(raw, c("Protein","protein")),
      gene          = .pick_col(raw, c("Gene","gene")),
      stringsAsFactors = FALSE
    )
  }

  # ============================================================
  # 模式 A：BSA 或 iRT（PSM + AlphaTims）
  # ============================================================
  if (table_type %in% c("bsa", "irt")) {
    psm <- .read_psm_dir(fragpipe_dir)

    # 过滤 pp_prob
    if (!all(is.na(psm$pp_prob))) {
      n0 <- nrow(psm)
      psm <- psm[is.na(psm$pp_prob) | psm$pp_prob >= min_pp_prob, ]
      message(sprintf(">>> pp_prob >= %.2f filter: %d -> %d PSMs",
                      min_pp_prob, n0, nrow(psm)))
    }

    # 按蛋白来源过滤（避免 BSA PSM 混入 iRT 表，反之亦然）
    if (table_type == "irt" && "protein" %in% colnames(psm)) {
      n0 <- nrow(psm)
      psm <- psm[grepl("iRT|Biognosys|iRTKit", psm$protein, ignore.case = TRUE), ]
      message(sprintf(">>> iRT protein filter: %d -> %d PSMs", n0, nrow(psm)))
    } else if (table_type == "bsa" && "protein" %in% colnames(psm)) {
      n0 <- nrow(psm)
      psm <- psm[grepl("ALBU_BOVIN|P02769", psm$protein, ignore.case = TRUE), ]
      message(sprintf(">>> BSA protein filter: %d -> %d PSMs", n0, nrow(psm)))
    }

    # 去重（sequence + charge + sample）
    psm <- psm[order(psm$pp_prob, decreasing = TRUE, na.last = TRUE), ]
    psm <- psm[!duplicated(paste(psm$sequence, psm$charge, psm$sample)), ]
    message(sprintf(">>> After dedup: %d unique peptide×sample", nrow(psm)))

    # 读 AlphaTims CSV
    at_df <- data.frame()
    if (!is.null(at_csv) && file.exists(at_csv)) {
      at_df <- read.csv(at_csv, stringsAsFactors = FALSE)
      # NA-safe detected 过滤
      at_df <- at_df[!is.na(at_df$detected) & at_df$detected == TRUE, ]
      message(sprintf(">>> AlphaTims (%s): %d detected rows", table_type, nrow(at_df)))
    } else {
      message(">>> AlphaTims CSV not found — EIC columns will be NA")
    }

    # 聚合 AlphaTims（跨所有样本取中位数）
    if (nrow(at_df) > 0) {
      at_agg <- stats::aggregate(
        cbind(mz_theoretical, apex_rt, apex_intensity, peak_area, snr, fwhm) ~
          sequence + charge,
        data  = at_df,
        FUN   = function(x) round(median(x, na.rm = TRUE), 4)
      )
      # 统计检出样本数
      key <- paste(at_df$sequence, at_df$charge)
      at_agg$n_detected <- as.integer(table(key)[
        paste(at_agg$sequence, at_agg$charge)])
      result <- merge(psm, at_agg, by = c("sequence", "charge"), all.x = TRUE)
    } else {
      result <- psm
      result[c("mz_theoretical","apex_rt","apex_intensity",
               "peak_area","snr","fwhm","n_detected")] <- NA_real_
    }

    # m/z 误差（保护零分母和 NA）
    if (!all(is.na(result$mz_observed)) && !all(is.na(result$mz_theoretical))) {
      denom <- result$mz_theoretical
      denom[!is.na(denom) & denom == 0] <- NA_real_
      result$mz_error_ppm <- round(
        (result$mz_observed - denom) / denom * 1e6, 2)
    } else {
      result$mz_error_ppm <- NA_real_
    }

    # 列顺序
    col_order <- c("sequence","mod_sequence","charge",
                   "mz_observed","mz_theoretical","mz_calculated","mz_error_ppm",
                   "rt_min","apex_rt",
                   "apex_intensity","peak_area","snr","fwhm",
                   "intensity_ms1","pp_prob","hyperscore",
                   "protein","gene","n_detected","sample")
    result <- result[, intersect(col_order, colnames(result))]
    result <- result[order(result$pp_prob, decreasing = TRUE, na.last = TRUE), ]

    label <- toupper(table_type)  # "BSA" or "IRT"

  # ============================================================
  # 模式 B：HeLa（combined_peptide.tsv）
  # ============================================================
  } else {
    # 尝试 combined_peptide.tsv
    cpep <- .find_file(fragpipe_dir, "^combined_peptide\\.(tsv|csv|txt)$")
    if (is.null(cpep)) {
      # 搜子目录
      for (sd in list.dirs(fragpipe_dir, recursive = FALSE, full.names = TRUE)) {
        cpep <- .find_file(sd, "^combined_peptide\\.(tsv|csv|txt)$")
        if (!is.null(cpep)) break
      }
    }
    if (is.null(cpep))
      stop("combined_peptide.tsv not found under: ", fragpipe_dir)

    message(sprintf(">>> Reading: %s", basename(cpep)))
    raw <- read.delim(cpep, stringsAsFactors = FALSE, check.names = FALSE)

    # 提取核心注释列
    id_cols <- c("Peptide","Modified Peptide","Charge",
                 "Prev AA","Next AA","Protein","Protein ID","Gene",
                 "Organism","Calculated M/Z","Calibrated Observed M/Z",
                 "Retention","Assigned Modifications","Is Unique",
                 "PeptideProphet Probability","Hyperscore","Number of Missed Cleavages")
    id_keep <- intersect(id_cols, colnames(raw))

    # 提取各样本的 MaxLFQ Intensity 列（per-experiment）
    lfq_cols <- grep("MaxLFQ Intensity$|\\bIntensity$", colnames(raw),
                     value = TRUE, ignore.case = TRUE)
    lfq_cols <- lfq_cols[!grepl("^(Total|Combined)", lfq_cols, ignore.case = TRUE)]

    result <- raw[, c(id_keep, lfq_cols), drop = FALSE]

    # 精确列名映射（避免链式 gsub 的歧义替换）
    col_rename <- c(
      "Peptide"                    = "sequence",
      "Modified Peptide"           = "mod_sequence",
      "Charge"                     = "charge",
      "Calibrated Observed M/Z"    = "mz_observed",
      "Calculated M/Z"             = "mz_calculated",
      "Retention"                  = "rt_min",
      "PeptideProphet Probability" = "pp_prob",
      "Hyperscore"                 = "hyperscore",
      "Number of Missed Cleavages" = "n_missed",
      "Is Unique"                  = "is_unique",
      "Assigned Modifications"     = "modifications"
    )
    colnames(result) <- ifelse(
      colnames(result) %in% names(col_rename),
      col_rename[colnames(result)],
      colnames(result)
    )

    # pp_prob 过滤
    if ("pp_prob" %in% colnames(result)) {
      result$pp_prob <- suppressWarnings(as.numeric(result$pp_prob))
      n0 <- nrow(result)
      result <- result[is.na(result$pp_prob) | result$pp_prob >= min_pp_prob, ]
      message(sprintf(">>> pp_prob >= %.2f filter: %d -> %d peptides",
                      min_pp_prob, n0, nrow(result)))
    }

    # m/z 误差
    if (all(c("mz_observed","mz_calculated") %in% colnames(result))) {
      result$mz_observed   <- suppressWarnings(as.numeric(result$mz_observed))
      result$mz_calculated <- suppressWarnings(as.numeric(result$mz_calculated))
      result$mz_error_ppm  <- round(
        (result$mz_observed - result$mz_calculated) /
          result$mz_calculated * 1e6, 2)
    }

    # 按 pp_prob 排序
    if ("pp_prob" %in% colnames(result))
      result <- result[order(result$pp_prob, decreasing = TRUE, na.last = TRUE), ]

    label <- "HeLa"
    message(sprintf(">>> HeLa peptide table: %d peptides, %d intensity columns",
                    nrow(result), length(lfq_cols)))
  }

  rownames(result) <- NULL

  # ============================================================
  # 保存 CSV
  # ============================================================
  out_csv <- file.path(output_dir,
                        sprintf("%s_%s_verification_table.csv",
                                project_name, label))
  write.csv(result, out_csv, row.names = FALSE)
  message(sprintf(">>> Saved CSV: %s", basename(out_csv)))

  # ============================================================
  # 保存 Excel（如果有 openxlsx）
  # ============================================================
  if (requireNamespace("openxlsx", quietly = TRUE)) {
    out_xlsx <- sub("\\.csv$", ".xlsx", out_csv)
    wb <- openxlsx::createWorkbook()
    sheet <- paste0(label, "_Verification")
    openxlsx::addWorksheet(wb, sheet)
    openxlsx::writeData(wb, sheet, result)
    openxlsx::freezePane(wb, sheet, firstRow = TRUE)

    # 格式：表头加粗 + 背景色
    hdr_style <- openxlsx::createStyle(
      fontColour = "#FFFFFF", fgFill = "#2C5F8A",
      halign = "center", textDecoration = "bold", wrapText = TRUE)
    openxlsx::addStyle(wb, sheet, hdr_style,
                        rows = 1, cols = seq_len(ncol(result)),
                        gridExpand = TRUE)

    # 低概率行（pp_prob < 0.95）标黄
    if ("pp_prob" %in% colnames(result)) {
      low_rows <- which(!is.na(result$pp_prob) & result$pp_prob < 0.95) + 1
      if (length(low_rows) > 0) {
        openxlsx::addStyle(wb, sheet,
          openxlsx::createStyle(fgFill = "#FFF2CC"),
          rows = low_rows, cols = seq_len(ncol(result)),
          gridExpand = TRUE)
      }
    }

    # 自动列宽
    openxlsx::setColWidths(wb, sheet,
                            cols = seq_len(ncol(result)), widths = "auto")
    openxlsx::saveWorkbook(wb, out_xlsx, overwrite = TRUE)
    message(sprintf(">>> Saved Excel: %s", basename(out_xlsx)))
  } else {
    message(">>> openxlsx not available; install it for Excel output.")
  }

  invisible(result)
}
