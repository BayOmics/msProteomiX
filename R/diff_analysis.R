# ==============================================================================
# msProteomiX — 差异分析模块
# ==============================================================================

#' 运行差异分析
#'
#' 支持 Limma 和 t-test 两种方法，统一接口替代原来的 6 个火山图脚本变体。
#'
#' @param ms_data MsDataSet 对象
#' @param group_info 分组信息 data.frame
#' @param method 差异分析方法: "limma" (推荐) 或 "ttest"
#' @param p_cutoff P 值阈值 (默认 0.05)
#' @param fc_cutoff log2 Fold Change 阈值 (默认 1，即 2 倍变化)
#' @param p_type P 值类型: "raw" 或 "adj" (BH 校正)
#' @param filter_pct 最少非缺失值比例 (默认 0，不过滤)
#' @param ref_group 对照组名称 (NULL 表示交互选择)
#' @param test_group 实验组名称 (NULL 表示交互选择)
#' @return data.frame，包含差异分析结果
#' @export
run_diff_analysis <- function(ms_data, group_info,
                               method = "limma",
                               p_cutoff = 0.05,
                               fc_cutoff = 1,
                               p_type = "raw",
                               filter_pct = 0,
                               ref_group = NULL,
                               test_group = NULL) {
  stopifnot(inherits(ms_data, "MsDataSet"))

  # 构建表达矩阵
  expr_matrix <- ms_data$proteins
  colnames(expr_matrix) <- ms_data$sample_names
  meta_data <- ms_data$protein_info

  # 匹配分组
  valid_samples <- intersect(ms_data$sample_names, group_info$sample_name)
  if (length(valid_samples) == 0) stop("No matching samples between data and group info.")

  grp_vec <- group_info$user_group[match(valid_samples, group_info$sample_name)]
  expr_sub <- expr_matrix[, valid_samples, drop = FALSE]

  # Log2 转换
  expr_sub[expr_sub == 0] <- NA
  expr_log2 <- log2(as.matrix(expr_sub))

  # 分组因子
  grps <- factor(grp_vec)
  levels(grps) <- make.names(levels(grps))
  lvl <- levels(grps)

  # 交互选择对照组和实验组
  if (is.null(ref_group) || is.null(test_group)) {
    if (!interactive()) stop("ref_group and test_group must be specified in non-interactive mode.")
    selected <- .interactive_contrast_selection(lvl)
    ref_group <- selected$ref
    test_group <- selected$test
  } else {
    ref_group <- make.names(ref_group)
    test_group <- make.names(test_group)
  }

  message(paste0(">>> Running: ", test_group, " vs ", ref_group, " (", method, ")"))

  # Pairwise 提取
  target_cols <- grps %in% c(ref_group, test_group)
  expr_pair <- expr_log2[, target_cols, drop = FALSE]
  grps_pair <- droplevels(grps[target_cols])

  # 过滤低质量行
  if (filter_pct > 0) {
    keep <- apply(expr_pair, 1, function(x) sum(!is.na(x)) / length(x) >= filter_pct)
    expr_pair <- expr_pair[keep, , drop = FALSE]
    meta_data <- meta_data[keep, , drop = FALSE]
  }

  # 至少一组有 >=2 个有效值
  keep_rows <- rep(FALSE, nrow(expr_pair))
  for (g in levels(grps_pair)) {
    keep_rows <- keep_rows | (rowSums(!is.na(expr_pair[, grps_pair == g, drop = FALSE])) >= 2)
  }
  expr_pair <- expr_pair[keep_rows, , drop = FALSE]
  meta_sub <- meta_data[keep_rows, , drop = FALSE]

  if (nrow(expr_pair) == 0) stop("No valid proteins for this contrast.")

  # 差异分析
  if (method == "limma") {
    result <- .run_limma(expr_pair, grps_pair, test_group, ref_group)
  } else if (method == "ttest") {
    result <- .run_ttest(expr_pair, grps_pair, test_group, ref_group)
  } else {
    stop("Unsupported method: ", method, ". Use 'limma' or 'ttest'.")
  }

  # 合并元数据
  df <- cbind(meta_sub, result)

  # 标记差异
  p_col <- if (p_type == "adj" && "adj.P.Val" %in% colnames(df)) "adj.P.Val" else "P.Value"
  df$diff <- "NO"
  df$diff[df$logFC > fc_cutoff  & df[[p_col]] < p_cutoff] <- "UP"
  df$diff[df$logFC < -fc_cutoff & df[[p_col]] < p_cutoff] <- "DOWN"
  df$diff[is.na(df$diff)] <- "NO"
  df$diff <- factor(df$diff, levels = c("UP", "DOWN", "NO"))

  # 添加属性
  attr(df, "contrast") <- paste0(test_group, "_vs_", ref_group)
  attr(df, "method") <- method
  attr(df, "p_type") <- p_type
  attr(df, "p_cutoff") <- p_cutoff
  attr(df, "fc_cutoff") <- fc_cutoff

  # 统计
  n_up <- sum(df$diff == "UP", na.rm = TRUE)
  n_down <- sum(df$diff == "DOWN", na.rm = TRUE)
  message(sprintf("  Results: UP = %d, DOWN = %d, Total = %d", n_up, n_down, nrow(df)))

  df
}


# ==============================================================================
# 内部函数
# ==============================================================================

#' Limma 差异分析
#' @keywords internal
.run_limma <- function(expr_mat, grps, test, ref) {
  if (!requireNamespace("limma", quietly = TRUE)) {
    stop("Please install limma: BiocManager::install('limma')")
  }

  design <- stats::model.matrix(~0 + grps)
  colnames(design) <- levels(grps)

  fit <- limma::lmFit(expr_mat, design)
  contrast_str <- paste(test, "-", ref)
  contrast.matrix <- limma::makeContrasts(contrasts = contrast_str, levels = design)
  fit2 <- limma::contrasts.fit(fit, contrast.matrix)
  fit2 <- limma::eBayes(fit2, trend = TRUE)
  limma::topTable(fit2, number = Inf, sort.by = "none")
}

#' t-test 差异分析
#' @keywords internal
.run_ttest <- function(expr_mat, grps, test, ref) {
  test_idx <- grps == test
  ref_idx <- grps == ref

  results <- data.frame(
    logFC   = numeric(nrow(expr_mat)),
    P.Value = numeric(nrow(expr_mat))
  )

  for (i in seq_len(nrow(expr_mat))) {
    x <- as.numeric(expr_mat[i, test_idx])
    y <- as.numeric(expr_mat[i, ref_idx])
    x <- x[!is.na(x)]; y <- y[!is.na(y)]

    results$logFC[i] <- mean(x, na.rm = TRUE) - mean(y, na.rm = TRUE)

    if (length(x) >= 2 && length(y) >= 2) {
      tt <- tryCatch(stats::t.test(x, y), error = function(e) NULL)
      results$P.Value[i] <- if (!is.null(tt)) tt$p.value else NA_real_
    } else {
      results$P.Value[i] <- NA_real_
    }
  }

  results$adj.P.Val <- stats::p.adjust(results$P.Value, method = "BH")
  results
}

#' 交互式对比组选择
#' @keywords internal
.interactive_contrast_selection <- function(lvl) {
  message("\n==========================================================")
  message("  Differential Analysis: Select Groups")
  message("==========================================================")
  message("Available groups:")
  for (i in seq_along(lvl)) cat(sprintf("  [%d] %s\n", i, lvl[i]))

  repeat {
    ref_idx <- suppressWarnings(as.numeric(
      readline(">>> Enter CONTROL group number (0 to exit): ")))
    if (!is.na(ref_idx) && (ref_idx == 0 || (ref_idx >= 1 && ref_idx <= length(lvl)))) break
  }
  if (ref_idx == 0) stop("User cancelled.")
  ref <- lvl[ref_idx]

  remain <- lvl[-ref_idx]
  message(paste0("\n  Control group: ", ref))
  message("  Select TREATMENT group:")
  for (i in seq_along(remain)) cat(sprintf("  [%d] %s\n", i, remain[i]))

  repeat {
    test_idx <- suppressWarnings(as.numeric(
      readline(">>> Enter TREATMENT group number: ")))
    if (!is.na(test_idx) && test_idx >= 1 && test_idx <= length(remain)) break
  }

  list(ref = ref, test = remain[test_idx])
}
