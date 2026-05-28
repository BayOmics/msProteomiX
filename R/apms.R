# ==============================================================================
# msProteomiX — AP-MS 亲和纯化质谱分析模块
# ==============================================================================

# ---- 内嵌 CRAPome/cRAP 常见 AP-MS 污染蛋白基因名 ----
# 来源: GPM cRAP + CRAPome 高频蛋白 (出现频率 > 20% 的 AP-MS 对照实验)
.crapome_genes <- c(
  # 角蛋白 (环境污染)
  "KRT1", "KRT2", "KRT3", "KRT4", "KRT5", "KRT6A", "KRT6B", "KRT6C",
  "KRT7", "KRT8", "KRT9", "KRT10", "KRT13", "KRT14", "KRT15", "KRT16",
  "KRT17", "KRT18", "KRT19", "KRT20", "KRT24", "KRT25", "KRT26",
  "KRT27", "KRT28", "KRT31", "KRT32", "KRT33A", "KRT33B", "KRT34",
  "KRT35", "KRT36", "KRT37", "KRT38", "KRT39", "KRT40",
  "KRT71", "KRT72", "KRT73", "KRT74", "KRT75", "KRT76", "KRT77",
  "KRT78", "KRT79", "KRT80", "KRT81", "KRT82", "KRT83", "KRT84",
  "KRT85", "KRT86",
  # 血清蛋白 (来自培养基/血清)
  "ALB", "TF", "A2M", "SERPINA1", "HP", "APOA1", "APOA2", "APOB",
  "APOE", "AHSG", "FGA", "FGB", "FGG", "FN1", "VTN", "SERPINC1",
  # 酶/标准品 (实验试剂)
  "PRSS1", "PRSS2", "PRSS3",   # trypsin
  # 核糖体/翻译因子 (高丰度背景)
  "RPS27A", "UBA52", "UBB", "UBC",  # ubiquitin
  # 分子伴侣 (常见非特异性结合)
  "HSPA5", "HSPA8", "HSPA9", "HSP90AA1", "HSP90AB1", "HSP90B1",
  "HSPD1", "HSPE1", "HSPA1A", "HSPA1B", "CCT2", "CCT3", "CCT4",
  "CCT5", "CCT6A", "CCT7", "CCT8", "TCP1",
  # 细胞骨架 (高丰度)
  "ACTB", "ACTG1", "ACTA1", "ACTA2", "ACTC1", "ACTG2",
  "TUBA1A", "TUBA1B", "TUBA1C", "TUBA3C", "TUBA3D", "TUBA3E",
  "TUBA4A", "TUBA8", "TUBB", "TUBB1", "TUBB2A", "TUBB2B", "TUBB3",
  "TUBB4A", "TUBB4B", "TUBB6", "TUBB8",
  "VIM", "DES", "GFAP", "LMNA", "LMNB1", "LMNB2",
  "MYH9", "MYH10", "MYH14", "MYL6", "MYL12A", "MYL12B",
  # 翻译/核糖体 (极高丰度)
  "EEF1A1", "EEF1A2", "EEF2",
  # 核仁/RNA 结合蛋白 (常见非特异性)
  "NCL", "NPM1", "HNRNPA1", "HNRNPA2B1", "HNRNPC", "HNRNPD",
  "HNRNPK", "HNRNPM", "HNRNPU", "PABPC1",
  # 线粒体 (高丰度背景)
  "VDAC1", "VDAC2", "VDAC3", "PHB", "PHB2",
  # 组蛋白 (核提取物常见)
  "H2AFZ", "H3F3A", "H3F3B", "HIST1H1C", "HIST1H2AB", "HIST1H4A",
  "H2AX", "H2BC1", "H3C1", "H4C1",
  # 其他高频 CRAPome 蛋白
  "FLNA", "FLNB", "PLEC", "DSP", "JUP", "PKM",
  "ENO1", "GAPDH", "LDHA", "LDHB", "PGK1", "TPI1", "ALDOA",
  "FASN", "PDIA3", "PDIA4", "PDIA6", "P4HB", "CALR", "CANX",
  "ANXA1", "ANXA2", "ANXA5", "ANXA6"
)


# ==============================================================================
# Bait 归一化
# ==============================================================================

#' AP-MS Bait 归一化
#'
#' 以 bait 蛋白的丰度为基准，对所有样本进行归一化。
#' 确保 bait 蛋白在不同条件/时间点的丰度一致，
#' 从而使 prey 蛋白的变化反映真实的结合动态。
#'
#' @param ms_data MsDataSet 对象
#' @param bait_gene bait 蛋白的基因名 (如 "UBASH3B", "CBL")
#' @return 归一化后的 MsDataSet 对象 (proteins 矩阵已更新)
#' @export
normalize_bait <- function(ms_data, bait_gene) {
  stopifnot(inherits(ms_data, "MsDataSet"))

  proteins <- ms_data$proteins
  pinfo    <- ms_data$protein_info

  # 查找 bait 蛋白
  gene_col <- if ("Gene" %in% colnames(pinfo)) pinfo$Gene else pinfo$Label_Name
  # 支持多基因名 (用分号分隔的情况, 取第一个)
  first_gene <- stringr::str_extract(gene_col, "^[^;]+")
  bait_idx <- which(first_gene == bait_gene)

  if (length(bait_idx) == 0) {
    stop(sprintf("\u274c \u672a\u627e\u5230 bait \u86cb\u767d: '%s'\n\u53ef\u7528\u57fa\u56e0\u540d\u793a\u4f8b: %s",
                 bait_gene, paste(head(first_gene[!is.na(first_gene)], 10), collapse = ", ")))
  }
  if (length(bait_idx) > 1) {
    message(sprintf("  \u26a0\ufe0f \u627e\u5230 %d \u4e2a\u5339\u914d bait \u86cb\u767d\uff0c\u4f7f\u7528\u7b2c\u4e00\u4e2a", length(bait_idx)))
    bait_idx <- bait_idx[1]
  }

  bait_values <- as.numeric(proteins[bait_idx, ])
  if (any(bait_values == 0 | is.na(bait_values))) {
    n_zero <- sum(bait_values == 0 | is.na(bait_values))
    warning(sprintf("Bait \u86cb\u767d\u5728 %d \u4e2a\u6837\u672c\u4e2d\u7f3a\u5931/\u4e3a0\uff0c\u5c06\u7528\u6700\u5c0f\u503c\u586b\u5145", n_zero))
    min_val <- min(bait_values[bait_values > 0 & !is.na(bait_values)], na.rm = TRUE)
    bait_values[bait_values == 0 | is.na(bait_values)] <- min_val * 0.1
  }

  # 计算归一化系数
  bait_mean <- mean(bait_values)
  norm_factors <- bait_mean / bait_values

  # 对每列乘以归一化系数
  proteins_norm <- as.data.frame(
    mapply(function(col, factor) col * factor,
           proteins, norm_factors, SIMPLIFY = FALSE)
  )
  colnames(proteins_norm) <- colnames(proteins)

  ms_data$proteins <- proteins_norm
  message(sprintf("  \u2705 Bait \u5f52\u4e00\u5316\u5b8c\u6210 (%s, mean=%.1f)", bait_gene, bait_mean))
  ms_data
}


# ==============================================================================
# 相对丰度计算
# ==============================================================================

#' 计算蛋白相对丰度
#'
#' 对每个蛋白计算各组的均值、标准差，以及相对丰度 (max=1 归一化)。
#' 返回 long-format data.frame，适合直接用于可视化。
#'
#' @param ms_data MsDataSet 对象
#' @param group_info 分组信息 (来自 set_groups/interactive_grouping)
#' @param min_valid 每组至少需要的有效值数目 (默认 2)
#' @return data.frame: Gene, Protein, Group, Mean, SD, RelAbundance, MaxGroup
#' @export
calc_relative_abundance <- function(ms_data, group_info, min_valid = 2) {
  stopifnot(inherits(ms_data, "MsDataSet"))

  proteins <- ms_data$proteins
  pinfo    <- ms_data$protein_info
  groups   <- unique(group_info$user_group)
  n_groups <- length(groups)

  # 样本-组对应
  sample_group <- setNames(group_info$user_group,
                           group_info$sample_name)

  # 计算每组均值和 SD
  mat_mean <- matrix(NA, nrow = nrow(proteins), ncol = n_groups)
  mat_sd   <- matrix(NA, nrow = nrow(proteins), ncol = n_groups)
  colnames(mat_mean) <- colnames(mat_sd) <- groups

  for (j in seq_along(groups)) {
    grp <- groups[j]
    grp_samples <- group_info$sample_name[group_info$user_group == grp]
    grp_samples <- intersect(grp_samples, colnames(proteins))

    if (length(grp_samples) == 0) next

    grp_data <- proteins[, grp_samples, drop = FALSE]

    for (i in seq_len(nrow(grp_data))) {
      vals <- as.numeric(grp_data[i, ])
      valid_vals <- vals[!is.na(vals) & vals > 0]

      if (length(valid_vals) >= min_valid) {
        mat_mean[i, j] <- mean(valid_vals)
        mat_sd[i, j]   <- if (length(valid_vals) > 1) sd(valid_vals) else 0
      } else {
        mat_mean[i, j] <- 0
        mat_sd[i, j]   <- 0
      }
    }
  }

  # 相对丰度 (每行的最大均值 = 1)
  row_max <- apply(mat_mean, 1, max, na.rm = TRUE)
  row_max[row_max == 0] <- 1  # 避免除 0
  mat_rel <- mat_mean / row_max

  # 最大丰度的组
  max_group <- apply(mat_mean, 1, function(x) which.max(x)[1])

  # 构建 long-format 结果
  gene_col <- if ("Gene" %in% colnames(pinfo)) {
    stringr::str_extract(pinfo$Gene, "^[^;]+")
  } else {
    pinfo$Label_Name
  }
  protein_col <- if ("Protein" %in% colnames(pinfo)) pinfo$Protein else gene_col

  result_list <- vector("list", n_groups)
  for (j in seq_along(groups)) {
    result_list[[j]] <- data.frame(
      Protein      = protein_col,
      Gene         = gene_col,
      Label_Name   = pinfo$Label_Name,
      Group        = groups[j],
      Mean         = mat_mean[, j],
      SD           = mat_sd[, j],
      RelAbundance = mat_rel[, j],
      MaxGroup     = max_group,
      stringsAsFactors = FALSE
    )
  }

  result <- do.call(rbind, result_list)
  # MaxGroup 是整数索引, 转为组名
  result$MaxGroup <- groups[result$MaxGroup]
  result
}


# ==============================================================================
# ANOVA 时序检验
# ==============================================================================

#' AP-MS 多组 ANOVA 检验
#'
#' 对每个蛋白做单因素 ANOVA 检验 (组间差异)。
#' 使用 Perseus 风格的缺失值填补 (downshift + narrowed distribution)。
#'
#' @param ms_data MsDataSet 对象
#' @param group_info 分组信息
#' @param p_cutoff p 值阈值 (默认 0.05)
#' @param p_adjust p 值校正方法 (默认 "BH")
#' @param impute 是否做缺失值填补 (默认 TRUE)
#' @param shift Perseus downshift 参数 (默认 1.8)
#' @param width Perseus width 参数 (默认 0.3)
#' @return data.frame 含 Protein, Gene, p.value, p.adj, q.significant 等列
#' @export
run_anova_timecourse <- function(ms_data, group_info,
                                  p_cutoff = 0.05,
                                  p_adjust = "BH",
                                  impute = TRUE,
                                  shift = 1.8, width = 0.3) {
  stopifnot(inherits(ms_data, "MsDataSet"))

  proteins <- ms_data$proteins
  pinfo    <- ms_data$protein_info

  # 匹配样本
  valid_samples <- intersect(colnames(proteins), group_info$sample_name)
  grp_vec <- group_info$user_group[match(valid_samples, group_info$sample_name)]
  expr_sub <- proteins[, valid_samples, drop = FALSE]

  # Log2 转换
  expr_sub[expr_sub == 0] <- NA
  if (isTRUE(ms_data$is_log2)) {
    expr_log2 <- as.matrix(expr_sub)
  } else {
    expr_log2 <- log2(as.matrix(expr_sub))
  }

  # 缺失值填补 (Perseus 风格)
  if (impute) {
    expr_log2 <- .perseus_impute(expr_log2, shift = shift, width = width)
  }

  # 对每行做 ANOVA
  p_values <- vapply(seq_len(nrow(expr_log2)), function(i) {
    df_row <- data.frame(
      value = expr_log2[i, ],
      group = grp_vec,
      stringsAsFactors = FALSE
    )
    # 至少需要 2 组有值
    if (length(unique(df_row$group[!is.na(df_row$value)])) < 2) return(NA_real_)
    tryCatch({
      fit <- aov(value ~ group, data = df_row)
      summary(fit)[[1]]$`Pr(>F)`[1]
    }, error = function(e) NA_real_)
  }, numeric(1))

  # 多重检验校正
  p_adj <- p.adjust(p_values, method = p_adjust)

  # 构建结果
  gene_col <- if ("Gene" %in% colnames(pinfo)) {
    stringr::str_extract(pinfo$Gene, "^[^;]+")
  } else {
    pinfo$Label_Name
  }
  protein_col <- if ("Protein" %in% colnames(pinfo)) pinfo$Protein else gene_col

  result <- data.frame(
    Protein    = protein_col,
    Gene       = gene_col,
    Label_Name = pinfo$Label_Name,
    p.value    = p_values,
    p.adj      = p_adj,
    significant = ifelse(!is.na(p_adj) & p_adj <= p_cutoff, "YES", "NO"),
    stringsAsFactors = FALSE
  )

  n_sig <- sum(result$significant == "YES", na.rm = TRUE)
  message(sprintf("  \u2705 ANOVA \u68c0\u9a8c\u5b8c\u6210: %d/%d \u86cb\u767d\u663e\u8457 (p.adj \u2264 %.2f)",
                  n_sig, nrow(result), p_cutoff))
  result
}


# ==============================================================================
# CRAPome 过滤
# ==============================================================================

#' CRAPome/cRAP 污染蛋白过滤
#'
#' 标记或过滤 AP-MS 实验中常见的非特异性结合蛋白。
#' 使用内嵌的 CRAPome + cRAP 蛋白列表 (角蛋白、血清蛋白、分子伴侣等)。
#'
#' @param ms_data MsDataSet 对象
#' @param action "flag" (默认) 或 "remove"。"flag" 添加 is_contaminant 列，
#'   "remove" 直接从数据中删除
#' @param extra_genes 额外的污染基因名向量 (添加到内置列表)
#' @return MsDataSet 对象 (protein_info 多一列 is_contaminant) 或过滤后的对象
#' @export
filter_crapome <- function(ms_data, action = "flag", extra_genes = NULL) {
  stopifnot(inherits(ms_data, "MsDataSet"))
  stopifnot(action %in% c("flag", "remove"))

  pinfo <- ms_data$protein_info
  contam_genes <- .crapome_genes
  if (!is.null(extra_genes)) {
    contam_genes <- unique(c(contam_genes, extra_genes))
  }

  # 匹配: Gene 列的第一个基因名
  gene_col <- if ("Gene" %in% colnames(pinfo)) pinfo$Gene else pinfo$Label_Name
  first_gene <- toupper(stringr::str_extract(gene_col, "^[^;]+"))
  contam_upper <- toupper(contam_genes)

  is_contam <- first_gene %in% contam_upper
  n_contam <- sum(is_contam, na.rm = TRUE)

  if (action == "flag") {
    ms_data$protein_info$is_contaminant <- is_contam
    message(sprintf("  \u2139\ufe0f \u6807\u8bb0 %d \u4e2a CRAPome \u6c61\u67d3\u86cb\u767d (%d \u4e2a\u5728\u5185\u7f6e\u5217\u8868\u4e2d)",
                    n_contam, length(contam_upper)))
  } else {
    if (n_contam > 0) {
      ms_data$proteins     <- ms_data$proteins[!is_contam, , drop = FALSE]
      ms_data$protein_info <- ms_data$protein_info[!is_contam, , drop = FALSE]
      message(sprintf("  \u2139\ufe0f \u79fb\u9664 %d \u4e2a CRAPome \u6c61\u67d3\u86cb\u767d\uff0c\u5269\u4f59 %d \u4e2a",
                      n_contam, nrow(ms_data$proteins)))
    } else {
      message("  \u2139\ufe0f \u672a\u53d1\u73b0 CRAPome \u6c61\u67d3\u86cb\u767d")
    }
  }

  ms_data
}


# ==============================================================================
# Mfuzz 时序聚类
# ==============================================================================

#' Mfuzz 软聚类
#'
#' 对显著蛋白做 Mfuzz 软聚类，按时序模式分组。
#' 需要安装 Bioconductor 包 Mfuzz。
#'
#' @param ms_data MsDataSet 对象
#' @param group_info 分组信息
#' @param sig_proteins 显著蛋白 ID 向量 (来自 ANOVA 结果)
#' @param n_clusters 聚类数目 (默认 4)
#' @param min_valid 每组最少有效值数 (默认 2)
#' @return list(eset=ExpressionSet, cl=Mfuzz聚类结果, membership=成员矩阵)
#' @export
run_mfuzz_cluster <- function(ms_data, group_info, sig_proteins,
                               n_clusters = 4, min_valid = 2) {
  if (!requireNamespace("Mfuzz", quietly = TRUE)) {
    stop("\u274c \u9700\u8981\u5b89\u88c5 Mfuzz \u5305\u3002\u8bf7\u8fd0\u884c:\n",
         "  BiocManager::install('Mfuzz')")
  }

  proteins <- ms_data$proteins
  pinfo    <- ms_data$protein_info
  groups   <- unique(group_info$user_group)
  n_groups <- length(groups)

  # 匹配蛋白 ID
  protein_col <- if ("Protein" %in% colnames(pinfo)) pinfo$Protein else pinfo$Label_Name
  keep_idx <- which(protein_col %in% sig_proteins)

  if (length(keep_idx) < n_clusters) {
    stop(sprintf("\u274c \u663e\u8457\u86cb\u767d\u6570 (%d) \u5c0f\u4e8e\u805a\u7c7b\u6570 (%d)",
                 length(keep_idx), n_clusters))
  }

  # 计算组均值矩阵
  mat_mean <- matrix(NA, nrow = length(keep_idx), ncol = n_groups)
  colnames(mat_mean) <- groups

  for (j in seq_along(groups)) {
    grp <- groups[j]
    grp_samples <- group_info$sample_name[group_info$user_group == grp]
    grp_samples <- intersect(grp_samples, colnames(proteins))
    if (length(grp_samples) == 0) next

    grp_data <- proteins[keep_idx, grp_samples, drop = FALSE]
    for (i in seq_along(keep_idx)) {
      vals <- as.numeric(grp_data[i, ])
      valid_vals <- vals[!is.na(vals) & vals > 0]
      if (length(valid_vals) >= min_valid) {
        mat_mean[i, j] <- mean(valid_vals)
      }
    }
  }

  # Log2 转换 + 行标准化
  mat_mean[mat_mean == 0] <- NA
  if (!isTRUE(ms_data$is_log2)) {
    mat_mean <- log2(mat_mean)
  }

  # Z-score 行标准化
  mat_z <- t(apply(mat_mean, 1, function(x) {
    if (all(is.na(x))) return(x)
    (x - mean(x, na.rm = TRUE)) / max(sd(x, na.rm = TRUE), 1e-10)
  }))
  colnames(mat_z) <- groups

  # 行名
  gene_col <- if ("Gene" %in% colnames(pinfo)) {
    stringr::str_extract(pinfo$Gene[keep_idx], "^[^;]+")
  } else {
    pinfo$Label_Name[keep_idx]
  }
  gene_col[is.na(gene_col)] <- protein_col[keep_idx][is.na(gene_col)]
  rownames(mat_z) <- make.unique(gene_col)

  # 去掉含 NA 的行
  valid <- complete.cases(mat_z)
  mat_z <- mat_z[valid, , drop = FALSE]

  if (nrow(mat_z) < n_clusters) {
    stop(sprintf("\u274c \u5b8c\u6574\u6570\u636e\u884c\u6570 (%d) \u5c0f\u4e8e\u805a\u7c7b\u6570 (%d)",
                 nrow(mat_z), n_clusters))
  }

  message(sprintf("  \u2139\ufe0f Mfuzz \u805a\u7c7b: %d \u86cb\u767d, %d \u7c7b", nrow(mat_z), n_clusters))

  # 构建 ExpressionSet
  eset <- Biobase::ExpressionSet(assayData = mat_z)

  # 估计 fuzzifier
  m <- Mfuzz::mestimate(eset)

  # 聚类
  cl <- Mfuzz::mfuzz(eset, c = n_clusters, m = m)

  # 成员矩阵
  membership <- cl$membership
  rownames(membership) <- rownames(mat_z)

  list(
    eset       = eset,
    cl         = cl,
    membership = membership,
    n_clusters = n_clusters,
    fuzzifier  = m
  )
}


# ==============================================================================
# 内部辅助函数
# ==============================================================================

#' Perseus 风格缺失值填补
#'
#' 从正态分布采样 (均值下移 shift*SD, 宽度缩窄为 width*SD)
#' @keywords internal
.perseus_impute <- function(mat, shift = 1.8, width = 0.3) {
  # 全局统计量
  all_vals <- as.vector(mat)
  all_vals <- all_vals[!is.na(all_vals)]
  if (length(all_vals) == 0) return(mat)

  global_mean <- mean(all_vals)
  global_sd   <- sd(all_vals)

  impute_mean <- global_mean - shift * global_sd
  impute_sd   <- width * global_sd

  n_na <- sum(is.na(mat))
  if (n_na > 0) {
    set.seed(1234)
    mat[is.na(mat)] <- rnorm(n_na, mean = impute_mean, sd = impute_sd)
    message(sprintf("  \u2139\ufe0f Perseus \u586b\u5145 %d \u4e2a\u7f3a\u5931\u503c (shift=%.1f, width=%.1f)",
                    n_na, shift, width))
  }
  mat
}
