# ==============================================================================
# msProteomiX - MaxQuant parser
# ==============================================================================

# ==============================================================================

#' 解析 MaxQuant 搜库结果
#'
#' 读取 MaxQuant 输出的 proteinGroups.txt 文件，返回 MsDataSet 对象。
#' 自动过滤 reverse, contaminant, only-identified-by-site 条目。
#' 定量列优先级: LFQ intensity > iBAQ > Intensity
#'
#' @param path 文件路径或包含 MaxQuant txt/ 结果的目录
#' @return MsDataSet 对象
#' @export
#' @examples
#' \dontrun{
#' ms <- read_ms_data("path/to/proteinGroups.txt", engine = "maxquant")
#' ms <- read_ms_data("path/to/txt/")  # 自动检测
#' }
parse_maxquant <- function(path) {
  # 定位文件
  pg_file <- NULL

  if (dir.exists(path)) {
    pg_file <- .find_file(path, "proteinGroups\\.txt$")
    # 搜索 txt/ 子目录
    if (is.null(pg_file)) {
      txt_dir <- file.path(path, "txt")
      if (dir.exists(txt_dir)) pg_file <- .find_file(txt_dir, "proteinGroups\\.txt$")
    }
    # 递归搜索子目录
    if (is.null(pg_file)) {
      subdirs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
      for (sd in subdirs) {
        pg_file <- .find_file(sd, "proteinGroups\\.txt$")
        if (!is.null(pg_file)) break
      }
    }
  } else {
    pg_file <- path
  }

  if (is.null(pg_file) || !file.exists(pg_file)) {
    stop("\u274c \u672a\u627e\u5230 MaxQuant \u7684 proteinGroups.txt \u6587\u4ef6\u3002")
  }

  message(sprintf(">>> \u8bfb\u53d6 MaxQuant proteinGroups: %s", basename(pg_file)))

  df <- .read_omics_file(pg_file)
  raw_cols <- colnames(df)

  # --- 过滤 reverse, contaminant, only-identified-by-site ---
  keep <- rep(TRUE, nrow(df))
  if ("Reverse" %in% raw_cols) {
    keep <- keep & (is.na(df$Reverse) | df$Reverse != "+")
  }
  if ("Potential contaminant" %in% raw_cols) {
    keep <- keep & (is.na(df$`Potential contaminant`) | df$`Potential contaminant` != "+")
  }
  if ("Only identified by site" %in% raw_cols) {
    keep <- keep & (is.na(df$`Only identified by site`) | df$`Only identified by site` != "+")
  }
  # 额外: Protein IDs 前缀过滤
  if ("Protein IDs" %in% raw_cols) {
    keep <- keep & !grepl("^REV__|^CON__", df$`Protein IDs`)
  }
  n_removed <- sum(!keep)
  if (n_removed > 0) {
    message(sprintf("  \u2139\ufe0f \u8fc7\u6ee4 %d \u4e2a reverse/contaminant/only-by-site \u6761\u76ee", n_removed))
  }
  df <- df[keep, ]

  # --- 提取定量列 ---
  quant_info <- .extract_mq_quant_cols(raw_cols)

  if (length(quant_info$cols) == 0) {
    stop("\u274c \u65e0\u6cd5\u8bc6\u522b MaxQuant \u5b9a\u91cf\u5217\u3002\n",
         "\u652f\u6301\u7684\u5217\u7c7b\u578b: LFQ intensity, iBAQ, Intensity")
  }

  quant_cols <- quant_info$cols
  sample_names <- quant_info$samples

  # --- 构建蛋白定量矩阵 ---
  proteins <- as.data.frame(lapply(df[, quant_cols, drop = FALSE], function(x) {
    x <- suppressWarnings(as.numeric(x))
    x[x == 0] <- NA  # MaxQuant 用 0 表示缺失
    x
  }))
  colnames(proteins) <- sample_names
  rownames(proteins) <- NULL

  # --- 过滤全 NA 行 ---
  valid_rows <- rowSums(!is.na(proteins)) > 0
  if (any(!valid_rows)) {
    message(sprintf("  \u2139\ufe0f \u79fb\u9664 %d \u4e2a\u5168\u7f3a\u5931\u86cb\u767d", sum(!valid_rows)))
    proteins <- proteins[valid_rows, , drop = FALSE]
    df <- df[valid_rows, , drop = FALSE]
  }

  # --- 构建蛋白注释信息 ---
  info_col_map <- c(
    "Protein IDs"          = "Protein",
    "Majority protein IDs" = "Protein ID",
    "Gene names"           = "Gene",
    "Protein names"        = "Description",
    "Fasta headers"        = "Fasta Header",
    "Score"                = "Score",
    "Number of proteins"   = "Num Proteins"
  )
  avail_info <- intersect(names(info_col_map), raw_cols)
  protein_info <- df[, avail_info, drop = FALSE]
  colnames(protein_info) <- info_col_map[avail_info]

  # 确保 Protein 列存在
  if (!("Protein" %in% colnames(protein_info)) && "Protein ID" %in% colnames(protein_info)) {
    protein_info$Protein <- protein_info$`Protein ID`
  }
  if (!("Protein ID" %in% colnames(protein_info)) && "Protein" %in% colnames(protein_info)) {
    protein_info$`Protein ID` <- protein_info$Protein
  }

  # Gene 列去除 trailing 分号
  if ("Gene" %in% colnames(protein_info)) {
    protein_info$Gene <- stringr::str_remove(protein_info$Gene, ";+$")
  }

  # 构建标签名
  protein_info$Label_Name <- .build_mq_label_name(df, raw_cols)

  new_MsDataSet(
    proteins     = proteins,
    peptides     = data.frame(),
    psms         = data.frame(),
    ions         = data.frame(),
    sample_names = sample_names,
    protein_info = protein_info,
    engine       = "maxquant",
    quant_type   = quant_info$type
  )
}


# ==============================================================================

# ==============================================================================

#' 提取 MaxQuant 定量列名和样本名
#'
#' MaxQuant proteinGroups.txt 的定量列格式为:
#'   "LFQ intensity SampleName" (前缀固定，样本名在后)
#'   "iBAQ SampleName" 或 "Intensity SampleName"
#'
#' @keywords internal
.extract_mq_quant_cols <- function(raw_cols) {
  # 优先级 1: LFQ intensity (列名格式: "LFQ intensity <sample>")
  idx <- grep("^LFQ intensity ", raw_cols, ignore.case = TRUE)
  if (length(idx) > 0) {
    cols <- raw_cols[idx]
    samples <- trimws(stringr::str_remove(cols, "(?i)^LFQ\\s+intensity\\s+"))
    return(list(cols = cols, samples = samples, type = "LFQ"))
  }

  # 优先级 2: iBAQ (列名格式: "iBAQ <sample>")
  idx <- grep("^iBAQ ", raw_cols, ignore.case = TRUE)
  # 排除 "iBAQ peptides" 等非定量列
  idx <- idx[!grepl("peptides", raw_cols[idx], ignore.case = TRUE)]
  if (length(idx) > 0) {
    cols <- raw_cols[idx]
    samples <- trimws(stringr::str_remove(cols, "(?i)^iBAQ\\s+"))
    return(list(cols = cols, samples = samples, type = "iBAQ"))
  }

  # 优先级 3: Intensity (列名格式: "Intensity <sample>")
  idx <- grep("^Intensity ", raw_cols, ignore.case = TRUE)
  if (length(idx) > 0) {
    cols <- raw_cols[idx]
    samples <- trimws(stringr::str_remove(cols, "(?i)^Intensity\\s+"))
    return(list(cols = cols, samples = samples, type = "Intensity"))
  }

  list(cols = character(), samples = character(), type = "Unknown")
}

#' 构建 MaxQuant 标签名
#' @keywords internal
.build_mq_label_name <- function(df, raw_cols) {
  if ("Gene names" %in% raw_cols) {
    labels <- df$`Gene names`
    # 多基因用分号分隔，取第一个
    labels <- stringr::str_extract(labels, "^[^;]+")
  } else if ("Protein names" %in% raw_cols) {
    labels <- df$`Protein names`
    labels <- stringr::str_extract(labels, "^[^;]+")
  } else if ("Majority protein IDs" %in% raw_cols) {
    labels <- df$`Majority protein IDs`
  } else if ("Protein IDs" %in% raw_cols) {
    labels <- df$`Protein IDs`
  } else {
    labels <- as.character(seq_len(nrow(df)))
  }
  labels[is.na(labels)] <- "Unknown"
  labels
}


# ==============================================================================
# Proteome Discoverer 内部辅助函数

