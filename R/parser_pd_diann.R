# ==============================================================================
# msProteomiX - Proteome Discoverer & DIA-NN parsers
# ==============================================================================

parse_pd <- function(path) {
  # 定位文件
  pd_file <- NULL

  if (dir.exists(path)) {
    pd_file <- .find_file(path, "_Proteins\\.(txt|tsv)$")
    if (is.null(pd_file)) {
      subdirs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
      for (sd in subdirs) {
        pd_file <- .find_file(sd, "_Proteins\\.(txt|tsv)$")
        if (!is.null(pd_file)) break
      }
    }
  } else {
    pd_file <- path
  }

  if (is.null(pd_file) || !file.exists(pd_file)) {
    stop("\u274c \u672a\u627e\u5230 Proteome Discoverer \u7684 *_Proteins.txt \u6587\u4ef6\u3002")
  }

  message(sprintf(">>> \u8bfb\u53d6 Proteome Discoverer: %s", basename(pd_file)))

  df <- .read_omics_file(pd_file)
  raw_cols <- colnames(df)

  # --- 过滤 contaminant ---
  if ("Contaminant" %in% raw_cols) {
    keep <- is.na(df$Contaminant) | df$Contaminant != "TRUE"
    df <- df[keep, ]
  }
  # Accession 前缀过滤
  acc_col <- intersect(c("Accession", "Master.Protein.Accessions"), raw_cols)[1]
  if (!is.na(acc_col)) {
    keep <- !grepl("^REV__|^CON__|^CONT_", df[[acc_col]])
    df <- df[keep, ]
  }

  # --- 提取定量列 ---
  quant_info <- .extract_pd_quant_cols(raw_cols)

  if (length(quant_info$cols) == 0) {
    stop("\u274c \u65e0\u6cd5\u8bc6\u522b PD \u5b9a\u91cf\u5217\u3002\n",
         "\u652f\u6301\u7684\u5217\u7c7b\u578b: Abundance, Abundances.Normalized, Abundances.Scaled")
  }

  quant_cols <- quant_info$cols
  sample_names <- quant_info$samples

  # --- 构建蛋白定量矩阵 ---
  proteins <- as.data.frame(lapply(df[, quant_cols, drop = FALSE], function(x) {
    x <- suppressWarnings(as.numeric(x))
    x[x == 0] <- NA
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
  # PD 标准列名 和 R-Friendly 列名映射
  info_col_map <- c(
    "Accession"                    = "Protein ID",
    "Master.Protein.Accessions"    = "Protein ID",
    "Description"                  = "Description",
    "Protein.Descriptions"         = "Description",
    "Gene Symbol"                  = "Gene",
    "Gene.Symbol"                  = "Gene",
    "Genes"                        = "Gene",
    "Coverage [%]"                 = "Coverage",
    "Coverage.in.Percent"          = "Coverage",
    "# Peptides"                   = "Num Peptides",
    "Number.of.Peptides"           = "Num Peptides",
    "# Unique Peptides"            = "Num Unique Peptides",
    "Number.of.Unique.Peptides"    = "Num Unique Peptides",
    "MW [kDa]"                     = "MW",
    "Molecular.Weight.in.kDa"      = "MW"
  )
  avail_info <- intersect(names(info_col_map), raw_cols)
  protein_info <- df[, avail_info, drop = FALSE]
  colnames(protein_info) <- info_col_map[avail_info]

  # 确保 Protein 和 Protein ID 列存在
  if (!("Protein" %in% colnames(protein_info))) {
    if ("Protein ID" %in% colnames(protein_info)) {
      protein_info$Protein <- protein_info$`Protein ID`
    } else if ("Gene" %in% colnames(protein_info)) {
      protein_info$Protein <- protein_info$Gene
    }
  }
  if (!("Protein ID" %in% colnames(protein_info)) && "Protein" %in% colnames(protein_info)) {
    protein_info$`Protein ID` <- protein_info$Protein
  }

  # Gene 列去除 trailing 分号
  if ("Gene" %in% colnames(protein_info)) {
    protein_info$Gene <- stringr::str_remove(protein_info$Gene, ";+$")
  }

  # 构建标签名
  if ("Gene" %in% colnames(protein_info)) {
    labels <- stringr::str_extract(protein_info$Gene, "^[^;]+")
    labels[is.na(labels)] <- protein_info$`Protein ID`[is.na(labels)]
  } else if ("Protein ID" %in% colnames(protein_info)) {
    labels <- protein_info$`Protein ID`
  } else {
    labels <- as.character(seq_len(nrow(protein_info)))
  }
  labels[is.na(labels)] <- "Unknown"
  protein_info$Label_Name <- labels

  new_MsDataSet(
    proteins     = proteins,
    peptides     = data.frame(),
    psms         = data.frame(),
    ions         = data.frame(),
    sample_names = sample_names,
    protein_info = protein_info,
    engine       = "pd",
    quant_type   = quant_info$type
  )
}


# ==============================================================================
# DIA-NN 解析器
# ==============================================================================

#' 解析 DIA-NN 搜库结果
#'
#' 读取 DIA-NN 输出的 report.pg_matrix.tsv 文件 (宽格式蛋白定量矩阵)。
#' 元数据列 (Protein.Group, Genes 等) 自动识别，剩余列作为样本定量列。
#' 样本名从完整文件路径中提取干净名称。
#'
#' @param path 文件路径或包含 DIA-NN 报告的目录
#' @return MsDataSet 对象
#' @export
#' @examples
#' \dontrun{
#' ms <- read_ms_data("path/to/report.pg_matrix.tsv", engine = "diann")
#' }
parse_diann <- function(path) {
  # 定位文件
  pg_file <- NULL

  if (dir.exists(path)) {
    pg_file <- .find_file(path, "report\\.pg_matrix\\.(tsv|csv)$")
    if (is.null(pg_file)) {
      # 也搜索 .pg_matrix 后缀的变体
      pg_file <- .find_file(path, "pg_matrix\\.(tsv|csv)$")
    }
    if (is.null(pg_file)) {
      subdirs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
      for (sd in subdirs) {
        pg_file <- .find_file(sd, "report\\.pg_matrix\\.(tsv|csv)$")
        if (!is.null(pg_file)) break
        pg_file <- .find_file(sd, "pg_matrix\\.(tsv|csv)$")
        if (!is.null(pg_file)) break
      }
    }
  } else {
    pg_file <- path
  }

  if (is.null(pg_file) || !file.exists(pg_file)) {
    stop("\u274c \u672a\u627e\u5230 DIA-NN \u7684 report.pg_matrix.tsv \u6587\u4ef6\u3002")
  }

  message(sprintf(">>> \u8bfb\u53d6 DIA-NN pg_matrix: %s", basename(pg_file)))

  df <- .read_omics_file(pg_file)
  raw_cols <- colnames(df)

  # --- 识别元数据列 vs 定量列 ---
  # DIA-NN 元数据列固定名称
  diann_meta_cols <- c("Protein.Group", "Protein.Ids", "Protein.Names",
                       "Genes", "First.Protein.Description",
                       "N.Sequences", "N.Proteotypic.Sequences")
  is_meta <- raw_cols %in% diann_meta_cols
  quant_cols <- raw_cols[!is_meta]

  if (length(quant_cols) == 0) {
    stop("\u274c \u65e0\u6cd5\u8bc6\u522b DIA-NN \u5b9a\u91cf\u5217\u3002")
  }

  # --- 提取样本名 (从文件路径中提取干净名称) ---
  sample_names <- .extract_diann_sample_names(quant_cols)

  # --- 过滤 contaminant / reverse ---
  if ("Protein.Group" %in% raw_cols) {
    keep <- !grepl("^REV__|^CON__|^CONT_|^con_", df$Protein.Group)
    df <- df[keep, ]
  }

  # --- 构建蛋白定量矩阵 ---
  proteins <- as.data.frame(lapply(df[, quant_cols, drop = FALSE], function(x) {
    x <- suppressWarnings(as.numeric(x))
    x[x == 0] <- NA  # DIA-NN 用 0 表示缺失
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
    "Protein.Group"              = "Protein",
    "Protein.Ids"                = "Protein ID",
    "Genes"                      = "Gene",
    "First.Protein.Description"  = "Description",
    "Protein.Names"              = "Entry Name"
  )
  avail_info <- intersect(names(info_col_map), raw_cols)
  protein_info <- df[, avail_info, drop = FALSE]
  colnames(protein_info) <- info_col_map[avail_info]

  # 确保 Protein / Protein ID 列
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
  if ("Gene" %in% colnames(protein_info)) {
    labels <- stringr::str_extract(protein_info$Gene, "^[^;]+")
    labels[is.na(labels)] <- if ("Protein" %in% colnames(protein_info)) {
      protein_info$Protein[is.na(labels)]
    } else {
      "Unknown"
    }
  } else if ("Protein" %in% colnames(protein_info)) {
    labels <- protein_info$Protein
  } else {
    labels <- as.character(seq_len(nrow(protein_info)))
  }
  labels[is.na(labels)] <- "Unknown"
  protein_info$Label_Name <- labels

  new_MsDataSet(
    proteins     = proteins,
    peptides     = data.frame(),
    psms         = data.frame(),
    ions         = data.frame(),
    sample_names = sample_names,
    protein_info = protein_info,
    engine       = "diann",
    quant_type   = "MaxLFQ"
  )
}


# ==============================================================================
# 内部辅助函数
# ==============================================================================

#' 在目录中查找匹配文件
#' @keywords internal

.extract_pd_quant_cols <- function(raw_cols) {
  # 优先级 1: Abundances Scaled (R-Friendly: Abundances.Scaled.)
  idx <- grep("^Abundances?\\.?\\s*Scaled", raw_cols, ignore.case = TRUE)
  if (length(idx) > 1) {
    cols <- raw_cols[idx]
    samples <- .clean_pd_sample_names(cols, "Scaled")
    return(list(cols = cols, samples = samples, type = "Abundance.Scaled"))
  }

  # 优先级 2: Abundances Normalized
  idx <- grep("^Abundances?\\.?\\s*Normalized", raw_cols, ignore.case = TRUE)
  if (length(idx) > 1) {
    cols <- raw_cols[idx]
    samples <- .clean_pd_sample_names(cols, "Normalized")
    return(list(cols = cols, samples = samples, type = "Abundance.Normalized"))
  }

  # 优先级 3: Abundance (raw)
  # 标准格式: "Abundance: F1: Sample" 或 R-Friendly: "Abundances..F1..Sample"
  idx <- grep("^Abundance", raw_cols, ignore.case = TRUE)
  # 排除已经匹配的 Normalized/Scaled/Grouped/Ratio/Count
  idx <- idx[!grepl("Normalized|Scaled|Grouped|Ratio|Count|CV|Variability",
                     raw_cols[idx], ignore.case = TRUE)]
  if (length(idx) > 1) {
    cols <- raw_cols[idx]
    samples <- .clean_pd_sample_names(cols, "")
    return(list(cols = cols, samples = samples, type = "Abundance"))
  }

  list(cols = character(), samples = character(), type = "Unknown")
}

#' 清理 PD 样本名
#' @keywords internal
.clean_pd_sample_names <- function(cols, quant_type) {
  # 去掉前缀 "Abundance: " / "Abundances Normalized: " 等
  samples <- cols
  # R-Friendly: "Abundances.Normalized.F1.Sample"
  if (any(grepl("\\.", samples))) {
    samples <- stringr::str_remove(samples, "(?i)^Abundances?\\.?")
    if (nchar(quant_type) > 0) {
      samples <- stringr::str_remove(samples, paste0("(?i)\\.?", quant_type, "\\.?"))
    }
    samples <- stringr::str_remove(samples, "^\\.")
  } else {
    # 标准: "Abundance: F1: Sample" 或 "Abundance Normalized: F1: Sample"
    samples <- stringr::str_remove(samples, "(?i)^Abundances?\\s*")
    if (nchar(quant_type) > 0) {
      samples <- stringr::str_remove(samples, paste0("(?i)", quant_type, "\\s*:?\\s*"))
    }
    samples <- stringr::str_remove(samples, "^:\\s*")
  }
  # 去掉 F1:, F2: 等 fraction 编号
  samples <- stringr::str_remove(samples, "(?i)^F\\d+\\.?:?\\s*")
  samples <- trimws(samples)
  # 如果还是太长或有路径, 取 basename
  if (any(nchar(samples) > 60)) {
    samples <- basename(samples)
  }
  samples
}


# ==============================================================================
# DIA-NN 内部辅助函数
# ==============================================================================

#' 从 DIA-NN 列名中提取干净样本名
#'
#' DIA-NN 的列名通常是完整文件路径，如:
#'   "\\\\server\\path\\DIA_Sample1.d"
#' 需要提取干净的样本名。
#'
#' @keywords internal
.extract_diann_sample_names <- function(quant_cols) {
  # 标准化路径分隔符 (DIA-NN 的列名常含 Windows 反斜杠路径)
  normalized <- gsub("\\\\", "/", quant_cols)
  # 取 basename (最后一段路径)
  samples <- basename(normalized)
  # 去掉文件扩展名
  samples <- stringr::str_remove(samples, "\\.(d|raw|mzML|wiff|dia)$")
  # 如果所有样本名都一样 (说明 basename 不够区分), 用倒数第二段
  if (length(unique(samples)) < length(samples)) {
    parents <- basename(dirname(normalized))
    samples <- paste0(parents, "/", samples)
  }
  samples
}


