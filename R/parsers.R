# ==============================================================================
# msProteomiX — 搜库引擎解析器
# ==============================================================================

#' 自动检测并解析质谱搜库结果
#'
#' 根据文件或目录的特征自动检测搜库引擎类型，并将结果解析为统一的
#' MsDataSet 对象。支持 FragPipe、MaxQuant、Proteome Discoverer 和 DIA-NN。
#'
#' @param path 文件路径或包含搜库结果的目录
#' @param engine 引擎名称。"auto" 表示自动检测。
#'   可选值: "auto", "fragpipe", "maxquant", "pd", "diann"
#' @return MsDataSet 对象
#' @export
#' @examples
#' \dontrun{
#' # 自动检测
#' ms <- read_ms_data("path/to/results/")
#'
#' # 指定引擎
#' ms <- read_ms_data("path/to/combined_protein.tsv", engine = "fragpipe")
#' }
read_ms_data <- function(path, engine = "auto") {
  if (!file.exists(path) && !dir.exists(path)) {
    stop("\u274c \u9519\u8bef: \u8def\u5f84\u4e0d\u5b58\u5728: ", path)
  }

  if (engine == "auto") {
    engine <- detect_engine(path)
    message(sprintf("\u2705 \u81ea\u52a8\u68c0\u6d4b\u5230\u641c\u5e93\u5f15\u64ce: %s", engine))
  }

  result <- switch(engine,
    fragpipe = parse_fragpipe(path),
    maxquant = parse_maxquant(path),
    pd       = parse_pd(path),
    diann    = parse_diann(path),
    stop("\u274c \u4e0d\u652f\u6301\u7684\u641c\u5e93\u5f15\u64ce: ", engine,
         "\n\u652f\u6301\u7684\u5f15\u64ce: fragpipe, maxquant, pd, diann")
  )

  message(sprintf("\u2705 \u6570\u636e\u8bfb\u53d6\u5b8c\u6210: %d \u4e2a\u86cb\u767d, %d \u4e2a\u6837\u54c1",
                  nrow(result$proteins), length(result$sample_names)))
  result
}


#' 自动检测搜库引擎类型
#'
#' 通过检查目录或文件的特征文件名来判断使用的搜库引擎。
#'
#' @param path 文件路径或目录
#' @return 引擎名称字符串
#' @export
detect_engine <- function(path) {
  # 如果是文件，检查文件名

  if (file.exists(path) && !.is_directory(path)) {
    fname <- basename(path)
    if (grepl("combined_protein", fname, ignore.case = TRUE)) return("fragpipe")
    if (grepl("proteinGroups", fname, ignore.case = TRUE))    return("maxquant")
    if (grepl("_Proteins", fname, ignore.case = TRUE))        return("pd")
    if (grepl("report\\.pg_matrix", fname, ignore.case = TRUE)) return("diann")
  }

  # 如果是目录，扫描子文件
  if (dir.exists(path)) {
    files <- list.files(path, recursive = TRUE, full.names = FALSE)

    if (any(grepl("combined_protein", files, ignore.case = TRUE))) return("fragpipe")
    if (any(grepl("proteinGroups\\.txt", files, ignore.case = TRUE))) return("maxquant")
    if (any(grepl("_Proteins\\.txt", files, ignore.case = TRUE)))    return("pd")
    if (any(grepl("report\\.pg_matrix", files, ignore.case = TRUE))) return("diann")

    # 递归搜索子目录
    subdirs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
    for (subdir in subdirs) {
      subfiles <- list.files(subdir, recursive = TRUE, full.names = FALSE)
      if (any(grepl("combined_protein", subfiles, ignore.case = TRUE))) return("fragpipe")
      if (any(grepl("proteinGroups\\.txt", subfiles, ignore.case = TRUE))) return("maxquant")
    }
  }

  stop("\u274c \u65e0\u6cd5\u81ea\u52a8\u8bc6\u522b\u641c\u5e93\u5f15\u64ce\u3002\n",
       "\u8bf7\u901a\u8fc7 engine \u53c2\u6570\u624b\u52a8\u6307\u5b9a: ",
       "'fragpipe', 'maxquant', 'pd', 'diann'")
}


# ==============================================================================
# FragPipe 解析器
# ==============================================================================

#' 解析 FragPipe 搜库结果
#'
#' 读取 FragPipe 输出的 combined_protein.tsv (及可选的 combined_peptide.tsv,
#' combined_ion.tsv, psm.tsv) 文件，返回 MsDataSet 对象。
#'
#' @param path 文件路径或包含 FragPipe 结果的目录
#' @return MsDataSet 对象
#' @export
parse_fragpipe <- function(path) {
  # 定位文件
  if (dir.exists(path)) {
    data_dir <- path
    # 搜索子目录中的文件
    file_prot <- .find_file(data_dir, "combined_protein.*\\.(tsv|csv|txt)")
    if (is.null(file_prot)) {
      # 尝试搜索子目录
      subdirs <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE)
      for (sd in subdirs) {
        file_prot <- .find_file(sd, "combined_protein.*\\.(tsv|csv|txt)")
        if (!is.null(file_prot)) { data_dir <- sd; break }
      }
    }
  } else {
    file_prot <- path
    data_dir <- dirname(path)
  }

  if (is.null(file_prot)) {
    stop("\u274c \u672a\u627e\u5230 FragPipe \u7684 combined_protein \u6587\u4ef6\u3002")
  }

  message(sprintf(">>> \u8bfb\u53d6 Protein \u8868: %s", basename(file_prot)))

  # 读取 Protein 表
  prot_df <- .read_omics_file(file_prot)

  # 过滤反库和污染
  if ("Protein" %in% colnames(prot_df)) {
    prot_df <- prot_df[!grepl("^REV_|^CONT_", prot_df$Protein), ]
  }

  # 提取定量列和样本名
  raw_cols <- colnames(prot_df)
  quant_info <- .extract_fp_quant_cols(raw_cols)
  quant_cols <- quant_info$cols
  sample_names <- quant_info$samples

  if (length(quant_cols) == 0) {
    stop("\u274c \u65e0\u6cd5\u8bc6\u522b\u5b9a\u91cf\u5217\u3002")
  }

  # 构建蛋白定量矩阵
  proteins <- as.data.frame(lapply(prot_df[, quant_cols, drop = FALSE], as.numeric))
  colnames(proteins) <- sample_names
  rownames(proteins) <- NULL

  # 构建蛋白注释信息
  info_cols <- intersect(c("Protein", "Protein ID", "Entry Name", "Gene",
                            "Organism", "Protein Length", "Coverage",
                            "Protein Probability", "Top Peptide Probability"),
                          raw_cols)
  protein_info <- prot_df[, info_cols, drop = FALSE]

  # 构建标签名 (用于火山图标注)
  protein_info$Label_Name <- .build_label_name(prot_df)

  # 读取可选文件
  file_pep <- .find_file(data_dir, "combined_peptide.*\\.(tsv|csv|txt)")
  file_ion <- .find_file(data_dir, "combined_ion.*\\.(tsv|csv|txt)")

  peptides <- data.frame()
  ions <- data.frame()

  if (!is.null(file_pep)) {
    message(sprintf(">>> \u8bfb\u53d6 Peptide \u8868: %s", basename(file_pep)))
    peptides <- .read_omics_file(file_pep)
  }
  if (!is.null(file_ion)) {
    message(sprintf(">>> \u8bfb\u53d6 Ion \u8868: %s", basename(file_ion)))
    ions <- .read_omics_file(file_ion)
  }

  # 读取 PSM 文件 (如果有)
  psms <- data.frame()
  psm_files <- list.files(data_dir, pattern = "^psm\\.tsv$", full.names = TRUE,
                           recursive = TRUE)
  if (length(psm_files) > 0) {
    message(sprintf(">>> \u8bfb\u53d6 PSM \u8868: %d \u4e2a\u6587\u4ef6", length(psm_files)))
    psm_list <- lapply(psm_files, .read_omics_file)
    psms <- do.call(rbind, psm_list)
  }

  new_MsDataSet(
    proteins     = proteins,
    peptides     = peptides,
    psms         = psms,
    ions         = ions,
    sample_names = sample_names,
    protein_info = protein_info,
    engine       = "fragpipe",
    quant_type   = quant_info$type
  )
}


# ==============================================================================
# MaxQuant 解析器 (Phase 2 占位)
# ==============================================================================

#' 解析 MaxQuant 搜库结果
#'
#' @param path 文件路径或目录
#' @return MsDataSet 对象
#' @export
parse_maxquant <- function(path) {
  stop("\u274c MaxQuant \u89e3\u6790\u5668\u5c06\u5728 Phase 2 \u5b9e\u73b0\u3002\n",
       "\u5f53\u524d\u7248\u672c\u652f\u6301: FragPipe")
}


# ==============================================================================
# Proteome Discoverer 解析器 (Phase 2 占位)
# ==============================================================================

#' 解析 Proteome Discoverer 搜库结果
#'
#' @param path 文件路径或目录
#' @return MsDataSet 对象
#' @export
parse_pd <- function(path) {
  stop("\u274c Proteome Discoverer \u89e3\u6790\u5668\u5c06\u5728 Phase 2 \u5b9e\u73b0\u3002\n",
       "\u5f53\u524d\u7248\u672c\u652f\u6301: FragPipe")
}


# ==============================================================================
# DIA-NN 解析器 (Phase 2 占位)
# ==============================================================================

#' 解析 DIA-NN 搜库结果
#'
#' @param path 文件路径或目录
#' @return MsDataSet 对象
#' @export
parse_diann <- function(path) {
  stop("\u274c DIA-NN \u89e3\u6790\u5668\u5c06\u5728 Phase 2 \u5b9e\u73b0\u3002\n",
       "\u5f53\u524d\u7248\u672c\u652f\u6301: FragPipe")
}


# ==============================================================================
# 内部辅助函数
# ==============================================================================

#' 在目录中查找匹配文件
#' @keywords internal
.find_file <- function(dir, pattern) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(f) > 0) f[1] else NULL
}

#' 读取 TSV/CSV 文件
#' @keywords internal
.read_omics_file <- function(fpath) {
  if (grepl("\\.csv$", fpath, ignore.case = TRUE)) {
    readr::read_csv(fpath, show_col_types = FALSE)
  } else {
    readr::read_delim(fpath, delim = "\t", escape_double = FALSE,
                      trim_ws = TRUE, show_col_types = FALSE)
  }
}

#' 提取 FragPipe 定量列名和样本名
#'
#' FragPipe combined_protein.tsv 的定量列格式为:
#'   "SampleName MaxLFQ Intensity" (样本名在前，后缀固定)
#'
#' @keywords internal
.extract_fp_quant_cols <- function(raw_cols) {
  # 优先级 1: MaxLFQ Intensity (列名格式: "<sample> MaxLFQ Intensity")
  idx <- grep("MaxLFQ Intensity$", raw_cols, ignore.case = TRUE)
  if (length(idx) > 0) {
    cols <- raw_cols[idx]
    samples <- trimws(stringr::str_remove(cols, "(?i)\\s*MaxLFQ\\s*Intensity$"))
    return(list(cols = cols, samples = samples, type = "LFQ"))
  }

  # 优先级 2: Intensity (非 MaxLFQ, 列名格式: "<sample> Intensity")
  idx <- grep("Intensity$", raw_cols, ignore.case = TRUE)
  idx <- idx[!grepl("MaxLFQ", raw_cols[idx], ignore.case = TRUE)]
  # 排除 "Total Intensity" 等汇总列
  idx <- idx[!grepl("^(Total|Combined)", raw_cols[idx], ignore.case = TRUE)]
  if (length(idx) > 0) {
    cols <- raw_cols[idx]
    samples <- trimws(stringr::str_remove(cols, "(?i)\\s*Intensity$"))
    return(list(cols = cols, samples = samples, type = "Intensity"))
  }

  # 优先级 3: Spectral Count (列名格式: "<sample> Spectral Count")
  idx <- grep("Spectral Count$", raw_cols, ignore.case = TRUE)
  idx <- idx[!grepl("(Unique|Total)", raw_cols[idx], ignore.case = TRUE)]
  if (length(idx) > 0) {
    cols <- raw_cols[idx]
    samples <- trimws(stringr::str_remove(cols, "(?i)\\s*Spectral Count$"))
    return(list(cols = cols, samples = samples, type = "SpectralCount"))
  }

  list(cols = character(), samples = character(), type = "Unknown")
}

#' 构建标签名 (用于火山图标注)
#' @keywords internal
.build_label_name <- function(df) {
  if ("Entry Name" %in% colnames(df)) {
    labels <- stringr::str_match(df$`Entry Name`, "\\|([^|]+)$")[, 2]
    labels[is.na(labels)] <- df$`Entry Name`[is.na(labels)]
  } else if ("Gene" %in% colnames(df)) {
    labels <- df$Gene
  } else if ("Protein" %in% colnames(df)) {
    labels <- df$Protein
  } else {
    labels <- seq_len(nrow(df))
  }
  labels <- stringr::str_remove(as.character(labels), "^(sp|tr)\\|")
  labels
}

#' 检查路径是否为目录
#' @keywords internal
.is_directory <- function(path) {
  file.info(path)$isdir
}
