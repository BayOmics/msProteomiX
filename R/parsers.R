# ==============================================================================
# msProteomiX — 搜库引擎解析器
# ==============================================================================

#' 自动检测并解析质谱搜库结果
#'
#' 根据文件或目录的特征自动检测搜库引擎类型，并将结果解析为统一的
#' MsDataSet 对象。支持 FragPipe、Spectronaut、MaxQuant、Proteome Discoverer 和 DIA-NN。
#'
#' @param path 文件路径或包含搜库结果的目录
#' @param engine 引擎名称。"auto" 表示自动检测。
#'   可选值: "auto", "fragpipe", "spectronaut", "maxquant", "pd", "diann"
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
    fragpipe    = parse_fragpipe(path),
    spectronaut = parse_spectronaut(path),
    maxquant    = parse_maxquant(path),
    pd          = parse_pd(path),
    diann       = parse_diann(path),
    stop("\u274c \u4e0d\u652f\u6301\u7684\u641c\u5e93\u5f15\u64ce: ", engine,
         "\n\u652f\u6301\u7684\u5f15\u64ce: fragpipe, spectronaut, maxquant, pd, diann")
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
  # 如果是文件，检查文件名和内容特征
  if (file.exists(path) && !.is_directory(path)) {
    fname <- basename(path)
    if (grepl("combined_protein", fname, ignore.case = TRUE)) return("fragpipe")
    if (grepl("proteinGroups", fname, ignore.case = TRUE))    return("maxquant")
    if (grepl("_Proteins", fname, ignore.case = TRUE))        return("pd")
    if (grepl("report\\.pg_matrix", fname, ignore.case = TRUE)) return("diann")

    # Spectronaut: 文件名含 _Report.tsv 或内容含 PG. 前缀列
    if (grepl("_Report\\.(tsv|csv)$", fname, ignore.case = TRUE) ||
        grepl("_ME_Report", fname, ignore.case = TRUE)) {
      if (.detect_spectronaut_header(path)) return("spectronaut")
    }
    # 兜底: 尝试读取表头检测 Spectronaut
    if (grepl("\\.(tsv|csv|txt)$", fname, ignore.case = TRUE)) {
      if (.detect_spectronaut_header(path)) return("spectronaut")
    }
  }

  # 如果是目录，扫描子文件
  if (dir.exists(path)) {
    files <- list.files(path, recursive = TRUE, full.names = FALSE)

    if (any(grepl("combined_protein", files, ignore.case = TRUE))) return("fragpipe")
    if (any(grepl("proteinGroups\\.txt", files, ignore.case = TRUE))) return("maxquant")
    if (any(grepl("_Proteins\\.txt", files, ignore.case = TRUE)))    return("pd")
    if (any(grepl("report\\.pg_matrix", files, ignore.case = TRUE))) return("diann")

    # Spectronaut: 目录中含 _Report.tsv 或 PG. 前缀的 TSV
    sn_candidates <- files[grepl("_Report\\.(tsv|csv)$", files, ignore.case = TRUE)]
    if (length(sn_candidates) > 0) {
      full_path <- file.path(path, sn_candidates[1])
      if (.detect_spectronaut_header(full_path)) return("spectronaut")
    }

    # 递归搜索子目录
    subdirs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
    for (subdir in subdirs) {
      subfiles <- list.files(subdir, recursive = TRUE, full.names = FALSE)
      if (any(grepl("combined_protein", subfiles, ignore.case = TRUE))) return("fragpipe")
      if (any(grepl("proteinGroups\\.txt", subfiles, ignore.case = TRUE))) return("maxquant")
      # Spectronaut in subdir
      sn_sub <- subfiles[grepl("_Report\\.(tsv|csv)$", subfiles, ignore.case = TRUE)]
      if (length(sn_sub) > 0) {
        full_path <- file.path(subdir, sn_sub[1])
        if (.detect_spectronaut_header(full_path)) return("spectronaut")
      }
    }
  }

  stop("\u274c \u65e0\u6cd5\u81ea\u52a8\u8bc6\u522b\u641c\u5e93\u5f15\u64ce\u3002\n",
       "\u8bf7\u901a\u8fc7 engine \u53c2\u6570\u624b\u52a8\u6307\u5b9a: ",
       "'fragpipe', 'spectronaut', 'maxquant', 'pd', 'diann'")
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
# Spectronaut 解析器
# ==============================================================================

#' 解析 Spectronaut 搜库结果
#'
#' 支持 Spectronaut Run Pivot (宽格式) 导出，自动识别两种常见格式：
#' \itemize{
#'   \item 格式A: 列名为干净样本名 (如 \code{PlasmaX-1}, \code{PlasmaX-2})
#'   \item 格式B: 列名含索引和后缀 (如 \code{[1] filename.d.PG.Quantity})
#' }
#'
#' @param path 文件路径或包含 Spectronaut 导出文件的目录
#' @return MsDataSet 对象
#' @export
#' @examples
#' \dontrun{
#' ms <- read_ms_data("path/to/report.tsv", engine = "spectronaut")
#' ms <- read_ms_data("path/to/wkdir/")  # 自动检测
#' }
parse_spectronaut <- function(path) {
  # 定位报告文件
  report_file <- NULL

  if (dir.exists(path)) {
    # 在目录中搜索 Spectronaut 报告文件
    candidates <- list.files(path, pattern = "_Report\\.(tsv|csv)$",
                              full.names = TRUE, ignore.case = TRUE)
    if (length(candidates) == 0) {
      # 尝试搜索所有 tsv/csv 文件并通过表头检测
      all_files <- list.files(path, pattern = "\\.(tsv|csv)$",
                               full.names = TRUE, ignore.case = TRUE)
      # 排除 IdentificationsOverview 和 group_info 等辅助文件
      all_files <- all_files[!grepl("(group_info|IdentificationsOverview)",
                                      basename(all_files), ignore.case = TRUE)]
      for (f in all_files) {
        if (.detect_spectronaut_header(f)) { report_file <- f; break }
      }
    } else {
      report_file <- candidates[1]
    }

    # 搜索子目录
    if (is.null(report_file)) {
      subdirs <- list.dirs(path, recursive = FALSE, full.names = TRUE)
      for (sd in subdirs) {
        sub_candidates <- list.files(sd, pattern = "_Report\\.(tsv|csv)$",
                                      full.names = TRUE, ignore.case = TRUE)
        if (length(sub_candidates) > 0) { report_file <- sub_candidates[1]; break }
      }
    }
  } else {
    report_file <- path
  }

  if (is.null(report_file) || !file.exists(report_file)) {
    stop("\u274c \u672a\u627e\u5230 Spectronaut \u62a5\u544a\u6587\u4ef6\u3002\n",
         "\u8bf7\u786e\u4fdd\u76ee\u5f55\u4e2d\u5305\u542b *_Report.tsv \u6216 *_Report.csv \u6587\u4ef6\u3002")
  }

  message(sprintf(">>> \u8bfb\u53d6 Spectronaut \u62a5\u544a: %s", basename(report_file)))

  # 读取数据
  df <- .read_omics_file(report_file)
  raw_cols <- colnames(df)

  # --- 识别 PG. 元数据列 ---
  pg_meta_patterns <- c("^PG\\.", "^PG\\.ProteinGroups$", "^PG\\.ProteinAccessions$",
                         "^PG\\.Genes$", "^PG\\.ProteinDescriptions$",
                         "^PG\\.ProteinNames$", "^PG\\.UniProtIds$",
                         "^PG\\.FastaFiles$", "^PG\\.Qvalue$",
                         "^PG\\.MolecularWeight$")
  is_pg_col <- grepl("^PG\\.", raw_cols)

  # --- 提取定量列和样本名 ---
  quant_info <- .extract_sn_quant_cols(raw_cols, is_pg_col)
  quant_cols <- quant_info$cols
  sample_names <- quant_info$samples

  if (length(quant_cols) == 0) {
    stop("\u274c \u65e0\u6cd5\u8bc6\u522b Spectronaut \u5b9a\u91cf\u5217\u3002\n",
         "\u8bf7\u786e\u4fdd\u5bfc\u51fa\u65f6\u9009\u62e9\u4e86 Run Pivot \u683c\u5f0f\u3002")
  }

  # --- 过滤反库和污染蛋白 ---
  if ("PG.ProteinGroups" %in% raw_cols) {
    keep <- !grepl("^REV_|^CONT_|^CON__|^con_", df$PG.ProteinGroups)
    df <- df[keep, ]
  } else if ("PG.ProteinAccessions" %in% raw_cols) {
    keep <- !grepl("^REV_|^CONT_|^CON__|^con_", df$PG.ProteinAccessions)
    df <- df[keep, ]
  }

  # --- 构建蛋白定量矩阵 ---
  proteins <- as.data.frame(lapply(df[, quant_cols, drop = FALSE], function(x) {
    # Spectronaut 用 NaN 表示缺失值
    x <- as.numeric(x)
    x[is.nan(x)] <- NA
    x
  }))
  colnames(proteins) <- sample_names
  rownames(proteins) <- NULL

  # --- 构建蛋白注释信息 ---
  info_col_map <- c(
    "PG.ProteinGroups"       = "Protein",
    "PG.ProteinAccessions"   = "Protein ID",
    "PG.Genes"               = "Gene",
    "PG.ProteinDescriptions" = "Description",
    "PG.ProteinNames"        = "Entry Name",
    "PG.UniProtIds"          = "UniProt ID",
    "PG.Qvalue"              = "Q-value",
    "PG.MolecularWeight"     = "MW"
  )
  avail_info <- intersect(names(info_col_map), raw_cols)
  protein_info <- df[, avail_info, drop = FALSE]
  colnames(protein_info) <- info_col_map[avail_info]

  # 构建标签名 (用于火山图标注)
  protein_info$Label_Name <- .build_sn_label_name(df, raw_cols)

  new_MsDataSet(
    proteins     = proteins,
    peptides     = data.frame(),
    psms         = data.frame(),
    ions         = data.frame(),
    sample_names = sample_names,
    protein_info = protein_info,
    engine       = "spectronaut",
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
       "\u5f53\u524d\u7248\u672c\u652f\u6301: FragPipe, Spectronaut")
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
       "\u5f53\u524d\u7248\u672c\u652f\u6301: FragPipe, Spectronaut")
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
       "\u5f53\u524d\u7248\u672c\u652f\u6301: FragPipe, Spectronaut")
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


# ==============================================================================
# Spectronaut 内部辅助函数
# ==============================================================================

#' 通过读取表头判断是否为 Spectronaut 文件
#' @keywords internal
.detect_spectronaut_header <- function(fpath) {
  tryCatch({
    header <- readLines(fpath, n = 1, warn = FALSE)
    # Spectronaut 表头特征: PG. 前缀列
    grepl("PG\\.ProteinGroups|PG\\.Genes|PG\\.ProteinAccessions", header) ||
    # 或者列名含 .PG.Quantity 后缀
    grepl("\\.PG\\.Quantity", header)
  }, error = function(e) FALSE)
}

#' 提取 Spectronaut 定量列名和样本名
#'
#' 支持两种 Spectronaut 宽格式导出:
#'   格式A: 列名为干净样本名 (非 PG. 开头), 如 "PlasmaX-1", "PlasmaX-2"
#'   格式B: 列名含索引+文件名+后缀, 如 "[1] filename.d.PG.Quantity"
#'
#' @keywords internal
.extract_sn_quant_cols <- function(raw_cols, is_pg_col) {
  # --- 策略 1: 检测 .PG.Quantity / .PG.Log2Quantity 后缀 (格式B) ---
  idx_qty <- grep("\\.PG\\.Quantity$", raw_cols)
  idx_log2 <- grep("\\.PG\\.Log2Quantity$", raw_cols)

  if (length(idx_qty) > 0) {
    cols <- raw_cols[idx_qty]
    # 样本名提取: 去掉 [N] 前缀 和 .PG.Quantity 后缀
    samples <- stringr::str_remove(cols, "^\\[\\d+\\]\\s*")
    samples <- stringr::str_remove(samples, "\\.PG\\.Quantity$")
    # 去掉 .d 后缀 (Bruker raw file extension)
    samples <- stringr::str_remove(samples, "\\.d$")
    return(list(cols = cols, samples = samples, type = "PG.Quantity"))
  }

  if (length(idx_log2) > 0) {
    cols <- raw_cols[idx_log2]
    samples <- stringr::str_remove(cols, "^\\[\\d+\\]\\s*")
    samples <- stringr::str_remove(samples, "\\.PG\\.Log2Quantity$")
    samples <- stringr::str_remove(samples, "\\.d$")
    return(list(cols = cols, samples = samples, type = "PG.Log2Quantity"))
  }

  # --- 策略 2: 非 PG. 列即为定量列 (格式A) ---
  non_pg_idx <- which(!is_pg_col)
  if (length(non_pg_idx) > 0) {
    cols <- raw_cols[non_pg_idx]
    samples <- cols  # 列名就是样本名
    return(list(cols = cols, samples = samples, type = "PG.Quantity"))
  }

  list(cols = character(), samples = character(), type = "Unknown")
}

#' 构建 Spectronaut 标签名
#' @keywords internal
.build_sn_label_name <- function(df, raw_cols) {
  if ("PG.Genes" %in% raw_cols) {
    labels <- df$PG.Genes
    # 多基因用分号分隔，取第一个
    labels <- stringr::str_extract(labels, "^[^;]+")
  } else if ("PG.ProteinNames" %in% raw_cols) {
    labels <- df$PG.ProteinNames
    labels <- stringr::str_extract(labels, "^[^;]+")
  } else if ("PG.ProteinGroups" %in% raw_cols) {
    labels <- df$PG.ProteinGroups
  } else if ("PG.ProteinAccessions" %in% raw_cols) {
    labels <- df$PG.ProteinAccessions
  } else {
    labels <- as.character(seq_len(nrow(df)))
  }
  labels[is.na(labels)] <- "Unknown"
  labels
}
