# ==============================================================================
# msProteomiX - Core parsers & shared helpers
# ==============================================================================

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
    sn_candidates <- sn_candidates[!grepl("Precise_Report", sn_candidates, ignore.case = TRUE)]
    if (length(sn_candidates) > 0) {
      full_path <- file.path(path, sn_candidates[1])
      if (.detect_spectronaut_header(full_path)) return("spectronaut")
    }

    # Spectronaut fallback: check any Report*.csv by header (catches "Report.csv")
    sn_report <- files[grepl("^Report\\.(tsv|csv)$", basename(files), ignore.case = TRUE)]
    sn_report <- sn_report[!grepl("Precise_Report", sn_report, ignore.case = TRUE)]
    if (length(sn_report) > 0) {
      full_path <- file.path(path, sn_report[1])
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

#' 在目录中查找匹配文件
#' @keywords internal
.find_file <- function(dir, pattern) {
  f <- list.files(dir, pattern = pattern, full.names = TRUE, ignore.case = TRUE)
  if (length(f) > 0) f[1] else NULL
}

#' 读取 TSV/CSV 文件
#'
#' 使用 base R 读取，避免 readr/vroom 在 Apple Silicon RStudio 上的内存问题
#' @keywords internal
.read_omics_file <- function(fpath) {
  if (grepl("\\.csv$", fpath, ignore.case = TRUE)) {
    utils::read.csv(fpath, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    utils::read.delim(fpath, sep = "\t", stringsAsFactors = FALSE,
                      check.names = FALSE, strip.white = TRUE)
  }
}

#' 提取 FragPipe 定量列名和样本名
#'
#' FragPipe combined_protein.tsv 的定量列格式为:
#'   "SampleName MaxLFQ Intensity" (样本名在前，后缀固定)
#'

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



