# ==============================================================================
# msProteomiX - Spectronaut parser
# ==============================================================================

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
    # 排除 Precise_Report (它是肽段级数据, 不是蛋白级报告)
    candidates <- candidates[!grepl("Precise_Report|Peptide_Report", basename(candidates), ignore.case = TRUE)]
    if (length(candidates) == 0) {
      # 尝试搜索所有 tsv/csv 文件并通过表头检测
      all_files <- list.files(path, pattern = "\\.(tsv|csv)$",
                               full.names = TRUE, ignore.case = TRUE)
      # 排除 IdentificationsOverview 和 group_info 等辅助文件
      all_files <- all_files[!grepl("(group_info|IdentificationsOverview|Precise_Report|Peptide_Report|Overview)",
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
        sub_candidates <- sub_candidates[!grepl("Precise_Report|Peptide_Report", basename(sub_candidates),
                                                  ignore.case = TRUE)]
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
    # Spectronaut 用 NaN 表示缺失值; "Filtered" 也可能出现
    x <- suppressWarnings(as.numeric(x))
    x[is.nan(x)] <- NA
    x
  }))
  colnames(proteins) <- sample_names
  rownames(proteins) <- NULL

  # --- 过滤全 NA 行 (所有样本均无定量值的蛋白) ---
  valid_rows <- rowSums(!is.na(proteins)) > 0
  if (any(!valid_rows)) {
    n_removed <- sum(!valid_rows)
    message(sprintf("  >>> Removed %d all-missing proteins", n_removed))
    proteins <- proteins[valid_rows, , drop = FALSE]
    df <- df[valid_rows, , drop = FALSE]
  }

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

  # --- C2 修复: 确保 Protein 列存在 (下游函数依赖) ---
  if (!("Protein" %in% colnames(protein_info))) {
    if ("Protein ID" %in% colnames(protein_info)) {
      protein_info$Protein <- protein_info$`Protein ID`
    } else if ("Entry Name" %in% colnames(protein_info)) {
      protein_info$Protein <- protein_info$`Entry Name`
    } else if ("Gene" %in% colnames(protein_info)) {
      protein_info$Protein <- protein_info$Gene
    }
  }
  if (!("Protein ID" %in% colnames(protein_info))) {
    if ("Protein" %in% colnames(protein_info)) {
      protein_info$`Protein ID` <- protein_info$Protein
    }
  }

  # --- H1 修复: Gene 列去除 trailing 分号 ---
  if ("Gene" %in% colnames(protein_info)) {
    protein_info$Gene <- stringr::str_remove(protein_info$Gene, ";+$")
  }

  # 构建标签名 (用于火山图标注)
  protein_info$Label_Name <- .build_sn_label_name(df, raw_cols)

  # --- 可选: 解析 Precise_Report.csv (肽段/前体水平数据, 用于 QC) ---
  psm_df <- .parse_sn_precise_report(report_file)

  new_MsDataSet(
    proteins     = proteins,
    peptides     = data.frame(),
    psms         = psm_df,
    ions         = data.frame(),
    sample_names = sample_names,
    protein_info = protein_info,
    engine       = "spectronaut",
    quant_type   = quant_info$type
  )
}


# ==============================================================================
# MaxQuant 解析器

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


# ==============================================================================
# Precise_Report 解析 (可选)
# ==============================================================================

#' 搜索 Precise_Report.csv 文件
#' @keywords internal
.sn_find_precise_report <- function(report_file) {
  report_dir <- dirname(report_file)
  candidates <- list.files(report_dir, pattern = "Precise_Report\\.(csv|tsv)$",
                             full.names = TRUE, ignore.case = TRUE)
  if (length(candidates) > 0) return(candidates[1])
  NULL
}


#' 解析 Spectronaut Precise_Report (肽段/前体水平, 可选)
#'
#' 读取 Precise_Report.csv 并标准化列名，用于 QC 函数。
#' 只保留需要的 8 列以节省内存。
#'
#' @param report_file Report.csv 的路径（用于定位同目录的 Precise_Report）
#' @return data.frame（标准化列名）或空 data.frame
#' @keywords internal
.parse_sn_precise_report <- function(report_file) {
  precise_file <- .sn_find_precise_report(report_file)

  if (is.null(precise_file)) {
    message("  >>> Precise_Report not found. PSM-level QC plots will be skipped.")
    return(data.frame())
  }

  message(sprintf(">>> Reading Precise_Report: %s", basename(precise_file)))

  # 先读取表头确定可用列
  header <- readLines(precise_file, n = 1, warn = FALSE)
  sep <- if (grepl("\t", header)) "\t" else ","
  all_cols <- unlist(strsplit(header, sep))

  # 定义需要的列 (Spectronaut → 标准名)
  col_map <- c(
    "PEP.StrippedSequence"       = "Peptide",
    "EG.ModifiedSequence"        = "Modified Peptide",
    "PEP.NrOfMissedCleavages"    = "Number of Missed Cleavages",
    "FG.Charge"                  = "Charge",
    "FG.PrecMz"                  = "Calibrated Observed M/Z",
    "FG.PrecMzCalibrated"        = "Calibrated Observed M/Z",
    "FG.Quantity"                = "Intensity",
    "R.Label"                    = "Spectrum File",
    "R.FileName"                 = "Run",
    "PG.ProteinAccessions"       = "Protein",
    "EG.ApexRT"                  = "Retention"
  )

  # 找到实际存在的列
  avail <- intersect(names(col_map), all_cols)
  if (length(avail) == 0) {
    message("  >>> Precise_Report columns not recognized. Skipping.")
    return(data.frame())
  }

  # 构建 colClasses: 只读需要的列, 其余设为 NULL
  col_classes <- rep("NULL", length(all_cols))
  names(col_classes) <- all_cols
  for (cn in avail) {
    col_classes[cn] <- NA  # NA = auto detect type
  }

  # 读取数据 (只加载需要的列)
  psm_df <- tryCatch({
    read.csv(precise_file, sep = sep, stringsAsFactors = FALSE,
             check.names = FALSE, colClasses = col_classes, quote = "\"")
  }, error = function(e) {
    message("  >>> Failed to read Precise_Report: ", e$message)
    return(data.frame())
  })

  if (nrow(psm_df) == 0) return(data.frame())

  # 标准化列名
  rename_map <- col_map[colnames(psm_df)]
  rename_map <- rename_map[!is.na(rename_map)]

  # 处理重复目标名 (FG.PrecMz 和 FG.PrecMzCalibrated 都映射到同一个名字)
  # 优先使用 FG.PrecMzCalibrated
  if ("FG.PrecMzCalibrated" %in% colnames(psm_df) && "FG.PrecMz" %in% colnames(psm_df)) {
    psm_df$FG.PrecMz <- NULL
    rename_map <- rename_map[names(rename_map) != "FG.PrecMz"]
  }

  for (old_name in names(rename_map)) {
    new_name <- rename_map[[old_name]]
    idx <- which(colnames(psm_df) == old_name)
    if (length(idx) == 1) colnames(psm_df)[idx] <- new_name
  }

  # NaN → NA
  for (j in seq_along(psm_df)) {
    if (is.numeric(psm_df[[j]])) {
      psm_df[[j]][is.nan(psm_df[[j]])] <- NA
    }
  }

  n_rows <- nrow(psm_df)
  n_files <- if ("Spectrum File" %in% colnames(psm_df)) {
    length(unique(psm_df[["Spectrum File"]]))
  } else {
    NA_integer_
  }
  message(sprintf("  >>> Loaded %d precursor records from %d runs", n_rows, n_files))

  psm_df
}
