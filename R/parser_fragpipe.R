# ==============================================================================
# msProteomiX - FragPipe parser
# ==============================================================================

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
    # 搜索子目录中的文件 (精确匹配 combined_protein.tsv，排除 template 等变体)
    file_prot <- .find_file(data_dir, "^combined_protein\\.(tsv|csv|txt)$")
    if (is.null(file_prot)) {
      # 尝试搜索子目录
      subdirs <- list.dirs(data_dir, recursive = FALSE, full.names = TRUE)
      for (sd in subdirs) {
        file_prot <- .find_file(sd, "^combined_protein\\.(tsv|csv|txt)$")
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
  proteins <- as.data.frame(lapply(prot_df[, quant_cols, drop = FALSE], function(x) suppressWarnings(as.numeric(x))))
  colnames(proteins) <- sample_names
  rownames(proteins) <- NULL

  # 智能 fallback: 如果 Intensity 列全为 0, 尝试 Spectral Count
  if (quant_info$type %in% c("LFQ", "Intensity")) {
    total_nonzero <- sum(rowSums(proteins, na.rm = TRUE) > 0)
    if (total_nonzero == 0) {
      message(">>> Intensity \u5217\u5168\u4e3a 0, \u5c1d\u8bd5\u4f7f\u7528 Spectral Count...")
      sc_idx <- grep("Spectral Count$", raw_cols, ignore.case = TRUE)
      sc_idx <- sc_idx[!grepl("(Unique|Total|Combined)", raw_cols[sc_idx], ignore.case = TRUE)]
      if (length(sc_idx) > 0) {
        sc_cols <- raw_cols[sc_idx]
        sc_samples <- trimws(stringr::str_remove(sc_cols, "(?i)\\s*Spectral Count$"))
        proteins <- as.data.frame(lapply(prot_df[, sc_cols, drop = FALSE], function(x) suppressWarnings(as.numeric(x))))
        colnames(proteins) <- sc_samples
        rownames(proteins) <- NULL
        sample_names <- sc_samples
        quant_info$type <- "SpectralCount"
        message(sprintf(">>> \u4f7f\u7528 Spectral Count (%d \u4e2a\u6837\u54c1)", length(sc_samples)))
      }
    }
  }

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
