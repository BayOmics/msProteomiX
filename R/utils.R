# ==============================================================================
# msProteomiX — 通用工具函数
# ==============================================================================

# NULL-coalescing operator (base R only has this from 4.4.0+)
`%||%` <- function(a, b) if (is.null(a)) b else a

#' @title msProteomiX 包内部环境
#' @description 用于在 RStudio Source 模式下跨脚本传递数据
#' @keywords internal
.msProteomiX_env <- new.env(parent = emptyenv())

#' 正则转义辅助函数
#'
#' 转义字符串中的正则特殊字符，保护样本名中的 +, ., (, ) 等
#'
#' @param s 要转义的字符串
#' @return 转义后的字符串
#' @export
#' @examples
#' escape_regex("Sample+1 (Rep.2)")
escape_regex <- function(s) {
  stringr::str_replace_all(s, "([.|()\\\\^{}+$*?]|\\[|\\])", "\\\\\\1")
}

#' 保存图表和数据的通用函数
#'
#' 将 ggplot 图表保存为 PDF，同时将源数据保存为 CSV
#'
#' @param plot_obj ggplot 对象
#' @param data_df 源数据 data.frame
#' @param project_name 项目名称 (用于文件名前缀)
#' @param suffix 文件名后缀 (如 "01_Protein_ID_Count")
#' @param output_dir 输出目录 (默认 "output")
#' @param width PDF 宽度 (英寸)
#' @param height PDF 高度 (英寸)
#' @return 不可见地返回保存的文件路径
#' @export
save_plot_and_data <- function(plot_obj, data_df, project_name, suffix,
                                output_dir = "output", width = 7, height = 6) {
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  base_name <- file.path(output_dir, paste0(project_name, "_", suffix))

  # 保存 CSV
  readr::write_csv(data_df, paste0(base_name, ".csv"))

  # 保存 PDF (仅 ggplot 对象)
  if ("ggplot" %in% class(plot_obj)) {
    ggplot2::ggsave(paste0(base_name, ".pdf"), plot_obj, width = width, height = height)
    print(plot_obj)
  }

  message(paste0("\u2705 \u5df2\u751f\u6210: ", suffix))
  Sys.sleep(0.3)


  invisible(paste0(base_name, c(".csv", ".pdf")))
}

#' msProteomiX 标准调色板
#'
#' 返回用于分组着色的标准调色板 (7色)
#'
#' @param n 需要的颜色数量 (最多7)
#' @return 颜色向量
#' @export
mspx_colors <- function(n = 7) {
  palette <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
               "#FF7F00", "#FFFF33", "#A65628")
  if (n > length(palette)) {
    palette <- rep(palette, length.out = n)
  }
  palette[1:n]
}

#' 获取 MsDataSet 中的样本名
#'
#' @param ms_data MsDataSet 对象
#' @return 样本名字符向量
#' @export
get_sample_names <- function(ms_data) {
  stopifnot(inherits(ms_data, "MsDataSet"))
  ms_data$sample_names
}

#' 创建 MsDataSet 对象
#'
#' @param proteins 蛋白水平 data.frame
#' @param peptides 肽段水平 data.frame (可选)
#' @param psms PSM 水平 data.frame (可选)
#' @param ions 离子水平 data.frame (可选)
#' @param sample_names 样本名向量
#' @param protein_info 蛋白注释信息 data.frame
#' @param engine 搜库引擎名称
#' @param quant_type 定量方式
#' @return MsDataSet 对象
#' @keywords internal
new_MsDataSet <- function(proteins = data.frame(),
                           peptides = data.frame(),
                           psms = data.frame(),
                           ions = data.frame(),
                           sample_names = character(),
                           protein_info = data.frame(),
                           engine = "unknown",
                           quant_type = "LFQ") {
  # Normalize tibble → data.frame to prevent subsetting gotchas
  protein_info <- as.data.frame(protein_info, stringsAsFactors = FALSE)

  # Detect if data is already log2-transformed (e.g. Spectronaut PG.Log2Quantity)
  is_log2 <- grepl("Log2", quant_type, ignore.case = TRUE)

  obj <- list(
    proteins     = proteins,
    peptides     = peptides,
    psms         = psms,
    ions         = ions,
    sample_names = sample_names,
    protein_info = protein_info,
    engine       = engine,
    quant_type   = quant_type,
    is_log2      = is_log2
  )
  class(obj) <- "MsDataSet"
  obj
}

#' 打印 MsDataSet 摘要
#' @param x MsDataSet 对象
#' @param ... 其他参数
#' @export
print.MsDataSet <- function(x, ...) {
  cat("== MsDataSet ==\n")
  cat(sprintf("  Engine:       %s\n", x$engine))
  cat(sprintf("  Quant Type:   %s%s\n", x$quant_type,
              if (isTRUE(x$is_log2)) " (log2)" else ""))
  cat(sprintf("  Samples:      %d\n", length(x$sample_names)))
  cat(sprintf("  Proteins:     %d\n", nrow(x$proteins)))
  if (nrow(x$peptides) > 0) cat(sprintf("  Peptides:     %d\n", nrow(x$peptides)))
  if (nrow(x$psms) > 0)     cat(sprintf("  PSMs:         %d\n", nrow(x$psms)))
  if (nrow(x$ions) > 0)     cat(sprintf("  Ions:         %d\n", nrow(x$ions)))
  cat(sprintf("  Sample Names: %s\n", paste(head(x$sample_names, 5), collapse = ", ")))
  if (length(x$sample_names) > 5) cat(sprintf("                ... and %d more\n", length(x$sample_names) - 5))
  invisible(x)
}

#' 确保输出目录存在
#' @keywords internal
ensure_output_dir <- function(dir = "output") {
  if (!dir.exists(dir)) dir.create(dir, recursive = TRUE)
}
