# ==============================================================================
# msProteomiX — 通用工具函数
# ==============================================================================

# NULL-coalescing operator (base R only has this from 4.4.0+)
`%||%` <- function(a, b) if (is.null(a)) b else a

#' @title msProteomiX 包内部环境
#' @description 用于在 RStudio Source 模式下跨脚本传递数据
#' @keywords internal
.msProteomiX_env <- new.env(parent = emptyenv())
#' 创建 msProteomiX 分析项目
#'
#' 在指定目录创建标准项目结构，并从包中复制所有分析脚本。
#' 项目目录结构:
#' \preformatted{
#'   <project_dir>/
#'   ├── scripts/     ← 分析脚本 (从包复制)
#'   ├── wkdir/       ← 搜库结果文件放这里
#'   └── output/      ← 分析结果自动保存在这里
#' }
#'
#' @param project_dir 项目目录路径 (如 "~/my_project")
#' @return 不可见地返回项目路径
#' @export
#' @examples
#' \dontrun{
#' # 创建项目
#' create_project("~/Desktop/HeLa_APMS_2025")
#' # 然后在 RStudio 中打开 scripts/ 里的脚本, Source 运行即可
#' }
create_project <- function(project_dir) {
  project_dir <- normalizePath(project_dir, mustWork = FALSE)

  # 创建目录结构
  dirs <- c(
    file.path(project_dir, "scripts"),
    file.path(project_dir, "wkdir"),
    file.path(project_dir, "output")
  )
  for (d in dirs) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
  }

  # 复制脚本 (始终覆盖, 确保最新版本)
  n_copied <- .copy_pkg_scripts(file.path(project_dir, "scripts"))

  message(sprintf("\n>>> Project created: %s", project_dir))
  message(sprintf("   scripts/  <- %d scripts (latest version)", n_copied))
  message("   wkdir/    <- Place search results here")
  message("   output/   <- Analysis output saved here")
  message("\n>>> Next: Open scripts/01 in RStudio and Source")

  invisible(project_dir)
}


#' 更新项目脚本到最新版本
#'
#' 将项目 scripts/ 目录中的脚本更新为当前安装的 msProteomiX 包版本。
#' 仅覆盖包提供的脚本，用户自定义脚本不会被删除。
#'
#' @param scripts_dir 脚本目录路径 (默认自动检测当前项目)
#' @return 不可见地返回更新的文件数
#' @export
#' @examples
#' \dontrun{
#' update_scripts()                    # 自动检测当前项目
#' update_scripts("~/project/scripts") # 指定路径
#' }
update_scripts <- function(scripts_dir = NULL) {
  # 自动检测
  if (is.null(scripts_dir)) {
    if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
      ctx <- tryCatch(rstudioapi::getSourceEditorContext(), error = function(e) NULL)
      if (!is.null(ctx) && nchar(ctx$path) > 0) {
        sd <- dirname(ctx$path)
        if (basename(sd) == "scripts") scripts_dir <- sd
      }
    }
    if (is.null(scripts_dir)) {
      # 尝试从当前工作目录推断
      candidate <- file.path(getwd(), "scripts")
      if (dir.exists(candidate)) scripts_dir <- candidate
    }
    if (is.null(scripts_dir)) {
      stop("Cannot auto-detect scripts directory. Please specify scripts_dir.")
    }
  }

  n <- .copy_pkg_scripts(scripts_dir)
  message(sprintf(">>> Updated %d scripts in: %s", n, scripts_dir))
  invisible(n)
}


#' @keywords internal
.copy_pkg_scripts <- function(dest_dir) {
  pkg_scripts <- system.file("scripts", package = "msProteomiX")
  if (nchar(pkg_scripts) == 0) {
    stop("Cannot find msProteomiX scripts. Is the package installed?")
  }
  if (!dir.exists(dest_dir)) dir.create(dest_dir, recursive = TRUE)

  script_files <- list.files(pkg_scripts, full.names = TRUE)
  n_copied <- 0
  for (f in script_files) {
    dest <- file.path(dest_dir, basename(f))
    file.copy(f, dest, overwrite = TRUE)
    n_copied <- n_copied + 1
  }
  n_copied
}


#' 自动定位工作目录到项目根目录
#'
#' 在 RStudio 中运行时，自动将工作目录切换到脚本所在文件夹的
#' **上一级** (即项目根目录)。这样 wkdir/ 和 output/ 都在项目根目录下。
#'
#' 项目目录结构应为:
#' \preformatted{
#'   project_root/     ← 工作目录设在这里
#'   ├── scripts/      ← 脚本在这里
#'   ├── wkdir/
#'   └── output/
#' }
#'
#' @return 不可见地返回工作目录路径
#' @export
setup_workdir <- function() {
  project_dir <- NULL
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    ctx <- tryCatch(rstudioapi::getSourceEditorContext(), error = function(e) NULL)
    if (!is.null(ctx) && nchar(ctx$path) > 0) {
      script_dir <- dirname(ctx$path)
      if (basename(script_dir) == "scripts") {
        project_dir <- dirname(script_dir)
      } else {
        project_dir <- script_dir
      }
      setwd(project_dir)
      message(sprintf(">>> Working directory: %s", project_dir))
    }
  }
  if (is.null(project_dir)) {
    project_dir <- getwd()
    message(sprintf(">>> Working directory: %s", project_dir))
  }

  # 自动检查脚本是否需要更新
  scripts_dir <- file.path(project_dir, "scripts")
  if (dir.exists(scripts_dir)) {
    .check_scripts_version(scripts_dir)
  }

  invisible(project_dir)
}


#' @keywords internal
.check_scripts_version <- function(scripts_dir) {
  pkg_scripts <- system.file("scripts", package = "msProteomiX")
  if (nchar(pkg_scripts) == 0) return(invisible(NULL))

  # 比较每个脚本文件的 MD5
  pkg_files <- list.files(pkg_scripts, full.names = TRUE)
  outdated <- character(0)

  for (f in pkg_files) {
    local_f <- file.path(scripts_dir, basename(f))
    if (!file.exists(local_f)) {
      outdated <- c(outdated, basename(f))
    } else {
      pkg_md5 <- tools::md5sum(f)
      local_md5 <- tools::md5sum(local_f)
      if (pkg_md5 != local_md5) {
        outdated <- c(outdated, basename(f))
      }
    }
  }

  if (length(outdated) > 0) {
    message(sprintf("  >>> %d script(s) outdated. Run update_scripts() to sync.",
                    length(outdated)))
  }
}

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
  utils::write.csv(data_df, paste0(base_name, ".csv"), row.names = FALSE)

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
