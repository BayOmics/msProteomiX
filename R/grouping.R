# ==============================================================================
# msProteomiX — 样品分组模块
# ==============================================================================

#' 交互式样品分组向导
#'
#' 在 RStudio Console 中通过交互式对话引导用户对样品进行分组。
#' 支持自动保存/读取分组信息文件，避免重复分组。
#'
#' @param sample_names 样品名字符向量
#' @param group_file 分组信息文件路径。如果文件已存在，会提示用户是否复用。
#'   默认为 NULL (不自动保存)。
#' @param context 当前分析上下文名称，用于显示提示信息 (如 "数据导入")
#' @return data.frame，包含 sample_name 和 user_group 两列
#' @export
#' @examples
#' \dontrun{
#' groups <- interactive_grouping(c("Sample1", "Sample2", "Sample3"),
#'                                 group_file = "wkdir/group_info.csv")
#' }
interactive_grouping <- function(sample_names, group_file = NULL, context = "Analysis") {
  if (!interactive()) {
    stop("Please run this script in RStudio using the Source button.")
  }

  final_groups <- NULL

  # --- 尝试读取已有分组文件 ---
  if (!is.null(group_file) && file.exists(group_file)) {
    old_groups <- tryCatch(
      readr::read_csv(group_file, show_col_types = FALSE),
      error = function(e) NULL
    )

    if (!is.null(old_groups) && "sample_name" %in% colnames(old_groups) &&
        "user_group" %in% colnames(old_groups)) {
      f_samps <- trimws(old_groups$sample_name)
      d_samps <- trimws(sample_names)
      common <- intersect(d_samps, f_samps)

      if (length(common) > 0) {
        cat("\n")
        message("==========================================================")
        message(paste0("  Detected existing group file: ", basename(group_file)))
        message(paste0("  Matched samples: ", length(common), " / ", length(d_samps)))
        message("----------------------------------------------------------")
        message("1. Use existing groups (recommended)")
        message("2. Re-group (overwrite)")
        message("==========================================================")
        flush.console()
        choice <- readline(">>> Please select (1/2): ")

        if (choice == "1") {
          ret_grps <- rep("Undefined", length(d_samps))
          match_idx <- match(d_samps, f_samps)
          valid <- !is.na(match_idx)
          ret_grps[valid] <- old_groups$user_group[match_idx[valid]]
          final_groups <- ret_grps
          message("  Successfully loaded existing group file!")
        }
      }
    }
  }

  # --- 新建分组 (交互式向导) ---
  if (is.null(final_groups)) {
    final_groups <- .run_grouping_wizard(sample_names, context)
  }

  # --- 构建结果 data.frame ---
  result <- data.frame(
    sample_name = sample_names,
    user_group  = final_groups,
    stringsAsFactors = FALSE
  )

  # --- 保存分组文件 ---
  if (!is.null(group_file)) {
    group_dir <- dirname(group_file)
    if (!dir.exists(group_dir)) dir.create(group_dir, recursive = TRUE)
    readr::write_csv(result, group_file)
    message(paste0("  Group info saved to: ", group_file))
  }

  # 过滤掉未分组的样本
  result_clean <- result[result$user_group != "Undefined", ]

  if (nrow(result_clean) == 0) {
    stop("Error: no valid group assignments.")
  }

  result_clean
}


#' 从文件读取分组信息
#'
#' 直接读取已保存的分组 CSV 文件，不进行交互。
#'
#' @param file_path 分组文件路径
#' @return data.frame，包含 sample_name 和 user_group
#' @export
read_group_info <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("Group file not found: ", file_path,
         "\nPlease run 01_data_import_and_grouping.R first.")
  }
  df <- readr::read_csv(file_path, show_col_types = FALSE)
  if (!all(c("sample_name", "user_group") %in% colnames(df))) {
    stop("Invalid group file format: missing sample_name or user_group column.")
  }
  df[df$user_group != "Undefined", ]
}


#' 编程式分组 (非交互)
#'
#' 直接指定分组，适用于脚本自动化场景。
#'
#' @param sample_names 样品名向量
#' @param groups 分组名向量 (与 sample_names 等长)
#' @return data.frame，包含 sample_name 和 user_group
#' @export
#' @examples
#' groups <- set_groups(
#'   c("S1", "S2", "S3", "S4"),
#'   c("Control", "Control", "Treatment", "Treatment")
#' )
set_groups <- function(sample_names, groups) {
  if (length(sample_names) != length(groups)) {
    stop("sample_names and groups must have the same length.")
  }
  data.frame(
    sample_name = sample_names,
    user_group  = groups,
    stringsAsFactors = FALSE
  )
}


#' 检查前置脚本是否已执行
#'
#' 检查 .msProteomiX_env 中是否已有数据和分组信息。
#' 如果未执行前置脚本，给出详细的中文操作提示。
#'
#' @param step 需要检查的步骤编号 (如 "01")
#' @return 不可见地返回 TRUE
#' @export
check_prerequisites <- function(step = "01") {
  env <- .msProteomiX_env

  has_data <- exists("ms_data", envir = env) && !is.null(env$ms_data)
  has_group <- exists("group_info", envir = env) && !is.null(env$group_info)

  if (!has_data || !has_group) {
    message("===================================================")
    message("  WARNING: Please run 01_data_import_and_grouping.R first!")
    message("")
    message("  How to do it:")
    message("  1. In RStudio, find inst/scripts/ folder in the Files panel")
    message("  2. Open 01_data_import_and_grouping.R")
    message("  3. Click the 'Source' button (top-right)")
    message("  4. After it finishes, come back to run this script")
    message("===================================================")
    stop("Prerequisites not met.")
  }
  invisible(TRUE)
}


# ==============================================================================
# 内部辅助函数
# ==============================================================================

#' 交互式分组向导 (内部)
#' @keywords internal
.run_grouping_wizard <- function(sample_names, context = "Analysis") {
  groups <- rep("Undefined", length(sample_names))

  repeat {
    cat("\n")
    message("==========================================================")
    message(paste0("  Sample Grouping - ", context))
    message("----------------------------------------------------------")
    print(data.frame(
      Index  = seq_along(sample_names),
      Sample = sample_names,
      Group  = groups
    ))
    message("----------------------------------------------------------")
    message("Instructions:")
    message("  - Enter sample index range, e.g.: 1-3 or 1,3,5")
    message("  - Then enter a group name, e.g.: Control")
    message("  - Enter 'y' to confirm grouping")
    message("==========================================================")
    flush.console()

    input_str <- readline(">>> Enter command: ")

    if (tolower(input_str) == "y") {
      if (any(groups == "Undefined")) {
        confirm <- readline("  Some samples are ungrouped and will be excluded. Continue? (y/n): ")
        if (tolower(confirm) == "y") break
      } else {
        message("  Grouping confirmed!")
        break
      }
    } else {
      # 解析序号
      idx <- suppressWarnings(tryCatch({
        parts <- unlist(strsplit(input_str, "[,\uff0c]"))
        vec <- c()
        for (p in parts) {
          p <- trimws(p)
          if (grepl("-", p)) {
            r <- as.numeric(unlist(strsplit(p, "-")))
            if (length(r) == 2 && !any(is.na(r))) vec <- c(vec, r[1]:r[2])
          } else if (nchar(p) > 0) {
            vec <- c(vec, as.numeric(p))
          }
        }
        unique(vec[!is.na(vec) & vec > 0 & vec <= length(sample_names)])
      }, error = function(e) NULL))

      if (length(idx) > 0) {
        cat(sprintf("  Selected %d samples.\n", length(idx)))
        nm <- readline(">>> Enter group name (e.g. Control): ")
        if (nchar(nm) > 0) groups[idx] <- nm
      } else {
        message("  Invalid input, please try again.")
        Sys.sleep(0.5)
      }
    }
  }
  groups
}
