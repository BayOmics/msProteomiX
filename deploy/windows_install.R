# ╔══════════════════════════════════════════════════════════════╗
# ║     msProteomiX — Windows 一键安装脚本                      ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【操作方法】                                                ║
# ║   1. 用 RStudio 打开本文件                                  ║
# ║   2. 点击右上角 "Source" 按钮                                ║
# ║   3. 等待安装完成 (约 10-30 分钟, 取决于网速)               ║
# ║   4. 看到 "INSTALL COMPLETE" 即为成功                       ║
# ║                                                              ║
# ║  【注意事项】                                                ║
# ║   - 首次运行需联网 (下载依赖包)                             ║
# ║   - 请勿关闭 RStudio 窗口                                   ║
# ║   - 如遇报错, 重启 RStudio 后再次 Source 即可               ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

# ==============================================================================
# 环境检测
# ==============================================================================

.installer_start_time <- Sys.time()

message("==================================================")
message("  msProteomiX Windows Installer")
message("==================================================\n")

# --- 检测 R 版本 ---
r_ver <- getRversion()
message(sprintf("  R version:  %s", r_ver))
if (r_ver < "4.1.0") {
  stop("R >= 4.1.0 is required. Please update R from: https://cran.r-project.org/bin/windows/base/")
}

# --- 检测操作系统 ---
message(sprintf("  OS:         %s", .Platform$OS.type))
message(sprintf("  Platform:   %s", R.version$platform))

# --- 检测 Rtools (Windows only) ---
if (.Platform$OS.type == "windows") {
  rtools_ok <- nzchar(Sys.which("make"))
  if (rtools_ok) {
    message(sprintf("  Rtools:     OK (%s)", Sys.which("make")))
  } else {
    message("  Rtools:     NOT FOUND")
    message("")
    message("  !! WARNING: Rtools is required to compile some packages.")
    message("  !! Download from: https://cran.r-project.org/bin/windows/Rtools/")
    message("  !! Install Rtools, restart RStudio, then run this script again.")
    message("")
    # NOTE: readline() returns "" silently in RStudio Source mode,
    #       so use utils::menu() which works in both interactive and Source mode.
    choice <- utils::menu(
      choices = c("Continue anyway", "Abort (install Rtools first)"),
      title = "Rtools not found. What would you like to do?"
    )
    if (choice != 1) {
      stop("Aborted. Please install Rtools first.")
    }
  }
}

message("")

# ==============================================================================
# 镜像配置 (国内用户加速)
# ==============================================================================

# CRAN: 清华镜像 (备选: 中科大 https://mirrors.ustc.edu.cn/CRAN/)
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# Bioconductor: 清华镜像
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")

message(">>> Mirror: CRAN (Tsinghua) + Bioconductor (Tsinghua)\n")

# ==============================================================================
# 辅助函数
# ==============================================================================

.safe_install_cran <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    v <- as.character(packageVersion(pkg))
    message(sprintf("  [OK]  %-20s %s", pkg, v))
    return(TRUE)
  }
  message(sprintf("  [..] Installing %-15s ...", pkg))
  tryCatch({
    install.packages(pkg, quiet = TRUE)
    if (requireNamespace(pkg, quietly = TRUE)) {
      v <- as.character(packageVersion(pkg))
      message(sprintf("  [OK]  %-20s %s", pkg, v))
      return(TRUE)
    } else {
      message(sprintf("  [!!]  %-20s install returned but package not loadable", pkg))
      return(FALSE)
    }
  }, error = function(e) {
    message(sprintf("  [!!]  %-20s FAILED: %s", pkg, e$message))
    return(FALSE)
  })
}

.safe_install_bioc <- function(pkg) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    v <- as.character(packageVersion(pkg))
    message(sprintf("  [OK]  %-20s %s", pkg, v))
    return(TRUE)
  }
  message(sprintf("  [..] Installing %-15s ... (Bioconductor, may take a while)", pkg))
  tryCatch({
    BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
    if (requireNamespace(pkg, quietly = TRUE)) {
      v <- as.character(packageVersion(pkg))
      message(sprintf("  [OK]  %-20s %s", pkg, v))
      return(TRUE)
    } else {
      message(sprintf("  [!!]  %-20s install returned but package not loadable", pkg))
      return(FALSE)
    }
  }, error = function(e) {
    message(sprintf("  [!!]  %-20s FAILED: %s", pkg, e$message))
    return(FALSE)
  })
}

# ==============================================================================
# Step 1: CRAN 包
# ==============================================================================

message(">>> [Step 1/4] CRAN packages\n")

cran_pkgs <- c(
  # --- Core dependencies (Imports) ---
  "dplyr", "magrittr", "rlang", "stringr",
  "ggplot2", "ggrepel", "reshape2", "scales", "tibble",
  # --- Suggests (needed by analysis scripts) ---
  # NOTE: "grid" is a base R package, no install needed
  "devtools", "ggvenn", "UpSetR", "gridExtra",
  "pheatmap", "igraph", "knitr", "rmarkdown", "base64enc"
)

cran_results <- vapply(cran_pkgs, .safe_install_cran, logical(1))
cran_ok <- sum(cran_results)
cran_fail <- sum(!cran_results)
message(sprintf("\n  CRAN: %d/%d OK", cran_ok, length(cran_pkgs)))
if (cran_fail > 0) {
  message(sprintf("  !! %d package(s) failed. Will retry after Bioconductor.\n", cran_fail))
}

# ==============================================================================
# Step 2: Bioconductor 包
# ==============================================================================

message("\n>>> [Step 2/4] Bioconductor packages\n")

# 安装 BiocManager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  message("  [..] Installing BiocManager ...")
  install.packages("BiocManager", quiet = TRUE)
}
message(sprintf("  [OK]  %-20s %s", "BiocManager", as.character(packageVersion("BiocManager"))))

bioc_pkgs <- c(
  "limma",              # 差异分析
  "Biostrings",         # 序列覆盖度
  "clusterProfiler",    # GO/KEGG 富集
  "org.Hs.eg.db",       # 人类基因注释 (大约 80 MB)
  "DOSE",               # clusterProfiler 依赖
  "ReactomePA"          # Reactome 通路富集
)

bioc_results <- vapply(bioc_pkgs, .safe_install_bioc, logical(1))
bioc_ok <- sum(bioc_results)
bioc_fail <- sum(!bioc_results)
message(sprintf("\n  Bioconductor (core): %d/%d OK", bioc_ok, length(bioc_pkgs)))

# --- Optional Bioconductor packages ---
# These are large or rarely needed; install with tryCatch so failures are non-fatal
message("\n  Installing optional packages (failures are non-fatal)...")

optional_bioc <- c(
  "org.Mm.eg.db",       # 小鼠基因注释 (仅分析小鼠数据时需要, ~80 MB)
  "STRINGdb"            # STRING PPI 网络 (需联网查询)
)
for (pkg in optional_bioc) {
  tryCatch(
    .safe_install_bioc(pkg),
    error = function(e) {
      message(sprintf("  [--]  %-20s skipped (optional): %s", pkg, e$message))
    }
  )
}

# (bioc_results moved above with optional packages)

# ==============================================================================
# Step 3: 安装 msProteomiX
# ==============================================================================

message("\n>>> [Step 3/4] Installing msProteomiX\n")

# 定位当前脚本所在目录
script_dir <- tryCatch({
  dirname(rstudioapi::getSourceEditorContext()$path)
}, error = function(e) {
  getwd()
})

# 查找本地 tar.gz
tar_files <- list.files(script_dir, pattern = "^msProteomiX.*\\.tar\\.gz$",
                        full.names = TRUE)

if (length(tar_files) > 0) {
  # 如有多个版本，取最新的
  tar_file <- sort(tar_files, decreasing = TRUE)[1]
  message(sprintf("  Found local package: %s", basename(tar_file)))
  message("  Installing from local file ...")
  tryCatch({
    install.packages(tar_file, repos = NULL, type = "source")
  }, error = function(e) {
    message(sprintf("  Local install failed: %s", e$message))
    message("  Trying GitHub fallback ...")
    if (requireNamespace("devtools", quietly = TRUE)) {
      devtools::install_github("BayOmics/msProteomiX", upgrade = "never")
    } else {
      stop("Cannot install msProteomiX: devtools not available.")
    }
  })
} else {
  # 没有本地包, 从 GitHub 安装
  message("  No local .tar.gz found. Installing from GitHub ...")
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::install_github("BayOmics/msProteomiX", upgrade = "never")
  } else {
    stop("Cannot install msProteomiX: neither local .tar.gz nor devtools available.")
  }
}

# ==============================================================================
# Step 4: 验证
# ==============================================================================

message("\n>>> [Step 4/4] Verification\n")

# --- 验证核心包 ---
all_ok <- TRUE

# 核心检查
critical_pkgs <- c("msProteomiX", "dplyr", "ggplot2", "limma", "clusterProfiler")
for (pkg in critical_pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    v <- as.character(packageVersion(pkg))
    message(sprintf("  [OK]  %-20s %s", pkg, v))
  } else {
    message(sprintf("  [!!]  %-20s NOT FOUND", pkg))
    all_ok <- FALSE
  }
}

# --- Chrome 检查 (非必需) ---
chrome_found <- FALSE
if (.Platform$OS.type == "windows") {
  chrome_paths <- c(
    file.path(Sys.getenv("ProgramFiles"), "Google/Chrome/Application/chrome.exe"),
    file.path(Sys.getenv("ProgramFiles(x86)"), "Google/Chrome/Application/chrome.exe"),
    file.path(Sys.getenv("LOCALAPPDATA"), "Google/Chrome/Application/chrome.exe")
  )
  for (p in chrome_paths) {
    if (nzchar(p) && file.exists(p)) {
      message(sprintf("  [OK]  %-20s %s", "Google Chrome", "found"))
      chrome_found <- TRUE
      break
    }
  }
} else {
  mac_chrome <- "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome"
  if (file.exists(mac_chrome)) {
    chrome_found <- TRUE
    message(sprintf("  [OK]  %-20s %s", "Google Chrome", "found"))
  }
}
if (!chrome_found) {
  message(sprintf("  [--]  %-20s not found (PDF reports will fallback to HTML)", "Google Chrome"))
}

# --- 总耗时 ---
elapsed <- round(as.numeric(difftime(Sys.time(), .installer_start_time, units = "mins")), 1)

# --- 最终结果 ---
message("")
message("==================================================")
if (all_ok) {
  pkg_ver <- as.character(packageVersion("msProteomiX"))
  message(sprintf("  INSTALL COMPLETE  -  msProteomiX v%s", pkg_ver))
  message(sprintf("  Time elapsed: %s min", elapsed))
  message("==================================================")
  message("")
  message("  Next steps:")
  message("  ------------------------------------------------")
  message("  1. In RStudio Console, type:")
  message("")
  message("     library(msProteomiX)")
  message("     create_project(\"D:/MyProject\")")
  message("")
  message("  2. Place search engine output in wkdir/ folder")
  message("  3. Open scripts/01 in RStudio, click Source")
  message("  ------------------------------------------------")
} else {
  message("  INSTALL INCOMPLETE  -  some packages failed")
  message(sprintf("  Time elapsed: %s min", elapsed))
  message("==================================================")
  message("")
  message("  Troubleshooting:")
  message("  - Restart RStudio and run this script again")
  message("  - Check internet connection")
  message("  - Make sure Rtools is installed (Windows)")
  message("  - Contact BayOmics tech support")
}
message("")
