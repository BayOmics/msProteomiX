# ╔══════════════════════════════════════════════════════════════╗
# ║           msProteomiX — 步骤 0: 安装指南                    ║
# ╠══════════════════════════════════════════════════════════════╣
# ║                                                              ║
# ║  【这是什么？】                                              ║
# ║   这个脚本会自动安装 msProteomiX 及其所有依赖包。           ║
# ║   只需运行一次，后续不需要再运行。                          ║
# ║                                                              ║
# ║  【操作步骤】                                                ║
# ║   1. 在 RStudio 中打开本文件                                ║
# ║   2. 点击右上角 "Source" 按钮                               ║
# ║   3. 等待安装完成 (可能需要几分钟)                          ║
# ║   4. 看到 "安装完成" 提示后即可关闭                         ║
# ║                                                              ║
# ║  【常见问题】                                                ║
# ║   Q: 安装失败怎么办？                                       ║
# ║   A: 检查网络连接，或尝试设置镜像:                          ║
# ║      options(repos = "https://cloud.r-project.org")         ║
# ║                                                              ║
# ╚══════════════════════════════════════════════════════════════╝

message("=== msProteomiX Installation ===\n")

# --- 1. 安装 CRAN 包 ---
cran_pkgs <- c(
  "devtools", "readr", "dplyr", "stringr", "ggplot2", "ggrepel",
  "reshape2", "scales", "tibble", "ggvenn", "UpSetR", "gridExtra"
)

message(">>> Installing CRAN packages...")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("  Installing:", pkg))
    tryCatch(
      install.packages(pkg, quiet = TRUE),
      error = function(e) message(paste("  FAILED:", pkg, "-", e$message))
    )
  } else {
    message(paste("  Already installed:", pkg))
  }
}

# --- 2. 安装 Bioconductor 包 ---
message("\n>>> Installing Bioconductor packages...")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c("limma", "Biostrings", "clusterProfiler", "org.Hs.eg.db")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("  Installing:", pkg))
    tryCatch(
      BiocManager::install(pkg, update = FALSE, ask = FALSE),
      error = function(e) message(paste("  FAILED:", pkg, "-", e$message))
    )
  } else {
    message(paste("  Already installed:", pkg))
  }
}

# --- 3. 安装 msProteomiX 本身 ---
message("\n>>> Installing msProteomiX...")
tryCatch({
  # 如果在包开发目录中，使用 devtools 安装
  pkg_dir <- getwd()
  if (file.exists(file.path(pkg_dir, "DESCRIPTION"))) {
    devtools::install(pkg_dir, dependencies = FALSE, upgrade = "never")
  } else {
    # 否则尝试从父目录安装
    parent <- dirname(dirname(pkg_dir))
    if (file.exists(file.path(parent, "DESCRIPTION"))) {
      devtools::install(parent, dependencies = FALSE, upgrade = "never")
    } else {
      message("  Please run this from within the msProteomiX project directory.")
    }
  }
}, error = function(e) {
  message(paste("  Installation error:", e$message))
})

# --- 4. 验证 ---
message("\n>>> Verifying installation...")
if (requireNamespace("msProteomiX", quietly = TRUE)) {
  message("\n========================================")
  message("  \u2705 \u5b89\u88c5\u5b8c\u6210!")
  message("========================================")
  message("")
  message("  \u4e0b\u4e00\u6b65: \u521b\u5efa\u5206\u6790\u9879\u76ee\u6587\u4ef6\u5939:")
  message('  msProteomiX::create_project("~/Desktop/\u6211\u7684\u86cb\u767d\u7ec4\u5b66\u9879\u76ee")')
  message("")
  message("  \u7136\u540e\u5728 RStudio \u4e2d\u6253\u5f00 scripts/ \u91cc\u7684\u811a\u672c\u8fd0\u884c")
} else {
  message("\n========================================")
  message("  \u26a0\ufe0f msProteomiX \u672a\u627e\u5230")
  message("  \u8bf7\u91cd\u542f RStudio \u540e\u518d\u6b21\u8fd0\u884c")
  message("========================================")
}
