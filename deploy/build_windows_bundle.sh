#!/usr/bin/env bash
# ==============================================================================
# msProteomiX — Windows 部署包构建脚本
# 在 Mac 上运行，生成可直接交付给市场部的部署包
# ==============================================================================
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PKG_DIR="$(dirname "$SCRIPT_DIR")"
BUNDLE_NAME="msProteomiX_Windows部署包"
BUNDLE_DIR="${PKG_DIR}/${BUNDLE_NAME}"
TIMESTAMP=$(date +%Y%m%d_%H%M)

echo "============================================"
echo "  msProteomiX Windows 部署包构建"
echo "============================================"
echo ""

# --- 1. 检查 R 环境 ---
echo ">>> [1/4] Checking R environment..."
if ! command -v Rscript &> /dev/null; then
    echo "ERROR: Rscript not found. Please install R first."
    exit 1
fi

R_VERSION=$(Rscript -e 'cat(R.version.string)')
echo "    R: ${R_VERSION}"

# --- 2. document + build ---
echo ""
echo ">>> [2/4] Building package..."
cd "$PKG_DIR"

PKG_TAR=$(Rscript -e '
  suppressMessages(devtools::document())
  pkg_path <- devtools::build(quiet = TRUE)
  cat(pkg_path)
')

PKG_VERSION=$(Rscript -e 'cat(as.character(read.dcf("DESCRIPTION")[,"Version"]))')

echo "    Package: $(basename "$PKG_TAR")"
echo "    Version: ${PKG_VERSION}"

# --- 3. 组装部署包 ---
echo ""
echo ">>> [3/4] Assembling deployment bundle..."

# 清理旧的部署包
rm -rf "$BUNDLE_DIR"
mkdir -p "$BUNDLE_DIR"

# 复制 tar.gz
cp "$PKG_TAR" "$BUNDLE_DIR/"

# 复制安装脚本
cp "${SCRIPT_DIR}/windows_install.R" "$BUNDLE_DIR/"

# 生成 README
cat > "${BUNDLE_DIR}/安装说明.txt" << 'README_EOF'
╔══════════════════════════════════════════════════════════════╗
║           msProteomiX Windows 安装说明                       ║
╠══════════════════════════════════════════════════════════════╣
║                                                              ║
║  【前置条件】(按顺序安装，已安装可跳过)                      ║
║                                                              ║
║   1. R (>= 4.1.0)                                           ║
║      https://cran.r-project.org/bin/windows/base/            ║
║                                                              ║
║   2. Rtools (与 R 版本匹配, 如 R 4.4 -> Rtools44)           ║
║      https://cran.r-project.org/bin/windows/Rtools/          ║
║                                                              ║
║   3. RStudio Desktop (免费版)                                ║
║      https://posit.co/download/rstudio-desktop/              ║
║                                                              ║
║   4. Google Chrome (可选, 用于生成 PDF 报告)                 ║
║      https://www.google.com/chrome/                          ║
║                                                              ║
║  【安装 msProteomiX】                                        ║
║                                                              ║
║   1. 打开 RStudio                                            ║
║   2. File -> Open File -> 选择 "windows_install.R"           ║
║   3. 点击右上角 "Source" 按钮                                ║
║   4. 等待安装完成 (约 10-30 分钟)                            ║
║   5. 看到 "INSTALL COMPLETE" 表示成功                        ║
║                                                              ║
║  【创建分析项目】                                            ║
║                                                              ║
║   安装成功后, 在 RStudio Console 输入:                       ║
║                                                              ║
║     library(msProteomiX)                                     ║
║     create_project("D:/MyProject")                           ║
║                                                              ║
║   然后将搜库结果放入 wkdir/ 文件夹,                         ║
║   打开 scripts/01 脚本并 Source 运行即可。                   ║
║                                                              ║
║  【遇到问题?】                                               ║
║   联系 BayOmics 技术团队                                     ║
║                                                              ║
╚══════════════════════════════════════════════════════════════╝
README_EOF

# 写入版本信息
cat > "${BUNDLE_DIR}/VERSION.txt" << VERSION_EOF
msProteomiX ${PKG_VERSION}
Built: $(date '+%Y-%m-%d %H:%M:%S')
Builder: $(whoami)@$(hostname)
R: ${R_VERSION}
VERSION_EOF

# --- 4. 统计 ---
echo ""
echo ">>> [4/4] Done!"
echo ""
echo "============================================"
echo "  Bundle created: ${BUNDLE_NAME}/"
echo "============================================"
echo ""
ls -lh "$BUNDLE_DIR/"
echo ""
TOTAL_SIZE=$(du -sh "$BUNDLE_DIR" | cut -f1)
echo "  Total size: ${TOTAL_SIZE}"
echo "  Location:   ${BUNDLE_DIR}"
echo ""
echo "  Next: Copy this folder to USB or share via network."
echo ""
