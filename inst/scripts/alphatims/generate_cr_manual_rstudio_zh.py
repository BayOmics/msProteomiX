#!/usr/bin/env python3
"""
generate_cr_manual_rstudio_zh.py — Build Chinese RStudio CR Manual.

Output:
    inst/scripts/alphatims/CR_Manual_RStudio_ZH.html
"""
import base64
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ASSETS_DIR = os.path.join(
    os.path.dirname(SCRIPT_DIR), "..", "..", "wkdir", "cr_manual_assets"
)
if not os.path.isdir(ASSETS_DIR):
    ASSETS_DIR = os.path.join(SCRIPT_DIR, "assets")

OUTPUT_HTML = os.path.join(SCRIPT_DIR, "CR_Manual_RStudio_ZH.html")


def img_base64(path):
    if not os.path.exists(path):
        print(f"  WARNING: Image not found: {path}")
        return ""
    with open(path, "rb") as f:
        data = base64.b64encode(f.read()).decode()
    ext = os.path.splitext(path)[1].lower().replace(".", "")
    if ext == "jpg":
        ext = "jpeg"
    return f"data:image/{ext};base64,{data}"


def img_tag(path, caption="", width="90%"):
    src = img_base64(path)
    if not src:
        return f'<p style="color:#999;">[图片未找到: {os.path.basename(path)}]</p>'
    return (
        f'<figure style="text-align:center;margin:20px auto;">'
        f'<img src="{src}" style="max-width:{width};border:1px solid #ddd;'
        f'border-radius:8px;box-shadow:0 2px 8px rgba(0,0,0,0.1);">'
        f'<figcaption style="color:#666;font-size:12px;margin-top:8px;">'
        f'{caption}</figcaption></figure>'
    )


CSS = """
<style>
  * { box-sizing: border-box; }
  body {
    font-family: "Microsoft YaHei", "PingFang SC", "Segoe UI", sans-serif;
    max-width: 900px; margin: 0 auto; padding: 30px 20px;
    background: #fafbfc; color: #333; line-height: 1.8;
  }
  h1 { color: #1a1a2e; border-bottom: 3px solid #2196F3; padding-bottom: 12px;
       font-size: 26px; }
  h2 { color: #16213e; margin-top: 40px; padding: 8px 0;
       border-bottom: 2px solid #e0e0e0; font-size: 21px; }
  h3 { color: #0f3460; margin-top: 25px; font-size: 17px; }
  h4 { color: #555; margin-top: 15px; }
  code { background: #f0f0f0; padding: 2px 6px; border-radius: 4px;
         font-family: "Cascadia Code", "Fira Code", Consolas, monospace; font-size: 13px; }
  pre { background: #1e1e2e; color: #cdd6f4; padding: 16px 20px;
        border-radius: 8px; overflow-x: auto; font-size: 13px;
        line-height: 1.5; border: 1px solid #313244; }
  pre code { background: none; padding: 0; color: inherit; }
  .note { background: #e8f5e9; border-left: 4px solid #4CAF50;
          padding: 12px 16px; margin: 15px 0; border-radius: 0 8px 8px 0; }
  .warning { background: #fff3e0; border-left: 4px solid #FF9800;
             padding: 12px 16px; margin: 15px 0; border-radius: 0 8px 8px 0; }
  .danger { background: #ffebee; border-left: 4px solid #f44336;
            padding: 12px 16px; margin: 15px 0; border-radius: 0 8px 8px 0; }
  .info { background: #e3f2fd; border-left: 4px solid #2196F3;
          padding: 12px 16px; margin: 15px 0; border-radius: 0 8px 8px 0; }
  table { border-collapse: collapse; margin: 15px 0; width: 100%; }
  th, td { border: 1px solid #ddd; padding: 10px 14px; text-align: left; }
  th { background: #2196F3; color: white; font-weight: 600; }
  tr:nth-child(even) { background: #f9f9f9; }
  .toc { background: #f5f5f5; padding: 20px; border-radius: 8px; margin: 20px 0; }
  .toc ul { list-style: none; padding-left: 0; }
  .toc ul ul { padding-left: 20px; }
  .toc a { text-decoration: none; color: #0f3460; }
  .toc a:hover { text-decoration: underline; }
  .step-num { display: inline-block; background: #2196F3; color: white;
              width: 28px; height: 28px; text-align: center; line-height: 28px;
              border-radius: 50%; font-weight: bold; margin-right: 8px; font-size: 14px; }
  .output-box { background: #263238; color: #aed581; padding: 12px 16px;
                border-radius: 8px; font-family: monospace; font-size: 12px;
                margin: 10px 0; white-space: pre-wrap; overflow-x: auto;
                border: 1px solid #37474f; }
  .gui-step { background: #fff; border: 1px solid #e0e0e0; border-radius: 8px;
              padding: 16px 20px; margin: 12px 0;
              box-shadow: 0 1px 4px rgba(0,0,0,0.06); }
  .gui-step h4 { margin-top: 0; color: #0f3460; }
  .badge-rec { display: inline-block; background: #2196F3; color: white;
               padding: 2px 8px; border-radius: 4px; font-size: 11px;
               font-weight: bold; margin-left: 6px; vertical-align: middle; }
  .file-tree { background: #f5f5f5; padding: 16px 20px; border-radius: 8px;
               font-family: monospace; font-size: 13px; margin: 10px 0;
               border: 1px solid #e0e0e0; line-height: 1.8; }
  .rstudio-tip { background: #f3e5f5; border-left: 4px solid #9C27B0;
                 padding: 12px 16px; margin: 15px 0; border-radius: 0 8px 8px 0; }
  footer { text-align: center; color: #999; font-size: 12px;
           margin-top: 50px; padding-top: 20px; border-top: 1px solid #eee; }
  @media print {
    body { background: white; }
    .gui-step { break-inside: avoid; }
    pre { white-space: pre-wrap; word-wrap: break-word; }
  }
</style>
"""


def build_html():
    princ_img = os.path.join(ASSETS_DIR, "cr_principle.png")
    real_data_dir = os.path.join(ASSETS_DIR, "real_data")
    fig_hela_grid = os.path.join(real_data_dir, "hela_heatmap_grid.png")
    fig_ddm2_grid = os.path.join(real_data_dir, "mb_ddm2_heatmap_grid.png")
    fig_ffhe2_grid = os.path.join(real_data_dir, "mb_ffhe2_heatmap_grid.png")
    fig_3sample_comp = os.path.join(real_data_dir, "3sample_heatmap_comparison.png")
    fig_cr_alif_comp = os.path.join(real_data_dir, "cr_alif_comparison.png")
    fig_summary_bar = os.path.join(real_data_dir, "sample_summary.png")
    fig_metrics_comp = os.path.join(real_data_dir, "metrics_comparison.png")

    html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>CR \u5206\u6790 - RStudio \u5feb\u901f\u5165\u95e8\u6307\u5357</title>
{CSS}
</head>
<body>

<h1>污染比 (CR) 分析
<span style="font-size:14px;color:#2196F3;font-weight:normal;margin-left:12px;">
RStudio + Windows 版</span></h1>
<p style="color:#666;font-size:14px;">
msProteomiX 多引擎蛋白质组学数据处理包 &mdash; timsTOF Pro 离子淌度-质谱 (IM-MS) 污染比分析模块<br>
本指南面向 <strong>Windows + RStudio</strong> 用户。
全部操作均可在 RStudio 内完成，无需命令行知识。
</p>

<!-- \u76ee\u5f55 -->
<div class="toc">
<strong>\u76ee\u5f55</strong>
<ul>
<li><a href="#sec1">1 \u6982\u8ff0</a></li>
<li><a href="#sec2">2 \u73af\u5883\u914d\u7f6e\uff08\u4ec5\u9700\u4e00\u6b21\uff0c\u7ea6 5 \u5206\u949f\uff09</a>
  <ul>
    <li><a href="#sec21">2.1 \u5b89\u88c5 R + RStudio</a></li>
    <li><a href="#sec22">2.2 \u5b89\u88c5 msProteomiX</a></li>
    <li><a href="#sec23">2.3 \u5b89\u88c5 Python + alphatims</a></li>
    <li><a href="#sec24">2.4 \u9a8c\u8bc1\u5b89\u88c5</a></li>
  </ul>
</li>
<li><a href="#sec3">3 \u5b8c\u6574\u5de5\u4f5c\u6d41\u7a0b\uff08\u5206\u6b65\u64cd\u4f5c\uff09</a>
  <ul>
    <li><a href="#sec31">3.1 \u521b\u5efa\u9879\u76ee</a></li>
    <li><a href="#sec32">3.2 \u4ece .d \u6587\u4ef6\u63d0\u53d6 CR \u6570\u636e</a></li>
    <li><a href="#sec33">3.3 \u4ea4\u4e92\u5f0f\u5206\u7ec4\u4e0e\u53ef\u89c6\u5316</a></li>
    <li><a href="#sec34">3.4 \u67e5\u770b\u7ed3\u679c</a></li>
  </ul>
</li>
<li><a href="#sec4">4 \u7ed3\u679c\u89e3\u8bfb</a>
  <ul>
    <li><a href="#sec41">4.1 IM-m/z \u5bc6\u5ea6\u56fe</a></li>
    <li><a href="#sec42">4.2 CR \u4e0e ALIF \u68af\u5ea6\u66f2\u7ebf</a></li>
    <li><a href="#sec43">4.3 \u8de8\u6837\u54c1\u5bf9\u6bd4</a></li>
    <li><a href="#sec44">4.4 \u6307\u6807\u9c81\u68d2\u6027\u5206\u6790</a></li>
  </ul>
</li>
<li><a href="#sec5">5 \u5e38\u89c1\u95ee\u9898</a></li>
<li><a href="#sec6">6 \u53c2\u6570\u53c2\u8003</a></li>
<li><a href="#secA">\u9644\u5f55\uff1aCR \u516c\u5f0f\u63a8\u5bfc</a></li>
</ul>
</div>

<!-- ================================================================ -->
<h2 id="sec1">1 \u6982\u8ff0</h2>

<h3>\u4ec0\u4e48\u662f CR\uff1f</h3>
<p>
CR\uff08Contamination Ratio\uff0c\u6c61\u67d3\u6bd4\uff09\u7528\u4e8e\u8bc4\u4f30\u5355\u7ec6\u80de/\u5fae\u91cf\u86cb\u767d\u8d28\u7ec4\u5b66\u6837\u54c1\u7684\u6e05\u6d01\u5ea6\u3002
\u5b83\u901a\u8fc7\u6d4b\u91cf timsTOF Pro \u6570\u636e IM-m/z \u7a7a\u95f4\u4e2d
<strong>\u5355\u7535\u8377\u6c61\u67d3\u7269</strong>\uff08\u805a\u5408\u7269\u3001\u8868\u9762\u6d3b\u6027\u5242\u3001\u5316\u5b66\u566a\u58f0\uff09\u4e0e
<strong>\u591a\u7535\u8377\u80bd\u6bb5\u4fe1\u53f7</strong>\u7684\u6bd4\u4f8b\u6765\u5b9e\u73b0\u3002
</p>

<div class="info">
<strong>核心原理：</strong>在捕获型离子淌度 (TIMS) 中，单电荷污染物的
1/K<sub>0</sub> 值显著高于同 m/z 的多电荷肽段。
利用一条分界线可将两者分开。
<br><br>
<strong>CR = &Sigma;I<sub>污染</sub> / &Sigma;I<sub>信号</sub> &times; 100%</strong>
&mdash; 数值越低越干净。
</div>

{img_tag(princ_img, "\u56fe 1\uff1aCR \u5206\u754c\u7ebf\u539f\u7406\u3002\u7ebf\u4e0a\u65b9\uff1a\u5355\u7535\u8377\u6c61\u67d3\u7269\uff1b\u7ebf\u4e0b\u65b9\uff1a\u591a\u7535\u8377\u80bd\u6bb5\u4fe1\u53f7\u3002")}

<h3>\u5de5\u4f5c\u6d41\u7a0b\u603b\u89c8</h3>
<p>\u5168\u90e8\u6d41\u7a0b\u5728 <strong>RStudio</strong> \u5185\u5b8c\u6210\uff1a</p>
<table>
<tr><th>\u6b65\u9aa4</th><th>\u5185\u5bb9</th><th>\u65f6\u95f4</th></tr>
<tr><td><span class="step-num">1</span></td>
    <td>\u521b\u5efa\u9879\u76ee\u5e76\u590d\u5236 .d \u6587\u4ef6</td><td>1 \u5206\u949f</td></tr>
<tr><td><span class="step-num">2</span></td>
    <td><code>extract_cr_data()</code> \u63d0\u53d6\u524d\u4f53\u79bb\u5b50\u5e76\u8ba1\u7b97 CR</td>
    <td>\u6bcf\u4e2a\u6587\u4ef6 2-5 \u5206\u949f</td></tr>
<tr><td><span class="step-num">3</span></td>
    <td>\u8fd0\u884c <code>33_contamination_ratio.R</code> &mdash; \u4ea4\u4e92\u5f0f\u5206\u7ec4 + \u81ea\u52a8\u7ed8\u56fe</td>
    <td>1 \u5206\u949f</td></tr>
<tr><td><span class="step-num">4</span></td>
    <td>\u5728 <code>output/</code> \u67e5\u770b HTML \u62a5\u544a</td><td>&mdash;</td></tr>
</table>

<!-- ================================================================ -->
<h2 id="sec2">2 \u73af\u5883\u914d\u7f6e\uff08\u4ec5\u9700\u4e00\u6b21\uff09</h2>

<p>\u4ee5\u4e0b\u6b65\u9aa4\u53ea\u9700\u5728 Windows \u7535\u8111\u4e0a\u6267\u884c\u4e00\u6b21\u3002</p>

<h3 id="sec21">2.1 \u5b89\u88c5 R + RStudio</h3>
<div class="gui-step">
<h4><span class="step-num">1</span> \u4e0b\u8f7d\u5e76\u5b89\u88c5</h4>
<ul>
<li><strong>R</strong> (4.3+)\uff1a<a href="https://cran.r-project.org/bin/windows/base/" target="_blank">
https://cran.r-project.org</a> &mdash; \u70b9\u51fb "Download R x.x.x for Windows"</li>
<li><strong>RStudio</strong>\uff1a<a href="https://posit.co/download/rstudio-desktop/" target="_blank">
https://posit.co/download/rstudio-desktop/</a> &mdash; \u4e0b\u8f7d Windows \u5b89\u88c5\u5305</li>
</ul>
<p>\u5168\u90e8\u4f7f\u7528\u9ed8\u8ba4\u8bbe\u7f6e\u5b89\u88c5\u5373\u53ef\u3002</p>
</div>

<h3 id="sec22">2.2 \u5b89\u88c5 msProteomiX</h3>
<div class="gui-step">
<h4><span class="step-num">2</span> \u6253\u5f00 RStudio\uff0c\u5728\u63a7\u5236\u53f0\u7c98\u8d34\u4ee5\u4e0b\u4ee3\u7801</h4>
<pre><code># \u5b89\u88c5 msProteomiX\uff08\u4e00\u884c\u5373\u53ef\uff09
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("BayOmics/msProteomiX")

# \u5b89\u88c5 reticulate\uff08\u7528\u4e8e R \u8c03\u7528 Python\uff09
install.packages("reticulate")</code></pre>

<div class="output-box"># \u9884\u671f\u8f93\u51fa\uff1a
* DONE (msProteomiX)
=== msProteomiX v0.3.0 ===</div>
</div>

<h3 id="sec23">2.3 \u5b89\u88c5 Python + alphatims</h3>
<div class="gui-step">
<h4><span class="step-num">3</span> \u65b9\u6848 A\uff1a\u4ece RStudio \u5185\u5b89\u88c5\uff08\u6700\u7b80\u5355\uff09</h4>
<pre><code># \u5728 R \u63a7\u5236\u53f0\u6267\u884c\uff1a
reticulate::install_miniconda()     # \u81ea\u52a8\u5b89\u88c5 Miniconda
reticulate::py_install("alphatims") # \u5b89\u88c5 alphatims</code></pre>

<div class="note">
<strong>\u63d0\u793a\uff1a</strong>\u8fd9\u4f1a\u5b89\u88c5\u4e00\u4e2a\u7531 R \u7ba1\u7406\u7684\u79c1\u6709 Python \u73af\u5883\u3002
\u4f60\u4e0d\u9700\u8981\u5355\u72ec\u5b89\u88c5 Anaconda\uff0c\u4e5f\u4e0d\u9700\u8981\u8f93\u5165\u4efb\u4f55 Python \u547d\u4ee4\u3002
</div>
</div>

<div class="gui-step">
<h4>\u65b9\u6848 B\uff1a\u4f7f\u7528\u5df2\u6709\u7684 Anaconda</h4>
<p>\u5982\u679c\u4f60\u7535\u8111\u4e0a\u5df2\u7ecf\u6709 Anaconda/Miniconda\uff0c\u53ea\u9700\u5728 Anaconda Prompt \u4e2d\u5b89\u88c5 alphatims\uff1a</p>
<pre><code># \u5728 Anaconda Prompt \u4e2d\u6267\u884c\uff08\u4e0d\u662f R\uff09\uff1a
pip install alphatims</code></pre>
<p>RStudio \u7684 <code>reticulate</code> \u4f1a\u81ea\u52a8\u68c0\u6d4b\u4f60\u7684 Anaconda \u5b89\u88c5\u3002</p>
</div>

<h3 id="sec24">2.4 \u9a8c\u8bc1\u5b89\u88c5</h3>
<div class="gui-step">
<h4><span class="step-num">4</span> \u5728 R \u63a7\u5236\u53f0\u8fd0\u884c</h4>
<pre><code>library(msProteomiX)    # \u5e94\u6253\u5370\u7248\u672c\u4fe1\u606f
library(reticulate)
py_config()             # \u5e94\u663e\u793a Python \u8def\u5f84
py_module_available("alphatims")  # \u5e94\u8fd4\u56de TRUE</code></pre>

<div class="output-box">=== msProteomiX v0.3.0 ===
[1] TRUE</div>

<div class="note">
\u5982\u679c <code>py_module_available("alphatims")</code> \u8fd4\u56de <code>FALSE</code>\uff0c
\u8bf7\u6267\u884c <code>reticulate::py_install("alphatims")</code>\uff0c\u7136\u540e\u91cd\u542f RStudio\u3002
</div>
</div>

<!-- ================================================================ -->
<h2 id="sec3">3 \u5b8c\u6574\u5de5\u4f5c\u6d41\u7a0b</h2>

<h3 id="sec31">3.1 \u521b\u5efa\u9879\u76ee</h3>

<div class="gui-step">
<h4><span class="step-num">1</span> \u521b\u5efa\u9879\u76ee\u6587\u4ef6\u5939</h4>
<pre><code>library(msProteomiX)
create_project("D:/Projects/CR_Analysis")
# \u8fd9\u4f1a\u521b\u5efa\uff1a
#   D:/Projects/CR_Analysis/
#     \u251c\u2500\u2500 scripts/      \uff08\u7528\u6237\u811a\u672c\uff0c\u5305\u542b 33_contamination_ratio.R\uff09
#     \u251c\u2500\u2500 wkdir/        \uff08\u6570\u636e\u76ee\u5f55\uff09
#     \u2514\u2500\u2500 output/       \uff08\u7ed3\u679c\u8f93\u51fa\uff09</code></pre>
</div>

<div class="gui-step">
<h4><span class="step-num">2</span> \u590d\u5236 .d \u6587\u4ef6</h4>
<p>\u5c06 timsTOF Pro \u539f\u59cb\u6570\u636e\u6587\u4ef6\u5939\uff08<code>.d</code>\uff09\u590d\u5236\u5230 <code>wkdir/</code>
\u6216\u4fdd\u7559\u5728\u539f\u59cb\u4f4d\u7f6e\uff08\u540e\u7eed\u6307\u5b9a\u8def\u5f84\u5373\u53ef\uff09\u3002</p>
<p>\u793a\u4f8b\u76ee\u5f55\u7ed3\u6784\uff1a</p>
<div class="file-tree">
D:/data/<br>
&nbsp;&nbsp;\u251c\u2500\u2500 Hela_10ng_Rep1.d/<br>
&nbsp;&nbsp;\u251c\u2500\u2500 Hela_10ng_Rep2.d/<br>
&nbsp;&nbsp;\u251c\u2500\u2500 MouseBrain_DDM_Rep1.d/<br>
&nbsp;&nbsp;\u251c\u2500\u2500 MouseBrain_DDM_Rep2.d/<br>
&nbsp;&nbsp;\u251c\u2500\u2500 MouseBrain_FFHE_Rep1.d/<br>
&nbsp;&nbsp;\u2514\u2500\u2500 MouseBrain_FFHE_Rep2.d/
</div>
</div>

<h3 id="sec32">3.2 \u4ece .d \u6587\u4ef6\u63d0\u53d6 CR \u6570\u636e</h3>

<div class="gui-step">
<h4><span class="step-num">3</span> \u5728 R \u63a7\u5236\u53f0\u6267\u884c</h4>
<pre><code>library(msProteomiX)

# \u5904\u7406\u67d0\u4e2a\u76ee\u5f55\u4e0b\u7684\u6240\u6709 .d \u6587\u4ef6\uff08\u81ea\u52a8\u53d1\u73b0\uff09
extract_cr_data(
  input_path   = "D:/data/",
  output_dir   = "wkdir/cr_output",
  n_slices     = 10,           # \u6bcf\u4e2a\u6587\u4ef6\u53d6 10 \u4e2a RT \u5207\u7247
  export_slices = TRUE          # \u5bfc\u51fa\u5bc6\u5ea6\u56fe\u6570\u636e
)</code></pre>

<div class="output-box">>>> Running extract_precursors.py ...
    Input:  D:/data/
    Output: D:/Projects/CR_Analysis/wkdir/cr_output
    Slices: 10
    Python: 3.10.12
    alphatims: 1.1.1
[INFO] Processing: Hela_10ng_Rep1.d
[INFO]   Loaded 2325916050 detector events
[INFO]   Slice 1/10: CR=6.4%, TIC(log2)=25.5
...
>>> Done! 6 summary file(s) in: wkdir/cr_output</div>

<div class="rstudio-tip">
<strong>RStudio \u63d0\u793a\uff1a</strong>\u6bcf\u4e2a .d \u6587\u4ef6\u5904\u7406\u9700\u8981 2-5 \u5206\u949f\uff08\u53d6\u51b3\u4e8e\u6587\u4ef6\u5927\u5c0f\uff09\u3002
\u63a7\u5236\u53f0\u4f1a\u663e\u793a\u8fdb\u5ea6\u4fe1\u606f\uff0c\u8bf7\u52ff\u5173\u95ed RStudio\u3002
</div>
</div>

<h3 id="sec33">3.3 \u4ea4\u4e92\u5f0f\u5206\u7ec4\u4e0e\u53ef\u89c6\u5316</h3>

<div class="gui-step">
<h4><span class="step-num">4</span> \u6253\u5f00\u5e76\u8fd0\u884c\u811a\u672c</h4>
<p>\u5728 RStudio \u4e2d\u6253\u5f00 <code>scripts/33_contamination_ratio.R</code>\uff0c\u70b9\u51fb\u53f3\u4e0a\u89d2
<strong>Source</strong> \u6309\u94ae\u3002</p>

<p>\u811a\u672c\u4f1a\u81ea\u52a8\u6267\u884c\u4ee5\u4e0b\u64cd\u4f5c\uff1a</p>
<ol>
<li><strong>\u81ea\u52a8\u53d1\u73b0</strong> <code>wkdir/cr_output/</code> \u4e2d\u7684\u6240\u6709\u7ed3\u679c\u6587\u4ef6</li>
<li><strong>\u5217\u51fa\u6837\u54c1</strong>\uff1a
<div class="output-box">=== Step 1: Auto-discovering CR results ===
>>> Found 6 sample(s):
    [1] Hela_10ng_Rep1
    [2] Hela_10ng_Rep2
    [3] MouseBrain_DDM_Rep1
    [4] MouseBrain_DDM_Rep2
    [5] MouseBrain_FFHE_Rep1
    [6] MouseBrain_FFHE_Rep2</div>
</li>
<li><strong>\u4ea4\u4e92\u5f0f\u5206\u7ec4\u5411\u5bfc</strong> &mdash; \u6309\u63d0\u793a\u8f93\u5165\u7f16\u53f7\u548c\u7ec4\u540d\uff1a
<div class="output-box">=== Step 2: Sample grouping ===
==========================================================
  Samples to group:
----------------------------------------------------------
  [1] Hela_10ng_Rep1        [2] Hela_10ng_Rep2
  [3] MouseBrain_DDM_Rep1   [4] MouseBrain_DDM_Rep2
  [5] MouseBrain_FFHE_Rep1  [6] MouseBrain_FFHE_Rep2
==========================================================
>>> \u8f93\u5165\u6837\u54c1\u7f16\u53f7 (e.g. 1-2,5): 1-2
>>> \u8f93\u5165\u7ec4\u540d: HeLa
>>> \u8f93\u5165\u6837\u54c1\u7f16\u53f7 (e.g. 1-2,5): 3-4
>>> \u8f93\u5165\u7ec4\u540d: DDM
>>> \u8f93\u5165\u6837\u54c1\u7f16\u53f7 (e.g. 1-2,5): 5-6
>>> \u8f93\u5165\u7ec4\u540d: FFHE
>>> \u5206\u7ec4\u5b8c\u6210? (y/n): y</div>
</li>
<li><strong>\u81ea\u52a8\u751f\u6210</strong>\u6240\u6709\u56fe\u8868\u548c HTML \u62a5\u544a</li>
</ol>

<div class="note">
<strong>\u7075\u6d3b\u5206\u7ec4\uff1a</strong>\u652f\u6301<strong>\u4efb\u610f\u6570\u91cf</strong>\u7684\u7ec4\uff08\u4e0d\u9650\u4e8e 3 \u7ec4\uff09\uff0c\u7ec4\u540d\u53ef\u81ea\u7531\u5b9a\u4e49\u3002
\u5206\u7ec4\u4fe1\u606f\u4fdd\u5b58\u5230 <code>cr_group_info.csv</code>\uff0c\u4e0b\u6b21\u8fd0\u884c\u65f6\u81ea\u52a8\u52a0\u8f7d\u3002
</div>
</div>

<h3 id="sec34">3.4 \u67e5\u770b\u7ed3\u679c</h3>

<div class="gui-step">
<h4><span class="step-num">5</span> \u68c0\u67e5 output/ \u6587\u4ef6\u5939</h4>
<div class="file-tree">
output/<br>
&nbsp;&nbsp;\u251c\u2500\u2500 CR_Analysis_CR_Gradient.png &nbsp;&nbsp; \u2190 CR vs RT \u66f2\u7ebf<br>
&nbsp;&nbsp;\u251c\u2500\u2500 CR_Analysis_CR_TIC.png &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; \u2190 \u603b\u79bb\u5b50\u6d41\u5f3a\u5ea6<br>
&nbsp;&nbsp;\u251c\u2500\u2500 CR_Analysis_CR_Heatmap_*.png &nbsp; \u2190 IM-m/z \u5bc6\u5ea6\u56fe<br>
&nbsp;&nbsp;\u2514\u2500\u2500 CR_Analysis_CR_Report.html &nbsp;&nbsp; \u2190 \u6c47\u603b HTML \u62a5\u544a
</div>
<p>\u53cc\u51fb <code>*_CR_Report.html</code> \u5373\u53ef\u5728\u6d4f\u89c8\u5668\u4e2d\u67e5\u770b\u3002</p>
</div>

<!-- ================================================================ -->
<h2 id="sec4">4 \u7ed3\u679c\u89e3\u8bfb</h2>

<h3 id="sec41">4.1 IM-m/z \u5bc6\u5ea6\u56fe</h3>
<p>\u6bcf\u5f20\u5bc6\u5ea6\u56fe\u5c55\u793a\u4e00\u4e2a RT \u5207\u7247\u7684\u79bb\u5b50\u5206\u5e03\u3002
<strong>\u767d\u8272\u5206\u754c\u7ebf</strong>\u5c06\u6c61\u67d3\u7269\uff08\u4e0a\u65b9\uff09\u4e0e\u80bd\u6bb5\u4fe1\u53f7\uff08\u4e0b\u65b9\uff09\u5206\u5f00\u3002
\u6bcf\u4e2a\u9762\u677f\u6807\u6ce8\u4e86 CR \u503c\u3002</p>

<h4>4.1.1 HeLa 10ng\uff08\u6e05\u6d01\u53c2\u8003\u6837\u54c1\uff09</h4>
{img_tag(fig_hela_grid, "\u56fe 2a\uff1aHeLa 10ng IM-m/z \u5bc6\u5ea6\u56fe\u3002CR \u8303\u56f4 6.4%\uff08\u68af\u5ea6\u4e2d\u6bb5\uff0c\u6700\u5e72\u51c0\uff09\u81f3 61.9%\uff08\u68af\u5ea6\u672b\u7aef\uff09\u3002")}

<h4>4.1.2 \u5c0f\u9f20\u8111\u7ec4\u7ec7 DDM2\uff08\u9ad8\u6c61\u67d3\uff09</h4>
{img_tag(fig_ddm2_grid, "\u56fe 2b\uff1a\u5c0f\u9f20\u8111 DDM2 \u5bc6\u5ea6\u56fe\u3002CR \u8303\u56f4 8.2% \u81f3 140.1%\u3002\u5206\u754c\u7ebf\u4e0a\u65b9\u53ef\u89c1\u660e\u663e\u7684\u6c61\u67d3\u7269\u6761\u5e26\u3002")}

<h4>4.1.3 \u5c0f\u9f20\u8111\u7ec4\u7ec7 FFHE2\uff08FFPE \u7ec4\u7ec7\uff09</h4>
{img_tag(fig_ffhe2_grid, "\u56fe 2c\uff1a\u5c0f\u9f20\u8111 FFHE2 \u5bc6\u5ea6\u56fe\u3002\u4e2d\u7b49\u6c61\u67d3\uff0c\u68af\u5ea6\u5c3e\u90e8\u5bcc\u96c6\u3002")}

<h3 id="sec42">4.2 CR \u4e0e ALIF \u68af\u5ea6\u66f2\u7ebf</h3>
<p>\u5c55\u793a\u5404\u6837\u54c1\u7c7b\u578b\u6c61\u67d3\u7a0b\u5ea6\u968f\u6db2\u76f8\u68af\u5ea6\u7684\u53d8\u5316\uff1a</p>
{img_tag(fig_cr_alif_comp, "\u56fe 3\uff1aCR\uff08\u4e0a\uff09\u548c ALIF\uff08\u4e0b\uff09\u968f\u68af\u5ea6\u53d8\u5316\u3002MB DDM2\uff08\u84dd\u8272\uff09\u5168\u7a0b\u6c61\u67d3\u504f\u9ad8\u3002")}

<div class="info">
<strong>\u5982\u4f55\u89e3\u8bfb\uff1a</strong>
<ul>
<li><strong>U \u5f62\u66f2\u7ebf</strong>\uff08HeLa\uff09\uff1a\u6b63\u5e38 &mdash; \u68af\u5ea6\u5f00\u59cb/\u7ed3\u675f\u65f6\u80bd\u6bb5\u4fe1\u53f7\u8f83\u5c11</li>
<li><strong>\u6574\u4f53\u504f\u9ad8</strong>\uff08DDM2\uff09\uff1a\u7cfb\u7edf\u6027\u6c61\u67d3 &mdash; \u6837\u54c1\u5236\u5907\u95ee\u9898</li>
<li><strong>CR &gt; 100%</strong>\uff1a\u6c61\u67d3\u7269\u4fe1\u53f7\u8d85\u8fc7\u80bd\u6bb5\u4fe1\u53f7\uff08\u4e25\u91cd\u6c61\u67d3\uff09</li>
</ul>
</div>

<h3 id="sec43">4.3 \u8de8\u6837\u54c1\u5bf9\u6bd4</h3>
<h4>3&times;3 \u5bc6\u5ea6\u56fe\u77e9\u9635</h4>
{img_tag(fig_3sample_comp, "\u56fe 4\uff1a3 \u4e2a\u6837\u54c1 \u00d7 3 \u4e2a RT \u5207\u7247\u3002\u7edf\u4e00\u8272\u6807\u3002MB DDM2 \u6301\u7eed\u9ad8\u6c61\u67d3\u3002")}

<h4>\u6837\u54c1\u6c34\u5e73\u6c47\u603b</h4>
{img_tag(fig_summary_bar, "\u56fe 5\uff1a\u5e73\u5747 CR \u548c ALIF\uff08\u542b\u6807\u51c6\u5dee\u8bef\u5dee\u7ebf\uff09\u3002")}

<h3 id="sec44">4.4 \u6307\u6807\u9c81\u68d2\u6027\u5206\u6790</h3>
{img_tag(fig_metrics_comp, "\u56fe 6\uff1a5 \u79cd\u5149\u8c31\u6e05\u6d01\u5ea6\u6307\u6807\u5bf9\u6bd4\u3002CR \u548c ALIF \u9ad8\u5ea6\u76f8\u5173 (r=0.99)\uff0c\u4f18\u4e8e\u5176\u4ed6\u66ff\u4ee3\u6307\u6807\u3002")}

<div class="note">
<strong>\u7ed3\u8bba\uff1a</strong>CR \u662f\u6700\u4f73\u5355\u4e00\u6307\u6807\uff0c\u52a8\u6001\u8303\u56f4\u8fbe 10 \u500d\u3002
ALIF\uff08\u6709\u754c [0, 100%]\uff09\u63a8\u8350\u4f5c\u4e3a\u8865\u5145\u6307\u6807\u3002
</div>

<!-- ================================================================ -->
<h2 id="sec5">5 \u5e38\u89c1\u95ee\u9898</h2>

<table>
<tr><th>\u95ee\u9898</th><th>\u539f\u56e0</th><th>\u89e3\u51b3\u65b9\u6cd5</th></tr>
<tr><td><code>reticulate::py_config()</code> \u65e0 Python</td>
    <td>\u672a\u5b89\u88c5 Python</td>
    <td>\u5728 R \u4e2d\u6267\u884c <code>reticulate::install_miniconda()</code></td></tr>
<tr><td><code>py_module_available("alphatims")</code> \u4e3a FALSE</td>
    <td>alphatims \u672a\u5b89\u88c5</td>
    <td>\u6267\u884c <code>reticulate::py_install("alphatims")</code></td></tr>
<tr><td><code>extract_cr_data()</code> \u65e0\u54cd\u5e94</td>
    <td>\u5927\u6587\u4ef6\uff08>10GB\uff09</td>
    <td>\u6b63\u5e38\u73b0\u8c61 &mdash; \u8bf7\u7b49\u5f85 5-10 \u5206\u949f\uff0c\u67e5\u770b\u63a7\u5236\u53f0\u8fdb\u5ea6</td></tr>
<tr><td>CR \u503c\u5168\u90e8\u5f88\u9ad8 (&gt;50%)</td>
    <td>\u6837\u54c1\u6c61\u67d3\u6216\u5206\u754c\u7ebf\u4e0d\u5f53</td>
    <td>\u68c0\u67e5\u6837\u54c1\u5236\u5907\uff1b\u8c03\u6574 <code>LINE_POINT1/2</code></td></tr>
<tr><td>\u627e\u4e0d\u5230\u5bc6\u5ea6\u56fe CSV</td>
    <td>\u63d0\u53d6\u65f6\u672a\u5bfc\u51fa\u5207\u7247</td>
    <td>\u91cd\u65b0\u8fd0\u884c <code>extract_cr_data()</code>\uff0c\u8bbe\u7f6e <code>export_slices = TRUE</code></td></tr>
<tr><td><code>Error: CR output directory not found</code></td>
    <td>\u8def\u5f84\u4e0d\u5339\u914d</td>
    <td>\u68c0\u67e5\u811a\u672c\u4e2d <code>CR_DIR</code> \u662f\u5426\u4e0e\u63d0\u53d6\u65f6\u7684 <code>output_dir</code> \u4e00\u81f4</td></tr>
<tr><td>RStudio \u63d0\u793a\u201c\u627e\u4e0d\u5230\u51fd\u6570\u201d</td>
    <td>\u672a\u52a0\u8f7d\u5305</td>
    <td>\u5148\u6267\u884c <code>library(msProteomiX)</code></td></tr>
</table>

<!-- ================================================================ -->
<h2 id="sec6">6 \u53c2\u6570\u53c2\u8003</h2>

<h3><code>extract_cr_data()</code> \u51fd\u6570\u53c2\u6570</h3>
<table>
<tr><th>\u53c2\u6570</th><th>\u9ed8\u8ba4\u503c</th><th>\u8bf4\u660e</th></tr>
<tr><td><code>input_path</code></td><td>\uff08\u5fc5\u586b\uff09</td>
    <td>.d \u6587\u4ef6\u5939\u3001.hdf \u6587\u4ef6\u6216\u5305\u542b\u5b83\u4eec\u7684\u76ee\u5f55</td></tr>
<tr><td><code>output_dir</code></td><td><code>"wkdir/cr_output"</code></td>
    <td>\u7ed3\u679c\u4fdd\u5b58\u8def\u5f84</td></tr>
<tr><td><code>n_slices</code></td><td><code>10</code></td>
    <td>\u7b49\u95f4\u8ddd RT \u5207\u7247\u6570</td></tr>
<tr><td><code>export_slices</code></td><td><code>TRUE</code></td>
    <td>\u662f\u5426\u5bfc\u51fa\u6bcf\u4e2a\u5207\u7247\u7684 CSV \u6587\u4ef6</td></tr>
<tr><td><code>heatmap_rts</code></td><td><code>NULL</code></td>
    <td>\u6307\u5b9a\u7528\u4e8e\u5bc6\u5ea6\u56fe\u7684 RT \u503c\uff08\u5206\u949f\uff09\uff0cNULL = \u5168\u90e8</td></tr>
<tr><td><code>line_points</code></td><td><code>"350,0.8,950,1.3"</code></td>
    <td>CR \u5206\u754c\u7ebf: mz1,im1,mz2,im2</td></tr>
<tr><td><code>python</code></td><td><code>NULL</code></td>
    <td>Python \u8def\u5f84\uff08\u81ea\u52a8\u68c0\u6d4b\uff09</td></tr>
</table>

<h3>\u811a\u672c <code>33_contamination_ratio.R</code> \u914d\u7f6e</h3>
<table>
<tr><th>\u53d8\u91cf</th><th>\u9ed8\u8ba4\u503c</th><th>\u8bf4\u660e</th></tr>
<tr><td><code>PROJECT_NAME</code></td><td><code>"CR_Analysis"</code></td>
    <td>\u8f93\u51fa\u6587\u4ef6\u540d\u524d\u7f00</td></tr>
<tr><td><code>CR_DIR</code></td><td><code>"wkdir/cr_output"</code></td>
    <td>\u5305\u542b <code>*_cr_summary.csv</code> \u7684\u76ee\u5f55</td></tr>
<tr><td><code>OUTPUT</code></td><td><code>"output"</code></td>
    <td>\u56fe\u8868\u548c\u62a5\u544a\u8f93\u51fa\u76ee\u5f55</td></tr>
<tr><td><code>LINE_POINT1</code></td><td><code>c(350, 0.8)</code></td>
    <td>\u5206\u754c\u7ebf\u7aef\u70b9 1: (m/z, 1/K<sub>0</sub>)</td></tr>
<tr><td><code>LINE_POINT2</code></td><td><code>c(950, 1.3)</code></td>
    <td>\u5206\u754c\u7ebf\u7aef\u70b9 2: (m/z, 1/K<sub>0</sub>)</td></tr>
</table>

<!-- ================================================================ -->
<h2 id="secA">\u9644\u5f55\uff1aCR \u516c\u5f0f\u63a8\u5bfc</h2>

<h3>\u7269\u7406\u57fa\u7840</h3>
<p>\u5728 TIMS \u4e2d\uff0c1/K<sub>0</sub> &prop; CCS / z\u3002\u5728\u76f8\u540c m/z \u4e0b\uff0c
\u5355\u7535\u8377\u79bb\u5b50 (z=1) \u7684 1/K<sub>0</sub> \u6bd4\u53cc\u7535\u8377\u79bb\u5b50 (z=2) \u9ad8\u7ea6 26%\u3002</p>

<h3>\u5206\u754c\u7ebf</h3>
<p style="text-align:center;">
\u7ecf\u8fc7 (350, 0.8) \u548c (950, 1.3) \u4e24\u70b9<br>
<strong>1/K<sub>0</sub> = 0.000833 &times; m/z + 0.5083</strong>
</p>

<h3>CR \u8ba1\u7b97</h3>
<ul>
<li>1/K<sub>0</sub> &ge; \u5206\u754c\u7ebf\u503c &rarr; <strong>\u6c61\u67d3\u7269</strong>\uff08\u5355\u7535\u8377\uff09</li>
<li>1/K<sub>0</sub> &lt; \u5206\u754c\u7ebf\u503c &rarr; <strong>\u4fe1\u53f7</strong>\uff08\u591a\u7535\u8377\uff09</li>
</ul>
<p style="text-align:center;font-size:16px;">
<strong>CR = &Sigma;I<sub>\u6c61\u67d3</sub> / &Sigma;I<sub>\u4fe1\u53f7</sub> &times; 100%</strong>
</p>

<footer>
msProteomiX v0.3.0 &mdash; BayOmics<br>
RStudio + Windows \u4e2d\u6587\u7248<br>
Generated by generate_cr_manual_rstudio_zh.py
</footer>

</body>
</html>"""

    with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Generated: {OUTPUT_HTML}")
    print(f"File size: {os.path.getsize(OUTPUT_HTML) / 1024:.1f} KB")


if __name__ == "__main__":
    build_html()
