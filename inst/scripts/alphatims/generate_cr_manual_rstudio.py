#!/usr/bin/env python3
"""
generate_cr_manual_rstudio.py — Build RStudio + Windows version of CR Manual.

This version targets users who work primarily in RStudio on Windows.
The original CR_Manual.html (advanced/cross-platform) is kept intact.

Usage:
    python generate_cr_manual_rstudio.py

Output:
    inst/scripts/alphatims/CR_Manual_RStudio.html
"""
import base64
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ASSETS_DIR = os.path.join(
    os.path.dirname(SCRIPT_DIR), "..", "..", "wkdir", "cr_manual_assets"
)
if not os.path.isdir(ASSETS_DIR):
    ASSETS_DIR = os.path.join(SCRIPT_DIR, "assets")

OUTPUT_HTML = os.path.join(SCRIPT_DIR, "CR_Manual_RStudio.html")


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
        return f'<p style="color:#999;">[Image not found: {os.path.basename(path)}]</p>'
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
    font-family: "Segoe UI", Roboto, "Helvetica Neue", sans-serif;
    max-width: 900px; margin: 0 auto; padding: 30px 20px;
    background: #fafbfc; color: #333; line-height: 1.7;
  }
  h1 { color: #1a1a2e; border-bottom: 3px solid #2196F3; padding-bottom: 12px;
       font-size: 28px; }
  h2 { color: #16213e; margin-top: 40px; padding: 8px 0;
       border-bottom: 2px solid #e0e0e0; font-size: 22px; }
  h3 { color: #0f3460; margin-top: 25px; font-size: 18px; }
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
  footer { text-align: center; color: #999; font-size: 12px;
           margin-top: 50px; padding-top: 20px; border-top: 1px solid #eee; }
  .rstudio-tip { background: #f3e5f5; border-left: 4px solid #9C27B0;
                 padding: 12px 16px; margin: 15px 0; border-radius: 0 8px 8px 0; }
</style>
"""


def build_html():
    # --- Image paths ---
    princ_img = os.path.join(ASSETS_DIR, "cr_principle.png")
    arch_img = os.path.join(ASSETS_DIR, "architecture.png")

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
<title>CR Analysis - RStudio Quick Start Guide</title>
{CSS}
</head>
<body>

<h1>Contamination Ratio (CR) Analysis
<span style="font-size:14px;color:#2196F3;font-weight:normal;margin-left:12px;">
RStudio + Windows Edition</span></h1>
<p style="color:#666;font-size:14px;">
msProteomiX &mdash; timsTOF Pro IM-MS Contamination Evaluation<br>
This guide is designed for <strong>Windows + RStudio</strong> users.
All operations can be completed within RStudio &mdash; no command-line knowledge required.
</p>

<!-- Table of Contents -->
<div class="toc">
<strong>Table of Contents</strong>
<ul>
<li><a href="#sec1">1 Overview</a></li>
<li><a href="#sec2">2 One-Time Setup (5 min)</a>
  <ul>
    <li><a href="#sec21">2.1 Install R + RStudio</a></li>
    <li><a href="#sec22">2.2 Install msProteomiX</a></li>
    <li><a href="#sec23">2.3 Install Python + alphatims</a></li>
    <li><a href="#sec24">2.4 Verify Everything</a></li>
  </ul>
</li>
<li><a href="#sec3">3 Complete Workflow (Step by Step)</a>
  <ul>
    <li><a href="#sec31">3.1 Create Project</a></li>
    <li><a href="#sec32">3.2 Extract CR Data from .d Files</a></li>
    <li><a href="#sec33">3.3 Interactive Grouping &amp; Visualization</a></li>
    <li><a href="#sec34">3.4 View Results</a></li>
  </ul>
</li>
<li><a href="#sec4">4 Understanding the Results</a>
  <ul>
    <li><a href="#sec41">4.1 IM-m/z Density Maps</a></li>
    <li><a href="#sec42">4.2 CR &amp; ALIF Gradient Profiles</a></li>
    <li><a href="#sec43">4.3 Cross-Sample Comparison</a></li>
    <li><a href="#sec44">4.4 Metric Robustness</a></li>
  </ul>
</li>
<li><a href="#sec5">5 Troubleshooting</a></li>
<li><a href="#sec6">6 Parameter Reference</a></li>
<li><a href="#secA">Appendix: CR Formula</a></li>
</ul>
</div>

<!-- ================================================================ -->
<h2 id="sec1">1 Overview</h2>

<h3>What is CR?</h3>
<p>
CR (Contamination Ratio) evaluates sample cleanliness by measuring
<strong>singly-charged contaminants</strong> vs <strong>multi-charged peptide signals</strong>
in the IM-m/z space of timsTOF Pro data.
</p>

<div class="info">
<strong>Key Insight:</strong> In trapped ion mobility, singly-charged contaminants
(polymers, detergents, noise) sit <strong>above</strong> the peptide corridor
in the IM-m/z plane. A dividing line separates these two populations.
<br><br>
<strong>CR = &Sigma;I<sub>above</sub> / &Sigma;I<sub>below</sub> &times; 100%</strong>
&mdash; lower is cleaner.
</div>

{img_tag(princ_img, "Fig. 1: CR dividing line in IM-m/z space. Above: singly-charged contaminants. Below: multi-charge peptide signal.")}

<h3>Workflow Overview</h3>
<p>The entire workflow runs in <strong>RStudio</strong>:</p>
<table>
<tr><th>Step</th><th>What Happens</th><th>Time</th></tr>
<tr><td><span class="step-num">1</span></td>
    <td>Create project &amp; copy .d files</td><td>1 min</td></tr>
<tr><td><span class="step-num">2</span></td>
    <td><code>extract_cr_data()</code> extracts precursors &amp; computes CR</td>
    <td>2-5 min/file</td></tr>
<tr><td><span class="step-num">3</span></td>
    <td>Run <code>33_contamination_ratio.R</code> &mdash; interactive grouping + plots</td>
    <td>1 min</td></tr>
<tr><td><span class="step-num">4</span></td>
    <td>View HTML report in <code>output/</code></td><td>&mdash;</td></tr>
</table>

<!-- ================================================================ -->
<h2 id="sec2">2 One-Time Setup (5 min)</h2>

<p>You only need to do this once on your Windows computer.</p>

<h3 id="sec21">2.1 Install R + RStudio</h3>
<div class="gui-step">
<h4><span class="step-num">1</span> Download and Install</h4>
<ul>
<li><strong>R</strong> (4.3+): <a href="https://cran.r-project.org/bin/windows/base/" target="_blank">
https://cran.r-project.org</a> &mdash; click "Download R x.x.x for Windows"</li>
<li><strong>RStudio</strong>: <a href="https://posit.co/download/rstudio-desktop/" target="_blank">
https://posit.co/download/rstudio-desktop/</a> &mdash; download the Windows installer</li>
</ul>
<p>Install both with default settings.</p>
</div>

<h3 id="sec22">2.2 Install msProteomiX</h3>
<div class="gui-step">
<h4><span class="step-num">2</span> Open RStudio, paste in Console</h4>
<pre><code># Install msProteomiX (one-liner)
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("BayOmics/msProteomiX")

# Install reticulate for Python integration
install.packages("reticulate")</code></pre>

<div class="output-box"># Expected output:
* DONE (msProteomiX)
=== msProteomiX v0.3.0 ===</div>
</div>

<h3 id="sec23">2.3 Install Python + alphatims</h3>
<div class="gui-step">
<h4><span class="step-num">3</span> Option A: Install from RStudio (easiest)</h4>
<pre><code># In R console:
reticulate::install_miniconda()    # Install Miniconda automatically
reticulate::py_install("alphatims") # Install alphatims</code></pre>

<div class="note">
<strong>Tip:</strong> This installs a private Python environment managed by R.
You don't need to install Anaconda or know any Python commands.
</div>
</div>

<div class="gui-step">
<h4>Option B: Use existing Anaconda (if you already have it)</h4>
<p>If you already have Anaconda/Miniconda installed, just install alphatims:</p>
<pre><code># In Anaconda Prompt (not R):
pip install alphatims</code></pre>
<p>RStudio's <code>reticulate</code> will auto-detect your Anaconda installation.</p>
</div>

<h3 id="sec24">2.4 Verify Everything</h3>
<div class="gui-step">
<h4><span class="step-num">4</span> Run in R Console</h4>
<pre><code>library(msProteomiX)    # Should print version banner
library(reticulate)
py_config()             # Should show Python path
py_module_available("alphatims")  # Should return TRUE</code></pre>

<div class="output-box">=== msProteomiX v0.3.0 ===
Multi-Engine Mass Spectrometry Proteomics Data Processing
[1] TRUE</div>

<div class="note">
If <code>py_module_available("alphatims")</code> returns <code>FALSE</code>,
run <code>reticulate::py_install("alphatims")</code> and restart RStudio.
</div>
</div>

<!-- ================================================================ -->
<h2 id="sec3">3 Complete Workflow (Step by Step)</h2>

<h3 id="sec31">3.1 Create Project</h3>

<div class="gui-step">
<h4><span class="step-num">1</span> Create a project folder</h4>
<pre><code>library(msProteomiX)
create_project("D:/Projects/CR_Analysis")
# This creates:
#   D:/Projects/CR_Analysis/
#     ├── scripts/      (user scripts including 33_contamination_ratio.R)
#     ├── wkdir/        (put your data here)
#     └── output/       (results go here)</code></pre>
</div>

<div class="gui-step">
<h4><span class="step-num">2</span> Copy your .d files</h4>
<p>Copy your timsTOF <code>.d</code> folders into <code>wkdir/</code>
or keep them in their original location (you'll specify the path later).</p>
<p>Example directory structure:</p>
<div class="file-tree">
D:/data/<br>
&nbsp;&nbsp;&#9500;&#9472; Hela_10ng_Rep1.d/<br>
&nbsp;&nbsp;&#9500;&#9472; Hela_10ng_Rep2.d/<br>
&nbsp;&nbsp;&#9500;&#9472; MouseBrain_DDM_Rep1.d/<br>
&nbsp;&nbsp;&#9500;&#9472; MouseBrain_DDM_Rep2.d/<br>
&nbsp;&nbsp;&#9500;&#9472; MouseBrain_FFHE_Rep1.d/<br>
&nbsp;&nbsp;&#9492;&#9472; MouseBrain_FFHE_Rep2.d/
</div>
</div>

<h3 id="sec32">3.2 Extract CR Data from .d Files</h3>

<div class="gui-step">
<h4><span class="step-num">3</span> Run extraction in R Console</h4>
<pre><code>library(msProteomiX)

# Process all .d files in a directory (auto-discovers files)
extract_cr_data(
  input_path   = "D:/data/",
  output_dir   = "wkdir/cr_output",
  n_slices     = 10,           # 10 RT slices per file
  export_slices = TRUE          # Export CSVs for density maps
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
<strong>RStudio Tip:</strong> This step may take 2-5 minutes per .d file
(depending on file size). You'll see progress messages in the Console.
Don't close RStudio while it's running.
</div>
</div>

<h3 id="sec33">3.3 Interactive Grouping &amp; Visualization</h3>

<div class="gui-step">
<h4><span class="step-num">4</span> Open and Source the script</h4>
<p>In RStudio, open <code>scripts/33_contamination_ratio.R</code>, then click
the <strong>Source</strong> button (top-right of the editor).</p>

<p>The script will:</p>
<ol>
<li><strong>Auto-discover</strong> all CR result files in <code>wkdir/cr_output/</code></li>
<li><strong>List your samples</strong> in the Console:
<div class="output-box">=== Step 1: Auto-discovering CR results ===
>>> Found 6 sample(s):
    [1] Hela_10ng_Rep1
    [2] Hela_10ng_Rep2
    [3] MouseBrain_DDM_Rep1
    [4] MouseBrain_DDM_Rep2
    [5] MouseBrain_FFHE_Rep1
    [6] MouseBrain_FFHE_Rep2</div>
</li>
<li><strong>Interactive grouping wizard</strong> &mdash; you assign groups at the prompt:
<div class="output-box">=== Step 2: Sample grouping ===
==========================================================
  Samples to group:
----------------------------------------------------------
  [1] Hela_10ng_Rep1        [2] Hela_10ng_Rep2
  [3] MouseBrain_DDM_Rep1   [4] MouseBrain_DDM_Rep2
  [5] MouseBrain_FFHE_Rep1  [6] MouseBrain_FFHE_Rep2
==========================================================
>>> Enter sample numbers (e.g. 1-2,5): 1-2
>>> Enter group name: HeLa
>>> Enter sample numbers (e.g. 1-2,5): 3-4
>>> Enter group name: DDM
>>> Enter sample numbers (e.g. 1-2,5): 5-6
>>> Enter group name: FFHE
>>> Grouping complete? (y/n): y</div>
</li>
<li><strong>Auto-generate</strong> all plots and HTML report</li>
</ol>

<div class="note">
<strong>Flexible:</strong> You can define <strong>any number of groups</strong>
with any names. The grouping is saved to <code>cr_group_info.csv</code> and
reused on subsequent runs.
</div>
</div>

<h3 id="sec34">3.4 View Results</h3>

<div class="gui-step">
<h4><span class="step-num">5</span> Check the output/ folder</h4>
<div class="file-tree">
output/<br>
&nbsp;&nbsp;&#9500;&#9472; CR_Analysis_CR_Gradient.png &nbsp;&nbsp; &larr; CR vs RT curve<br>
&nbsp;&nbsp;&#9500;&#9472; CR_Analysis_CR_TIC.png &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &larr; TIC intensity<br>
&nbsp;&nbsp;&#9500;&#9472; CR_Analysis_CR_Heatmap_*.png &nbsp; &larr; IM-m/z density maps<br>
&nbsp;&nbsp;&#9492;&#9472; CR_Analysis_CR_Report.html &nbsp;&nbsp; &larr; Combined HTML report
</div>
<p>Double-click the <code>*_CR_Report.html</code> file to view in your browser.</p>
</div>

<!-- ================================================================ -->
<h2 id="sec4">4 Understanding the Results</h2>

<h3 id="sec41">4.1 IM-m/z Density Maps</h3>
<p>Each density map shows ion distribution in the IM-m/z plane for one RT slice.
The <strong>white dividing line</strong> separates contaminants (above) from peptide signal (below).
CR is annotated on each panel.</p>

<h4>4.1.1 HeLa 10ng (clean reference)</h4>
{img_tag(fig_hela_grid, "Fig. 2a: HeLa 10ng IM-m/z density maps across LC gradient. CR ranges from 6.4% (mid-gradient, cleanest) to 61.9% (late gradient).")}

<h4>4.1.2 Mouse Brain DDM2 (high contamination)</h4>
{img_tag(fig_ddm2_grid, "Fig. 2b: Mouse Brain DDM2 density maps. CR ranges from 8.2% to 140.1%. Prominent contaminant bands above the dividing line.")}

<h4>4.1.3 Mouse Brain FFHE2 (FFPE tissue)</h4>
{img_tag(fig_ffhe2_grid, "Fig. 2c: Mouse Brain FFHE2 density maps. Moderate contamination with tail-gradient enrichment.")}

<h3 id="sec42">4.2 CR &amp; ALIF Gradient Profiles</h3>
<p>These plots show how contamination changes across the LC gradient for each sample type.</p>
{img_tag(fig_cr_alif_comp, "Fig. 3: CR (upper) and ALIF (lower) profiles across LC gradient. MB DDM2 (blue) shows elevated contamination throughout.")}

<div class="info">
<strong>Interpreting the profiles:</strong>
<ul>
<li><strong>U-shaped curve</strong> (HeLa): Normal &mdash; gradient start/end have less peptide signal</li>
<li><strong>Elevated baseline</strong> (DDM2): Systematic contamination &mdash; sample prep issue</li>
<li><strong>CR &gt; 100%</strong>: Contaminant signal exceeds peptide signal (severely dirty)</li>
</ul>
</div>

<h3 id="sec43">4.3 Cross-Sample Comparison</h3>
<h4>3&times;3 Heatmap Matrix</h4>
{img_tag(fig_3sample_comp, "Fig. 4: 3 samples x 3 RT slices (early, mid, late). Same color scale. MB DDM2 shows consistently higher contamination.")}

<h4>Summary Statistics</h4>
{img_tag(fig_summary_bar, "Fig. 5: Mean CR and ALIF with SD error bars.")}

<h3 id="sec44">4.4 Metric Robustness</h3>
{img_tag(fig_metrics_comp, "Fig. 6: Comparison of 5 spectral cleanliness metrics. CR and ALIF are highly correlated (r=0.99) and outperform alternatives.")}

<div class="note">
<strong>Conclusion:</strong> CR is the best single indicator with 10&times; dynamic range.
ALIF (bounded [0, 100%]) is recommended as a complementary metric.
</div>

<!-- ================================================================ -->
<h2 id="sec5">5 Troubleshooting</h2>

<table>
<tr><th>Problem</th><th>Cause</th><th>Solution</th></tr>
<tr><td><code>reticulate::py_config()</code> shows no Python</td>
    <td>No Python installed</td>
    <td>Run <code>reticulate::install_miniconda()</code> in R</td></tr>
<tr><td><code>py_module_available("alphatims")</code> is FALSE</td>
    <td>alphatims not installed</td>
    <td>Run <code>reticulate::py_install("alphatims")</code></td></tr>
<tr><td><code>extract_cr_data()</code> hangs</td>
    <td>Large .d file (>10GB)</td>
    <td>Normal &mdash; wait 5-10 min. Check Console for progress.</td></tr>
<tr><td>CR values all very high (&gt;50%)</td>
    <td>Possible contamination or wrong dividing line</td>
    <td>Check sample prep; try adjusting <code>LINE_POINT1/2</code></td></tr>
<tr><td>No heatmap slice CSVs found</td>
    <td><code>export_slices = FALSE</code> in extract step</td>
    <td>Re-run <code>extract_cr_data()</code> with <code>export_slices = TRUE</code></td></tr>
<tr><td><code>Error: CR output directory not found</code></td>
    <td>Path mismatch</td>
    <td>Check <code>CR_DIR</code> in script matches <code>output_dir</code> in extract</td></tr>
<tr><td>RStudio says "function not found"</td>
    <td>Package not loaded</td>
    <td>Run <code>library(msProteomiX)</code> first</td></tr>
</table>

<!-- ================================================================ -->
<h2 id="sec6">6 Parameter Reference</h2>

<h3><code>extract_cr_data()</code></h3>
<table>
<tr><th>Parameter</th><th>Default</th><th>Description</th></tr>
<tr><td><code>input_path</code></td><td>(required)</td>
    <td>Path to .d folder, .hdf file, or a directory containing them</td></tr>
<tr><td><code>output_dir</code></td><td><code>"wkdir/cr_output"</code></td>
    <td>Where to save results</td></tr>
<tr><td><code>n_slices</code></td><td><code>10</code></td>
    <td>Number of equally-spaced RT slices</td></tr>
<tr><td><code>export_slices</code></td><td><code>TRUE</code></td>
    <td>Export per-slice CSV files for density maps</td></tr>
<tr><td><code>heatmap_rts</code></td><td><code>NULL</code></td>
    <td>Specific RT values (min) for heatmap. NULL = all slices.</td></tr>
<tr><td><code>line_points</code></td><td><code>"350,0.8,950,1.3"</code></td>
    <td>CR dividing line: mz1,im1,mz2,im2</td></tr>
<tr><td><code>python</code></td><td><code>NULL</code></td>
    <td>Python path (auto-detected via reticulate)</td></tr>
</table>

<h3>Script <code>33_contamination_ratio.R</code></h3>
<table>
<tr><th>Variable</th><th>Default</th><th>Description</th></tr>
<tr><td><code>PROJECT_NAME</code></td><td><code>"CR_Analysis"</code></td>
    <td>Prefix for output filenames</td></tr>
<tr><td><code>CR_DIR</code></td><td><code>"wkdir/cr_output"</code></td>
    <td>Directory with <code>*_cr_summary.csv</code> files</td></tr>
<tr><td><code>OUTPUT</code></td><td><code>"output"</code></td>
    <td>Output directory for plots and report</td></tr>
<tr><td><code>LINE_POINT1</code></td><td><code>c(350, 0.8)</code></td>
    <td>Dividing line point 1: (m/z, 1/K<sub>0</sub>)</td></tr>
<tr><td><code>LINE_POINT2</code></td><td><code>c(950, 1.3)</code></td>
    <td>Dividing line point 2: (m/z, 1/K<sub>0</sub>)</td></tr>
</table>

<!-- ================================================================ -->
<h2 id="secA">Appendix: CR Formula</h2>

<h3>Physical Basis</h3>
<p>In TIMS, 1/K<sub>0</sub> &prop; CCS / z. At the same m/z, singly-charged
ions (z=1) have ~26% higher 1/K<sub>0</sub> than doubly-charged ions (z=2).</p>

<h3>Dividing Line</h3>
<p style="text-align:center;">
(350, 0.8) &rarr; (950, 1.3)<br>
<strong>1/K<sub>0</sub> = 0.000833 &times; m/z + 0.5083</strong>
</p>

<h3>CR Computation</h3>
<ul>
<li>1/K<sub>0</sub> &ge; line value &rarr; <strong>Contaminant</strong></li>
<li>1/K<sub>0</sub> &lt; line value &rarr; <strong>Signal</strong></li>
</ul>
<p style="text-align:center;font-size:16px;">
<strong>CR = &Sigma;I<sub>contam</sub> / &Sigma;I<sub>signal</sub> &times; 100%</strong>
</p>

<footer>
msProteomiX v0.3.0 &mdash; BayOmics<br>
RStudio + Windows Edition<br>
Generated by generate_cr_manual_rstudio.py
</footer>

</body>
</html>"""

    with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Generated: {OUTPUT_HTML}")
    print(f"File size: {os.path.getsize(OUTPUT_HTML) / 1024:.1f} KB")


if __name__ == "__main__":
    build_html()
