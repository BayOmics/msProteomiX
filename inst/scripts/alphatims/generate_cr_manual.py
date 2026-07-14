#!/usr/bin/env python3
"""
generate_cr_manual.py — Build the self-contained HTML CR operation manual.

All images are embedded as base64. Run this script to regenerate CR_Manual.html.

Usage:
    python generate_cr_manual.py

Output:
    inst/scripts/alphatims/CR_Manual.html
"""
import base64
import os
import sys

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ASSETS_DIR = os.path.join(os.path.dirname(SCRIPT_DIR), "..", "..", "wkdir", "cr_manual_assets")
# Also check relative to script for installed package
if not os.path.isdir(ASSETS_DIR):
    ASSETS_DIR = os.path.join(SCRIPT_DIR, "assets")

OUTPUT_HTML = os.path.join(SCRIPT_DIR, "CR_Manual.html")


def img_base64(path):
    """Read an image file and return a base64 data URI string."""
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
    """Create an HTML img tag with caption."""
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
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
    max-width: 900px; margin: 0 auto; padding: 30px 20px;
    background: #fafbfc; color: #333; line-height: 1.7;
  }
  h1 { color: #1a1a2e; border-bottom: 3px solid #4CAF50; padding-bottom: 12px;
       font-size: 28px; }
  h2 { color: #16213e; margin-top: 40px; padding: 8px 0;
       border-bottom: 2px solid #e0e0e0; font-size: 22px; }
  h3 { color: #0f3460; margin-top: 25px; font-size: 18px; }
  h4 { color: #555; margin-top: 15px; }
  code { background: #f0f0f0; padding: 2px 6px; border-radius: 4px;
         font-family: "SF Mono", "Fira Code", Consolas, monospace; font-size: 13px; }
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
  th { background: #4CAF50; color: white; font-weight: 600; }
  tr:nth-child(even) { background: #f9f9f9; }
  .toc { background: #f5f5f5; padding: 20px; border-radius: 8px;
         margin: 20px 0; }
  .toc ul { list-style: none; padding-left: 0; }
  .toc ul ul { padding-left: 20px; }
  .toc a { text-decoration: none; color: #0f3460; }
  .toc a:hover { text-decoration: underline; }
  .step-num { display: inline-block; background: #4CAF50; color: white;
              width: 28px; height: 28px; text-align: center; line-height: 28px;
              border-radius: 50%; font-weight: bold; margin-right: 8px;
              font-size: 14px; }
  .output-box { background: #263238; color: #aed581; padding: 12px 16px;
                border-radius: 8px; font-family: monospace; font-size: 12px;
                margin: 10px 0; white-space: pre-wrap; overflow-x: auto;
                border: 1px solid #37474f; }
  .gui-step { background: #fff; border: 1px solid #e0e0e0; border-radius: 8px;
              padding: 16px 20px; margin: 12px 0;
              box-shadow: 0 1px 4px rgba(0,0,0,0.06); }
  .gui-step h4 { margin-top: 0; color: #0f3460; }
  .badge-rec { display: inline-block; background: #4CAF50; color: white;
               padding: 2px 8px; border-radius: 4px; font-size: 11px;
               font-weight: bold; margin-left: 6px; vertical-align: middle; }
  .file-tree { background: #f5f5f5; padding: 16px 20px; border-radius: 8px;
               font-family: monospace; font-size: 13px; margin: 10px 0;
               border: 1px solid #e0e0e0; line-height: 1.8; }
  footer { text-align: center; color: #999; font-size: 12px;
           margin-top: 50px; padding-top: 20px; border-top: 1px solid #eee; }
</style>
"""


def build_html():
    """Build the complete HTML manual."""

    # Resolve image paths
    arch_img = os.path.join(ASSETS_DIR, "architecture.png")
    princ_img = os.path.join(ASSETS_DIR, "cr_principle.png")

    # E2E test output images (fallback if available)
    plots_dir = os.path.join(
        os.path.dirname(SCRIPT_DIR), "..", "..",
        "wkdir", "cr_output", "plots"
    )

    # Real data figures
    real_data_dir = os.path.join(ASSETS_DIR, "real_data")
    fig2f_img = os.path.join(real_data_dir, "hela10ng_cr_summary.png")
    fig2g_img = os.path.join(real_data_dir, "hela10ng_heatmap_grid.png")
    # 3-sample comparison figures
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
<title>Contamination Ratio (CR) Analysis - Operation Manual</title>
{CSS}
</head>
<body>

<h1>Contamination Ratio (CR) Analysis</h1>
<p style="color:#666;font-size:14px;">
msProteomiX &mdash; timsTOF Pro IM-MS Contamination Evaluation<br>
Reference: "Multimodal single cell-resolved spatial proteomics reveal pancreatic tumor heterogeneity" (Fig. 2f, 2g)
</p>

<!-- Table of Contents -->
<div class="toc">
<strong>Table of Contents</strong>
<ul>
<li><a href="#sec1">1 Overview</a></li>
<li><a href="#sec2">2 Environment Setup</a>
  <ul>
    <li><a href="#sec21">2.1 Install Python (Anaconda)</a></li>
    <li><a href="#sec22">2.2 Install alphatims</a></li>
    <li><a href="#sec23">2.3 Install R + msProteomiX</a></li>
    <li><a href="#sec24">2.4 Verify Installation</a></li>
  </ul>
</li>
<li><a href="#sec3">3 Data Preparation</a>
  <ul>
    <li><a href="#sec31">3.1 timsTOF .d Files</a></li>
    <li><a href="#sec32">3.2 Experiment Design</a></li>
    <li><a href="#sec33">3.3 AlphaTims GUI (.d &rarr; .hdf, Windows) <span style='color:#4CAF50;font-size:11px;'>&star; Recommended</span></a></li>
    <li><a href="#sec34">3.4 CLI Conversion (Alternative)</a></li>
    <li><a href="#sec35">3.5 Transfer to macOS</a></li>
  </ul>
</li>
<li><a href="#sec4">4 Extract Precursors &amp; Compute CR</a>
  <ul>
    <li><a href="#sec41">4.1 From RStudio (Recommended) &star;</a></li>
    <li><a href="#sec42">4.2 From Command Line (Advanced)</a></li>
    <li><a href="#sec43">4.3 R Function Reference</a></li>
    <li><a href="#sec44">4.4 CLI Parameter Reference</a></li>
    <li><a href="#sec45">4.5 Output Structure</a></li>
  </ul>
</li>
<li><a href="#sec5">5 Visualization (R)</a>
  <ul>
    <li><a href="#sec51">5.1 Import CR Results</a></li>
    <li><a href="#sec52">5.2 Group Assignment</a></li>
    <li><a href="#sec53">5.3 CR Gradient (Fig. 2f upper)</a></li>
    <li><a href="#sec54">5.4 TIC Intensity (Fig. 2f lower)</a></li>
    <li><a href="#sec55">5.5 IM-m/z Heatmap (Fig. 2g)</a></li>
    <li><a href="#sec56">5.6 HTML Report</a></li>
  </ul>
</li>
<li><a href="#sec6">6 Advanced</a></li>
<li><a href="#secA">Appendix A: Parameter Reference</a></li>
<li><a href="#secB">Appendix B: Troubleshooting</a></li>
<li><a href="#secC">Appendix C: CR Formula Derivation</a></li>
</ul>
</div>

<!-- ================================================================ -->
<h2 id="sec1">1 Overview</h2>

<h3>1.1 What is Contamination Ratio (CR)?</h3>
<p>
CR (Contamination Ratio) evaluates single-cell / micro-proteomics sample cleanliness
by measuring the proportion of <strong>singly-charged contaminants</strong> relative to
<strong>multi-charged peptide signals</strong> in the IM-m/z (ion mobility vs mass-to-charge)
space of timsTOF Pro data.
</p>

<div class="info">
<strong>Key Insight:</strong> In trapped ion mobility spectrometry, singly-charged contaminants
(polymers, detergents, chemical noise) have <strong>higher 1/K<sub>0</sub></strong> than
multi-charged peptides at the same m/z, because 1/K<sub>0</sub> = CCS/z and z=1 yields a higher ratio.
A dividing line in the IM-m/z plane separates these two populations.
</div>

<h3>1.2 Principle</h3>
{img_tag(princ_img, "Fig. 1: CR dividing line in IM-m/z space. Above the line: singly-charged contaminants. Below: multi-charge peptide signal.")}

<p>The CR is computed as:</p>
<p style="text-align:center;font-size:18px;">
<strong>CR = &Sigma;I<sub>above</sub> / &Sigma;I<sub>below</sub> &times; 100%</strong>
</p>
<p>where I<sub>above</sub> is the total intensity of ions above the dividing line (contaminants)
and I<sub>below</sub> is the total intensity of ions below (signal).</p>

<h3>1.3 Pipeline Architecture</h3>
{img_tag(arch_img, "Fig. 2: CR analysis pipeline architecture. Python (alphatims) extracts data, R (msProteomiX) visualizes.")}

<!-- ================================================================ -->
<h2 id="sec2">2 Environment Setup</h2>

<h3 id="sec21">2.1 Install Python</h3>
<p><span class="step-num">1</span> Download and install
<a href="https://www.anaconda.com/download" target="_blank">Anaconda</a> (Python 3.10+).</p>

<h3 id="sec22">2.2 Install alphatims</h3>
<p><span class="step-num">2</span> Open terminal (or Anaconda Prompt on Windows) and run:</p>
<pre><code>pip install alphatims</code></pre>

<div class="output-box">Successfully installed alphatims-1.1.1 alphabase-1.9.0 alpharaw-0.6.1 ...</div>

<div class="note">
<strong>Alternative (from R):</strong> If you prefer to install alphatims from within RStudio:
<pre><code>reticulate::py_install("alphatims")</code></pre>
This uses the Python environment managed by <code>reticulate</code>.
</div>

<h3 id="sec23">2.3 Install R + msProteomiX</h3>
<p><span class="step-num">3</span> Install <a href="https://cran.r-project.org/" target="_blank">R</a>
(4.3+) and <a href="https://posit.co/download/rstudio-desktop/" target="_blank">RStudio</a>.</p>
<p><span class="step-num">4</span> Install msProteomiX:</p>
<pre><code># In R console (one-liner, no devtools needed)
if (!requireNamespace("remotes")) install.packages("remotes")
remotes::install_github("BayOmics/msProteomiX")

# Also install reticulate for Python integration
install.packages("reticulate")</code></pre>

<h3 id="sec24">2.4 Verify Installation</h3>
<p><span class="step-num">5</span> Verify both tools:</p>
<pre><code># Python
python -c "import alphatims; print(alphatims.__version__)"

# R
Rscript -e "library(msProteomiX); cat('OK\\n')"</code></pre>

<div class="output-box">1.1.1
=== msProteomiX v0.3.0 ===
OK</div>

<!-- ================================================================ -->
<h2 id="sec3">3 Data Preparation</h2>

<h3 id="sec31">3.1 timsTOF Pro .d Files</h3>
<p>Prepare your timsTOF Pro raw data files (<code>.d</code> folders).
Each <code>.d</code> folder contains <code>analysis.tdf</code> and <code>analysis.tdf_bin</code>.</p>

<div class="note">
<strong>Naming Convention:</strong> Use descriptive sample names, e.g.:
<code>Hela_10ng_Rep1.d</code>, <code>MouseBrain_DDM_Rep2.d</code>
</div>

<h3 id="sec32">3.2 Experiment Design</h3>
<p>Recommended design for CR evaluation:</p>
<table>
<tr><th>Group</th><th>Description</th><th>Replicates</th></tr>
<tr><td>HeLa</td><td>Standard reference (clean sample)</td><td>&ge;3</td></tr>
<tr><td>W/o Ext. Wash</td><td>Target tissue, standard wash</td><td>&ge;3</td></tr>
<tr><td>W/ Ext. Wash</td><td>Target tissue, extended wash</td><td>&ge;3</td></tr>
</table>

<div class="danger">
<strong>Important for macOS Users:</strong> On macOS, the Bruker native DLL is not available.
AlphaTims will <strong>estimate</strong> m/z and mobility values from metadata, which can introduce
significant errors (&gt;2x in m/z). This makes CR calculation unreliable on macOS with .d files.
<br><br>
<strong>Solution:</strong> Convert <code>.d</code> to <code>.hdf</code> on <strong>Windows</strong> first,
then transfer <code>.hdf</code> files to macOS for analysis.
</div>

<h3 id="sec33">3.3 AlphaTims GUI Conversion (Windows) <span class="badge-rec">Recommended</span></h3>

<p>AlphaTims provides a <strong>stand-alone graphical user interface (GUI)</strong> that requires
no Python/CLI knowledge. This is the recommended method for Windows users.</p>

<div class="gui-step">
<h4><span class="step-num">1</span> Download &amp; Install AlphaTims GUI</h4>
<p>Download the latest <strong>one-click installer</strong> from
<a href="https://github.com/MannLabs/alphatims/releases/latest" target="_blank">
https://github.com/MannLabs/alphatims/releases/latest</a></p>
<ul>
<li>Choose <code>alphatims-X.Y.Z-windows-amd64.exe</code></li>
<li>Double-click the installer. If you see "Windows protected your PC", click
<strong>More info</strong> &rarr; <strong>Run anyway</strong></li>
<li>Select <strong>"Install for me only"</strong> (recommended, no admin needed)</li>
<li>Follow the wizard and finish installation</li>
</ul>
<div class="note">
<strong>Tip:</strong> No Python or Anaconda installation is required for the GUI version.
It is a fully self-contained application.
</div>
</div>

<div class="gui-step">
<h4><span class="step-num">2</span> Launch AlphaTims GUI</h4>
<p>Double-click the <strong>AlphaTims</strong> desktop icon, or find it in the Start menu.
The GUI opens in your default web browser (Chrome or Firefox recommended).</p>
<p>You will see the main interface with:</p>
<ul>
<li><strong>Upload File</strong> &mdash; input field for the .d folder path</li>
<li><strong>UPLOAD</strong> button &mdash; load the dataset</li>
<li><strong>Download GUI manual / Test data / Citation</strong> buttons</li>
</ul>
</div>

<div class="gui-step">
<h4><span class="step-num">3</span> Load a .d File</h4>
<p>In the "Upload File" text box, enter the full path to your <code>.d</code> folder:</p>
<pre><code>D:\\data\\Hela_10ng_Rep1.d</code></pre>
<p>Click <strong>UPLOAD</strong>. A loading spinner will appear while the data is being indexed.
For a typical 80-min run (~2 billion detector events), loading takes about 1&ndash;3 minutes.</p>
<div class="info">
<strong>Tip:</strong> You can also drag the <code>.d</code> folder path from Windows Explorer
into the text box.
</div>
</div>

<div class="gui-step">
<h4><span class="step-num">4</span> Explore the Data (Optional)</h4>
<p>After loading, AlphaTims displays an interactive IM-m/z heatmap.
You can adjust the <strong>Parameters</strong> panel:</p>
<ul>
<li><strong>Frame / RT</strong> slider &mdash; scroll through retention time</li>
<li><strong>Scan / IM</strong> slider &mdash; filter by ion mobility range</li>
<li><strong>TOF / m/z</strong> slider &mdash; filter by m/z range</li>
<li><strong>Show MS1 ions / MS2 ions</strong> checkboxes &mdash; toggle MS1 precursors</li>
<li><strong>Select axis for plots</strong> &mdash; switch between heatmap, spectrum, XIC, mobilogram</li>
</ul>
<p>This is useful for visually inspecting the IM-m/z distribution and confirming
the presence of singly-charged contaminant bands.</p>
</div>

<div class="gui-step">
<h4><span class="step-num">5</span> Save as HDF &mdash; Key Step!</h4>
<p>Scroll down to the <strong>"Export data"</strong> card. Find the <strong>HDF</strong> section:</p>
<ol>
<li>In the <strong>"Specify a path to save all data as a portable .hdf file"</strong> text box,
enter the output path, e.g.:<br>
<code>D:\\data\\Hela_10ng_Rep1</code> (the <code>.hdf</code> extension is added automatically)</li>
<li>Optionally check <strong>compress</strong> if you need smaller files
(saves ~50% space but 3&ndash;6x slower to load)</li>
<li>Click <strong>SAVE AS HDF</strong></li>
<li>Wait for the spinner to complete. A success message will appear.</li>
</ol>
<div class="warning">
<strong>Important:</strong> The HDF file preserves Bruker-calibrated m/z and mobility values.
This is essential for accurate CR calculation on macOS.
HDF files are roughly <strong>2x the size</strong> of the original .d folder.
</div>
</div>

<div class="gui-step">
<h4><span class="step-num">6</span> Repeat for All Samples</h4>
<p>For each <code>.d</code> file in your experiment:</p>
<ol>
<li>Type the new path in the <strong>Upload File</strong> box</li>
<li>Click <strong>UPLOAD</strong> to load the new dataset</li>
<li>Click <strong>SAVE AS HDF</strong> to export</li>
</ol>
<p>All <code>.hdf</code> files will be saved in the same directory as the original <code>.d</code> folders.</p>
</div>

<h3 id="sec34">3.4 CLI Conversion (Alternative)</h3>
<p>If you prefer command-line tools, alphatims also provides a CLI for batch conversion:</p>

<h4>3.4.1 Single File (Windows PowerShell / CMD)</h4>
<pre><code>alphatims export hdf "D:\\data\\sample.d"
# Output: D:\\data\\sample.hdf</code></pre>

<h4>3.4.2 Batch Conversion (Windows PowerShell)</h4>
<pre><code>Get-ChildItem -Path "D:\\data\\" -Filter "*.d" -Directory | ForEach-Object {{
    alphatims export hdf $_.FullName
    Write-Host "Converted: $($_.Name)"
}}</code></pre>

<div class="note">
The CLI requires a Python environment with <code>pip install alphatims</code>.
For batch processing of many files, the CLI is more efficient than the GUI.
</div>

<h3 id="sec35">3.5 Transfer to macOS</h3>
<p>Copy the <code>.hdf</code> files from Windows to your Mac via:</p>
<ul>
<li><strong>USB drive</strong> &mdash; fastest for large files</li>
<li><strong>Network share / SMB</strong> &mdash; convenient for shared lab storage</li>
<li><strong>Cloud storage</strong> (OneDrive, Google Drive) &mdash; for remote access</li>
</ul>
<p>Then use the <code>.hdf</code> files as input in all subsequent steps.</p>
<div class="info">
<strong>Note:</strong> If you run the full pipeline on Windows, no transfer is needed.
The Python <code>extract_precursors.py</code> and R visualization both work directly on Windows.
</div>

<!-- ================================================================ -->
<h2 id="sec4">4 Extract Precursors &amp; Compute CR</h2>

<h3 id="sec41">4.1 From RStudio (Recommended) &star;</h3>
<p>msProteomiX provides <code>extract_cr_data()</code> &mdash; a one-line R function
that calls the Python extraction script automatically. <strong>No command-line
knowledge needed.</strong></p>

<pre><code>library(msProteomiX)

# Single file (.d on Windows, .hdf on macOS)
extract_cr_data("D:/data/sample.d")

# Batch: all .d/.hdf files in a directory
extract_cr_data("D:/data/")

# Custom output directory and RT slices
extract_cr_data("sample.hdf",
                output_dir = "wkdir/cr_output",
                n_slices = 10,
                heatmap_rts = c(20, 40, 60))</code></pre>

<div class="output-box">>>> Running extract_precursors.py ...
    Input:  D:/data/sample.d
    Output: D:/data/cr_output
    Slices: 10
    Python: 3.10.12
    alphatims: 1.1.1
13:47:32 [INFO] Loaded in 138.0s: 2325916050 detector events
13:47:33 [INFO] Selected 10 slices...
>>> Done! 1 summary file(s) in: D:/data/cr_output</div>

<div class="note">
<strong>How it works:</strong> <code>extract_cr_data()</code> uses the <code>reticulate</code>
package to find your Python installation and verify that <code>alphatims</code>
is available. The actual extraction runs via the bundled
<code>extract_precursors.py</code> script &mdash; no manual Python command needed.
</div>

<div class="danger">
<strong>macOS Users:</strong> On macOS, <code>.d</code> files cannot be read accurately
(Bruker DLL unavailable). The function will warn you automatically.
Convert <code>.d</code> &rarr; <code>.hdf</code> on <strong>Windows</strong> first
(see Section 3.3), then use <code>.hdf</code> files on macOS.
</div>

<h3 id="sec42">4.2 From Command Line (Advanced)</h3>
<p>For power users or batch automation, the Python CLI is also available:</p>

<h4>4.2.1 Single File</h4>
<pre><code>python extract_precursors.py "D:\\data\\sample.d" --output_dir ./cr_output</code></pre>

<h4>4.2.2 Batch Processing</h4>
<pre><code>python extract_precursors.py "D:\\data\\" --output_dir ./cr_output</code></pre>

<h4>4.2.3 Export Heatmap Data</h4>
<pre><code>python extract_precursors.py sample.d \\
    --output_dir ./cr_output \\
    --export_slices \\
    --heatmap_rts 32,48</code></pre>

<h3 id="sec43">4.3 R Function Reference</h3>
<table>
<tr><th>Parameter</th><th>Default</th><th>Description</th></tr>
<tr><td><code>input_path</code></td><td>(required)</td><td><code>.d</code> folder, <code>.hdf</code> file, or directory</td></tr>
<tr><td><code>output_dir</code></td><td><code>"wkdir/cr_output"</code></td><td>Output directory</td></tr>
<tr><td><code>n_slices</code></td><td><code>10</code></td><td>Number of equally-spaced RT slices</td></tr>
<tr><td><code>export_slices</code></td><td><code>TRUE</code></td><td>Export per-slice CSV files</td></tr>
<tr><td><code>heatmap_rts</code></td><td><code>NULL</code></td><td>Numeric vector of RT values (min)</td></tr>
<tr><td><code>line_points</code></td><td><code>"350,0.8,950,1.3"</code></td><td>CR dividing line: mz1,im1,mz2,im2</td></tr>
<tr><td><code>python</code></td><td><code>NULL</code></td><td>Path to Python executable (auto-detect)</td></tr>
</table>

<h3 id="sec44">4.4 CLI Parameter Reference</h3>
<table>
<tr><th>Parameter</th><th>Default</th><th>Description</th></tr>
<tr><td><code>input</code></td><td>(required)</td><td><code>.d</code> folder, <code>.hdf</code> file, or directory</td></tr>
<tr><td><code>--output_dir</code></td><td><code>./cr_output</code></td><td>Output directory</td></tr>
<tr><td><code>--n_slices</code></td><td><code>10</code></td><td>Number of equally-spaced RT slices</td></tr>
<tr><td><code>--line_points</code></td><td><code>350,0.8,950,1.3</code></td><td>CR dividing line: mz1,im1,mz2,im2</td></tr>
<tr><td><code>--export_slices</code></td><td><code>False</code></td><td>Export per-slice CSV files</td></tr>
<tr><td><code>--heatmap_rts</code></td><td>None</td><td>RT values for heatmap (e.g. <code>32,48</code>)</td></tr>
</table>

<h3 id="sec45">4.5 Output Structure</h3>
<div class="file-tree">
cr_output/<br>
&nbsp;&nbsp;&#9500;&#9472; sample1_cr_summary.csv &nbsp;&nbsp;&nbsp; &larr; CR summary per slice<br>
&nbsp;&nbsp;&#9500;&#9472; sample2_cr_summary.csv<br>
&nbsp;&nbsp;&#9500;&#9472; all_cr_summary.csv &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; &larr; Combined (batch mode)<br>
&nbsp;&nbsp;&#9500;&#9472; sample1/<br>
&nbsp;&nbsp;&#9474;&nbsp;&nbsp;&#9500;&#9472; sample1_slice_rt0.0.csv &nbsp; &larr; Per-slice data<br>
&nbsp;&nbsp;&#9474;&nbsp;&nbsp;&#9500;&#9472; sample1_slice_rt8.9.csv<br>
&nbsp;&nbsp;&#9474;&nbsp;&nbsp;&#9492;&#9472; ...<br>
&nbsp;&nbsp;&#9492;&#9472; sample2/<br>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&#9492;&#9472; ...
</div>

<!-- ================================================================ -->
<h2 id="sec5">5 Visualization (R)</h2>

<p>Open RStudio and run the user script <code>33_contamination_ratio.R</code>,
or use the functions below interactively.</p>

<h3 id="sec51">5.1 Import CR Results</h3>
<pre><code>library(msProteomiX)

cr_data &lt;- import_cr_results(
  cr_dir    = "wkdir/cr_output",
  group_map = list(
    "HeLa"          = c("Hela"),
    "W/o Ext. Wash" = c("DDM"),
    "W/ Ext. Wash"  = c("FFHE")
  )
)</code></pre>

<div class="output-box">&gt;&gt;&gt; Loading CR results: 9 files from wkdir/cr_output
    9 samples, 3 groups, 90 data points</div>

<h3 id="sec52">5.2 Group Assignment</h3>
<p>The <code>group_map</code> parameter maps sample name patterns to group labels.
Each value is a character vector of patterns to match (case-insensitive):</p>
<pre><code>group_map = list(
  "GroupLabel" = c("pattern1", "pattern2"),
  "Another"   = c("pattern3")
)</code></pre>

<div class="note">
Samples not matching any pattern are assigned to group <code>"Other"</code>.
</div>

<h3 id="sec53">5.3 Real Data: IM-m/z Heatmaps</h3>
<p>The following figures are generated from <strong>real timsTOF Pro data</strong> (3 sample types,
10 RT slices each, ~500K-1.5M events per slice).</p>

<h4>5.3.1 HeLa 10ng (Picoliter frit, clean reference)</h4>
{img_tag(fig_hela_grid, "Fig. 3a: HeLa 10ng IM-m/z heatmaps across LC gradient. CR ranges from 6.4% (mid-gradient) to 61.9% (late gradient). Clear peptide bands visible in clean slices.")}

<h4>5.3.2 Mouse Brain DDM2 (high contamination)</h4>
{img_tag(fig_ddm2_grid, "Fig. 3b: Mouse Brain DDM2 IM-m/z heatmaps. CR ranges from 8.2% to 140.1%. Prominent vertical stripe artifacts visible at high m/z. Polymer contamination dominates mid-gradient.")}

<h4>5.3.3 Mouse Brain FFHE2 (FFPE tissue)</h4>
{img_tag(fig_ffhe2_grid, "Fig. 3c: Mouse Brain FFHE2 IM-m/z heatmaps. CR ranges from 8.8% to 57.8%. Moderate contamination with tail-gradient enrichment.")}

<h3 id="sec54">5.4 Cross-Sample Comparison</h3>
<h4>5.4.1 3-Sample Heatmap Matrix (3 &times; 3)</h4>
{img_tag(fig_3sample_comp, "Fig. 4: 3 samples (HeLa, MB DDM2, MB FFHE2) &times; 3 RT slices (early, mid, late). Same color scale. CR and ALIF annotated per panel. MB DDM2 shows consistently higher contamination.")}

<h4>5.4.2 CR &amp; ALIF across LC Gradient</h4>
{img_tag(fig_cr_alif_comp, "Fig. 5: CR (upper) and ALIF (lower) profiles across LC gradient for 3 samples. MB DDM2 (blue) shows elevated contamination throughout. HeLa (green) is cleanest mid-gradient. ALIF is bounded [0, 100%] unlike CR which can exceed 100%.")}

<h4>5.4.3 Sample-Level Summary</h4>
{img_tag(fig_summary_bar, "Fig. 6: Mean CR and ALIF with SD error bars. MB DDM2 has the highest contamination (mean CR ~54%). HeLa and FFHE2 are comparable (~24-27%).")}

<h3 id="sec55">5.5 Metric Robustness Analysis</h3>
<p>We compared 5 candidate spectral cleanliness metrics on real data.</p>
{img_tag(fig_metrics_comp, "Fig. 7: Comparison of CR, ALIF, Spectral Entropy, SCI, and BSI. Correlation matrix shows CR and ALIF are highly correlated (r=0.99). Entropy and SCI have poor dynamic range. BSI fails (wrong direction).")}

<div class="note">
<strong>Conclusion:</strong> CR remains the best single indicator with 10&times; dynamic range
between clean and contaminated spectra. ALIF (bounded [0, 100%]) is recommended
as a complementary metric. Spectral Entropy and SCI lack discrimination power
for this specific task.
</div>

<h3 id="sec56">5.6 HTML Report</h3>
<pre><code>generate_cr_report(
  cr_data,
  heatmap_csvs = slice_csvs,   # vector of slice CSV paths
  output_dir   = "output",
  project_name = "MyProject"
)
# Output: output/MyProject_CR_Report.html</code></pre>

<!-- ================================================================ -->
<h2 id="sec6">6 Advanced</h2>

<h3>6.1 Custom Dividing Line</h3>
<p>The default dividing line passes through (350, 0.8) and (950, 1.3) in the
(m/z, 1/K<sub>0</sub>) space. To customize:</p>
<pre><code># Python: use --line_points
python extract_precursors.py sample.d --line_points 400,0.85,900,1.25

# R: pass line_point1/line_point2
plot_cr_heatmap(csv, line_point1 = c(400, 0.85),
                     line_point2 = c(900, 1.25))</code></pre>

<div class="info">
<strong>Tip:</strong> The diaPASEF paper (Meier et al., 2020) uses
<code>1/K<sub>0</sub> = 0.0009 &times; m/z + 0.48</code> as the dividing line,
which is very similar to the default parameters.
</div>

<h3>6.2 Adjust Slice Count</h3>
<pre><code># Use 20 equally-spaced slices instead of 10
python extract_precursors.py sample.d --n_slices 20</code></pre>

<h3>6.3 Complete User Script Template</h3>
<pre><code>library(msProteomiX)
setup_workdir()

# Configuration
PROJECT  &lt;- "CR_Analysis"
CR_DIR   &lt;- "wkdir/cr_output"
OUTPUT   &lt;- "output"
GROUP_MAP &lt;- list(
  "HeLa"          = c("Hela"),
  "W/o Ext. Wash" = c("DDM"),
  "W/ Ext. Wash"  = c("FFHE")
)

# Import
cr_data &lt;- import_cr_results(CR_DIR, group_map = GROUP_MAP)

# Plot
plot_cr_gradient(cr_data, output_dir = OUTPUT, project_name = PROJECT)
plot_cr_tic(cr_data, output_dir = OUTPUT, project_name = PROJECT)

# Heatmap (if slice CSVs exist)
slice_csvs &lt;- list.files(CR_DIR, pattern = "slice_rt",
                          recursive = TRUE, full.names = TRUE)
for (csv in slice_csvs) {{
  plot_cr_heatmap(csv, output_dir = OUTPUT, project_name = PROJECT)
}}

# Report
generate_cr_report(cr_data, heatmap_csvs = slice_csvs,
                   output_dir = OUTPUT, project_name = PROJECT)</code></pre>

<!-- ================================================================ -->
<h2 id="secA">Appendix A: Parameter Reference</h2>

<h3>Python Script Parameters</h3>
<table>
<tr><th>Parameter</th><th>Type</th><th>Default</th><th>Description</th></tr>
<tr><td><code>input</code></td><td>path</td><td>required</td><td>Input .d folder, .hdf file, or directory</td></tr>
<tr><td><code>--output_dir</code></td><td>path</td><td>./cr_output</td><td>Output directory</td></tr>
<tr><td><code>--n_slices</code></td><td>int</td><td>10</td><td>Number of RT slices (equally spaced)</td></tr>
<tr><td><code>--line_points</code></td><td>str</td><td>350,0.8,950,1.3</td><td>mz1,im1,mz2,im2 for dividing line</td></tr>
<tr><td><code>--export_slices</code></td><td>flag</td><td>False</td><td>Export per-slice CSV files</td></tr>
<tr><td><code>--heatmap_rts</code></td><td>str</td><td>None</td><td>Comma-separated RT values (min)</td></tr>
</table>

<h3>R Function Parameters</h3>
<table>
<tr><th>Function</th><th>Key Parameters</th></tr>
<tr><td><code>import_cr_results()</code></td><td><code>cr_dir</code>, <code>pattern</code>, <code>group_map</code></td></tr>
<tr><td><code>plot_cr_gradient()</code></td><td><code>cr_data</code>, <code>sd_factor</code>, <code>xlim</code>, <code>ylim</code></td></tr>
<tr><td><code>plot_cr_tic()</code></td><td><code>cr_data</code>, <code>group_col</code></td></tr>
<tr><td><code>plot_cr_heatmap()</code></td><td><code>slice_csv</code>, <code>line_point1</code>, <code>line_point2</code>, <code>log_intensity</code></td></tr>
<tr><td><code>generate_cr_report()</code></td><td><code>cr_data</code>, <code>heatmap_csvs</code>, <code>sd_factor</code></td></tr>
</table>

<!-- ================================================================ -->
<h2 id="secB">Appendix B: Troubleshooting</h2>

<table>
<tr><th>Problem</th><th>Cause</th><th>Solution</th></tr>
<tr><td>CR values &gt;100%</td>
    <td>macOS m/z estimation error</td>
    <td>Convert .d to .hdf on Windows first</td></tr>
<tr><td><code>ModuleNotFoundError: alphatims</code></td>
    <td>alphatims not installed</td>
    <td><code>pip install alphatims</code></td></tr>
<tr><td><code>No precursors in frame 0</code></td>
    <td>First frame is empty (normal)</td>
    <td>This is expected; slice 1 will be skipped</td></tr>
<tr><td>Heatmap looks like a rectangle</td>
    <td>macOS estimated m/z (no IM-m/z correlation)</td>
    <td>Use .hdf converted on Windows</td></tr>
<tr><td><code>Bruker DLL not available</code></td>
    <td>macOS or Linux without Bruker SDK</td>
    <td>Expected on macOS; use .hdf files for accuracy</td></tr>
<tr><td>No CR summary files found</td>
    <td>Wrong directory path</td>
    <td>Check <code>cr_dir</code> points to the Python output</td></tr>
</table>

<!-- ================================================================ -->
<h2 id="secC">Appendix C: CR Formula Derivation</h2>

<h3>Physical Basis</h3>
<p>In trapped ion mobility spectrometry (TIMS), the inverse reduced mobility
1/K<sub>0</sub> relates to the collision cross section (CCS) and ion charge:</p>

<p style="text-align:center;font-size:16px;">
1/K<sub>0</sub> &prop; CCS / z
</p>

<p>For the same m/z value:</p>
<ul>
<li><strong>Singly-charged ion</strong> (z=1): mass &asymp; m/z, CCS &prop; mass<sup>2/3</sup>,
    so 1/K<sub>0</sub> &prop; (m/z)<sup>2/3</sup></li>
<li><strong>Doubly-charged ion</strong> (z=2): mass &asymp; 2&times;m/z, CCS &prop; (2&times;m/z)<sup>2/3</sup>,
    so 1/K<sub>0</sub> &prop; (2&times;m/z)<sup>2/3</sup>/2 &asymp; 0.79 &times; (m/z)<sup>2/3</sup></li>
</ul>

<p>Therefore, <strong>singly-charged ions have ~26% higher 1/K<sub>0</sub></strong> than doubly-charged
ions at the same m/z. This separation enables the dividing line approach.</p>

<h3>Dividing Line</h3>
<p>The line passes through two empirically-determined points in the (m/z, 1/K<sub>0</sub>) plane:</p>
<p style="text-align:center;">
(350, 0.8) and (950, 1.3)
</p>
<p>Slope: k = (1.3 &minus; 0.8) / (950 &minus; 350) = <strong>0.000833</strong><br>
Intercept: b = 0.8 &minus; 0.000833 &times; 350 = <strong>0.5083</strong></p>
<p style="text-align:center;font-size:16px;">
<strong>1/K<sub>0</sub> = 0.000833 &times; m/z + 0.5083</strong>
</p>

<h3>CR Computation</h3>
<p>For each RT slice, all MS1 detector events are classified:</p>
<ul>
<li>1/K<sub>0</sub> &ge; line value &rarr; <strong>Contaminant</strong> (singly-charged)</li>
<li>1/K<sub>0</sub> &lt; line value &rarr; <strong>Signal</strong> (multiply-charged)</li>
</ul>
<p style="text-align:center;font-size:16px;">
<strong>CR = &Sigma;I<sub>contam</sub> / &Sigma;I<sub>signal</sub> &times; 100%</strong>
</p>

<footer>
msProteomiX v0.3.0 &mdash; BayOmics<br>
Generated by generate_cr_manual.py
</footer>

</body>
</html>"""

    with open(OUTPUT_HTML, "w", encoding="utf-8") as f:
        f.write(html)
    print(f"Generated: {OUTPUT_HTML}")
    print(f"File size: {os.path.getsize(OUTPUT_HTML) / 1024:.1f} KB")


if __name__ == "__main__":
    build_html()
