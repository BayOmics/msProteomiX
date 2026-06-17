# ==============================================================================
# msProteomiX — Plasma QC Evaluation Report Generator
# ==============================================================================
# Pipeline: R ggplot -> PNG (tempdir) -> base64 -> HTML -> Chrome PDF
# Pitfall #17: Use paste0 + gsub, NOT sprintf for Chinese Unicode
# Pitfall #18: Use tempdir() for PNG, NOT paths with Chinese characters
# ==============================================================================


#' Find Chrome executable path
#' @return Path to Chrome or NULL
#' @keywords internal
.find_chrome <- function() {
  candidates <- c(
    # macOS
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    # Windows (standard install locations) — only add if env var is set
    if (nzchar(Sys.getenv("ProgramFiles")))
      file.path(Sys.getenv("ProgramFiles"), "Google/Chrome/Application/chrome.exe"),
    if (nzchar(Sys.getenv("ProgramFiles(x86)")))
      file.path(Sys.getenv("ProgramFiles(x86)"), "Google/Chrome/Application/chrome.exe"),
    if (nzchar(Sys.getenv("LOCALAPPDATA")))
      file.path(Sys.getenv("LOCALAPPDATA"), "Google/Chrome/Application/chrome.exe"),
    # Linux
    Sys.which("google-chrome"),
    Sys.which("chromium-browser"),
    Sys.which("chromium"),
    Sys.which("chrome")
  )
  for (p in candidates) {
    if (nzchar(p) && file.exists(p)) return(p)
  }
  NULL
}


#' Encode PNG file as base64 img tag
#' @param path Path to PNG file
#' @param width CSS width (default "85%")
#' @param caption Figure caption
#' @return HTML string
#' @keywords internal
.img_to_b64 <- function(path, width = "85%", caption = "") {
  if (!file.exists(path)) return("")
  b64 <- base64enc::base64encode(path)
  paste0(
    '<div class="figure">',
    sprintf('<img src="data:image/png;base64,%s" style="max-width:%s">', b64, width),
    if (nzchar(caption)) sprintf('<div class="caption">%s</div>', caption) else "",
    '</div>'
  )
}


#' Generate plasma proteomics quality evaluation report
#'
#' Creates a self-contained HTML report (optionally converted to PDF via
#' Chrome headless) evaluating blood cell contamination in plasma samples.
#'
#' @param ms MsDataSet object (human plasma data)
#' @param ci_result data.frame from \code{calc_contamination_index()}
#' @param ms_fbs Optional MsDataSet for FBS control data (BOVIN FASTA search).
#'   If provided, FBS summary statistics are included in the report.
#' @param project Character: project name displayed in report header.
#' @param instrument Character: instrument name (default "Thermo Orbitrap").
#' @param search_mode Character: search strategy, e.g. "directDIA", "DDA".
#' @param lang Character: "en" for English, "zh" for Chinese (default "en").
#' @param output_dir Character: output directory for the report.
#' @param chrome_path Character or NULL: path to Chrome/Chromium executable.
#'   NULL triggers auto-detection. If Chrome is unavailable, only HTML is generated.
#'
#' @return Path to the generated report (PDF if Chrome available, otherwise HTML).
#'
#' @details
#' The report includes:
#' \enumerate{
#'   \item Data overview (samples, protein groups, identification counts)
#'   \item CI methodology explanation (PLT / RBC / PBMC panels)
#'   \item CI results with Pass/Flag/Fail judgments
#'   \item Protein coverage analysis
#'   \item FBS control summary (if provided)
#'   \item QC metrics (peptide length, charge state, etc.)
#'   \item Conclusions and recommendations
#' }
#'
#' @export
generate_plasma_report <- function(ms, ci_result,
                                    ms_fbs = NULL,
                                    project = "Plasma QC Evaluation",
                                    instrument = "Thermo Orbitrap",
                                    search_mode = "directDIA",
                                    lang = "en",
                                    output_dir = "output",
                                    chrome_path = NULL) {
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    stop("Package 'base64enc' is required. Install with: install.packages('base64enc')")
  }

  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  # Use tempdir for PNG to avoid non-ASCII path issues (Pitfall #18)
  png_dir <- file.path(tempdir(), "msProteomiX_plasma_report_png")
  dir.create(png_dir, recursive = TRUE, showWarnings = FALSE)

  message(">>> [Report] Generating plots ...")

  # ============================================================
  # Module B: Generate Plots
  # ============================================================

  n_prot <- nrow(ms$proteins)
  n_sample <- ncol(ms$proteins)
  samples <- colnames(ms$proteins)
  prot_info <- ms$protein_info

  # B1: CI Barplot
  jcol <- c(Pass = "#27ae60", Flag = "#f39c12", Fail = "#e74c3c", "N/A" = "#bdc3c7")

  ci_long <- rbind(
    data.frame(Sample = gsub("\\.raw$", "", ci_result$Sample), Panel = "CI_Platelet (Gao2026)",
               Value = ci_result$CI_PLT, Judgment = ci_result$J_PLT, stringsAsFactors = FALSE),
    data.frame(Sample = gsub("\\.raw$", "", ci_result$Sample), Panel = "CI_Erythrocyte (Gao2026)",
               Value = ci_result$CI_RBC, Judgment = ci_result$J_RBC, stringsAsFactors = FALSE),
    data.frame(Sample = gsub("\\.raw$", "", ci_result$Sample), Panel = "CI_PBMC (Korff2025)",
               Value = ci_result$CI_PBMC, Judgment = ci_result$J_PBMC, stringsAsFactors = FALSE)
  )
  ci_long$Panel <- factor(ci_long$Panel,
    levels = c("CI_Platelet (Gao2026)", "CI_Erythrocyte (Gao2026)", "CI_PBMC (Korff2025)"))

  thresh_df <- data.frame(
    Panel = factor(c(rep("CI_Platelet (Gao2026)", 2), "CI_Erythrocyte (Gao2026)",
                     rep("CI_PBMC (Korff2025)", 2)), levels = levels(ci_long$Panel)),
    yval = c(0.009, 0.05, 0.1, 0.005, 0.02),
    Label = c("Pass 0.009", "Fail 0.05", "Fail 0.1", "Pass 0.005", "Fail 0.02")
  )

  p_ci <- ggplot2::ggplot(ci_long, ggplot2::aes(x = .data$Sample, y = .data$Value, fill = .data$Judgment)) +
    ggplot2::geom_col(width = 0.6) +
    ggplot2::facet_wrap(~ .data$Panel, scales = "free_y", ncol = 1) +
    ggplot2::geom_hline(data = thresh_df, ggplot2::aes(yintercept = .data$yval),
                        linetype = "dashed", color = "grey40", linewidth = 0.5) +
    ggplot2::geom_text(data = thresh_df,
                       ggplot2::aes(x = Inf, y = .data$yval, label = .data$Label),
                       hjust = 1.1, vjust = -0.5, size = 2.8, color = "grey40",
                       inherit.aes = FALSE) +
    ggplot2::scale_fill_manual(values = jcol, name = "Judgment") +
    ggplot2::labs(title = "Contamination Index (CI) Evaluation",
                  x = NULL, y = "Contamination Index") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                   strip.text = ggplot2::element_text(face = "bold"))
  ggplot2::ggsave(file.path(png_dir, "fig1_ci_barplot.png"), p_ci, width = 8, height = 8, dpi = 300)
  message("    fig1 (CI barplot) OK")

  # B2: Protein Rank
  mat_raw <- if (isTRUE(ms$is_log2)) 2^as.matrix(ms$proteins) else as.matrix(ms$proteins)
  mat_raw[is.na(mat_raw)] <- 0
  med_intensity <- apply(mat_raw, 1, stats::median, na.rm = TRUE)
  rank_df <- data.frame(Rank = seq_along(sort(med_intensity, decreasing = TRUE)),
                         log2I = log2(sort(med_intensity, decreasing = TRUE) + 1))
  p_rank <- ggplot2::ggplot(rank_df, ggplot2::aes(x = .data$Rank, y = .data$log2I)) +
    ggplot2::geom_point(size = 0.5, alpha = 0.5, color = "#2c3e50") +
    ggplot2::labs(title = "Protein Rank by Median Intensity",
                  x = "Protein Rank", y = "log2(Median Intensity + 1)") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"))
  ggplot2::ggsave(file.path(png_dir, "fig2_protein_rank.png"), p_rank, width = 8, height = 5, dpi = 300)
  message("    fig2 (protein rank) OK")

  # B3: Intensity distribution boxplot
  mat_log2 <- if (isTRUE(ms$is_log2)) as.matrix(ms$proteins) else log2(as.matrix(ms$proteins) + 1)
  int_long <- data.frame(
    Sample = rep(colnames(mat_log2), each = nrow(mat_log2)),
    log2I = as.vector(mat_log2),
    stringsAsFactors = FALSE
  )
  int_long$Sample <- gsub("\\.raw$", "", int_long$Sample)
  int_long <- int_long[!is.na(int_long$log2I) & is.finite(int_long$log2I), ]

  p_box <- ggplot2::ggplot(int_long, ggplot2::aes(x = .data$Sample, y = .data$log2I, fill = .data$Sample)) +
    ggplot2::geom_boxplot(outlier.size = 0.3, alpha = 0.8) +
    ggplot2::labs(title = "Intensity Distribution per Sample",
                  x = NULL, y = "log2(Intensity)") +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
                   axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
                   legend.position = "none")
  ggplot2::ggsave(file.path(png_dir, "fig3_intensity_boxplot.png"), p_box, width = 8, height = 5, dpi = 300)
  message("    fig3 (intensity boxplot) OK")

  # ============================================================
  # Module C: Build HTML
  # ============================================================
  message(">>> [Report] Building HTML ...")

  samples_clean <- gsub("\\.raw$", "", samples)

  # CI results table rows
  ci_rows <- ""
  for (i in seq_len(nrow(ci_result))) {
    s <- gsub("\\.raw$", "", ci_result$Sample[i])
    ci_rows <- paste0(ci_rows, sprintf(
      '<tr><td><b>%s</b></td><td>%s</td><td>%.4f</td><td>%.4f</td><td>%.4f</td><td>%s</td><td>%s</td><td>%s</td><td><b>%s</b></td></tr>\n',
      s, format(ci_result$PGs[i], big.mark = ","),
      ci_result$CI_PLT[i], ci_result$CI_RBC[i], ci_result$CI_PBMC[i],
      ci_result$J_PLT[i], ci_result$J_RBC[i], ci_result$J_PBMC[i],
      ci_result$Overall[i]
    ))
  }

  # Match info
  match_info <- attr(ci_result, "match_info")
  match_text <- ""
  if (!is.null(match_info)) {
    for (nm in names(match_info)) {
      mi <- match_info[[nm]]
      match_text <- paste0(match_text, sprintf("%s: %d/%d matched. ", nm, mi$matched, mi$total))
    }
  }

  # FBS section
  fbs_section <- ""
  if (!is.null(ms_fbs)) {
    fbs_n <- nrow(ms_fbs$proteins)
    fbs_section <- paste0(
      '<h2>5. FBS Control</h2>',
      '<p>FBS was searched independently with <b>BOVIN FASTA</b>, identifying <b>',
      format(fbs_n, big.mark = ","), ' bovine protein groups</b>.</p>',
      '<div class="note">CI marker panels (PLT / RBC / PBMC) are human-specific ',
      'UniProt IDs and are NOT applicable to bovine proteins.</div>'
    )
  }

  # CSS
  css <- '
<style>
  @page { margin: 2cm 2.5cm; size: A4; }
  body { font-family: "Helvetica Neue", "PingFang SC", "Noto Sans SC", sans-serif;
         font-size: 11pt; line-height: 1.6; color: #333; max-width: 800px; margin: 0 auto; padding: 20px; }
  h1 { font-size: 20pt; border-bottom: 3px solid #2c3e50; padding-bottom: 8px; color: #2c3e50; }
  h2 { font-size: 15pt; color: #2c3e50; border-bottom: 1px solid #bdc3c7; padding-bottom: 4px;
       margin-top: 30px; }
  table { border-collapse: collapse; width: 100%; margin: 10px 0; font-size: 10pt; }
  th { background: #2c3e50; color: white; padding: 8px 10px; text-align: center; }
  td { border: 1px solid #ddd; padding: 6px 10px; }
  tr:nth-child(even) { background: #f9f9f9; }
  .figure { margin: 15px 0; page-break-inside: avoid; text-align: center; }
  .caption { font-size: 9pt; color: #7f8c8d; margin-top: 4px; text-align: center; }
  .note { background: #eaf2f8; border-left: 4px solid #3498db; padding: 10px 15px; margin: 10px 0; font-size: 10pt; }
  .header-meta { color: #7f8c8d; font-size: 10pt; }
  code { background: #f4f4f4; padding: 2px 5px; border-radius: 3px; font-size: 9.5pt; }
  hr { border: none; border-top: 1px solid #ddd; margin: 25px 0; }
</style>
'

  # Build HTML body
  html <- paste0(
    '<!DOCTYPE html><html><head><meta charset="utf-8">',
    '<title>', project, ' - Plasma QC Report</title>', css, '</head><body>',

    '<h1>', project, '</h1>',
    '<p class="header-meta">Plasma Proteomics Quality Evaluation Report</p>',
    '<p class="header-meta">Generated: ', Sys.Date(), ' | Engine: ',
    ms$engine, ' | Mode: ', search_mode, ' | Instrument: ', instrument, '</p>',
    '<hr>',

    # Section 1: Overview
    '<h2>1. Data Overview</h2>',
    '<table>',
    '<tr><th>Metric</th><th>Value</th></tr>',
    '<tr><td>Search Engine</td><td>', ms$engine, '</td></tr>',
    '<tr><td>Search Mode</td><td>', search_mode, '</td></tr>',
    '<tr><td>Samples</td><td>', n_sample, '</td></tr>',
    '<tr><td>Protein Groups</td><td>', format(n_prot, big.mark = ","), '</td></tr>',
    '</table>',
    '<hr>',

    # Section 2: CI Methodology
    '<h2>2. Contamination Assessment Methodology</h2>',
    '<p>Contamination Index (CI) = SUM(marker protein intensities) / SUM(all protein intensities).</p>',
    '<table>',
    '<tr><th>Panel</th><th>Source</th><th>Markers</th><th>Pass</th><th>Fail</th></tr>',
    '<tr><td>PLT (Platelet)</td><td>Gao 2026</td><td>30 UniProt</td>',
    '<td>&le;0.009</td><td>&gt;0.05</td></tr>',
    '<tr><td>RBC (Erythrocyte)</td><td>Gao 2026</td><td>30 UniProt</td>',
    '<td>&le;0.05</td><td>&gt;0.1</td></tr>',
    '<tr><td>PBMC</td><td>Korff 2025</td><td>5 Gene</td>',
    '<td>&lt;0.005</td><td>&gt;0.02</td></tr>',
    '</table>',
    '<div class="note">', match_text, '</div>',
    '<hr>',

    # Section 3: CI Results
    '<h2>3. Contamination Index Results</h2>',
    .img_to_b64(file.path(png_dir, "fig1_ci_barplot.png"), "90%", "Figure 1: CI Barplot"),
    '<table>',
    '<tr><th>Sample</th><th>PGs</th><th>CI_PLT</th><th>CI_RBC</th><th>CI_PBMC</th>',
    '<th>J_PLT</th><th>J_RBC</th><th>J_PBMC</th><th>Overall</th></tr>',
    ci_rows,
    '</table>',
    '<hr>',

    # Section 4: Intensity
    '<h2>4. Protein Quantification</h2>',
    .img_to_b64(file.path(png_dir, "fig3_intensity_boxplot.png"), "85%", "Figure 2: Intensity Distribution"),
    .img_to_b64(file.path(png_dir, "fig2_protein_rank.png"), "85%", "Figure 3: Protein Rank"),
    '<hr>',

    # Section 5: FBS (optional)
    fbs_section,
    if (nzchar(fbs_section)) '<hr>' else '',

    # Section 6: Conclusions
    '<h2>', if (is.null(ms_fbs)) '5' else '6', '. Conclusions</h2>',
    '<ol>',
    '<li>Contamination assessment completed for ', n_sample, ' samples.</li>',
    '<li>Overall results: ',
    sum(ci_result$Overall == "Pass"), ' Pass, ',
    sum(ci_result$Overall == "Flag"), ' Flag, ',
    sum(ci_result$Overall == "Fail"), ' Fail.</li>',
    '</ol>',
    '<hr>',

    # Section 7: Methodology
    '<h2>', if (is.null(ms_fbs)) '6' else '7', '. References</h2>',
    '<ol>',
    '<li>Geyer PE et al. <i>EMBO Mol Med</i> 2019; 11(11): e10427.</li>',
    '<li>Gao et al. <i>EMBO Mol Med</i> 2026 (NP-specific panels).</li>',
    '<li>Korff et al. <i>bioRxiv</i> 2025 (PBMC panel).</li>',
    '</ol>',
    '<p><em>Report generated by msProteomiX v',
    as.character(utils::packageVersion("msProteomiX")), '</em></p>',

    '</body></html>'
  )

  # ============================================================
  # Module D: Chinese translation (if lang == "zh")
  # ============================================================
  if (lang == "zh") {
    message(">>> [Report] Applying Chinese translation ...")
    # Using gsub translation layer (Pitfall #17: avoid sprintf with CJK)
    html <- gsub("Plasma Proteomics Quality Evaluation Report",
                 "\u8840\u6d46\u86cb\u767d\u7ec4\u5b66\u6570\u636e\u8d28\u91cf\u8bc4\u4f30\u62a5\u544a", html, fixed = TRUE)
    html <- gsub(">1. Data Overview<", ">1. \u6570\u636e\u6982\u89c8<", html, fixed = TRUE)
    html <- gsub(">Search Engine<", ">\u641c\u5e93\u5f15\u64ce<", html, fixed = TRUE)
    html <- gsub(">Search Mode<", ">\u641c\u7d22\u6a21\u5f0f<", html, fixed = TRUE)
    html <- gsub(">Samples<", ">\u6837\u672c\u6570<", html, fixed = TRUE)
    html <- gsub(">Protein Groups<", ">\u86cb\u767d\u7ec4\u6570<", html, fixed = TRUE)
    html <- gsub(">2. Contamination Assessment Methodology<",
                 ">2. \u6c61\u67d3\u8bc4\u4f30\u65b9\u6cd5<", html, fixed = TRUE)
    html <- gsub("Contamination Index (CI) = SUM(marker protein intensities) / SUM(all protein intensities).",
                 "\u6c61\u67d3\u6307\u6570 (CI) = SUM(\u6807\u8bb0\u86cb\u767d\u5f3a\u5ea6) / SUM(\u5168\u90e8\u86cb\u767d\u5f3a\u5ea6)\u3002", html, fixed = TRUE)
    html <- gsub("<th>Panel</th><th>Source</th><th>Markers</th><th>Pass</th><th>Fail</th>",
                 "<th>\u9762\u677f</th><th>\u6765\u6e90</th><th>\u6807\u8bb0\u86cb\u767d</th><th>Pass</th><th>Fail</th>", html, fixed = TRUE)
    html <- gsub(">3. Contamination Index Results<",
                 ">3. \u6c61\u67d3\u6307\u6570\u7ed3\u679c<", html, fixed = TRUE)
    html <- gsub(">4. Protein Quantification<",
                 ">4. \u86cb\u767d\u5b9a\u91cf\u5206\u6790<", html, fixed = TRUE)
    html <- gsub(">5. FBS Control<", ">5. FBS \u5bf9\u7167\u5206\u6790<", html, fixed = TRUE)
    html <- gsub("FBS was searched independently with",
                 "FBS \u4f7f\u7528", html, fixed = TRUE)
    html <- gsub("bovine protein groups",
                 "\u4e2a\u725b\u6e90\u86cb\u767d\u7ec4", html, fixed = TRUE)
    html <- gsub("CI marker panels (PLT / RBC / PBMC) are human-specific UniProt IDs and are NOT applicable to bovine proteins.",
                 "CI \u6807\u8bb0\u86cb\u767d (PLT / RBC / PBMC) \u5747\u4e3a\u4eba\u6e90 UniProt ID\uff0c\u4e0d\u9002\u7528\u4e8e\u725b\u6e90\u86cb\u767d\u3002", html, fixed = TRUE)
    html <- gsub(">Conclusions<", ">\u7ed3\u8bba<", html, fixed = TRUE)
    html <- gsub("Contamination assessment completed for",
                 "\u5df2\u5b8c\u6210", html, fixed = TRUE)
    html <- gsub("samples.", "\u4e2a\u6837\u672c\u7684\u6c61\u67d3\u8bc4\u4f30\u3002", html, fixed = TRUE)
    html <- gsub("Overall results:", "\u7efc\u5408\u7ed3\u679c:", html, fixed = TRUE)
    html <- gsub(">References<", ">\u53c2\u8003\u6587\u732e<", html, fixed = TRUE)
    html <- gsub("Figure 1:", "\u56fe 1:", html, fixed = TRUE)
    html <- gsub("Figure 2:", "\u56fe 2:", html, fixed = TRUE)
    html <- gsub("Figure 3:", "\u56fe 3:", html, fixed = TRUE)
    html <- gsub(">Metric</th><th>Value<", ">\u6307\u6807</th><th>\u6570\u503c<", html, fixed = TRUE)
    html <- gsub("Generated:", "\u751f\u6210\u65e5\u671f:", html, fixed = TRUE)
    html <- gsub("Report generated by", "\u62a5\u544a\u7531", html, fixed = TRUE)
  }

  # ============================================================
  # Module E: Output (HTML + PDF)
  # ============================================================
  lang_suffix <- if (lang == "zh") "_zh" else ""
  html_path <- file.path(output_dir, paste0("plasma_qc_report", lang_suffix, ".html"))
  writeLines(html, html_path, useBytes = TRUE)
  message(sprintf(">>> HTML report saved: %s", html_path))

  # Try Chrome PDF
  if (is.null(chrome_path)) chrome_path <- .find_chrome()

  if (!is.null(chrome_path)) {
    pdf_path <- file.path(output_dir, paste0("plasma_qc_report", lang_suffix, ".pdf"))
    tmp_pdf <- file.path(tempdir(), "plasma_qc_report.pdf")

    cmd <- sprintf('"%s" --headless --disable-gpu --no-sandbox --print-to-pdf="%s" "%s"',
                   chrome_path, tmp_pdf, html_path)
    res <- system(cmd, intern = FALSE, ignore.stdout = TRUE, ignore.stderr = TRUE)

    if (res == 0 && file.exists(tmp_pdf)) {
      file.copy(tmp_pdf, pdf_path, overwrite = TRUE)
      file.remove(tmp_pdf)
      message(sprintf(">>> PDF report saved: %s", pdf_path))
      # Cleanup PNG temp dir
      unlink(png_dir, recursive = TRUE)
      return(invisible(pdf_path))
    } else {
      message(">>> Chrome PDF conversion failed. HTML report available.")
    }
  } else {
    message(">>> Chrome not found. HTML-only report generated.")
    message("    To generate PDF, install Google Chrome or pass chrome_path argument.")
  }

  # Cleanup
  unlink(png_dir, recursive = TRUE)
  invisible(html_path)
}
