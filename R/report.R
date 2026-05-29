# ==============================================================================
# msProteomiX — Automated Report Generation
# ==============================================================================

#' Generate analysis report
#'
#' Scans the output directory for analysis results and generates
#' an HTML or PDF report with all figures and summary tables.
#'
#' @param output_dir Output directory containing analysis results
#' @param project_name Project name
#' @param format Report format: "html" (default) or "pdf"
#' @param open_report Whether to open the report after generation (default TRUE)
#' @return Path to generated report
#' @export
generate_report <- function(output_dir = "output",
                             project_name = "Project",
                             format = "html",
                             open_report = TRUE) {
  if (!requireNamespace("rmarkdown", quietly = TRUE)) {
    stop("Please install rmarkdown: install.packages('rmarkdown')")
  }
  if (!requireNamespace("knitr", quietly = TRUE)) {
    stop("Please install knitr: install.packages('knitr')")
  }

  if (!dir.exists(output_dir)) {
    stop("Output directory not found: ", output_dir)
  }

  # Find all PDF files in output
  pdf_files <- list.files(output_dir, pattern = "\\.pdf$", full.names = TRUE)
  csv_files <- list.files(output_dir, pattern = "\\.csv$", full.names = TRUE)

  if (length(pdf_files) == 0) {
    stop("No analysis results found in: ", output_dir,
         "\nPlease run analysis scripts first.")
  }

  # Locate template
  template_path <- system.file("report_template.Rmd", package = "msProteomiX")
  if (template_path == "") {
    # Fallback: generate in-memory
    template_path <- .generate_rmd_template(output_dir, project_name,
                                             pdf_files, csv_files)
  }

  # Output format
  out_format <- if (format == "pdf") "pdf_document" else "html_document"
  report_ext <- if (format == "pdf") ".pdf" else ".html"
  report_file <- file.path(output_dir,
                            paste0(project_name, "_Report", report_ext))

  message(sprintf(">>> Generating %s report...", toupper(format)))

  tryCatch({
    rmarkdown::render(
      input = template_path,
      output_format = out_format,
      output_file = basename(report_file),
      output_dir = output_dir,
      params = list(
        project_name = project_name,
        output_dir = normalizePath(output_dir)
      ),
      quiet = TRUE
    )

    message(sprintf(">>> Report generated: %s", report_file))

    if (open_report && interactive()) {
      utils::browseURL(report_file)
    }

    invisible(report_file)

  }, error = function(e) {
    message(">>> Report generation failed: ", e$message)
    message(">>> Falling back to simple report...")
    .generate_simple_report(output_dir, project_name, pdf_files, csv_files)
  })
}


#' Generate RMarkdown template dynamically
#' @keywords internal
.generate_rmd_template <- function(output_dir, project_name,
                                    pdf_files, csv_files) {

  # Classify files by analysis type
  sections <- .classify_output_files(pdf_files)

  # Build Rmd content
  rmd_lines <- c(
    "---",
    sprintf("title: \"%s - Proteomics Analysis Report\"", project_name),
    sprintf("date: \"`r Sys.Date()`\""),
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_float: true",
    "    theme: flatly",
    "    code_folding: hide",
    "params:",
    sprintf("  project_name: \"%s\"", project_name),
    sprintf("  output_dir: \"%s\"", normalizePath(output_dir)),
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE,",
    "                      fig.align = 'center', out.width = '90%')",
    "output_dir <- params$output_dir",
    "```",
    "",
    sprintf("# %s", project_name),
    "",
    sprintf("Analysis report generated on `r Sys.Date()`."),
    ""
  )

  # Add sections
  section_order <- c(
    "Identification" = "Identification & Quantification",
    "Venn"           = "Venn Diagram",
    "PCA"            = "PCA Analysis",
    "CV"             = "Coefficient of Variation",
    "Corr"           = "Correlation Analysis",
    "QC"             = "Quality Control",
    "Volcano"        = "Volcano Plot",
    "Diff"           = "Differential Expression",
    "Heatmap"        = "Heatmap",
    "GO"             = "GO Enrichment",
    "KEGG"           = "KEGG Enrichment",
    "GSEA"           = "GSEA Analysis",
    "Marker"         = "Marker Gene Expression",
    "PPI"            = "Protein-Protein Interaction",
    "Reactome"       = "Reactome Enrichment"
  )

  for (key in names(section_order)) {
    if (key %in% names(sections)) {
      files <- sections[[key]]
      rmd_lines <- c(rmd_lines,
        sprintf("## %s", section_order[key]),
        ""
      )
      for (f in files) {
        fname <- basename(f)
        # Convert PDF to PNG for HTML embedding
        rmd_lines <- c(rmd_lines,
          sprintf("### %s", sub("^.*?_", "", sub("\\.pdf$", "", fname))),
          "",
          sprintf("```{r, out.width='100%%'}"),
          sprintf("f <- file.path(output_dir, \"%s\")", fname),
          "if (file.exists(f)) {",
          "  # Try to include PDF as image",
          sprintf("  knitr::include_graphics(f)"),
          "}",
          "```",
          ""
        )
      }
    }
  }

  # Summary tables
  rmd_lines <- c(rmd_lines,
    "## Summary Tables",
    "",
    "```{r}",
    "csv_files <- list.files(output_dir, pattern = '\\\\.csv$', full.names = TRUE)",
    "csv_files <- csv_files[!grepl('group_info', csv_files)]",
    "if (length(csv_files) > 0) {",
    "  cat('\\n')",
    "  for (f in csv_files) {",
    "    cat(sprintf('- %s\\n', basename(f)))",
    "  }",
    "}",
    "```",
    "",
    "---",
    "",
    "*Report generated by msProteomiX*"
  )

  # Write temp Rmd file
  rmd_path <- file.path(output_dir, paste0(project_name, "_report.Rmd"))
  writeLines(rmd_lines, rmd_path)
  rmd_path
}


#' Classify output files by analysis type
#' @keywords internal
.classify_output_files <- function(pdf_files) {
  basenames <- basename(pdf_files)

  patterns <- c(
    "Identification" = "Identification",
    "Venn"           = "Venn",
    "PCA"            = "PCA",
    "CV"             = "CV_",
    "Corr"           = "Corr",
    "QC"             = "QC_",
    "Volcano"        = "Volcano",
    "Diff"           = "Diff_",
    "Heatmap"        = "Heatmap",
    "GO"             = "GO_",
    "KEGG"           = "KEGG_",
    "GSEA"           = "GSEA",
    "Marker"         = "Marker",
    "PPI"            = "PPI_",
    "Reactome"       = "Reactome"
  )

  sections <- list()
  for (key in names(patterns)) {
    matched <- pdf_files[grepl(patterns[key], basenames, ignore.case = TRUE)]
    if (length(matched) > 0) sections[[key]] <- matched
  }
  sections
}


#' Generate a simple HTML report (fallback)
#' @keywords internal
.generate_simple_report <- function(output_dir, project_name,
                                     pdf_files, csv_files) {
  report_file <- file.path(output_dir,
                            paste0(project_name, "_Report.html"))

  html_lines <- c(
    "<!DOCTYPE html>",
    "<html><head>",
    sprintf("<title>%s - Analysis Report</title>", project_name),
    "<style>",
    "body { font-family: 'Segoe UI', sans-serif; max-width: 900px; margin: 0 auto; padding: 20px; }",
    "h1 { color: #2c3e50; border-bottom: 2px solid #3498db; }",
    "h2 { color: #34495e; }",
    ".file-list { background: #f8f9fa; padding: 15px; border-radius: 8px; }",
    ".file-list a { display: block; padding: 4px 0; color: #2980b9; }",
    "</style>",
    "</head><body>",
    sprintf("<h1>%s - Proteomics Analysis Report</h1>", project_name),
    sprintf("<p>Generated: %s</p>", Sys.Date()),
    "<h2>Analysis Results</h2>",
    "<div class='file-list'>"
  )

  for (f in sort(basename(pdf_files))) {
    html_lines <- c(html_lines,
      sprintf("<a href='%s'>%s</a>", f, sub("\\.pdf$", "", f)))
  }

  html_lines <- c(html_lines,
    "</div>",
    "<h2>Data Tables</h2>",
    "<div class='file-list'>"
  )

  for (f in sort(basename(csv_files))) {
    html_lines <- c(html_lines,
      sprintf("<a href='%s'>%s</a>", f, f))
  }

  html_lines <- c(html_lines,
    "</div>",
    "<hr><p><em>Generated by msProteomiX</em></p>",
    "</body></html>"
  )

  writeLines(html_lines, report_file)
  message(sprintf(">>> Simple report generated: %s", report_file))

  if (interactive()) utils::browseURL(report_file)
  invisible(report_file)
}
