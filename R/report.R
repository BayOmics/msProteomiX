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
  if (!dir.exists(output_dir)) {
    stop("Output directory not found: ", output_dir)
  }

  # Find all PNG files in output (saved alongside PDFs by save_plot_and_data)
  png_files <- list.files(output_dir, pattern = "\\.png$", full.names = TRUE)
  csv_files <- list.files(output_dir, pattern = "\\.csv$", full.names = TRUE)

  if (length(png_files) == 0) {
    # Fallback: check for PDF files (older runs without PNG)
    pdf_files <- list.files(output_dir, pattern = "\\.pdf$", full.names = TRUE)
    if (length(pdf_files) == 0) {
      stop("No analysis results found in: ", output_dir,
           "\nPlease run analysis scripts first.")
    }
    message(">>> No PNG files found. Generating link-only report...")
    return(.generate_simple_report(output_dir, project_name, pdf_files, csv_files))
  }

  # Use rmarkdown if available, otherwise simple HTML
  if (requireNamespace("rmarkdown", quietly = TRUE) &&
      requireNamespace("knitr", quietly = TRUE)) {

    template_path <- .generate_rmd_template(output_dir, project_name,
                                             png_files, csv_files)

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
        output_dir = normalizePath(output_dir),
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

      return(invisible(report_file))

    }, error = function(e) {
      message(">>> RMarkdown failed: ", e$message)
      message(">>> Falling back to simple report...")
    })
  }

  # Fallback: simple HTML with embedded images
  .generate_simple_report_with_images(output_dir, project_name,
                                       png_files, csv_files, open_report)
}


#' Generate RMarkdown template dynamically
#' @keywords internal
.generate_rmd_template <- function(output_dir, project_name,
                                    png_files, csv_files) {

  # Classify files by analysis type
  sections <- .classify_output_files(png_files)
  abs_output <- normalizePath(output_dir)

  # Build Rmd content
  rmd_lines <- c(
    "---",
    sprintf("title: \"%s - Proteomics Analysis Report\"", project_name),
    "date: \"`r Sys.Date()`\"",
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_float: true",
    "    theme: flatly",
    "    self_contained: true",
    "params:",
    sprintf("  project_name: \"%s\"", project_name),
    sprintf("  output_dir: \"%s\"", abs_output),
    "---",
    "",
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE,",
    "                      fig.align = 'center', out.width = '100%')",
    "output_dir <- params$output_dir",
    "```",
    "",
    sprintf("Analysis report for **%s**, generated on `r Sys.Date()`.", project_name),
    "",
    "---",
    ""
  )

  # Add sections
  section_order <- c(
    "Identification" = "Identification & Quantification",
    "Venn"           = "Venn Diagram",
    "UpSet"          = "UpSet Diagram",
    "PCA"            = "PCA Analysis",
    "CV"             = "Coefficient of Variation",
    "Corr"           = "Correlation Analysis",
    "QC"             = "Quality Control",
    "Volcano"        = "Volcano Plot",
    "Diff"           = "Differential Expression",
    "Heatmap"        = "Heatmap",
    "TopHeatmap"     = "Top Protein Heatmap",
    "Coverage"       = "Sequence Coverage",
    "GO"             = "GO Enrichment",
    "KEGG"           = "KEGG Enrichment",
    "Reactome"       = "Reactome Enrichment",
    "GSEA"           = "GSEA Analysis",
    "Marker"         = "Marker Gene Expression",
    "PPI"            = "Protein Network"
  )

  for (key in names(section_order)) {
    if (key %in% names(sections)) {
      files <- sections[[key]]
      rmd_lines <- c(rmd_lines,
        sprintf("## %s {.tabset}", section_order[key]),
        ""
      )
      for (f in files) {
        fname <- basename(f)
        # Clean display name: remove project prefix and extension
        display_name <- sub("^.*?_", "", sub("\\.png$", "", fname))
        display_name <- gsub("_", " ", display_name)
        abs_path <- normalizePath(f)

        rmd_lines <- c(rmd_lines,
          sprintf("### %s", display_name),
          "",
          sprintf("```{r}"),
          sprintf("knitr::include_graphics(\"%s\")", abs_path),
          "```",
          ""
        )
      }
    }
  }

  # Summary tables section
  rmd_lines <- c(rmd_lines,
    "---",
    "",
    "## Data Tables",
    "",
    "The following data files are available in the output directory:",
    "",
    "```{r results='asis'}",
    sprintf("csv_files <- list.files('%s', pattern = '[.]csv$', full.names = FALSE)", abs_output),
    "csv_files <- csv_files[!grepl('group_info', csv_files)]",
    "if (length(csv_files) > 0) {",
    "  for (f in sort(csv_files)) {",
    "    cat(sprintf('- `%s`\\n', f))",
    "  }",
    "}",
    "```",
    "",
    "---",
    "",
    "*Report generated by [msProteomiX](https://github.com/BayOmics/msProteomiX)*"
  )

  # Write Rmd file
  rmd_path <- file.path(output_dir, paste0(project_name, "_report.Rmd"))
  writeLines(rmd_lines, rmd_path)
  rmd_path
}


#' Classify output files by analysis type
#' @keywords internal
.classify_output_files <- function(image_files) {
  basenames <- basename(image_files)

  patterns <- c(
    "Identification" = "Identification",
    "Venn"           = "Venn",
    "UpSet"          = "UpSet",
    "PCA"            = "PCA",
    "CV"             = "CV_",
    "Corr"           = "Corr",
    "QC"             = "QC_",
    "Volcano"        = "Volcano",
    "Diff"           = "Diff_",
    "Heatmap"        = "(?<!Top)Heatmap",
    "TopHeatmap"     = "TopHeatmap",
    "Coverage"       = "Coverage",
    "GO"             = "GO_",
    "KEGG"           = "KEGG_",
    "Reactome"       = "Reactome",
    "GSEA"           = "GSEA",
    "Marker"         = "Marker",
    "PPI"            = "PPI_"
  )

  sections <- list()
  for (key in names(patterns)) {
    matched <- image_files[grepl(patterns[key], basenames, ignore.case = TRUE, perl = TRUE)]
    if (length(matched) > 0) sections[[key]] <- matched
  }
  sections
}


#' Generate simple HTML report with embedded images (no rmarkdown needed)
#' @keywords internal
.generate_simple_report_with_images <- function(output_dir, project_name,
                                                 png_files, csv_files,
                                                 open_report = TRUE) {
  report_file <- file.path(output_dir,
                            paste0(project_name, "_Report.html"))

  sections <- .classify_output_files(png_files)

  section_order <- c(
    "Identification" = "Identification & Quantification",
    "Venn"           = "Venn Diagram",
    "UpSet"          = "UpSet Diagram",
    "PCA"            = "PCA Analysis",
    "CV"             = "Coefficient of Variation",
    "Corr"           = "Correlation Analysis",
    "QC"             = "Quality Control",
    "Volcano"        = "Volcano Plot",
    "Diff"           = "Differential Expression",
    "Heatmap"        = "Heatmap",
    "TopHeatmap"     = "Top Protein Heatmap",
    "Coverage"       = "Sequence Coverage",
    "GO"             = "GO Enrichment",
    "KEGG"           = "KEGG Enrichment",
    "Reactome"       = "Reactome Enrichment",
    "GSEA"           = "GSEA Analysis",
    "Marker"         = "Marker Gene Expression",
    "PPI"            = "Protein Network"
  )

  html <- c(
    "<!DOCTYPE html>",
    "<html><head>",
    sprintf("<title>%s - Analysis Report</title>", project_name),
    "<meta charset='utf-8'>",
    "<style>",
    "body { font-family: -apple-system, 'Segoe UI', Roboto, sans-serif; max-width: 1000px; margin: 0 auto; padding: 30px; background: #f5f5f5; }",
    ".container { background: white; padding: 40px; border-radius: 12px; box-shadow: 0 2px 10px rgba(0,0,0,0.1); }",
    "h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }",
    "h2 { color: #34495e; margin-top: 40px; border-left: 4px solid #3498db; padding-left: 12px; }",
    ".figure { text-align: center; margin: 20px 0; }",
    ".figure img { max-width: 100%; border: 1px solid #eee; border-radius: 8px; }",
    ".figure .caption { color: #7f8c8d; font-size: 0.9em; margin-top: 8px; }",
    ".file-list { background: #f8f9fa; padding: 15px 20px; border-radius: 8px; }",
    ".file-list li { padding: 3px 0; }",
    ".footer { text-align: center; color: #95a5a6; margin-top: 30px; font-size: 0.85em; }",
    "</style>",
    "</head><body>",
    "<div class='container'>",
    sprintf("<h1>%s</h1>", project_name),
    sprintf("<p>Proteomics Analysis Report &mdash; %s</p>", Sys.Date())
  )

  for (key in names(section_order)) {
    if (key %in% names(sections)) {
      html <- c(html, sprintf("<h2>%s</h2>", section_order[key]))
      for (f in sections[[key]]) {
        fname <- basename(f)
        display <- sub("^.*?_", "", sub("\\.png$", "", fname))
        display <- gsub("_", " ", display)
        html <- c(html,
          "<div class='figure'>",
          sprintf("<img src='%s' alt='%s'>", fname, display),
          sprintf("<div class='caption'>%s</div>", display),
          "</div>"
        )
      }
    }
  }

  # CSV list
  if (length(csv_files) > 0) {
    html <- c(html, "<h2>Data Tables</h2>", "<div class='file-list'><ul>")
    for (f in sort(basename(csv_files))) {
      html <- c(html, sprintf("<li><a href='%s'>%s</a></li>", f, f))
    }
    html <- c(html, "</ul></div>")
  }

  html <- c(html,
    sprintf("<div class='footer'>Generated by msProteomiX v%s</div>",
            utils::packageVersion("msProteomiX")),
    "</div></body></html>"
  )

  writeLines(html, report_file)
  message(sprintf(">>> Report generated: %s", report_file))

  if (open_report && interactive()) utils::browseURL(report_file)
  invisible(report_file)
}


#' Generate a simple link-only HTML report (fallback for no PNG)
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
    sprintf("<h1>%s - Analysis Report</h1>", project_name),
    sprintf("<p>Generated: %s</p>", Sys.Date()),
    "<h2>Analysis Results (PDF)</h2>",
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
  invisible(report_file)
}
