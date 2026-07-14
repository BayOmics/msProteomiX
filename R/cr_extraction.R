# ==============================================================================
# msProteomiX - CR data extraction via reticulate + alphatims
# ==============================================================================
# Provides extract_cr_data() which wraps the bundled Python script
# extract_precursors.py, calling it through reticulate so users
# can stay entirely in RStudio.
# ==============================================================================

#' Extract CR data from timsTOF .d/.hdf files
#'
#' Calls the bundled Python script \code{extract_precursors.py} via
#' \code{reticulate} to extract MS1 precursors from timsTOF data and
#' compute the Contamination Ratio (CR) per RT slice.
#'
#' @param input_path Path to a \code{.d} folder, \code{.hdf} file,
#'   or a directory containing multiple such files.
#' @param output_dir Output directory for results
#'   (default: \code{"wkdir/cr_output"}).
#' @param n_slices Number of equally-spaced RT slices (default: 10).
#' @param export_slices If \code{TRUE}, export per-slice CSV files
#'   for heatmap visualization (default: \code{TRUE}).
#' @param heatmap_rts Optional numeric vector of RT values (minutes)
#'   for heatmap export. \code{NULL} means use all slices.
#' @param line_points CR dividing line definition as
#'   \code{"mz1,im1,mz2,im2"} (default: \code{"350,0.8,950,1.3"}).
#' @param python Path to a Python executable. If \code{NULL},
#'   uses the \code{reticulate} default.
#'
#' @details
#' This function requires:
#' \itemize{
#'   \item A working Python installation (3.8+)
#'   \item The \code{alphatims} Python package
#'     (\code{pip install alphatims})
#'   \item The R package \code{reticulate}
#' }
#'
#' On macOS, \code{.d} files cannot be read accurately because the
#' Bruker native DLL is only available on Windows. Convert \code{.d}
#' to \code{.hdf} on Windows first (see the CR Manual for details),
#' then transfer the \code{.hdf} files to macOS for analysis.
#'
#' @return Invisibly returns the path to the output directory.
#'
#' @examples
#' \dontrun{
#' # Single file
#' extract_cr_data("D:/data/sample.d")
#'
#' # Batch processing
#' extract_cr_data("D:/data/")
#'
#' # Custom RT slices for heatmap
#' extract_cr_data("sample.hdf", heatmap_rts = c(20, 40, 60))
#' }
#'
#' @export
extract_cr_data <- function(input_path,
                            output_dir = "wkdir/cr_output",
                            n_slices = 10,
                            export_slices = TRUE,
                            heatmap_rts = NULL,
                            line_points = "350,0.8,950,1.3",
                            python = NULL) {

  # --- 1. Check reticulate availability ---
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop(
      "Package 'reticulate' is required for extract_cr_data().\n",
      "Install it with: install.packages(\"reticulate\")",
      call. = FALSE
    )
  }

  # --- 2. Configure Python ---
  if (!is.null(python)) {
    reticulate::use_python(python, required = TRUE)
  }

  # --- 3. macOS .d warning ---
  if (.Platform$OS.type != "windows" && grepl("\\.d$", input_path)) {
    warning(
      "On macOS/Linux, .d files cannot be read accurately ",
      "(Bruker DLL unavailable).\n",
      "Convert .d to .hdf on Windows first, then use the .hdf file here.\n",
      "See the CR Manual (Section 3.3) for instructions.",
      call. = FALSE
    )
  }

  # --- 4. Verify Python + alphatims ---
  .check_python_alphatims()

  # --- 5. Locate the bundled script ---
  script_path <- system.file("scripts", "alphatims", "extract_precursors.py",
                             package = "msProteomiX")
  if (script_path == "") {
    stop(
      "Cannot find extract_precursors.py in the installed package.\n",
      "Try reinstalling: remotes::install_github(\"BayOmics/msProteomiX\")",
      call. = FALSE
    )
  }

  # --- 6. Resolve paths ---
  input_path <- normalizePath(input_path, mustWork = TRUE)
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  output_dir <- normalizePath(output_dir)

  # --- 7. Build command ---
  args <- c(
    script_path,
    input_path,
    "--output_dir", output_dir,
    "--n_slices", as.character(n_slices),
    "--line_points", line_points
  )
  if (export_slices) {
    args <- c(args, "--export_slices")
  }
  if (!is.null(heatmap_rts)) {
    args <- c(args, "--heatmap_rts", paste(heatmap_rts, collapse = ","))
  }

  # --- 8. Run via reticulate ---
  message(">>> Running extract_precursors.py ...")
  message("    Input:  ", input_path)
  message("    Output: ", output_dir)
  message("    Slices: ", n_slices)

  # Use system2 with the Python executable found by reticulate
  py_exe <- reticulate::py_config()$python
  ret <- system2(py_exe, args = args, stdout = "", stderr = "")

  if (ret != 0) {
    stop(
      "Python script failed (exit code ", ret, ").\n",
      "Check the console output above for details.\n",
      "Common causes:\n",
      "  - alphatims not installed: pip install alphatims\n",
      "  - .d file on macOS: convert to .hdf on Windows first\n",
      "  - Invalid input path",
      call. = FALSE
    )
  }

  # --- 9. Report results ---
  csv_files <- list.files(output_dir, pattern = "cr_summary\\.csv$",
                          full.names = TRUE)
  message(">>> Done! ", length(csv_files), " summary file(s) in: ", output_dir)

  invisible(output_dir)
}


#' Check Python and alphatims availability
#'
#' @keywords internal
.check_python_alphatims <- function() {
  # Check Python is accessible
  tryCatch({
    py_ver <- reticulate::py_config()$version
    message("    Python: ", py_ver)
  }, error = function(e) {
    stop(
      "No Python installation found.\n\n",
      "Please install Python (3.8+). Recommended options:\n",
      "  1. Anaconda: https://www.anaconda.com/download\n",
      "  2. Python.org: https://www.python.org/downloads/\n\n",
      "After installing, restart RStudio and try again.",
      call. = FALSE
    )
  })

  # Check alphatims is installed
  has_alphatims <- reticulate::py_module_available("alphatims")
  if (!has_alphatims) {
    stop(
      "Python package 'alphatims' is not installed.\n\n",
      "Install it by running in a terminal:\n",
      "  pip install alphatims\n\n",
      "Or from R:\n",
      "  reticulate::py_install(\"alphatims\")\n\n",
      "Then restart RStudio and try again.",
      call. = FALSE
    )
  }

  at_ver <- tryCatch({
    at <- reticulate::import("alphatims")
    at$`__version__`
  }, error = function(e) "unknown")
  message("    alphatims: ", at_ver)
}
