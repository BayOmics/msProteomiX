# ==============================================================================
# msProteomiX — Plasma Contamination Index (CI) Calculation
# ==============================================================================
# References:
#   - Geyer PE et al., EMBO Mol Med 2019; 11(11): e10427  (CI formula + panels)
#   - Gao et al., EMBO Mol Med 2026 (NP-specific 30-protein panels)
#   - Korff et al., bioRxiv 2025 (PBMC panel)
# ==============================================================================


#' Load marker panel definitions
#'
#' @param panels "default" to use built-in panels, or path to custom CSV
#' @return data.frame with columns: Panel, Source, UniProt_ID, Gene, Description, Match_By
#' @keywords internal
.load_marker_panels <- function(panels = "default") {
  if (identical(panels, "default") || identical(panels, "all")) {
    csv_path <- system.file("extdata", "plasma_marker_panels.csv",
                            package = "msProteomiX")
    if (!nzchar(csv_path) || !file.exists(csv_path)) {
      stop("Built-in marker panel file not found. ",
           "Please reinstall msProteomiX or provide a custom CSV path.")
    }
  } else if (is.character(panels) && length(panels) == 1 && file.exists(panels)) {
    csv_path <- panels
  } else {
    stop("'panels' must be \"default\" or a path to a CSV file.")
  }
  df <- read.csv(csv_path, stringsAsFactors = FALSE)
  required <- c("Panel", "UniProt_ID", "Gene", "Match_By")
  missing <- setdiff(required, colnames(df))
  if (length(missing) > 0) {
    stop("Panel CSV missing required columns: ", paste(missing, collapse = ", "))
  }
  df
}


#' Match marker proteins by UniProt accession
#'
#' Handles semicolon-separated multi-accession entries in Protein ID column.
#'
#' @param ids Character vector of UniProt IDs to match
#' @param pinfo data.frame with "Protein ID" column
#' @return Integer vector of row indices
#' @keywords internal
.match_markers_uniprot <- function(ids, pinfo) {
  acc_col <- as.character(pinfo[["Protein ID"]])
  rows <- integer(0)
  for (i in seq_along(acc_col)) {
    accs <- trimws(unlist(strsplit(acc_col[i], ";")))
    if (any(accs %in% ids)) rows <- c(rows, i)
  }
  rows
}


#' Match marker proteins by gene name
#'
#' @param genes Character vector of gene names to match
#' @param pinfo data.frame with "Gene" column
#' @return Integer vector of row indices
#' @keywords internal
.match_markers_gene <- function(genes, pinfo) {
  gene_col <- as.character(pinfo[["Gene"]])
  rows <- integer(0)
  for (i in seq_along(gene_col)) {
    gs <- trimws(unlist(strsplit(gene_col[i], ";")))
    if (any(gs %in% genes)) rows <- c(rows, i)
  }
  rows
}


#' Calculate CI per sample
#'
#' CI = SUM(marker intensities) / SUM(all intensities)
#'
#' @param mat_raw Raw intensity matrix (non-log2)
#' @param marker_rows Integer vector of marker row indices
#' @return Numeric vector of CI values (one per sample/column)
#' @keywords internal
.calc_ci_per_sample <- function(mat_raw, marker_rows) {
  if (length(marker_rows) == 0) return(rep(NA_real_, ncol(mat_raw)))
  sapply(seq_len(ncol(mat_raw)), function(j) {
    vals <- mat_raw[, j]
    denom <- sum(vals, na.rm = TRUE)
    if (denom == 0) return(NA_real_)
    sum(vals[marker_rows], na.rm = TRUE) / denom
  })
}


#' Judge CI value as Pass/Flag/Fail
#'
#' @param value CI value
#' @param pass_threshold Upper bound for Pass
#' @param fail_threshold Lower bound for Fail
#' @return Character: "Pass", "Flag", "Fail", or "N/A"
#' @keywords internal
.judge_ci <- function(value, pass_threshold, fail_threshold) {
  if (is.na(value)) return("N/A")
  if (value <= pass_threshold) return("Pass")
  if (value > fail_threshold) return("Fail")
  "Flag"
}


#' Calculate Contamination Index (CI) for plasma proteomics
#'
#' Evaluates blood cell contamination in plasma proteomics samples using
#' marker protein panels. Based on the CI formula from Geyer et al. 2019:
#' CI = SUM(marker intensities) / SUM(all protein intensities).
#'
#' Default panels:
#' \itemize{
#'   \item \strong{PLT_Gao2026}: 30 platelet-specific proteins (Gao et al. 2026, NP-optimized)
#'   \item \strong{RBC_Gao2026}: 30 erythrocyte-specific proteins (Gao et al. 2026)
#'   \item \strong{PBMC_Korff2025}: 5 PBMC marker genes (Korff et al. 2025)
#' }
#'
#' Additionally, PLT_Geyer2019 (29 proteins) is computed as a reference
#' but does NOT participate in the Overall judgment.
#'
#' @param ms MsDataSet object from \code{parse_spectronaut()} or \code{parse_fragpipe()}.
#' @param panels Character: "default" to use built-in panels, or a path to a
#'   custom CSV file with columns: Panel, UniProt_ID, Gene, Match_By.
#' @param thresholds Named list of Pass/Fail thresholds per panel category.
#'   Default: \code{list(PLT = c(0.009, 0.05), RBC = c(0.05, 0.1), PBMC = c(0.005, 0.02))}.
#'   Each element is \code{c(pass_upper, fail_lower)}.
#' @param min_match Integer: minimum number of matched markers required to
#'   compute CI for a panel. Default 2 for PBMC (5-protein panel), 5 for others.
#'
#' @return A data.frame with one row per sample and columns:
#'   \describe{
#'     \item{Sample}{Sample name}
#'     \item{PGs}{Number of detected protein groups}
#'     \item{CI_PLT}{PLT contamination index (Gao2026 panel)}
#'     \item{CI_PLT_Geyer}{PLT contamination index (Geyer2019, reference only)}
#'     \item{CI_RBC}{RBC contamination index (Gao2026 panel)}
#'     \item{CI_PBMC}{PBMC contamination index (Korff2025 panel)}
#'     \item{J_PLT, J_RBC, J_PBMC}{Judgment: "Pass", "Flag", or "Fail"}
#'     \item{Overall}{"Pass", "Flag", or "Fail" (worst of PLT + RBC + PBMC)}
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ms <- parse_spectronaut("path/to/data")
#' ci <- calc_contamination_index(ms)
#' print(ci)
#' }
calc_contamination_index <- function(ms,
                                      panels = "default",
                                      thresholds = NULL,
                                      min_match = NULL) {
  # Validate input
  if (!inherits(ms, "MsDataSet")) {
    stop("'ms' must be an MsDataSet object from parse_spectronaut() or parse_fragpipe().")
  }

  # Default thresholds
  if (is.null(thresholds)) {
    thresholds <- list(
      PLT  = c(0.009, 0.05),
      RBC  = c(0.05,  0.1),
      PBMC = c(0.005, 0.02)
    )
  }

  # Default minimum match counts
  if (is.null(min_match)) {
    min_match <- list(PLT = 5, RBC = 5, PBMC = 2)
  } else if (is.numeric(min_match) && length(min_match) == 1) {
    min_match <- list(PLT = min_match, RBC = min_match, PBMC = min_match)
  }

  # Get raw intensity matrix (ensure matrix, not data.frame)
  prot_mat <- as.matrix(ms$proteins)
  if (isTRUE(ms$is_log2)) {
    message("    Data is log2-transformed. Converting to raw intensities for CI calculation.")
    prot_mat <- 2^prot_mat
  }
  prot_mat[is.na(prot_mat)] <- 0
  prot_info <- ms$protein_info
  samples <- colnames(prot_mat)

  # Load panels
  panel_df <- .load_marker_panels(panels)
  panel_names <- unique(panel_df$Panel)

  # --- Match markers ---
  message(">>> Matching marker panels ...")

  .match_panel <- function(panel_name) {
    sub_df <- panel_df[panel_df$Panel == panel_name, ]
    if (nrow(sub_df) == 0) return(list(rows = integer(0), total = 0L))
    match_by <- sub_df$Match_By[1]
    if (match_by == "gene") {
      genes <- sub_df$Gene[nzchar(sub_df$Gene)]
      rows <- if ("Gene" %in% colnames(prot_info) && length(genes) > 0) {
        .match_markers_gene(genes, prot_info)
      } else integer(0)
      # Fallback: if Gene matching failed, try UniProt ID
      if (length(rows) == 0) {
        ids <- sub_df$UniProt_ID[nzchar(sub_df$UniProt_ID)]
        if (length(ids) > 0 && "Protein ID" %in% colnames(prot_info)) {
          rows <- .match_markers_uniprot(ids, prot_info)
          if (length(rows) > 0) message(sprintf("      %s: Gene column missing/empty, using UniProt fallback", panel_name))
        }
      }
    } else {
      ids <- sub_df$UniProt_ID[nzchar(sub_df$UniProt_ID)]
      rows <- .match_markers_uniprot(ids, prot_info)
    }
    list(rows = rows, total = nrow(sub_df))
  }

  # Match default panels
  plt_gao   <- .match_panel("PLT_Gao2026")
  rbc_gao   <- .match_panel("RBC_Gao2026")
  pbmc      <- .match_panel("PBMC_Korff2025")
  plt_geyer <- .match_panel("PLT_Geyer2019")

  message(sprintf("    PLT  Gao2026:   %d/%d matched", length(plt_gao$rows), plt_gao$total))
  message(sprintf("    PLT  Geyer2019: %d/%d matched (reference)", length(plt_geyer$rows), plt_geyer$total))
  message(sprintf("    RBC  Gao2026:   %d/%d matched", length(rbc_gao$rows), rbc_gao$total))
  message(sprintf("    PBMC Korff2025: %d/%d matched", length(pbmc$rows), pbmc$total))

  # Check minimum match
  plt_ok  <- length(plt_gao$rows) >= min_match$PLT
  rbc_ok  <- length(rbc_gao$rows) >= min_match$RBC
  pbmc_ok <- length(pbmc$rows)    >= min_match$PBMC

  if (!plt_ok) message("    WARNING: PLT markers below minimum (", min_match$PLT, "). CI_PLT will be NA.")
  if (!rbc_ok) message("    WARNING: RBC markers below minimum (", min_match$RBC, "). CI_RBC will be NA.")
  if (!pbmc_ok) message("    WARNING: PBMC markers below minimum (", min_match$PBMC, "). CI_PBMC will be NA.")

  # --- Calculate CI ---
  ci_plt_gao   <- if (plt_ok)  .calc_ci_per_sample(prot_mat, plt_gao$rows)   else rep(NA_real_, length(samples))
  ci_rbc_gao   <- if (rbc_ok)  .calc_ci_per_sample(prot_mat, rbc_gao$rows)   else rep(NA_real_, length(samples))
  ci_pbmc      <- if (pbmc_ok) .calc_ci_per_sample(prot_mat, pbmc$rows)      else rep(NA_real_, length(samples))
  ci_plt_geyer <- .calc_ci_per_sample(prot_mat, plt_geyer$rows)

  # PGs per sample
  pgs <- sapply(seq_len(ncol(ms$proteins)), function(j) sum(!is.na(ms$proteins[, j])))

  # --- Build result ---
  ci_df <- data.frame(
    Sample       = samples,
    PGs          = pgs,
    CI_PLT       = ci_plt_gao,
    CI_PLT_Geyer = ci_plt_geyer,
    CI_RBC       = ci_rbc_gao,
    CI_PBMC      = ci_pbmc,
    stringsAsFactors = FALSE,
    row.names = NULL
  )

  # Judgments
  plt_thresh  <- thresholds$PLT
  rbc_thresh  <- thresholds$RBC
  pbmc_thresh <- thresholds$PBMC

  ci_df$J_PLT  <- vapply(ci_df$CI_PLT,  .judge_ci, character(1), plt_thresh[1],  plt_thresh[2])
  ci_df$J_RBC  <- vapply(ci_df$CI_RBC,  .judge_ci, character(1), rbc_thresh[1],  rbc_thresh[2])
  ci_df$J_PBMC <- vapply(ci_df$CI_PBMC, .judge_ci, character(1), pbmc_thresh[1], pbmc_thresh[2])

  # Overall: worst of PLT + RBC + PBMC (Geyer2019 excluded)
  ci_df$Overall <- vapply(seq_len(nrow(ci_df)), function(i) {
    js <- c(ci_df$J_PLT[i], ci_df$J_RBC[i], ci_df$J_PBMC[i])
    if (any(js == "Fail", na.rm = TRUE)) return("Fail")
    if (any(js == "Flag", na.rm = TRUE)) return("Flag")
    if (all(js == "N/A")) return("N/A")
    "Pass"
  }, character(1))

  # Report match info as attribute for downstream use
  attr(ci_df, "match_info") <- list(
    PLT_Gao2026   = list(matched = length(plt_gao$rows), total = plt_gao$total),
    PLT_Geyer2019 = list(matched = length(plt_geyer$rows), total = plt_geyer$total),
    RBC_Gao2026   = list(matched = length(rbc_gao$rows), total = rbc_gao$total),
    PBMC_Korff2025 = list(matched = length(pbmc$rows), total = pbmc$total)
  )

  message(">>> CI calculation complete.")
  ci_df
}
