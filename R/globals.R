# ==============================================================================
# msProteomiX - Global variable declarations (suppress R CMD check NOTEs)
# ==============================================================================
# These variables are used in dplyr/ggplot2 NSE (non-standard evaluation)
# contexts and are not actual global variables.

utils::globalVariables(c(
  # visualization_barplot.R
  "mean_val", "sd_val",
  # visualization_pca.R
  "PC1", "PC2", "Group", "ShortLabel",
  # visualization_volcano.R
  "logFC",
  # visualization_apms.R
  "Value", "DisplayLabel", "Label", "MaxGroup", "Mean", "SD",
  "RelAbundance", "ClusterLabel", "Timepoint",
  # visualization_cv.R
  "CV", "MedianCV", "group", "label",
  # visualization_correlation.R
  "Sample1", "Sample2", "R2",
  # coverage.R
  "Coverage", "MedianCov",
  # enrichment.R
  "Description", "Count",
  # diff_analysis.R
  "diff",
  # grouping.R / general
  "user_group", "sample_name"
))
