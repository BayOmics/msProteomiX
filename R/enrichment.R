# ==============================================================================
# msProteomiX — GO 富集分析
# ==============================================================================

#' 运行 GO 富集分析 (Cellular Component)
#'
#' @param ms_data MsDataSet 对象 (需含 peptides 数据)
#' @param group_info 分组信息 data.frame
#' @param org_db 物种注释包名 (默认 "org.Hs.eg.db" 人类)
#' @param ont GO Ontology: "CC" (Cellular Component), "BP", "MF"
#' @param top_n 每组保留前 N 个结果
#' @return data.frame
#' @export
run_go_enrichment <- function(ms_data, group_info,
                               org_db = "org.Hs.eg.db",
                               ont = "CC",
                               top_n = 10) {
  stopifnot(inherits(ms_data, "MsDataSet"))
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Please install clusterProfiler: BiocManager::install('clusterProfiler')")
  }
  if (!requireNamespace(org_db, quietly = TRUE)) {
    stop(paste("Please install", org_db, ": BiocManager::install('", org_db, "')"))
  }

  pep_df <- ms_data$peptides
  if (nrow(pep_df) == 0) stop("No peptide data available.")

  gene_col <- grep("Mapped Gene|Gene Name|Gene", colnames(pep_df),
                    ignore.case = TRUE, value = TRUE)[1]
  if (is.na(gene_col)) stop("Cannot find Gene column in peptide data.")

  pep_raw_cols <- colnames(pep_df)
  unique_groups <- unique(group_info$user_group)
  gocc_results <- list()

  for (grp in unique_groups) {
    message(paste("  Analyzing group:", grp))
    grp_samples <- group_info$sample_name[group_info$user_group == grp]
    grp_cols <- c()
    for (s in grp_samples) {
      s_safe <- escape_regex(s)
      pat <- paste0("^", s_safe, ".*(Intensity|Spectral Count)$")
      hits <- grep(pat, pep_raw_cols, ignore.case = TRUE, value = TRUE)
      hits <- hits[!grepl("(Unique|Total)", hits, ignore.case = TRUE)]
      grp_cols <- c(grp_cols, hits)
    }
    grp_cols <- unique(grp_cols)
    if (length(grp_cols) == 0) next

    is_detected <- rowSums(pep_df[, grp_cols, drop = FALSE] > 0, na.rm = TRUE) > 0
    detected_genes_str <- pep_df[[gene_col]][is_detected]
    detected_genes <- unlist(strsplit(as.character(detected_genes_str), "[;\\|,]"))
    detected_genes <- unique(trimws(detected_genes))
    detected_genes <- detected_genes[detected_genes != "" & !is.na(detected_genes)]

    if (length(detected_genes) < 10) {
      message(paste("  Too few genes (<10) for group", grp, "- skipping."))
      next
    }

    ego <- tryCatch({
      clusterProfiler::enrichGO(
        gene          = detected_genes,
        OrgDb         = get(org_db),
        keyType       = "SYMBOL",
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2
      )
    }, error = function(e) NULL)

    if (!is.null(ego) && nrow(ego) > 0) {
      top_res <- ego@result %>%
        dplyr::arrange(p.adjust) %>%
        utils::head(top_n)
      top_res$Group <- grp
      gocc_results[[grp]] <- top_res
    }
  }

  if (length(gocc_results) == 0) {
    message("  No significant GO enrichment results.")
    return(data.frame())
  }
  do.call(rbind, gocc_results)
}


#' 绘制 GO 富集气泡图
#'
#' @param go_df run_go_enrichment() 返回的 data.frame
#' @param output_dir 输出目录
#' @param project_name 项目名
#' @return ggplot 对象
#' @export
plot_go_bubble <- function(go_df,
                            output_dir = "output",
                            project_name = "Project") {
  ensure_output_dir(output_dir)

  if (is.null(go_df) || nrow(go_df) == 0) {
    message("  No GO data to plot.")
    return(invisible(NULL))
  }

  # 自动换行
  go_df$Description <- stringr::str_wrap(go_df$Description, width = 50)
  go_df$Description <- factor(go_df$Description,
                                levels = unique(go_df$Description[order(go_df$Count)]))

  p <- ggplot2::ggplot(go_df, ggplot2::aes(x = Count, y = Description)) +
    ggplot2::geom_point(ggplot2::aes(size = Count, color = p.adjust)) +
    ggplot2::scale_color_gradient(low = "red", high = "blue") +
    ggplot2::facet_grid(Group ~ ., scales = "free_y", space = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "GO Cellular Component Enrichment",
                  x = "Gene Count", y = NULL,
                  color = "p.adjust", size = "Count") +
    ggplot2::theme(
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      axis.text.y = ggplot2::element_text(size = 9)
    )

  plot_height <- 2 + (nrow(go_df) * 0.35)
  if (plot_height < 5) plot_height <- 5

  save_plot_and_data(p, go_df, project_name, "GO_Enrichment",
                     output_dir = output_dir, height = plot_height, width = 9)
  p
}
