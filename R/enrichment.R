# ==============================================================================
# msProteomiX — GO 富集分析
# ==============================================================================

#' 运行 GO 富集分析 (Cellular Component)
#'
#' @param ms_data MsDataSet 对象 (使用 protein_info 和 proteins 矩阵提取基因)
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

  # 使用 protein_info + protein matrix 提取各组基因 (兼容所有引擎)
  gene_list_by_group <- .extract_group_genes(ms_data, group_info)
  if (length(gene_list_by_group) == 0) {
    message("  No gene data available for GO enrichment.")
    return(data.frame())
  }

  gocc_results <- list()
  for (grp in names(gene_list_by_group)) {
    message(paste("  GO enrichment for group:", grp))
    detected_genes <- gene_list_by_group[[grp]]

    if (length(detected_genes) < 10) {
      message(paste("  Too few genes (<10) for group", grp, "- skipping."))
      next
    }

    ego <- tryCatch({
      clusterProfiler::enrichGO(
        gene          = detected_genes,
        OrgDb         = loadNamespace(org_db)[[org_db]],
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


#' Plot GO enrichment bar chart
#'
#' Horizontal bar chart with facet by group, complementary to plot_go_bubble().
#'
#' @param go_df run_go_enrichment() result data.frame
#' @param output_dir Output directory
#' @param project_name Project name
#' @return ggplot object
#' @export
plot_go_bar <- function(go_df,
                        output_dir = "output",
                        project_name = "Project") {
  ensure_output_dir(output_dir)

  if (is.null(go_df) || nrow(go_df) == 0) {
    message("  No GO data to plot.")
    return(invisible(NULL))
  }

  go_df$Description <- stringr::str_trunc(go_df$Description, 40)
  go_df$Description <- factor(go_df$Description,
                              levels = unique(go_df$Description[order(go_df$Count)]))

  p <- ggplot2::ggplot(go_df,
                        ggplot2::aes(x = stats::reorder(Description, Count),
                                    y = Count, fill = p.adjust)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::scale_fill_gradient(low = "red", high = "blue") +
    ggplot2::facet_wrap(~ Group, scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "GO Cellular Component Enrichment",
                  x = NULL, y = "Gene Count",
                  fill = "p.adjust") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      axis.text.y = ggplot2::element_text(size = 9)
    )

  n_groups <- length(unique(go_df$Group))
  plot_height <- 6 + n_groups * 2

  save_plot_and_data(p, go_df, project_name, "GO_Enrichment_Bar",
                     output_dir = output_dir, height = plot_height, width = 8)
  p
}


# ==============================================================================
# KEGG Pathway Enrichment
# ==============================================================================

#' Run KEGG pathway enrichment analysis (offline)
#'
#' Uses local org.db annotation database for KEGG pathway-gene mappings.
#' No internet connection required.
#'
#' @param ms_data MsDataSet object (used when diff_result is NULL)
#' @param group_info Group info data.frame (used when diff_result is NULL)
#' @param diff_result Optional. Output from run_diff_analysis(). If provided,
#'   enrichment is performed on differentially expressed proteins only.
#' @param org_db Annotation package name (default "org.Hs.eg.db")
#' @param organism KEGG organism code (default "hsa" for human)
#' @param top_n Number of top results to keep per group/set
#' @return data.frame with enrichment results
#' @export
run_kegg_enrichment <- function(ms_data = NULL, group_info = NULL,
                                diff_result = NULL,
                                org_db = "org.Hs.eg.db",
                                organism = "hsa",
                                top_n = 15) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Please install clusterProfiler: BiocManager::install('clusterProfiler')")
  }
  if (!requireNamespace(org_db, quietly = TRUE)) {
    stop(paste("Please install", org_db, ": BiocManager::install('", org_db, "')"))
  }

  # Build local KEGG TERM2GENE from org.db
  org_obj <- loadNamespace(org_db)[[org_db]]
  kegg_keys <- tryCatch(
    AnnotationDbi::keys(org_obj, keytype = "PATH"),
    error = function(e) character(0)
  )

  if (length(kegg_keys) == 0) {
    message(">>> No KEGG pathway data found in ", org_db, ". Check annotation package.")
    return(data.frame())
  }

  message(sprintf(">>> Loading KEGG pathways from local %s (%d pathways, offline)...",
                  org_db, length(kegg_keys)))
  path_data <- AnnotationDbi::select(org_obj,
                                      keys = kegg_keys,
                                      columns = c("ENTREZID", "PATH"),
                                      keytype = "PATH")
  path_data <- path_data[!is.na(path_data$ENTREZID) & !is.na(path_data$PATH), ]

  # Format TERM2GENE: add organism prefix for readability
  prefix <- paste0(organism, ":")
  term2gene <- data.frame(
    term = paste0(prefix, path_data$PATH),
    gene = path_data$ENTREZID,
    stringsAsFactors = FALSE
  )

  # Build TERM2NAME if possible (pathway ID -> description)
  term2name <- NULL

  # Helper function for single enrichment
  .do_kegg_enricher <- function(entrez_ids, term2gene, term2name, top_n) {
    ekegg <- tryCatch({
      clusterProfiler::enricher(
        gene          = entrez_ids,
        TERM2GENE     = term2gene,
        TERM2NAME     = term2name,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2
      )
    }, error = function(e) {
      message(">>> KEGG enrichment failed: ", e$message)
      NULL
    })

    if (is.null(ekegg) || nrow(ekegg) == 0) return(NULL)

    ekegg@result %>%
      dplyr::arrange(p.adjust) %>%
      utils::head(top_n)
  }

  if (!is.null(diff_result)) {
    # --- Mode 1: Enrichment on differentially expressed proteins ---
    gene_list <- .extract_diff_genes(diff_result)
    if (length(gene_list) < 5) {
      message(">>> Too few differentially expressed genes (<5) for KEGG enrichment.")
      return(data.frame())
    }

    entrez_ids <- .symbol_to_entrez(gene_list, org_db)
    if (length(entrez_ids) < 5) {
      message(">>> Too few mapped Entrez IDs (<5). Check gene symbols and org_db.")
      return(data.frame())
    }

    result <- .do_kegg_enricher(entrez_ids, term2gene, term2name, top_n)
    if (is.null(result)) {
      message(">>> No significant KEGG enrichment results.")
      return(data.frame())
    }
    result$Group <- attr(diff_result, "contrast") %||% "DiffExpr"
    return(result)

  } else {
    # --- Mode 2: Enrichment by group ---
    if (is.null(ms_data) || is.null(group_info)) {
      stop("Provide either diff_result, or both ms_data and group_info.")
    }
    stopifnot(inherits(ms_data, "MsDataSet"))

    gene_list_by_group <- .extract_group_genes(ms_data, group_info)
    kegg_results <- list()

    for (grp in names(gene_list_by_group)) {
      genes <- gene_list_by_group[[grp]]
      if (length(genes) < 10) {
        message(paste("  Too few genes (<10) for group", grp, "- skipping."))
        next
      }

      entrez_ids <- .symbol_to_entrez(genes, org_db)
      if (length(entrez_ids) < 5) next

      message(paste("  KEGG enrichment for group:", grp))
      top_res <- .do_kegg_enricher(entrez_ids, term2gene, term2name, top_n)
      if (!is.null(top_res)) {
        top_res$Group <- grp
        kegg_results[[grp]] <- top_res
      }
    }

    if (length(kegg_results) == 0) {
      message(">>> No significant KEGG enrichment results.")
      return(data.frame())
    }
    do.call(rbind, kegg_results)
  }
}


#' Plot KEGG enrichment bar chart
#'
#' @param kegg_df run_kegg_enrichment() result data.frame
#' @param output_dir Output directory
#' @param project_name Project name
#' @return ggplot object
#' @export
plot_kegg_bar <- function(kegg_df,
                          output_dir = "output",
                          project_name = "Project") {
  ensure_output_dir(output_dir)

  if (is.null(kegg_df) || nrow(kegg_df) == 0) {
    message("  No KEGG data to plot.")
    return(invisible(NULL))
  }

  kegg_df$Description <- stringr::str_wrap(kegg_df$Description, width = 50)
  kegg_df$Description <- factor(kegg_df$Description,
                                levels = unique(kegg_df$Description[order(kegg_df$Count)]))

  p <- ggplot2::ggplot(kegg_df, ggplot2::aes(x = Count, y = Description)) +
    ggplot2::geom_point(ggplot2::aes(size = Count, color = p.adjust)) +
    ggplot2::scale_color_gradient(low = "#E64B35", high = "#4DBBD5") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "KEGG Pathway Enrichment",
                  x = "Gene Count", y = NULL,
                  color = "p.adjust", size = "Count") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      axis.text.y = ggplot2::element_text(size = 9)
    )

  # Add facet if multiple groups
  if (length(unique(kegg_df$Group)) > 1) {
    p <- p + ggplot2::facet_grid(Group ~ ., scales = "free_y", space = "free_y")
  }

  plot_height <- 2 + (nrow(kegg_df) * 0.35)
  if (plot_height < 5) plot_height <- 5

  save_plot_and_data(p, kegg_df, project_name, "KEGG_Enrichment",
                     output_dir = output_dir, height = plot_height, width = 10)
  p
}


# ==============================================================================
# Reactome Pathway Enrichment
# ==============================================================================

#' Run Reactome pathway enrichment analysis
#'
#' @param ms_data MsDataSet object (used when diff_result is NULL)
#' @param group_info Group info data.frame (used when diff_result is NULL)
#' @param diff_result Optional. Output from run_diff_analysis().
#' @param org_db Annotation package name (default "org.Hs.eg.db")
#' @param organism Reactome organism (default "human")
#' @param top_n Number of top results to keep
#' @return data.frame with enrichment results
#' @export
run_reactome_enrichment <- function(ms_data = NULL, group_info = NULL,
                                    diff_result = NULL,
                                    org_db = "org.Hs.eg.db",
                                    organism = "human",
                                    top_n = 20) {
  if (!requireNamespace("ReactomePA", quietly = TRUE)) {
    stop("Please install ReactomePA: BiocManager::install('ReactomePA')")
  }
  if (!requireNamespace("reactome.db", quietly = TRUE)) {
    message(">>> NOTE: reactome.db required. Install: BiocManager::install('reactome.db')")
  }
  if (!requireNamespace(org_db, quietly = TRUE)) {
    stop(paste("Please install", org_db, ": BiocManager::install('", org_db, "')"))
  }

  if (!is.null(diff_result)) {
    # --- Mode 1: Differentially expressed proteins ---
    gene_list <- .extract_diff_genes(diff_result)
    if (length(gene_list) < 5) {
      message(">>> Too few DE genes (<5) for Reactome enrichment.")
      return(data.frame())
    }

    entrez_ids <- .symbol_to_entrez(gene_list, org_db)
    if (length(entrez_ids) < 5) {
      message(">>> Too few mapped Entrez IDs (<5).")
      return(data.frame())
    }

    epa <- tryCatch({
      ReactomePA::enrichPathway(
        gene         = entrez_ids,
        organism     = organism,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2
      )
    }, error = function(e) {
      message(">>> Reactome enrichment failed: ", e$message)
      NULL
    })

    if (is.null(epa) || nrow(epa) == 0) {
      message(">>> No significant Reactome enrichment results.")
      return(data.frame())
    }

    result <- epa@result %>%
      dplyr::arrange(p.adjust) %>%
      utils::head(top_n)
    result$Group <- attr(diff_result, "contrast") %||% "DiffExpr"
    return(result)

  } else {
    # --- Mode 2: By group ---
    if (is.null(ms_data) || is.null(group_info)) {
      stop("Provide either diff_result, or both ms_data and group_info.")
    }
    stopifnot(inherits(ms_data, "MsDataSet"))

    gene_list_by_group <- .extract_group_genes(ms_data, group_info)
    react_results <- list()

    for (grp in names(gene_list_by_group)) {
      genes <- gene_list_by_group[[grp]]
      if (length(genes) < 10) {
        message(paste("  Too few genes (<10) for group", grp, "- skipping."))
        next
      }

      entrez_ids <- .symbol_to_entrez(genes, org_db)
      if (length(entrez_ids) < 5) next

      message(paste("  Reactome enrichment for group:", grp))
      epa <- tryCatch({
        ReactomePA::enrichPathway(
          gene         = entrez_ids,
          organism     = organism,
          pAdjustMethod = "BH",
          pvalueCutoff = 0.05,
          qvalueCutoff = 0.2
        )
      }, error = function(e) NULL)

      if (!is.null(epa) && nrow(epa) > 0) {
        top_res <- epa@result %>%
          dplyr::arrange(p.adjust) %>%
          utils::head(top_n)
        top_res$Group <- grp
        react_results[[grp]] <- top_res
      }
    }

    if (length(react_results) == 0) {
      message(">>> No significant Reactome enrichment results.")
      return(data.frame())
    }
    do.call(rbind, react_results)
  }
}


#' Plot Reactome enrichment bar chart
#'
#' @param reactome_df run_reactome_enrichment() result data.frame
#' @param output_dir Output directory
#' @param project_name Project name
#' @return ggplot object
#' @export
plot_reactome_bar <- function(reactome_df,
                              output_dir = "output",
                              project_name = "Project") {
  ensure_output_dir(output_dir)

  if (is.null(reactome_df) || nrow(reactome_df) == 0) {
    message("  No Reactome data to plot.")
    return(invisible(NULL))
  }

  reactome_df$Description <- stringr::str_wrap(reactome_df$Description, width = 55)
  reactome_df$Description <- factor(reactome_df$Description,
                                    levels = unique(reactome_df$Description[
                                      order(reactome_df$Count)]))

  p <- ggplot2::ggplot(reactome_df,
                        ggplot2::aes(x = Count, y = Description)) +
    ggplot2::geom_point(ggplot2::aes(size = Count, color = p.adjust)) +
    ggplot2::scale_color_gradient(low = "#7E6148", high = "#B09C85") +
    ggplot2::theme_bw() +
    ggplot2::labs(title = "Reactome Pathway Enrichment",
                  x = "Gene Count", y = NULL,
                  color = "p.adjust", size = "Count") +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      strip.background = ggplot2::element_rect(fill = "grey90"),
      strip.text = ggplot2::element_text(face = "bold", size = 10),
      axis.text.y = ggplot2::element_text(size = 9)
    )

  if (length(unique(reactome_df$Group)) > 1) {
    p <- p + ggplot2::facet_grid(Group ~ ., scales = "free_y", space = "free_y")
  }

  plot_height <- 2 + (nrow(reactome_df) * 0.35)
  if (plot_height < 5) plot_height <- 5

  save_plot_and_data(p, reactome_df, project_name, "Reactome_Enrichment",
                     output_dir = output_dir, height = plot_height, width = 11)
  p
}


# ==============================================================================
# Shared internal helpers for enrichment
# ==============================================================================

#' Extract gene symbols from differential analysis result
#' @keywords internal
.extract_diff_genes <- function(diff_result) {
  # Get genes from UP and DOWN proteins
  diff_df <- diff_result[diff_result$diff %in% c("UP", "DOWN"), ]
  if (nrow(diff_df) == 0) return(character())

  # Use Gene column for enrichment (Label_Name = UniProt entry name, not gene symbol)
  gene_col <- if ("Gene" %in% colnames(diff_df)) "Gene"
              else if ("Label_Name" %in% colnames(diff_df)) "Label_Name"
              else NULL

  if (is.null(gene_col)) return(character())

  genes <- as.character(diff_df[[gene_col]])
  genes <- unlist(strsplit(genes, "[;|,]"))
  genes <- unique(trimws(genes))
  genes[genes != "" & !is.na(genes)]
}


#' Extract gene symbols from MsDataSet by group
#' @keywords internal
.extract_group_genes <- function(ms_data, group_info) {
  prot_mat <- ms_data$proteins
  pinfo <- ms_data$protein_info

  # Use Gene column for enrichment (Label_Name = UniProt entry name, not gene symbol)
  gene_col <- if ("Gene" %in% colnames(pinfo)) "Gene"
              else if ("Label_Name" %in% colnames(pinfo)) "Label_Name"
              else NULL

  if (is.null(gene_col)) {
    message(">>> Cannot find gene column in protein_info.")
    return(list())
  }

  unique_groups <- unique(group_info$user_group)
  result <- list()

  for (grp in unique_groups) {
    grp_samples <- group_info$sample_name[group_info$user_group == grp]
    grp_samples <- intersect(grp_samples, colnames(prot_mat))
    if (length(grp_samples) == 0) next

    is_detected <- rowSums(prot_mat[, grp_samples, drop = FALSE] > 0,
                           na.rm = TRUE) > 0
    detected_genes_str <- pinfo[[gene_col]][is_detected]
    detected_genes <- unlist(strsplit(as.character(detected_genes_str),
                                     "[;|,]"))
    detected_genes <- unique(trimws(detected_genes))
    detected_genes <- detected_genes[detected_genes != "" &
                                     !is.na(detected_genes)]
    result[[grp]] <- detected_genes
  }
  result
}


#' Convert gene symbols to Entrez IDs
#' @keywords internal
.symbol_to_entrez <- function(gene_symbols, org_db) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) return(character())

  id_map <- tryCatch({
    clusterProfiler::bitr(
      gene_symbols,
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = loadNamespace(org_db)[[org_db]]
    )
  }, error = function(e) {
    message(">>> Gene ID conversion failed: ", e$message)
    data.frame(SYMBOL = character(), ENTREZID = character())
  })

  if (nrow(id_map) == 0) return(character())
  unique(id_map$ENTREZID)
}


# ==============================================================================
# GSEA (Gene Set Enrichment Analysis)
# ==============================================================================

#' Run GSEA pre-ranked analysis (offline by default)
#'
#' Uses all proteins' logFC from diff_result (not just significant ones)
#' to perform gene set enrichment analysis.
#'
#' By default uses msigdbr package for offline gene sets (no internet needed).
#' Supported gene_sets: "kegg" (default), "hallmark", "reactome", "go_bp",
#' or "kegg_online" (requires internet via gseKEGG).
#'
#' @param diff_result data.frame from run_diff_analysis()
#' @param org_db Annotation database (default "org.Hs.eg.db")
#' @param organism KEGG organism code (default "hsa") or msigdbr species
#' @param gene_sets Gene set collection: "kegg" (offline), "hallmark",
#'   "reactome", "go_bp", or "kegg_online" (requires internet)
#' @param pvalue_cutoff P-value cutoff (default 0.05)
#' @param min_gs_size Minimum gene set size (default 10)
#' @param max_gs_size Maximum gene set size (default 500)
#' @return gseaResult object or NULL
#' @export
run_gsea_analysis <- function(diff_result,
                               org_db = "org.Hs.eg.db",
                               organism = "hsa",
                               gene_sets = "kegg",
                               pvalue_cutoff = 0.05,
                               min_gs_size = 10,
                               max_gs_size = 500) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("Please install clusterProfiler: BiocManager::install('clusterProfiler')")
  }

  # --- Build ranked gene list ---
  gene_col <- if ("Gene" %in% colnames(diff_result)) "Gene"
              else if ("Label_Name" %in% colnames(diff_result)) "Label_Name"
              else NULL
  if (is.null(gene_col)) {
    message(">>> Cannot find Gene column in diff_result.")
    return(NULL)
  }

  gene_symbols <- as.character(diff_result[[gene_col]])
  logfc_values <- diff_result$logFC

  valid <- !is.na(gene_symbols) & gene_symbols != "" & !is.na(logfc_values)
  gene_symbols <- gene_symbols[valid]
  logfc_values <- logfc_values[valid]
  first_gene <- sub(";.*$", "", gene_symbols)

  # --- Online gseKEGG path (requires internet) ---
  if (gene_sets == "kegg_online") {
    return(.run_gsea_online(diff_result, first_gene, logfc_values,
                            org_db, organism, pvalue_cutoff,
                            min_gs_size, max_gs_size))
  }

  # --- Offline path: gseGO (truly offline, uses local org.db) ---
  if (gene_sets %in% c("go_bp", "go_cc", "go_mf")) {
    if (!requireNamespace(org_db, quietly = TRUE)) {
      stop(paste("Please install", org_db, ": BiocManager::install('", org_db, "')"))
    }

    ont <- switch(gene_sets,
      "go_bp" = "BP", "go_cc" = "CC", "go_mf" = "MF", "BP"
    )

    message(sprintf(">>> Converting gene symbols to Entrez IDs..."))
    id_map <- tryCatch({
      clusterProfiler::bitr(
        unique(first_gene),
        fromType = "SYMBOL",
        toType   = "ENTREZID",
        OrgDb    = loadNamespace(org_db)[[org_db]]
      )
    }, error = function(e) {
      message(">>> Gene ID conversion failed: ", e$message)
      return(NULL)
    })

    if (is.null(id_map) || nrow(id_map) == 0) {
      message(">>> No genes could be mapped. GSEA aborted.")
      return(NULL)
    }

    gene_df <- data.frame(SYMBOL = first_gene, logFC = logfc_values,
                           stringsAsFactors = FALSE)
    merged <- merge(gene_df, id_map, by = "SYMBOL")
    merged <- merged[order(abs(merged$logFC), decreasing = TRUE), ]
    merged <- merged[!duplicated(merged$ENTREZID), ]

    gene_list <- merged$logFC
    names(gene_list) <- merged$ENTREZID
    gene_list <- sort(gene_list, decreasing = TRUE)

    if (length(gene_list) < 15) {
      message(sprintf(">>> Too few mapped genes (%d).", length(gene_list)))
      return(NULL)
    }

    message(sprintf(">>> Running gseGO (%s) with %d ranked genes (offline)...",
                    ont, length(gene_list)))

    gsea_result <- tryCatch({
      clusterProfiler::gseGO(
        geneList      = gene_list,
        OrgDb         = loadNamespace(org_db)[[org_db]],
        ont           = ont,
        minGSSize     = min_gs_size,
        maxGSSize     = max_gs_size,
        pvalueCutoff  = pvalue_cutoff,
        pAdjustMethod = "BH",
        verbose       = FALSE
      )
    }, error = function(e) {
      message(">>> gseGO failed: ", e$message)
      return(NULL)
    })

    if (is.null(gsea_result) || nrow(gsea_result) == 0) {
      message(">>> No significant GSEA results found.")
      return(gsea_result)
    }

    contrast <- attr(diff_result, "contrast")
    if (!is.null(contrast)) attr(gsea_result, "contrast") <- contrast

    n_up <- sum(gsea_result@result$NES > 0)
    n_down <- sum(gsea_result@result$NES < 0)
    message(sprintf(">>> GSEA complete: %d enriched (%d activated, %d suppressed)",
                    nrow(gsea_result), n_up, n_down))
    return(gsea_result)
  }

  # --- msigdbr path (kegg/hallmark/reactome, offline after first cache) ---
  if (gene_sets %in% c("kegg", "hallmark", "reactome")) {
    if (!requireNamespace("msigdbr", quietly = TRUE)) {
      stop("Please install msigdbr for ", gene_sets, " GSEA: install.packages('msigdbr')\n",
           "Or use gene_sets = 'go_bp' (fully offline) or 'kegg_online' (requires internet).")
    }

    species_map <- c("hsa" = "Homo sapiens", "mmu" = "Mus musculus",
                     "rno" = "Rattus norvegicus")
    species <- species_map[organism]
    if (is.na(species)) species <- "Homo sapiens"

    gs_config <- switch(gene_sets,
      "kegg"     = list(collection = "C2", subcollection = "CP:KEGG_MEDICUS",
                         label = "KEGG"),
      "hallmark" = list(collection = "H",  subcollection = NULL,
                         label = "Hallmark"),
      "reactome" = list(collection = "C2", subcollection = "CP:REACTOME",
                         label = "Reactome")
    )

    message(sprintf(">>> Loading %s gene sets (via msigdbr)...", gs_config$label))

    msig_df <- tryCatch({
      if (!is.null(gs_config$subcollection)) {
        msigdbr::msigdbr(species = species,
                         collection = gs_config$collection,
                         subcollection = gs_config$subcollection)
      } else {
        msigdbr::msigdbr(species = species,
                         collection = gs_config$collection)
      }
    }, error = function(e) {
      message(">>> msigdbr failed: ", e$message)
      message(">>> Try gene_sets='go_bp' (fully offline) or gene_sets='kegg_online'.")
      return(NULL)
    })

    if (is.null(msig_df) || nrow(msig_df) == 0) {
      message(">>> No gene sets found. Check organism/gene_sets settings.")
      return(NULL)
    }

    term2gene <- data.frame(
      gs_name     = msig_df$gs_name,
      gene_symbol = msig_df$gene_symbol,
      stringsAsFactors = FALSE
    )

    gene_df <- data.frame(SYMBOL = first_gene, logFC = logfc_values,
                           stringsAsFactors = FALSE)
    gene_df <- gene_df[order(abs(gene_df$logFC), decreasing = TRUE), ]
    gene_df <- gene_df[!duplicated(gene_df$SYMBOL), ]

    gene_list <- gene_df$logFC
    names(gene_list) <- gene_df$SYMBOL
    gene_list <- sort(gene_list, decreasing = TRUE)

    if (length(gene_list) < 15) {
      message(sprintf(">>> Too few genes (%d).", length(gene_list)))
      return(NULL)
    }

    message(sprintf(">>> Running GSEA (%s) with %d ranked genes...",
                    gs_config$label, length(gene_list)))

    gsea_result <- tryCatch({
      clusterProfiler::GSEA(
        geneList      = gene_list,
        TERM2GENE     = term2gene,
        minGSSize     = min_gs_size,
        maxGSSize     = max_gs_size,
        pvalueCutoff  = pvalue_cutoff,
        pAdjustMethod = "BH",
        verbose       = FALSE
      )
    }, error = function(e) {
      message(">>> GSEA failed: ", e$message)
      return(NULL)
    })

    if (is.null(gsea_result) || nrow(gsea_result) == 0) {
      message(">>> No significant GSEA results found.")
      return(gsea_result)
    }

    # Clean up pathway names
    gsea_result@result$Description <- gsub(
      "^KEGG_MEDICUS_|^KEGG_LEGACY_|^KEGG_|^HALLMARK_|^REACTOME_|^GOBP_", "",
      gsea_result@result$Description)
    gsea_result@result$Description <- gsub("_", " ", gsea_result@result$Description)
  } else {
    message(sprintf(">>> Unknown gene_sets: '%s'. Use 'go_bp', 'kegg', 'hallmark', 'reactome', or 'kegg_online'.",
                    gene_sets))
    return(NULL)
  }

  contrast <- attr(diff_result, "contrast")
  if (!is.null(contrast)) attr(gsea_result, "contrast") <- contrast

  n_up <- sum(gsea_result@result$NES > 0)
  n_down <- sum(gsea_result@result$NES < 0)
  message(sprintf(">>> GSEA complete: %d enriched (%d activated, %d suppressed)",
                  nrow(gsea_result), n_up, n_down))
  gsea_result
}


#' Run GSEA online via gseKEGG (internal helper)
#' @keywords internal
.run_gsea_online <- function(diff_result, first_gene, logfc_values,
                              org_db, organism, pvalue_cutoff,
                              min_gs_size, max_gs_size) {
  if (!requireNamespace(org_db, quietly = TRUE)) {
    stop(paste("Please install", org_db))
  }

  message(">>> Converting gene symbols to Entrez IDs...")
  id_map <- tryCatch({
    clusterProfiler::bitr(
      unique(first_gene),
      fromType = "SYMBOL",
      toType   = "ENTREZID",
      OrgDb    = loadNamespace(org_db)[[org_db]]
    )
  }, error = function(e) {
    message(">>> Gene ID conversion failed: ", e$message)
    return(NULL)
  })

  if (is.null(id_map) || nrow(id_map) == 0) {
    message(">>> No genes could be mapped. GSEA aborted.")
    return(NULL)
  }

  gene_df <- data.frame(SYMBOL = first_gene, logFC = logfc_values,
                         stringsAsFactors = FALSE)
  merged <- merge(gene_df, id_map, by = "SYMBOL")
  merged <- merged[order(abs(merged$logFC), decreasing = TRUE), ]
  merged <- merged[!duplicated(merged$ENTREZID), ]

  gene_list <- merged$logFC
  names(gene_list) <- merged$ENTREZID
  gene_list <- sort(gene_list, decreasing = TRUE)

  if (length(gene_list) < 15) {
    message(sprintf(">>> Too few mapped genes (%d).", length(gene_list)))
    return(NULL)
  }

  message(sprintf(">>> Running gseKEGG with %d ranked genes (online)...",
                  length(gene_list)))

  gsea_result <- tryCatch({
    clusterProfiler::gseKEGG(
      geneList      = gene_list,
      organism      = organism,
      minGSSize     = min_gs_size,
      maxGSSize     = max_gs_size,
      pvalueCutoff  = pvalue_cutoff,
      pAdjustMethod = "BH",
      verbose       = FALSE
    )
  }, error = function(e) {
    message(">>> gseKEGG failed: ", e$message)
    return(NULL)
  })

  if (is.null(gsea_result) || nrow(gsea_result) == 0) {
    message(">>> No significant GSEA results found.")
    return(gsea_result)
  }

  contrast <- attr(diff_result, "contrast")
  if (!is.null(contrast)) attr(gsea_result, "contrast") <- contrast

  n_up <- sum(gsea_result@result$NES > 0)
  n_down <- sum(gsea_result@result$NES < 0)
  message(sprintf(">>> GSEA complete: %d enriched (%d activated, %d suppressed)",
                  nrow(gsea_result), n_up, n_down))
  gsea_result
}


#' Plot GSEA results
#'
#' Generates a dot plot of GSEA results.
#'
#' @param gsea_result gseaResult object from run_gsea_analysis()
#' @param top_n Number of top pathways to show (default 20)
#' @param output_dir Output directory
#' @param project_name Project name
#' @return ggplot object or NULL
#' @export
plot_gsea_result <- function(gsea_result,
                              top_n = 20,
                              output_dir = "output",
                              project_name = "Project") {
  if (is.null(gsea_result) || nrow(gsea_result) == 0) {
    message(">>> No GSEA data to plot.")
    return(invisible(NULL))
  }
  ensure_output_dir(output_dir)

  res_df <- gsea_result@result
  res_df <- res_df[order(res_df$p.adjust), ]
  res_df <- utils::head(res_df, top_n)

  # Dot plot: x = NES, y = pathway, size = setSize, color = p.adjust
  res_df$Description <- factor(res_df$Description,
                                levels = rev(res_df$Description))

  p <- ggplot2::ggplot(res_df, ggplot2::aes(
    x = NES, y = Description, size = setSize, color = p.adjust
  )) +
    ggplot2::geom_point() +
    ggplot2::scale_color_gradient(low = "#E74C3C", high = "#3498DB",
                                  name = "p.adjust") +
    ggplot2::scale_size_continuous(range = c(3, 8), name = "Gene Set Size") +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::labs(
      title = "GSEA - KEGG Pathways",
      x = "Normalized Enrichment Score (NES)",
      y = NULL
    ) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold"),
      axis.text.y = ggplot2::element_text(size = 10)
    )

  contrast <- attr(gsea_result, "contrast")
  suffix <- if (!is.null(contrast)) paste0("GSEA_KEGG_", contrast) else "GSEA_KEGG"

  save_plot_and_data(p, res_df, project_name, suffix,
                     output_dir = output_dir, width = 12, height = 8)
  p
}


#' Plot classic GSEA enrichment running score plots
#'
#' Generates the classic GSEA visualization with:
#' - Running enrichment score curve
#' - Gene hit barcode positions
#' - Ranked list metric
#'
#' @param gsea_result gseaResult object from run_gsea_analysis()
#' @param top_n Number of top pathways to plot (default 3, each activated and suppressed)
#' @param output_dir Output directory
#' @param project_name Project name
#' @return List of ggplot objects (invisible)
#' @export
plot_gsea_enrichment <- function(gsea_result,
                                  top_n = 3,
                                  output_dir = "output",
                                  project_name = "Project") {
  if (is.null(gsea_result) || nrow(gsea_result) == 0) {
    message(">>> No GSEA data to plot.")
    return(invisible(NULL))
  }

  if (!requireNamespace("enrichplot", quietly = TRUE)) {
    stop("Please install enrichplot: BiocManager::install('enrichplot')")
  }

  ensure_output_dir(output_dir)
  res_df <- gsea_result@result

  contrast <- attr(gsea_result, "contrast")
  suffix_base <- if (!is.null(contrast)) paste0("GSEA_ES_", contrast) else "GSEA_ES"

  plots <- list()

  # Top activated (NES > 0) and suppressed (NES < 0)
  activated <- res_df[res_df$NES > 0, ]
  activated <- activated[order(activated$NES, decreasing = TRUE), ]
  suppressed <- res_df[res_df$NES < 0, ]
  suppressed <- suppressed[order(suppressed$NES), ]

  selected <- rbind(
    utils::head(activated, top_n),
    utils::head(suppressed, top_n)
  )

  if (nrow(selected) == 0) {
    message(">>> No pathways to plot.")
    return(invisible(NULL))
  }

  # Find indices of selected pathways in the result
  pathway_ids <- selected$ID
  pathway_idx <- which(gsea_result@result$ID %in% pathway_ids)

  if (length(pathway_idx) == 0) {
    message(">>> Pathway index mismatch.")
    return(invisible(NULL))
  }

  # Generate combined enrichment plot
  tryCatch({
    p <- enrichplot::gseaplot2(
      gsea_result,
      geneSetID = pathway_idx,
      pvalue_table = TRUE,
      ES_geom = "line",
      base_size = 11
    )

    suffix <- paste0(suffix_base, "_top", length(pathway_idx))

    # Save using ggsave for patchwork objects
    base_name <- file.path(output_dir,
                            paste0(project_name, "_", suffix))

    plot_height <- 4 + length(pathway_idx) * 0.3
    ggplot2::ggsave(paste0(base_name, ".pdf"), p,
                    width = 10, height = plot_height)
    ggplot2::ggsave(paste0(base_name, ".png"), p,
                    width = 10, height = plot_height, dpi = 300, bg = "white")

    # Save pathway info CSV
    utils::write.csv(selected[, c("ID", "Description", "NES",
                                   "pvalue", "p.adjust", "setSize")],
                     paste0(base_name, ".csv"), row.names = FALSE)

    message(paste0("\u2705 \u5df2\u751f\u6210: ", suffix))
    Sys.sleep(0.3)
    plots[["combined"]] <- p

  }, error = function(e) {
    message(">>> Combined enrichment plot failed: ", e$message)
  })

  # Also generate individual plots for top 2 activated + 2 suppressed
  top_individual <- rbind(
    utils::head(activated, min(2, nrow(activated))),
    utils::head(suppressed, min(2, nrow(suppressed)))
  )

  for (i in seq_len(nrow(top_individual))) {
    pw_id <- top_individual$ID[i]
    pw_name <- top_individual$Description[i]
    pw_idx <- which(gsea_result@result$ID == pw_id)

    if (length(pw_idx) == 0) next

    tryCatch({
      p_single <- enrichplot::gseaplot2(
        gsea_result,
        geneSetID = pw_idx,
        title = pw_name,
        pvalue_table = TRUE,
        ES_geom = "line",
        base_size = 11
      )

      # Clean pathway name for filename
      clean_name <- gsub("[^A-Za-z0-9_]", "_", pw_name)
      clean_name <- gsub("_+", "_", clean_name)
      clean_name <- substr(clean_name, 1, 60)
      single_suffix <- paste0(suffix_base, "_", clean_name)
      single_base <- file.path(output_dir,
                                paste0(project_name, "_", single_suffix))

      ggplot2::ggsave(paste0(single_base, ".pdf"), p_single,
                      width = 8, height = 6)
      ggplot2::ggsave(paste0(single_base, ".png"), p_single,
                      width = 8, height = 6, dpi = 300, bg = "white")

      message(paste0("\u2705 \u5df2\u751f\u6210: ", single_suffix))
      plots[[pw_id]] <- p_single

    }, error = function(e) {
      message(sprintf(">>> Plot for %s failed: %s", pw_name, e$message))
    })
  }

  invisible(plots)
}
