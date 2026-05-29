# ==============================================================================
# msProteomiX — PPI (Protein-Protein Interaction) Network Analysis
# ==============================================================================

#' Run protein network analysis
#'
#' Two modes available:
#' \itemize{
#'   \item mode = "correlation" (default, offline): Build co-expression network
#'     from protein abundance correlation. No internet needed.
#'   \item mode = "string" (requires internet): Query STRING database for
#'     protein-protein interactions.
#' }
#'
#' @param diff_result Output from run_diff_analysis()
#' @param ms_data MsDataSet object (required for mode="correlation")
#' @param mode Analysis mode: "correlation" (offline) or "string" (online)
#' @param cor_threshold Correlation threshold for mode="correlation" (default 0.8)
#' @param species STRING species NCBI taxonomy ID (default 9606 for human)
#' @param score_threshold Combined score threshold 0-1000 (default 700)
#' @param network_type Network type: "full" or "physical" (default "full")
#' @return List with components: nodes, edges, n_mapped, n_edges
#' @export
run_ppi_network <- function(diff_result,
                            ms_data = NULL,
                            mode = "correlation",
                            cor_threshold = 0.8,
                            species = 9606,
                            score_threshold = 700,
                            network_type = "full") {

  # --- Offline: correlation-based co-expression network ---
  if (mode == "correlation") {
    return(.run_correlation_network(diff_result, ms_data, cor_threshold))
  }

  # --- Online: STRINGdb ---
  if (!requireNamespace("STRINGdb", quietly = TRUE)) {
    stop("Please install STRINGdb: BiocManager::install('STRINGdb')")
  }

  # Extract DE genes
  diff_df <- diff_result[diff_result$diff %in% c("UP", "DOWN"), ]
  if (nrow(diff_df) == 0) {
    message(">>> No differentially expressed proteins for PPI analysis.")
    return(NULL)
  }

  gene_col <- if ("Label_Name" %in% colnames(diff_df)) "Label_Name"
              else if ("Gene" %in% colnames(diff_df)) "Gene"
              else NULL

  if (is.null(gene_col)) {
    stop("Cannot find gene column in diff_result.")
  }

  genes <- as.character(diff_df[[gene_col]])
  genes_clean <- unique(trimws(unlist(strsplit(genes, "[;|,]"))))
  genes_clean <- genes_clean[genes_clean != "" & !is.na(genes_clean)]

  if (length(genes_clean) < 3) {
    message(">>> Too few genes (<3) for PPI analysis.")
    return(NULL)
  }

  message(sprintf(">>> Querying STRING database (species=%d, score>=%d)...",
                  species, score_threshold))
  message(sprintf("    Input: %d unique gene symbols", length(genes_clean)))

  # Initialize STRINGdb
  string_db <- tryCatch({
    STRINGdb::STRINGdb$new(
      version      = "12.0",
      species      = species,
      score_threshold = score_threshold,
      network_type = network_type
    )
  }, error = function(e) {
    message(">>> Failed to connect to STRING database: ", e$message)
    message("    Check your internet connection and try again.")
    return(NULL)
  })

  if (is.null(string_db)) return(NULL)

  # Map gene symbols to STRING IDs
  gene_df <- data.frame(gene = genes_clean, stringsAsFactors = FALSE)
  mapped <- tryCatch({
    string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)
  }, error = function(e) {
    message(">>> Gene mapping failed: ", e$message)
    return(NULL)
  })

  if (is.null(mapped) || nrow(mapped) == 0) {
    message(">>> No genes could be mapped to STRING IDs.")
    return(NULL)
  }

  message(sprintf("    Mapped: %d / %d genes", nrow(mapped), length(genes_clean)))

  # Get interactions
  interactions <- tryCatch({
    string_db$get_interactions(mapped$STRING_id)
  }, error = function(e) {
    message(">>> Failed to retrieve interactions: ", e$message)
    return(NULL)
  })

  if (is.null(interactions) || nrow(interactions) == 0) {
    message(">>> No interactions found above the score threshold.")
    return(NULL)
  }

  message(sprintf("    Interactions found: %d edges", nrow(interactions)))

  # Build node information
  # Map STRING IDs back to gene symbols
  id_to_gene <- stats::setNames(mapped$gene, mapped$STRING_id)

  nodes <- data.frame(
    gene      = mapped$gene,
    STRING_id = mapped$STRING_id,
    stringsAsFactors = FALSE
  )

  # Add diff status
  diff_status <- stats::setNames(as.character(diff_df$diff), diff_df[[gene_col]])
  diff_fc <- stats::setNames(diff_df$logFC, diff_df[[gene_col]])

  nodes$diff <- sapply(nodes$gene, function(g) {
    if (g %in% names(diff_status)) diff_status[g] else "NO"
  })
  nodes$logFC <- sapply(nodes$gene, function(g) {
    if (g %in% names(diff_fc)) diff_fc[g] else 0
  })

  # Build edges with gene names
  edges <- data.frame(
    from           = id_to_gene[interactions$from],
    to             = id_to_gene[interactions$to],
    combined_score = interactions$combined_score,
    stringsAsFactors = FALSE
  )
  # Remove edges with unmapped nodes
  edges <- edges[!is.na(edges$from) & !is.na(edges$to), ]

  list(
    nodes     = nodes,
    edges     = edges,
    string_db = string_db,
    mapped    = mapped,
    n_mapped  = nrow(mapped),
    n_edges   = nrow(edges)
  )
}


#' Plot PPI network using igraph
#'
#' Visualize protein-protein interaction network with nodes colored by
#' differential expression status and sized by connectivity degree.
#'
#' @param ppi_result Output from run_ppi_network()
#' @param diff_result Output from run_diff_analysis() (for coloring)
#' @param layout Network layout algorithm (default "fruchterman.reingold")
#' @param vertex_label_size Label font size (default 0.7)
#' @param output_dir Output directory
#' @param project_name Project name
#' @return igraph object (invisible)
#' @export
plot_ppi_network <- function(ppi_result,
                             diff_result = NULL,
                             layout = "fruchterman.reingold",
                             vertex_label_size = 0.7,
                             output_dir = "output",
                             project_name = "Project") {
  ensure_output_dir(output_dir)

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("Please install igraph: install.packages('igraph')")
  }

  if (is.null(ppi_result) || ppi_result$n_edges == 0) {
    message(">>> No PPI data to plot.")
    return(invisible(NULL))
  }

  edges <- ppi_result$edges
  nodes <- ppi_result$nodes

  # Build igraph object
  g <- igraph::graph_from_data_frame(
    d        = edges[, c("from", "to"), drop = FALSE],
    directed = FALSE,
    vertices = nodes[, c("gene"), drop = FALSE]
  )

  # Remove isolated vertices (no connections)
  isolated <- igraph::degree(g) == 0
  if (any(isolated)) {
    g <- igraph::delete_vertices(g, which(isolated))
    nodes <- nodes[!isolated, , drop = FALSE]
  }

  if (igraph::vcount(g) == 0) {
    message(">>> All nodes are isolated. No network to plot.")
    return(invisible(NULL))
  }

  # Node colors based on diff status
  node_colors <- ifelse(nodes$diff == "UP", "#E64B35",
                        ifelse(nodes$diff == "DOWN", "#3C5488", "#B8B8B8"))

  # Node size based on degree (scaled)
  deg <- igraph::degree(g)
  node_sizes <- 3 + (deg / max(deg, 1)) * 8

  # Edge width scaled by score
  edge_widths <- 0.3 + (edges$combined_score / 1000) * 1.5

  # Layout
  set.seed(42)
  if (layout == "fruchterman.reingold") {
    lo <- igraph::layout_with_fr(g)
  } else if (layout == "kamada.kawai") {
    lo <- igraph::layout_with_kk(g)
  } else if (layout == "circle") {
    lo <- igraph::layout_in_circle(g)
  } else {
    lo <- igraph::layout_with_fr(g)
  }

  contrast <- attr(diff_result, "contrast") %||% "PPI"
  fname_base <- file.path(output_dir,
                          paste0(project_name, "_PPI_Network_", contrast))

  # --- Save igraph network plot ---
  grDevices::pdf(paste0(fname_base, ".pdf"), width = 10, height = 10)

  graphics::par(mar = c(1, 1, 3, 1))
  igraph::plot.igraph(
    g,
    layout           = lo,
    vertex.color     = node_colors,
    vertex.size      = node_sizes,
    vertex.label     = igraph::V(g)$name,
    vertex.label.cex = vertex_label_size,
    vertex.label.color = "black",
    vertex.frame.color = "grey50",
    edge.color       = "grey70",
    edge.width       = edge_widths,
    main             = paste("PPI Network:", contrast)
  )

  # Legend
  graphics::legend("bottomleft",
                   legend = c("Up-regulated", "Down-regulated", "Not significant"),
                   fill   = c("#E64B35", "#3C5488", "#B8B8B8"),
                   bty    = "n", cex = 0.8)

  grDevices::dev.off()

  # --- Save PNG (for HTML report) ---
  grDevices::png(paste0(fname_base, ".png"),
                 width = 10, height = 10, units = "in", res = 300)
  graphics::par(mar = c(1, 1, 3, 1))
  igraph::plot.igraph(
    g,
    layout           = lo,
    vertex.color     = node_colors,
    vertex.size      = node_sizes,
    vertex.label     = igraph::V(g)$name,
    vertex.label.cex = vertex_label_size,
    vertex.label.color = "black",
    vertex.frame.color = "grey50",
    edge.color       = "grey70",
    edge.width       = edge_widths,
    main             = paste("Protein Network:", contrast)
  )
  graphics::legend("bottomleft",
                   legend = c("Up-regulated", "Down-regulated", "Not significant"),
                   fill   = c("#E64B35", "#3C5488", "#B8B8B8"),
                   bty    = "n", cex = 0.8)
  grDevices::dev.off()

  # --- Save edge list as CSV ---
  utils::write.csv(edges, paste0(fname_base, "_edges.csv"), row.names = FALSE)

  # --- Save node info as CSV ---
  node_out <- nodes
  node_out$degree <- deg[match(nodes$gene, names(deg))]
  utils::write.csv(node_out, paste0(fname_base, "_nodes.csv"), row.names = FALSE)

  # --- Print to screen ---
  graphics::par(mar = c(1, 1, 3, 1))
  igraph::plot.igraph(
    g,
    layout           = lo,
    vertex.color     = node_colors,
    vertex.size      = node_sizes,
    vertex.label     = igraph::V(g)$name,
    vertex.label.cex = vertex_label_size,
    vertex.label.color = "black",
    vertex.frame.color = "grey50",
    edge.color       = "grey70",
    edge.width       = edge_widths,
    main             = paste("PPI Network:", contrast)
  )
  graphics::legend("bottomleft",
                   legend = c("Up-regulated", "Down-regulated", "Not significant"),
                   fill   = c("#E64B35", "#3C5488", "#B8B8B8"),
                   bty    = "n", cex = 0.8)

  message(sprintf(">>> PPI network saved: %s", basename(fname_base)))
  message(sprintf("    Nodes: %d, Edges: %d",
                  igraph::vcount(g), igraph::ecount(g)))

  invisible(g)
}


# ==============================================================================
# Offline correlation-based co-expression network
# ==============================================================================

#' Build co-expression network from protein abundance correlation (internal)
#' @keywords internal
.run_correlation_network <- function(diff_result, ms_data, cor_threshold = 0.8) {
  if (is.null(ms_data)) {
    stop("ms_data is required for mode='correlation'. ",
         "Pass the MsDataSet object to run_ppi_network().")
  }
  stopifnot(inherits(ms_data, "MsDataSet"))

  # Extract DE genes
  diff_df <- diff_result[diff_result$diff %in% c("UP", "DOWN"), ]
  if (nrow(diff_df) == 0) {
    message(">>> No differentially expressed proteins for network analysis.")
    return(NULL)
  }

  gene_col <- if ("Gene" %in% colnames(diff_df)) "Gene"
              else if ("Label_Name" %in% colnames(diff_df)) "Label_Name"
              else NULL
  if (is.null(gene_col)) stop("Cannot find gene column in diff_result.")

  genes <- as.character(diff_df[[gene_col]])
  genes_clean <- unique(trimws(sub(";.*$", "", genes)))
  genes_clean <- genes_clean[genes_clean != "" & !is.na(genes_clean)]

  if (length(genes_clean) < 3) {
    message(">>> Too few DE genes (<3) for network analysis.")
    return(NULL)
  }

  # Get expression matrix for DE proteins
  protein_genes <- sub(";.*$", "", ms_data$protein_info$Gene)
  idx <- which(protein_genes %in% genes_clean)

  if (length(idx) < 3) {
    message(">>> Too few proteins matched in expression matrix.")
    return(NULL)
  }

  # Limit to manageable size
  if (length(idx) > 200) {
    # Take top 200 by |logFC|
    match_idx <- match(protein_genes[idx], genes_clean)
    fc_vals <- diff_df$logFC[match(protein_genes[idx], sub(";.*$", "", diff_df[[gene_col]]))]
    fc_vals[is.na(fc_vals)] <- 0
    fc_order <- order(abs(fc_vals), decreasing = TRUE)
    idx <- idx[fc_order[seq_len(200)]]
  }

  abund <- as.matrix(ms_data$proteins[idx, , drop = FALSE])
  rownames(abund) <- protein_genes[idx]

  # Remove rows with all NA
  valid_rows <- rowSums(!is.na(abund)) >= 3
  abund <- abund[valid_rows, , drop = FALSE]

  if (nrow(abund) < 3) {
    message(">>> Too few proteins with sufficient data for correlation.")
    return(NULL)
  }

  message(sprintf(">>> Building co-expression network (offline, %d proteins, |r| >= %.2f)...",
                  nrow(abund), cor_threshold))

  # Compute pairwise Pearson correlation
  cor_mat <- stats::cor(t(abund), use = "pairwise.complete.obs", method = "pearson")
  diag(cor_mat) <- 0

  # Build edges from high correlations
  edge_list <- which(abs(cor_mat) >= cor_threshold & upper.tri(cor_mat), arr.ind = TRUE)

  if (nrow(edge_list) == 0) {
    message(sprintf(">>> No edges above threshold |r| >= %.2f. Try lowering cor_threshold.", cor_threshold))
    return(NULL)
  }

  gene_names <- rownames(cor_mat)
  edges <- data.frame(
    from           = gene_names[edge_list[, 1]],
    to             = gene_names[edge_list[, 2]],
    combined_score = round(abs(cor_mat[edge_list]) * 1000),
    correlation    = round(cor_mat[edge_list], 3),
    stringsAsFactors = FALSE
  )

  # Build nodes
  all_nodes <- unique(c(edges$from, edges$to))
  diff_status <- stats::setNames(as.character(diff_df$diff), sub(";.*$", "", diff_df[[gene_col]]))
  diff_fc <- stats::setNames(diff_df$logFC, sub(";.*$", "", diff_df[[gene_col]]))

  nodes <- data.frame(
    gene      = all_nodes,
    STRING_id = all_nodes,  # compatibility with plot_ppi_network
    diff      = sapply(all_nodes, function(g) if (g %in% names(diff_status)) diff_status[g] else "NO"),
    logFC     = sapply(all_nodes, function(g) if (g %in% names(diff_fc)) diff_fc[g] else 0),
    stringsAsFactors = FALSE
  )

  message(sprintf(">>> Co-expression network: %d nodes, %d edges (offline)",
                  nrow(nodes), nrow(edges)))

  list(
    nodes     = nodes,
    edges     = edges,
    string_db = NULL,
    mapped    = nodes,
    n_mapped  = nrow(nodes),
    n_edges   = nrow(edges)
  )
}
