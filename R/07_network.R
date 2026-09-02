protein_symbol_matches <- function(feature_data, symbol) {
  vapply(strsplit(ifelse(is.na(feature_data$PG.Genes), "", feature_data$PG.Genes), ";", fixed = TRUE),
         function(z) symbol %in% trimws(z), logical(1))
}

aggregate_target_protein <- function(protein_layer, symbol) {
  index <- protein_symbol_matches(protein_layer$feature_data, symbol)
  if (!any(index)) return(NULL)
  block <- protein_layer$values[index, , drop = FALSE]
  values <- if (nrow(block) == 1L) as.numeric(block[1, ]) else colMeans(block, na.rm = TRUE)
  data.frame(
    Sample = colnames(block),
    Canonical_animal_ID = canonical_sample_id(colnames(block)),
    Group = sample_group_from_name(colnames(block)),
    Sex = ifelse(substr(sample_group_from_name(colnames(block)), 1, 1) == "F", "Female", "Male"),
    Protein = symbol,
    Protein_log2_abundance = values,
    stringsAsFactors = FALSE
  )
}

formal_metabolite_union <- function(layer_results, sex, cfg) {
  ids <- cfg$comparisons$comparison_id[cfg$comparisons$sex == sex]
  sort(unique(unlist(lapply(ids, function(id) layer_results[[id]]$formal$Feature_ID))))
}

correlate_target_to_metabolites <- function(protein_data, metabolite_layer, feature_ids, sex, cfg) {
  if (is.null(protein_data) || !length(feature_ids)) return(data.frame())
  sample_meta <- metabolite_layer$sample_data
  sample_meta <- sample_meta[sample_meta$Sex == sex, , drop = FALSE]
  joined <- dplyr::inner_join(
    protein_data[protein_data$Sex == sex, , drop = FALSE],
    sample_meta[, c("Sample", "Canonical_animal_ID", "Group")],
    by = c("Canonical_animal_ID", "Group"), suffix = c("_protein", "_metabolite")
  )
  if (!nrow(joined)) return(data.frame())
  feature_ids <- intersect(feature_ids, rownames(metabolite_layer$values))
  annotation <- metabolite_layer$feature_data
  annotation <- annotation[match(feature_ids, annotation$feature_id), , drop = FALSE]
  rows <- vector("list", length(feature_ids))

  for (i in seq_along(feature_ids)) {
    metabolite_values <- metabolite_layer$values[feature_ids[i], joined$Sample_metabolite]
    protein_values <- joined$Protein_log2_abundance
    ok <- is.finite(metabolite_values) & is.finite(protein_values)
    n <- sum(ok)
    estimate <- p_value <- NA_real_
    if (n >= cfg$network$minimum_complete_pairs &&
        stats::sd(metabolite_values[ok]) > 0 && stats::sd(protein_values[ok]) > 0) {
      test <- suppressWarnings(stats::cor.test(
        protein_values[ok], metabolite_values[ok],
        method = cfg$network$correlation_method, exact = FALSE
      ))
      estimate <- unname(test$estimate)
      p_value <- test$p.value
    }
    rows[[i]] <- data.frame(
      Protein = unique(protein_data$Protein)[1],
      Metabolite_feature_ID = feature_ids[i],
      Metabolite_name = annotation$name[i],
      KEGG_ID = annotation$id_kegg[i],
      Sex = sex,
      Complete_pairs = n,
      Correlation = estimate,
      P.Value = p_value,
      stringsAsFactors = FALSE
    )
  }
  out <- dplyr::bind_rows(rows)
  out$adj.P.Val <- stats::p.adjust(out$P.Value, method = "BH")
  out$Formal_edge <- is.finite(out$adj.P.Val) & out$adj.P.Val < cfg$network$fdr_cutoff &
    is.finite(out$Correlation) & abs(out$Correlation) >= cfg$network$minimum_abs_correlation
  out
}

plot_network <- function(edges, title) {
  edges <- edges[edges$Formal_edge, , drop = FALSE]
  if (!nrow(edges)) {
    return(ggplot2::ggplot() +
      ggplot2::annotate("text", x = 0, y = 0, label = "No edge passed the formal network rule") +
      ggplot2::theme_void() + ggplot2::labs(title = title))
  }
  edge_table <- data.frame(
    from = edges$Protein,
    to = ifelse(is.na(edges$Metabolite_name), edges$Metabolite_feature_ID,
                sub(";.*$", "", edges$Metabolite_name)),
    Correlation = edges$Correlation,
    stringsAsFactors = FALSE
  )
  graph <- igraph::graph_from_data_frame(edge_table, directed = FALSE)
  vertices <- igraph::as_data_frame(graph, what = "vertices")
  vertices$Type <- ifelse(vertices$name %in% unique(edges$Protein), "Protein", "Metabolite")
  igraph::V(graph)$Type <- vertices$Type[match(igraph::V(graph)$name, vertices$name)]
  set.seed(20260828L)
  ggraph::ggraph(graph, layout = "fr") +
    ggraph::geom_edge_link(ggplot2::aes(color = .data$Correlation,
                                        width = abs(.data$Correlation)), alpha = 0.75) +
    ggraph::geom_node_point(ggplot2::aes(color = .data$Type, shape = .data$Type), size = 4) +
    ggraph::geom_node_text(ggplot2::aes(label = .data$name), repel = TRUE, size = 3) +
    ggraph::scale_edge_colour_gradient2(low = "#2D6FA3", mid = "#EEEEEE", high = "#C23B33",
                                        midpoint = 0, limits = c(-1, 1)) +
    ggraph::scale_edge_width(range = c(0.4, 1.6)) +
    ggplot2::scale_color_manual(values = c(Protein = "#7B3294", Metabolite = "#008837")) +
    ggplot2::labs(title = title, color = NULL, shape = NULL,
                  edge_color = "Pearson r", edge_width = "|r|") +
    ggraph::theme_graph(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"))
}

run_network_all <- function(layers, all_results, cfg) {
  if (!isTRUE(cfg$network$enabled)) return(list(edges = data.frame(), summary = data.frame()))
  output_dir <- ensure_dir(file.path(cfg$paths$output_dir, "network"))
  protein_layer <- layers$Cardiac_Proteome
  edge_tables <- list()
  summary_rows <- list()

  for (metabolite_layer_id in setdiff(names(layers), "Cardiac_Proteome")) {
    metabolite_layer <- layers[[metabolite_layer_id]]
    for (sex in c("Female", "Male")) {
      feature_ids <- formal_metabolite_union(all_results[[metabolite_layer_id]], sex, cfg)
      for (symbol in cfg$network$protein_symbols) {
        protein_data <- aggregate_target_protein(protein_layer, symbol)
        edges <- correlate_target_to_metabolites(protein_data, metabolite_layer,
                                                 feature_ids, sex, cfg)
        if (nrow(edges)) {
          edges$Metabolite_layer <- metabolite_layer_id
          edge_tables[[length(edge_tables) + 1L]] <- edges
        }
        key <- paste(metabolite_layer_id, sex, symbol, sep = "_")
        write_csv_utf8(edges, file.path(output_dir, paste0(key, "_all_correlations.csv")))
        formal_edges <- if (nrow(edges)) edges[edges$Formal_edge, , drop = FALSE] else edges
        write_csv_utf8(formal_edges, file.path(output_dir, paste0(key, "_formal_edges.csv")))
        plot <- plot_network(edges, paste(metabolite_layer_id, sex, symbol, sep = " | "))
        save_plot_formats(plot, file.path(output_dir, paste0(key, "_network")), cfg, 7.5, 6.0)
        summary_rows[[length(summary_rows) + 1L]] <- data.frame(
          Metabolite_layer = metabolite_layer_id,
          Sex = sex,
          Protein = symbol,
          Candidate_metabolites = length(feature_ids),
          Tested_correlations = if (nrow(edges)) sum(is.finite(edges$P.Value)) else 0L,
          Formal_edges = if (nrow(edges)) sum(edges$Formal_edge, na.rm = TRUE) else 0L,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  list(edges = dplyr::bind_rows(edge_tables), summary = dplyr::bind_rows(summary_rows))
}
