split_semicolon_ids <- function(x) {
  ids <- trim_na(unlist(strsplit(paste(x, collapse = ";"), ";", fixed = TRUE)))
  sort(unique(ids[!is.na(ids)]))
}

ora_from_mapping <- function(query, universe, mapping, min_size = 3L, max_size = 500L,
                             fdr_cutoff = 0.05) {
  assert_columns(mapping, c("Feature_ID", "Term_ID", "Term_name"), "enrichment mapping")
  mapping$Feature_ID <- trim_na(mapping$Feature_ID)
  mapping$Term_ID <- trim_na(mapping$Term_ID)
  mapping <- unique(mapping[!is.na(mapping$Feature_ID) & !is.na(mapping$Term_ID), , drop = FALSE])
  universe <- intersect(sort(unique(universe)), unique(mapping$Feature_ID))
  query <- intersect(sort(unique(query)), universe)
  term_members <- split(mapping$Feature_ID, mapping$Term_ID)
  term_members <- lapply(term_members, function(x) intersect(unique(x), universe))
  sizes <- vapply(term_members, length, integer(1))
  eligible <- names(sizes)[sizes >= min_size & sizes <= max_size]
  term_members <- term_members[eligible]
  if (!length(term_members)) return(data.frame())

  N <- length(universe)
  n <- length(query)
  term_names <- tapply(mapping$Term_name, mapping$Term_ID, function(z) trim_na(z)[1])
  rows <- lapply(names(term_members), function(term_id) {
    members <- term_members[[term_id]]
    hits <- intersect(query, members)
    M <- length(members)
    k <- length(hits)
    p <- if (n == 0L || k == 0L) 1 else stats::phyper(k - 1L, M, N - M, n, lower.tail = FALSE)
    data.frame(
      Term_ID = term_id,
      Term_name = unname(term_names[term_id]),
      Query_hits = k,
      Query_size = n,
      Background_term_size = M,
      Background_size = N,
      P.Value = p,
      Matching_IDs = if (length(hits)) paste(sort(hits), collapse = ";") else NA_character_,
      stringsAsFactors = FALSE
    )
  })
  result <- dplyr::bind_rows(rows)
  result$adj.P.Val <- stats::p.adjust(result$P.Value, method = "BH")
  result$GeneRatio <- result$Query_hits / pmax(result$Query_size, 1L)
  result$BgRatio <- result$Background_term_size / pmax(result$Background_size, 1L)
  result$Formal_enrichment <- result$adj.P.Val < fdr_cutoff & result$Query_hits > 0L
  result[order(result$adj.P.Val, result$P.Value, -result$Query_hits), , drop = FALSE]
}

build_go_mapping <- function(symbols, ontology, cache_dir) {
  ensure_dir(cache_dir)
  cache_file <- file.path(cache_dir, paste0("GO_", ontology, "_symbol_term_mapping.tsv"))
  if (file.exists(cache_file)) return(readr::read_tsv(cache_file, show_col_types = FALSE))
  raw <- AnnotationDbi::select(
    org.Mm.eg.db::org.Mm.eg.db,
    keys = unique(symbols), keytype = "SYMBOL", columns = c("GO", "ONTOLOGY")
  )
  raw <- raw[raw$ONTOLOGY == ontology & !is.na(raw$GO), c("SYMBOL", "GO", "ONTOLOGY")]
  term_info <- AnnotationDbi::select(
    GO.db::GO.db,
    keys = unique(raw$GO), keytype = "GOID", columns = c("TERM", "ONTOLOGY")
  )
  mapping <- dplyr::left_join(raw, term_info, by = c("GO" = "GOID", "ONTOLOGY" = "ONTOLOGY"))
  mapping <- unique(data.frame(
    Feature_ID = mapping$SYMBOL,
    Term_ID = mapping$GO,
    Term_name = mapping$TERM,
    Ontology = ontology,
    stringsAsFactors = FALSE
  ))
  write_tsv(mapping, cache_file)
  mapping
}

kegg_cache_files <- function(cache_dir) {
  list(
    pathway_names = file.path(cache_dir, "KEGG_pathway_names.tsv"),
    gene_mapping = file.path(cache_dir, "KEGG_mmu_gene_pathway.tsv"),
    compound_mapping = file.path(cache_dir, "KEGG_compound_pathway.tsv")
  )
}

download_kegg_cache <- function(cache_dir) {
  ensure_dir(cache_dir)
  files <- kegg_cache_files(cache_dir)
  pathway_list <- KEGGREST::keggList("pathway", "mmu")
  pathway_names <- data.frame(
    Term_ID = sub("^path:", "", names(pathway_list)),
    Term_name = sub(" - Mus musculus \\(house mouse\\)$", "", unname(pathway_list)),
    stringsAsFactors = FALSE
  )
  pathway_names$Map_ID <- sub("^mmu", "map", pathway_names$Term_ID)
  write_tsv(pathway_names, files$pathway_names)

  gene_links <- KEGGREST::keggLink("pathway", "mmu")
  gene_mapping <- data.frame(
    Feature_ID = sub("^mmu:", "", names(gene_links)),
    Term_ID = sub("^path:", "", unname(gene_links)),
    stringsAsFactors = FALSE
  )
  gene_mapping <- dplyr::left_join(gene_mapping, pathway_names[, c("Term_ID", "Term_name")], by = "Term_ID")
  write_tsv(gene_mapping, files$gene_mapping)

  compound_links <- KEGGREST::keggLink("pathway", "compound")
  compound_mapping <- data.frame(
    Feature_ID = sub("^cpd:", "", names(compound_links)),
    Map_ID = sub("^path:", "", unname(compound_links)),
    stringsAsFactors = FALSE
  )
  compound_mapping <- dplyr::inner_join(
    compound_mapping,
    pathway_names[, c("Map_ID", "Term_ID", "Term_name")],
    by = "Map_ID"
  )
  compound_mapping <- unique(compound_mapping[, c("Feature_ID", "Term_ID", "Term_name")])
  write_tsv(compound_mapping, files$compound_mapping)
  metadata <- data.frame(
    Source = "KEGG REST via KEGGREST",
    Organism = "mmu",
    Retrieved_UTC = format(Sys.time(), tz = "UTC", usetz = TRUE),
    Pathway_terms = nrow(pathway_names),
    Gene_term_edges = nrow(gene_mapping),
    Compound_term_edges = nrow(compound_mapping),
    stringsAsFactors = FALSE
  )
  write_tsv(metadata, file.path(cache_dir, "KEGG_cache_metadata.tsv"))
  invisible(files)
}

get_kegg_mappings <- function(cache_dir, cfg) {
  files <- kegg_cache_files(cache_dir)
  if (!all(file.exists(unlist(files)))) {
    if (!isTRUE(cfg$enrichment$use_kegg_online_if_cache_missing)) {
      abortf("KEGG cache is missing and online retrieval is disabled: %s", cache_dir)
    }
    message("KEGG cache is absent; retrieving official mappings once through KEGGREST.")
    download_kegg_cache(cache_dir)
  }
  list(
    pathway_names = readr::read_tsv(files$pathway_names, show_col_types = FALSE),
    gene_mapping = readr::read_tsv(files$gene_mapping, show_col_types = FALSE),
    compound_mapping = readr::read_tsv(files$compound_mapping, show_col_types = FALSE)
  )
}

protein_symbol_vectors <- function(layer, result) {
  universe <- split_semicolon_ids(layer$feature_data$PG.Genes)
  query <- split_semicolon_ids(result$formal$PG.Genes)
  list(universe = universe, query = query)
}

protein_entrez_vectors <- function(symbol_vectors) {
  mapping <- AnnotationDbi::select(
    org.Mm.eg.db::org.Mm.eg.db,
    keys = symbol_vectors$universe, keytype = "SYMBOL", columns = "ENTREZID"
  )
  mapping <- unique(mapping[!is.na(mapping$ENTREZID), c("SYMBOL", "ENTREZID")])
  list(
    universe = unique(mapping$ENTREZID),
    query = unique(mapping$ENTREZID[mapping$SYMBOL %in% symbol_vectors$query]),
    symbol_entrez_map = mapping
  )
}

metabolite_kegg_vectors <- function(layer, result) {
  list(
    universe = split_semicolon_ids(layer$feature_data$id_kegg),
    query = split_semicolon_ids(result$formal$id_kegg)
  )
}

plot_enrichment_dot <- function(result, title, cfg) {
  result <- result[result$Query_hits > 0L, , drop = FALSE]
  result <- head(result[order(result$adj.P.Val, -result$Query_hits), , drop = FALSE],
                 cfg$enrichment$top_n_plot)
  if (!nrow(result)) {
    return(ggplot2::ggplot() + ggplot2::annotate("text", x = 0, y = 0,
      label = "No enriched term with at least one query hit") + ggplot2::theme_void() +
      ggplot2::labs(title = title))
  }
  result$Term_name <- factor(result$Term_name, levels = rev(result$Term_name))
  ggplot2::ggplot(result, ggplot2::aes(x = .data$GeneRatio, y = .data$Term_name)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$Query_hits, color = .data$adj.P.Val)) +
    ggplot2::scale_color_viridis_c(option = "C", direction = -1, trans = "reverse") +
    ggplot2::labs(title = title, x = "Query ratio", y = NULL, size = "Hits", color = "BH-FDR") +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   panel.grid.minor = ggplot2::element_blank())
}

run_enrichment_all <- function(layers, all_results, cfg) {
  if (!isTRUE(cfg$enrichment$enabled)) return(list(results = list(), summary = data.frame()))
  output_dir <- ensure_dir(file.path(cfg$paths$output_dir, "enrichment"))
  cache_dir <- ensure_dir(file.path(output_dir, "mapping_cache"))
  kegg <- tryCatch(
    get_kegg_mappings(cache_dir, cfg),
    error = function(e) {
      if (isTRUE(cfg$project$strict)) {
        abortf("KEGG mappings could not be obtained in strict mode: %s", conditionMessage(e))
      }
      warning(sprintf("KEGG enrichment skipped because mappings could not be obtained: %s", conditionMessage(e)))
      NULL
    }
  )
  results <- list()
  summary_rows <- list()

  for (layer_id in names(layers)) {
    layer <- layers[[layer_id]]
    layer_dir <- ensure_dir(file.path(output_dir, layer_id))
    for (comparison_id in names(all_results[[layer_id]])) {
      comparison_result <- all_results[[layer_id]][[comparison_id]]
      comparison_dir <- ensure_dir(file.path(layer_dir, comparison_id))
      if (layer_id == "Cardiac_Proteome") {
        symbols <- protein_symbol_vectors(layer, comparison_result)
        for (ontology in cfg$enrichment$ontology) {
          mapping <- build_go_mapping(symbols$universe, ontology, cache_dir)
          ora <- ora_from_mapping(symbols$query, symbols$universe, mapping,
                                  cfg$enrichment$minimum_term_size,
                                  cfg$enrichment$maximum_term_size,
                                  cfg$enrichment$fdr_cutoff_for_formal)
          key <- paste(layer_id, comparison_id, paste0("GO_", ontology), sep = "::")
          results[[key]] <- ora
          write_csv_utf8(ora, file.path(comparison_dir, paste0("GO_", ontology, "_full_universe.csv")))
          plot <- plot_enrichment_dot(ora, paste(layer_id, comparison_id, paste0("GO-", ontology)), cfg)
          save_plot_formats(plot, file.path(comparison_dir, paste0("GO_", ontology, "_dotplot")), cfg,
                            width = 8.2, height = 5.6)
          summary_rows[[length(summary_rows) + 1L]] <- data.frame(
            Layer = layer_id, Comparison = comparison_id, Database = paste0("GO_", ontology),
            Query_size = length(symbols$query), Eligible_terms = nrow(ora),
            Formal_terms = sum(ora$Formal_enrichment, na.rm = TRUE), stringsAsFactors = FALSE
          )
        }
        if (!is.null(kegg)) {
          entrez <- protein_entrez_vectors(symbols)
          ora <- ora_from_mapping(entrez$query, entrez$universe, kegg$gene_mapping,
                                  cfg$enrichment$minimum_term_size,
                                  cfg$enrichment$maximum_term_size,
                                  cfg$enrichment$fdr_cutoff_for_formal)
          key <- paste(layer_id, comparison_id, "KEGG", sep = "::")
          results[[key]] <- ora
          write_csv_utf8(ora, file.path(comparison_dir, "KEGG_full_universe.csv"))
          plot <- plot_enrichment_dot(ora, paste(layer_id, comparison_id, "KEGG"), cfg)
          save_plot_formats(plot, file.path(comparison_dir, "KEGG_dotplot"), cfg, 8.2, 5.6)
          summary_rows[[length(summary_rows) + 1L]] <- data.frame(
            Layer = layer_id, Comparison = comparison_id, Database = "KEGG",
            Query_size = length(entrez$query), Eligible_terms = nrow(ora),
            Formal_terms = sum(ora$Formal_enrichment, na.rm = TRUE), stringsAsFactors = FALSE
          )
        }
      } else if (!is.null(kegg)) {
        ids <- metabolite_kegg_vectors(layer, comparison_result)
        ora <- ora_from_mapping(ids$query, ids$universe, kegg$compound_mapping,
                                cfg$enrichment$minimum_term_size,
                                cfg$enrichment$maximum_term_size,
                                cfg$enrichment$fdr_cutoff_for_formal)
        key <- paste(layer_id, comparison_id, "KEGG", sep = "::")
        results[[key]] <- ora
        write_csv_utf8(ora, file.path(comparison_dir, "KEGG_full_universe.csv"))
        plot <- plot_enrichment_dot(ora, paste(layer_id, comparison_id, "KEGG"), cfg)
        save_plot_formats(plot, file.path(comparison_dir, "KEGG_dotplot"), cfg, 8.2, 5.6)
        summary_rows[[length(summary_rows) + 1L]] <- data.frame(
          Layer = layer_id, Comparison = comparison_id, Database = "KEGG",
          Query_size = length(ids$query), Eligible_terms = nrow(ora),
          Formal_terms = sum(ora$Formal_enrichment, na.rm = TRUE), stringsAsFactors = FALSE
        )
      }
    }
  }
  list(results = results, summary = dplyr::bind_rows(summary_rows))
}
