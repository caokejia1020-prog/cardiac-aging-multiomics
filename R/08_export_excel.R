excel_styles <- function() {
  list(
    title = openxlsx::createStyle(fontSize = 14, textDecoration = "bold", fontColour = "#FFFFFF",
                                  fgFill = "#17365D", halign = "left", valign = "center"),
    header = openxlsx::createStyle(textDecoration = "bold", fontColour = "#FFFFFF",
                                   fgFill = "#4472C4", halign = "center", valign = "center",
                                   border = "Bottom", borderColour = "#D9E2F3", wrapText = TRUE),
    formal = openxlsx::createStyle(fgFill = "#E2F0D9"),
    warning = openxlsx::createStyle(fgFill = "#FFF2CC"),
    integer = openxlsx::createStyle(numFmt = "0"),
    number = openxlsx::createStyle(numFmt = "0.0000"),
    pvalue = openxlsx::createStyle(numFmt = "0.000E+00")
  )
}

write_dataframe_sheet <- function(wb, sheet, data, title = NULL, styles = excel_styles()) {
  openxlsx::addWorksheet(wb, sheet, gridLines = FALSE)
  start_row <- 1L
  if (!is.null(title)) {
    openxlsx::writeData(wb, sheet, title, startRow = 1L, startCol = 1L)
    openxlsx::addStyle(wb, sheet, styles$title, rows = 1L, cols = 1L,
                       gridExpand = TRUE, stack = TRUE)
    openxlsx::mergeCells(wb, sheet, cols = 1:max(1L, ncol(data)), rows = 1L)
    openxlsx::setRowHeights(wb, sheet, rows = 1L, heights = 24)
    start_row <- 3L
  }
  if (!ncol(data)) data <- data.frame(Message = "No rows produced", stringsAsFactors = FALSE)
  openxlsx::writeData(wb, sheet, data, startRow = start_row, startCol = 1L,
                      withFilter = nrow(data) > 0L, keepNA = FALSE)
  openxlsx::addStyle(wb, sheet, styles$header, rows = start_row,
                     cols = seq_len(ncol(data)), gridExpand = TRUE, stack = TRUE)
  openxlsx::freezePane(wb, sheet, firstActiveRow = start_row + 1L, firstActiveCol = 1L)
  widths <- pmin(32, pmax(10, nchar(names(data)) + 2))
  openxlsx::setColWidths(wb, sheet, cols = seq_len(ncol(data)), widths = widths)
  if (nrow(data) > 0L) {
    pcols <- grep("(^P\\.Value$|adj\\.P\\.Val|FDR|pvalue|P_value)", names(data), ignore.case = TRUE)
    ncols <- which(vapply(data, is.numeric, logical(1)))
    if (length(ncols)) openxlsx::addStyle(wb, sheet, styles$number,
                                          rows = (start_row + 1L):(start_row + nrow(data)),
                                          cols = ncols, gridExpand = TRUE, stack = TRUE)
    if (length(pcols)) openxlsx::addStyle(wb, sheet, styles$pvalue,
                                          rows = (start_row + 1L):(start_row + nrow(data)),
                                          cols = pcols, gridExpand = TRUE, stack = TRUE)
  }
  invisible(start_row)
}

analysis_parameter_table <- function(cfg) {
  data.frame(
    Parameter = c(
      "Release", "Random seed", "Model", "Interaction model allowed", "Paired analysis",
      "eBayes trend", "eBayes robust", "eBayes proportion",
      "Formal P column", "Formal P cutoff", "Formal minimum |log2FC|",
      "Formal minimum fold change", "Minimum finite observations per group",
      "Protein input scale", "Protein additional normalization",
      "Protein additional imputation", "Protein additional transform",
      "Metabolite input scale", "Metabolite normalization",
      "Metabolite imputation", "Metabolite log2 transform", "Metabolite pseudocount",
      "Correlation method", "Network minimum |r|", "Network BH-FDR cutoff"
    ),
    Value = as.character(c(
      cfg$project$release, cfg$project$seed, cfg$analysis$model,
      cfg$analysis$allow_interaction, cfg$analysis$paired,
      cfg$analysis$ebayes$trend, cfg$analysis$ebayes$robust,
      cfg$analysis$ebayes$proportion, cfg$analysis$formal$p_column,
      cfg$analysis$formal$p_cutoff, cfg$analysis$formal$minimum_abs_log2FC,
      2^cfg$analysis$formal$minimum_abs_log2FC,
      cfg$analysis$minimum_non_missing_per_group,
      cfg$preprocessing$protein$input_scale,
      cfg$preprocessing$protein$additional_normalization,
      cfg$preprocessing$protein$additional_imputation,
      cfg$preprocessing$protein$additional_transform,
      cfg$preprocessing$metabolite$input_scale,
      cfg$preprocessing$metabolite$normalization,
      cfg$preprocessing$metabolite$imputation,
      cfg$preprocessing$metabolite$log2_transform,
      cfg$preprocessing$metabolite$pseudocount,
      cfg$network$correlation_method,
      cfg$network$minimum_abs_correlation,
      cfg$network$fdr_cutoff
    )),
    stringsAsFactors = FALSE
  )
}

export_differential_csvs <- function(all_results, cfg) {
  base <- ensure_dir(file.path(cfg$paths$output_dir, "tables", "differential"))
  for (layer_id in names(all_results)) {
    layer_dir <- ensure_dir(file.path(base, layer_id))
    for (comparison_id in names(all_results[[layer_id]])) {
      x <- all_results[[layer_id]][[comparison_id]]
      write_csv_utf8(x$all, file.path(layer_dir, paste0(comparison_id, "_all.csv")))
      write_csv_utf8(x$formal, file.path(layer_dir, paste0(comparison_id, "_formal_BH.csv")))
    }
  }
}

write_layer_workbook <- function(layer_id, layer, layer_results, cfg) {
  output_dir <- file.path(cfg$paths$output_dir, "excel")
  ensure_dir(output_dir)
  path <- file.path(output_dir, paste0(layer_id, "_limma_results.xlsx"))
  wb <- openxlsx::createWorkbook(creator = "R reproducibility workflow")
  styles <- excel_styles()
  summary <- dplyr::bind_rows(lapply(layer_results, `[[`, "summary"))
  write_dataframe_sheet(wb, "README", data.frame(
    Item = c("Layer", "Formal model", "Formal rule", "Direction"),
    Description = c(
      layer_id,
      cfg$analysis$model,
      "adj.P.Val < 0.05 and |logFC| >= log2(1.5)",
      "logFC = first displayed age group minus second displayed age group"
    ), stringsAsFactors = FALSE
  ), title = paste(layer_id, "complete limma results"), styles = styles)
  write_dataframe_sheet(wb, "Summary", summary, styles = styles)
  write_dataframe_sheet(wb, "Sample_metadata", layer$sample_data, styles = styles)
  write_dataframe_sheet(wb, "Feature_annotation", layer$feature_data, styles = styles)
  write_dataframe_sheet(wb, "Import_QC", layer$source_qc, styles = styles)
  write_dataframe_sheet(wb, "Preprocess_QC", layer$preprocessing_qc, styles = styles)
  used <- c("README", "Summary", "Sample_metadata", "Feature_annotation", "Import_QC", "Preprocess_QC")
  for (comparison_id in names(layer_results)) {
    x <- layer_results[[comparison_id]]
    all_name <- safe_sheet_name(paste0("All_", comparison_id), used)
    used <- c(used, all_name)
    formal_name <- safe_sheet_name(paste0("Formal_", comparison_id), used)
    used <- c(used, formal_name)
    write_dataframe_sheet(wb, all_name, x$all, styles = styles)
    write_dataframe_sheet(wb, formal_name, x$formal,
                          title = "Formal differential table: BH-FDR and fold-change rule",
                          styles = styles)
  }
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  path
}

write_intersection_workbook <- function(intersections, cfg) {
  output_dir <- file.path(cfg$paths$output_dir, "excel")
  ensure_dir(output_dir)
  path <- file.path(output_dir, "Shared_unique_and_UpSet_inputs.xlsx")
  wb <- openxlsx::createWorkbook(creator = "R reproducibility workflow")
  styles <- excel_styles()
  summary <- list()
  used <- character()
  for (layer_id in names(intersections)) {
    for (sex in c("Female", "Male")) {
      x <- intersections[[layer_id]][[sex]]
      counts_name <- safe_sheet_name(paste0(substr(layer_id, 1, 15), "_", sex, "_counts"), used)
      used <- c(used, counts_name)
      membership_name <- safe_sheet_name(paste0(substr(layer_id, 1, 12), "_", sex, "_members"), used)
      used <- c(used, membership_name)
      write_dataframe_sheet(wb, counts_name, x$counts, styles = styles)
      write_dataframe_sheet(wb, membership_name, x$membership, styles = styles)
      summary[[length(summary) + 1L]] <- cbind(Layer = layer_id, x$set_sizes)
    }
  }
  overview <- dplyr::bind_rows(summary)
  overview_name <- safe_sheet_name("Set_sizes", used)
  write_dataframe_sheet(wb, overview_name, overview,
                        title = "Formal differential sets used by the four UpSet panels",
                        styles = styles)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  path
}

write_enrichment_network_workbook <- function(enrichment, network, cfg) {
  output_dir <- file.path(cfg$paths$output_dir, "excel")
  ensure_dir(output_dir)
  path <- file.path(output_dir, "Enrichment_and_network_results.xlsx")
  wb <- openxlsx::createWorkbook(creator = "R reproducibility workflow")
  styles <- excel_styles()
  write_dataframe_sheet(wb, "Enrichment_summary", enrichment$summary, styles = styles)
  all_enrichment <- dplyr::bind_rows(lapply(names(enrichment$results), function(key) {
    x <- enrichment$results[[key]]
    if (!nrow(x)) return(NULL)
    cbind(Analysis = key, x)
  }))
  formal_enrichment <- if (nrow(all_enrichment)) {
    all_enrichment[all_enrichment$Formal_enrichment, , drop = FALSE]
  } else all_enrichment
  write_dataframe_sheet(wb, "All_enrichment_terms", all_enrichment, styles = styles)
  write_dataframe_sheet(wb, "Formal_enrichment", formal_enrichment, styles = styles)
  write_dataframe_sheet(wb, "Network_summary", network$summary, styles = styles)
  formal_edges <- if (nrow(network$edges)) network$edges[network$edges$Formal_edge, , drop = FALSE] else network$edges
  write_dataframe_sheet(wb, "All_network_tests", network$edges, styles = styles)
  write_dataframe_sheet(wb, "Formal_network_edges", formal_edges, styles = styles)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  path
}

export_figure_source_data <- function(all_results, intersections, cfg) {
  base <- ensure_dir(file.path(cfg$paths$output_dir, "figure_source_data"))
  panel_no <- 0L
  for (layer_id in cfg$layers$layer_id) {
    for (comparison_id in cfg$comparisons$comparison_id) {
      panel_no <- panel_no + 1L
      data <- all_results[[layer_id]][[comparison_id]]$all
      write_csv_utf8(data, file.path(base, sprintf("P%02d_source_data.csv", panel_no)))
    }
  }
  for (i in seq_along(cfg$layers$layer_id)) {
    panel_no <- cfg$expected$number_of_differential_panels + i
    layer_id <- cfg$layers$layer_id[i]
    counts <- dplyr::bind_rows(
      intersections[[layer_id]]$Female$counts,
      intersections[[layer_id]]$Male$counts
    )
    membership <- dplyr::bind_rows(
      cbind(Sex = rep("Female", nrow(intersections[[layer_id]]$Female$membership)),
            intersections[[layer_id]]$Female$membership),
      cbind(Sex = rep("Male", nrow(intersections[[layer_id]]$Male$membership)),
            intersections[[layer_id]]$Male$membership)
    )
    write_csv_utf8(counts, file.path(base, sprintf("P%02d_intersection_counts.csv", panel_no)))
    write_csv_utf8(membership, file.path(base, sprintf("P%02d_membership.csv", panel_no)))
  }
}

write_master_workbook <- function(layers, all_results, intersections, panel_manifest,
                                  enrichment, network, package_versions,
                                  model_specifications, cfg) {
  output_dir <- file.path(cfg$paths$output_dir, "excel")
  ensure_dir(output_dir)
  path <- file.path(output_dir, "Fuda_28Panel_complete_results.xlsx")
  wb <- openxlsx::createWorkbook(creator = "R reproducibility workflow")
  styles <- excel_styles()
  write_dataframe_sheet(wb, "README", data.frame(
    Item = c("Run entry", "Formal model", "Interaction analysis", "Formal criterion",
             "Main panel count", "Full tables", "Formal tables", "Reproducibility records"),
    Value = c("Rscript run_all.R", cfg$analysis$model, "Not performed",
              "BH-FDR < 0.05 and |log2FC| >= log2(1.5)",
              cfg$expected$number_of_main_panels,
              "Layer workbooks and outputs/tables/differential/*_all.csv",
              "Layer workbooks and outputs/tables/differential/*_formal_BH.csv",
              "outputs/reproducibility/sessionInfo.txt and package_versions_actual.csv"),
    stringsAsFactors = FALSE
  ), title = "Fuda 28-panel complete R analysis delivery", styles = styles)
  write_dataframe_sheet(wb, "Parameters", analysis_parameter_table(cfg), styles = styles)
  write_dataframe_sheet(wb, "Comparisons", cfg$comparisons, styles = styles)
  write_dataframe_sheet(wb, "Layers", cfg$layers, styles = styles)
  write_dataframe_sheet(wb, "Differential_summary", collect_differential_summary(all_results), styles = styles)
  write_dataframe_sheet(wb, "Panel_manifest", panel_manifest, styles = styles)
  write_dataframe_sheet(wb, "Model_specification", model_specifications, styles = styles)
  write_dataframe_sheet(wb, "Enrichment_summary", enrichment$summary, styles = styles)
  write_dataframe_sheet(wb, "Network_summary", network$summary, styles = styles)
  target_versions <- utils::read.delim("config/package_versions_target.tsv", check.names = FALSE,
                                       stringsAsFactors = FALSE)
  write_dataframe_sheet(wb, "Package_versions_target", target_versions, styles = styles)
  write_dataframe_sheet(wb, "Package_versions", package_versions, styles = styles)
  import_qc <- dplyr::bind_rows(lapply(layers, `[[`, "source_qc"))
  preprocessing_qc <- dplyr::bind_rows(lapply(layers, `[[`, "preprocessing_qc"))
  write_dataframe_sheet(wb, "Import_QC", import_qc, styles = styles)
  write_dataframe_sheet(wb, "Preprocess_QC", preprocessing_qc, styles = styles)
  openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  path
}

export_all_results <- function(layers, all_results, intersections, panel_manifest,
                               enrichment, network, package_versions,
                               model_specifications, cfg) {
  export_differential_csvs(all_results, cfg)
  export_figure_source_data(all_results, intersections, cfg)
  layer_books <- vapply(names(layers), function(layer_id) {
    write_layer_workbook(layer_id, layers[[layer_id]], all_results[[layer_id]], cfg)
  }, character(1))
  intersection_book <- write_intersection_workbook(intersections, cfg)
  enrichment_book <- write_enrichment_network_workbook(enrichment, network, cfg)
  master_book <- write_master_workbook(layers, all_results, intersections, panel_manifest,
                                       enrichment, network, package_versions,
                                       model_specifications, cfg)
  c(layer_books, intersection_book, enrichment_book, master_book)
}
