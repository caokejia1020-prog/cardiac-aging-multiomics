validate_sample_mapping <- function(layers) {
  rows <- list()
  for (layer_id in names(layers)) {
    table_group <- table(layers[[layer_id]]$sample_data$Group)
    groups <- c("FL", "FM", "FH", "ML", "MM", "MH")
    counts <- as.integer(table_group[groups])
    counts[is.na(counts)] <- 0L
    rows[[length(rows) + 1L]] <- data.frame(
      Check = paste0(layer_id, ": six groups each contain six samples"),
      Passed = all(counts == 6L),
      Detail = paste(paste0(groups, "=", counts), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }
  dplyr::bind_rows(rows)
}

validate_limma_tables <- function(all_results, cfg) {
  rows <- list()
  for (layer_id in names(all_results)) {
    for (comparison_id in names(all_results[[layer_id]])) {
      result <- all_results[[layer_id]][[comparison_id]]$all
      tested <- result[result$Testable & is.finite(result$P.Value), , drop = FALSE]
      recalculated_bh <- stats::p.adjust(tested$P.Value, method = "BH")
      bh_delta <- if (nrow(tested)) max(abs(recalculated_bh - tested$adj.P.Val), na.rm = TRUE) else 0
      fc_delta <- if (nrow(tested)) {
        max(abs((tested$Mean_first_log2 - tested$Mean_second_log2) - tested$logFC), na.rm = TRUE)
      } else 0
      formal_expected <- classify_differential(result, cfg$analysis$formal)
      formal_match <- identical(as.character(formal_expected), as.character(result$Formal_call))
      rows[[length(rows) + 1L]] <- data.frame(
        Check = paste(layer_id, comparison_id, sep = ": "),
        Passed = is.finite(bh_delta) && bh_delta < 1e-12 &&
          is.finite(fc_delta) && fc_delta < 1e-10 && formal_match,
        Detail = sprintf("max_BH_delta=%.3g; max_logFC_delta=%.3g; formal_match=%s",
                         bh_delta, fc_delta, formal_match),
        stringsAsFactors = FALSE
      )
    }
  }
  dplyr::bind_rows(rows)
}

validate_panel_outputs <- function(panel_manifest, cfg) {
  expected_formats <- cfg$figures$formats
  rows <- list()
  for (i in seq_len(nrow(panel_manifest))) {
    primary <- panel_manifest$Primary_file[i]
    stem <- sub("\\.png$", "", primary)
    files <- paste0(stem, ".", expected_formats)
    rows[[length(rows) + 1L]] <- data.frame(
      Check = paste0(panel_manifest$Panel[i], ": all figure formats"),
      Passed = all(file.exists(files) & file.info(files)$size > 0),
      Detail = paste(basename(files), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }
  count_check <- data.frame(
    Check = "Exactly 28 independent main panels",
    Passed = nrow(panel_manifest) == cfg$expected$number_of_main_panels &&
      !anyDuplicated(panel_manifest$Panel),
    Detail = sprintf("observed=%d; expected=%d", nrow(panel_manifest),
                     cfg$expected$number_of_main_panels),
    stringsAsFactors = FALSE
  )
  expected <- utils::read.csv("config/panel_manifest_expected.csv", check.names = FALSE,
                              stringsAsFactors = FALSE)
  observed_key <- panel_manifest[, c("Panel", "Panel_type", "Layer", "Comparison")]
  mapping_check <- data.frame(
    Check = "Panel mapping matches the locked 28-panel manifest",
    Passed = identical(observed_key, expected),
    Detail = "config/panel_manifest_expected.csv",
    stringsAsFactors = FALSE
  )
  dplyr::bind_rows(count_check, mapping_check, rows)
}

validate_model_specifications <- function(model_specifications) {
  data.frame(
    Check = "All 24 formal models contain zero interaction terms",
    Passed = nrow(model_specifications) == 24L &&
      all(model_specifications$Interaction_terms == 0L) &&
      all(model_specifications$Design_columns == "First;Second") &&
      all(model_specifications$Contrast == "First - Second"),
    Detail = sprintf("models=%d; maximum_interaction_terms=%d", nrow(model_specifications),
                     max(model_specifications$Interaction_terms)),
    stringsAsFactors = FALSE
  )
}

run_validation <- function(layers, all_results, panel_manifest, model_specifications, cfg) {
  qc <- dplyr::bind_rows(
    validate_sample_mapping(layers),
    validate_limma_tables(all_results, cfg),
    validate_panel_outputs(panel_manifest, cfg),
    validate_model_specifications(model_specifications)
  )
  qc_dir <- ensure_dir(file.path(cfg$paths$output_dir, "QC"))
  write_csv_utf8(qc, file.path(qc_dir, "QC_gate_results.csv"))
  summary <- data.frame(
    Checks = nrow(qc), Passed = sum(qc$Passed), Failed = sum(!qc$Passed),
    Overall_status = if (all(qc$Passed)) "PASS" else "FAIL",
    stringsAsFactors = FALSE
  )
  write_csv_utf8(summary, file.path(qc_dir, "QC_summary.csv"))
  if (!all(qc$Passed) && isTRUE(cfg$project$strict)) {
    failed <- qc$Check[!qc$Passed]
    abortf("QC failed: %s", paste(failed, collapse = " | "))
  }
  list(gates = qc, summary = summary)
}
