classify_differential <- function(result, rule) {
  p <- result[[rule$p_column]]
  hit <- is.finite(p) & p < rule$p_cutoff &
    is.finite(result$logFC) & abs(result$logFC) >= rule$minimum_abs_log2FC
  direction <- rep("Not_significant", nrow(result))
  direction[hit & result$logFC > 0] <- "Up"
  direction[hit & result$logFC < 0] <- "Down"
  factor(direction, levels = c("Up", "Down", "Not_significant"))
}

run_limma_comparison <- function(layer, comparison, cfg) {
  first_samples <- layer$sample_data$Sample[layer$sample_data$Group == comparison$first_code]
  second_samples <- layer$sample_data$Sample[layer$sample_data$Group == comparison$second_code]
  if (!length(first_samples) || !length(second_samples)) {
    abortf("Missing samples for %s in %s", comparison$comparison_id, layer$layer_id)
  }
  selected <- c(first_samples, second_samples)
  expression <- layer$values[, selected, drop = FALSE]
  n_first <- rowSums(is.finite(expression[, first_samples, drop = FALSE]))
  n_second <- rowSums(is.finite(expression[, second_samples, drop = FALSE]))
  keep <- n_first >= cfg$analysis$minimum_non_missing_per_group &
    n_second >= cfg$analysis$minimum_non_missing_per_group
  if (!any(keep)) abortf("No testable features for %s in %s", comparison$comparison_id, layer$layer_id)

  group <- factor(c(rep("First", length(first_samples)), rep("Second", length(second_samples))),
                  levels = c("First", "Second"))
  design <- stats::model.matrix(~ 0 + group)
  colnames(design) <- c("First", "Second")
  contrast <- limma::makeContrasts(First_minus_Second = First - Second, levels = design)
  fit <- limma::lmFit(expression[keep, , drop = FALSE], design)
  fit <- limma::contrasts.fit(fit, contrast)
  fit <- limma::eBayes(
    fit,
    trend = cfg$analysis$ebayes$trend,
    robust = cfg$analysis$ebayes$robust,
    proportion = cfg$analysis$ebayes$proportion
  )
  tab <- limma::topTable(
    fit,
    coef = "First_minus_Second",
    number = Inf,
    adjust.method = "BH",
    p.value = 1,
    sort.by = "none"
  )
  tab$Feature_ID <- rownames(tab)

  all_ids <- rownames(layer$values)
  result <- data.frame(
    Feature_ID = all_ids,
    Layer = layer$layer_id,
    Comparison = comparison$comparison_id,
    Sex = comparison$sex,
    First_group = comparison$first_label,
    Second_group = comparison$second_label,
    First_group_source_code = comparison$first_code,
    Second_group_source_code = comparison$second_code,
    N_first = n_first,
    N_second = n_second,
    Mean_first_log2 = rowMeans(expression[, first_samples, drop = FALSE], na.rm = TRUE),
    Mean_second_log2 = rowMeans(expression[, second_samples, drop = FALSE], na.rm = TRUE),
    SD_first_log2 = apply(expression[, first_samples, drop = FALSE], 1, stats::sd, na.rm = TRUE),
    SD_second_log2 = apply(expression[, second_samples, drop = FALSE], 1, stats::sd, na.rm = TRUE),
    Testable = keep,
    Exclusion_reason = ifelse(keep, NA_character_,
                              sprintf("fewer than %d finite observations in at least one group",
                                      cfg$analysis$minimum_non_missing_per_group)),
    stringsAsFactors = FALSE
  )
  result$Mean_first_log2[!is.finite(result$Mean_first_log2)] <- NA_real_
  result$Mean_second_log2[!is.finite(result$Mean_second_log2)] <- NA_real_
  result$SD_first_log2[!is.finite(result$SD_first_log2)] <- NA_real_
  result$SD_second_log2[!is.finite(result$SD_second_log2)] <- NA_real_

  stat_columns <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")
  for (column in stat_columns) result[[column]] <- NA_real_
  idx <- match(tab$Feature_ID, result$Feature_ID)
  for (column in stat_columns) result[[column]][idx] <- tab[[column]]
  result$FC_first_over_second <- 2^result$logFC
  result$Formal_call <- classify_differential(result, cfg$analysis$formal)

  annotation <- layer$feature_data
  names(annotation)[names(annotation) == layer$feature_id_column] <- "Feature_ID"
  result <- dplyr::left_join(annotation, result, by = "Feature_ID")
  formal <- result[result$Formal_call %in% c("Up", "Down"), , drop = FALSE]
  formal <- formal[order(formal$adj.P.Val, -abs(formal$logFC)), , drop = FALSE]

  summary <- data.frame(
    Layer = layer$layer_id,
    Comparison = comparison$comparison_id,
    Sex = comparison$sex,
    First_group = comparison$first_label,
    Second_group = comparison$second_label,
    Total_features = nrow(result),
    Testable_features = sum(result$Testable),
    Formal_up = sum(result$Formal_call == "Up", na.rm = TRUE),
    Formal_down = sum(result$Formal_call == "Down", na.rm = TRUE),
    Formal_total = nrow(formal),
    Formal_P_column = cfg$analysis$formal$p_column,
    Formal_P_cutoff = cfg$analysis$formal$p_cutoff,
    Minimum_abs_log2FC = cfg$analysis$formal$minimum_abs_log2FC,
    Model = cfg$analysis$model,
    stringsAsFactors = FALSE
  )

  list(all = result, formal = formal, summary = summary, design = design, contrast = contrast)
}

run_layer_limma <- function(layer, cfg) {
  results <- vector("list", nrow(cfg$comparisons))
  names(results) <- cfg$comparisons$comparison_id
  for (i in seq_len(nrow(cfg$comparisons))) {
    comparison <- as.list(cfg$comparisons[i, , drop = FALSE])
    messagef("limma: %s / %s", layer$layer_id, comparison$comparison_id)
    results[[comparison$comparison_id]] <- run_limma_comparison(layer, comparison, cfg)
  }
  results
}

run_all_limma <- function(layers, cfg) {
  out <- lapply(layers, run_layer_limma, cfg = cfg)
  names(out) <- names(layers)
  out
}

collect_differential_summary <- function(all_results) {
  dplyr::bind_rows(unlist(lapply(all_results, function(layer_results) {
    lapply(layer_results, `[[`, "summary")
  }), recursive = FALSE))
}

collect_model_specifications <- function(all_results) {
  rows <- list()
  for (layer_id in names(all_results)) {
    for (comparison_id in names(all_results[[layer_id]])) {
      x <- all_results[[layer_id]][[comparison_id]]
      rows[[length(rows) + 1L]] <- data.frame(
        Layer = layer_id,
        Comparison = comparison_id,
        Design_columns = paste(colnames(x$design), collapse = ";"),
        Contrast = "First - Second",
        Interaction_terms = 0L,
        stringsAsFactors = FALSE
      )
    }
  }
  dplyr::bind_rows(rows)
}
