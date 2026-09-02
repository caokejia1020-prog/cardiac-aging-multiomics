preprocess_layer <- function(layer, cfg) {
  values <- layer$values
  feature_data <- layer$feature_data
  before <- data.frame(
    Layer = layer$layer_id,
    Features_imported = nrow(values),
    Samples_imported = ncol(values),
    Missing_or_nonfinite_before = sum(!is.finite(values)),
    Nonpositive_before = sum(values <= 0, na.rm = TRUE),
    stringsAsFactors = FALSE
  )

  if (layer$layer_id == "Cardiac_Proteome") {
    if (!identical(layer$scale, "log2")) {
      abortf("Unsupported protein input scale: %s", layer$scale)
    }
    values[!is.finite(values)] <- NA_real_
  } else {
    values[!is.finite(values)] <- NA_real_
    if (isTRUE(cfg$preprocessing$metabolite$zero_as_missing)) values[values <= 0] <- NA_real_
    if (isTRUE(cfg$preprocessing$metabolite$log2_transform)) {
      finite_positive <- is.finite(values) & values > 0
      values[finite_positive] <- log2(values[finite_positive])
    }
  }

  if (anyDuplicated(rownames(values))) abortf("Duplicated matrix row names after preprocessing: %s", layer$layer_id)
  if (anyDuplicated(layer$sample_data$Sample)) abortf("Duplicated sample IDs: %s", layer$layer_id)
  if (!identical(colnames(values), layer$sample_data$Sample)) {
    abortf("Sample metadata and matrix columns are not in identical order: %s", layer$layer_id)
  }
  if (nrow(feature_data) != nrow(values)) abortf("Feature annotation mismatch: %s", layer$layer_id)

  after <- data.frame(
    Layer = layer$layer_id,
    Features_after_preprocessing = nrow(values),
    Samples_after_preprocessing = ncol(values),
    Missing_after = sum(!is.finite(values)),
    stringsAsFactors = FALSE
  )
  layer$feature_data <- feature_data
  layer$values <- values
  layer$scale <- "log2"
  layer$preprocessing_qc <- cbind(before, after[, -1, drop = FALSE])
  layer
}

preprocess_all_layers <- function(layers, cfg) {
  out <- lapply(layers, preprocess_layer, cfg = cfg)
  names(out) <- names(layers)
  out
}
