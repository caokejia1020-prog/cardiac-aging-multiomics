protein_filename_for_comparison <- function(source_sheet) {
  paste0(gsub("-", "_vs_", source_sheet, fixed = TRUE), ".result.xlsx")
}

read_excel_all <- function(path, sheet) {
  assert_file(path)
  out <- readxl::read_excel(path, sheet = sheet, .name_repair = "minimal")
  out <- as.data.frame(out, check.names = FALSE)
  if (ncol(out) <= 1L) {
    out <- openxlsx::read.xlsx(path, sheet = sheet, check.names = FALSE,
                              skipEmptyRows = TRUE, skipEmptyCols = FALSE)
    out <- as.data.frame(out, check.names = FALSE)
  }
  if (!nrow(out)) abortf("No data rows found in %s [%s]", path, sheet)
  out
}

append_sample_vectors <- function(sample_store, feature_ids, values, sample_names,
                                  layer_id, tolerance = 1e-10) {
  for (j in seq_along(sample_names)) {
    sample_name <- sample_names[j]
    incoming <- setNames(values[, j], feature_ids)
    if (is.null(sample_store[[sample_name]])) {
      sample_store[[sample_name]] <- incoming
      next
    }
    current <- sample_store[[sample_name]]
    common <- intersect(names(current), names(incoming))
    comparable <- common[is.finite(current[common]) & is.finite(incoming[common])]
    if (length(comparable)) {
      delta <- max(abs(current[comparable] - incoming[comparable]), na.rm = TRUE)
      if (is.finite(delta) && delta > tolerance) {
        abortf("Repeated sample values disagree in %s for %s (maximum difference %.6g)",
               layer_id, sample_name, delta)
      }
    }
    missing_names <- setdiff(names(incoming), names(current))
    if (length(missing_names)) current[missing_names] <- incoming[missing_names]
    fillable <- intersect(names(current)[!is.finite(current)], names(incoming)[is.finite(incoming)])
    if (length(fillable)) current[fillable] <- incoming[fillable]
    sample_store[[sample_name]] <- current
  }
  sample_store
}

assemble_sample_matrix <- function(sample_store, feature_ids, layer_id) {
  if (!length(sample_store)) abortf("No samples were imported for %s", layer_id)
  matrix_out <- vapply(sample_store, function(z) unname(z[feature_ids]), numeric(length(feature_ids)))
  if (is.null(dim(matrix_out))) matrix_out <- matrix(matrix_out, ncol = 1L)
  rownames(matrix_out) <- feature_ids
  colnames(matrix_out) <- names(sample_store)
  matrix_out
}

import_protein_pairwise <- function(cfg) {
  layer_id <- "Cardiac_Proteome"
  input_dir <- cfg$paths$protein_pairwise_dir
  sample_store <- list()
  annotations <- list()
  qc <- list()
  feature_order <- NULL

  for (i in seq_len(nrow(cfg$comparisons))) {
    comp <- cfg$comparisons[i, ]
    path <- file.path(input_dir, protein_filename_for_comparison(comp$source_sheet))
    dat <- read_excel_all(path, "ALL")
    assert_columns(dat, c("Accession", "PG.Genes", "PG.ProteinDescriptions"), basename(path))
    dat$Accession <- trim_na(dat$Accession)
    if (anyNA(dat$Accession) || anyDuplicated(dat$Accession)) {
      abortf("Protein feature IDs must be non-missing and unique in %s", basename(path))
    }
    sample_cols <- grep("^[FM][LMH]H[._]", names(dat), value = TRUE)
    if (length(sample_cols) != 12L) {
      abortf("Expected 12 protein sample columns in %s, found %d", basename(path), length(sample_cols))
    }
    observed_groups <- sort(unique(sample_group_from_name(sample_cols)))
    expected_groups <- sort(c(comp$first_code, comp$second_code))
    if (!identical(observed_groups, expected_groups)) {
      abortf("Group mismatch in %s: observed %s, expected %s", basename(path),
             paste(observed_groups, collapse = ","), paste(expected_groups, collapse = ","))
    }
    values <- as_numeric_matrix(dat[, sample_cols, drop = FALSE], basename(path))
    if (is.null(feature_order)) feature_order <- dat$Accession
    sample_store <- append_sample_vectors(sample_store, dat$Accession, values, sample_cols, layer_id)
    annotations[[length(annotations) + 1L]] <- dat[, c("Accession", "PG.Genes", "PG.ProteinDescriptions")]
    qc[[length(qc) + 1L]] <- data.frame(
      Layer = layer_id,
      Source_file = basename(path),
      Source_sheet = "ALL",
      Feature_rows = nrow(dat),
      Sample_columns = length(sample_cols),
      Duplicate_feature_IDs = sum(duplicated(dat$Accession)),
      Missing_quantitative_cells = sum(!is.finite(values)),
      stringsAsFactors = FALSE
    )
  }

  annotation_all <- dplyr::bind_rows(annotations)
  annotation <- annotation_all |>
    dplyr::group_by(.data$Accession) |>
    dplyr::summarise(
      PG.Genes = collapse_unique(.data$PG.Genes),
      PG.ProteinDescriptions = collapse_unique(.data$PG.ProteinDescriptions),
      .groups = "drop"
    )
  feature_order <- unique(c(feature_order, annotation$Accession))
  annotation <- annotation[match(feature_order, annotation$Accession), , drop = FALSE]
  matrix_out <- assemble_sample_matrix(sample_store, feature_order, layer_id)
  metadata <- sample_metadata_from_names(colnames(matrix_out), layer_id)
  group_rank <- match(metadata$Group, c("FL", "FM", "FH", "ML", "MM", "MH"))
  ord <- order(group_rank, metadata$Sample)
  metadata <- metadata[ord, , drop = FALSE]
  matrix_out <- matrix_out[, metadata$Sample, drop = FALSE]

  list(
    layer_id = layer_id,
    feature_id_column = "Accession",
    label_column = "PG.Genes",
    feature_data = annotation,
    sample_data = metadata,
    values = matrix_out,
    scale = "log2",
    source_qc = dplyr::bind_rows(qc)
  )
}

collapse_metabolite_sheet <- function(dat, sample_cols, context) {
  assert_columns(dat, c("peak_name", "mz", "rt", "name", "id_kegg", "Ion mode"), context)
  dat$peak_name <- trim_na(dat$peak_name)
  dat$`Ion mode` <- toupper(trim_na(dat$`Ion mode`))
  if (anyNA(dat$peak_name) || anyNA(dat$`Ion mode`)) {
    abortf("Missing peak_name or Ion mode in %s", context)
  }
  dat$feature_id <- paste(dat$peak_name, dat$`Ion mode`, sep = "_")
  values <- as_numeric_matrix(dat[, sample_cols, drop = FALSE], context)
  unique_ids <- unique(dat$feature_id)
  collapsed_values <- matrix(NA_real_, nrow = length(unique_ids), ncol = length(sample_cols),
                             dimnames = list(unique_ids, sample_cols))
  annotation_rows <- vector("list", length(unique_ids))

  for (i in seq_along(unique_ids)) {
    idx <- which(dat$feature_id == unique_ids[i])
    block <- values[idx, , drop = FALSE]
    if (length(idx) > 1L) {
      for (j in seq_len(ncol(block))) {
        finite <- block[is.finite(block[, j]), j]
        if (length(finite) > 1L && diff(range(finite)) > 1e-8) {
          abortf("Duplicate annotation rows have inconsistent abundances for %s in %s",
                 unique_ids[i], context)
        }
      }
    }
    collapsed_values[i, ] <- apply(block, 2, function(z) {
      z <- z[is.finite(z)]
      if (length(z)) z[1] else NA_real_
    })
    annotation_rows[[i]] <- data.frame(
      feature_id = unique_ids[i],
      peak_name = dat$peak_name[idx[1]],
      Ion_mode = dat$`Ion mode`[idx[1]],
      mz = suppressWarnings(as.numeric(dat$mz[idx[1]])),
      rt = suppressWarnings(as.numeric(dat$rt[idx[1]])),
      name = collapse_unique(dat$name[idx]),
      id_kegg = collapse_unique(dat$id_kegg[idx]),
      Candidate_annotation_rows = length(idx),
      stringsAsFactors = FALSE
    )
  }
  list(
    feature_data = dplyr::bind_rows(annotation_rows),
    values = collapsed_values,
    duplicate_rows_collapsed = nrow(dat) - length(unique_ids)
  )
}

import_metabolite_workbook <- function(path, layer_id, cfg) {
  sample_store <- list()
  annotations <- list()
  qc <- list()
  feature_order <- character()

  for (i in seq_len(nrow(cfg$comparisons))) {
    comp <- cfg$comparisons[i, ]
    dat <- read_excel_all(path, comp$source_sheet)
    sample_cols <- grep("^[FM][LMH][HSU][._]", names(dat), value = TRUE)
    if (length(sample_cols) != 12L) {
      abortf("Expected 12 metabolite sample columns in %s [%s], found %d",
             basename(path), comp$source_sheet, length(sample_cols))
    }
    observed_groups <- sort(unique(sample_group_from_name(sample_cols)))
    expected_groups <- sort(c(comp$first_code, comp$second_code))
    if (!identical(observed_groups, expected_groups)) {
      abortf("Group mismatch in %s [%s]", basename(path), comp$source_sheet)
    }
    collapsed <- collapse_metabolite_sheet(dat, sample_cols,
                                           sprintf("%s [%s]", basename(path), comp$source_sheet))
    fdat <- collapsed$feature_data
    values <- collapsed$values
    feature_order <- unique(c(feature_order, fdat$feature_id))
    sample_store <- append_sample_vectors(sample_store, fdat$feature_id, values,
                                          sample_cols, layer_id, tolerance = 1e-8)
    annotations[[length(annotations) + 1L]] <- fdat
    qc[[length(qc) + 1L]] <- data.frame(
      Layer = layer_id,
      Source_file = basename(path),
      Source_sheet = comp$source_sheet,
      Feature_rows = nrow(dat),
      Quantitative_features = nrow(fdat),
      Duplicate_annotation_rows_collapsed = collapsed$duplicate_rows_collapsed,
      Sample_columns = length(sample_cols),
      Missing_or_nonfinite_cells = sum(!is.finite(values)),
      Nonpositive_cells = sum(values <= 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }

  annotation_all <- dplyr::bind_rows(annotations)
  annotation <- annotation_all |>
    dplyr::group_by(.data$feature_id) |>
    dplyr::summarise(
      peak_name = dplyr::first(stats::na.omit(.data$peak_name)),
      Ion_mode = dplyr::first(stats::na.omit(.data$Ion_mode)),
      mz = dplyr::first(stats::na.omit(.data$mz)),
      rt = dplyr::first(stats::na.omit(.data$rt)),
      name = collapse_unique(.data$name),
      id_kegg = collapse_unique(.data$id_kegg),
      Candidate_annotation_rows = max(.data$Candidate_annotation_rows, na.rm = TRUE),
      .groups = "drop"
    )
  annotation <- annotation[match(feature_order, annotation$feature_id), , drop = FALSE]
  matrix_out <- assemble_sample_matrix(sample_store, feature_order, layer_id)
  metadata <- sample_metadata_from_names(colnames(matrix_out), layer_id)
  group_rank <- match(metadata$Group, c("FL", "FM", "FH", "ML", "MM", "MH"))
  ord <- order(group_rank, metadata$Sample)
  metadata <- metadata[ord, , drop = FALSE]
  matrix_out <- matrix_out[, metadata$Sample, drop = FALSE]

  list(
    layer_id = layer_id,
    feature_id_column = "feature_id",
    label_column = "name",
    feature_data = annotation,
    sample_data = metadata,
    values = matrix_out,
    scale = "linear",
    source_qc = dplyr::bind_rows(qc)
  )
}

import_all_layers <- function(cfg) {
  result <- list(Cardiac_Proteome = import_protein_pairwise(cfg))
  for (layer_id in names(cfg$paths$metabolite_files)) {
    result[[layer_id]] <- import_metabolite_workbook(
      cfg$paths$metabolite_files[[layer_id]], layer_id, cfg
    )
  }
  if (!identical(sort(names(result)), sort(cfg$layers$layer_id))) {
    abortf("Imported layers do not match configured layers")
  }
  result[cfg$layers$layer_id]
}
