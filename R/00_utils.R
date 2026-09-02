options(stringsAsFactors = FALSE)

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L) y else x

abortf <- function(fmt, ...) {
  stop(sprintf(fmt, ...), call. = FALSE)
}

messagef <- function(fmt, ...) {
  message(sprintf(fmt, ...))
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  normalizePath(path, winslash = "/", mustWork = TRUE)
}

assert_file <- function(path) {
  if (!file.exists(path)) abortf("Required input file is missing: %s", path)
  invisible(path)
}

assert_columns <- function(x, columns, context = "data frame") {
  missing <- setdiff(columns, names(x))
  if (length(missing)) {
    abortf("%s is missing required columns: %s", context, paste(missing, collapse = ", "))
  }
  invisible(TRUE)
}

trim_na <- function(x) {
  x <- trimws(as.character(x))
  x[x %in% c("", "NA", "N/A", "None", "null")] <- NA_character_
  x
}

collapse_unique <- function(x, separator = ";") {
  values <- trim_na(unlist(strsplit(paste(x, collapse = separator), separator, fixed = TRUE)))
  values <- sort(unique(values[!is.na(values)]))
  if (!length(values)) NA_character_ else paste(values, collapse = separator)
}

as_numeric_matrix <- function(x, context = "matrix") {
  out <- as.matrix(data.frame(lapply(x, function(z) suppressWarnings(as.numeric(z))),
                              check.names = FALSE))
  storage.mode(out) <- "double"
  if (!ncol(out) || !nrow(out)) abortf("%s is empty", context)
  out
}

safe_neg_log10 <- function(p) {
  p <- suppressWarnings(as.numeric(p))
  -log10(pmax(p, .Machine$double.xmin, na.rm = FALSE))
}

safe_sheet_name <- function(x, used = character()) {
  x <- gsub("[\\\\/:?*\\[\\]]", "_", x)
  x <- substr(x, 1L, 31L)
  candidate <- x
  i <- 1L
  while (candidate %in% used) {
    suffix <- paste0("_", i)
    candidate <- paste0(substr(x, 1L, 31L - nchar(suffix)), suffix)
    i <- i + 1L
  }
  candidate
}

write_tsv <- function(x, path) {
  ensure_dir(dirname(path))
  if (is.data.frame(x) && ncol(x) == 0L) {
    x <- data.frame(Message = "No rows produced under the configured rule", stringsAsFactors = FALSE)
  }
  readr::write_tsv(x, path, na = "")
  invisible(path)
}

write_csv_utf8 <- function(x, path) {
  ensure_dir(dirname(path))
  if (is.data.frame(x) && ncol(x) == 0L) {
    x <- data.frame(Message = "No rows produced under the configured rule", stringsAsFactors = FALSE)
  }
  readr::write_excel_csv(x, path, na = "")
  invisible(path)
}

canonical_sample_id <- function(x) {
  x <- gsub("_", ".", as.character(x), fixed = TRUE)
  token <- sub("^[FM][LMH][HSU]\\.?", "", x)
  token <- sub("^[FM][LMH]H\\.?", "", token)
  token <- sub("^A", "A", token)
  token
}

sample_group_from_name <- function(x) {
  x <- toupper(gsub("[._]", "", as.character(x)))
  prefix <- substr(x, 1L, 2L)
  prefix[!prefix %in% c("FL", "FM", "FH", "ML", "MM", "MH")] <- NA_character_
  prefix
}

sample_metadata_from_names <- function(sample_names, layer_id) {
  group <- sample_group_from_name(sample_names)
  if (anyNA(group)) {
    abortf("Cannot infer group for samples in %s: %s", layer_id,
           paste(sample_names[is.na(group)], collapse = ", "))
  }
  data.frame(
    Sample = sample_names,
    Canonical_animal_ID = canonical_sample_id(sample_names),
    Group = group,
    Sex = ifelse(substr(group, 1L, 1L) == "F", "Female", "Male"),
    Age_code_original = substr(group, 2L, 2L),
    Age_display = unname(c(L = "E", M = "M", H = "L")[substr(group, 2L, 2L)]),
    Layer = layer_id,
    stringsAsFactors = FALSE
  )
}

hash_manifest <- function(root, output_file) {
  files <- list.files(root, recursive = TRUE, full.names = TRUE, all.files = FALSE)
  files <- files[file.info(files)$isdir %in% FALSE]
  files <- setdiff(files, output_file)
  if (!length(files)) return(invisible(NULL))
  hashes <- tools::md5sum(files)
  manifest <- data.frame(
    md5 = unname(hashes),
    file = substring(normalizePath(names(hashes), winslash = "/"), nchar(normalizePath(root, winslash = "/")) + 2L),
    stringsAsFactors = FALSE
  )
  write.table(manifest, output_file, sep = "\t", quote = FALSE, row.names = FALSE)
  invisible(manifest)
}

verify_input_md5 <- function(project_root, manifest_file = "INPUT_MD5SUMS.tsv", strict = TRUE) {
  path <- file.path(project_root, manifest_file)
  assert_file(path)
  manifest <- utils::read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
  assert_columns(manifest, c("md5", "file"), manifest_file)
  input_paths <- file.path(project_root, manifest$file)
  missing <- !file.exists(input_paths)
  actual <- rep(NA_character_, length(input_paths))
  actual[!missing] <- unname(tools::md5sum(input_paths[!missing]))
  passed <- !missing & tolower(actual) == tolower(manifest$md5)
  verification <- data.frame(File = manifest$file, Expected_MD5 = manifest$md5,
                             Actual_MD5 = actual, Passed = passed, stringsAsFactors = FALSE)
  if (any(!passed) && isTRUE(strict)) {
    abortf("Input MD5 verification failed: %s", paste(manifest$file[!passed], collapse = ", "))
  }
  verification
}

verify_target_environment <- function(target_file = "config/package_versions_target.tsv",
                                      strict = TRUE) {
  assert_file(target_file)
  targets <- utils::read.delim(target_file, check.names = FALSE, stringsAsFactors = FALSE)
  assert_columns(targets, c("Package", "Version", "Repository_or_component"), target_file)
  actual <- vapply(targets$Package, function(package) {
    if (identical(package, "R")) return(as.character(getRversion()))
    if (identical(package, "Bioconductor")) {
      if (!requireNamespace("BiocManager", quietly = TRUE)) return(NA_character_)
      return(as.character(BiocManager::version()))
    }
    if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
    as.character(utils::packageVersion(package))
  }, character(1))
  verification <- data.frame(
    Package = targets$Package,
    Expected_version = targets$Version,
    Actual_version = actual,
    Passed = !is.na(actual) & actual == targets$Version,
    stringsAsFactors = FALSE
  )
  if (any(!verification$Passed) && isTRUE(strict)) {
    failed <- sprintf(
      "%s expected %s, observed %s",
      verification$Package[!verification$Passed],
      verification$Expected_version[!verification$Passed],
      ifelse(is.na(verification$Actual_version[!verification$Passed]),
             "not installed", verification$Actual_version[!verification$Passed])
    )
    abortf("Target environment verification failed: %s", paste(failed, collapse = " | "))
  }
  verification
}

write_session_records <- function(output_dir, required_packages) {
  ensure_dir(file.path(output_dir, "reproducibility"))
  capture.output(sessionInfo(), file = file.path(output_dir, "reproducibility", "sessionInfo.txt"))
  capture.output(sessioninfo::session_info(),
                 file = file.path(output_dir, "reproducibility", "session_info_detailed.txt"))
  versions <- data.frame(
    Package = required_packages,
    Version = vapply(required_packages, function(pkg) as.character(utils::packageVersion(pkg)), character(1)),
    stringsAsFactors = FALSE
  )
  write_csv_utf8(versions, file.path(output_dir, "reproducibility", "package_versions_actual.csv"))
  versions
}

assert_no_interaction_model <- function(cfg, project_root) {
  if (!identical(cfg$analysis$allow_interaction, FALSE)) {
    abortf("Interaction modelling is forbidden by this delivery specification.")
  }
  r_files <- file.path(project_root, "R", "03_limma_pairwise.R")
  code <- paste(unlist(lapply(r_files, readLines, warn = FALSE)), collapse = "\n")
  forbidden <- c("Sex:Age", "Age:Sex", "Sex\\*Age", "Age\\*Sex", "sex_by_age_interaction")
  hit <- vapply(forbidden, function(pattern) grepl(pattern, code, perl = TRUE), logical(1))
  if (any(hit)) abortf("Forbidden interaction-model token detected: %s", paste(forbidden[hit], collapse = ", "))
  invisible(TRUE)
}
