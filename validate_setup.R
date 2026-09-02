#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
project_root <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE))
} else normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(project_root)

r_files <- c(
  "run_all.R", "install_packages.R", "validate_setup.R",
  list.files("R", pattern = "\\.R$", full.names = TRUE),
  "config/config.R"
)
for (file in r_files) {
  tryCatch(parse(file = file), error = function(e) {
    stop(sprintf("R syntax error in %s: %s", file, conditionMessage(e)), call. = FALSE)
  })
}

model_code <- paste(readLines("R/03_limma_pairwise.R", warn = FALSE), collapse = "\n")
forbidden <- c(
  paste0("Sex", ":", "Age"), paste0("Age", ":", "Sex"),
  paste0("Sex", "*", "Age"), paste0("Age", "*", "Sex")
)
hits <- forbidden[vapply(forbidden, function(x) grepl(x, model_code, fixed = TRUE), logical(1))]
if (length(hits)) {
  stop(sprintf("Forbidden interaction token found: %s", paste(hits, collapse = ", ")),
       call. = FALSE)
}

required_model_tokens <- c(
  "lmFit", "makeContrasts", "contrasts.fit", "eBayes", "topTable",
  'adjust.method = "BH"'
)
missing_model_tokens <- required_model_tokens[
  !vapply(required_model_tokens, function(x) grepl(x, model_code, fixed = TRUE), logical(1))
]
if (length(missing_model_tokens)) {
  stop(sprintf("Required limma token missing: %s", paste(missing_model_tokens, collapse = ", ")),
       call. = FALSE)
}

source("config/config.R", local = FALSE)
source("R/01_import_validate.R", local = FALSE)
protein_input_names <- protein_filename_for_comparison(cfg$comparisons$source_sheet)
expected_protein_names <- paste0(
  c("FL_vs_FM", "FM_vs_FH", "FL_vs_FH", "ML_vs_MM", "MM_vs_MH", "ML_vs_MH"),
  ".result.xlsx"
)
if (!identical(protein_input_names, expected_protein_names)) {
  stop("Protein importer filenames do not match the six supplied workbooks", call. = FALSE)
}

required_inputs <- c(
  file.path(
    "input/protein_pairwise",
    protein_input_names
  ),
  "input/metabolite/heart_data_unfiltered.xlsx",
  "input/metabolite/serum_data_unfiltered.xlsx",
  "input/metabolite/urine_data_unfiltered.xlsx"
)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop(sprintf("Missing input files: %s", paste(missing_inputs, collapse = ", ")),
       call. = FALSE)
}

checksum_manifest <- utils::read.delim("INPUT_MD5SUMS.tsv", check.names = FALSE,
                                       stringsAsFactors = FALSE)
if (!identical(names(checksum_manifest), c("md5", "file"))) {
  stop("INPUT_MD5SUMS.tsv must contain exactly: md5, file", call. = FALSE)
}
if (!setequal(checksum_manifest$file, required_inputs)) {
  stop("INPUT_MD5SUMS.tsv does not match the nine required inputs", call. = FALSE)
}
actual_md5 <- unname(tools::md5sum(checksum_manifest$file))
if (!all(tolower(actual_md5) == tolower(checksum_manifest$md5))) {
  bad <- checksum_manifest$file[tolower(actual_md5) != tolower(checksum_manifest$md5)]
  stop(sprintf("Input checksum mismatch: %s", paste(bad, collapse = ", ")), call. = FALSE)
}

panel_manifest <- utils::read.csv("config/panel_manifest_expected.csv", check.names = FALSE,
                                  stringsAsFactors = FALSE)
if (nrow(panel_manifest) != 28L || anyDuplicated(panel_manifest$Panel)) {
  stop("The expected panel manifest must contain 28 unique panels", call. = FALSE)
}
if (sum(panel_manifest$Panel_type == "Differential_volcano") != 24L ||
    sum(panel_manifest$Panel_type == "Shared_unique_UpSet") != 4L) {
  stop("The expected panel manifest must contain 24 differential and 4 UpSet panels",
       call. = FALSE)
}

source("R/00_utils.R", local = FALSE)
verify_target_environment(strict = TRUE)

message(sprintf(
  "Pre-run validation passed: %d R files, 9 verified inputs, 28 fixed panels, no interaction model.",
  length(r_files)
))
