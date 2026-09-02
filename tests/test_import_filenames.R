#!/usr/bin/env Rscript
# Run from the package root: Rscript tests/test_import_filenames.R
source("config/config.R")
source("R/01_import_validate.R")
expected <- c("FL_vs_FM.result.xlsx", "FM_vs_FH.result.xlsx", "FL_vs_FH.result.xlsx",
              "ML_vs_MM.result.xlsx", "MM_vs_MH.result.xlsx", "ML_vs_MH.result.xlsx")
actual <- protein_filename_for_comparison(cfg$comparisons$source_sheet)
stopifnot(identical(actual, expected))
stopifnot(all(file.exists(file.path(cfg$paths$protein_pairwise_dir, actual))))
for (i in seq_along(expected)) {
  stopifnot(identical(protein_filename_for_comparison(cfg$comparisons$source_sheet[i]),
                      expected[i]))
}
manifest <- read.delim("INPUT_MD5SUMS.tsv", stringsAsFactors = FALSE)
stopifnot(all(unname(tools::md5sum(manifest$file)) == tolower(manifest$md5)))
all_scripts <- c("install_packages.R", "validate_setup.R", "run_all.R",
                 list.files("R", pattern = "\\.R$", full.names = TRUE),
                 "config/config.R")
invisible(lapply(all_scripts, function(path) parse(file = path)))
message("PASS: six importer paths, nine original input hashes, and R syntax.")
