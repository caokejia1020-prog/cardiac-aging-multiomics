#!/usr/bin/env Rscript
source("R/08_export_excel.R")
source("config/config.R")
source("R/00_utils.R")
cfg_test <- cfg
cfg_test$layers <- cfg$layers[1, , drop = FALSE]
cfg_test$paths$output_dir <- tempfile("fuda-empty-export-")
layer <- cfg_test$layers$layer_id[1]
empty <- data.frame(Feature_ID = character(), stringsAsFactors = FALSE)
counts <- data.frame(Pattern = "000", Count = 0L)
sets <- setNames(list(list(Female = list(membership = empty, counts = counts),
                           Male = list(membership = empty, counts = counts))), layer)
comparisons <- setNames(lapply(cfg_test$comparisons$comparison_id,
                              function(id) list(all = empty)),
                         cfg_test$comparisons$comparison_id)
all_results <- setNames(list(comparisons), layer)
export_figure_source_data(all_results, sets, cfg_test)
observed <- read.csv(file.path(cfg_test$paths$output_dir, "figure_source_data/P25_membership.csv"))
stopifnot(nrow(observed) == 0L, identical(names(observed), c("Sex", "Feature_ID")))
message("PASS: zero-row female/male membership exports without adding any feature.")
