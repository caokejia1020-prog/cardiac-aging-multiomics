#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
project_root <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE))
} else {
  normalizePath(".", winslash = "/", mustWork = TRUE)
}
setwd(project_root)

source(file.path("config", "config.R"), local = FALSE)
source(file.path("R", "00_utils.R"), local = FALSE)

required_packages <- c(
  "limma", "readxl", "openxlsx", "dplyr", "tidyr", "purrr", "stringr",
  "readr", "tibble", "ggplot2", "ggrepel", "patchwork", "igraph", "ggraph",
  "scales", "svglite", "sessioninfo", "AnnotationDbi", "org.Mm.eg.db", "GO.db",
  "KEGGREST", "statmod", "BiocManager"
)
missing_packages <- required_packages[!vapply(required_packages, requireNamespace,
                                              quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages)) {
  stop(sprintf("Missing R packages: %s\nRun: Rscript install_packages.R",
               paste(missing_packages, collapse = ", ")), call. = FALSE)
}

for (script in list.files("R", pattern = "^[0-9]{2}_.+\\.R$", full.names = TRUE)) {
  if (!grepl("00_utils\\.R$", script)) source(script, local = FALSE)
}

set.seed(cfg$project$seed)
ensure_dir(cfg$paths$output_dir)
assert_no_interaction_model(cfg, project_root)
runtime_version_verification <- verify_target_environment(strict = cfg$project$strict)
ensure_dir(file.path(cfg$paths$output_dir, "QC"))
write_csv_utf8(runtime_version_verification,
               file.path(cfg$paths$output_dir, "QC", "runtime_version_verification.csv"))
input_checksum_verification <- verify_input_md5(project_root, strict = cfg$project$strict)
write_csv_utf8(input_checksum_verification,
               file.path(cfg$paths$output_dir, "QC", "input_checksum_verification.csv"))

message("[1/9] Importing and validating all four layers")
layers <- import_all_layers(cfg)

message("[2/9] Applying declared preprocessing")
layers <- preprocess_all_layers(layers, cfg)
ensure_dir(file.path(cfg$paths$output_dir, "r_objects"))
saveRDS(layers, file.path(cfg$paths$output_dir, "r_objects", "processed_layers.rds"))

message("[3/9] Running 24 sex-stratified two-group limma contrasts")
all_results <- run_all_limma(layers, cfg)
saveRDS(all_results, file.path(cfg$paths$output_dir, "r_objects", "limma_results.rds"))
model_specifications <- collect_model_specifications(all_results)

message("[4/9] Building shared/unique sets")
intersections <- build_all_intersections(all_results, cfg)
saveRDS(intersections, file.path(cfg$paths$output_dir, "r_objects", "intersections.rds"))

message("[5/9] Rendering all 28 independent main panels")
differential_manifest <- plot_all_differential_panels(all_results, cfg)
upset_manifest <- plot_all_upsets(intersections, cfg)
panel_manifest <- dplyr::bind_rows(differential_manifest, upset_manifest)
write_csv_utf8(panel_manifest, file.path(cfg$paths$output_dir, "figures", "panel_manifest.csv"))

message("[6/9] Running enrichment with the complete eligible term universe")
enrichment <- run_enrichment_all(layers, all_results, cfg)
saveRDS(enrichment, file.path(cfg$paths$output_dir, "r_objects", "enrichment_results.rds"))

message("[7/9] Running target protein-metabolite correlation networks")
network <- run_network_all(layers, all_results, cfg)
saveRDS(network, file.path(cfg$paths$output_dir, "r_objects", "network_results.rds"))

message("[8/9] Exporting full tables, formal tables, figure inputs and Excel workbooks")
package_versions <- write_session_records(cfg$paths$output_dir, required_packages)
workbooks <- export_all_results(layers, all_results, intersections, panel_manifest,
                                enrichment, network, package_versions, model_specifications, cfg)
write_csv_utf8(data.frame(Workbook = workbooks),
               file.path(cfg$paths$output_dir, "excel", "workbook_manifest.csv"))

message("[9/9] Running reproducibility and output QC gates")
qc <- run_validation(layers, all_results, panel_manifest, model_specifications, cfg)
saveRDS(list(config = cfg, qc = qc$summary),
        file.path(cfg$paths$output_dir, "r_objects", "run_summary.rds"))
hash_manifest(cfg$paths$output_dir,
              file.path(cfg$paths$output_dir, "reproducibility", "OUTPUT_MD5SUMS.tsv"))

message(sprintf("Completed successfully: %s", normalizePath(cfg$paths$output_dir, winslash = "/")))
