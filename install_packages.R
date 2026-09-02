#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
project_root <- if (length(script_arg)) {
  dirname(normalizePath(sub("^--file=", "", script_arg[1]), winslash = "/", mustWork = TRUE))
} else normalizePath(".", winslash = "/", mustWork = TRUE)
setwd(project_root)

target_file <- "config/package_versions_target.tsv"
if (!file.exists(target_file)) stop(sprintf("Missing version file: %s", target_file), call. = FALSE)
targets <- utils::read.delim(target_file, check.names = FALSE, stringsAsFactors = FALSE)
required_columns <- c("Package", "Version", "Repository_or_component")
if (!all(required_columns %in% names(targets))) {
  stop(sprintf("%s must contain: %s", target_file, paste(required_columns, collapse = ", ")),
       call. = FALSE)
}

target_version <- function(package) {
  value <- targets$Version[targets$Package == package]
  if (length(value) != 1L) stop(sprintf("Missing or duplicated target: %s", package), call. = FALSE)
  value
}

r_target <- target_version("R")
if (!identical(as.character(getRversion()), r_target)) {
  stop(sprintf("R %s is required; current R is %s", r_target, getRversion()), call. = FALSE)
}

options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 1200)
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", dependencies = c("Depends", "Imports", "LinkingTo"))
}

cran_targets <- targets[targets$Repository_or_component == "CRAN", , drop = FALSE]
for (i in seq_len(nrow(cran_targets))) {
  package <- cran_targets$Package[i]
  version <- cran_targets$Version[i]
  installed <- requireNamespace(package, quietly = TRUE) &&
    identical(as.character(utils::packageVersion(package)), version)
  if (!installed) {
    message(sprintf("Installing %s %s from CRAN archive", package, version))
    remotes::install_version(
      package = package,
      version = version,
      repos = getOption("repos")[["CRAN"]],
      dependencies = c("Depends", "Imports", "LinkingTo"),
      upgrade = "never",
      quiet = FALSE
    )
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
bioc_target <- target_version("Bioconductor")
if (!identical(as.character(BiocManager::version()), bioc_target)) {
  BiocManager::install(version = bioc_target, ask = FALSE, update = FALSE)
}

# A Bioconductor release does not pin each package's patch version.
# Install the archived limma source explicitly; do not substitute 3.62.2.
bioc_targets <- targets[targets$Repository_or_component == "Bioconductor" &
                         targets$Package != "Bioconductor", , drop = FALSE]
for (i in seq_len(nrow(bioc_targets))) {
  package <- bioc_targets$Package[i]
  version <- bioc_targets$Version[i]
  installed <- requireNamespace(package, quietly = TRUE) &&
    identical(as.character(utils::packageVersion(package)), version)
  if (!installed) {
    message(sprintf("Installing %s %s from Bioconductor %s", package, version, bioc_target))
    if (identical(package, "limma")) {
      source_url <- sprintf(
        "https://bioconductor.org/packages/%s/bioc/src/contrib/Archive/limma/limma_%s.tar.gz",
        bioc_target, version
      )
      archive <- file.path("vendor", sprintf("limma_%s.tar.gz", version))
      if (!file.exists(archive)) {
        archive <- tempfile(pattern = "limma_", fileext = ".tar.gz")
        utils::download.file(source_url, archive, mode = "wb", quiet = FALSE)
      }
      expected_md5 <- "60f71c513c6724401b8d114b8d669e08"
      if (!identical(unname(tools::md5sum(archive)), expected_md5)) {
        stop("Archived limma source checksum mismatch", call. = FALSE)
      }
      utils::install.packages(archive, repos = NULL, type = "source")
    } else {
      BiocManager::install(package, version = bioc_target, ask = FALSE, update = FALSE)
    }
  }
}

observed <- vapply(targets$Package, function(package) {
  if (identical(package, "R")) return(as.character(getRversion()))
  if (identical(package, "Bioconductor")) return(as.character(BiocManager::version()))
  if (!requireNamespace(package, quietly = TRUE)) return(NA_character_)
  as.character(utils::packageVersion(package))
}, character(1))
passed <- !is.na(observed) & observed == targets$Version
if (any(!passed)) {
  failed <- sprintf(
    "%s expected %s, observed %s",
    targets$Package[!passed], targets$Version[!passed],
    ifelse(is.na(observed[!passed]), "not installed", observed[!passed])
  )
  stop(sprintf("Package installation did not reach the target environment: %s",
               paste(failed, collapse = " | ")), call. = FALSE)
}

message("Target package environment installed. Run: Rscript validate_setup.R")
