# Cardiac Aging Multi-omics

Analysis code and processed quantitative data supporting the cardiac proteomic and multi-compartment metabolomic analyses of naturally long-lived mice during late-life aging.

## Overview

This repository contains the R code and processed quantitative data used for the analyses reported in our study of age- and sex-associated molecular remodeling during late-life cardiac aging.

The analyses include:

- sex-stratified pairwise differential analyses across early-aged (E, 29 months), middle-aged (M, 32 months), and late-aged (L, 35 months) mice;
- cardiac proteomic analyses;
- cardiac, serum, and urinary metabolomic analyses;
- multiple-testing-controlled differential analyses and age-by-sex interaction analyses;
- overlap and trajectory analyses;
- functional enrichment analyses;
- exploratory protein-metabolite correlation analyses;
- SIRT4-centered exploratory association analyses; and
- generation of corresponding statistical outputs and figures.

Cardiac proteomics and cardiac, serum, and urinary metabolomics were analyzed within their respective data layers. Cross-compartment comparisons and protein-metabolite association analyses were subsequently performed as exploratory analyses. The repository does not implement formal joint multi-omics integration.

## Repository structure

- `R/` – R scripts for statistical analyses and figure generation.
- `config/` – configuration files used by the analysis pipeline.
- `input/` – processed quantitative data used as input for reproducible downstream analyses.
- `outputs/` – analysis outputs and enrichment-related results.
- `tests/` – scripts and files for testing the analysis workflow.
- `validation/` – validation resources for checking reproducibility and expected outputs.
- `vendor/` – supporting resources required by the workflow.
- `CODEBOOK.md` – description of variables and data organization.
- `README_CN.md` – detailed Chinese-language instructions.
- `INPUT_MD5SUMS.tsv` – checksums for input files.
- `PACKAGE_SHA256SUMS.tsv` – checksums for files included in the analysis package.

## Requirements

The analyses were performed in R. Package dependencies and the expected software environment are documented in:

- `install_packages.R`
- `sessionInfo_expected.txt`

Required R packages can be installed by running:

```r
source("install_packages.R")
```

The local analysis environment can then be checked using:

```r
source("validate_setup.R")
```

## Running the analyses

After installing the required packages and confirming the analysis environment, the complete downstream analysis workflow can be initiated using:

```r
source("run_all.R")
```

Input file organization and variable definitions are described in `input/README_inputs.md` and `CODEBOOK.md`.

## Data availability

Processed quantitative data required to reproduce the downstream statistical analyses are included in the `input/` directory.

Raw mass spectrometry files are not included in this GitHub repository. Raw proteomic and metabolomic data will be deposited in appropriate public repositories and their accession numbers will be added here upon completion of data deposition.

## Scope

This repository is intended to reproduce the downstream statistical analyses and figure generation performed using processed quantitative proteomic and metabolomic data. It does not reproduce raw mass spectrometry data processing, peptide/protein identification, or primary metabolomic feature extraction and annotation.

## License

The code in this repository is provided under the MIT License.
